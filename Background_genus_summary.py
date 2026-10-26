#!/usr/bin/env python3
"""
Summarise samtools-coverage style TSVs annotated with GTDB taxonomy, to
profile the composition of a background mock community. Species-level
output is a simple raw-count summary; GENUS-level output is expressed as
relative abundance against the total reads/bases that mapped to ANY genus
(target + non-target) in that sample's own coverage TSV - not against the
raw nanostat sequencing total, and not just relative to the target list.

NB on input columns: in this pipeline's TSVs the taxonomy column is
literally headed ``domain``, but it actually holds the *whole*
semicolon-joined GTDB lineage (e.g.
``Bacteria;p__Bacillota_I;c__Bacilli_A;o__...;f__...;g__...;s__...``),
not just the domain-level tag. Species/genus are pulled out of that
string by searching for the ``s__``/``g__`` tag, not by position - the
lineage strings here start with a bare ``Bacteria``, no leading ``d__``,
so position-based splitting isn't safe. Pass --taxonomy-col to point at
a different column if the upstream format ever changes.

For every matching coverage TSV in a directory it will:

  1. Drop every row with ``numreads == 0``  (typically ~99.9% of rows).
  2. Sum ``numreads`` and ``covbases`` per target species (see
     TARGET_SPECIES - 14 species across 9 genera) purely as raw counts -
     no percentage column, this level isn't part of the relative-
     abundance feature below. GTDB alphabet suffixes are collapsed by
     default, i.e. ``Bacillus cereus_A`` and ``Bacillus cereus_B`` are
     both counted as ``Bacillus cereus``.
  3. Independently, sum ``numreads`` per target genus (see TARGET_GENERA),
     read straight off the GTDB ``g__`` rank. This is NOT simply the
     species totals rolled up: it also picks up reads from any off-target
     species within a target genus (e.g. an un-listed Listeria species,
     or a row that only resolved to genus level), so it's a genuine
     genus-level community profile.
  4. Also per target genus, sum a DEPTH-WEIGHTED bases figure:
     ``sum(meandepth * (endpos - startpos + 1))`` across every contig
     assigned to that genus. This is deliberately NOT ``covbases``:
     ``covbases`` is breadth of coverage (positions covered >=1x), capped
     at the reference length, so on ONT-length reads a handful of long
     reads can blanket a small genome at >=1x almost immediately - after
     that, 10x more data looks identical to 10,000x more data. depth *
     length reconstructs actual mapped bases (samtools coverage's default
     startpos/endpos already span the whole reference contig, so this is
     just depth * contig length) and isn't sensitive to individual read
     length the way covbases is.
  5. Looks up, per sample (basename + Timepoint), the total reads and
     total bases reported by a separate nanostat TSV (--nanostat; requires
     basename, Timepoint, num_of_reads and num_of_bases columns - any
     other columns, e.g. treatment, are read but not currently propagated
     into the output). These are carried into qc_per_file.tsv for
     reference/context only - they are NOT used as the relative-abundance
     denominator (see point 6).
  6. Expresses every target genus as relative abundance (%) against the
     total reads/bases that mapped to ANY genus (target + non-target) in
     that sample's own coverage TSV:
       relative_abundance_reads = genus total_reads        / reads_all_taxa
       relative_abundance_bases = genus total_mapped_bases / bases_all_taxa
     where reads_all_taxa / bases_all_taxa are that sample's summed reads
     / depth-weighted mapped bases across every genus row in the coverage
     TSV (target and non-target combined) - NOT the nanostat sequencing
     total, and NOT just the target-genus sum. Reads and bases are
     normalised separately since they can diverge (genome size / depth
     completeness differ by organism).
  7. Adds one synthetic genus row per sample so the picture is complete:
       - "false_positive": reads/bases that mapped, but to a genus
         outside TARGET_GENERA (that sample's all-taxa total minus what
         target genera accounted for), normalised against the same
         reads_all_taxa/bases_all_taxa denominator as every other genus.
     For a given sample, relative_abundance_reads across every target
     genus + false_positive sums to 100% (0% if that sample had zero
     all-taxa reads/bases). Same for relative_abundance_bases. This uses
     the same depth-weighted bases as every other genus row.
  8. Join everything to the sample metadata on (basename, Timepoint).

Species-level output carries no relative-abundance/false-positive
figures - this feature is genus-level only.

Files are processed in parallel, one worker per file.

Usage
-----
    python script.py --tsv-dir DIR --metadata meta.tsv \
        --nanostat nanostat.tsv --outdir results

    # more/less parallelism (default: one worker per CPU, capped at file count)
    python script.py -d DIR -m meta.tsv -n nanostat.tsv -o results --workers 16

    # keep every surviving row, not just the summaries
    python script.py -d DIR -m meta.tsv -n nanostat.tsv -o results --keep-rows

    # if the taxonomy column is ever renamed away from 'domain'
    python script.py -d DIR -m meta.tsv -n nanostat.tsv -o results --taxonomy-col gtdb_taxonomy

Author: generated for background mock-community composition analysis.
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd


TARGET_SPECIES = [
    "Listeria monocytogenes",
    "Listeria innocua",
    "Listeria ivanovii",
    "Listeria welshimeri",
    "Bacillus cereus",
    "Bacillus subtilis",
    "Escherichia coli",
    "Pseudomonas aeruginosa",
    "Proteus mirabilis",
    "Rhodococcus equi",
    "Staphylococcus aureus",
    "Yersinia enterocolitica",
    "Yersinia frederiksenii",
    "Enterococcus faecalis",
]

TARGET_GENERA = [
    "Listeria",
    "Bacillus",
    "Escherichia",
    "Pseudomonas",
    "Proteus",
    "Rhodococcus",
    "Staphylococcus",
    "Yersinia",
    "Enterococcus",
]


# TARGET_GENERA, so relative_abundance_reads/_bases sums to 100% per sample.
FALSE_POSITIVE_LABEL = "false_positive"


REQUIRED_BASE_COLS = ["numreads", "covbases", "meandepth", "startpos", "endpos"]


def _rank_value(tax: pd.Series, prefix: str) -> pd.Series:
    """Value following `prefix` (e.g. ``g__Listeria`` -> ``Listeria``).

    Found by searching for the prefix anywhere in the semicolon-joined
    lineage string, NOT by position: these lineage strings start with a
    bare ``Bacteria`` (no leading ``d__`` tag), so "last/first token"
    slicing isn't reliable once ranks are missing or reordered.
    """
    pattern = rf"(?:^|;)\s*{re.escape(prefix)}([^;]*)"
    return tax.astype("string").str.extract(pattern, expand=False).str.strip()


def _rank_tag(tax: pd.Series, prefix: str) -> pd.Series:
    """Whole tagged segment INCLUDING the prefix (e.g. ``s__Listeria monocytogenes_B``),
    kept for the raw/canonical audit trail so the mapping is easy to eyeball."""
    pattern = rf"(?:^|;)\s*({re.escape(prefix)}[^;]*)"
    return tax.astype("string").str.extract(pattern, expand=False).str.strip()


def canonical_species(tax: pd.Series, keep_genus_suffix: bool = False) -> pd.Series:
    """``...;g__Listeria;s__Listeria monocytogenes_B`` -> ``Listeria monocytogenes``.

    GTDB appends ``_A``/``_B``/... to genus and/or species epithet to mark
    polyphyletic splits. Stripping them counts e.g. all of
    ``Bacillus cereus``, ``Bacillus cereus_A``, ``Bacillus cereus_B`` as the
    single target species ``Bacillus cereus``. Vectorised, so it stays cheap
    even when a lot of rows survive the zero-read filter.
    """
    sp = _rank_value(tax, "s__")
    if keep_genus_suffix:
        # only the trailing epithet suffix: Listeria_A booriae keeps its _A
        pattern = r"_[A-Z]+$"
    else:
        # every token: Listeria_A booriae -> Listeria booriae
        pattern = r"_[A-Z]+(?=\s|$)"
    return sp.str.replace(pattern, "", regex=True).fillna("")


def canonical_genus(tax: pd.Series, keep_genus_suffix: bool = False) -> pd.Series:
    """``...;g__Listeria_A;s__...`` -> ``Listeria`` (independent of species).

    Same suffix-collapsing behaviour as canonical_species, controlled by the
    same --keep-genus-suffix flag, so the two stay consistent with each
    other.
    """
    g = _rank_value(tax, "g__")
    if keep_genus_suffix:
        return g.fillna("")
    return g.str.replace(r"_[A-Z]+$", "", regex=True).fillna("")


def parse_sample_id(filename: str, known_timepoints: list[str] | None = None):
    """Pull (basename, timepoint) out of e.g. ``barcode01_N_filt_0h_cov_lis_annot_st.tsv``
    or ``barcode01_N_filt_0h_cov_lis_annot.tsv`` (the trailing ``_st`` is optional,
    in case the MLST-annotation step is no longer run upstream now that ST
    isn't used here).

    Primary rule: basename is the leading token, timepoint is everything after
    ``_filt_``. Falls back to matching against the timepoints present in the
    metadata (longest first, so ``24h_2`` wins over ``24h``).
    """
    stem = filename
    if stem.endswith(".gz"):
        stem = stem[:-3]
    stem = re.sub(r"\.tsv$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"_cov_lis_annot(?:_st)?$", "", stem, flags=re.IGNORECASE)

    basename = stem.split("_")[0]

    timepoint = None
    m = re.search(r"_filt_(?P<tp>.+)$", stem)
    if m:
        timepoint = m.group("tp")

    if known_timepoints and timepoint not in known_timepoints:
        for cand in sorted(known_timepoints, key=len, reverse=True):
            if re.search(rf"(?:^|_){re.escape(str(cand))}(?:_|$)", stem):
                timepoint = cand
                break
    return basename, timepoint


def natural_key(tp) -> tuple:
    """Sort 0h < 4h < 12h < 24h < 24h_2 instead of lexicographically."""
    s = str(tp)
    m = re.match(r"^(\d+(?:\.\d+)?)", s)
    return (float(m.group(1)) if m else float("inf"), s)



AWK = shutil.which("awk") or shutil.which("gawk") or shutil.which("mawk")


def _read_nonzero(path: Path, use_awk: bool, taxonomy_col: str) -> tuple[pd.DataFrame, int]:
    """Read a TSV, keeping only rows with numreads != 0.

    When awk is available and the file is uncompressed, the zero rows are
    discarded before Python ever sees them. On the test data that is ~3x
    faster than parsing the whole file with pandas, because >99% of the bytes
    never get converted into Python/NumPy objects.

    Returns (dataframe_of_kept_rows, n_rows_in_file).
    """
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")

    needed = REQUIRED_BASE_COLS + [taxonomy_col]
    missing = [c for c in needed if c not in header]
    if missing:
        raise ValueError(f"{path.name}: missing required column(s): {', '.join(missing)}")

    if use_awk and AWK and path.suffix != ".gz":
        col = header.index("numreads") + 1                 # awk is 1-based
        prog = 'NR==1 || ($c+0)!=0 { print } END { print NR-1 > "/dev/stderr" }'
        with open(path, "rb") as fh:
            proc = subprocess.run(
                [AWK, "-F", "\t", "-v", f"c={col}", prog],
                stdin=fh, capture_output=True, text=True, check=True,
            )
        # awk counted the rows during the same pass - no second read of 35 MB
        n_total = int(proc.stderr.strip().splitlines()[-1])
        df = pd.read_csv(io.StringIO(proc.stdout), sep="\t")
    else:
        df = pd.read_csv(path, sep="\t")
        n_total = len(df)

    for col in ("numreads", "covbases", "meandepth", "startpos", "endpos"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df = df.loc[df["numreads"] != 0].copy()
    return df, n_total


def process_file(args) -> dict:
    """Worker entry point. Must be top-level so it is picklable."""
    path_str, known_timepoints, keep_genus_suffix, use_awk, keep_rows, taxonomy_col = args
    path = Path(path_str)
    t0 = time.perf_counter()

    try:
        df, n_total = _read_nonzero(path, use_awk, taxonomy_col)
    except Exception as exc:                                  # noqa: BLE001
        return {"file": path.name, "error": f"{type(exc).__name__}: {exc}"}

    basename, timepoint = parse_sample_id(path.name, known_timepoints)

    df["species"] = canonical_species(df[taxonomy_col], keep_genus_suffix)
    df["genus"] = canonical_genus(df[taxonomy_col], keep_genus_suffix)
    df["basename"] = basename
    df["Timepoint"] = timepoint
    df["source_file"] = path.name

    # Depth-weighted mapped bases per row: meandepth * region length. NOT
    # covbases (breadth, capped at reference length) - see module docstring
    # for why that badly undercounts deep targets on ONT-length reads.
    # clip(lower=0) guards against a malformed endpos < startpos row quietly
    # turning into a negative contribution to a sample's total.
    region_len = (df["endpos"] - df["startpos"] + 1).clip(lower=0)
    df["mapped_bases"] = df["meandepth"] * region_len

    # ---- species totals (targets only; raw counts, no % - see docstring) -- #
    sp = (
        df.loc[df["species"].isin(TARGET_SPECIES)]
        .groupby("species", as_index=False)
        .agg(total_reads=("numreads", "sum"),
             total_covbases=("covbases", "sum"),
             n_contigs=("numreads", "size"))
    )

    # ---- genus totals (targets only) --------------------------------------- #
    # Independent of the species filter above: this also catches reads
    # assigned to an off-target species within a target genus, or a row only
    # resolved to genus level, so it is NOT simply the species totals rolled
    # up - it's a real genus-level profile.
    gn = (
        df.loc[df["genus"].isin(TARGET_GENERA)]
        .groupby("genus", as_index=False)
        .agg(total_reads=("numreads", "sum"),
             total_covbases=("covbases", "sum"),
             total_mapped_bases=("mapped_bases", "sum"),
             n_contigs=("numreads", "size"))
    )

    # ---- audit: which raw GTDB labels folded into which canonical name --- #
    df["raw_species_label"] = _rank_tag(df[taxonomy_col], "s__").fillna("<no species tag>")
    df["raw_genus_label"] = _rank_tag(df[taxonomy_col], "g__").fillna("<no genus tag>")

    species_audit = (
        df.groupby(["raw_species_label", "species"], as_index=False, dropna=False)
        .agg(n_contigs=("numreads", "size"), total_reads=("numreads", "sum"))
    )
    genus_audit = (
        df.groupby(["raw_genus_label", "genus"], as_index=False, dropna=False)
        .agg(n_contigs=("numreads", "size"), total_reads=("numreads", "sum"))
    )

    qc = {
        "file": path.name,
        "basename": basename,
        "Timepoint": timepoint,
        "rows_in_file": n_total,
        "rows_kept": len(df),
        "rows_dropped_zero_reads": n_total - len(df),
        "reads_all_taxa": int(df["numreads"].sum()),
        "reads_target_species": int(sp["total_reads"].sum()) if len(sp) else 0,
        "reads_target_genus": int(gn["total_reads"].sum()) if len(gn) else 0,
        "bases_all_taxa": float(df["mapped_bases"].sum()),
        "bases_target_genus": float(gn["total_mapped_bases"].sum()) if len(gn) else 0.0,
        # Raw covbases (breadth) sums, alongside the depth-weighted bases
        # above - used specifically for the breadth-based false-positive
        # metric in false_positive_summary.tsv (see build_false_positive_
        # summary()): a handful of reads spread across a large off-target
        # genome usually rounds meandepth to 0.00 in the coverage TSV,
        # zeroing the depth-weighted figure, but covbases stays non-zero.
        "covbases_all_taxa": int(df["covbases"].sum()),
        "covbases_target_genus": int(gn["total_covbases"].sum()) if len(gn) else 0,
        "seconds": round(time.perf_counter() - t0, 3),
    }

    for frame in (sp, gn, species_audit, genus_audit):
        frame.insert(0, "Timepoint", timepoint)
        frame.insert(0, "basename", basename)

    return {
        "file": path.name,
        "error": None,
        "qc": qc,
        "species": sp,
        "genus": gn,
        "species_audit": species_audit,
        "genus_audit": genus_audit,
        "rows": df if keep_rows else None,
    }



def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str)
    if meta.shape[1] < 2:                       # not really tab separated
        meta = pd.read_csv(path, sep=r"\t|\s{2,}", engine="python", dtype=str)
    meta.columns = [c.strip() for c in meta.columns]
    for col in ("basename", "Timepoint"):
        if col not in meta.columns:
            raise SystemExit(
                f"metadata is missing a '{col}' column; found: {list(meta.columns)}"
            )
    for c in meta.columns:
        meta[c] = meta[c].str.strip()
    dupes = meta.duplicated(["basename", "Timepoint"], keep=False)
    if dupes.any():
        print("WARNING: duplicate (basename, Timepoint) rows in metadata:",
              file=sys.stderr)
        print(meta.loc[dupes].to_string(index=False), file=sys.stderr)
    return meta


def load_nanostat(path: Path) -> pd.DataFrame:
    """Per-sample sequencing totals - the reads/bases that went into the
    alignment in the first place - used as the relative-abundance
    denominator.

    Requires basename, Timepoint, num_of_reads and num_of_bases columns,
    joined onto the coverage-derived results the same way metadata is - on
    (basename, Timepoint). Any other columns (e.g. treatment) are read and
    validated for duplicates but not currently propagated into the output
    tables.
    """
    ns = pd.read_csv(path, sep="\t", dtype=str)
    if ns.shape[1] < 2:                       # not really tab separated
        ns = pd.read_csv(path, sep=r"\t|\s{2,}", engine="python", dtype=str)
    ns.columns = [c.strip() for c in ns.columns]
    required = ["basename", "Timepoint", "num_of_reads", "num_of_bases"]
    missing = [c for c in required if c not in ns.columns]
    if missing:
        raise SystemExit(
            f"nanostat file is missing required column(s): {', '.join(missing)}; "
            f"found: {list(ns.columns)}"
        )
    for c in ns.columns:
        ns[c] = ns[c].str.strip()

    dupes = ns.duplicated(["basename", "Timepoint"], keep=False)
    if dupes.any():
        print("WARNING: duplicate (basename, Timepoint) rows in nanostat file:",
              file=sys.stderr)
        print(ns.loc[dupes].to_string(index=False), file=sys.stderr)

    out = ns[["basename", "Timepoint", "num_of_reads", "num_of_bases"]].copy()
    for c in ("num_of_reads", "num_of_bases"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    bad = out[["num_of_reads", "num_of_bases"]].isna().any(axis=1)
    if bad.any():
        print(f"WARNING: {int(bad.sum())} nanostat row(s) have a non-numeric "
              f"num_of_reads/num_of_bases and will read as NaN in the "
              f"qc_per_file.tsv reference columns (relative abundance is "
              f"unaffected, since it no longer depends on nanostat):",
              file=sys.stderr)
        print(out.loc[bad, ["basename", "Timepoint"]].to_string(index=False),
              file=sys.stderr)
    return out


def attach_metadata(df: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    out = df.merge(meta, on=["basename", "Timepoint"], how="left", validate="m:1")
    meta_cols = [c for c in meta.columns if c not in ("basename", "Timepoint")]
    front = ["basename", "Timepoint"] + meta_cols
    return out[front + [c for c in out.columns if c not in front]]




def complete_grid(long: pd.DataFrame, key: str, targets: list[str],
                  samples: pd.DataFrame) -> pd.DataFrame:
    """Emit an explicit zero row for every sample x target that had no hits.

    Without this, downstream stats silently treat 'absent' as 'missing'.
    """
    grid = samples.merge(pd.DataFrame({key: targets}), how="cross")
    out = grid.merge(long, on=["basename", "Timepoint", key], how="left")
    for c in ("total_reads", "total_covbases", "n_contigs"):
        if c in out.columns:
            out[c] = out[c].fillna(0).astype("int64")
    # mapped bases = meandepth * length; keep exact (unrounded) float values
    if "total_mapped_bases" in out.columns:
        out["total_mapped_bases"] = out["total_mapped_bases"].fillna(0).astype("float64")
    return out


def add_relative_abundance(genus_df: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    """Add relative_abundance_reads / relative_abundance_bases (%) to every
    TARGET_GENERA row, normalised against reads_all_taxa/bases_all_taxa -
    the total reads/depth-weighted mapped bases that mapped to ANY genus
    (target + non-target) for that sample, straight from that sample's own
    coverage TSV - NOT against the nanostat sequencing total, and NOT
    against just the target-genus sum. 0% for a sample with zero all-taxa
    reads/bases, rather than a divide-by-zero.

    `qc` must already have reads_all_taxa/bases_all_taxa (see process_file()).
    """
    out = genus_df.merge(qc[["basename", "Timepoint", "reads_all_taxa", "bases_all_taxa"]],
                         on=["basename", "Timepoint"], how="left")
    for total_col, metric_col, pct_col in (
        ("reads_all_taxa", "total_reads", "relative_abundance_reads"),
        ("bases_all_taxa", "total_mapped_bases", "relative_abundance_bases"),
    ):
        denom = out[total_col].where(out[total_col] > 0)
        out[pct_col] = ((out[metric_col] / denom).fillna(0) * 100).round(4)
    return out.drop(columns=["reads_all_taxa", "bases_all_taxa"])


def append_offtarget_rows(genus_df: pd.DataFrame, qc: pd.DataFrame) -> pd.DataFrame:
    """Add one synthetic genus row per sample - FALSE_POSITIVE_LABEL - so
    relative_abundance_reads (and _bases) across every genus row for a
    sample sums to 100% (0% for a sample with zero all-taxa reads/bases).
    Reads/bases come straight from qc (see main(), where
    false_positive_reads/bases are derived as reads_all_taxa -
    reads_target_genus / bases_all_taxa - bases_target_genus) - this
    function just reshapes those columns into a long-format genus row,
    normalised against the same reads_all_taxa/bases_all_taxa denominator
    as every other genus (see add_relative_abundance).
    """
    r = qc[["basename", "Timepoint", "false_positive_reads", "false_positive_bases",
            "reads_all_taxa", "bases_all_taxa"]].rename(
        columns={"false_positive_reads": "total_reads",
                 "false_positive_bases": "total_mapped_bases"}).copy()
    r["total_reads"] = r["total_reads"].fillna(0).round().astype("int64")
    r["total_mapped_bases"] = r["total_mapped_bases"].fillna(0).astype("float64")
    r["genus"] = FALSE_POSITIVE_LABEL
    # total_covbases (breadth) doesn't have a meaningful summed value across
    # mixed/absent taxa the way depth-weighted bases does - see module
    # docstring on why covbases was dropped as the bases metric.
    r["total_covbases"] = 0
    r["n_contigs"] = 0
    denom_reads = r["reads_all_taxa"].where(r["reads_all_taxa"] > 0)
    denom_bases = r["bases_all_taxa"].where(r["bases_all_taxa"] > 0)
    r["relative_abundance_reads"] = ((r["total_reads"] / denom_reads).fillna(0) * 100).round(4)
    r["relative_abundance_bases"] = ((r["total_mapped_bases"] / denom_bases).fillna(0) * 100).round(4)
    r = r.drop(columns=["reads_all_taxa", "bases_all_taxa"])

    cols = [c for c in genus_df.columns if c in r.columns]
    return pd.concat([genus_df, r[cols]], ignore_index=True)


def build_false_positive_summary(qc: pd.DataFrame) -> pd.DataFrame:
    """Per-sample false-positive detail (reads/bases mapped to a genus
    outside TARGET_GENERA), reported with BOTH bases methods side by side
    so they can be compared directly:

      - *_depth columns: the pipeline's normal depth-weighted bases
        (meandepth * length). Matches the "false_positive" row in
        genus_summary_long.tsv exactly - left as-is. A handful of reads
        spread thinly across a large off-target genome usually rounds
        meandepth to 0.00 in the coverage TSV (2 dp), so this figure can
        legitimately read as 0/near-0 even when reads clearly mapped.
      - *_covbases columns: the same false-positive reads/bases but using
        raw covbases (breadth: positions covered >=1x) instead of the
        depth-weighted figure. This doesn't collapse to zero at very low
        depth, since even a single read covers real bases - at the cost of
        the breadth-capping behaviour the depth-weighted method exists to
        avoid for well-covered TARGET genera (not usually a concern for
        off-target hits, which rarely saturate a genome).

    Each relative-abundance column is normalised against the all-genus
    total computed the SAME way (reads against reads_all_taxa, *_depth
    against bases_all_taxa, *_covbases against covbases_all_taxa), so every
    percentage here is apples-to-apples within its own column. 0% for a
    sample with a zero denominator, rather than a divide-by-zero.

    `qc` must already have false_positive_reads/false_positive_bases/
    false_positive_bases_covbases (see main()).
    """
    out = qc[["basename", "Timepoint",
              "reads_all_taxa", "reads_target_genus", "false_positive_reads",
              "bases_all_taxa", "bases_target_genus", "false_positive_bases",
              "covbases_all_taxa", "covbases_target_genus",
              "false_positive_bases_covbases"]].copy()
    out = out.rename(columns={"false_positive_bases": "false_positive_bases_depth"})

    denom_reads = out["reads_all_taxa"].where(out["reads_all_taxa"] > 0)
    denom_bases_depth = out["bases_all_taxa"].where(out["bases_all_taxa"] > 0)
    denom_bases_cov = out["covbases_all_taxa"].where(out["covbases_all_taxa"] > 0)

    out["relative_abundance_reads"] = (
        (out["false_positive_reads"] / denom_reads).fillna(0) * 100
    ).round(4)
    out["relative_abundance_bases_depth"] = (
        (out["false_positive_bases_depth"] / denom_bases_depth).fillna(0) * 100
    ).round(4)
    out["relative_abundance_bases_covbases"] = (
        (out["false_positive_bases_covbases"] / denom_bases_cov).fillna(0) * 100
    ).round(4)

    return out


def to_wide(long: pd.DataFrame, key: str, value: str,
            id_cols: list[str]) -> pd.DataFrame:
    wide = long.pivot_table(index=id_cols, columns=key, values=value,
                            aggfunc="sum", fill_value=0).reset_index()
    wide.columns.name = None
    return wide



def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("-d", "--tsv-dir", required=True, type=Path,
                    help="directory containing the annotated coverage TSVs")
    ap.add_argument("-m", "--metadata", required=True, type=Path,
                    help="tab-separated metadata with basename + Timepoint columns")
    ap.add_argument("-n", "--nanostat", required=True, type=Path,
                    help="tab-separated per-sample sequencing totals with "
                         "basename, Timepoint, num_of_reads and num_of_bases "
                         "columns - carried into qc_per_file.tsv for "
                         "reference only; genus-level relative abundance is "
                         "normalised against that sample's own all-genus "
                         "mapped totals instead (see docstring)")
    ap.add_argument("-o", "--outdir", default=Path("mock_community_summary"), type=Path,
                    help="output directory (default: mock_community_summary)")
    ap.add_argument("-p", "--pattern", default="*.tsv*",
                    help="glob for input files (default: %(default)s)")
    ap.add_argument("--taxonomy-col", default="taxonomy",
                    help="column holding the full semicolon-joined GTDB lineage "
                         "string (default: %(default)s). In this pipeline's TSVs "
                         "the column literally named 'domain' holds the whole "
                         "d;p;c;o;f;g;s lineage, not just the domain rank - "
                         "override this if that ever changes")
    ap.add_argument("-w", "--workers", type=int, default=0,
                    help="parallel workers; 0 = one per CPU (default)")
    ap.add_argument("--executor", choices=("process", "thread"), default="process",
                    help="process pool (default) or thread pool for small runs")
    ap.add_argument("--keep-rows", action="store_true",
                    help="also write every surviving non-zero row to disk")
    ap.add_argument("--keep-genus-suffix", action="store_true",
                    help="do NOT collapse GTDB genus-split suffixes (e.g. "
                         "'Listeria_A' stays distinct from 'Listeria') - this "
                         "affects both the genus-level rollup and the genus "
                         "token inside species names; the species epithet "
                         "suffix (e.g. 'cereus_A') is always collapsed")
    ap.add_argument("--no-awk", action="store_true",
                    help="disable the awk pre-filter and parse everything in pandas")
    args = ap.parse_args()

    files = sorted(
        f for f in args.tsv_dir.glob(args.pattern)
        if f.is_file() and not f.name.startswith(".")
    )
    if not files:
        raise SystemExit(f"no files matching {args.pattern!r} in {args.tsv_dir}")

    meta = load_metadata(args.metadata)
    nanostat = load_nanostat(args.nanostat)
    known_tps = sorted(meta["Timepoint"].dropna().unique().tolist())

    n_workers = args.workers or min(len(files), os.cpu_count() or 1)
    n_workers = max(1, min(n_workers, len(files)))
    use_awk = (not args.no_awk) and AWK is not None

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"{len(files)} file(s) | {n_workers} {args.executor} worker(s) | "
          f"awk pre-filter: {'on' if use_awk else 'off'} | "
          f"taxonomy col: {args.taxonomy_col!r} | "
          f"nanostat: {args.nanostat}", file=sys.stderr)

    payloads = [
        (str(f), known_tps, args.keep_genus_suffix, use_awk, args.keep_rows, args.taxonomy_col)
        for f in files
    ]

    Pool = ProcessPoolExecutor if args.executor == "process" else ThreadPoolExecutor
    results, errors = [], []
    t0 = time.perf_counter()
    with Pool(max_workers=n_workers) as pool:
        futures = {pool.submit(process_file, p): p[0] for p in payloads}
        for i, fut in enumerate(as_completed(futures), 1):
            res = fut.result()
            if res.get("error"):
                errors.append((res["file"], res["error"]))
                print(f"  [{i}/{len(files)}] FAILED {res['file']}: {res['error']}",
                      file=sys.stderr)
            else:
                results.append(res)
                q = res["qc"]
                print(f"  [{i}/{len(files)}] {res['file']}: "
                      f"{q['rows_kept']:,}/{q['rows_in_file']:,} rows kept "
                      f"({q['seconds']}s)", file=sys.stderr)
    elapsed = time.perf_counter() - t0

    if not results:
        raise SystemExit("every file failed; nothing to write")

    qc = pd.DataFrame([r["qc"] for r in results])
    qc = qc.sort_values(["basename", "Timepoint"], key=lambda s:
                        s.map(natural_key) if s.name == "Timepoint" else s)

    qc = qc.merge(nanostat, on=["basename", "Timepoint"], how="left")
    missing_nanostat = qc["num_of_reads"].isna()
    if missing_nanostat.any():
        print(f"WARNING: {int(missing_nanostat.sum())} sample(s) had no nanostat "
              f"match (joined on basename + Timepoint) - num_of_reads/num_of_bases "
              f"will read as blank for these in qc_per_file.tsv (relative abundance "
              f"is unaffected, since it no longer depends on nanostat):", file=sys.stderr)
        for _, row in qc.loc[missing_nanostat, ["basename", "Timepoint"]].iterrows():
            print(f"  basename={row['basename']!r} Timepoint={row['Timepoint']!r}",
                  file=sys.stderr)

    # false_positive = mapped, but to a genus outside TARGET_GENERA. Always
    # >= 0 by construction (reads_target_genus/bases_target_genus are a
    # subset sum of reads_all_taxa/bases_all_taxa), so no negative-value
    # check is needed here.
    qc["false_positive_reads"] = qc["reads_all_taxa"] - qc["reads_target_genus"]
    qc["false_positive_bases"] = qc["bases_all_taxa"] - qc["bases_target_genus"]
    # Breadth-based (covbases) version of the same false-positive bases
    # figure - see build_false_positive_summary() / false_positive_summary.tsv.
    qc["false_positive_bases_covbases"] = qc["covbases_all_taxa"] - qc["covbases_target_genus"]

    samples = (
        qc[["basename", "Timepoint"]].drop_duplicates()
        .sort_values("Timepoint", key=lambda s: s.map(natural_key), kind="stable")
        .reset_index(drop=True)
    )

    species_long = pd.concat([r["species"] for r in results], ignore_index=True)
    genus_long = pd.concat([r["genus"] for r in results], ignore_index=True)

    species_audit = (
        pd.concat([r["species_audit"] for r in results], ignore_index=True)
        .groupby(["raw_species_label", "species"], as_index=False)
        .agg(n_contigs=("n_contigs", "sum"), total_reads=("total_reads", "sum"))
        .sort_values(["species", "raw_species_label"])
    )
    genus_audit = (
        pd.concat([r["genus_audit"] for r in results], ignore_index=True)
        .groupby(["raw_genus_label", "genus"], as_index=False)
        .agg(n_contigs=("n_contigs", "sum"), total_reads=("total_reads", "sum"))
        .sort_values(["genus", "raw_genus_label"])
    )

    species_long = complete_grid(species_long, "species", TARGET_SPECIES, samples)
    genus_long = complete_grid(genus_long, "genus", TARGET_GENERA, samples)

    # relative abundance against that sample's own all-genus mapped totals,
    # plus a false_positive row so it reconciles to 100% - genus only.
    genus_long = add_relative_abundance(genus_long, qc)
    genus_long = append_offtarget_rows(genus_long, qc)

    recon = (genus_long.groupby(["basename", "Timepoint"], as_index=False)
             .agg(sum_pct_reads=("relative_abundance_reads", "sum"),
                  sum_pct_bases=("relative_abundance_bases", "sum")))
    bad_recon = recon[~(recon["sum_pct_reads"].round(2).isin([0.0, 100.0]) &
                        recon["sum_pct_bases"].round(2).isin([0.0, 100.0]))]
    if len(bad_recon):
        print(f"WARNING: {len(bad_recon)} sample(s) did NOT reconcile to 100% "
              f"across target genera + false_positive - check for data issues:",
              file=sys.stderr)
        print(bad_recon.to_string(index=False), file=sys.stderr)

    false_positive_summary = build_false_positive_summary(qc)

    species_long = attach_metadata(species_long, meta)
    genus_long = attach_metadata(genus_long, meta)
    qc = attach_metadata(qc, meta)
    false_positive_summary = attach_metadata(false_positive_summary, meta)

    meta_keys = set(zip(meta["basename"], meta["Timepoint"]))
    unmatched = [(b, t) for b, t in zip(samples["basename"], samples["Timepoint"])
                 if (b, t) not in meta_keys]
    if unmatched:
        print(f"WARNING: {len(unmatched)} sample(s) had no metadata match "
              f"(joined on basename + Timepoint):", file=sys.stderr)
        for b, t in unmatched:
            print(f"  basename={b!r} Timepoint={t!r}", file=sys.stderr)

    metric_cols = {"total_reads", "total_covbases", "total_mapped_bases", "n_contigs",
                   "relative_abundance_reads", "relative_abundance_bases"}
    id_cols_species = [c for c in species_long.columns if c not in metric_cols | {"species"}]
    id_cols_genus = [c for c in genus_long.columns if c not in metric_cols | {"genus"}]

    written = []

    def write(df: pd.DataFrame, name: str):
        p = args.outdir / name
        df.to_csv(p, sep="\t", index=False)
        written.append(p)

    
    write(species_long.sort_values(["Timepoint", "basename", "species"],
                                   key=lambda s: s.map(natural_key)
                                   if s.name == "Timepoint" else s),
          "species_summary_long.tsv")
    write(to_wide(species_long, "species", "total_reads", id_cols_species),
          "species_reads_wide.tsv")
    write(to_wide(species_long, "species", "total_covbases", id_cols_species),
          "species_covbases_wide.tsv")

    
    # genus_long includes a "false_positive" row alongside TARGET_GENERA
    # (see append_offtarget_rows), so every wide table below naturally gets
    # it as an extra column - a built-in reconciliation check: the
    # relative-abundance wide tables' rows should sum to ~100 each.
    write(genus_long.sort_values(["Timepoint", "basename", "genus"],
                                 key=lambda s: s.map(natural_key)
                                 if s.name == "Timepoint" else s),
          "genus_summary_long.tsv")
    write(to_wide(genus_long, "genus", "total_reads", id_cols_genus),
          "genus_reads_wide.tsv")
    write(to_wide(genus_long, "genus", "total_covbases", id_cols_genus),
          "genus_covbases_wide.tsv")
    write(to_wide(genus_long, "genus", "total_mapped_bases", id_cols_genus),
          "genus_mapped_bases_wide.tsv")
    write(to_wide(genus_long, "genus", "relative_abundance_reads", id_cols_genus),
          "genus_relative_abundance_reads_wide.tsv")
    write(to_wide(genus_long, "genus", "relative_abundance_bases", id_cols_genus),
          "genus_relative_abundance_bases_wide.tsv")

    write(qc, "qc_per_file.tsv")
    write(false_positive_summary.sort_values(["Timepoint", "basename"],
                                             key=lambda s: s.map(natural_key)
                                             if s.name == "Timepoint" else s),
          "false_positive_summary.tsv")
    write(species_audit, "species_label_audit.tsv")
    write(genus_audit, "genus_label_audit.tsv")

    def safe_pct(total: pd.Series) -> pd.Series:
        s = total.sum()
        if not s:
            return total.astype(float) * 0.0
        return (total / s * 100).round(4)

    grand_sp = (species_long.groupby("species", as_index=False)
                .agg(total_reads=("total_reads", "sum"),
                     total_covbases=("total_covbases", "sum")))
    grand_sp["pct_reads"] = safe_pct(grand_sp["total_reads"])
    grand_sp["pct_covbases"] = safe_pct(grand_sp["total_covbases"])

    # Denominator is each sample's own all-genus mapped totals (reads_all_taxa/
    # bases_all_taxa), computed directly from the coverage TSVs - unlike the
    # old nanostat-based denominator, every sample has this by construction,
    # so no scoping/filtering step is needed here.
    grand_gen = (genus_long.groupby("genus", as_index=False)
                 .agg(total_reads=("total_reads", "sum"),
                      total_covbases=("total_covbases", "sum"),
                      total_mapped_bases=("total_mapped_bases", "sum")))
    grand_reads_all_taxa = qc["reads_all_taxa"].sum()
    grand_bases_all_taxa = qc["bases_all_taxa"].sum()
    grand_gen["relative_abundance_reads"] = (
        (grand_gen["total_reads"] / grand_reads_all_taxa * 100).round(4)
        if grand_reads_all_taxa else 0.0
    )
    grand_gen["relative_abundance_bases"] = (
        (grand_gen["total_mapped_bases"] / grand_bases_all_taxa * 100).round(4)
        if grand_bases_all_taxa else 0.0
    )

    write(grand_sp, "grand_totals_species.tsv")
    write(grand_gen, "grand_totals_genus.tsv")

    if args.keep_rows:
        rows = pd.concat([r["rows"] for r in results if r["rows"] is not None],
                         ignore_index=True)
        rows = attach_metadata(rows, meta)
        p = args.outdir / "nonzero_rows.tsv.gz"
        rows.to_csv(p, sep="\t", index=False, compression="gzip")
        written.append(p)

    print(f"\nProcessed {len(results)} file(s) in {elapsed:.2f}s "
          f"({qc['rows_in_file'].sum():,} rows read, "
          f"{qc['rows_kept'].sum():,} kept, "
          f"{qc['rows_dropped_zero_reads'].sum():,} dropped as zero-read)",
          file=sys.stderr)
    if errors:
        print(f"{len(errors)} file(s) failed - see messages above", file=sys.stderr)
    print("\nWrote:", file=sys.stderr)
    for p in written:
        print(f"  {p}", file=sys.stderr)

    print("\n== Grand total reads / covbases per species (% of target-species total) ==")
    print(grand_sp.to_string(index=False))
    print("\n== Grand total reads / mapped bases per genus, incl. false_positive "
          "(% of total reads/bases mapped to any genus) ==")
    print(grand_gen.to_string(index=False))

    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())