#!/usr/bin/env python3
"""
listeria_cov_summary.py
=======================

Summarise samtools-coverage style TSVs annotated with GTDB taxonomy and MLST
ST, expressing Listeria composition as relative abundance with a
DELIBERATELY MIXED denominator, by row:
  - TARGET rows (the four species, the four STs) are % of the total
    whole-community reads/bases for that sample (reads/bases mapped to ALL
    genera, supplied externally via --genus-totals) - the SAME denominator
    at both levels, so a target ST's % is guaranteed to sit at or below its
    species' Lm % on the same plot.
  - false_positive (species AND ST level) and untyped (ST level only) are
    instead % of the Listeria-only pool they were drawn from - all-Listeria
    genus reads for species-level false_positive, monocytogenes-only reads
    for ST-level false_positive/untyped - since these are typing/
    contamination-quality metrics that are more meaningful scoped to the
    Listeria reads actually recovered, not diluted by the whole community.
Not the target list, not the raw nanostat sequencing total.

NB: ST9 and ST35 are combined and reported as ST35 everywhere in this script
(aggregation, targets, output) - a row typed ST9 has its ST rewritten to 35
immediately after normalisation, before anything else touches it.

For every ``*_cov_lis_annot_st.tsv`` file in a directory it will:

  1. Drop every row with ``numreads == 0``  (typically ~99.9% of rows).
  2. Determine, per row, whether it is a genuine SPECIES-level target:
       - species is innocua, ivanovii, welshimeri OR monocytogenes: target
         on species match alone. ST is NOT consulted for species-level
         classification - a monocytogenes row with an untyped or
         off-target ST still counts as a "Listeria monocytogenes" target,
         same as any other monocytogenes row (ST only matters for the
         SEPARATE ST-level table in step 3).
     Species-level false_positive = every OTHER (non-target) Listeria
     species (i.e. genus-wide Listeria minus the species-level target sum).
     Non-Listeria genera never enter this table at all (they're neither a
     target nor a false_positive here).
  3. Independently, determine ST-level targets: ST in {1, 35, 70, 121}
     (post-merge). Because ST is only ever called on monocytogenes rows in
     this pipeline, ST-level composition is scoped to monocytogenes rows
     only - innocua/ivanovii/welshimeri never appear in the ST table at
     all, since they never receive an ST call. Within that monocytogenes
     subset:
       - target ST      -> the four target-ST rows
       - false_positive  -> monocytogenes rows with a called ST that is NOT
                            in the target set (e.g. ST8, ST6, ...)
       - untyped         -> monocytogenes rows with no ST call at all
  4. For both levels, also sum a DEPTH-WEIGHTED bases figure per target:
     ``sum(meandepth * (endpos - startpos + 1))``, NOT ``covbases`` -
     covbases is breadth of coverage, capped at reference length, which on
     ONT-length reads badly undercounts genuinely deep targets (a handful of
     long reads can blanket a small genome at >=1x almost immediately).
     depth * length reconstructs actual mapped bases.
  5. Look up, per sample (basename + Timepoint), the total reads and bases
     that went into the alignment in the first place, from a nanostat TSV
     (--nanostat; requires basename, Timepoint, num_of_reads, num_of_bases).
     These are carried through to qc_per_file.tsv as context only - they do
     NOT feed any relative-abundance calculation.
  6. Look up, per sample, the whole-community totals - reads/bases mapped
     to ALL genera present, not just Listeria - from a separate TSV
     (--genus-totals; requires basename, Timepoint, reads_all_taxa,
     bases_all_taxa). This is supplied externally rather than summed from
     the coverage TSV itself, since that TSV may only contain reads that
     mapped to a curated multi-genus reference panel, not literally every
     genus in the sample (the coverage TSV's own all-taxa sum is still kept,
     as reads_all_taxa_covfile/bases_all_taxa_covfile in qc_per_file.tsv,
     purely as a sanity check against the externally supplied total).
     Express every species-level AND ST-level TARGET row as relative
     abundance (%) against THIS reads_all_taxa/bases_all_taxa total. The
     false_positive (species and ST level) and untyped (ST level) synthetic
     rows instead use reads_all_listeria/bases_all_listeria (species) or
     reads_all_monocytogenes/bases_all_monocytogenes (ST) - see the mixed-
     denominator note at the top of this docstring. Because target rows
     share one denominator across both levels, and target-ST reads are a
     subset of Lm-species reads are a subset of genus-wide Listeria reads,
     target percentages nest correctly: an individual target ST's % is
     always <= its species' target %. false_positive/untyped percentages
     are NOT part of that nesting relationship and should not be summed
     together with target percentages - each row's % is only meaningful
     against its own stated denominator (see qc_per_file.tsv for the raw
     reads_all_taxa/reads_all_listeria/reads_all_monocytogenes columns).
     0% for a sample with no genus-totals match (see the WARNING for those).
  7. Join everything to the sample metadata on (basename, Timepoint).

Files are processed in parallel, one worker per file.

Usage
-----
    python listeria_cov_summary.py --tsv-dir DIR --metadata meta.tsv \
        --nanostat nanostat.tsv --outdir results

    # more/less parallelism (default: one worker per CPU, capped at file count)
    python listeria_cov_summary.py -d DIR -m meta.tsv -n nanostat.tsv -o results --workers 16

    # keep every surviving row, not just the summaries
    python listeria_cov_summary.py -d DIR -m meta.tsv -n nanostat.tsv -o results --keep-rows

Author: generated for Listeria sponge-swab timecourse analysis.
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

# --------------------------------------------------------------------------- #
# Targets
# --------------------------------------------------------------------------- #

TARGET_SPECIES = [
    "Listeria innocua",
    "Listeria ivanovii",
    "Listeria monocytogenes",
    "Listeria welshimeri",
]

# Species-level classification does NOT gate monocytogenes on ST - every
# monocytogenes row counts toward the species-level "Listeria monocytogenes"
# target regardless of ST (untyped and off-target-ST included). ST_GATED_SPECIES
# is kept only to identify the monocytogenes subset for the SEPARATE ST-level
# table (is_monocytogenes, below), where ST typing still matters.
ST_GATED_SPECIES = "Listeria monocytogenes"

# ST9 and ST35 are combined and reported as ST35 - applied to ST_norm right
# after normalisation, before anything downstream ever sees a "9".
ST_MERGE_MAP = {"9": "35"}

TARGET_STS = ["1", "35", "70", "121"]          # stored bare, displayed as ST1 ... (post-merge, no separate "9")

# Synthetic labels added alongside the target lists so
# relative_abundance_reads/_bases sums to 100% per sample at each level
# independently. Species level gets false_positive only (every row is
# always "typed" to a species). ST level gets false_positive AND untyped,
# since a monocytogenes row can map with no ST call at all.
FALSE_POSITIVE_LABEL = "false_positive"
UNTYPED_LABEL = "untyped"

# meandepth/startpos/endpos are needed to reconstruct depth-weighted mapped
# bases (see module docstring).
NEEDED_COLS = ["numreads", "covbases", "meandepth", "startpos", "endpos", "gtdb_taxonomy", "ST"]

# --------------------------------------------------------------------------- #
# Normalisation helpers
# --------------------------------------------------------------------------- #

_NULL_ST = {"", "-", "na", "n/a", "nan", "none", "null", "<na>", "*", "?"}


def canonical_species(tax: pd.Series, keep_genus_suffix: bool = False) -> pd.Series:
    """``d__Bacteria;...;s__Listeria monocytogenes_B`` -> ``Listeria monocytogenes``.

    GTDB appends ``_A``/``_B``/... to genus and/or species epithet to mark
    polyphyletic splits. Stripping them counts the four target species as
    single biological units, which is what was asked for. Vectorised, so it
    stays cheap even when a lot of rows survive the zero-read filter.
    """
    last = (
        tax.astype("string")
        .str.split(";").str[-1]
        .str.strip()
        .str.replace(r"^[a-z]__", "", regex=True)
    )
    if keep_genus_suffix:
        # only the trailing epithet suffix: Listeria_A booriae keeps its _A
        pattern = r"_[A-Z]+$"
    else:
        # every token: Listeria_A booriae -> Listeria booriae
        pattern = r"_[A-Z]+(?=\s|$)"
    return last.str.replace(pattern, "", regex=True).fillna("")


def canonical_st(value: pd.Series) -> pd.Series:
    """Normalise an ST column to bare number strings, with <NA> for untyped."""
    s = value.astype("string").str.strip()
    s = s.str.replace(r"^ST[_\-\s]*", "", regex=True, flags=re.IGNORECASE)
    s = s.str.replace(r"\.0+$", "", regex=True)   # 121.0 -> 121 if read as float
    return s.mask(s.str.lower().isin(_NULL_ST), pd.NA)


def parse_sample_id(filename: str, known_timepoints: list[str] | None = None):
    """Pull (basename, timepoint) out of e.g. ``barcode01_N_filt_0h_cov_lis_annot_st.tsv``.

    Primary rule: basename is the leading token, timepoint is everything after
    ``_filt_``. Falls back to matching against the timepoints present in the
    metadata (longest first, so ``24h_2`` wins over ``24h``).
    """
    stem = filename
    if stem.endswith(".gz"):
        stem = stem[:-3]
    stem = re.sub(r"\.tsv$", "", stem, flags=re.IGNORECASE)
    stem = re.sub(r"_cov_lis_annot_st$", "", stem, flags=re.IGNORECASE)

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


# --------------------------------------------------------------------------- #
# Per-file worker
# --------------------------------------------------------------------------- #

AWK = shutil.which("awk") or shutil.which("gawk") or shutil.which("mawk")


def _read_nonzero(path: Path, use_awk: bool) -> tuple[pd.DataFrame, int]:
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

    missing = [c for c in NEEDED_COLS if c not in header]
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
        df = pd.read_csv(io.StringIO(proc.stdout), sep="\t", dtype={"ST": "string"})
    else:
        df = pd.read_csv(path, sep="\t", dtype={"ST": "string"})
        n_total = len(df)

    for col in ("numreads", "covbases", "meandepth", "startpos", "endpos"):
        df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    df = df.loc[df["numreads"] != 0].copy()
    return df, n_total


def process_file(args) -> dict:
    """Worker entry point. Must be top-level so it is picklable."""
    path_str, known_timepoints, keep_genus_suffix, use_awk, keep_rows = args
    path = Path(path_str)
    t0 = time.perf_counter()

    try:
        df, n_total = _read_nonzero(path, use_awk)
    except Exception as exc:                                  # noqa: BLE001
        return {"file": path.name, "error": f"{type(exc).__name__}: {exc}"}

    basename, timepoint = parse_sample_id(path.name, known_timepoints)

    df["species"] = canonical_species(df["gtdb_taxonomy"], keep_genus_suffix)
    df["ST_norm"] = canonical_st(df["ST"]).replace(ST_MERGE_MAP)
    df["basename"] = basename
    df["Timepoint"] = timepoint
    df["source_file"] = path.name

    # Depth-weighted mapped bases per row: meandepth * region length. NOT
    # covbases (breadth, capped at reference length) - see module docstring.
    region_len = (df["endpos"] - df["startpos"] + 1).clip(lower=0)
    df["mapped_bases"] = df["meandepth"] * region_len

    # ---- genus-wide Listeria mask ------------------------------------------ #
    # First whitespace token of the canonical species name, with any GTDB
    # suffix stripped, so this is robust to --keep-genus-suffix as well
    # (e.g. "Listeria_A booriae" -> genus token "Listeria"). This is the
    # SPECIES-LEVEL denominator: every Listeria species, target or not.
    genus_token = df["species"].str.split(" ").str[0].str.replace(r"_[A-Z]+$", "", regex=True)
    is_listeria = genus_token == "Listeria"

    # ---- species-level target mask ---------------------------------------- #
    # innocua/ivanovii/welshimeri/monocytogenes: target on species match
    # alone. ST is NOT consulted here - an untyped or off-target-ST
    # monocytogenes row still counts as a "Listeria monocytogenes"
    # species-level target, same as any other monocytogenes row. (ST
    # continues to gate the SEPARATE ST-level table below, via
    # is_monocytogenes/target_st_mask - this change only affects the
    # species-level target/false_positive split.)
    is_target_species = df["species"].isin(TARGET_SPECIES)
    is_nontarget_listeria = is_listeria & ~is_target_species

    # ---- species totals (targets only, per is_target_species above) ------- #
    sp = (
        df.loc[is_target_species]
        .groupby("species", as_index=False)
        .agg(total_reads=("numreads", "sum"),
             total_covbases=("covbases", "sum"),
             total_mapped_bases=("mapped_bases", "sum"),
             n_contigs=("numreads", "size"))
    )

    # ---- ST-level masks (monocytogenes-only universe) ---------------------- #
    # ST is only ever called on monocytogenes rows in this pipeline, so the
    # ST-level denominator is "every monocytogenes row", NOT genus-wide
    # Listeria - innocua/ivanovii/welshimeri never appear in the ST table.
    is_monocytogenes = df["species"] == ST_GATED_SPECIES
    target_st_mask = is_monocytogenes & df["ST_norm"].isin(TARGET_STS)
    nontarget_st_mask = is_monocytogenes & df["ST_norm"].notna() & ~df["ST_norm"].isin(TARGET_STS)
    untyped_st_mask = is_monocytogenes & df["ST_norm"].isna()

    # ---- ST totals (targets only), plus a species breakdown -------------- #
    st = (
        df.loc[target_st_mask]
        .groupby("ST_norm", as_index=False)
        .agg(total_reads=("numreads", "sum"),
             total_covbases=("covbases", "sum"),
             total_mapped_bases=("mapped_bases", "sum"),
             n_contigs=("numreads", "size"))
        .rename(columns={"ST_norm": "ST"})
    )
    st_by_sp = (
        df.loc[target_st_mask]
        .groupby(["ST_norm", "species"], as_index=False)
        .agg(total_reads=("numreads", "sum"),
             total_covbases=("covbases", "sum"),
             n_contigs=("numreads", "size"))
        .rename(columns={"ST_norm": "ST"})
    )

    # ---- audit: which raw GTDB labels folded into which canonical name --- #
    # raw_label must be a real column: pandas 2.x drops an inline Series
    # grouper from the result (FutureWarning), pandas 3.x keeps it.
    df["raw_label"] = (
        df["gtdb_taxonomy"].astype("string")
        .str.split(";").str[-1].str.strip()
        .fillna("<no taxonomy>")
    )
    audit = (
        df.groupby(["raw_label", "species"], as_index=False, dropna=False)
        .agg(n_contigs=("numreads", "size"), total_reads=("numreads", "sum"))
    )

    # ---- ST audit: raw ST value (pre-merge) -> canonical ST (post-merge) - #
    # exists mainly so the ST9 -> ST35 merge is checkable against real data.
    df["raw_st_label"] = df["ST"].astype("string").fillna("<no ST>")
    st_audit = (
        df.groupby(["raw_st_label", "ST_norm"], as_index=False, dropna=False)
        .agg(n_contigs=("numreads", "size"), total_reads=("numreads", "sum"))
        .rename(columns={"ST_norm": "ST"})
    )

    qc = {
        "file": path.name,
        "basename": basename,
        "Timepoint": timepoint,
        "rows_in_file": n_total,
        "rows_kept": len(df),
        "rows_dropped_zero_reads": n_total - len(df),
        # Sum of every row in THIS coverage TSV (every genus it happened to
        # contain). Named "_covfile" to distinguish it from the externally
        # supplied "reads_all_taxa"/"bases_all_taxa" (see --genus-totals in
        # main()), which is the TRUE whole-community total and the actual
        # relative-abundance denominator - the coverage TSV may only reflect
        # reads that mapped to a curated reference panel, not literally every
        # genus in the sample, so the two are not expected to match exactly.
        # Kept here only for the sanity check in main() (in-file sum should
        # never exceed the externally supplied total).
        "reads_all_taxa_covfile": int(df["numreads"].sum()),
        "bases_all_taxa_covfile": int(round(df["mapped_bases"].sum())),
        # genus-wide Listeria and monocytogenes-only totals - informational
        # context only now (see main(): both species- and ST-level relative
        # abundance are expressed against the externally supplied
        # reads_all_taxa/bases_all_taxa instead, so that ST% is guaranteed to
        # nest under species% under the Listeria genus%, all against the same
        # whole-community denominator).
        "reads_all_listeria": int(df.loc[is_listeria, "numreads"].sum()),
        "bases_all_listeria": int(round(df.loc[is_listeria, "mapped_bases"].sum())),
        "reads_target_species": int(sp["total_reads"].sum()) if len(sp) else 0,
        "bases_target_species": int(round(sp["total_mapped_bases"].sum())) if len(sp) else 0,
        "reads_nontarget_species": int(df.loc[is_nontarget_listeria, "numreads"].sum()),
        "bases_nontarget_species": int(round(df.loc[is_nontarget_listeria, "mapped_bases"].sum())),
        "reads_all_monocytogenes": int(df.loc[is_monocytogenes, "numreads"].sum()),
        "bases_all_monocytogenes": int(round(df.loc[is_monocytogenes, "mapped_bases"].sum())),
        "reads_target_st": int(st["total_reads"].sum()) if len(st) else 0,
        "bases_target_st": int(round(st["total_mapped_bases"].sum())) if len(st) else 0,
        "reads_nontarget_st": int(df.loc[nontarget_st_mask, "numreads"].sum()),
        "bases_nontarget_st": int(round(df.loc[nontarget_st_mask, "mapped_bases"].sum())),
        "reads_untyped_st": int(df.loc[untyped_st_mask, "numreads"].sum()),
        "bases_untyped_st": int(round(df.loc[untyped_st_mask, "mapped_bases"].sum())),
        "rows_without_ST": int(df["ST_norm"].isna().sum()),
        "seconds": round(time.perf_counter() - t0, 3),
    }

    for frame in (sp, st, st_by_sp, audit, st_audit):
        frame.insert(0, "Timepoint", timepoint)
        frame.insert(0, "basename", basename)

    return {
        "file": path.name,
        "error": None,
        "qc": qc,
        "species": sp,
        "st": st,
        "st_by_species": st_by_sp,
        "audit": audit,
        "st_audit": st_audit,
        "rows": df if keep_rows else None,
    }


# --------------------------------------------------------------------------- #
# Metadata
# --------------------------------------------------------------------------- #

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
    alignment in the first place. Carried through to qc_per_file.tsv as
    context only; does NOT feed relative abundance (see load_genus_totals
    for the actual denominator source). Requires basename, Timepoint,
    num_of_reads, num_of_bases, joined the same way metadata is - on
    (basename, Timepoint).
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
              f"num_of_reads/num_of_bases and will read as NaN in "
              f"qc_per_file.tsv (context only - does not affect relative "
              f"abundance):", file=sys.stderr)
        print(out.loc[bad, ["basename", "Timepoint"]].to_string(index=False),
              file=sys.stderr)
    return out


def load_genus_totals(path: Path) -> pd.DataFrame:
    """Per-sample whole-community totals - reads/bases mapped to ALL genera
    (not just Listeria) - the actual relative-abundance denominator for BOTH
    species-level and ST-level output, so the two levels share one common
    whole and nest properly (target-ST reads <= Lm-species reads <= all-taxa
    reads, all as % of the same total). Supplied externally rather than
    summed from the coverage TSV itself, since that TSV may only contain
    reads that mapped to a curated multi-genus reference panel, not
    literally every genus present in the sample. Requires basename,
    Timepoint, reads_all_taxa, bases_all_taxa, joined the same way metadata
    and nanostat are - on (basename, Timepoint).
    """
    gt = pd.read_csv(path, sep="\t", dtype=str)
    if gt.shape[1] < 2:                       # not really tab separated
        gt = pd.read_csv(path, sep=r"\t|\s{2,}", engine="python", dtype=str)
    gt.columns = [c.strip() for c in gt.columns]
    required = ["basename", "Timepoint", "reads_all_taxa", "bases_all_taxa"]
    missing = [c for c in required if c not in gt.columns]
    if missing:
        raise SystemExit(
            f"genus-totals file is missing required column(s): {', '.join(missing)}; "
            f"found: {list(gt.columns)}"
        )
    for c in gt.columns:
        gt[c] = gt[c].str.strip()

    dupes = gt.duplicated(["basename", "Timepoint"], keep=False)
    if dupes.any():
        print("WARNING: duplicate (basename, Timepoint) rows in genus-totals file:",
              file=sys.stderr)
        print(gt.loc[dupes].to_string(index=False), file=sys.stderr)

    out = gt[["basename", "Timepoint", "reads_all_taxa", "bases_all_taxa"]].copy()
    for c in ("reads_all_taxa", "bases_all_taxa"):
        out[c] = pd.to_numeric(out[c], errors="coerce")
    bad = out[["reads_all_taxa", "bases_all_taxa"]].isna().any(axis=1)
    if bad.any():
        print(f"WARNING: {int(bad.sum())} genus-totals row(s) have a non-numeric "
              f"reads_all_taxa/bases_all_taxa and will read as NaN (samples "
              f"will show 0% relative abundance):", file=sys.stderr)
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


# --------------------------------------------------------------------------- #
# Grid completion + wide pivots
# --------------------------------------------------------------------------- #

def complete_grid(long: pd.DataFrame, key: str, targets: list[str],
                  samples: pd.DataFrame) -> pd.DataFrame:
    """Emit an explicit zero row for every sample x target that had no hits.

    Without this, downstream stats silently treat 'absent' as 'missing'.
    """
    grid = samples.merge(pd.DataFrame({key: targets}), how="cross")
    out = grid.merge(long, on=["basename", "Timepoint", key], how="left")
    for c in ("total_reads", "total_covbases", "total_mapped_bases", "n_contigs"):
        if c in out.columns:
            out[c] = out[c].fillna(0).astype("int64")
    return out


def to_wide(long: pd.DataFrame, key: str, value: str,
            id_cols: list[str]) -> pd.DataFrame:
    wide = long.pivot_table(index=id_cols, columns=key, values=value,
                            aggfunc="sum", fill_value=0).reset_index()
    wide.columns.name = None
    return wide


def add_relative_abundance(target_df: pd.DataFrame, qc: pd.DataFrame,
                           reads_denom_col: str, bases_denom_col: str) -> pd.DataFrame:
    """Add relative_abundance_reads / relative_abundance_bases (%) to every
    target row (species- or ST-level - generic over total_reads/
    total_mapped_bases, whichever key column the caller used), normalised
    against `reads_denom_col`/`bases_denom_col` from `qc` - the per-sample
    Listeria-relative denominator for that level (genus-wide Listeria totals
    for species, monocytogenes-only totals for ST - see main()), NOT against
    nanostat and NOT against the target-list sum. 0% for a sample with no
    reads in that denominator.
    """
    out = target_df.merge(qc[["basename", "Timepoint", reads_denom_col, bases_denom_col]],
                          on=["basename", "Timepoint"], how="left")
    for total_col, metric_col, pct_col in (
        (reads_denom_col, "total_reads", "relative_abundance_reads"),
        (bases_denom_col, "total_mapped_bases", "relative_abundance_bases"),
    ):
        denom = out[total_col].where(out[total_col] > 0)
        out[pct_col] = ((out[metric_col] / denom).fillna(0) * 100).round(4)
    return out.drop(columns=[reads_denom_col, bases_denom_col])


def append_offtarget_rows(target_df: pd.DataFrame, qc: pd.DataFrame, key_col: str,
                          rows_spec: list[tuple[str, str, str]],
                          reads_denom_col: str, bases_denom_col: str) -> pd.DataFrame:
    """Add one synthetic row per sample per entry in `rows_spec`, under
    `key_col` ("species" or "ST"), so relative_abundance_reads (and _bases)
    across every row for a sample sums to 100% (0% for a sample with no
    reads in that level's denominator).

    `rows_spec` is a list of (label, reads_col, bases_col), pulled straight
    from `qc` (see main()) - species level passes one entry
    (false_positive = non-target Listeria); ST level passes two
    (false_positive = off-target-ST monocytogenes, untyped = no ST call).
    Percentages use the same `reads_denom_col`/`bases_denom_col` as the real
    target rows for that level, so the whole table reconciles to 100%.
    """
    rows = []
    for label, reads_col, bases_col in rows_spec:
        r = qc[["basename", "Timepoint", reads_col, bases_col,
                reads_denom_col, bases_denom_col]].rename(
            columns={reads_col: "total_reads", bases_col: "total_mapped_bases"}).copy()
        r["total_reads"] = r["total_reads"].fillna(0).round().astype("int64")
        r["total_mapped_bases"] = r["total_mapped_bases"].fillna(0).round().astype("int64")
        r[key_col] = label
        # total_covbases (breadth) has no meaningful summed value across
        # mixed/absent taxa - see module docstring on the depth-weighted
        # bases metric used for the percentages instead.
        r["total_covbases"] = 0
        r["n_contigs"] = 0
        denom_reads = r[reads_denom_col].where(r[reads_denom_col] > 0)
        denom_bases = r[bases_denom_col].where(r[bases_denom_col] > 0)
        r["relative_abundance_reads"] = ((r["total_reads"] / denom_reads).fillna(0) * 100).round(4)
        r["relative_abundance_bases"] = ((r["total_mapped_bases"] / denom_bases).fillna(0) * 100).round(4)
        rows.append(r.drop(columns=[reads_denom_col, bases_denom_col]))

    extra = pd.concat(rows, ignore_index=True)
    cols = [c for c in target_df.columns if c in extra.columns]
    return pd.concat([target_df, extra[cols]], ignore_index=True)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

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
                         "columns - carried through to qc_per_file.tsv as "
                         "context only; does not feed relative abundance "
                         "(see --genus-totals for that)")
    ap.add_argument("-g", "--genus-totals", required=True, type=Path,
                    help="tab-separated per-sample whole-community totals "
                         "with basename, Timepoint, reads_all_taxa and "
                         "bases_all_taxa columns - reads/bases mapped to ALL "
                         "genera, not just Listeria. This is the relative-"
                         "abundance denominator for BOTH species- and "
                         "ST-level output, so the two share one common whole")
    ap.add_argument("-o", "--outdir", default=Path("listeria_summary"), type=Path,
                    help="output directory (default: listeria_summary)")
    ap.add_argument("-p", "--pattern", default="*_cov_lis_cleaned_annot_st.tsv*",
                    help="glob for input files (default: %(default)s)")
    ap.add_argument("-w", "--workers", type=int, default=0,
                    help="parallel workers; 0 = one per CPU (default)")
    ap.add_argument("--executor", choices=("process", "thread"), default="process",
                    help="process pool (default) or thread pool for small runs")
    ap.add_argument("--keep-rows", action="store_true",
                    help="also write every surviving non-zero row to disk")
    ap.add_argument("--keep-genus-suffix", action="store_true",
                    help="do NOT collapse GTDB genus suffixes (Listeria_A stays "
                         "distinct from Listeria); species epithet suffixes such "
                         "as monocytogenes_B are always collapsed")
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
    genus_totals = load_genus_totals(args.genus_totals)
    known_tps = sorted(meta["Timepoint"].dropna().unique().tolist())

    n_workers = args.workers or min(len(files), os.cpu_count() or 1)
    n_workers = max(1, min(n_workers, len(files)))
    use_awk = (not args.no_awk) and AWK is not None

    args.outdir.mkdir(parents=True, exist_ok=True)

    print(f"{len(files)} file(s) | {n_workers} {args.executor} worker(s) | "
          f"awk pre-filter: {'on' if use_awk else 'off'} | "
          f"nanostat: {args.nanostat} | genus-totals: {args.genus_totals}",
          file=sys.stderr)

    payloads = [
        (str(f), known_tps, args.keep_genus_suffix, use_awk, args.keep_rows)
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

    # nanostat is joined in as context only (see qc_per_file.tsv) - it no
    # longer feeds false_positive/unmapped or any relative-abundance calc.
    qc = qc.merge(nanostat, on=["basename", "Timepoint"], how="left")
    missing_nanostat = qc["num_of_reads"].isna()
    if missing_nanostat.any():
        print(f"NOTE: {int(missing_nanostat.sum())} sample(s) had no nanostat "
              f"match (joined on basename + Timepoint) - num_of_reads/"
              f"num_of_bases will be blank for these in qc_per_file.tsv; "
              f"relative abundance is unaffected since it no longer depends "
              f"on nanostat:", file=sys.stderr)
        for _, row in qc.loc[missing_nanostat, ["basename", "Timepoint"]].iterrows():
            print(f"  basename={row['basename']!r} Timepoint={row['Timepoint']!r}",
                  file=sys.stderr)

    # genus_totals supplies reads_all_taxa/bases_all_taxa - the actual
    # relative-abundance denominator for BOTH species- and ST-level output
    # (see add_relative_abundance()/append_offtarget_rows() calls below).
    # A missing match leaves reads_all_taxa/bases_all_taxa as NaN, which the
    # existing "denom.where(> 0)" logic in add_relative_abundance() already
    # turns into 0% (same forgiving behavior nanostat used to have).
    qc = qc.merge(genus_totals, on=["basename", "Timepoint"], how="left")
    missing_genus_totals = qc["reads_all_taxa"].isna()
    if missing_genus_totals.any():
        print(f"WARNING: {int(missing_genus_totals.sum())} sample(s) had no "
              f"genus-totals match (joined on basename + Timepoint) - "
              f"relative abundance will read as 0% for these:", file=sys.stderr)
        for _, row in qc.loc[missing_genus_totals, ["basename", "Timepoint"]].iterrows():
            print(f"  basename={row['basename']!r} Timepoint={row['Timepoint']!r}",
                  file=sys.stderr)

    # Sanity check: the coverage TSV's own all-taxa sum should never exceed
    # the externally supplied whole-community total - if it does, the two
    # inputs likely aren't from the same run/sample.
    over = qc[(qc["reads_all_taxa_covfile"] > qc["reads_all_taxa"].fillna(float("inf"))) |
              (qc["bases_all_taxa_covfile"] > qc["bases_all_taxa"].fillna(float("inf")))]
    if len(over):
        print(f"WARNING: {len(over)} sample(s) have MORE reads/bases in the "
              f"coverage TSV itself than the externally supplied "
              f"reads_all_taxa/bases_all_taxa total - check the two inputs "
              f"are from the same run:", file=sys.stderr)
        print(over[["basename", "Timepoint", "reads_all_taxa_covfile", "reads_all_taxa",
                   "bases_all_taxa_covfile", "bases_all_taxa"]].to_string(index=False),
              file=sys.stderr)

    samples = (
        qc[["basename", "Timepoint"]].drop_duplicates()
        .sort_values("Timepoint", key=lambda s: s.map(natural_key), kind="stable")
        .reset_index(drop=True)
    )

    species_long = pd.concat([r["species"] for r in results], ignore_index=True)
    st_long = pd.concat([r["st"] for r in results], ignore_index=True)
    st_sp_long = pd.concat([r["st_by_species"] for r in results], ignore_index=True)
    audit = (
        pd.concat([r["audit"] for r in results], ignore_index=True)
        .groupby(["raw_label", "species"], as_index=False)
        .agg(n_contigs=("n_contigs", "sum"), total_reads=("total_reads", "sum"))
        .sort_values(["species", "raw_label"])
    )
    st_audit = (
        pd.concat([r["st_audit"] for r in results], ignore_index=True)
        .groupby(["raw_st_label", "ST"], as_index=False, dropna=False)
        .agg(n_contigs=("n_contigs", "sum"), total_reads=("total_reads", "sum"))
        .sort_values(["ST", "raw_st_label"])
    )

    species_long = complete_grid(species_long, "species", TARGET_SPECIES, samples)
    st_long = complete_grid(st_long, "ST", TARGET_STS, samples)
    st_long["ST"] = "ST" + st_long["ST"].astype(str)
    if not st_sp_long.empty:
        st_sp_long["ST"] = "ST" + st_sp_long["ST"].astype(str)

    # Species level: TARGET rows (innocua/ivanovii/monocytogenes/welshimeri)
    # are relative abundance against the whole-community total
    # (reads_all_taxa/bases_all_taxa, from --genus-totals) - this is what
    # makes them nest under the genus/ST comparisons. The "Listeria
    # monocytogenes" target row now includes ALL monocytogenes reads
    # regardless of ST (untyped and off-target-ST included) - see
    # is_target_species in process_file(). false_positive is deliberately
    # kept on its own previous denominator (reads_all_listeria/
    # bases_all_listeria - all Listeria genus reads) instead, since it's a
    # typing/contamination-quality metric that's more meaningful scoped to
    # the Listeria reads actually recovered, not diluted by the whole
    # community; it now represents only OTHER (non-target) Listeria species,
    # since off-target-ST/untyped monocytogenes reads are counted as
    # monocytogenes targets, not false_positive, under this change. This
    # means the rows in species_summary_long.tsv do NOT all share one
    # denominator - see the per-row (not per-sample-sum) sanity check below.
    species_long = add_relative_abundance(
        species_long, qc, "reads_all_taxa", "bases_all_taxa")
    species_long = append_offtarget_rows(
        species_long, qc, key_col="species",
        rows_spec=[(FALSE_POSITIVE_LABEL, "reads_nontarget_species", "bases_nontarget_species")],
        reads_denom_col="reads_all_listeria", bases_denom_col="bases_all_listeria")

    # ST level: TARGET ST rows are relative abundance against the SAME
    # whole-community total as species-level targets, so target-ST% nests
    # under species-level Lm%. false_positive AND untyped are both kept on
    # their own previous denominator (reads_all_monocytogenes/
    # bases_all_monocytogenes) instead - "of the monocytogenes reads I
    # recovered, what fraction typed off-target / didn't type at all" - a
    # typing-quality question that's diluted into near-meaninglessness if
    # measured against the whole community.
    st_long = add_relative_abundance(
        st_long, qc, "reads_all_taxa", "bases_all_taxa")
    st_long = append_offtarget_rows(
        st_long, qc, key_col="ST",
        rows_spec=[(FALSE_POSITIVE_LABEL, "reads_nontarget_st", "bases_nontarget_st"),
                   (UNTYPED_LABEL, "reads_untyped_st", "bases_untyped_st")],
        reads_denom_col="reads_all_monocytogenes", bases_denom_col="bases_all_monocytogenes")

    # Rows in species_long/st_long now deliberately mix denominators (target
    # rows: whole-community; false_positive/untyped: Listeria-genus or
    # monocytogenes-only), so summing percentages ACROSS rows in a sample no
    # longer means anything - that check is gone. What must still hold is a
    # per-ROW bound: no single relative_abundance value can be outside
    # [0, 100], since every row is still numerator <= its own denominator by
    # construction.
    for name, df_ in (("species", species_long), ("ST", st_long)):
        bad_recon = df_[(df_["relative_abundance_reads"].round(2) > 100.0) |
                        (df_["relative_abundance_reads"].round(2) < 0.0) |
                        (df_["relative_abundance_bases"].round(2) > 100.0) |
                        (df_["relative_abundance_bases"].round(2) < 0.0)]
        if len(bad_recon):
            print(f"WARNING: {len(bad_recon)} row(s) have a relative_abundance "
                  f"value outside [0, 100] at {name} level - this should be "
                  f"impossible by construction; check the genus-totals file "
                  f"is from the same run as the coverage TSVs:", file=sys.stderr)
            print(bad_recon.to_string(index=False), file=sys.stderr)

    species_long = attach_metadata(species_long, meta)
    st_long = attach_metadata(st_long, meta)
    st_sp_long = attach_metadata(st_sp_long, meta)
    qc = attach_metadata(qc, meta)

    meta_keys = set(zip(meta["basename"], meta["Timepoint"]))
    unmatched = [(b, t) for b, t in zip(samples["basename"], samples["Timepoint"])
                 if (b, t) not in meta_keys]
    if unmatched:
        print(f"WARNING: {len(unmatched)} sample(s) had no metadata match "
              f"(joined on basename + Timepoint):", file=sys.stderr)
        for b, t in unmatched:
            print(f"  basename={b!r} Timepoint={t!r}", file=sys.stderr)

    id_cols = [c for c in species_long.columns
               if c not in ("species", "ST", "total_reads", "total_covbases",
                           "total_mapped_bases", "n_contigs",
                           "relative_abundance_reads", "relative_abundance_bases")]

    written = []

    def write(df: pd.DataFrame, name: str):
        p = args.outdir / name
        df.to_csv(p, sep="\t", index=False)
        written.append(p)

    write(species_long.sort_values(["Timepoint", "basename", "species"],
                                   key=lambda s: s.map(natural_key)
                                   if s.name == "Timepoint" else s),
          "species_summary_long.tsv")
    write(to_wide(species_long, "species", "total_reads", id_cols),
          "species_reads_wide.tsv")
    write(to_wide(species_long, "species", "total_covbases", id_cols),
          "species_covbases_wide.tsv")
    write(to_wide(species_long, "species", "total_mapped_bases", id_cols),
          "species_mapped_bases_wide.tsv")
    write(to_wide(species_long, "species", "relative_abundance_reads", id_cols),
          "species_relative_abundance_reads_wide.tsv")
    write(to_wide(species_long, "species", "relative_abundance_bases", id_cols),
          "species_relative_abundance_bases_wide.tsv")

    write(st_long.sort_values(["Timepoint", "basename", "ST"],
                              key=lambda s: s.map(natural_key)
                              if s.name == "Timepoint" else s),
          "st_summary_long.tsv")
    write(to_wide(st_long, "ST", "total_reads", id_cols), "st_reads_wide.tsv")
    write(to_wide(st_long, "ST", "total_covbases", id_cols), "st_covbases_wide.tsv")
    write(to_wide(st_long, "ST", "total_mapped_bases", id_cols), "st_mapped_bases_wide.tsv")
    write(to_wide(st_long, "ST", "relative_abundance_reads", id_cols),
          "st_relative_abundance_reads_wide.tsv")
    write(to_wide(st_long, "ST", "relative_abundance_bases", id_cols),
          "st_relative_abundance_bases_wide.tsv")
    if not st_sp_long.empty:
        write(st_sp_long, "st_by_species_long.tsv")

    write(qc, "qc_per_file.tsv")
    write(audit, "species_label_audit.tsv")
    write(st_audit, "st_label_audit.tsv")

    # Grand-total denominators mirror the per-sample mixed design above:
    # target rows against the whole-community grand total, false_positive/
    # untyped rows against their own (Listeria-genus or monocytogenes-only)
    # grand total, so a sample's contribution to each row uses the same
    # denominator basis it used per-sample.
    grand_reads_taxa = qc["reads_all_taxa"].sum()
    grand_bases_taxa = qc["bases_all_taxa"].sum()
    grand_reads_listeria = qc["reads_all_listeria"].sum()
    grand_bases_listeria = qc["bases_all_listeria"].sum()
    grand_reads_mono = qc["reads_all_monocytogenes"].sum()
    grand_bases_mono = qc["bases_all_monocytogenes"].sum()

    def grand_totals(df_: pd.DataFrame, key: str,
                     grand_reads: float, grand_bases: float) -> pd.DataFrame:
        g = (df_.groupby(key, as_index=False)
             .agg(total_reads=("total_reads", "sum"),
                  total_covbases=("total_covbases", "sum"),
                  total_mapped_bases=("total_mapped_bases", "sum")))
        g["relative_abundance_reads"] = (
            (g["total_reads"] / grand_reads * 100).round(4) if grand_reads else 0.0
        )
        g["relative_abundance_bases"] = (
            (g["total_mapped_bases"] / grand_bases * 100).round(4) if grand_bases else 0.0
        )
        return g

    is_sp_offtarget = species_long["species"] == FALSE_POSITIVE_LABEL
    grand_sp = pd.concat([
        grand_totals(species_long[~is_sp_offtarget], "species", grand_reads_taxa, grand_bases_taxa),
        grand_totals(species_long[is_sp_offtarget], "species", grand_reads_listeria, grand_bases_listeria),
    ], ignore_index=True)

    is_st_offtarget = st_long["ST"].isin([FALSE_POSITIVE_LABEL, UNTYPED_LABEL])
    grand_st = pd.concat([
        grand_totals(st_long[~is_st_offtarget], "ST", grand_reads_taxa, grand_bases_taxa),
        grand_totals(st_long[is_st_offtarget], "ST", grand_reads_mono, grand_bases_mono),
    ], ignore_index=True)

    write(grand_sp, "grand_totals_species.tsv")
    write(grand_st, "grand_totals_st.tsv")

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

    print("\n== Grand total reads / mapped bases per species (targets: % of "
          "whole-community; false_positive: % of all-Listeria) ==")
    print(grand_sp.to_string(index=False))
    print("\n== Grand total reads / mapped bases per ST (targets: % of "
          "whole-community; false_positive/untyped: % of monocytogenes-only) ==")
    print(grand_st.to_string(index=False))

    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
