#!/usr/bin/env python3
"""
add_taxa.py

Annotate samtools-coverage-style files with GTDB taxonomy, matching the
assembly accession embedded in each coverage file's `rname` column to the
accession in a GTDB-style taxonomy file. The taxonomy is added AS IS (the
lineage string from the taxonomy file, not split into ranks, no NA filler).

Accession matching:
  - taxonomy accessions carry a prefix:  RS_GCF_051534235.1 / GB_GCA_950287535.1
  - coverage rnames carry a suffix:      GCF_006385165.1~NZ_VFJA01000002.1
  The RS_/GB_ prefix and the ~contig suffix are stripped, then the bare
  accession is matched. If the exact versioned accession is missing, the
  version number (.N) is ignored as a fallback.

Usage:
    python3 add_taxa.py bac_taxa.tsv coverage_dir/
    python3 add_taxa.py bac_taxa.tsv coverage_dir/ -o annotated_dir/
    python3 add_taxa.py bac_taxa.tsv coverage_dir/ --recursive

Input assumptions:
    - taxonomy file: tab-separated, no header; first column is the accession,
      the remaining column(s) hold the taxonomy (e.g. one column
      "d__Bacteria;p__...;s__..." or several tab-separated columns).
    - coverage files: samtools-coverage-style, first line is a "#rname ..." header.
"""

import sys
import csv
import argparse
from pathlib import Path

ANNOTATED_SUFFIX = "_ann.tsv"


def strip_source_prefix(accession):
    """Drop GTDB's RS_/GB_ source-database prefix, if present."""
    if accession.startswith("RS_") or accession.startswith("GB_"):
        return accession[3:]
    return accession


def base_accession(accession):
    """Strip a trailing .N version suffix, e.g. GCF_000123456.2 -> GCF_000123456."""
    prefix, sep, suffix = accession.rpartition(".")
    if sep and suffix.isdigit():
        return prefix
    return accession


def load_taxonomy(path):
    """
    Returns two dicts mapping accession -> taxonomy string (as in the file):
      exact_map (versioned accession) and base_map (version-agnostic fallback).
    """
    exact_map = {}
    base_map = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if not row or not row[0].strip():
                continue
            accession = strip_source_prefix(row[0].strip())
            # Everything after the accession, kept as is. If the lineage is
            # spread over several tab columns, join them with ";".
            taxonomy = ";".join(v.strip() for v in row[1:] if v.strip())
            exact_map[accession] = taxonomy
            base_map.setdefault(base_accession(accession), taxonomy)
    return exact_map, base_map


def find_coverage_files(input_dir, pattern, taxonomy_path, recursive=False):
    """Files in input_dir matching pattern, minus prior outputs and the taxonomy file itself."""
    glob_fn = input_dir.rglob if recursive else input_dir.glob
    return sorted(
        p for p in glob_fn(pattern)
        if p.is_file()
        and not p.name.endswith(ANNOTATED_SUFFIX)
        and p.resolve() != taxonomy_path
    )


def annotate_file(coverage_path, exact_map, base_map, output_path):
    """Annotate a single coverage file. Returns (n_rows, n_exact, n_fallback, missing_set)."""
    n_rows = n_exact = n_fallback = 0
    missing = set()

    with open(coverage_path) as fin, open(output_path, "w") as fout:
        for i, line in enumerate(fin):
            fields = line.rstrip("\n").split()
            if not fields:
                continue

            if i == 0 and fields[0].startswith("#"):
                fout.write("\t".join(fields + ["accession", "taxonomy"]) + "\n")
                continue

            n_rows += 1
            accession = fields[0].split("~", 1)[0]

            taxonomy = exact_map.get(accession)
            if taxonomy is not None:
                n_exact += 1
            else:
                taxonomy = base_map.get(base_accession(accession))
                if taxonomy is not None:
                    n_fallback += 1
                else:
                    missing.add(accession)
                    taxonomy = ""  # no match: leave empty rather than writing NA

            fout.write("\t".join(fields + [accession, taxonomy]) + "\n")

    return n_rows, n_exact, n_fallback, missing


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("taxonomy_file", help="taxonomy TSV: accession + taxonomy")
    ap.add_argument("input_dir", type=Path, help="directory containing samtools coverage-style .tsv files")
    ap.add_argument("-o", "--output-dir", type=Path, default=None,
                    help="directory to write annotated files to (default: same as input_dir)")
    ap.add_argument("--pattern", default="*.tsv",
                    help="glob pattern for coverage files inside input_dir (default: *.tsv)")
    ap.add_argument("--recursive", action="store_true",
                    help="also search subdirectories of input_dir")
    args = ap.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"[error] input directory not found: {args.input_dir}")

    output_dir = args.output_dir or args.input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    exact_map, base_map = load_taxonomy(args.taxonomy_file)
    print(f"[info] loaded {len(exact_map)} taxonomy records from {args.taxonomy_file}", file=sys.stderr)

    taxonomy_path = Path(args.taxonomy_file).resolve()
    coverage_files = find_coverage_files(args.input_dir, args.pattern, taxonomy_path, args.recursive)
    if not coverage_files:
        sys.exit(f"[error] no matching coverage files found in {args.input_dir} (pattern: {args.pattern})")

    print(f"[info] found {len(coverage_files)} coverage file(s) in {args.input_dir}", file=sys.stderr)

    total_rows = total_exact = total_fallback = 0
    all_missing = set()

    for coverage_path in coverage_files:
        output_path = output_dir / (coverage_path.stem + ANNOTATED_SUFFIX)
        n_rows, n_exact, n_fallback, missing = annotate_file(coverage_path, exact_map, base_map, output_path)

        print(f"[info] {coverage_path.name}: matched {n_exact + n_fallback}/{n_rows} rows "
              f"({n_exact} exact, {n_fallback} version-agnostic fallback) -> {output_path.name}",
              file=sys.stderr)
        if missing:
            example = ", ".join(sorted(missing)[:5])
            print(f"[warn] {coverage_path.name}: {len(missing)} distinct accession(s) had no taxonomy match, "
                  f"e.g.: {example}", file=sys.stderr)

        total_rows += n_rows
        total_exact += n_exact
        total_fallback += n_fallback
        all_missing |= missing

    print(f"[info] done: {len(coverage_files)} file(s) written to {output_dir}, "
          f"{total_exact + total_fallback}/{total_rows} rows matched overall "
          f"({total_exact} exact, {total_fallback} fallback); "
          f"{len(all_missing)} distinct accession(s) never matched",
          file=sys.stderr)


if __name__ == "__main__":
    main()