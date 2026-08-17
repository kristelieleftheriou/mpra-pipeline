#!/usr/bin/env python3
"""
build_count_matrix.py -- combine per-sample counts into cross-sample matrices.

Replaces the earlier aggregate_counts.py / aggregate_starcode_counts.py pair:
pass --suffix to switch between raw and starcode-collapsed inputs.

Inputs, per sample directory:
  <name>_barcode_counts[.<suffix>.collapsed].tsv   barcode <TAB> count
  qc_summary.json

Outputs, in --outdir:
  <name>_count_matrix[.<suffix>].tsv   rows = barcodes, cols = samples
  qc_summary_all.tsv                   one row per sample
"""
import argparse
import glob
import json
import os
import sys
from collections import defaultdict


def read_counts(path):
    counts = {}
    if not os.path.exists(path):
        return counts
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    counts[parts[0]] = int(parts[1])
                except ValueError:
                    continue
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--samples", nargs="+", default=None,
                    help="sample subdirectory names, in the order you want "
                         "the columns (default: all, sorted)")
    ap.add_argument("--suffix", default=None,
                    help="e.g. starcode_d2 to read the collapsed files")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    samples = args.samples or sorted(
        d for d in os.listdir(args.indir)
        if os.path.isdir(os.path.join(args.indir, d)))
    if not samples:
        sys.exit(f"No sample directories in {args.indir}")
    print(f"[matrix] {len(samples)} sample(s)")

    tail = f".{args.suffix}.collapsed.tsv" if args.suffix else ".tsv"

    # Which count files exist (one per marker in the design) -- take the union.
    names = set()
    for sample in samples:
        for path in glob.glob(os.path.join(args.indir, sample,
                                           f"*_barcode_counts{tail}")):
            names.add(os.path.basename(path).split("_barcode_counts")[0])
    if not names:
        sys.exit(f"No *_barcode_counts{tail} files found under {args.indir}")

    for name in sorted(names):
        table = defaultdict(lambda: [0] * len(samples))
        for i, sample in enumerate(samples):
            path = os.path.join(args.indir, sample,
                                f"{name}_barcode_counts{tail}")
            if not os.path.exists(path):
                print(f"[matrix] WARN: missing {path}", file=sys.stderr)
                continue
            for barcode, count in read_counts(path).items():
                table[barcode][i] = count

        stem = f"{name}_count_matrix" + (f".{args.suffix}" if args.suffix else "")
        out = os.path.join(args.outdir, stem + ".tsv")
        with open(out, "w") as fh:
            fh.write("barcode\t" + "\t".join(samples) + "\n")
            for barcode, counts in sorted(table.items(), key=lambda kv: -sum(kv[1])):
                fh.write(barcode + "\t" + "\t".join(map(str, counts)) + "\n")
        print(f"[matrix] wrote {out}  ({len(table):,} barcodes)")

    # ---------------- QC summary (independent of --suffix) ----------------
    summaries = []
    for sample in samples:
        path = os.path.join(args.indir, sample, "qc_summary.json")
        if os.path.exists(path):
            with open(path) as fh:
                summaries.append(json.load(fh))
        else:
            print(f"[matrix] WARN: missing {path}", file=sys.stderr)

    if summaries:
        first = ["sample", "r1", "total_reads"]
        cols = first + sorted(
            {k for d in summaries for k in d} - set(first))
        out = os.path.join(args.outdir, "qc_summary_all.tsv")
        with open(out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for d in summaries:
                fh.write("\t".join(str(d.get(c, "")) for c in cols) + "\n")
        print(f"[matrix] wrote {out}")


if __name__ == "__main__":
    main()
