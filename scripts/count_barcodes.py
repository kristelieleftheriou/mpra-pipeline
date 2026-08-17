#!/usr/bin/env python3
"""
count_barcodes.py -- count barcodes in one R1 FASTQ.

Read layout (positions fixed, taken from the design YAML):

    [anchor][barcode][marker][rest]

Each read falls into exactly one category:

  <marker_name>   anchor matched and the marker matched -> barcode counted
                  (one category per entry in `markers`, typically 'library')
  other_marker    anchor matched, no marker matched
  no_anchor       anchor did not match
  too_short       read shorter than anchor + barcode + marker

Outputs, under <outdir>/<sample>/:
  <marker>_barcode_counts.tsv   barcode <TAB> count, sorted by count
  qc_summary.json               read counts and percentages per category
  top_unassigned.fa             most abundant unassigned 60-mers (BLAST these)
"""
import argparse
import gzip
import json
import os
from collections import Counter

import yaml


def n_mismatches(observed: str, expected: str) -> int:
    """Mismatches between equal-length strings; N is a wildcard."""
    return sum(1 for a, b in zip(observed, expected)
               if a != b and a != "N" and b != "N")


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def iter_reads(path):
    with open_fastq(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                yield line.rstrip("\n")


class Design:
    def __init__(self, d):
        self.anchor = d["anchor"].upper()
        self.anchor_max_mm = int(d.get("anchor_max_mm", 2))
        self.bc_len = int(d["barcode_length"])
        self.markers = {k: v.upper() for k, v in (d.get("markers") or {}).items()}
        self.marker_max_mm = int(d.get("marker_max_mm", 0))
        self.top_unassigned = int(d.get("top_unassigned", 100))
        self.bc_start = len(self.anchor)
        self.bc_end = self.bc_start + self.bc_len
        self.marker_len = max((len(m) for m in self.markers.values()), default=0)
        self.min_len = self.bc_end + self.marker_len

    def classify(self, seq):
        """Return (category, payload). Payload is a barcode or a 60-mer."""
        if len(seq) < self.min_len:
            return "too_short", None
        if n_mismatches(seq[:len(self.anchor)], self.anchor) > self.anchor_max_mm:
            return "no_anchor", seq[:60]
        barcode = seq[self.bc_start:self.bc_end]
        observed = seq[self.bc_end:self.bc_end + self.marker_len]
        for name, marker in self.markers.items():
            if not marker or n_mismatches(observed[:len(marker)], marker) <= self.marker_max_mm:
                return name, barcode
        return "other_marker", seq[:60]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--r1", required=True, help="R1 FASTQ (.gz ok)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True,
                    help="results go to <outdir>/<sample>/")
    ap.add_argument("--design", required=True, help="design YAML")
    args = ap.parse_args()

    with open(args.design) as fh:
        design = Design(yaml.safe_load(fh))

    sample_dir = os.path.join(args.outdir, args.sample)
    os.makedirs(sample_dir, exist_ok=True)

    counts = {name: Counter() for name in design.markers} or {"library": Counter()}
    unassigned = Counter()
    n_cat = Counter()
    total = 0

    for seq in iter_reads(args.r1):
        total += 1
        cat, payload = design.classify(seq)
        n_cat[cat] += 1
        if cat in counts:
            counts[cat][payload] += 1
        elif payload is not None:
            unassigned[payload] += 1

    for name, counter in counts.items():
        with open(os.path.join(sample_dir, f"{name}_barcode_counts.tsv"), "w") as fh:
            fh.write("barcode\tcount\n")
            for bc, c in counter.most_common():
                fh.write(f"{bc}\t{c}\n")

    with open(os.path.join(sample_dir, "top_unassigned.fa"), "w") as fh:
        for i, (seq, c) in enumerate(unassigned.most_common(design.top_unassigned), 1):
            fh.write(f">{args.sample}_unassigned{i:03d}_n{c}\n{seq}\n")

    pct = lambda x: round(100.0 * x / total, 3) if total else 0.0
    summary = {"sample": args.sample,
               "r1": os.path.basename(args.r1),
               "total_reads": total}
    for name, counter in counts.items():
        summary[f"{name}_reads"] = sum(counter.values())
        summary[f"{name}_unique_barcodes"] = len(counter)
        summary[f"pct_{name}"] = pct(sum(counter.values()))
    for cat in ("other_marker", "no_anchor", "too_short"):
        summary[f"{cat}_reads"] = n_cat[cat]
        summary[f"pct_{cat}"] = pct(n_cat[cat])
    summary["pct_assigned"] = round(
        sum(summary[f"pct_{name}"] for name in counts), 3)

    with open(os.path.join(sample_dir, "qc_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"[count_barcodes] {args.sample}: {total:,} reads, "
          f"{summary['pct_assigned']}% assigned")
    for name in counts:
        print(f"    {name:<12} {summary[f'{name}_reads']:>12,}  "
              f"({summary[f'pct_{name}']}%)  "
              f"unique={summary[f'{name}_unique_barcodes']:,}")
    for cat in ("other_marker", "no_anchor", "too_short"):
        print(f"    {cat:<12} {n_cat[cat]:>12,}  ({summary[f'pct_{cat}']}%)")


if __name__ == "__main__":
    main()
