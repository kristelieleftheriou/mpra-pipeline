#!/usr/bin/env python3
"""
qc_pba.py

Build a clean, deduplicated promoter-barcode association (PBA) table
from raw (promoter, barcode, n_reads) counts.

For each barcode:
  - Drop if total reads < --min-reads.
  - If it maps to multiple promoters, keep only if the dominant promoter's
    fraction of reads is >= --dominance; otherwise flag as ambiguous.

Reference entries whose name starts with --empty-prefix (default "Empty")
are tracked separately in the summary as empty-vector / empty-intron
assignments. They still appear in the clean PBA table just like any
other promoter, so the user gets the BC -> Empty_* assignment.

Outputs:
  --out-pba              promoter, barcode, n_reads_dominant, n_reads_total, dominance
  --out-ambiguous        barcodes that failed the dominance test
  --out-promoter-stats   per-promoter coverage
  --out-summary          one-page QC summary (incl. empty/promoter split)
"""
import argparse
import collections
import csv
import statistics


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True,
                    help="TSV with columns promoter, barcode, n_reads.")
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("--dominance", type=float, default=0.9,
                    help="Minimum fraction of reads on the dominant promoter "
                         "for a barcode to be kept (when >1 promoter is seen).")
    ap.add_argument("--n-promoters-expected", type=int, default=1000)
    ap.add_argument("--empty-prefix", default="Empty",
                    help="Reference names starting with this prefix are "
                         "treated as empty-vector references in the summary "
                         "(case-sensitive). Default: 'Empty'.")
    ap.add_argument("--out-pba",            required=True)
    ap.add_argument("--out-ambiguous",      required=True)
    ap.add_argument("--out-promoter-stats", required=True)
    ap.add_argument("--out-summary",        required=True)
    args = ap.parse_args()

    def is_empty(ref):
        return ref.startswith(args.empty_prefix)

    # ---- Group counts by barcode
    bc_to_prom: dict[str, collections.Counter] = collections.defaultdict(
        collections.Counter)
    with open(args.counts) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            bc_to_prom[row["barcode"]][row["promoter"]] += int(row["n_reads"])

    n_bc_total = len(bc_to_prom)
    n_bc_lowreads = 0
    n_bc_ambiguous = 0
    pba_rows = []
    amb_rows = []

    for bc, prom_counts in bc_to_prom.items():
        total = sum(prom_counts.values())
        if total < args.min_reads:
            n_bc_lowreads += 1
            continue
        top_prom, top_n = prom_counts.most_common(1)[0]
        dom = top_n / total
        if dom < args.dominance and len(prom_counts) > 1:
            n_bc_ambiguous += 1
            promoters_str = ";".join(f"{p}:{c}"
                                     for p, c in prom_counts.most_common())
            amb_rows.append([bc, total, len(prom_counts), top_prom,
                             top_n, f"{dom:.3f}", promoters_str])
            continue
        pba_rows.append([top_prom, bc, top_n, total, f"{dom:.3f}"])

    n_bc_kept = len(pba_rows)

    # ---- Clean PBA table
    with open(args.out_pba, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["promoter", "barcode", "n_reads_dominant",
                    "n_reads_total", "dominance"])
        for row in sorted(pba_rows, key=lambda r: (r[0], -r[2])):
            w.writerow(row)

    # ---- Ambiguous barcodes
    with open(args.out_ambiguous, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["barcode", "n_reads_total", "n_promoters",
                    "top_promoter", "top_n_reads", "dominance",
                    "all_promoters"])
        for row in sorted(amb_rows, key=lambda r: -r[1]):
            w.writerow(row)

    # ---- Per-promoter stats
    prom_to_reads: dict[str, list[int]] = collections.defaultdict(list)
    for prom, _bc, n_dom, _n_tot, _dom in pba_rows:
        prom_to_reads[prom].append(n_dom)

    with open(args.out_promoter_stats, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["promoter", "category", "n_barcodes", "total_reads",
                    "median_reads_per_bc", "max_reads_per_bc"])
        for prom, ns in sorted(prom_to_reads.items()):
            cat = "empty" if is_empty(prom) else "promoter"
            w.writerow([prom, cat, len(ns), sum(ns),
                        f"{statistics.median(ns):.1f}", max(ns)])

    # ---- Summary
    bc_kept_promoter = sum(1 for r in pba_rows if not is_empty(r[0]))
    bc_kept_empty    = sum(1 for r in pba_rows if     is_empty(r[0]))
    reads_promoter   = sum(r[3] for r in pba_rows if not is_empty(r[0]))  # n_reads_total
    reads_empty      = sum(r[3] for r in pba_rows if     is_empty(r[0]))

    promoters_only = {p: ns for p, ns in prom_to_reads.items() if not is_empty(p)}
    empties_only   = {p: ns for p, ns in prom_to_reads.items() if     is_empty(p)}

    with open(args.out_summary, "w") as fh:
        # Settings
        fh.write(f"min_reads_threshold\t{args.min_reads}\n")
        fh.write(f"dominance_threshold\t{args.dominance}\n")
        fh.write(f"empty_prefix\t{args.empty_prefix}\n")

        # Barcode-level outcome
        fh.write(f"barcodes_observed\t{n_bc_total}\n")
        fh.write(f"barcodes_dropped_low_reads\t{n_bc_lowreads}\n")
        fh.write(f"barcodes_dropped_ambiguous\t{n_bc_ambiguous}\n")
        fh.write(f"barcodes_kept_total\t{n_bc_kept}\n")
        fh.write(f"barcodes_kept_promoter\t{bc_kept_promoter}\n")
        fh.write(f"barcodes_kept_empty\t{bc_kept_empty}\n")

        # Read-level totals (across kept barcodes only)
        fh.write(f"reads_in_kept_promoter_bcs\t{reads_promoter}\n")
        fh.write(f"reads_in_kept_empty_bcs\t{reads_empty}\n")
        if reads_promoter + reads_empty > 0:
            frac_empty = reads_empty / (reads_promoter + reads_empty)
            fh.write(f"empty_read_fraction\t{frac_empty:.4f}\n")

        # Promoter coverage
        fh.write(f"promoters_detected\t{len(promoters_only)}\n")
        fh.write(f"promoters_expected\t{args.n_promoters_expected}\n")
        fh.write(f"promoters_missing\t"
                 f"{args.n_promoters_expected - len(promoters_only)}\n")
        fh.write(f"empty_refs_detected\t{len(empties_only)}\n")
        for ref, ns in sorted(empties_only.items()):
            fh.write(f"empty_ref_{ref}_n_barcodes\t{len(ns)}\n")
            fh.write(f"empty_ref_{ref}_total_reads\t{sum(ns)}\n")

        # Per-promoter BC distribution (excluding empty refs)
        if promoters_only:
            bcs_per_prom = sorted(len(v) for v in promoters_only.values())
            fh.write(f"barcodes_per_promoter_min\t{bcs_per_prom[0]}\n")
            fh.write(f"barcodes_per_promoter_p10\t"
                     f"{bcs_per_prom[max(0, len(bcs_per_prom) // 10)]}\n")
            fh.write(f"barcodes_per_promoter_median\t"
                     f"{statistics.median(bcs_per_prom):.1f}\n")
            fh.write(f"barcodes_per_promoter_p90\t"
                     f"{bcs_per_prom[min(len(bcs_per_prom) - 1, (9 * len(bcs_per_prom)) // 10)]}\n")
            fh.write(f"barcodes_per_promoter_max\t{bcs_per_prom[-1]}\n")


if __name__ == "__main__":
    main()

