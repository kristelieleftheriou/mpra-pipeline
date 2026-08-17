#!/usr/bin/env python3
"""
extract_pba.py

Build a raw promoter-barcode association (PBA) table by joining:
  - bowtie2 BAM of promoter-side reads (one alignment per read ID)
  - FASTQ of barcode-side reads (read IDs match the BAM)

The barcode read is structured (post-cutadapt) as:
    [BC_rc (12 bp)] [revcomp(universal)] [revcomp(promoter ...)]

We take the first --bc-len bases as the barcode (in revcomp orientation
relative to the canonical top strand), check that the next
--univ-anchor-len bases match the start of revcomp(universal) (allowing
--univ-anchor-mm mismatches), and optionally reverse-complement the
barcode so it is reported in canonical (top-strand) orientation.
"""
import argparse
import collections
import csv
import gzip
import sys

import pysam


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def hamming_le(a: str, b: str, max_mm: int) -> bool:
    """True if Hamming distance between equal-length a, b is <= max_mm."""
    if len(a) != len(b):
        return False
    mm = 0
    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > max_mm:
                return False
    return True


def iter_fastq(path: str):
    """Yield (read_id, seq) from a (possibly gzipped) FASTQ."""
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()  # +
            fh.readline()  # qual
            rid = h[1:].split()[0]
            yield rid, seq


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--barcode-fq", required=True)
    ap.add_argument("--bc-len", type=int, default=12)
    ap.add_argument("--mapq-min", type=int, default=20)
    ap.add_argument("--rc-barcode", action="store_true",
                    help="Reverse-complement the extracted barcode so it is "
                         "reported in canonical (top-strand) orientation. "
                         "Recommended.")
    ap.add_argument("--universal", required=True,
                    help="Universal sequence in canonical (top-strand) "
                         "orientation. The first --univ-anchor-len bases of "
                         "its reverse complement are required to follow the "
                         "barcode on the barcode read.")
    ap.add_argument("--univ-anchor-len", type=int, default=10)
    ap.add_argument("--univ-anchor-mm", type=int, default=1)
    ap.add_argument("--out-pba",    required=True,
                    help="Per-read TSV: read_id, promoter, barcode, mapq")
    ap.add_argument("--out-counts", required=True,
                    help="Aggregated TSV: promoter, barcode, n_reads")
    ap.add_argument("--out-stats",  required=True)
    args = ap.parse_args()

    universal = args.universal.upper()
    univ_rc = rc(universal)
    anchor = univ_rc[:args.univ_anchor_len]
    print(f"[INFO] universal:           {universal}", file=sys.stderr)
    print(f"[INFO] revcomp(universal):  {univ_rc}", file=sys.stderr)
    print(f"[INFO] anchor on barcode read (after BC): {anchor} "
          f"(allow {args.univ_anchor_mm} mm)", file=sys.stderr)

    # ---- 1. Read barcode FASTQ into dict keyed by read ID
    bc_map: dict[str, str] = {}
    n_bc = 0
    n_short = 0
    n_with_n = 0
    n_anchor_fail = 0
    need_len = args.bc_len + args.univ_anchor_len

    for rid, seq in iter_fastq(args.barcode_fq):
        n_bc += 1
        if len(seq) < need_len:
            n_short += 1
            continue
        bc = seq[:args.bc_len].upper()
        if "N" in bc:
            n_with_n += 1
            continue
        observed = seq[args.bc_len:args.bc_len + args.univ_anchor_len].upper()
        if not hamming_le(observed, anchor, args.univ_anchor_mm):
            n_anchor_fail += 1
            continue
        if args.rc_barcode:
            bc = rc(bc)
        bc_map[rid] = bc

    print(f"[INFO] barcode reads: total={n_bc} short={n_short} "
          f"N_in_BC={n_with_n} anchor_fail={n_anchor_fail} "
          f"usable={len(bc_map)}", file=sys.stderr)

    # ---- 2. Walk BAM, write per-read PBA, aggregate counts
    n_aln = 0
    n_unmapped = 0
    n_secondary = 0
    n_lowmapq = 0
    n_no_bc = 0
    n_kept = 0
    counts: collections.Counter = collections.Counter()

    with pysam.AlignmentFile(args.bam, "rb") as bam, \
            open(args.out_pba, "w", newline="") as outf:
        w = csv.writer(outf, delimiter="\t")
        w.writerow(["read_id", "promoter", "barcode", "mapq"])

        for read in bam:
            n_aln += 1
            if read.is_unmapped:
                n_unmapped += 1
                continue
            if read.is_secondary or read.is_supplementary:
                n_secondary += 1
                continue
            if read.mapping_quality < args.mapq_min:
                n_lowmapq += 1
                continue
            rid = read.query_name
            bc = bc_map.get(rid)
            if bc is None:
                n_no_bc += 1
                continue
            promoter = read.reference_name
            w.writerow([rid, promoter, bc, read.mapping_quality])
            counts[(promoter, bc)] += 1
            n_kept += 1

    # ---- 3. Aggregated (promoter, barcode, n_reads) table
    with open(args.out_counts, "w", newline="") as outf:
        w = csv.writer(outf, delimiter="\t")
        w.writerow(["promoter", "barcode", "n_reads"])
        for (p, b), c in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
            w.writerow([p, b, c])

    # ---- 4. Stats
    with open(args.out_stats, "w") as outf:
        outf.write(f"barcode_reads_total\t{n_bc}\n")
        outf.write(f"barcode_reads_too_short\t{n_short}\n")
        outf.write(f"barcode_reads_with_N_in_bc\t{n_with_n}\n")
        outf.write(f"barcode_reads_anchor_fail\t{n_anchor_fail}\n")
        outf.write(f"barcode_reads_usable\t{len(bc_map)}\n")
        outf.write(f"alignments_total\t{n_aln}\n")
        outf.write(f"alignments_unmapped\t{n_unmapped}\n")
        outf.write(f"alignments_secondary_supplementary\t{n_secondary}\n")
        outf.write(f"alignments_below_mapq\t{n_lowmapq}\n")
        outf.write(f"alignments_no_matching_barcode\t{n_no_bc}\n")
        outf.write(f"alignments_kept\t{n_kept}\n")
        outf.write(f"unique_promoter_barcode_pairs\t{len(counts)}\n")

    print(f"[INFO] kept {n_kept} read events; "
          f"{len(counts)} unique (promoter, barcode) pairs.", file=sys.stderr)


if __name__ == "__main__":
    main()

