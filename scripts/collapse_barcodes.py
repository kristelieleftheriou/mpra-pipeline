#!/usr/bin/env python3
"""
collapse_barcodes.py -- error-correct barcodes with starcode.

Runs `starcode -d <dist> --print-clusters` on every *_barcode_counts.tsv found
in each per-sample directory and writes, next to the input:

  <name>_barcode_counts.starcode_d<dist>.tsv   centroid <TAB> count <TAB> members
  <name>_barcode_counts.starcode_d<dist>.collapsed.tsv   barcode <TAB> count

Note: clustering is done per sample, so a centroid in one sample is not
guaranteed to be the centroid in another. For 1-2 bp Illumina errors on
well-covered barcodes this is not an issue; if you need exact cross-sample
identity, snap the centroids to your PBA barcodes downstream.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile


def run_starcode(binary, pairs, dist, threads):
    """Cluster (sequence, count) pairs; return [(centroid, count, members)]."""
    if not pairs:
        return []
    with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
        in_path = fh.name
        for seq, count in pairs:
            fh.write(f"{seq}\t{count}\n")
    out_path = in_path + ".out"
    try:
        subprocess.run([binary, "-d", str(dist), "--print-clusters",
                        "-t", str(threads), "-i", in_path, "-o", out_path,
                        "--quiet"], check=True)
        clusters = []
        with open(out_path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 3:
                    clusters.append((parts[0], int(parts[1]), parts[2].split(",")))
        return clusters
    finally:
        for path in (in_path, out_path):
            if os.path.exists(path):
                os.unlink(path)


def read_counts(path):
    pairs = []
    with open(path) as fh:
        fh.readline()  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    pairs.append((parts[0], int(parts[1])))
                except ValueError:
                    continue
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True,
                    help="directory of per-sample subdirectories")
    ap.add_argument("--samples", nargs="+", default=None)
    ap.add_argument("--dist", type=int, default=2)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--starcode", default="starcode")
    args = ap.parse_args()

    binary = shutil.which(args.starcode)
    if binary is None:
        sys.exit(f"ERROR: starcode not found ('{args.starcode}'). "
                 f"Activate the starcode env or pass --starcode /full/path.")

    samples = args.samples or sorted(
        d for d in os.listdir(args.indir)
        if os.path.isdir(os.path.join(args.indir, d)))
    if not samples:
        sys.exit(f"No sample directories in {args.indir}")

    for sample in samples:
        sample_dir = os.path.join(args.indir, sample)
        inputs = sorted(glob.glob(os.path.join(sample_dir, "*_barcode_counts.tsv")))
        if not inputs:
            print(f"[collapse] WARN: nothing to do in {sample_dir}", file=sys.stderr)
            continue
        for path in inputs:
            base = os.path.basename(path)[:-len(".tsv")]
            pairs = read_counts(path)
            clusters = run_starcode(binary, pairs, args.dist, args.threads)

            verbose = os.path.join(sample_dir, f"{base}.starcode_d{args.dist}.tsv")
            collapsed = os.path.join(sample_dir,
                                     f"{base}.starcode_d{args.dist}.collapsed.tsv")
            with open(verbose, "w") as fh:
                fh.write("centroid\tcount\tmembers\n")
                for centroid, count, members in clusters:
                    fh.write(f"{centroid}\t{count}\t{','.join(members)}\n")
            with open(collapsed, "w") as fh:
                fh.write("barcode\tcount\n")
                for centroid, count, _ in clusters:
                    fh.write(f"{centroid}\t{count}\n")

            print(f"[collapse] {sample}/{base}: {len(pairs):,} -> "
                  f"{len(clusters):,} centroids "
                  f"({100 * (len(pairs) - len(clusters)) / max(1, len(pairs)):.1f}% reduction)")


if __name__ == "__main__":
    main()
