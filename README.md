# mpra-pipeline

Analysis pipeline for a barcoded promoter library assayed as a massively
parallel reporter assay (MPRA). Two stages:

1. **PBA** — promoter–barcode association. Paired-end sequencing of the plasmid
   library links each 12 bp barcode to the promoter it sits next to.
2. **Counting and activity** — barcode counting in RNA and plasmid (pDNA)
   samples, then RNA/DNA activity per barcode summarised to promoters. If the
   sample sheet contains two conditions, the report also tests for differences
   between them.

```
FASTQ ──┬─ PBA:  cutadapt → bowtie2 → extract_pba → qc_pba ──→ pba_clean.tsv ──┐
        │                                                                     │
        └─ counts: count_barcodes → starcode → build_count_matrix ────────────┴─→ mpra_activity.Rmd
```

## Layout

```
config/     config.sh (paths + parameters), design.yaml (read structure),
            samples_pba.tsv, samples_mpra.tsv
scripts/    Python, each runnable standalone with --help
workflow/   SLURM submission scripts, numbered in run order
analysis/   R Markdown report
```

Everything site-specific lives in `config/`. The scripts take no hard-coded
paths, so a new library normally needs only a new `design.yaml` and sample
sheet.

## Setup

```bash
conda env create -f environment.yml                     # main env
conda create -n starcode -c bioconda starcode           # starcode only
```

Then edit `config/config.sh` (paths, conda env names) and
`config/design.yaml` (anchor, barcode length, marker).

## Running

All commands are run from the repository root; logs land in `logs/`.

```bash
# --- stage 1: promoter-barcode association ---
sbatch workflow/00_build_index.sh                 # once per reference FASTA
sbatch --array=1-2 workflow/01_pba.slurm          # one task per row of samples_pba.tsv

# --- stage 2: barcode counts ---
sbatch --array=1-6 workflow/02_count_barcodes.slurm
sbatch workflow/03_collapse_starcode.slurm        # error-correct, d=2
sbatch workflow/04_build_matrix.slurm             # cross-sample matrices

# --- stage 3: analysis ---
Rscript -e 'rmarkdown::render("analysis/mpra_activity.Rmd")'
```

Array sizes must match the sample sheets:
`--array=1-$(($(wc -l < config/samples_mpra.tsv) - 1))`.

## Key outputs

| File | What it is |
| --- | --- |
| `<sample>.pba_clean.tsv` | barcode → promoter, after read-count and dominance filters |
| `<sample>.qc_summary.tsv` | PBA QC: barcodes kept/dropped, promoters detected, empty-vector fraction |
| `library_count_matrix.starcode_d2.tsv` | barcodes × samples, error-corrected |
| `qc_summary_all.tsv` | reads per category per sample |
| `promoter_activity.tsv` | per-promoter activity in each condition |
| `differential_activity.tsv` | between-condition test with adjusted p-values |

## Method notes

- **Orientation.** PBA reads arrive in both orientations, so cutadapt is run
  twice with the anchors swapped and the two passes concatenated. Barcode
  orientation between the PBA table and the count matrix is auto-detected in
  the R report rather than assumed.
- **`bowtie2 --local`.** Empty-vector reference entries are shorter than the
  trimmed reads; local alignment soft-clips the overhang. Real promoter reads
  still align near end-to-end and win on score.
- **Barcode assignment.** A barcode is kept only if it has ≥ `MIN_READS_PER_BC`
  reads and ≥ `DOMINANCE` of them fall on a single promoter.
- **Error correction.** starcode (`-d 2`) collapses 1–2 bp sequencing errors.
  Clustering is per sample, so for exact cross-sample identity snap centroids
  to the PBA barcodes downstream.
- **Normalization.** Counts are converted to CPM over the retained barcodes;
  the pDNA aliquots are averaged into one input reference. A single shared
  pseudocount is used across samples so the floor cancels when conditions are
  subtracted.
- **Filters.** Barcodes need ≥ `pdna_min` reads in each pDNA input; promoters
  need ≥ `min_bc_per_prom` barcodes.
