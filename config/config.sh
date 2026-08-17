#!/bin/bash
# config.sh -- the ONLY file you should need to edit per project.
# Sourced by every script in workflow/.

# ---------------------------------------------------------------- paths
PROJECT_DIR="/path/to/project"
REPO_DIR="${PROJECT_DIR}/mpra-pipeline"           # this repository

# PBA (promoter-barcode association) run
PBA_FASTQ_DIR="${PROJECT_DIR}/PBA/fastq"
PBA_OUT_DIR="${PROJECT_DIR}/PBA/results"
PBA_REF_FASTA="${PROJECT_DIR}/PBA/reference/library_reference.fasta"
PBA_INDEX="${PROJECT_DIR}/PBA/reference/library"  # bowtie2 index prefix

# Barcode-counting run (RNA and plasmid samples)
MPRA_FASTQ_DIR="${PROJECT_DIR}/MPRA/fastq"
MPRA_OUT_DIR="${PROJECT_DIR}/MPRA/results"

# ---------------------------------------------------------------- env
CONDA_ENV="mpra-pipeline"
STARCODE_ENV="starcode"           # separate env that provides `starcode`

# ---------------------------------------------------------------- design
DESIGN="${REPO_DIR}/config/design.yaml"
SAMPLES_PBA="${REPO_DIR}/config/samples_pba.tsv"
SAMPLES_MPRA="${REPO_DIR}/config/samples_mpra.tsv"

# ---------------------------------------------------------------- PBA params
UNIVERSAL="ATCAGCCCTGGGAAGGTGCATGTGCCATAGGGATAACAGGGTAAT"
BARCODE_LEN=12
UNIV_ANCHOR_LEN=10       # bp of revcomp(universal) required after the barcode
UNIV_ANCHOR_MM=1
MIN_PROMOTER_LEN=50
CUTADAPT_ERR=0.15
MAPQ_MIN=20
MIN_READS_PER_BC=3       # min reads to call a barcode -> promoter link
DOMINANCE=0.90           # min fraction of reads on the winning promoter
N_PROMOTERS_EXPECTED=1000
EMPTY_PREFIX="Empty"     # reference names for empty-vector controls

TRUSEQ_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRUSEQ_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ---------------------------------------------------------------- count params
STARCODE_DIST=2          # Levenshtein distance for barcode error correction
