#!/bin/bash
#SBATCH --job-name=mpra_count
#SBATCH --output=logs/02_count_%A_%a.out
#SBATCH --error=logs/02_count_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# Count barcodes in each MPRA/RMCE sample, one array task per row of
# samples_mpra.tsv.
#
#   sbatch --array=1-$(($(wc -l < config/samples_mpra.tsv) - 1)) workflow/02_count_barcodes.slurm

set -euo pipefail
mkdir -p logs
source "$(dirname "$0")/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

ROW=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))          # +1 skips the header line
SAMPLE=$(cut -f1 "${SAMPLES_MPRA}" | sed -n "${ROW}p")
[[ -n "${SAMPLE}" ]] || { echo "ERROR: no sample on line ${ROW}"; exit 1; }

R1="${MPRA_FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
[[ -s "${R1}" ]] || { echo "ERROR: missing ${R1}"; exit 1; }

python3 "${REPO_DIR}/scripts/count_barcodes.py" \
    --r1 "${R1}" \
    --sample "${SAMPLE}" \
    --outdir "${MPRA_OUT_DIR}/per_sample" \
    --design "${DESIGN}"
