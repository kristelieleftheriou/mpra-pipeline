#!/bin/bash
#SBATCH --job-name=mpra_matrix
#SBATCH --output=logs/04_matrix_%j.out
#SBATCH --error=logs/04_matrix_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#
# Build cross-sample count matrices (raw and starcode-collapsed) plus the
# combined QC table. Column order follows samples_mpra.tsv.
#
#   sbatch workflow/04_build_matrix.slurm

set -euo pipefail
mkdir -p logs
source "$(dirname "$0")/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

mapfile -t SAMPLES < <(tail -n +2 "${SAMPLES_MPRA}" | cut -f1)

python3 "${REPO_DIR}/scripts/build_count_matrix.py" \
    --indir  "${MPRA_OUT_DIR}/per_sample" \
    --outdir "${MPRA_OUT_DIR}/aggregated" \
    --samples "${SAMPLES[@]}"

python3 "${REPO_DIR}/scripts/build_count_matrix.py" \
    --indir  "${MPRA_OUT_DIR}/per_sample" \
    --outdir "${MPRA_OUT_DIR}/aggregated" \
    --suffix "starcode_d${STARCODE_DIST}" \
    --samples "${SAMPLES[@]}"
