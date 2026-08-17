#!/bin/bash
#SBATCH --job-name=mpra_starcode
#SBATCH --output=logs/03_starcode_%j.out
#SBATCH --error=logs/03_starcode_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#
# Error-correct barcodes across all samples at once. Uses the starcode env,
# not the main analysis env.
#
#   sbatch workflow/03_collapse_starcode.slurm

set -euo pipefail
mkdir -p logs
source "$(dirname "$0")/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${STARCODE_ENV}"

python3 "${REPO_DIR}/scripts/collapse_barcodes.py" \
    --indir   "${MPRA_OUT_DIR}/per_sample" \
    --dist    "${STARCODE_DIST}" \
    --threads "${SLURM_CPUS_PER_TASK:-4}"
