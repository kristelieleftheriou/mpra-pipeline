#!/bin/bash
#SBATCH --job-name=bt2_build
#SBATCH --output=logs/00_build_index_%j.out
#SBATCH --error=logs/00_build_index_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#
# Build the bowtie2 index of the promoter reference. Re-run whenever the
# reference FASTA changes (e.g. after adding Empty_* control entries).
#
#   sbatch workflow/00_build_index.sh

set -euo pipefail
mkdir -p logs
source "$(dirname "$0")/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

bowtie2-build --threads "${SLURM_CPUS_PER_TASK:-4}" -f "${PBA_REF_FASTA}" "${PBA_INDEX}"

echo "Reference entries:"
bowtie2-inspect -n "${PBA_INDEX}" | head
echo "Total: $(bowtie2-inspect -n "${PBA_INDEX}" | wc -l)"
