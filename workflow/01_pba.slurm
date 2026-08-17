#!/bin/bash
#SBATCH --job-name=pba
#SBATCH --output=logs/01_pba_%A_%a.out
#SBATCH --error=logs/01_pba_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
# Promoter-barcode association, one array task per row of samples_pba.tsv.
#
#   sbatch --array=1-$(($(wc -l < config/samples_pba.tsv) - 1)) workflow/01_pba.slurm
#
# Steps: cutadapt (both read orientations) -> bowtie2 --local -> extract_pba.py
#        -> qc_pba.py. Final table: <sample>.pba_clean.tsv

set -euo pipefail
mkdir -p logs
source "$(dirname "$0")/../config/config.sh"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

THREADS=${SLURM_CPUS_PER_TASK:-8}
ROW=$(( ${SLURM_ARRAY_TASK_ID:-1} + 1 ))          # +1 skips the header line
read -r SAMPLE FWD_TRIM REV_TRIM < <(sed -n "${ROW}p" "${SAMPLES_PBA}")
[[ -n "${SAMPLE:-}" ]] || { echo "ERROR: no sample on line ${ROW} of ${SAMPLES_PBA}"; exit 1; }

R1="${PBA_FASTQ_DIR}/${SAMPLE}_R1.fastq.gz"
R2="${PBA_FASTQ_DIR}/${SAMPLE}_R2.fastq.gz"
OUT="${PBA_OUT_DIR}/${SAMPLE}"
WORK="${OUT}/work"
mkdir -p "${OUT}/logs" "${OUT}/bam" "${OUT}/pba" "${WORK}"

for f in "${R1}" "${R2}" "${PBA_INDEX}.1.bt2"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing ${f}"; exit 1; }
done

echo "[$(date)] sample=${SAMPLE}  fwd=${FWD_TRIM}  rev=${REV_TRIM}"

# ---------------------------------------------------------------- 1. cutadapt
# Reads arrive in either orientation, so trim twice with the anchors swapped
# and concatenate. --discard-untrimmed keeps only correctly structured pairs.
MIN_BARCODE_LEN=$(( BARCODE_LEN + UNIV_ANCHOR_LEN ))

cutadapt -j "${THREADS}" \
    -g "^${FWD_TRIM}" -G "^${REV_TRIM}" \
    -a "${TRUSEQ_R1}" -A "${TRUSEQ_R2}" \
    -e "${CUTADAPT_ERR}" --no-indels --discard-untrimmed --pair-filter=any \
    --nextseq-trim=20 -q 20 -m "${MIN_PROMOTER_LEN}:${MIN_BARCODE_LEN}" \
    -o "${WORK}/A_promoter.fq.gz" -p "${WORK}/A_barcode.fq.gz" \
    "${R1}" "${R2}" > "${OUT}/logs/cutadapt_A.log"

cutadapt -j "${THREADS}" \
    -g "^${REV_TRIM}" -G "^${FWD_TRIM}" \
    -a "${TRUSEQ_R1}" -A "${TRUSEQ_R2}" \
    -e "${CUTADAPT_ERR}" --no-indels --discard-untrimmed --pair-filter=any \
    --nextseq-trim=20 -q 20 -m "${MIN_BARCODE_LEN}:${MIN_PROMOTER_LEN}" \
    -o "${WORK}/B_barcode.fq.gz" -p "${WORK}/B_promoter.fq.gz" \
    "${R1}" "${R2}" > "${OUT}/logs/cutadapt_B.log"

cat "${WORK}/A_promoter.fq.gz" "${WORK}/B_promoter.fq.gz" > "${WORK}/promoter.fq.gz"
cat "${WORK}/A_barcode.fq.gz"  "${WORK}/B_barcode.fq.gz"  > "${WORK}/barcode.fq.gz"

{
    echo -e "raw_pairs\t$(zcat "${R1}" | awk 'END{print NR/4}')"
    echo -e "orientation_A\t$(zcat "${WORK}/A_promoter.fq.gz" | awk 'END{print NR/4}')"
    echo -e "orientation_B\t$(zcat "${WORK}/B_promoter.fq.gz" | awk 'END{print NR/4}')"
} | tee "${OUT}/logs/read_counts.tsv"

# ---------------------------------------------------------------- 2. bowtie2
# --local so that short Empty_* references still align despite read overhang;
# full-length promoter alignments win on score.
bowtie2 -x "${PBA_INDEX}" -U "${WORK}/promoter.fq.gz" \
    -p "${THREADS}" --local --no-unal 2> "${OUT}/logs/bowtie2.log" \
    | samtools sort -@ "${THREADS}" -o "${OUT}/bam/${SAMPLE}.bam" -
samtools index "${OUT}/bam/${SAMPLE}.bam"
samtools flagstat "${OUT}/bam/${SAMPLE}.bam" > "${OUT}/logs/flagstat.txt"

# ---------------------------------------------------------------- 3. join
python3 "${REPO_DIR}/scripts/extract_pba.py" \
    --bam "${OUT}/bam/${SAMPLE}.bam" \
    --barcode-fq "${WORK}/barcode.fq.gz" \
    --bc-len "${BARCODE_LEN}" --mapq-min "${MAPQ_MIN}" --rc-barcode \
    --universal "${UNIVERSAL}" \
    --univ-anchor-len "${UNIV_ANCHOR_LEN}" --univ-anchor-mm "${UNIV_ANCHOR_MM}" \
    --out-pba    "${OUT}/pba/${SAMPLE}.pba_raw.tsv" \
    --out-counts "${OUT}/pba/${SAMPLE}.pba_counts.tsv" \
    --out-stats  "${OUT}/pba/${SAMPLE}.pba_stats.tsv"

# ---------------------------------------------------------------- 4. QC
python3 "${REPO_DIR}/scripts/qc_pba.py" \
    --counts "${OUT}/pba/${SAMPLE}.pba_counts.tsv" \
    --min-reads "${MIN_READS_PER_BC}" --dominance "${DOMINANCE}" \
    --n-promoters-expected "${N_PROMOTERS_EXPECTED}" \
    --empty-prefix "${EMPTY_PREFIX}" \
    --out-pba            "${OUT}/pba/${SAMPLE}.pba_clean.tsv" \
    --out-ambiguous      "${OUT}/pba/${SAMPLE}.ambiguous_barcodes.tsv" \
    --out-promoter-stats "${OUT}/pba/${SAMPLE}.promoter_stats.tsv" \
    --out-summary        "${OUT}/pba/${SAMPLE}.qc_summary.tsv"
cat "${OUT}/pba/${SAMPLE}.qc_summary.tsv"

rm -rf "${WORK}"
echo "[$(date)] DONE ${SAMPLE} -> ${OUT}/pba/${SAMPLE}.pba_clean.tsv"
