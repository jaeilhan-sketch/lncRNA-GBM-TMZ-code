#!/usr/bin/env bash
# ============================================================================
# Phase 2.2: Adapter Trimming and Quality Filtering (Trim Galore)
#
# Processes paired-end FASTQ files.
# Naming convention: {sample}_R1.fastq.gz, {sample}_R2.fastq.gz
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
RAW_DIR="${PROJECT_DIR}/data/raw"
TRIMMED_DIR="${PROJECT_DIR}/data/processed/trimmed"
FASTQC_TRIM_DIR="${PROJECT_DIR}/data/processed/fastqc/trimmed"
REPORT_DIR="${PROJECT_DIR}/results/reports"
THREADS=4
PARALLEL_JOBS=4

mkdir -p "${TRIMMED_DIR}" "${FASTQC_TRIM_DIR}"

echo "============================================"
echo " Phase 2.2: Adapter Trimming (Trim Galore)"
echo "============================================"

# ── Collect sample list ──
# Expects paired-end files: {sample}_R1.fastq.gz / {sample}_R2.fastq.gz
# Or: {sample}_1.fastq.gz / {sample}_2.fastq.gz
SAMPLE_LIST=$(find "${RAW_DIR}" -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" | \
    sed 's/_R1\.fastq\.gz//' | sed 's/_1\.fastq\.gz//' | sort -u)

N_SAMPLES=$(echo "${SAMPLE_LIST}" | wc -l)
echo "  Found ${N_SAMPLES} paired-end samples."

# ── Run Trim Galore ──
run_trim_galore() {
    local sample_prefix="$1"
    local sample_name=$(basename "${sample_prefix}")

    # Detect naming pattern
    if [[ -f "${sample_prefix}_R1.fastq.gz" ]]; then
        R1="${sample_prefix}_R1.fastq.gz"
        R2="${sample_prefix}_R2.fastq.gz"
    elif [[ -f "${sample_prefix}_1.fastq.gz" ]]; then
        R1="${sample_prefix}_1.fastq.gz"
        R2="${sample_prefix}_2.fastq.gz"
    else
        echo "  WARN: Cannot find paired files for ${sample_name}. Skipping."
        return
    fi

    echo "  Trimming: ${sample_name}"
    trim_galore \
        --paired \
        --quality 20 \
        --length 36 \
        --cores "${THREADS}" \
        --fastqc \
        --fastqc_args "--outdir ${FASTQC_TRIM_DIR}" \
        -o "${TRIMMED_DIR}" \
        "${R1}" "${R2}" \
        2>&1 | tee "${TRIMMED_DIR}/${sample_name}_trim.log"
}

export -f run_trim_galore
export TRIMMED_DIR FASTQC_TRIM_DIR THREADS

echo "[1/2] Running Trim Galore on all samples..."
echo "${SAMPLE_LIST}" | parallel -j "${PARALLEL_JOBS}" run_trim_galore {}

# ── MultiQC for trimmed data ──
echo "[2/2] Generating post-trimming MultiQC report..."
multiqc "${FASTQC_TRIM_DIR}" "${TRIMMED_DIR}" \
    -o "${REPORT_DIR}" \
    -n "multiqc_trimmed" \
    --force

echo ""
echo "  Trimmed files: ${TRIMMED_DIR}/"
echo "  Report: ${REPORT_DIR}/multiqc_trimmed.html"
echo "============================================"
