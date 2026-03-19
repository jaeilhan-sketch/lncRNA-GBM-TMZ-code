#!/usr/bin/env bash
# ============================================================================
# Phase 2.1: Raw FASTQ Quality Control with FastQC + MultiQC
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
RAW_DIR="${PROJECT_DIR}/data/raw"
FASTQC_DIR="${PROJECT_DIR}/data/processed/fastqc/raw"
REPORT_DIR="${PROJECT_DIR}/results/reports"
THREADS=8

mkdir -p "${FASTQC_DIR}" "${REPORT_DIR}"

echo "============================================"
echo " Phase 2.1: Raw Data QC (FastQC)"
echo "============================================"

# ── FastQC on all raw FASTQ files ──
echo "[1/2] Running FastQC..."
find "${RAW_DIR}" -name "*.fastq.gz" -type f | \
    parallel -j "${THREADS}" fastqc -o "${FASTQC_DIR}" --threads 2 {}

echo "  FastQC complete: ${FASTQC_DIR}/"

# ── MultiQC aggregate report ──
echo "[2/2] Generating MultiQC report..."
multiqc "${FASTQC_DIR}" \
    -o "${REPORT_DIR}" \
    -n "multiqc_raw_fastqc" \
    --force

echo ""
echo "  Report: ${REPORT_DIR}/multiqc_raw_fastqc.html"
echo "  Review the report before proceeding to trimming."
echo "============================================"
