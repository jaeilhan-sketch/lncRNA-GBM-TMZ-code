#!/usr/bin/env bash
# ============================================================================
# Phase 1: Download controlled-access FASTQ files from GDC
#
# Prerequisites:
#   - dbGaP approved access (eRA Commons)
#   - GDC user token (download from GDC portal after login)
#   - gdc-client installed (pip install gdc-client)
#
# Usage:
#   1. Place your GDC token at: data/raw/gdc-user-token.txt
#   2. Run: bash scripts/01_data_acquisition/03_download_fastq.sh
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
RAW_DIR="${PROJECT_DIR}/data/raw"
CLINICAL_DIR="${PROJECT_DIR}/data/clinical"
TOKEN_FILE="${RAW_DIR}/gdc-user-token.txt"

echo "============================================"
echo " Download TCGA-GBM FASTQ files"
echo "============================================"

# ── Check token ──
if [[ ! -f "${TOKEN_FILE}" ]]; then
    echo "ERROR: GDC token not found at ${TOKEN_FILE}"
    echo ""
    echo "Steps to get your token:"
    echo "  1. Go to https://portal.gdc.cancer.gov/"
    echo "  2. Login with eRA Commons credentials"
    echo "  3. Click your username → 'Download Token'"
    echo "  4. Save to: ${TOKEN_FILE}"
    exit 1
fi

# ── Generate manifest for selected samples ──
echo "[1/3] Generating download manifest from selected samples..."

SAMPLE_INFO="${CLINICAL_DIR}/sample_info_final.tsv"
if [[ ! -f "${SAMPLE_INFO}" ]]; then
    echo "ERROR: sample_info_final.tsv not found."
    echo "Run 01_query_gdc_metadata.py and 02_download_molecular_annotations.py first."
    exit 1
fi

# Extract file_ids from sample_info
MANIFEST="${RAW_DIR}/gdc_download_manifest.txt"
echo -e "id\tfilename\tmd5\tsize\tstate" > "${MANIFEST}"
tail -n +2 "${SAMPLE_INFO}" | cut -f1 | while read -r file_id; do
    echo -e "${file_id}\t\t\t\t"
done >> "${MANIFEST}"

N_FILES=$(( $(wc -l < "${MANIFEST}") - 1 ))
echo "  Files to download: ${N_FILES}"

# ── Download with gdc-client ──
echo "[2/3] Downloading FASTQ files via gdc-client..."
echo "  This may take several hours depending on bandwidth."
echo ""

cd "${RAW_DIR}"
gdc-client download \
    -m "${MANIFEST}" \
    -t "${TOKEN_FILE}" \
    -n 4 \
    --retry-amount 3 \
    --log-file "${RAW_DIR}/gdc_download.log"

# ── Organize downloaded files ──
echo "[3/3] Organizing downloaded files..."

# gdc-client creates subdirectories per file_id
# Move FASTQ files to flat structure with sample IDs
find "${RAW_DIR}" -name "*.fastq.gz" -type f | while read -r fq; do
    fname=$(basename "${fq}")
    if [[ ! -f "${RAW_DIR}/${fname}" ]]; then
        ln -s "${fq}" "${RAW_DIR}/${fname}"
    fi
done

echo ""
echo "============================================"
echo " Download complete!"
echo " FASTQ location: ${RAW_DIR}/"
echo " Log: ${RAW_DIR}/gdc_download.log"
echo " Next: Phase 2 - Preprocessing"
echo "============================================"
