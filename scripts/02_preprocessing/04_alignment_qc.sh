#!/usr/bin/env bash
# ============================================================================
# Phase 2.4: Alignment Quality Control
#
# Uses RSeQC and samtools for post-alignment QC metrics.
# Generates summary statistics and identifies samples to exclude.
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
ALIGNED_DIR="${PROJECT_DIR}/data/processed/aligned"
REPORT_DIR="${PROJECT_DIR}/results/reports"
QC_DIR="${ALIGNED_DIR}/qc"
BED_FILE="${PROJECT_DIR}/data/reference/annotation/gencode.v44.annotation.bed"

mkdir -p "${QC_DIR}"

echo "============================================"
echo " Phase 2.4: Alignment QC (RSeQC + samtools)"
echo "============================================"

# ── Thresholds ──
MIN_UNIQUELY_MAPPED_PCT=70
MAX_RRNA_PCT=5

# ── Per-sample QC ──
STATS_FILE="${QC_DIR}/alignment_stats_summary.tsv"
echo -e "sample\ttotal_reads\tuniquely_mapped\tuniquely_mapped_pct\tmulti_mapped\tunmapped" \
    > "${STATS_FILE}"

for bam in "${ALIGNED_DIR}"/*Aligned.sortedByCoord.out.bam; do
    sample=$(basename "${bam}" | sed 's/_Aligned.sortedByCoord.out.bam//')
    echo "  Processing: ${sample}"

    # ── samtools flagstat ──
    samtools flagstat "${bam}" > "${QC_DIR}/${sample}_flagstat.txt"

    # ── Parse STAR Log.final.out for detailed stats ──
    STAR_LOG="${ALIGNED_DIR}/${sample}_Log.final.out"
    if [[ -f "${STAR_LOG}" ]]; then
        total=$(grep "Number of input reads" "${STAR_LOG}" | awk '{print $NF}')
        unique=$(grep "Uniquely mapped reads number" "${STAR_LOG}" | awk '{print $NF}')
        unique_pct=$(grep "Uniquely mapped reads %" "${STAR_LOG}" | awk '{print $NF}' | tr -d '%')
        multi=$(grep "Number of reads mapped to multiple loci" "${STAR_LOG}" | awk '{print $NF}')
        unmapped=$(grep "Number of reads unmapped: too many mismatches" "${STAR_LOG}" | awk '{print $NF}')

        echo -e "${sample}\t${total}\t${unique}\t${unique_pct}\t${multi}\t${unmapped}" \
            >> "${STATS_FILE}"
    fi

    # ── RSeQC: infer_experiment (strandedness) ──
    if [[ -f "${BED_FILE}" ]] && [[ ! -f "${QC_DIR}/${sample}_infer_experiment.txt" ]]; then
        infer_experiment.py -i "${bam}" -r "${BED_FILE}" \
            > "${QC_DIR}/${sample}_infer_experiment.txt" 2>&1
    fi

    # ── RSeQC: read_distribution ──
    if [[ -f "${BED_FILE}" ]] && [[ ! -f "${QC_DIR}/${sample}_read_distribution.txt" ]]; then
        read_distribution.py -i "${bam}" -r "${BED_FILE}" \
            > "${QC_DIR}/${sample}_read_distribution.txt" 2>&1
    fi

done

# ── Identify QC failures ──
echo ""
echo "  Checking for QC failures..."
FAIL_FILE="${QC_DIR}/qc_failed_samples.txt"
> "${FAIL_FILE}"

tail -n +2 "${STATS_FILE}" | while IFS=$'\t' read -r sample total unique unique_pct multi unmapped; do
    pct_int=$(echo "${unique_pct}" | cut -d. -f1)
    if [[ "${pct_int}" -lt "${MIN_UNIQUELY_MAPPED_PCT}" ]]; then
        echo "${sample}  uniquely_mapped_pct=${unique_pct}%" >> "${FAIL_FILE}"
    fi
done

if [[ -s "${FAIL_FILE}" ]]; then
    echo "  WARNING: QC failures detected:"
    cat "${FAIL_FILE}"
else
    echo "  All samples passed alignment QC."
fi

# ── MultiQC for alignment ──
echo ""
echo "  Generating alignment MultiQC report..."
multiqc "${ALIGNED_DIR}" "${QC_DIR}" \
    -o "${REPORT_DIR}" \
    -n "multiqc_alignment" \
    --force

echo ""
echo "============================================"
echo " Alignment QC complete."
echo " Stats:  ${STATS_FILE}"
echo " Report: ${REPORT_DIR}/multiqc_alignment.html"
echo " Fails:  ${FAIL_FILE}"
echo "============================================"
