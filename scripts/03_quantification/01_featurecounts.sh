#!/usr/bin/env bash
# ============================================================================
# Phase 3.1: Gene-level read counting with featureCounts
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
ALIGNED_DIR="${PROJECT_DIR}/data/processed/aligned"
COUNTS_DIR="${PROJECT_DIR}/data/processed/counts"
GTF="${PROJECT_DIR}/data/reference/annotation/gencode.v44.annotation.gtf"
THREADS=8

mkdir -p "${COUNTS_DIR}"

echo "============================================"
echo " Phase 3.1: featureCounts"
echo "============================================"

# ── Collect BAM files (exclude QC failures) ──
QC_FAIL="${ALIGNED_DIR}/qc/qc_failed_samples.txt"
BAM_LIST=()

for bam in "${ALIGNED_DIR}"/*Aligned.sortedByCoord.out.bam; do
    sample=$(basename "${bam}" | sed 's/_Aligned.sortedByCoord.out.bam//')

    # Skip QC failures
    if [[ -f "${QC_FAIL}" ]] && grep -q "${sample}" "${QC_FAIL}"; then
        echo "  SKIP (QC fail): ${sample}"
        continue
    fi

    BAM_LIST+=("${bam}")
done

echo "  BAM files to count: ${#BAM_LIST[@]}"

# ── Run featureCounts ──
# -s 2: reversely stranded (common for Illumina TruSeq stranded)
# Verify strandedness from infer_experiment.py results before running
echo ""
echo "  Running featureCounts..."
featureCounts \
    -a "${GTF}" \
    -o "${COUNTS_DIR}/raw_counts_all.txt" \
    -T "${THREADS}" \
    -p --countReadPairs \
    -s 2 \
    -t exon \
    -g gene_id \
    --extraAttributes gene_name,gene_type \
    -B -C \
    "${BAM_LIST[@]}" \
    2>&1 | tee "${COUNTS_DIR}/featurecounts.log"

echo ""
echo "  Raw counts: ${COUNTS_DIR}/raw_counts_all.txt"
echo "  Summary:    ${COUNTS_DIR}/raw_counts_all.txt.summary"
echo "============================================"
