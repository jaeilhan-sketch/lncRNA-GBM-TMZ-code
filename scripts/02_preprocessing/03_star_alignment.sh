#!/usr/bin/env bash
# ============================================================================
# Phase 2.3: STAR Alignment
#
# Aligns trimmed reads to GRCh38 using STAR 2-pass mode.
# Outputs sorted BAM files.
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
TRIMMED_DIR="${PROJECT_DIR}/data/processed/trimmed"
ALIGNED_DIR="${PROJECT_DIR}/data/processed/aligned"
STAR_INDEX="${PROJECT_DIR}/data/reference/index/star"
GTF="${PROJECT_DIR}/data/reference/annotation/gencode.v44.annotation.gtf"
THREADS=16

mkdir -p "${ALIGNED_DIR}"

echo "============================================"
echo " Phase 2.3: STAR Alignment"
echo "============================================"

# ── Verify STAR index ──
if [[ ! -f "${STAR_INDEX}/SA" ]]; then
    echo "ERROR: STAR index not found at ${STAR_INDEX}/"
    echo "Run scripts/00_setup/setup_environment.sh first."
    exit 1
fi

# ── Collect trimmed sample list ──
# Trim Galore outputs: {sample}_R1_val_1.fq.gz, {sample}_R2_val_2.fq.gz
SAMPLE_LIST=$(find "${TRIMMED_DIR}" -name "*_val_1.fq.gz" | \
    sed 's/_R1_val_1\.fq\.gz//' | sed 's/_1_val_1\.fq\.gz//' | sort -u)

N_SAMPLES=$(echo "${SAMPLE_LIST}" | wc -l)
echo "  Found ${N_SAMPLES} trimmed samples to align."

# ── STAR alignment per sample ──
for sample_prefix in ${SAMPLE_LIST}; do
    sample_name=$(basename "${sample_prefix}")

    # Detect trimmed file naming
    if [[ -f "${sample_prefix}_R1_val_1.fq.gz" ]]; then
        R1="${sample_prefix}_R1_val_1.fq.gz"
        R2="${sample_prefix}_R2_val_2.fq.gz"
    elif [[ -f "${sample_prefix}_1_val_1.fq.gz" ]]; then
        R1="${sample_prefix}_1_val_1.fq.gz"
        R2="${sample_prefix}_2_val_2.fq.gz"
    else
        echo "  WARN: Trimmed files not found for ${sample_name}. Skipping."
        continue
    fi

    OUTPREFIX="${ALIGNED_DIR}/${sample_name}_"

    # Skip if already aligned
    if [[ -f "${OUTPREFIX}Aligned.sortedByCoord.out.bam" ]]; then
        echo "  SKIP (exists): ${sample_name}"
        continue
    fi

    echo "  Aligning: ${sample_name}"

    STAR \
        --runMode alignReads \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${R1}" "${R2}" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outSAMstrandField intronMotif \
        --outSAMattributes NH HI AS NM MD \
        --runThreadN "${THREADS}" \
        --outFileNamePrefix "${OUTPREFIX}" \
        2>&1 | tee "${ALIGNED_DIR}/${sample_name}_star.log"

    # Index BAM
    samtools index -@ 4 "${OUTPREFIX}Aligned.sortedByCoord.out.bam"

done

echo ""
echo "============================================"
echo " STAR alignment complete."
echo " BAM files: ${ALIGNED_DIR}/"
echo " Next: Phase 2.4 - Alignment QC"
echo "============================================"
