#!/usr/bin/env bash
# ============================================================================
# Disk-Efficient Pipeline Runner
#
# 전체 파이프라인을 샘플 단위로 실행하며, 각 단계 완료 후 중간 파일을 삭제.
# 디스크 사용량을 최소화하면서도 재현성은 스크립트로 보장.
#
# 디스크 사용 패턴:
#   Raw FASTQ (1 sample)   : ~20 GB    ← 처리 후 유지 (controlled access 재다운로드 어려움)
#   Trimmed FASTQ (1 sample): ~16 GB   ← alignment 후 삭제
#   BAM (1 sample)          : ~15 GB   ← featureCounts 후 삭제
#   동시 최대 사용: ~51 GB/sample + Raw FASTQ 전체
#
# 최종 보존 파일:
#   - Raw FASTQ (전체)           ← 재현성의 핵심 (약 600-1200 GB)
#   - featureCounts 결과          ← < 1 GB
#   - QC reports (FastQC/MultiQC) ← < 1 GB
#   - STAR Log.final.out          ← alignment 통계 (< 1 MB/sample)
#   - R 분석 결과 전체            ← < 5 GB
#
# Usage:
#   bash scripts/run_pipeline_batched.sh [--keep-bam] [--keep-trimmed]
#   bash scripts/run_pipeline_batched.sh --dry-run
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
RAW_DIR="${PROJECT_DIR}/data/raw"
TRIMMED_DIR="${PROJECT_DIR}/data/processed/trimmed"
ALIGNED_DIR="${PROJECT_DIR}/data/processed/aligned"
COUNTS_DIR="${PROJECT_DIR}/data/processed/counts"
FASTQC_DIR="${PROJECT_DIR}/data/processed/fastqc"
REPORT_DIR="${PROJECT_DIR}/results/reports"
STAR_INDEX="${PROJECT_DIR}/data/reference/index/star"
GTF="${PROJECT_DIR}/data/reference/annotation/gencode.v44.annotation.gtf"

THREADS=16
TRIM_THREADS=4

# ── Parse options ──
KEEP_BAM=false
KEEP_TRIMMED=false
DRY_RUN=false

for arg in "$@"; do
    case "$arg" in
        --keep-bam)     KEEP_BAM=true ;;
        --keep-trimmed) KEEP_TRIMMED=true ;;
        --dry-run)      DRY_RUN=true ;;
    esac
done

mkdir -p "${TRIMMED_DIR}" "${ALIGNED_DIR}" "${ALIGNED_DIR}/qc" \
         "${COUNTS_DIR}" "${FASTQC_DIR}/raw" "${FASTQC_DIR}/trimmed" \
         "${REPORT_DIR}" "${ALIGNED_DIR}/logs" "${ALIGNED_DIR}/star_gene_counts"

# ============================================================================
# Utility functions
# ============================================================================
disk_usage() {
    local used avail
    read -r used avail <<< $(df --output=used,avail /media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68 | tail -1)
    local used_gb=$(( used / 1048576 ))
    local avail_gb=$(( avail / 1048576 ))
    echo "  [DISK] Used: ${used_gb} GB | Available: ${avail_gb} GB"
}

log_step() {
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  $1"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    disk_usage
}

# ============================================================================
# Discover samples
# ============================================================================
SAMPLE_LIST=$(find "${RAW_DIR}" -maxdepth 2 -name "*_R1.fastq.gz" -o -name "*_1.fastq.gz" 2>/dev/null | \
    sed 's/_R1\.fastq\.gz//' | sed 's/_1\.fastq\.gz//' | sort -u)

N_SAMPLES=$(echo "${SAMPLE_LIST}" | grep -c '.' || true)

echo "╔══════════════════════════════════════════════════╗"
echo "║  lncRNA-GBM-TMZ: Disk-Efficient Pipeline        ║"
echo "╠══════════════════════════════════════════════════╣"
echo "║  Samples found : ${N_SAMPLES}                            "
echo "║  Keep BAM      : ${KEEP_BAM}                         "
echo "║  Keep trimmed  : ${KEEP_TRIMMED}                      "
echo "║  Dry run       : ${DRY_RUN}                           "
echo "╚══════════════════════════════════════════════════╝"
disk_usage

if [[ "${DRY_RUN}" == true ]]; then
    echo ""
    echo "DRY RUN — listing samples that would be processed:"
    echo "${SAMPLE_LIST}" | while read -r sp; do echo "  $(basename "$sp")"; done
    exit 0
fi

if [[ "${N_SAMPLES}" -eq 0 ]]; then
    echo "ERROR: No FASTQ files found in ${RAW_DIR}"
    exit 1
fi

# ============================================================================
# Phase 2-3: Per-sample processing loop
# ============================================================================
PROCESSED=0
FAILED=0
FAILED_SAMPLES=""

# Collect BAM paths for featureCounts (if keeping BAMs) or per-sample gene counts
STAR_GENE_COUNTS_DIR="${ALIGNED_DIR}/star_gene_counts"

for sample_prefix in ${SAMPLE_LIST}; do
    sample_name=$(basename "${sample_prefix}")
    PROCESSED=$((PROCESSED + 1))

    log_step "[${PROCESSED}/${N_SAMPLES}] Processing: ${sample_name}"

    # Detect R1/R2 naming
    if [[ -f "${sample_prefix}_R1.fastq.gz" ]]; then
        R1="${sample_prefix}_R1.fastq.gz"
        R2="${sample_prefix}_R2.fastq.gz"
    elif [[ -f "${sample_prefix}_1.fastq.gz" ]]; then
        R1="${sample_prefix}_1.fastq.gz"
        R2="${sample_prefix}_2.fastq.gz"
    else
        echo "  WARN: Paired files not found. Skipping."
        FAILED=$((FAILED + 1))
        FAILED_SAMPLES="${FAILED_SAMPLES} ${sample_name}"
        continue
    fi

    # Skip if STAR gene counts already exist (resume support)
    if [[ -f "${STAR_GENE_COUNTS_DIR}/${sample_name}_ReadsPerGene.out.tab" ]]; then
        echo "  SKIP (already processed): ${sample_name}"
        continue
    fi

    # ────────────────────────────────────────────────
    # Step 1: FastQC (raw)
    # ────────────────────────────────────────────────
    echo "  [1/5] FastQC (raw)..."
    fastqc -o "${FASTQC_DIR}/raw/" --threads 2 "${R1}" "${R2}" 2>/dev/null || true

    # ────────────────────────────────────────────────
    # Step 2: Trim Galore
    # ────────────────────────────────────────────────
    echo "  [2/5] Trim Galore..."
    trim_galore \
        --paired --quality 20 --length 36 \
        --cores "${TRIM_THREADS}" \
        --fastqc --fastqc_args "--outdir ${FASTQC_DIR}/trimmed" \
        -o "${TRIMMED_DIR}" \
        "${R1}" "${R2}" \
        > "${TRIMMED_DIR}/${sample_name}_trim.log" 2>&1

    # Detect trimmed file names
    if [[ -f "${TRIMMED_DIR}/$(basename "${R1}" .fastq.gz)_val_1.fq.gz" ]]; then
        TRIM_R1="${TRIMMED_DIR}/$(basename "${R1}" .fastq.gz)_val_1.fq.gz"
        TRIM_R2="${TRIMMED_DIR}/$(basename "${R2}" .fastq.gz)_val_2.fq.gz"
    else
        echo "  ERROR: Trimmed files not found for ${sample_name}"
        FAILED=$((FAILED + 1))
        FAILED_SAMPLES="${FAILED_SAMPLES} ${sample_name}"
        continue
    fi

    # ────────────────────────────────────────────────
    # Step 3: STAR Alignment
    # ────────────────────────────────────────────────
    echo "  [3/5] STAR alignment..."
    OUTPREFIX="${ALIGNED_DIR}/${sample_name}_"

    STAR \
        --runMode alignReads \
        --genomeDir "${STAR_INDEX}" \
        --readFilesIn "${TRIM_R1}" "${TRIM_R2}" \
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
        > "${ALIGNED_DIR}/logs/${sample_name}_star.log" 2>&1

    BAM_FILE="${OUTPREFIX}Aligned.sortedByCoord.out.bam"

    if [[ ! -f "${BAM_FILE}" ]]; then
        echo "  ERROR: STAR alignment failed for ${sample_name}"
        FAILED=$((FAILED + 1))
        FAILED_SAMPLES="${FAILED_SAMPLES} ${sample_name}"
        # Cleanup partial STAR output
        rm -f "${OUTPREFIX}"*.out* "${OUTPREFIX}"*._STARtmp 2>/dev/null || true
        rm -rf "${OUTPREFIX}_STARtmp" 2>/dev/null || true
        continue
    fi

    # ── DELETE trimmed FASTQ (no longer needed) ──
    if [[ "${KEEP_TRIMMED}" == false ]]; then
        echo "  [CLEANUP] Removing trimmed FASTQ..."
        rm -f "${TRIM_R1}" "${TRIM_R2}"
        rm -f "${TRIMMED_DIR}/${sample_name}"*_trimming_report.txt 2>/dev/null || true
    fi

    # ────────────────────────────────────────────────
    # Step 4: Quick alignment QC (from STAR Log)
    # ────────────────────────────────────────────────
    echo "  [4/5] Alignment QC..."
    STAR_LOG="${OUTPREFIX}Log.final.out"

    # Save STAR log permanently (tiny file, ~2 KB)
    cp "${STAR_LOG}" "${ALIGNED_DIR}/logs/${sample_name}_Log.final.out" 2>/dev/null || true

    # Save STAR gene counts permanently (tiny file, ~1 MB)
    cp "${OUTPREFIX}ReadsPerGene.out.tab" \
       "${STAR_GENE_COUNTS_DIR}/${sample_name}_ReadsPerGene.out.tab" 2>/dev/null || true

    # Check uniquely mapped rate
    if [[ -f "${STAR_LOG}" ]]; then
        UNIQUE_PCT=$(grep "Uniquely mapped reads %" "${STAR_LOG}" | awk '{print $NF}' | tr -d '%')
        echo "  Uniquely mapped: ${UNIQUE_PCT}%"

        # Flag if below threshold
        PCT_INT=$(echo "${UNIQUE_PCT}" | cut -d. -f1)
        if [[ "${PCT_INT}" -lt 70 ]]; then
            echo "  WARNING: Low mapping rate (${UNIQUE_PCT}%) — flagged for review"
            echo "${sample_name}  uniquely_mapped=${UNIQUE_PCT}%" \
                >> "${ALIGNED_DIR}/qc/qc_flagged_samples.txt"
        fi
    fi

    # ────────────────────────────────────────────────
    # Step 5: Per-sample featureCounts → DELETE BAM
    # ────────────────────────────────────────────────
    echo "  [5/5] featureCounts (per-sample)..."

    featureCounts \
        -a "${GTF}" \
        -o "${COUNTS_DIR}/${sample_name}_counts.txt" \
        -T 8 \
        -p --countReadPairs \
        -s 2 \
        -t exon \
        -g gene_id \
        --extraAttributes gene_name,gene_type \
        -B -C \
        "${BAM_FILE}" \
        > "${COUNTS_DIR}/${sample_name}_featurecounts.log" 2>&1

    # ── DELETE BAM + BAI (no longer needed) ──
    if [[ "${KEEP_BAM}" == false ]]; then
        echo "  [CLEANUP] Removing BAM file..."
        rm -f "${BAM_FILE}" "${BAM_FILE}.bai"
        # Also remove STAR temp files
        rm -f "${OUTPREFIX}Aligned.toTranscriptome.out.bam" 2>/dev/null || true
        rm -f "${OUTPREFIX}SJ.out.tab" 2>/dev/null || true
        rm -rf "${OUTPREFIX}_STARgenome" "${OUTPREFIX}_STARpass1" 2>/dev/null || true
    fi

    echo "  Done: ${sample_name}"
    disk_usage
done


# ============================================================================
# Phase 3.1: Merge per-sample counts into a single matrix
# ============================================================================
log_step "Merging per-sample counts into combined matrix"

python3 - << 'MERGE_SCRIPT'
import os
import pandas as pd
from pathlib import Path

counts_dir = Path("/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ/data/processed/counts")
count_files = sorted(counts_dir.glob("*_counts.txt"))

if not count_files:
    print("ERROR: No per-sample count files found.")
    exit(1)

print(f"  Merging {len(count_files)} count files...")

# Read first file for gene annotations
first = pd.read_csv(count_files[0], sep="\t", comment="#", header=0)
meta_cols = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
extra_cols = [c for c in first.columns if c in ["gene_name", "gene_type"]]
annot = first[meta_cols + extra_cols].copy()

# Merge all count columns
merged = annot.copy()
for f in count_files:
    sample_name = f.stem.replace("_counts", "")
    df = pd.read_csv(f, sep="\t", comment="#", header=0)
    # Last column is the count
    count_col = [c for c in df.columns if c not in meta_cols + extra_cols][0]
    merged[sample_name] = df[count_col].values

out_path = counts_dir / "raw_counts_all.txt"
merged.to_csv(out_path, sep="\t", index=False)
print(f"  Saved: {out_path}")
print(f"  Genes: {len(merged)}, Samples: {len(count_files)}")

# Cleanup per-sample count files (keep merged only)
for f in count_files:
    f.unlink()
    summary_f = f.with_suffix(".txt.summary")
    if summary_f.exists():
        summary_f.unlink()
    log_f = counts_dir / f"{f.stem.replace('_counts', '')}_featurecounts.log"
    if log_f.exists():
        log_f.unlink()

print("  Per-sample count files cleaned up.")
MERGE_SCRIPT

# ============================================================================
# Final cleanup: remove empty trimmed directory
# ============================================================================
if [[ "${KEEP_TRIMMED}" == false ]]; then
    rmdir "${TRIMMED_DIR}" 2>/dev/null || true
fi

# ============================================================================
# Generate aggregate QC reports
# ============================================================================
log_step "Generating aggregate QC reports"

echo "  MultiQC: raw FastQC..."
multiqc "${FASTQC_DIR}/raw/" -o "${REPORT_DIR}" -n "multiqc_raw_fastqc" --force --quiet 2>/dev/null || true

echo "  MultiQC: alignment logs..."
multiqc "${ALIGNED_DIR}/logs/" -o "${REPORT_DIR}" -n "multiqc_alignment" --force --quiet 2>/dev/null || true

# ============================================================================
# Summary
# ============================================================================
echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║  Pipeline Complete                               ║"
echo "╠══════════════════════════════════════════════════╣"
echo "║  Processed : ${PROCESSED} samples                       "
echo "║  Failed    : ${FAILED} samples                          "
echo "╚══════════════════════════════════════════════════╝"
disk_usage

if [[ "${FAILED}" -gt 0 ]]; then
    echo ""
    echo "  Failed samples:${FAILED_SAMPLES}"
    echo "  Details: ${ALIGNED_DIR}/qc/qc_flagged_samples.txt"
fi

echo ""
echo "  Preserved files:"
echo "    ├─ Raw FASTQ:       ${RAW_DIR}/"
echo "    ├─ STAR logs:       ${ALIGNED_DIR}/logs/"
echo "    ├─ STAR gene counts:${STAR_GENE_COUNTS_DIR}/"
echo "    ├─ Merged counts:   ${COUNTS_DIR}/raw_counts_all.txt"
echo "    ├─ FastQC reports:  ${FASTQC_DIR}/"
echo "    └─ MultiQC reports: ${REPORT_DIR}/"
echo ""
echo "  Next: Run the R analysis scripts:"
echo "    Rscript scripts/03_quantification/02_separate_lncrna_mrna.R"
echo "    Rscript scripts/04_differential_expression/01_deseq2_analysis.R"
echo "    Rscript scripts/05_functional_analysis/01_coexpression_analysis.R"
echo "    Rscript scripts/05_functional_analysis/02_pathway_enrichment.R"
echo "    Rscript scripts/05_functional_analysis/03_wgcna_network.R"
echo "    Rscript scripts/06_survival_analysis/01_survival_analysis.R"
echo "    Rscript scripts/07_validation/01_external_validation.R"
echo "    Rscript scripts/08_figures/01_generate_main_figures.R"
