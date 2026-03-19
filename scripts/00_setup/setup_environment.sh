#!/usr/bin/env bash
# ============================================================================
# Phase 0: Environment Setup
# - Create conda environment
# - Download reference genome and annotation
# - Build STAR genome index
# ============================================================================
set -euo pipefail

PROJECT_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
REF_DIR="${PROJECT_DIR}/data/reference"
THREADS=16

echo "============================================"
echo " Phase 0: Environment Setup"
echo "============================================"

# ── 1. Create conda environment ──
echo "[1/4] Creating conda environment..."
if conda env list | grep -q "lncrna-gbm-tmz"; then
    echo "  Environment 'lncrna-gbm-tmz' already exists. Skipping."
else
    mamba env create -f "${PROJECT_DIR}/envs/environment.yaml"
fi

echo "  Activate with: conda activate lncrna-gbm-tmz"

# ── 2. Download reference genome ──
echo "[2/4] Downloading GRCh38 reference genome..."
GENOME_FA="${REF_DIR}/genome/GRCh38.primary_assembly.genome.fa"
if [[ ! -f "${GENOME_FA}" ]]; then
    wget -P "${REF_DIR}/genome/" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
    gunzip "${REF_DIR}/genome/GRCh38.primary_assembly.genome.fa.gz"
    echo "  Genome downloaded and decompressed."
else
    echo "  Genome FASTA already exists. Skipping."
fi

# ── 3. Download GENCODE annotation ──
echo "[3/4] Downloading GENCODE v44 annotation..."
GTF_FILE="${REF_DIR}/annotation/gencode.v44.annotation.gtf"
if [[ ! -f "${GTF_FILE}" ]]; then
    wget -P "${REF_DIR}/annotation/" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
    gunzip "${REF_DIR}/annotation/gencode.v44.annotation.gtf.gz"
    echo "  Annotation downloaded and decompressed."
else
    echo "  GTF annotation already exists. Skipping."
fi

# ── 4. Build STAR genome index ──
echo "[4/4] Building STAR genome index..."
STAR_INDEX="${REF_DIR}/index/star"
if [[ -f "${STAR_INDEX}/SA" ]]; then
    echo "  STAR index already exists. Skipping."
else
    mkdir -p "${STAR_INDEX}"
    STAR --runMode genomeGenerate \
         --genomeDir "${STAR_INDEX}" \
         --genomeFastaFiles "${GENOME_FA}" \
         --sjdbGTFfile "${GTF_FILE}" \
         --sjdbOverhang 100 \
         --runThreadN "${THREADS}" \
         --limitGenomeGenerateRAM 40000000000
    echo "  STAR index built successfully."
fi

echo ""
echo "============================================"
echo " Setup complete!"
echo " Next: Phase 1 - Data Acquisition"
echo "============================================"
