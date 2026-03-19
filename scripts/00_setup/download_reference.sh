#!/usr/bin/env bash
# ============================================================================
# Download additional reference files (optional)
# - lncRNA-specific annotation
# - BED files for RSeQC
# ============================================================================
set -euo pipefail

REF_DIR="/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ/data/reference"

# ── lncRNA annotation (GENCODE lncRNA-only GTF) ──
echo "Downloading GENCODE v44 lncRNA annotation..."
if [[ ! -f "${REF_DIR}/annotation/gencode.v44.long_noncoding_RNAs.gtf" ]]; then
    wget -P "${REF_DIR}/annotation/" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.long_noncoding_RNAs.gtf.gz"
    gunzip "${REF_DIR}/annotation/gencode.v44.long_noncoding_RNAs.gtf.gz"
fi

# ── BED12 for RSeQC ──
echo "Creating BED12 from GTF for RSeQC..."
if [[ ! -f "${REF_DIR}/annotation/gencode.v44.annotation.bed" ]]; then
    # Requires gtfToGenePred and genePredToBed from UCSC utilities
    # Alternative: use pybedtools or RSeQC's gtf2bed script
    echo "  NOTE: Convert GTF to BED12 manually if UCSC tools are not installed."
    echo "  Command: gtfToGenePred gencode.v44.annotation.gtf /dev/stdout | genePredToBed /dev/stdin gencode.v44.annotation.bed"
fi

echo "Reference download complete."
