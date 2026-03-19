#!/usr/bin/env Rscript
# ============================================================================
# Phase 3.2: Separate lncRNA and mRNA count matrices
#
# Input:  data/processed/counts/raw_counts_all.txt (featureCounts output)
# Output: data/processed/counts/counts_lncrna.tsv
#         data/processed/counts/counts_mrna.tsv
#         data/processed/counts/gene_annotation.tsv
#         results/reports/count_distribution_summary.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(edgeR)
})

set.seed(42)

# ── Paths ──
project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
counts_dir  <- file.path(project_dir, "data/processed/counts")
ref_dir     <- file.path(project_dir, "data/reference/annotation")
results_dir <- file.path(project_dir, "results")

# ── lncRNA biotypes to include ──
lncrna_biotypes <- c(
  "lincRNA", "antisense", "sense_intronic", "sense_overlapping",
  "processed_transcript", "bidirectional_promoter_lncRNA", "macro_lncRNA",
  "3prime_overlapping_ncRNA", "non_coding"
)

cat("============================================\n")
cat(" Phase 3.2: Separate lncRNA / mRNA counts\n")
cat("============================================\n")

# ── 1. Load featureCounts output ──
cat("[1/5] Loading featureCounts output...\n")
fc <- fread(file.path(counts_dir, "raw_counts_all.txt"), skip = 1)

# Extract annotation columns
annot_cols <- c("Geneid", "gene_name", "gene_type", "Chr", "Start", "End",
                "Strand", "Length")
annot_cols_present <- intersect(annot_cols, colnames(fc))

# Count columns are BAM file paths → clean to sample names
bam_cols <- setdiff(colnames(fc), c(annot_cols_present))
sample_names <- gsub("_Aligned\\.sortedByCoord\\.out\\.bam$", "", basename(bam_cols))
sample_names <- gsub("^.*\\/", "", sample_names)
setnames(fc, bam_cols, sample_names)

cat(sprintf("  Total genes: %d\n", nrow(fc)))
cat(sprintf("  Total samples: %d\n", length(sample_names)))

# ── 2. Build gene annotation table ──
cat("[2/5] Building gene annotation table...\n")
gene_annot <- fc[, ..annot_cols_present]
setnames(gene_annot, "Geneid", "gene_id")

# Strip ENSG version suffix for compatibility
gene_annot[, gene_id_noversion := sub("\\.\\d+$", "", gene_id)]

fwrite(gene_annot, file.path(counts_dir, "gene_annotation.tsv"), sep = "\t")
cat(sprintf("  Gene annotation saved: %d genes\n", nrow(gene_annot)))

# ── 3. Separate lncRNA and mRNA ──
cat("[3/5] Separating lncRNA and mRNA...\n")

# Count matrix (genes x samples)
count_mat <- as.matrix(fc[, ..sample_names])
rownames(count_mat) <- fc$Geneid

# lncRNA
is_lncrna <- gene_annot$gene_type %in% lncrna_biotypes
counts_lncrna <- count_mat[is_lncrna, , drop = FALSE]
cat(sprintf("  lncRNA genes: %d\n", nrow(counts_lncrna)))

# mRNA (protein_coding)
is_mrna <- gene_annot$gene_type == "protein_coding"
counts_mrna <- count_mat[is_mrna, , drop = FALSE]
cat(sprintf("  mRNA genes: %d\n", nrow(counts_mrna)))

# ── 4. Low-expression filtering ──
cat("[4/5] Filtering low-expression genes...\n")

# --- lncRNA: relaxed filtering ---
# Keep if CPM > 0.5 in at least 30% of samples
min_samples_lncrna <- ceiling(ncol(counts_lncrna) * 0.3)
cpm_lncrna <- cpm(counts_lncrna)
keep_lncrna <- rowSums(cpm_lncrna > 0.5) >= min_samples_lncrna
counts_lncrna_filt <- counts_lncrna[keep_lncrna, , drop = FALSE]
cat(sprintf("  lncRNA after filtering: %d (removed %d)\n",
            nrow(counts_lncrna_filt),
            sum(!keep_lncrna)))

# --- mRNA: standard filtering ---
min_samples_mrna <- ceiling(ncol(counts_mrna) * 0.5)
cpm_mrna <- cpm(counts_mrna)
keep_mrna <- rowSums(cpm_mrna > 1) >= min_samples_mrna
counts_mrna_filt <- counts_mrna[keep_mrna, , drop = FALSE]
cat(sprintf("  mRNA after filtering: %d (removed %d)\n",
            nrow(counts_mrna_filt),
            sum(!keep_mrna)))

# ── 5. Save ──
cat("[5/5] Saving count matrices...\n")

# Raw (unfiltered)
write.table(counts_lncrna, file.path(counts_dir, "counts_lncrna_raw.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(counts_mrna, file.path(counts_dir, "counts_mrna_raw.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Filtered
write.table(counts_lncrna_filt, file.path(counts_dir, "counts_lncrna.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(counts_mrna_filt, file.path(counts_dir, "counts_mrna.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

# Summary table
summary_df <- data.frame(
  Category = c("Total genes", "lncRNA (raw)", "lncRNA (filtered)",
                "mRNA (raw)", "mRNA (filtered)"),
  Count = c(nrow(count_mat), nrow(counts_lncrna), nrow(counts_lncrna_filt),
            nrow(counts_mrna), nrow(counts_mrna_filt))
)
fwrite(summary_df, file.path(results_dir, "reports/count_distribution_summary.tsv"),
       sep = "\t")

cat("\n============================================\n")
cat(" Count separation complete.\n")
cat(sprintf(" lncRNA: %s\n", file.path(counts_dir, "counts_lncrna.tsv")))
cat(sprintf(" mRNA:   %s\n", file.path(counts_dir, "counts_mrna.tsv")))
cat("============================================\n")
