#!/usr/bin/env Rscript
# ============================================================================
# Phase 4: Differential Expression Analysis with DESeq2
#
# Performs DE analysis for lncRNA and mRNA separately.
# Includes batch effect assessment, SVA correction, and visualization.
#
# Input:
#   data/processed/counts/counts_lncrna.tsv
#   data/processed/counts/counts_mrna.tsv
#   data/clinical/sample_info_final.tsv
#
# Output:
#   results/tables/de_lncrna_results.tsv
#   results/tables/de_mrna_results.tsv
#   results/tables/de_lncrna_significant.tsv
#   results/figures/volcano_lncrna.pdf
#   results/figures/heatmap_top_lncrna.pdf
#   results/figures/pca_plot.pdf
#   results/figures/ma_plot.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(sva)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(EnhancedVolcano)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
})

set.seed(42)

# ── Paths ──
project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
counts_dir  <- file.path(project_dir, "data/processed/counts")
clinical_dir <- file.path(project_dir, "data/clinical")
results_dir <- file.path(project_dir, "results")
fig_dir     <- file.path(results_dir, "figures")
table_dir   <- file.path(results_dir, "tables")
report_dir  <- file.path(results_dir, "reports")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# ── Parameters ──
ALPHA <- 0.05
LFC_THRESHOLD <- 1.0
LFC_RELAXED <- 0.585

# ── Publication theme ──
theme_pub <- theme_classic() +
  theme(
    text = element_text(size = 12, family = "sans"),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 11)
  )

cat("============================================\n")
cat(" Phase 4: Differential Expression Analysis\n")
cat("============================================\n")

# ============================================================================
# 1. Load data
# ============================================================================
cat("[1/7] Loading data...\n")

# Count matrices
counts_lncrna <- as.matrix(read.table(
  file.path(counts_dir, "counts_lncrna.tsv"), header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
))

counts_mrna <- as.matrix(read.table(
  file.path(counts_dir, "counts_mrna.tsv"), header = TRUE, row.names = 1,
  sep = "\t", check.names = FALSE
))

# Gene annotation
gene_annot <- fread(file.path(counts_dir, "gene_annotation.tsv"))

# Clinical / sample info
sample_info <- fread(file.path(clinical_dir, "sample_info_final.tsv"))

# Match sample order
common_samples <- intersect(colnames(counts_lncrna), sample_info$submitter_id)
if (length(common_samples) == 0) {
  # Try matching by other ID columns
  for (id_col in c("sample_id", "case_id", "submitter_id")) {
    common_samples <- intersect(colnames(counts_lncrna), sample_info[[id_col]])
    if (length(common_samples) > 0) break
  }
}

cat(sprintf("  Matched samples: %d\n", length(common_samples)))

# Subset and reorder
sample_info <- sample_info[match(common_samples, sample_info$submitter_id), ]
counts_lncrna <- counts_lncrna[, common_samples]
counts_mrna <- counts_mrna[, common_samples]

# Prepare factors
sample_info$tmz_response <- factor(sample_info$tmz_response,
                                    levels = c("NonResponder", "Responder"))
sample_info$gender <- factor(sample_info$gender)
sample_info$age_group <- factor(sample_info$age_group)
sample_info$MGMT_status <- factor(sample_info$MGMT_status)

cat(sprintf("  Responder: %d, NonResponder: %d\n",
            sum(sample_info$tmz_response == "Responder"),
            sum(sample_info$tmz_response == "NonResponder")))


# ============================================================================
# 2. DESeq2 - lncRNA analysis
# ============================================================================
cat("[2/7] Running DESeq2 for lncRNA...\n")

col_data <- as.data.frame(sample_info)
rownames(col_data) <- col_data$submitter_id

# Build design based on available covariates
# NOTE: MGMT_status has 25 "Unknown" values out of 94 samples.
# Including it would drop 25 samples. Instead, use gender as covariate
# and let SVA capture other confounders (including MGMT-related variance).
has_gender <- "gender" %in% colnames(col_data) && !all(is.na(col_data$gender))

if (has_gender) {
  design_formula <- ~ gender + tmz_response
} else {
  design_formula <- ~ tmz_response
}

cat(sprintf("  Design: %s\n", deparse(design_formula)))

# Create DESeqDataSet
dds_lncrna <- DESeqDataSetFromMatrix(
  countData = counts_lncrna,
  colData = col_data,
  design = design_formula
)

# Run DESeq2
dds_lncrna <- DESeq(dds_lncrna)

# Extract results
res_lncrna <- results(dds_lncrna,
                       contrast = c("tmz_response", "Responder", "NonResponder"),
                       alpha = ALPHA)

# Add gene names
res_lncrna_df <- as.data.frame(res_lncrna)
res_lncrna_df$gene_id <- rownames(res_lncrna_df)
res_lncrna_df$gene_id_noversion <- sub("\\.\\d+$", "", res_lncrna_df$gene_id)
res_lncrna_df <- merge(res_lncrna_df,
                        gene_annot[, c("gene_id", "gene_name", "gene_type")],
                        by = "gene_id", all.x = TRUE)

# Sort by padj
res_lncrna_df <- res_lncrna_df[order(res_lncrna_df$padj), ]

# Significant DE-lncRNAs
sig_lncrna <- res_lncrna_df[!is.na(res_lncrna_df$padj) &
                              res_lncrna_df$padj < ALPHA &
                              abs(res_lncrna_df$log2FoldChange) > LFC_THRESHOLD, ]

cat(sprintf("  Total DE-lncRNA (padj < %.2f, |log2FC| > %.1f): %d\n",
            ALPHA, LFC_THRESHOLD, nrow(sig_lncrna)))
cat(sprintf("    Upregulated in Responder: %d\n",
            sum(sig_lncrna$log2FoldChange > 0)))
cat(sprintf("    Downregulated in Responder: %d\n",
            sum(sig_lncrna$log2FoldChange < 0)))


# ============================================================================
# 3. DESeq2 - mRNA analysis (for comparison)
# ============================================================================
cat("[3/7] Running DESeq2 for mRNA...\n")

dds_mrna <- DESeqDataSetFromMatrix(
  countData = counts_mrna,
  colData = col_data,
  design = design_formula
)
dds_mrna <- DESeq(dds_mrna)
res_mrna <- results(dds_mrna,
                     contrast = c("tmz_response", "Responder", "NonResponder"),
                     alpha = ALPHA)

res_mrna_df <- as.data.frame(res_mrna)
res_mrna_df$gene_id <- rownames(res_mrna_df)
res_mrna_df <- merge(res_mrna_df,
                      gene_annot[, c("gene_id", "gene_name", "gene_type")],
                      by = "gene_id", all.x = TRUE)
res_mrna_df <- res_mrna_df[order(res_mrna_df$padj), ]

sig_mrna <- res_mrna_df[!is.na(res_mrna_df$padj) &
                          res_mrna_df$padj < ALPHA &
                          abs(res_mrna_df$log2FoldChange) > LFC_THRESHOLD, ]
cat(sprintf("  Total DE-mRNA: %d\n", nrow(sig_mrna)))


# ============================================================================
# 4. Batch effect assessment with PCA
# ============================================================================
cat("[4/7] Assessing batch effects (PCA)...\n")

vsd_lncrna <- vst(dds_lncrna, blind = TRUE)

# PCA plot
pca_data <- plotPCA(vsd_lncrna, intgroup = "tmz_response", returnData = TRUE)
pct_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = tmz_response)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("NonResponder" = "#B2182B", "Responder" = "#2166AC"),
                     name = "TMZ Response") +
  labs(
    x = paste0("PC1 (", pct_var[1], "% variance)"),
    y = paste0("PC2 (", pct_var[2], "% variance)"),
    title = "PCA: lncRNA Expression by TMZ Response"
  ) +
  theme_pub

ggsave(file.path(fig_dir, "pca_lncrna_tmz_response.pdf"), p_pca,
       width = 180, height = 140, units = "mm", dpi = 300)
ggsave(file.path(fig_dir, "pca_lncrna_tmz_response.tiff"), p_pca,
       width = 180, height = 140, units = "mm", dpi = 300, compression = "lzw")
cat("  PCA plot saved.\n")


# ============================================================================
# 5. SVA batch correction (if needed)
# ============================================================================
cat("[5/7] Surrogate Variable Analysis...\n")

mod <- model.matrix(design_formula, data = col_data)
mod0 <- model.matrix(~ 1, data = col_data)

# Normalized counts for SVA
norm_counts <- counts(dds_lncrna, normalized = TRUE)

n_sv <- num.sv(norm_counts, mod, method = "be")
cat(sprintf("  Estimated surrogate variables: %d\n", n_sv))

if (n_sv > 0) {
  svobj <- sva(norm_counts, mod, mod0, n.sv = n_sv)

  # Add SVs to design and re-run
  col_data_sv <- col_data
  for (i in seq_len(n_sv)) {
    col_data_sv[[paste0("SV", i)]] <- svobj$sv[, i]
  }

  sv_terms <- paste0("SV", seq_len(n_sv), collapse = " + ")
  design_sv <- as.formula(paste0(deparse(design_formula), " + ", sv_terms))

  cat(sprintf("  Re-running DESeq2 with SVA: %s\n", deparse(design_sv)))

  dds_lncrna_sv <- DESeqDataSetFromMatrix(
    countData = counts_lncrna,
    colData = col_data_sv,
    design = design_sv
  )
  dds_lncrna_sv <- DESeq(dds_lncrna_sv)

  res_lncrna_sv <- results(dds_lncrna_sv,
                            contrast = c("tmz_response", "Responder", "NonResponder"),
                            alpha = ALPHA)

  res_lncrna_sv_df <- as.data.frame(res_lncrna_sv)
  res_lncrna_sv_df$gene_id <- rownames(res_lncrna_sv_df)
  res_lncrna_sv_df <- merge(res_lncrna_sv_df,
                              gene_annot[, c("gene_id", "gene_name", "gene_type")],
                              by = "gene_id", all.x = TRUE)

  sig_lncrna_sv <- res_lncrna_sv_df[!is.na(res_lncrna_sv_df$padj) &
                                      res_lncrna_sv_df$padj < ALPHA &
                                      abs(res_lncrna_sv_df$log2FoldChange) > LFC_THRESHOLD, ]

  cat(sprintf("  DE-lncRNA after SVA correction: %d\n", nrow(sig_lncrna_sv)))

  # Use SVA-corrected results as primary
  res_lncrna_df <- res_lncrna_sv_df[order(res_lncrna_sv_df$padj), ]
  sig_lncrna <- sig_lncrna_sv
}


# If no SVA was performed, set res_lncrna_final
if (!exists("res_lncrna_sv")) {
  res_lncrna_final <- res_lncrna
}

# Also apply relaxed threshold for additional candidates
sig_lncrna_relaxed <- res_lncrna_df[!is.na(res_lncrna_df$padj) &
                                      res_lncrna_df$padj < ALPHA &
                                      abs(res_lncrna_df$log2FoldChange) > LFC_RELAXED, ]
cat(sprintf("  Relaxed DE-lncRNA (|log2FC| > %.3f): %d\n", LFC_RELAXED, nrow(sig_lncrna_relaxed)))

# ============================================================================
# 6. Visualization
# ============================================================================
cat("[6/7] Generating figures...\n")

# ── Volcano Plot ──
p_volcano <- EnhancedVolcano(
  res_lncrna_df,
  lab = res_lncrna_df$gene_name,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = ALPHA,
  FCcutoff = LFC_THRESHOLD,
  pointSize = 2.0,
  labSize = 3.5,
  colAlpha = 0.6,
  title = "DE-lncRNAs: Responder vs Non-Responder",
  subtitle = paste0("padj < ", ALPHA, ", |log2FC| > ", LFC_THRESHOLD),
  col = c("grey30", "#4575B4", "#D73027", "#FFD700"),
  legendPosition = "right"
)

ggsave(file.path(fig_dir, "volcano_lncrna.pdf"), p_volcano,
       width = 200, height = 180, units = "mm", dpi = 300)
ggsave(file.path(fig_dir, "volcano_lncrna.tiff"), p_volcano,
       width = 200, height = 180, units = "mm", dpi = 300, compression = "lzw")
cat("  Volcano plot saved.\n")

# ── MA Plot ──  (use the latest results object)
res_lncrna_final <- if (exists("res_lncrna_sv")) res_lncrna_sv else res_lncrna
pdf(file.path(fig_dir, "ma_plot_lncrna.pdf"), width = 7, height = 5)
plotMA(res_lncrna_final, ylim = c(-5, 5), main = "MA Plot: lncRNA DE Analysis")
dev.off()
cat("  MA plot saved.\n")

# ── Heatmap (top 50 DE-lncRNAs) ──
if (nrow(sig_lncrna) > 0) {
  top_n <- min(50, nrow(sig_lncrna))
  top_genes <- sig_lncrna$gene_id[seq_len(top_n)]

  # Use VST-transformed values
  vsd_mat <- assay(vsd_lncrna)
  heatmap_mat <- vsd_mat[top_genes, ]

  # Z-score normalization per gene
  heatmap_scaled <- t(scale(t(heatmap_mat)))

  # Annotation
  ha <- HeatmapAnnotation(
    Response = col_data$tmz_response,
    col = list(Response = c("Responder" = "#2166AC", "NonResponder" = "#B2182B")),
    annotation_name_side = "left"
  )

  # Gene labels
  gene_labels <- sig_lncrna$gene_name[seq_len(top_n)]
  gene_labels[is.na(gene_labels)] <- sig_lncrna$gene_id[which(is.na(gene_labels))]

  ht <- Heatmap(
    heatmap_scaled,
    name = "Z-score",
    top_annotation = ha,
    row_labels = gene_labels,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    col = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    column_title = paste0("Top ", top_n, " DE-lncRNAs")
  )

  pdf(file.path(fig_dir, "heatmap_top_lncrna.pdf"), width = 10, height = 12)
  draw(ht)
  dev.off()

  tiff(file.path(fig_dir, "heatmap_top_lncrna.tiff"), width = 10, height = 12,
       units = "in", res = 300, compression = "lzw")
  draw(ht)
  dev.off()

  cat("  Heatmap saved.\n")
}


# ============================================================================
# 7. Save results
# ============================================================================
cat("[7/7] Saving result tables...\n")

# Full results
fwrite(res_lncrna_df, file.path(table_dir, "de_lncrna_results_full.tsv"), sep = "\t")
fwrite(res_mrna_df, file.path(table_dir, "de_mrna_results_full.tsv"), sep = "\t")

# Significant only
fwrite(sig_lncrna, file.path(table_dir, "de_lncrna_significant.tsv"), sep = "\t")
fwrite(sig_lncrna_relaxed, file.path(table_dir, "de_lncrna_significant_relaxed.tsv"), sep = "\t")
fwrite(sig_mrna, file.path(table_dir, "de_mrna_significant.tsv"), sep = "\t")

# Save DESeq2 objects for downstream use
save(dds_lncrna, vsd_lncrna, res_lncrna_final, res_lncrna_df, sig_lncrna,
     dds_mrna, res_mrna, res_mrna_df, sig_mrna,
     col_data, gene_annot,
     file = file.path(results_dir, "de_analysis_objects.RData"))

# Summary
cat(sprintf("\n  Summary:\n"))
cat(sprintf("  ├─ DE-lncRNA (|log2FC| > %.1f, padj < %.2f): %d\n",
            LFC_THRESHOLD, ALPHA, nrow(sig_lncrna)))
cat(sprintf("  │  ├─ Up in Responder:   %d\n",
            sum(sig_lncrna$log2FoldChange > 0)))
cat(sprintf("  │  └─ Down in Responder: %d\n",
            sum(sig_lncrna$log2FoldChange < 0)))
cat(sprintf("  └─ DE-mRNA (reference):  %d\n", nrow(sig_mrna)))

cat("\n============================================\n")
cat(" Phase 4 complete.\n")
cat("============================================\n")

# Session info
writeLines(capture.output(sessionInfo()),
           file.path(report_dir, "deseq2_session_info.txt"))
