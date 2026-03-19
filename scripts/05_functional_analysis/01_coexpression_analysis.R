#!/usr/bin/env Rscript
# ============================================================================
# Phase 5.1: lncRNA-mRNA Co-expression Analysis (Guilt-by-Association)
#
# Identifies mRNAs co-expressed with DE-lncRNAs to infer function.
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

# ── Parameters ──
COR_METHOD   <- "spearman"
COR_CUTOFF   <- 0.6
PVAL_CUTOFF  <- 0.01
TOP_LNCRNA_N <- 50

cat("============================================\n")
cat(" Phase 5.1: Co-expression Analysis\n")
cat("============================================\n")

# ── Load DE results ──
cat("[1/4] Loading data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))

# Get VST-normalized expression
vsd_mrna <- vst(dds_mrna, blind = TRUE)
expr_lncrna <- assay(vsd_lncrna)
expr_mrna <- assay(vsd_mrna)

# Focus on significant DE-lncRNAs
top_lncrna <- head(sig_lncrna$gene_id, TOP_LNCRNA_N)
expr_lncrna_top <- expr_lncrna[top_lncrna, , drop = FALSE]

cat(sprintf("  DE-lncRNAs for co-expression: %d\n", length(top_lncrna)))
cat(sprintf("  mRNAs available: %d\n", nrow(expr_mrna)))

# ── Compute correlations ──
cat("[2/4] Computing lncRNA-mRNA correlations...\n")

results_list <- list()

for (i in seq_along(top_lncrna)) {
  lnc_id <- top_lncrna[i]
  lnc_expr <- as.numeric(expr_lncrna_top[lnc_id, ])

  # Correlation with all mRNAs
  cor_vals <- apply(expr_mrna, 1, function(mrna_expr) {
    ct <- cor.test(lnc_expr, as.numeric(mrna_expr), method = COR_METHOD)
    c(cor = ct$estimate, pval = ct$p.value)
  })

  cor_df <- data.frame(
    lncrna_id = lnc_id,
    mrna_id = colnames(cor_vals),
    correlation = cor_vals["cor.rho", ],
    pvalue = cor_vals["pval", ],
    stringsAsFactors = FALSE
  )
  # Handle Pearson naming
  if (is.null(cor_df$correlation)) {
    cor_df$correlation <- cor_vals["cor.cor", ]
  }

  # BH correction per lncRNA
  cor_df$padj <- p.adjust(cor_df$pvalue, method = "BH")

  # Filter significant correlations
  cor_df_sig <- cor_df[abs(cor_df$correlation) >= COR_CUTOFF &
                        cor_df$padj < PVAL_CUTOFF, ]

  results_list[[lnc_id]] <- cor_df_sig

  if (i %% 10 == 0) cat(sprintf("  Processed %d / %d lncRNAs\n", i, length(top_lncrna)))
}

coexpr_all <- rbindlist(results_list)

# Add gene names
coexpr_all <- merge(coexpr_all,
                     gene_annot[, c("gene_id", "gene_name")],
                     by.x = "lncrna_id", by.y = "gene_id", all.x = TRUE)
setnames(coexpr_all, "gene_name", "lncrna_name")

coexpr_all <- merge(coexpr_all,
                     gene_annot[, c("gene_id", "gene_name")],
                     by.x = "mrna_id", by.y = "gene_id", all.x = TRUE)
setnames(coexpr_all, "gene_name", "mrna_name")

cat(sprintf("  Total significant lncRNA-mRNA pairs: %d\n", nrow(coexpr_all)))
cat(sprintf("  Positive correlations: %d\n", sum(coexpr_all$correlation > 0)))
cat(sprintf("  Negative correlations: %d\n", sum(coexpr_all$correlation < 0)))

# ── Save co-expression results ──
cat("[3/4] Saving results...\n")
fwrite(coexpr_all, file.path(table_dir, "coexpression_lncrna_mrna.tsv"), sep = "\t")

# ── Per-lncRNA co-expressed gene lists (for enrichment) ──
cat("[4/4] Generating per-lncRNA gene lists for enrichment...\n")

gene_lists_dir <- file.path(results_dir, "coexpr_gene_lists")
dir.create(gene_lists_dir, showWarnings = FALSE)

for (lnc_id in unique(coexpr_all$lncrna_id)) {
  genes <- coexpr_all[lncrna_id == lnc_id, mrna_id]
  lnc_name <- coexpr_all[lncrna_id == lnc_id, lncrna_name][1]
  if (is.na(lnc_name)) lnc_name <- lnc_id

  writeLines(genes, file.path(gene_lists_dir,
                               paste0(gsub("[^A-Za-z0-9]", "_", lnc_name), "_coexpr_genes.txt")))
}

cat(sprintf("  Gene lists saved for %d lncRNAs\n", length(unique(coexpr_all$lncrna_id))))

cat("\n============================================\n")
cat(" Co-expression analysis complete.\n")
cat("============================================\n")
