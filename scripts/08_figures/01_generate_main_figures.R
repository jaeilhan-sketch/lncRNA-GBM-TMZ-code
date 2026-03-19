#!/usr/bin/env Rscript
# ============================================================================
# Phase 8: Generate Publication-Quality Main Figures
#
# Fig 1: Study design (created externally, e.g., BioRender)
# Fig 2: DE-lncRNA landscape (Volcano + Heatmap composite)
# Fig 3: Top candidate lncRNA detail (expression + survival)
# Fig 4: Functional analysis (GO/KEGG/GSEA)
# Fig 5: Regulatory network summary
# Fig 6: Prognostic model + external validation
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(ComplexHeatmap)
  library(circlize)
  library(EnhancedVolcano)
  library(survminer)
  library(survival)
  library(clusterProfiler)
  library(enrichplot)
  library(data.table)
  library(gridExtra)
  library(RColorBrewer)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

# ── Publication theme ──
theme_pub <- theme_classic(base_size = 12) +
  theme(
    text = element_text(family = "sans"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 11, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# Color palette
colors <- list(
  responder = "#2166AC",
  non_responder = "#B2182B",
  up = "#D73027",
  down = "#4575B4",
  ns = "#999999",
  high = "#B2182B",
  low = "#2166AC"
)

save_fig <- function(plot, name, w, h) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), plot,
         width = w, height = h, units = "mm", dpi = 300)
  ggsave(file.path(fig_dir, paste0(name, ".tiff")), plot,
         width = w, height = h, units = "mm", dpi = 300, compression = "lzw")
}

cat("============================================\n")
cat(" Phase 8: Generate Publication Figures\n")
cat("============================================\n")

# ── Load all results ──
cat("[1/6] Loading analysis results...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
load(file.path(results_dir, "survival_objects.RData"))

if (file.exists(file.path(results_dir, "enrichment_objects.RData")))
  load(file.path(results_dir, "enrichment_objects.RData"))
if (file.exists(file.path(results_dir, "wgcna_objects.RData")))
  load(file.path(results_dir, "wgcna_objects.RData"))


# ============================================================================
# Figure 2: DE-lncRNA Landscape
# ============================================================================
cat("[2/6] Figure 2: DE-lncRNA landscape...\n")

# ── 2A: Volcano plot ──
res_plot <- res_lncrna_df[!is.na(res_lncrna_df$padj), ]
res_plot$significance <- "NS"
res_plot$significance[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Up"
res_plot$significance[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Down"
res_plot$significance <- factor(res_plot$significance, levels = c("NS", "Up", "Down"))

# Label top genes
top_up <- head(res_plot[res_plot$significance == "Up", ][order(res_plot$padj[res_plot$significance == "Up"]), ], 10)
top_down <- head(res_plot[res_plot$significance == "Down", ][order(res_plot$padj[res_plot$significance == "Down"]), ], 10)
res_plot$label <- ""
res_plot$label[res_plot$gene_id %in% c(top_up$gene_id, top_down$gene_id)] <-
  res_plot$gene_name[res_plot$gene_id %in% c(top_up$gene_id, top_down$gene_id)]

fig2a <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj),
                               color = significance)) +
  geom_point(size = 1.2, alpha = 0.6) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                   segment.size = 0.3) +
  scale_color_manual(values = c("NS" = colors$ns, "Up" = colors$up, "Down" = colors$down),
                     labels = c(paste0("NS (", sum(res_plot$significance == "NS"), ")"),
                                paste0("Up (", sum(res_plot$significance == "Up"), ")"),
                                paste0("Down (", sum(res_plot$significance == "Down"), ")"))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(x = expression(log[2]~"Fold Change (Responder/Non-Responder)"),
       y = expression(-log[10]~"adjusted p-value"),
       color = "Significance",
       title = "A") +
  theme_pub +
  theme(legend.position = c(0.85, 0.85))

save_fig(fig2a, "Fig2A_volcano", 180, 150)

# ── 2B: Heatmap (top 30 DE-lncRNAs) ──
if (nrow(sig_lncrna) > 0) {
  top_n <- min(30, nrow(sig_lncrna))
  top_ids <- sig_lncrna$gene_id[seq_len(top_n)]
  vsd_mat <- assay(vsd_lncrna)
  hm_mat <- t(scale(t(vsd_mat[top_ids, ])))

  # Column annotation
  ha_col <- HeatmapAnnotation(
    Response = col_data$tmz_response,
    col = list(Response = c("Responder" = colors$responder,
                             "NonResponder" = colors$non_responder)),
    annotation_name_side = "left",
    simple_anno_size = unit(4, "mm")
  )

  gene_labels <- sig_lncrna$gene_name[seq_len(top_n)]
  gene_labels[is.na(gene_labels)] <- sig_lncrna$gene_id[which(is.na(gene_labels[seq_len(top_n)]))]

  ht <- Heatmap(
    hm_mat,
    name = "Z-score",
    top_annotation = ha_col,
    row_labels = gene_labels,
    row_names_gp = gpar(fontsize = 8),
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    col = colorRamp2(c(-2, 0, 2), c(colors$down, "white", colors$up)),
    column_title = "B",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Z-score",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 9)
    )
  )

  pdf(file.path(fig_dir, "Fig2B_heatmap.pdf"), width = 8, height = 10)
  draw(ht)
  dev.off()

  tiff(file.path(fig_dir, "Fig2B_heatmap.tiff"), width = 8, height = 10,
       units = "in", res = 300, compression = "lzw")
  draw(ht)
  dev.off()
}

cat("  Figure 2 saved.\n")


# ============================================================================
# Figure 3: Top Candidate Detail
# ============================================================================
cat("[3/6] Figure 3: Top candidate lncRNA detail...\n")

if (nrow(sig_lncrna) > 0 && nrow(km_df) > 0) {
  # Pick the top candidate (best survival + DE)
  best <- km_df[km_df$km_pvalue == min(km_df$km_pvalue, na.rm = TRUE), ][1, ]
  gene_id <- best$gene_id
  gene_name <- best$gene_name

  expr_vals <- as.numeric(assay(vsd_lncrna)[gene_id, ])
  surv_data_plot <- col_data
  surv_data_plot$expression <- expr_vals

  # ── 3A: Boxplot ──
  fig3a <- ggplot(surv_data_plot, aes(x = tmz_response, y = expression,
                                       fill = tmz_response)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
    geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("NonResponder" = colors$non_responder,
                                  "Responder" = colors$responder)) +
    labs(x = "", y = "Normalized Expression (VST)",
         title = paste0("A: ", gene_name)) +
    theme_pub +
    theme(legend.position = "none")

  # ── 3B: KM curve ──
  surv_data_plot$os_time <- as.numeric(surv_data_plot$os_days)
  surv_data_plot$os_status <- as.numeric(surv_data_plot$os_status)
  surv_data_plot <- surv_data_plot[!is.na(surv_data_plot$os_time), ]
  surv_data_plot$group <- ifelse(surv_data_plot$expression > median(surv_data_plot$expression),
                                  "High", "Low")

  fit <- survfit(Surv(os_time, os_status) ~ group, data = surv_data_plot)
  fig3b <- ggsurvplot(
    fit, data = surv_data_plot,
    pval = TRUE, risk.table = TRUE,
    palette = c(colors$low, colors$high),
    title = paste0("B: ", gene_name, " - Overall Survival"),
    xlab = "Time (days)", ylab = "Survival Probability",
    legend.labs = c("High expression", "Low expression"),
    ggtheme = theme_pub
  )

  pdf(file.path(fig_dir, "Fig3_top_candidate.pdf"), width = 12, height = 6)
  grid.arrange(
    fig3a, fig3b$plot,
    ncol = 2, widths = c(1, 1.5)
  )
  dev.off()

  cat(sprintf("  Figure 3 saved: %s\n", gene_name))
}


# ============================================================================
# Figure 4: Functional Analysis
# ============================================================================
cat("[4/6] Figure 4: Functional analysis...\n")

if (exists("ego_bp") && !is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  # ── 4A: GO BP Dotplot ──
  fig4a <- dotplot(ego_bp, showCategory = 15,
                    title = "A: GO Biological Process") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))
  save_fig(fig4a, "Fig4A_GO_BP", 200, 220)
}

if (exists("ekegg") && !is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  # ── 4B: KEGG Dotplot ──
  fig4b <- dotplot(ekegg, showCategory = 15,
                    title = "B: KEGG Pathway") +
    theme_pub +
    theme(axis.text.y = element_text(size = 9))
  save_fig(fig4b, "Fig4B_KEGG", 200, 200)
}

if (exists("gsea_go") && !is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  # ── 4C: GSEA plot ──
  fig4c <- gseaplot2(gsea_go, geneSetID = 1:3,
                      title = "C: Gene Set Enrichment Analysis")
  ggsave(file.path(fig_dir, "Fig4C_GSEA.pdf"), fig4c,
         width = 200, height = 150, units = "mm", dpi = 300)
}

cat("  Figure 4 saved.\n")


# ============================================================================
# Figure 5: Network (placeholder for Cytoscape)
# ============================================================================
cat("[5/6] Figure 5: Network visualization...\n")
cat("  NOTE: Import cytoscape_edges.tsv and cytoscape_nodes.tsv\n")
cat("  into Cytoscape for publication-quality network figure.\n")


# ============================================================================
# Figure 6: Prognostic Model
# ============================================================================
cat("[6/6] Figure 6: Prognostic model...\n")
cat("  Figures generated in Phase 6 (lasso_cox_risk_km.pdf, time_dependent_roc.pdf).\n")
cat("  Combine into composite Figure 6 for manuscript.\n")

cat("\n============================================\n")
cat(" All main figures generated.\n")
cat(sprintf(" Figure directory: %s\n", fig_dir))
cat("============================================\n")

dir.create(file.path(results_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "reports/figures_session_info.txt"))
