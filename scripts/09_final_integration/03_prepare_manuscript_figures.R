#!/usr/bin/env Rscript
# ============================================================================
# Prepare manuscript figures for Neuro-Oncology Advances submission
# Convert PDF figures to TIFF 600dpi and organize in manuscript/figures/
# ============================================================================

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
fig_src     <- file.path(project_dir, "results/figures")
fig_dst     <- file.path(project_dir, "manuscript/figures")
sup_dst     <- file.path(project_dir, "manuscript/supplementary")

dir.create(fig_dst, recursive = TRUE, showWarnings = FALSE)
dir.create(sup_dst, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat(" Preparing Manuscript Figures for NOA\n")
cat("============================================\n\n")

# ── Define figure mapping: source PDF → destination TIFF ──
main_figures <- list(
  list(src = "Fig1_study_design.pdf",         dst = "Figure1.tiff", label = "Study Design"),
  list(src = "Fig2A_volcano.pdf",             dst = "Figure2A.tiff", label = "Volcano Plot"),
  list(src = "Fig2B_heatmap.pdf",             dst = "Figure2B.tiff", label = "Heatmap"),
  list(src = "Fig3_top_candidate.pdf",        dst = "Figure3.tiff",  label = "Top Candidate"),
  list(src = "Fig4A_GO_BP.pdf",               dst = "Figure4A.tiff", label = "GO BP Dotplot"),
  list(src = "Fig4B_KEGG.pdf",                dst = "Figure4B.tiff", label = "KEGG Dotplot"),
  list(src = "Fig4C_GSEA.pdf",                dst = "Figure4C.tiff", label = "GSEA Plot"),
  list(src = "Fig5_cerna_network.pdf",        dst = "Figure5A.tiff", label = "ceRNA Network"),
  list(src = "Fig5B_cerna_by_lncrna.pdf",     dst = "Figure5B.tiff", label = "ceRNA per lncRNA"),
  list(src = "Fig6A_forest_plot.pdf",         dst = "Figure6.tiff",  label = "Forest Plot"),
  list(src = "Fig7_candidate_evidence_heatmap.pdf", dst = "Figure7.tiff", label = "Evidence Heatmap")
)

supp_figures <- list(
  list(src = "pca_lncrna_tmz_response.pdf",       dst = "FigureS1_PCA.tiff",            label = "PCA"),
  list(src = "ma_plot_lncrna.pdf",                 dst = "FigureS2_MA_plot.tiff",         label = "MA Plot"),
  list(src = "wgcna_soft_threshold.pdf",           dst = "FigureS3_WGCNA_threshold.tiff", label = "WGCNA Threshold"),
  list(src = "wgcna_sample_dendrogram.pdf",        dst = "FigureS4_WGCNA_sample.tiff",    label = "WGCNA Sample"),
  list(src = "wgcna_module_dendrogram.pdf",        dst = "FigureS5_WGCNA_module.tiff",    label = "WGCNA Module"),
  list(src = "wgcna_module_trait_heatmap.pdf",     dst = "FigureS6_WGCNA_trait.tiff",     label = "WGCNA Trait"),
  list(src = "lasso_cox_cv_plot.pdf",              dst = "FigureS7_LASSO_CV.tiff",        label = "LASSO CV"),
  list(src = "time_dependent_roc.pdf",             dst = "FigureS8_ROC.tiff",             label = "Time-Dependent ROC"),
  list(src = "cgga_km_log2_H19.pdf",              dst = "FigureS9_CGGA_H19.tiff",        label = "CGGA H19 KM"),
  list(src = "cgga_km_log2_CYP1B1_AS1.pdf",      dst = "FigureS10_CGGA_CYP1B1AS1.tiff", label = "CGGA CYP1B1-AS1 KM"),
  list(src = "candidate_prioritization_heatmap.pdf", dst = "FigureS11_ranking.tiff",      label = "Initial Ranking")
)

# ── Function to convert PDF to TIFF using ghostscript ──
pdf_to_tiff <- function(src_path, dst_path, dpi = 600) {
  if (!file.exists(src_path)) {
    cat(sprintf("  [SKIP] Source not found: %s\n", basename(src_path)))
    return(FALSE)
  }

  cmd <- sprintf(
    'gs -dNOPAUSE -dBATCH -sDEVICE=tiff24nc -r%d -sOutputFile="%s" "%s" 2>/dev/null',
    dpi, dst_path, src_path
  )
  ret <- system(cmd)

  if (ret == 0 && file.exists(dst_path)) {
    size_mb <- round(file.info(dst_path)$size / 1024 / 1024, 1)
    cat(sprintf("  [OK] %s (%.1f MB)\n", basename(dst_path), size_mb))
    return(TRUE)
  } else {
    cat(sprintf("  [FAIL] %s\n", basename(dst_path)))
    return(FALSE)
  }
}

# Also copy PDFs (vector format, preferred by many journals)
pdf_copy <- function(src_path, dst_path_tiff) {
  dst_pdf <- sub("\\.tiff$", ".pdf", dst_path_tiff)
  if (file.exists(src_path)) {
    file.copy(src_path, dst_pdf, overwrite = TRUE)
  }
}

# ── Process Main Figures ──
cat("── Main Figures ──\n")
for (fig in main_figures) {
  src_path <- file.path(fig_src, fig$src)
  dst_path <- file.path(fig_dst, fig$dst)
  cat(sprintf("  Converting: %s → %s (%s)\n", fig$src, fig$dst, fig$label))
  pdf_to_tiff(src_path, dst_path, dpi = 600)
  pdf_copy(src_path, dst_path)
}

# ── Process Supplementary Figures ──
cat("\n── Supplementary Figures ──\n")
for (fig in supp_figures) {
  src_path <- file.path(fig_src, fig$src)
  dst_path <- file.path(sup_dst, fig$dst)
  cat(sprintf("  Converting: %s → %s (%s)\n", fig$src, fig$dst, fig$label))
  pdf_to_tiff(src_path, dst_path, dpi = 600)
  pdf_copy(src_path, dst_path)
}

# ── Copy individual KM curves to supplementary ──
cat("\n── Individual KM Curves (Supplementary) ──\n")
km_dir <- file.path(fig_src, "km_curves")
if (dir.exists(km_dir)) {
  km_files <- list.files(km_dir, pattern = "\\.pdf$", full.names = TRUE)
  for (km_f in km_files) {
    dst_name <- sub("^km_", "FigureS_KM_", basename(km_f))
    dst_name <- sub("\\.pdf$", ".tiff", dst_name)
    cat(sprintf("  Converting: %s\n", basename(km_f)))
    pdf_to_tiff(km_f, file.path(sup_dst, dst_name), dpi = 600)
    file.copy(km_f, file.path(sup_dst, sub("\\.tiff$", ".pdf", dst_name)), overwrite = TRUE)
  }
}

# ── Summary ──
cat("\n============================================\n")
main_count <- length(list.files(fig_dst, pattern = "\\.tiff$"))
supp_count <- length(list.files(sup_dst, pattern = "\\.tiff$"))
cat(sprintf(" Main figures: %d TIFF files in manuscript/figures/\n", main_count))
cat(sprintf(" Supplementary: %d TIFF files in manuscript/supplementary/\n", supp_count))
cat(sprintf(" PDF copies also created for vector format submission\n"))
cat("============================================\n")
