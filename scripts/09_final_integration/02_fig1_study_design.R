#!/usr/bin/env Rscript
# ============================================================================
# Figure 1: Study Design Overview
#
# Professional study design flowchart using ggplot2 + grid
# Output: results/figures/Fig1_study_design.pdf
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(gridExtra)
})

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
fig_dir     <- file.path(project_dir, "results/figures")

cat("============================================\n")
cat(" Generating Figure 1: Study Design\n")
cat("============================================\n")

# ============================================================================
# Helper functions
# ============================================================================

# Draw a rounded rectangle box with text
draw_box <- function(x, y, w, h, label, fill = "#FFFFFF", border = "#333333",
                     text_col = "#000000", font_size = 9, font_face = "plain",
                     line_width = 1.5, radius = unit(3, "mm")) {
  roundrectGrob(
    x = unit(x, "npc"), y = unit(y, "npc"),
    width = unit(w, "npc"), height = unit(h, "npc"),
    r = radius,
    gp = gpar(fill = fill, col = border, lwd = line_width),
    name = paste0("box_", gsub(" ", "_", substr(label, 1, 10)))
  )
}

draw_text <- function(x, y, label, font_size = 9, font_face = "plain",
                      text_col = "#000000", hjust = 0.5, vjust = 0.5,
                      lineheight = 1.1) {
  textGrob(
    label = label,
    x = unit(x, "npc"), y = unit(y, "npc"),
    hjust = hjust, vjust = vjust,
    gp = gpar(fontsize = font_size, fontface = font_face, col = text_col,
              lineheight = lineheight)
  )
}

# Draw arrow between two points
draw_arrow <- function(x0, y0, x1, y1, col = "#555555", lwd = 1.8,
                       arrow_size = unit(2.5, "mm")) {
  linesGrob(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    arrow = arrow(length = arrow_size, type = "closed"),
    gp = gpar(col = col, lwd = lwd, fill = col)
  )
}

# Draw simple line (no arrow)
draw_line <- function(x0, y0, x1, y1, col = "#555555", lwd = 1.5, lty = 1) {
  linesGrob(
    x = unit(c(x0, x1), "npc"),
    y = unit(c(y0, y1), "npc"),
    gp = gpar(col = col, lwd = lwd, lty = lty)
  )
}

# ============================================================================
# Color palette
# ============================================================================
col_header    <- "#2C3E50"   # Dark blue-gray (section headers)
col_data      <- "#3498DB"   # Blue (data acquisition)
col_data_lt   <- "#D6EAF8"   # Light blue
col_cohort    <- "#27AE60"   # Green (cohort)
col_cohort_lt <- "#D5F5E3"   # Light green
col_resp      <- "#E74C3C"   # Red (non-responder)
col_resp_lt   <- "#FADBD8"   # Light red
col_good      <- "#2ECC71"   # Bright green (responder)
col_good_lt   <- "#D4EFDF"   # Light bright green
col_analysis  <- "#8E44AD"   # Purple (analysis)
col_anal_lt   <- "#E8DAEF"   # Light purple
col_result    <- "#E67E22"   # Orange (results)
col_result_lt <- "#FDEBD0"   # Light orange
col_valid     <- "#16A085"   # Teal (validation)
col_valid_lt  <- "#D1F2EB"   # Light teal

# ============================================================================
# Build the figure using grid
# ============================================================================

cairo_pdf(file.path(fig_dir, "Fig1_study_design.pdf"), width = 14, height = 10)

grid.newpage()

# Background
grid.rect(gp = gpar(fill = "#FAFAFA", col = NA))

# ── Title ──
grid.draw(draw_text(0.50, 0.97, "Figure 1. Study Design Overview",
                    font_size = 16, font_face = "bold", text_col = col_header))
grid.draw(draw_text(0.50, 0.945,
                    "lncRNA Signatures Associated with Temozolomide Response in Glioblastoma",
                    font_size = 10, text_col = "#666666"))

# ============================================================================
# ROW 1: Data Source (y ~ 0.87)
# ============================================================================
y1 <- 0.87
# Box: TCGA-GBM
grid.draw(draw_box(0.18, y1, 0.30, 0.055, "", fill = col_data_lt, border = col_data))
grid.draw(draw_text(0.18, y1 + 0.005, "TCGA-GBM (GDC Portal)",
                    font_size = 11, font_face = "bold", text_col = col_data))
grid.draw(draw_text(0.18, y1 - 0.015, "RNA-seq STAR-Counts  |  GENCODE v36",
                    font_size = 8, text_col = "#555555"))

# Arrow
grid.draw(draw_arrow(0.33, y1, 0.39, y1))

# Box: Inclusion/Exclusion
grid.draw(draw_box(0.54, y1, 0.28, 0.055, "", fill = "#FFF9E6", border = "#F1C40F"))
grid.draw(draw_text(0.54, y1 + 0.013, "Inclusion / Exclusion Criteria",
                    font_size = 10, font_face = "bold", text_col = "#7D6608"))
grid.draw(draw_text(0.54, y1 - 0.012,
                    "IDH-wt  |  Primary GBM  |  TMZ-treated",
                    font_size = 8, text_col = "#555555"))

# Arrow
grid.draw(draw_arrow(0.68, y1, 0.74, y1))

# Box: Final Cohort
grid.draw(draw_box(0.85, y1, 0.22, 0.055, "", fill = col_cohort_lt, border = col_cohort))
grid.draw(draw_text(0.85, y1 + 0.008, "Final Cohort",
                    font_size = 11, font_face = "bold", text_col = col_cohort))
grid.draw(draw_text(0.85, y1 - 0.013, "n = 94 patients",
                    font_size = 9, text_col = "#555555"))

# ============================================================================
# ROW 2: Patient Classification (y ~ 0.77)
# ============================================================================
y2 <- 0.77

# Arrow from cohort down
grid.draw(draw_arrow(0.85, y1 - 0.028, 0.85, y2 + 0.033))

# Branching arrows
grid.draw(draw_line(0.85, y2 + 0.033, 0.72, y2 + 0.033, col = "#555555"))
grid.draw(draw_arrow(0.72, y2 + 0.033, 0.72, y2 + 0.005))
grid.draw(draw_line(0.85, y2 + 0.033, 0.96, y2 + 0.033, col = "#555555"))
grid.draw(draw_arrow(0.96, y2 + 0.033, 0.96, y2 + 0.005))

# Box: Responder
grid.draw(draw_box(0.72, y2 - 0.025, 0.18, 0.055, "", fill = col_good_lt, border = col_good))
grid.draw(draw_text(0.72, y2 - 0.010, "Responder (n = 72)",
                    font_size = 10, font_face = "bold", text_col = "#1B7A3D"))
grid.draw(draw_text(0.72, y2 - 0.038, "PFS \u2265 6 months",
                    font_size = 8, text_col = "#555555"))

# Box: Non-responder
grid.draw(draw_box(0.96, y2 - 0.025, 0.18, 0.055, "", fill = col_resp_lt, border = col_resp))
grid.draw(draw_text(0.96, y2 - 0.010, "Non-responder (n = 22)",
                    font_size = 10, font_face = "bold", text_col = "#922B21"))
grid.draw(draw_text(0.96, y2 - 0.038, "PFS < 6 months",
                    font_size = 8, text_col = "#555555"))

# TMZ response label
grid.draw(draw_text(0.85, y2 + 0.05, "TMZ Response Classification",
                    font_size = 9, font_face = "italic", text_col = "#888888"))

# ============================================================================
# ROW 3-6: Analysis Pipeline (left column, y ~ 0.63 to 0.08)
# ============================================================================

# Phase labels on the far left
phase_x <- 0.04

# ── Phase A: Differential Expression ──
ya <- 0.63
grid.draw(draw_box(0.22, ya, 0.36, 0.075, "", fill = col_anal_lt, border = col_analysis,
                   line_width = 2))
grid.draw(draw_text(0.22, ya + 0.018, "Phase 1: Differential Expression Analysis",
                    font_size = 11, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, ya - 0.005, "DESeq2  +  SVA (1 surrogate variable)",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, ya - 0.022, "Design: ~ gender + SV1 + tmz_response",
                    font_size = 8, text_col = "#777777"))

# Result box for Phase A
grid.draw(draw_box(0.22, ya - 0.070, 0.30, 0.045, "", fill = "#FFFFFF", border = col_analysis,
                   line_width = 1))
grid.draw(draw_text(0.22, ya - 0.060, "16 DE-lncRNAs (14 down / 2 up)",
                    font_size = 9, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, ya - 0.078, "|log2FC| > 1.0,  padj < 0.05",
                    font_size = 7.5, text_col = "#888888"))

# Arrow down
grid.draw(draw_arrow(0.22, ya - 0.095, 0.22, ya - 0.12))

# ── Phase B: Functional Analysis ──
yb <- 0.46
grid.draw(draw_box(0.22, yb, 0.36, 0.085, "", fill = col_anal_lt, border = col_analysis,
                   line_width = 2))
grid.draw(draw_text(0.22, yb + 0.025, "Phase 2: Functional Analysis",
                    font_size = 11, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, yb + 0.005,
                    "Co-expression: Spearman |\u03C1| \u2265 0.6, padj < 0.01",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, yb - 0.012, "GO / KEGG (ORA)  |  GSEA  |  WGCNA",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, yb - 0.028, "ceRNA Network (ENCORI / starBase)",
                    font_size = 8.5, text_col = "#555555"))

# Result box
grid.draw(draw_box(0.22, yb - 0.080, 0.34, 0.050, "", fill = "#FFFFFF", border = col_analysis,
                   line_width = 1))
grid.draw(draw_text(0.22, yb - 0.065, "1,619 co-expression pairs  |  114 GO terms",
                    font_size = 8.5, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, yb - 0.080, "109,710 ceRNA triplets  |  11 validated (DNM3OS)",
                    font_size = 8, text_col = "#888888"))
grid.draw(draw_text(0.22, yb - 0.095, "15 KEGG pathways  |  6 WGCNA modules",
                    font_size = 8, text_col = "#888888"))

# Arrow down
grid.draw(draw_arrow(0.22, yb - 0.108, 0.22, yb - 0.135))

# ── Phase C: Survival Analysis ──
yc <- 0.27
grid.draw(draw_box(0.22, yc, 0.36, 0.085, "", fill = col_anal_lt, border = col_analysis,
                   line_width = 2))
grid.draw(draw_text(0.22, yc + 0.025, "Phase 3: Survival & Prognostic Analysis",
                    font_size = 11, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, yc + 0.005, "Kaplan-Meier (log-rank)  |  optimal cutpoint",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, yc - 0.012, "Multivariate Cox (age + gender + lncRNA)",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, yc - 0.028, "LASSO-Cox  |  Bootstrap (1,000 iter)",
                    font_size = 8.5, text_col = "#555555"))

# Result box
grid.draw(draw_box(0.22, yc - 0.077, 0.34, 0.045, "", fill = "#FFFFFF", border = col_analysis,
                   line_width = 1))
grid.draw(draw_text(0.22, yc - 0.065, "H19: HR=1.14, p=0.007  |  Bootstrap 74.1%",
                    font_size = 8.5, font_face = "bold", text_col = col_analysis))
grid.draw(draw_text(0.22, yc - 0.082, "LINC01936: HR=1.48, p=0.012  |  Bootstrap 61.0%",
                    font_size = 8, text_col = "#888888"))

# Arrow down
grid.draw(draw_arrow(0.22, yc - 0.102, 0.22, yc - 0.125))

# ── Phase D: External Validation ──
yd <- 0.095
grid.draw(draw_box(0.22, yd, 0.36, 0.065, "", fill = col_valid_lt, border = col_valid,
                   line_width = 2))
grid.draw(draw_text(0.22, yd + 0.015, "Phase 4: External Validation (CGGA)",
                    font_size = 11, font_face = "bold", text_col = col_valid))
grid.draw(draw_text(0.22, yd - 0.005, "CGGA mRNAseq_693 (n = 92, TMZ-treated GBM)",
                    font_size = 8.5, text_col = "#555555"))
grid.draw(draw_text(0.22, yd - 0.020,
                    "H19: HR=1.06, p=0.094 (consistent direction)",
                    font_size = 8.5, text_col = "#555555"))

# ============================================================================
# RIGHT PANEL: Key Findings Summary
# ============================================================================

# Background box for right panel
grid.draw(roundrectGrob(
  x = unit(0.72, "npc"), y = unit(0.37, "npc"),
  width = unit(0.48, "npc"), height = unit(0.60, "npc"),
  r = unit(4, "mm"),
  gp = gpar(fill = col_result_lt, col = col_result, lwd = 2.5)
))

grid.draw(draw_text(0.72, 0.645, "Key Findings",
                    font_size = 14, font_face = "bold", text_col = col_result))

# Divider line
grid.draw(draw_line(0.52, 0.625, 0.92, 0.625, col = col_result, lwd = 1))

# ── Finding 1: H19 ──
yf1 <- 0.59
grid.draw(draw_box(0.72, yf1, 0.42, 0.055, "", fill = "#FFFFFF", border = "#E67E22",
                   line_width = 1.5))
grid.draw(draw_text(0.545, yf1 + 0.007, "\u2460", font_size = 12, text_col = col_result,
                    hjust = 0))
grid.draw(draw_text(0.57, yf1 + 0.007, "H19 — Primary Candidate (Score: 7.43)",
                    font_size = 10, font_face = "bold", text_col = "#333333", hjust = 0))
grid.draw(draw_text(0.57, yf1 - 0.014, "Independent prognostic factor  |  CGGA validated",
                    font_size = 8, text_col = "#666666", hjust = 0))

# ── Finding 2: LINC01936 ──
yf2 <- 0.52
grid.draw(draw_box(0.72, yf2, 0.42, 0.055, "", fill = "#FFFFFF", border = "#E67E22",
                   line_width = 1.5))
grid.draw(draw_text(0.545, yf2 + 0.007, "\u2461", font_size = 12, text_col = col_result,
                    hjust = 0))
grid.draw(draw_text(0.57, yf2 + 0.007, "LINC01936 — Novel in GBM (Score: 4.78)",
                    font_size = 10, font_face = "bold", text_col = "#333333", hjust = 0))
grid.draw(draw_text(0.57, yf2 - 0.014, "First report in GBM  |  Highest Cox HR = 1.48",
                    font_size = 8, text_col = "#666666", hjust = 0))

# ── Finding 3: DNM3OS ceRNA ──
yf3 <- 0.45
grid.draw(draw_box(0.72, yf3, 0.42, 0.055, "", fill = "#FFFFFF", border = "#E67E22",
                   line_width = 1.5))
grid.draw(draw_text(0.545, yf3 + 0.007, "\u2462", font_size = 12, text_col = col_result,
                    hjust = 0))
grid.draw(draw_text(0.57, yf3 + 0.007, "DNM3OS ceRNA Axis (Score: 5.64)",
                    font_size = 10, font_face = "bold", text_col = "#333333", hjust = 0))
grid.draw(draw_text(0.57, yf3 - 0.014, "11 validated triplets  |  6 miRNAs \u2192 5 target mRNAs",
                    font_size = 8, text_col = "#666666", hjust = 0))

# ── Summary Statistics ──
ys <- 0.35
grid.draw(draw_line(0.52, ys + 0.035, 0.92, ys + 0.035, col = "#E67E22", lwd = 0.8,
                    lty = 2))
grid.draw(draw_text(0.72, ys + 0.015, "Summary Statistics",
                    font_size = 10, font_face = "bold.italic", text_col = "#7D5A29"))

# Stats grid
grid.draw(draw_text(0.56, ys - 0.015, "94", font_size = 16, font_face = "bold",
                    text_col = col_data))
grid.draw(draw_text(0.56, ys - 0.038, "Patients", font_size = 7.5, text_col = "#666666"))

grid.draw(draw_text(0.66, ys - 0.015, "16", font_size = 16, font_face = "bold",
                    text_col = col_analysis))
grid.draw(draw_text(0.66, ys - 0.038, "DE-lncRNAs", font_size = 7.5, text_col = "#666666"))

grid.draw(draw_text(0.76, ys - 0.015, "1,619", font_size = 16, font_face = "bold",
                    text_col = col_analysis))
grid.draw(draw_text(0.76, ys - 0.038, "Co-expr pairs", font_size = 7.5, text_col = "#666666"))

grid.draw(draw_text(0.87, ys - 0.015, "2", font_size = 16, font_face = "bold",
                    text_col = col_result))
grid.draw(draw_text(0.87, ys - 0.038, "Prognostic", font_size = 7.5, text_col = "#666666"))

# ── Pathway highlights ──
yp <- 0.24
grid.draw(draw_line(0.52, yp + 0.035, 0.92, yp + 0.035, col = "#E67E22", lwd = 0.8,
                    lty = 2))
grid.draw(draw_text(0.72, yp + 0.015, "Key Pathways",
                    font_size = 10, font_face = "bold.italic", text_col = "#7D5A29"))

pathways <- c(
  "\u2022 Ubiquitin-mediated proteolysis",
  "\u2022 Endocytosis / Autophagy / Lysosome",
  "\u2022 Cellular senescence",
  "\u2022 RAS/MAPK signaling (via RASA2)",
  "\u2022 Glutamine metabolism (via GLS)"
)
for (i in seq_along(pathways)) {
  grid.draw(draw_text(0.56, yp - 0.005 - (i-1)*0.020, pathways[i],
                      font_size = 8.5, text_col = "#555555", hjust = 0))
}

# ── Figures reference ──
yref <- 0.095
grid.draw(draw_line(0.52, yref + 0.035, 0.92, yref + 0.035, col = "#E67E22", lwd = 0.8,
                    lty = 2))
grid.draw(draw_text(0.72, yref + 0.015, "Associated Figures",
                    font_size = 10, font_face = "bold.italic", text_col = "#7D5A29"))

figs <- c(
  "Fig 2: Volcano plot & Heatmap (DE-lncRNAs)",
  "Fig 3: Top candidate boxplot & KM curve",
  "Fig 4: GO / KEGG / GSEA functional enrichment",
  "Fig 5: ceRNA network (DNM3OS axis)",
  "Fig 6: Forest plot (TCGA vs CGGA validation)",
  "Fig 7: Multi-evidence candidate heatmap"
)
for (i in seq_along(figs)) {
  grid.draw(draw_text(0.56, yref - 0.005 - (i-1)*0.017, figs[i],
                      font_size = 7.5, text_col = "#666666", hjust = 0))
}

# ============================================================================
# Connecting arrows from left pipeline to right panel
# ============================================================================

# Dashed arrow from Phase A result to Key Findings
grid.draw(draw_line(0.40, ya - 0.070, 0.50, 0.59, col = "#CCCCCC", lwd = 1, lty = 3))

# Dashed arrow from Phase B to Key Findings
grid.draw(draw_line(0.40, yb - 0.080, 0.50, 0.45, col = "#CCCCCC", lwd = 1, lty = 3))

# Dashed arrow from Phase C to Key Findings
grid.draw(draw_line(0.40, yc - 0.077, 0.50, 0.52, col = "#CCCCCC", lwd = 1, lty = 3))

# ============================================================================
# Phase labels (left side)
# ============================================================================
grid.draw(draw_text(phase_x, ya, "A", font_size = 20, font_face = "bold",
                    text_col = "#CCCCCC"))
grid.draw(draw_text(phase_x, yb, "B", font_size = 20, font_face = "bold",
                    text_col = "#CCCCCC"))
grid.draw(draw_text(phase_x, yc, "C", font_size = 20, font_face = "bold",
                    text_col = "#CCCCCC"))
grid.draw(draw_text(phase_x, yd, "D", font_size = 20, font_face = "bold",
                    text_col = "#CCCCCC"))

# ============================================================================
# Arrow from Row 1 cohort to Pipeline
# ============================================================================
# Vertical arrow from cohort box down to Phase A
grid.draw(draw_line(0.85, y2 - 0.055, 0.50, y2 - 0.055, col = "#AAAAAA", lwd = 1, lty = 2))
grid.draw(draw_line(0.50, y2 - 0.055, 0.22, y2 - 0.055, col = "#AAAAAA", lwd = 1, lty = 2))
grid.draw(draw_arrow(0.22, y2 - 0.055, 0.22, ya + 0.040, col = "#555555"))

dev.off()

cat(sprintf("\n  Figure saved: %s\n", file.path(fig_dir, "Fig1_study_design.pdf")))
cat("============================================\n")
cat(" Figure 1 generation complete.\n")
cat("============================================\n")
