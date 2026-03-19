#!/usr/bin/env Rscript
# ============================================================================
# Phase 5.3: WGCNA Co-expression Network Analysis
#
# Builds weighted co-expression network using lncRNA + mRNA jointly.
# Identifies modules correlated with TMZ response.
# ============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

set.seed(42)
allowWGCNAThreads(nThreads = 8)

# Fix known WGCNA/stats cor() conflict in newer R versions
temp_cor <- cor
cor <- WGCNA::cor

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

cat("============================================\n")
cat(" Phase 5.3: WGCNA Network Analysis\n")
cat("============================================\n")

# ── Load data ──
cat("[1/6] Loading data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))

vsd_mrna <- vst(dds_mrna, blind = TRUE)

# Combine top variable lncRNAs + mRNAs
expr_lncrna <- assay(vsd_lncrna)
expr_mrna <- assay(vsd_mrna)

# Select top 5000 most variable mRNAs + all expressed lncRNAs
mrna_var <- apply(expr_mrna, 1, var)
top_mrna <- names(sort(mrna_var, decreasing = TRUE))[1:5000]

expr_combined <- rbind(expr_lncrna, expr_mrna[top_mrna, ])
cat(sprintf("  Combined matrix: %d genes x %d samples\n",
            nrow(expr_combined), ncol(expr_combined)))

# Transpose for WGCNA (samples x genes)
datExpr <- t(expr_combined)

# ── Check for outlier samples ──
cat("[2/6] Checking for outlier samples...\n")
sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(fig_dir, "wgcna_sample_dendrogram.pdf"), width = 12, height = 6)
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub = "", xlab = "", cex.lab = 1.2, cex.axis = 1, cex.main = 1.5)
dev.off()

# Remove outliers (if any — use manual cutoff if needed)
# For now, keep all samples
goodSamples <- rep(TRUE, nrow(datExpr))
datExpr <- datExpr[goodSamples, ]

# ── Soft-thresholding power ──
cat("[3/6] Choosing soft-thresholding power...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                          networkType = "signed", verbose = 0)

# Select power where R² > 0.8
r2_threshold <- 0.8
best_power <- sft$powerEstimate
if (is.na(best_power)) {
  # Fallback: pick first power with R² > threshold
  idx <- which(sft$fitIndices$SFT.R.sq > r2_threshold)[1]
  best_power <- if (!is.na(idx)) powers[idx] else 6
}
cat(sprintf("  Selected soft power: %d\n", best_power))

# Plot scale-free topology fit
pdf(file.path(fig_dir, "wgcna_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R²",
     main = "Scale independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = r2_threshold, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# ── Build network and detect modules ──
cat("[4/6] Building network and detecting modules...\n")
net <- blockwiseModules(
  datExpr,
  power = best_power,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3,
  maxBlockSize = 10000
)

moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)
cat(sprintf("  Modules detected: %d\n", length(unique(moduleLabels)) - 1))

# Module dendrogram
pdf(file.path(fig_dir, "wgcna_module_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# ── Module-trait correlation ──
cat("[5/6] Correlating modules with TMZ response...\n")

# Trait data
sample_ids <- rownames(datExpr)
trait_data <- col_data[match(sample_ids, rownames(col_data)), ]

# Numeric traits
traits <- data.frame(
  tmz_response = as.numeric(trait_data$tmz_response == "Responder"),
  row.names = sample_ids
)
if ("MGMT_status" %in% colnames(trait_data)) {
  traits$MGMT_methylated <- as.numeric(grepl("methyl|yes|positive",
                                              trait_data$MGMT_status, ignore.case = TRUE))
}

# Module eigengenes
MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs)

# Correlation
moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPval <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Heatmap
pdf(file.path(fig_dir, "wgcna_module_trait_heatmap.pdf"), width = 8, height = 10)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPval, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)
dev.off()

# ── Extract lncRNAs in TMZ-response modules ──
cat("[6/6] Extracting lncRNAs in TMZ-response modules...\n")

# Find modules significantly correlated with TMZ response
sig_modules <- names(MEs)[moduleTraitPval[, "tmz_response"] < 0.05]

# Gene-module assignments
gene_module <- data.frame(
  gene_id = colnames(datExpr),
  module_color = moduleColors,
  module_label = moduleLabels,
  stringsAsFactors = FALSE
)

# Annotate lncRNA vs mRNA
gene_module$gene_id_noversion <- sub("\\.\\d+$", "", gene_module$gene_id)
gene_module <- merge(gene_module,
                      gene_annot[, c("gene_id", "gene_name", "gene_type")],
                      by = "gene_id", all.x = TRUE)

# lncRNAs in significant modules (GDC GENCODE v36 uses unified "lncRNA" biotype)
lncrna_biotypes <- c("lncRNA", "lincRNA", "antisense", "sense_intronic",
                      "sense_overlapping", "processed_transcript",
                      "bidirectional_promoter_lncRNA")
lncrna_in_modules <- gene_module[gene_module$gene_type %in% lncrna_biotypes &
                                   paste0("ME", gene_module$module_color) %in% sig_modules, ]

cat(sprintf("  lncRNAs in TMZ-response modules: %d\n", nrow(lncrna_in_modules)))

# Save
fwrite(gene_module, file.path(table_dir, "wgcna_gene_modules.tsv"), sep = "\t")
fwrite(lncrna_in_modules, file.path(table_dir, "wgcna_lncrna_tmz_modules.tsv"), sep = "\t")

# Save WGCNA objects
save(net, MEs, moduleColors, moduleLabels, moduleTraitCor, moduleTraitPval,
     gene_module, lncrna_in_modules,
     file = file.path(results_dir, "wgcna_objects.RData"))

# Restore base cor
cor <- temp_cor

cat("\n============================================\n")
cat(" WGCNA analysis complete.\n")
cat("============================================\n")
