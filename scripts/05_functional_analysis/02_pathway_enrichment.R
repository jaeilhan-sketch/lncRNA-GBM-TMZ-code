#!/usr/bin/env Rscript
# ============================================================================
# Phase 5.2: GO & KEGG Pathway Enrichment + GSEA
#
# Performs enrichment analysis on co-expressed mRNAs of each DE-lncRNA.
# Also runs GSEA on the full ranked gene list.
# ============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ggplot2)
  library(data.table)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

cat("============================================\n")
cat(" Phase 5.2: Pathway Enrichment Analysis\n")
cat("============================================\n")

# ── Load data ──
cat("[1/5] Loading data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
coexpr <- fread(file.path(table_dir, "coexpression_lncrna_mrna.tsv"))

# ============================================================================
# 1. ORA: Co-expressed genes of all DE-lncRNAs (pooled)
# ============================================================================
cat("[2/5] Running GO/KEGG ORA on co-expressed genes...\n")

# Unique co-expressed mRNAs
coexpr_genes <- unique(coexpr$mrna_id)
coexpr_genes_noversion <- sub("\\.\\d+$", "", coexpr_genes)

# Convert to Entrez IDs
gene_mapping <- bitr(coexpr_genes_noversion,
                      fromType = "ENSEMBL",
                      toType = c("ENTREZID", "SYMBOL"),
                      OrgDb = org.Hs.eg.db)
cat(sprintf("  Mapped genes: %d / %d\n", nrow(gene_mapping), length(coexpr_genes)))

# ── GO Biological Process ──
ego_bp <- enrichGO(
  gene = gene_mapping$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1,
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
  fwrite(as.data.frame(ego_bp), file.path(table_dir, "enrichment_GO_BP.tsv"), sep = "\t")

  p_go <- dotplot(ego_bp, showCategory = 20, title = "GO Biological Process") +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(fig_dir, "enrichment_GO_BP_dotplot.pdf"), p_go,
         width = 200, height = 250, units = "mm", dpi = 300)

  cat(sprintf("  GO BP significant terms: %d\n", nrow(as.data.frame(ego_bp))))
}

# ── GO Molecular Function ──
ego_mf <- enrichGO(
  gene = gene_mapping$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

if (!is.null(ego_mf) && nrow(as.data.frame(ego_mf)) > 0) {
  fwrite(as.data.frame(ego_mf), file.path(table_dir, "enrichment_GO_MF.tsv"), sep = "\t")
}

# ── KEGG Pathway ──
ekegg <- enrichKEGG(
  gene = gene_mapping$ENTREZID,
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  minGSSize = 10
)

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
  fwrite(as.data.frame(ekegg), file.path(table_dir, "enrichment_KEGG.tsv"), sep = "\t")

  p_kegg <- dotplot(ekegg, showCategory = 20, title = "KEGG Pathway Enrichment") +
    theme(axis.text.y = element_text(size = 9))
  ggsave(file.path(fig_dir, "enrichment_KEGG_dotplot.pdf"), p_kegg,
         width = 200, height = 220, units = "mm", dpi = 300)

  cat(sprintf("  KEGG significant pathways: %d\n", nrow(as.data.frame(ekegg))))
}


# ============================================================================
# 2. GSEA: Full ranked gene list (mRNA)
# ============================================================================
cat("[3/5] Running GSEA...\n")

# Rank genes by DESeq2 stat (mRNA results)
ranked_genes <- res_mrna_df[!is.na(res_mrna_df$stat), ]
ranked_genes$gene_id_noversion <- sub("\\.\\d+$", "", ranked_genes$gene_id)

# Map to Entrez
ranked_mapping <- bitr(ranked_genes$gene_id_noversion,
                        fromType = "ENSEMBL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)
ranked_genes <- merge(ranked_genes, ranked_mapping,
                       by.x = "gene_id_noversion", by.y = "ENSEMBL")

# Create named vector (deduplicate)
gene_list <- ranked_genes$stat
names(gene_list) <- ranked_genes$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]

# ── GSEA GO ──
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  fwrite(as.data.frame(gsea_go), file.path(table_dir, "gsea_GO_BP.tsv"), sep = "\t")

  # Top enrichment plots
  for (i in seq_len(min(5, nrow(as.data.frame(gsea_go))))) {
    p <- gseaplot2(gsea_go, geneSetID = i, title = gsea_go@result$Description[i])
    ggsave(file.path(fig_dir, paste0("gsea_GO_BP_", i, ".pdf")), p,
           width = 180, height = 120, units = "mm", dpi = 300)
  }
  cat(sprintf("  GSEA GO BP significant terms: %d\n", nrow(as.data.frame(gsea_go))))
}

# ── GSEA KEGG ──
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  fwrite(as.data.frame(gsea_kegg), file.path(table_dir, "gsea_KEGG.tsv"), sep = "\t")
  cat(sprintf("  GSEA KEGG significant pathways: %d\n", nrow(as.data.frame(gsea_kegg))))
}


# ============================================================================
# 3. TMZ-specific pathway analysis
# ============================================================================
cat("[4/5] TMZ-specific pathway enrichment...\n")

# Define TMZ-related gene sets
tmz_genesets <- list(
  "DNA_damage_repair" = c("MGMT", "MLH1", "MSH2", "MSH6", "PMS2",
                           "BRCA1", "BRCA2", "ATM", "ATR", "CHEK1", "CHEK2",
                           "RAD51", "XRCC1", "PARP1", "PARP2"),
  "Apoptosis" = c("BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BIM",
                   "CASP3", "CASP8", "CASP9", "APAF1", "CYCS", "TP53"),
  "Drug_efflux" = c("ABCB1", "ABCC1", "ABCC2", "ABCG2", "ABCC3", "ABCC4"),
  "Cell_cycle" = c("CDKN1A", "CDKN2A", "CDK4", "CDK6", "RB1",
                    "CCND1", "CCNE1", "E2F1")
)

# Check enrichment of co-expressed genes in TMZ-related pathways
tmz_enrichment <- list()
for (gs_name in names(tmz_genesets)) {
  gs_symbols <- tmz_genesets[[gs_name]]
  gs_entrez <- bitr(gs_symbols, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)

  # Overlap with co-expressed genes
  overlap <- intersect(gene_mapping$SYMBOL,
                        bitr(gene_mapping$ENTREZID, fromType = "ENTREZID",
                             toType = "SYMBOL", OrgDb = org.Hs.eg.db)$SYMBOL)
  overlap_tmz <- intersect(overlap, gs_symbols)

  tmz_enrichment[[gs_name]] <- data.frame(
    pathway = gs_name,
    total_genes_in_set = length(gs_symbols),
    coexpr_overlap = length(overlap_tmz),
    overlap_genes = paste(overlap_tmz, collapse = ", ")
  )
}

tmz_df <- rbindlist(tmz_enrichment)
fwrite(tmz_df, file.path(table_dir, "tmz_pathway_overlap.tsv"), sep = "\t")


# ============================================================================
# 4. Save enrichment summary
# ============================================================================
cat("[5/5] Saving summary...\n")

save(ego_bp, ego_mf, ekegg, gsea_go, gsea_kegg,
     file = file.path(results_dir, "enrichment_objects.RData"))

cat("\n============================================\n")
cat(" Pathway enrichment analysis complete.\n")
cat("============================================\n")
