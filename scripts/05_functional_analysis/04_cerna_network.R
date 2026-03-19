#!/usr/bin/env Rscript
# ============================================================================
# Phase 5.4: ceRNA (competing endogenous RNA) Network Analysis
#
# Constructs lncRNA-miRNA-mRNA regulatory network based on shared MREs.
# Uses starBase/ENCORI and TargetScan predictions.
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

cat("============================================\n")
cat(" Phase 5.4: ceRNA Network Analysis\n")
cat("============================================\n")

# ── Load DE-lncRNA results ──
cat("[1/5] Loading DE-lncRNA data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
coexpr <- fread(file.path(table_dir, "coexpression_lncrna_mrna.tsv"))

sig_lncrna_names <- unique(sig_lncrna$gene_name[!is.na(sig_lncrna$gene_name)])
cat(sprintf("  DE-lncRNAs with gene names: %d\n", length(sig_lncrna_names)))

# ============================================================================
# 2. Query miRNA-lncRNA interactions (ENCORI/starBase)
# ============================================================================
cat("[2/5] Querying miRNA-lncRNA interactions...\n")
cat("  NOTE: This step requires pre-downloaded interaction databases.\n")
cat("  Download from:\n")
cat("    - ENCORI: https://rnasysu.com/encori/\n")
cat("    - miRcode: http://www.mircode.org/\n")
cat("    - LncBase v3: https://diana.e-ce.uth.gr/lncbasev3\n")
cat("\n")

# Expected input files (user needs to download these)
mirna_lncrna_file <- file.path(project_dir, "data/reference/mirna_lncrna_interactions.tsv")
mirna_mrna_file <- file.path(project_dir, "data/reference/mirna_mrna_interactions.tsv")

# Template for the expected format
template_lncrna <- data.table(
  miRNA = character(),
  lncRNA = character(),
  lncRNA_ensembl = character(),
  interaction_type = character(),
  source_db = character(),
  clip_support = logical()
)

template_mrna <- data.table(
  miRNA = character(),
  mRNA = character(),
  mRNA_ensembl = character(),
  interaction_type = character(),
  source_db = character(),
  target_score = numeric()
)

if (!file.exists(mirna_lncrna_file)) {
  cat("  Creating template file for miRNA-lncRNA interactions...\n")
  fwrite(template_lncrna, mirna_lncrna_file, sep = "\t")
  cat(sprintf("  Template: %s\n", mirna_lncrna_file))
  cat("  Please populate with data from ENCORI/miRcode.\n")
}

if (!file.exists(mirna_mrna_file)) {
  cat("  Creating template file for miRNA-mRNA interactions...\n")
  fwrite(template_mrna, mirna_mrna_file, sep = "\t")
  cat(sprintf("  Template: %s\n", mirna_mrna_file))
  cat("  Please populate with data from TargetScan/miRDB.\n")
}

# ============================================================================
# 3. Build ceRNA network (if interaction data available)
# ============================================================================
cat("[3/5] Building ceRNA network...\n")

if (file.exists(mirna_lncrna_file) && file.info(mirna_lncrna_file)$size > 100 &&
    file.exists(mirna_mrna_file) && file.info(mirna_mrna_file)$size > 100) {

  mirna_lncrna <- fread(mirna_lncrna_file)
  mirna_mrna <- fread(mirna_mrna_file)

  # Filter for DE-lncRNAs
  mirna_lncrna_de <- mirna_lncrna[lncRNA %in% sig_lncrna_names]

  # Find shared miRNAs between lncRNA and mRNA
  shared_mirna <- intersect(mirna_lncrna_de$miRNA, mirna_mrna$miRNA)
  cat(sprintf("  Shared miRNAs: %d\n", length(shared_mirna)))

  # Build ceRNA triplets: lncRNA - miRNA - mRNA
  cerna_triplets <- merge(
    mirna_lncrna_de[miRNA %in% shared_mirna, .(miRNA, lncRNA)],
    mirna_mrna[miRNA %in% shared_mirna, .(miRNA, mRNA)],
    by = "miRNA",
    allow.cartesian = TRUE
  )

  # Filter: require co-expression support (positive correlation)
  # lncRNA↑ → miRNA sponged → mRNA↑ (same direction)
  coexpr_pairs <- paste(coexpr$lncrna_name, coexpr$mrna_name, sep = "_")
  cerna_triplets$coexpr_support <- paste(cerna_triplets$lncRNA,
                                          cerna_triplets$mRNA, sep = "_") %in% coexpr_pairs

  cerna_supported <- cerna_triplets[coexpr_support == TRUE]
  cat(sprintf("  ceRNA triplets (co-expression supported): %d\n", nrow(cerna_supported)))

  # Save
  fwrite(cerna_triplets, file.path(table_dir, "cerna_triplets_all.tsv"), sep = "\t")
  fwrite(cerna_supported, file.path(table_dir, "cerna_triplets_supported.tsv"), sep = "\t")

} else {
  cat("  Skipping ceRNA network (interaction databases not yet populated).\n")
  cat("  To complete this step:\n")
  cat("    1. Download miRNA-lncRNA interactions from ENCORI\n")
  cat("    2. Download miRNA-mRNA interactions from TargetScan\n")
  cat("    3. Place in data/reference/ and re-run this script\n")
}

# ============================================================================
# 4. Generate network files for Cytoscape
# ============================================================================
cat("[4/5] Generating Cytoscape-compatible network files...\n")

# From co-expression data: lncRNA → mRNA edges
network_edges <- data.table(
  source = coexpr$lncrna_name,
  target = coexpr$mrna_name,
  interaction = ifelse(coexpr$correlation > 0, "positive_coexpr", "negative_coexpr"),
  weight = abs(coexpr$correlation),
  edge_type = "coexpression"
)

# Node attributes
lncrna_nodes <- unique(data.table(
  id = coexpr$lncrna_name,
  node_type = "lncRNA",
  DE_status = "significant"
))

mrna_nodes <- unique(data.table(
  id = coexpr$mrna_name,
  node_type = "mRNA",
  DE_status = "co-expressed"
))

network_nodes <- rbind(lncrna_nodes, mrna_nodes)

fwrite(network_edges, file.path(table_dir, "cytoscape_edges.tsv"), sep = "\t")
fwrite(network_nodes, file.path(table_dir, "cytoscape_nodes.tsv"), sep = "\t")

cat(sprintf("  Network: %d nodes, %d edges\n",
            nrow(network_nodes), nrow(network_edges)))

# ============================================================================
# 5. Summary
# ============================================================================
cat("[5/5] Summary...\n")

cat("\n============================================\n")
cat(" ceRNA network analysis complete.\n")
cat(" Cytoscape files:\n")
cat(sprintf("   Edges: %s\n", file.path(table_dir, "cytoscape_edges.tsv")))
cat(sprintf("   Nodes: %s\n", file.path(table_dir, "cytoscape_nodes.tsv")))
cat("============================================\n")
