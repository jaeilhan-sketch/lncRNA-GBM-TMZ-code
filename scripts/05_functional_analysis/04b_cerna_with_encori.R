#!/usr/bin/env Rscript
# ============================================================================
# Phase 5.4b: ceRNA Network with ENCORI Data
#
# Uses downloaded ENCORI miRNA-lncRNA/mRNA interactions to build
# ceRNA (lncRNA-miRNA-mRNA) regulatory network.
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
ref_dir     <- file.path(project_dir, "data/reference")

cat("============================================\n")
cat(" Phase 5.4b: ceRNA Network (ENCORI-based)\n")
cat("============================================\n")

# ── Load data ──
cat("[1/5] Loading data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
coexpr <- fread(file.path(table_dir, "coexpression_lncrna_mrna.tsv"))

# ── Load ENCORI data (with fill=TRUE for inconsistent column counts) ──
mirna_lncrna <- fread(file.path(ref_dir, "mirna_lncrna_interactions.tsv"),
                       fill = TRUE, header = TRUE)
mirna_mrna <- fread(file.path(ref_dir, "mirna_mrna_interactions.tsv"),
                     fill = TRUE, header = TRUE)
cerna_encori <- fread(file.path(ref_dir, "encori_cerna_pairs.tsv"),
                       fill = TRUE, header = TRUE)

cat(sprintf("  miRNA-lncRNA interactions: %d\n", nrow(mirna_lncrna)))
cat(sprintf("  miRNA-mRNA interactions: %d\n", nrow(mirna_mrna)))
cat(sprintf("  ENCORI ceRNA pairs: %d\n", nrow(cerna_encori)))
cat(sprintf("  miRNA-lncRNA columns: %s\n", paste(colnames(mirna_lncrna), collapse=", ")))

# ── Standardize column names ──
# Find the right columns by checking for known patterns
find_col <- function(dt, patterns) {
  for (p in patterns) {
    matches <- grep(p, colnames(dt), value = TRUE, ignore.case = TRUE)
    if (length(matches) > 0) return(matches[1])
  }
  return(NULL)
}

mirna_col <- find_col(mirna_lncrna, c("^miRNAname$", "miRNAname"))
lncrna_col <- find_col(mirna_lncrna, c("^geneName$", "geneName"))
mrna_name_col <- find_col(mirna_mrna, c("^geneName$", "geneName"))

# If column detection failed, try by position
if (is.null(mirna_col)) mirna_col <- colnames(mirna_lncrna)[2]
if (is.null(lncrna_col)) lncrna_col <- colnames(mirna_lncrna)[4]
if (is.null(mrna_name_col)) mrna_name_col <- colnames(mirna_mrna)[4]

cat(sprintf("  miRNA col: %s, lncRNA col: %s, mRNA col: %s\n",
            mirna_col, lncrna_col, mrna_name_col))

# DE-lncRNA names
sig_lncrna_names <- unique(sig_lncrna$gene_name[!is.na(sig_lncrna$gene_name)])
cat(sprintf("  DE-lncRNAs: %d\n", length(sig_lncrna_names)))

# ============================================================================
# 2. Build ceRNA triplets from miRNA interactions
# ============================================================================
cat("\n[2/5] Building ceRNA triplets...\n")

# Filter miRNA-lncRNA for DE-lncRNAs
mirna_lnc_de <- mirna_lncrna[get(lncrna_col) %in% sig_lncrna_names]
cat(sprintf("  miRNA-lncRNA pairs (DE-lncRNAs): %d\n", nrow(mirna_lnc_de)))

# Get miRNAs that target our DE-lncRNAs
shared_mirnas <- unique(mirna_lnc_de[[mirna_col]])
cat(sprintf("  miRNAs targeting DE-lncRNAs: %d\n", length(shared_mirnas)))

# Filter miRNA-mRNA for these shared miRNAs
mirna_col_mrna <- find_col(mirna_mrna, c("^miRNAname$", "miRNAname"))
if (is.null(mirna_col_mrna)) mirna_col_mrna <- colnames(mirna_mrna)[2]
mirna_mrna_shared <- mirna_mrna[get(mirna_col_mrna) %in% shared_mirnas]
cat(sprintf("  miRNA-mRNA pairs (shared miRNAs): %d\n", nrow(mirna_mrna_shared)))

# Build triplets: lncRNA - miRNA - mRNA
triplets <- merge(
  mirna_lnc_de[, .(miRNA = get(mirna_col), lncRNA = get(lncrna_col))],
  mirna_mrna_shared[, .(miRNA = get(mirna_col_mrna), mRNA = get(mrna_name_col))],
  by = "miRNA",
  allow.cartesian = TRUE
)

# Remove duplicates
triplets <- unique(triplets)
cat(sprintf("  Total ceRNA triplets: %d\n", nrow(triplets)))
cat(sprintf("  Unique lncRNAs in network: %d\n", length(unique(triplets$lncRNA))))
cat(sprintf("  Unique miRNAs in network: %d\n", length(unique(triplets$miRNA))))
cat(sprintf("  Unique mRNAs in network: %d\n", length(unique(triplets$mRNA))))


# ============================================================================
# 3. Filter by co-expression support
# ============================================================================
cat("\n[3/5] Filtering by co-expression support...\n")

# Create key for co-expression pairs
coexpr_keys <- paste(coexpr$lncrna_name, coexpr$mrna_name, sep = "|")

# Check co-expression support for triplets
triplets$coexpr_key <- paste(triplets$lncRNA, triplets$mRNA, sep = "|")
triplets$coexpr_support <- triplets$coexpr_key %in% coexpr_keys

supported <- triplets[coexpr_support == TRUE]
cat(sprintf("  Co-expression supported triplets: %d\n", nrow(supported)))

if (nrow(supported) > 0) {
  cat("  Supported lncRNAs:\n")
  for (lnc in unique(supported$lncRNA)) {
    n_mirna <- length(unique(supported[lncRNA == lnc, miRNA]))
    n_mrna <- length(unique(supported[lncRNA == lnc, mRNA]))
    cat(sprintf("    %s: %d miRNAs, %d mRNAs\n", lnc, n_mirna, n_mrna))
  }
}


# ============================================================================
# 4. Add ENCORI ceRNA pairs as additional evidence
# ============================================================================
cat("\n[4/5] Integrating ENCORI ceRNA pairs...\n")

# ENCORI ceRNA pairs have statistical support
cerna_lnc_col <- find_col(cerna_encori, c("^ceRNAname$", "ceRNAname"))
cerna_mrna_col <- find_col(cerna_encori, c("^geneName$", "geneName"))
if (is.null(cerna_lnc_col)) cerna_lnc_col <- colnames(cerna_encori)[5]
if (is.null(cerna_mrna_col)) cerna_mrna_col <- colnames(cerna_encori)[2]

if (!is.null(cerna_lnc_col) && !is.null(cerna_mrna_col)) {
  encori_pairs <- cerna_encori[get(cerna_lnc_col) %in% sig_lncrna_names,
                                .(lncRNA = get(cerna_lnc_col),
                                  mRNA = get(cerna_mrna_col))]
  encori_pairs <- unique(encori_pairs)
  cat(sprintf("  ENCORI validated ceRNA pairs: %d\n", nrow(encori_pairs)))
}


# ============================================================================
# 5. Save results and generate network files
# ============================================================================
cat("\n[5/5] Saving ceRNA network...\n")

# Save all triplets
fwrite(triplets, file.path(table_dir, "cerna_triplets_all.tsv"), sep = "\t")

# Save supported triplets
if (nrow(supported) > 0) {
  fwrite(supported, file.path(table_dir, "cerna_triplets_supported.tsv"), sep = "\t")
}

# Generate enhanced Cytoscape network files
# Edges: lncRNA-miRNA, miRNA-mRNA
edges_lnc_mirna <- unique(triplets[, .(source = lncRNA, target = miRNA,
                                         interaction = "sponge",
                                         edge_type = "lncRNA-miRNA")])
edges_mirna_mrna <- unique(triplets[, .(source = miRNA, target = mRNA,
                                          interaction = "target",
                                          edge_type = "miRNA-mRNA")])

# Also add co-expression edges
edges_coexpr <- data.table(
  source = coexpr$lncrna_name,
  target = coexpr$mrna_name,
  interaction = ifelse(coexpr$correlation > 0, "positive_coexpr", "negative_coexpr"),
  edge_type = "coexpression"
)

all_edges <- rbind(edges_lnc_mirna, edges_mirna_mrna, edges_coexpr)

# Node attributes
lncrna_nodes <- unique(data.table(
  id = unique(c(triplets$lncRNA, coexpr$lncrna_name)),
  node_type = "lncRNA"
))
mirna_nodes <- unique(data.table(
  id = unique(triplets$miRNA),
  node_type = "miRNA"
))
mrna_nodes <- unique(data.table(
  id = unique(c(triplets$mRNA, coexpr$mrna_name)),
  node_type = "mRNA"
))
all_nodes <- rbind(lncrna_nodes, mirna_nodes, mrna_nodes)
all_nodes <- all_nodes[!duplicated(id)]

fwrite(all_edges, file.path(table_dir, "cytoscape_cerna_edges.tsv"), sep = "\t")
fwrite(all_nodes, file.path(table_dir, "cytoscape_cerna_nodes.tsv"), sep = "\t")

# Summary statistics
cat(sprintf("\n  ceRNA Network Summary:\n"))
cat(sprintf("  ├─ Total triplets: %d\n", nrow(triplets)))
cat(sprintf("  ├─ Co-expression supported: %d\n", nrow(supported)))
cat(sprintf("  ├─ Network nodes: %d\n", nrow(all_nodes)))
cat(sprintf("  │  ├─ lncRNA: %d\n", nrow(lncrna_nodes)))
cat(sprintf("  │  ├─ miRNA: %d\n", nrow(mirna_nodes)))
cat(sprintf("  │  └─ mRNA: %d\n", nrow(mrna_nodes)))
cat(sprintf("  └─ Network edges: %d\n", nrow(all_edges)))

# Save summary
cerna_summary <- data.frame(
  Metric = c("Total_triplets", "Supported_triplets", "Unique_lncRNA",
             "Unique_miRNA", "Unique_mRNA", "Network_nodes", "Network_edges"),
  Value = c(nrow(triplets), nrow(supported), length(unique(triplets$lncRNA)),
            length(unique(triplets$miRNA)), length(unique(triplets$mRNA)),
            nrow(all_nodes), nrow(all_edges))
)
fwrite(cerna_summary, file.path(table_dir, "cerna_network_summary.tsv"), sep = "\t")

cat("\n============================================\n")
cat(" ceRNA network analysis complete.\n")
cat("============================================\n")
