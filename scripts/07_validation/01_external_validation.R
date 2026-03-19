#!/usr/bin/env Rscript
# ============================================================================
# Phase 7: Validation & Candidate Prioritization
#
# - Internal cross-validation (bootstrap survival analysis)
# - Candidate prioritization scoring
# - Comprehensive summary of all findings
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(DESeq2)
  library(ggplot2)
  library(glmnet)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
valid_dir   <- file.path(results_dir, "validation")
dir.create(valid_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Phase 7: Validation & Candidate Prioritization\n")
cat("============================================\n")

# ── Load key results ──
cat("[1/4] Loading key findings...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
load(file.path(results_dir, "survival_objects.RData"))
coexpr <- fread(file.path(table_dir, "coexpression_lncrna_mrna.tsv"))

# Prepare expression data
expr_lncrna <- assay(vsd_lncrna)
surv_data <- as.data.frame(col_data)
surv_data$os_time <- as.numeric(surv_data$os_days)
surv_data$os_status <- as.numeric(surv_data$os_status)
valid_idx <- !is.na(surv_data$os_time) & !is.na(surv_data$os_status) & surv_data$os_time > 0
surv_data <- surv_data[valid_idx, ]
expr_lncrna <- expr_lncrna[, rownames(surv_data)]


# ============================================================================
# 2. Bootstrap internal validation of survival associations
# ============================================================================
cat("[2/4] Bootstrap validation of survival results...\n")

n_boot <- 1000
sig_gene_ids <- sig_lncrna$gene_id[sig_lncrna$gene_id %in% rownames(expr_lncrna)]

boot_results <- list()
for (gene_id in sig_gene_ids) {
  gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene_id][1]
  expr_vals <- as.numeric(expr_lncrna[gene_id, rownames(surv_data)])

  boot_pvals <- numeric(n_boot)
  boot_hr <- numeric(n_boot)

  for (b in seq_len(n_boot)) {
    idx <- sample(nrow(surv_data), replace = TRUE)
    boot_dat <- surv_data[idx, ]
    boot_dat$expr <- expr_vals[idx]

    tryCatch({
      cox_fit <- coxph(Surv(os_time, os_status) ~ expr, data = boot_dat)
      s <- summary(cox_fit)
      boot_pvals[b] <- s$coefficients[, "Pr(>|z|)"]
      boot_hr[b] <- s$coefficients[, "exp(coef)"]
    }, error = function(e) {
      boot_pvals[b] <<- NA
      boot_hr[b] <<- NA
    })
  }

  boot_sig_rate <- mean(boot_pvals < 0.05, na.rm = TRUE)
  boot_results[[gene_id]] <- data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    boot_sig_rate = boot_sig_rate,
    median_HR = median(boot_hr, na.rm = TRUE),
    HR_95_lower = quantile(boot_hr, 0.025, na.rm = TRUE),
    HR_95_upper = quantile(boot_hr, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

boot_df <- rbindlist(boot_results)
boot_df <- boot_df[order(-boot_df$boot_sig_rate), ]
fwrite(boot_df, file.path(valid_dir, "bootstrap_validation.tsv"), sep = "\t")

cat("  Bootstrap survival validation (1000 iterations):\n")
for (i in seq_len(nrow(boot_df))) {
  cat(sprintf("    %s: sig_rate=%.1f%%, HR=%.2f [%.2f-%.2f]\n",
              boot_df$gene_name[i], boot_df$boot_sig_rate[i] * 100,
              boot_df$median_HR[i], boot_df$HR_95_lower[i], boot_df$HR_95_upper[i]))
}


# ============================================================================
# 3. Comprehensive candidate prioritization
# ============================================================================
cat("\n[3/4] Comprehensive candidate prioritization...\n")

# Merge all evidence layers
priority <- sig_lncrna[, c("gene_id", "gene_name", "baseMean", "log2FoldChange", "padj")]
colnames(priority)[4:5] <- c("log2FC", "DE_padj")

# Add KM survival results
priority <- merge(priority, km_df[, c("gene_id", "km_pvalue")],
                   by = "gene_id", all.x = TRUE)

# Add univariate Cox results
priority <- merge(priority, cox_uni_df[, c("gene_id", "cox_pvalue", "HR")],
                   by = "gene_id", all.x = TRUE)

# Add bootstrap validation
priority <- merge(priority, boot_df[, c("gene_id", "boot_sig_rate", "median_HR")],
                   by = "gene_id", all.x = TRUE)

# Add co-expression info (number of co-expressed mRNAs)
coexpr_count <- coexpr[, .N, by = lncrna_id]
setnames(coexpr_count, c("lncrna_id", "N"), c("gene_id", "n_coexpr_mrna"))
priority <- merge(priority, coexpr_count, by = "gene_id", all.x = TRUE)
priority$n_coexpr_mrna[is.na(priority$n_coexpr_mrna)] <- 0

# Literature support (known lncRNAs in GBM/cancer)
known_lncrna <- c("H19", "MALAT1", "HOTAIR", "MEG3", "NEAT1", "GAS5",
                   "CRNDE", "HULC", "TUG1", "UCA1", "SNHG1", "DNM3OS",
                   "LINC00707", "ZFPM2-AS1", "LINC01936")
priority$literature_support <- priority$gene_name %in% known_lncrna

# Compute priority score
priority$score_DE <- -log10(priority$DE_padj) / max(-log10(priority$DE_padj), na.rm = TRUE) * 3
priority$score_survival <- ifelse(!is.na(priority$km_pvalue),
                                    -log10(priority$km_pvalue) / 3, 0)
priority$score_effect <- pmin(abs(priority$log2FC) / 3, 1) * 2
priority$score_coexpr <- pmin(priority$n_coexpr_mrna / 200, 1) * 2
priority$score_bootstrap <- ifelse(!is.na(priority$boot_sig_rate),
                                     priority$boot_sig_rate * 2, 0)
priority$score_literature <- as.numeric(priority$literature_support) * 1

priority$total_score <- priority$score_DE + priority$score_survival +
                          priority$score_effect + priority$score_coexpr +
                          priority$score_bootstrap + priority$score_literature

priority <- priority[order(-priority$total_score), ]

# Save full ranking
fwrite(priority, file.path(valid_dir, "candidate_priority_ranking.tsv"), sep = "\t")

cat("\n  Candidate Priority Ranking:\n")
cat("  ══════════════════════════════════════════════════\n")
for (i in seq_len(nrow(priority))) {
  cat(sprintf("  %2d. %-15s | log2FC=%+5.2f | DE_p=%.1e | KM_p=%s | boot=%.0f%% | coexpr=%d | score=%.1f%s\n",
              i,
              priority$gene_name[i],
              priority$log2FC[i],
              priority$DE_padj[i],
              ifelse(is.na(priority$km_pvalue[i]), "  NA  ",
                     sprintf("%.3f", priority$km_pvalue[i])),
              ifelse(is.na(priority$boot_sig_rate[i]), 0,
                     priority$boot_sig_rate[i] * 100),
              priority$n_coexpr_mrna[i],
              priority$total_score[i],
              ifelse(priority$literature_support[i], " *", "")))
}
cat("  ══════════════════════════════════════════════════\n")
cat("  * = literature support\n")


# ============================================================================
# 4. Generate validation summary figure
# ============================================================================
cat("\n[4/4] Generating validation summary figure...\n")

# Heatmap-style summary of all evidence
plot_data <- as.data.table(priority[, c("gene_name", "score_DE", "score_survival",
                            "score_effect", "score_coexpr", "score_bootstrap",
                            "score_literature")])
plot_data_long <- melt(plot_data, id.vars = "gene_name",
                        variable.name = "Evidence", value.name = "Score")
plot_data_long$Evidence <- gsub("score_", "", plot_data_long$Evidence)
plot_data_long$Evidence <- factor(plot_data_long$Evidence,
                                    levels = c("DE", "effect", "survival",
                                               "bootstrap", "coexpr", "literature"))
plot_data_long$gene_name <- factor(plot_data_long$gene_name,
                                     levels = rev(priority$gene_name))

p_priority <- ggplot(plot_data_long, aes(x = Evidence, y = gene_name, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Score") +
  theme_classic() +
  labs(title = "Candidate lncRNA Prioritization",
       x = "Evidence Category", y = "") +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, face = "bold"))

ggsave(file.path(fig_dir, "candidate_prioritization_heatmap.pdf"), p_priority,
       width = 180, height = 200, units = "mm", dpi = 300)

cat("  Priority heatmap saved.\n")

cat("\n============================================\n")
cat(" Validation & prioritization complete.\n")
cat("============================================\n")
