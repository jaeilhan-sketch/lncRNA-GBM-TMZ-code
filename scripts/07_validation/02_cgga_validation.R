#!/usr/bin/env Rscript
# ============================================================================
# Phase 7b: CGGA External Validation
#
# Validates key DE-lncRNAs (H19, LINC00707, CYP1B1-AS1) in CGGA dataset.
# - Differential expression validation
# - Survival analysis validation
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(gridExtra)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
valid_dir   <- file.path(results_dir, "validation")
cgga_dir    <- file.path(project_dir, "data/external/cgga")
dir.create(valid_dir, showWarnings = FALSE, recursive = TRUE)

# Publication theme
theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

cat("============================================\n")
cat(" Phase 7b: CGGA External Validation\n")
cat("============================================\n")

# ── Load TCGA results ──
cat("[1/5] Loading TCGA results...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))
load(file.path(results_dir, "survival_objects.RData"))

# ── Load CGGA data ──
cat("[2/5] Loading CGGA data...\n")

# mRNAseq_693 (larger cohort)
cgga_expr <- fread(file.path(cgga_dir, "CGGA.mRNAseq_693.RSEM-genes.20200506.txt"))
cgga_clin <- fread(file.path(cgga_dir, "CGGA.mRNAseq_693_clinical.20200506.txt"))

# Standardize column names
setnames(cgga_clin, old = c("Censor (alive=0; dead=1)",
                              "Chemo_status (TMZ treated=1;un-treated=0)",
                              "Radio_status (treated=1;un-treated=0)"),
         new = c("Censor", "Chemo_TMZ", "Radio"), skip_absent = TRUE)

cat(sprintf("  CGGA mRNAseq_693: %d samples, %d genes\n",
            ncol(cgga_expr) - 1, nrow(cgga_expr)))

# Filter GBM patients
cgga_gbm <- cgga_clin[Histology == "GBM"]
cat(sprintf("  GBM patients: %d\n", nrow(cgga_gbm)))

# Exclude IDH-mutant
cgga_gbm <- cgga_gbm[IDH_mutation_status == "Wildtype" | is.na(IDH_mutation_status)]
cat(sprintf("  After IDH-wt filter: %d\n", nrow(cgga_gbm)))

# Identify TMZ-treated vs untreated
cgga_gbm$tmz_status <- ifelse(cgga_gbm$Chemo_TMZ == 1, "TMZ_treated", "No_TMZ")
cat(sprintf("  TMZ treated: %d, No TMZ: %d, NA: %d\n",
            sum(cgga_gbm$tmz_status == "TMZ_treated", na.rm = TRUE),
            sum(cgga_gbm$tmz_status == "No_TMZ", na.rm = TRUE),
            sum(is.na(cgga_gbm$Chemo_TMZ))))

# ── Match expression data ──
gene_col <- colnames(cgga_expr)[1]
cgga_samples <- intersect(cgga_gbm$CGGA_ID, colnames(cgga_expr))
cat(sprintf("  Matched GBM samples with expression: %d\n", length(cgga_samples)))
cgga_gbm <- cgga_gbm[CGGA_ID %in% cgga_samples]


# ============================================================================
# 3. Validate expression direction of candidate lncRNAs
# ============================================================================
cat("\n[3/5] Validating lncRNA expression in CGGA...\n")

# lncRNAs available in CGGA
target_lncrna <- c("H19", "LINC00707", "CYP1B1-AS1")
available_lncrna <- target_lncrna[target_lncrna %in% cgga_expr[[gene_col]]]
cat(sprintf("  Available for validation: %s\n", paste(available_lncrna, collapse = ", ")))

validation_results <- list()

for (gene_name in available_lncrna) {
  # Get expression for this gene
  gene_row <- cgga_expr[get(gene_col) == gene_name]
  if (nrow(gene_row) == 0) next

  # Extract expression values for GBM samples
  expr_vals <- as.numeric(gene_row[1, cgga_samples, with = FALSE])
  names(expr_vals) <- cgga_samples

  # ── TMZ-treated subgroup: compare high vs low expression with survival ──
  tmz_treated <- cgga_gbm[Chemo_TMZ == 1 & !is.na(OS) & OS > 0 & !is.na(Censor)]

  if (nrow(tmz_treated) >= 20) {
    tmz_expr <- expr_vals[tmz_treated$CGGA_ID]
    tmz_treated$expression <- tmz_expr
    tmz_treated <- tmz_treated[!is.na(expression) & expression > 0]

    if (nrow(tmz_treated) >= 20) {
      # Cox regression
      cox_fit <- coxph(Surv(OS, Censor) ~ expression, data = tmz_treated)
      s <- summary(cox_fit)

      # KM with median split
      tmz_treated$group <- ifelse(tmz_treated$expression > median(tmz_treated$expression),
                                   "High", "Low")
      fit_km <- survfit(Surv(OS, Censor) ~ group, data = tmz_treated)
      diff_km <- survdiff(Surv(OS, Censor) ~ group, data = tmz_treated)
      km_pval <- 1 - pchisq(diff_km$chisq, df = 1)

      validation_results[[gene_name]] <- data.frame(
        gene = gene_name,
        dataset = "CGGA_693_GBM_TMZ",
        n_samples = nrow(tmz_treated),
        HR = s$coefficients[, "exp(coef)"],
        HR_lower = s$conf.int[, "lower .95"],
        HR_upper = s$conf.int[, "upper .95"],
        cox_pvalue = s$coefficients[, "Pr(>|z|)"],
        km_pvalue = km_pval,
        validated = s$coefficients[, "Pr(>|z|)"] < 0.05 | km_pval < 0.05,
        stringsAsFactors = FALSE
      )

      cat(sprintf("  %s (TMZ-treated GBM, n=%d): HR=%.2f, Cox_p=%.4f, KM_p=%.4f %s\n",
                  gene_name, nrow(tmz_treated),
                  s$coefficients[, "exp(coef)"],
                  s$coefficients[, "Pr(>|z|)"],
                  km_pval,
                  ifelse(s$coefficients[, "Pr(>|z|)"] < 0.05 | km_pval < 0.05,
                         "<-- VALIDATED", "")))

      # Save KM plot
      p_km <- ggsurvplot(
        fit_km, data = tmz_treated,
        pval = TRUE, pval.method = TRUE,
        risk.table = TRUE,
        palette = c("#B2182B", "#2166AC"),
        title = paste0("CGGA Validation: ", gene_name, " (TMZ-treated GBM, n=",
                       nrow(tmz_treated), ")"),
        xlab = "Time (days)", ylab = "Overall Survival",
        legend.labs = c("High", "Low"),
        ggtheme = theme_pub
      )

      pdf(file.path(fig_dir, paste0("cgga_km_", gene_name, ".pdf")), width = 7, height = 6)
      print(p_km)
      dev.off()
    }
  }

  # ── All GBM: expression comparison by OS outcome ──
  all_gbm <- cgga_gbm[!is.na(OS) & OS > 0 & !is.na(Censor)]
  all_gbm$expression <- expr_vals[all_gbm$CGGA_ID]
  all_gbm <- all_gbm[!is.na(expression) & expression > 0]

  if (nrow(all_gbm) >= 20) {
    cox_all <- coxph(Surv(OS, Censor) ~ expression, data = all_gbm)
    s_all <- summary(cox_all)

    all_gbm$group <- ifelse(all_gbm$expression > median(all_gbm$expression), "High", "Low")
    fit_all <- survfit(Surv(OS, Censor) ~ group, data = all_gbm)
    diff_all <- survdiff(Surv(OS, Censor) ~ group, data = all_gbm)
    km_pval_all <- 1 - pchisq(diff_all$chisq, df = 1)

    cat(sprintf("  %s (All GBM, n=%d): HR=%.2f, Cox_p=%.4f, KM_p=%.4f\n",
                gene_name, nrow(all_gbm),
                s_all$coefficients[, "exp(coef)"],
                s_all$coefficients[, "Pr(>|z|)"],
                km_pval_all))

    # Save KM plot for all GBM
    p_km_all <- ggsurvplot(
      fit_all, data = all_gbm,
      pval = TRUE, pval.method = TRUE,
      risk.table = TRUE,
      palette = c("#B2182B", "#2166AC"),
      title = paste0("CGGA All GBM: ", gene_name, " (n=", nrow(all_gbm), ")"),
      xlab = "Time (days)", ylab = "Overall Survival",
      legend.labs = c("High", "Low"),
      ggtheme = theme_pub
    )

    pdf(file.path(fig_dir, paste0("cgga_km_all_gbm_", gene_name, ".pdf")),
        width = 7, height = 6)
    print(p_km_all)
    dev.off()
  }
}


# ============================================================================
# 4. Also validate with mRNAseq_325 dataset
# ============================================================================
cat("\n[4/5] Validating in CGGA mRNAseq_325...\n")

cgga_expr2 <- fread(file.path(cgga_dir, "CGGA.mRNAseq_325.RSEM-genes.20200506.txt"))
cgga_clin2 <- fread(file.path(cgga_dir, "CGGA.mRNAseq_325_clinical.20200506.txt"))
setnames(cgga_clin2, old = c("Censor (alive=0; dead=1)",
                               "Chemo_status (TMZ treated=1;un-treated=0)",
                               "Radio_status (treated=1;un-treated=0)"),
         new = c("Censor", "Chemo_TMZ", "Radio"), skip_absent = TRUE)

cgga_gbm2 <- cgga_clin2[Histology == "GBM" &
                            (IDH_mutation_status == "Wildtype" | is.na(IDH_mutation_status))]
cgga_samples2 <- intersect(cgga_gbm2$CGGA_ID, colnames(cgga_expr2))
cgga_gbm2 <- cgga_gbm2[CGGA_ID %in% cgga_samples2]

cat(sprintf("  mRNAseq_325 IDH-wt GBM: %d\n", nrow(cgga_gbm2)))

for (gene_name in available_lncrna) {
  gene_row2 <- cgga_expr2[get(colnames(cgga_expr2)[1]) == gene_name]
  if (nrow(gene_row2) == 0) next

  expr_vals2 <- as.numeric(gene_row2[1, cgga_samples2, with = FALSE])
  val_data2 <- cgga_gbm2[!is.na(OS) & OS > 0 & !is.na(Censor)]
  val_data2$expression <- expr_vals2[match(val_data2$CGGA_ID, cgga_samples2)]
  val_data2 <- val_data2[!is.na(expression) & expression > 0]

  if (nrow(val_data2) >= 10) {
    cox2 <- coxph(Surv(OS, Censor) ~ expression, data = val_data2)
    s2 <- summary(cox2)

    val_data2$group <- ifelse(val_data2$expression > median(val_data2$expression), "High", "Low")
    diff2 <- survdiff(Surv(OS, Censor) ~ group, data = val_data2)
    km_pval2 <- 1 - pchisq(diff2$chisq, df = 1)

    validation_results[[paste0(gene_name, "_325")]] <- data.frame(
      gene = gene_name,
      dataset = "CGGA_325_GBM",
      n_samples = nrow(val_data2),
      HR = s2$coefficients[, "exp(coef)"],
      HR_lower = s2$conf.int[, "lower .95"],
      HR_upper = s2$conf.int[, "upper .95"],
      cox_pvalue = s2$coefficients[, "Pr(>|z|)"],
      km_pvalue = km_pval2,
      validated = s2$coefficients[, "Pr(>|z|)"] < 0.05 | km_pval2 < 0.05,
      stringsAsFactors = FALSE
    )

    cat(sprintf("  %s (325 GBM, n=%d): HR=%.2f, Cox_p=%.4f, KM_p=%.4f\n",
                gene_name, nrow(val_data2),
                s2$coefficients[, "exp(coef)"],
                s2$coefficients[, "Pr(>|z|)"],
                km_pval2))
  }
}


# ============================================================================
# 5. Compile and save validation results
# ============================================================================
cat("\n[5/5] Compiling validation results...\n")

if (length(validation_results) > 0) {
  val_df <- rbindlist(validation_results)
  fwrite(val_df, file.path(valid_dir, "cgga_validation_results.tsv"), sep = "\t")

  cat("\n  CGGA Validation Summary:\n")
  cat("  ══════════════════════════════════════════════\n")
  for (i in seq_len(nrow(val_df))) {
    cat(sprintf("  %s (%s, n=%d): HR=%.2f [%.2f-%.2f], Cox_p=%.4f, KM_p=%.4f %s\n",
                val_df$gene[i], val_df$dataset[i], val_df$n_samples[i],
                val_df$HR[i], val_df$HR_lower[i], val_df$HR_upper[i],
                val_df$cox_pvalue[i], val_df$km_pvalue[i],
                ifelse(val_df$validated[i], "VALIDATED", "")))
  }
  cat("  ══════════════════════════════════════════════\n")
} else {
  cat("  No validation results obtained.\n")
}

cat("\n============================================\n")
cat(" CGGA external validation complete.\n")
cat("============================================\n")
