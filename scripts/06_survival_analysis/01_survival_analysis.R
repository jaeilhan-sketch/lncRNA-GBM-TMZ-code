#!/usr/bin/env Rscript
# ============================================================================
# Phase 6: Survival Analysis
#
# - Kaplan-Meier survival curves for each DE-lncRNA
# - Log-rank test
# - Univariate and multivariate Cox regression
# - Optimal cutpoint determination
# - LASSO-Cox prognostic signature
# - Time-dependent ROC analysis
# ============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(DESeq2)
  library(glmnet)
  library(timeROC)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

set.seed(42)

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")

cat("============================================\n")
cat(" Phase 6: Survival Analysis\n")
cat("============================================\n")

# ── Load data ──
cat("[1/6] Loading data...\n")
load(file.path(results_dir, "de_analysis_objects.RData"))

# Prepare survival data
surv_data <- as.data.frame(col_data)
surv_data$os_time <- as.numeric(surv_data$os_days)
surv_data$os_status <- as.numeric(surv_data$os_status)
surv_data$pfs_time <- as.numeric(surv_data$pfs_days)

# VST-normalized lncRNA expression
expr_lncrna <- assay(vsd_lncrna)

# Remove samples with missing survival data
valid_idx <- !is.na(surv_data$os_time) & !is.na(surv_data$os_status) &
             surv_data$os_time > 0
surv_data <- surv_data[valid_idx, ]
expr_lncrna <- expr_lncrna[, rownames(surv_data)]

cat(sprintf("  Samples with valid OS data: %d\n", nrow(surv_data)))


# ============================================================================
# 2. Kaplan-Meier for each DE-lncRNA
# ============================================================================
cat("[2/6] Kaplan-Meier survival analysis for DE-lncRNAs...\n")

km_results <- list()
km_dir <- file.path(fig_dir, "km_curves")
dir.create(km_dir, showWarnings = FALSE)

for (i in seq_len(nrow(sig_lncrna))) {
  gene_id <- sig_lncrna$gene_id[i]
  gene_name <- sig_lncrna$gene_name[i]
  if (is.na(gene_name)) gene_name <- gene_id

  if (!(gene_id %in% rownames(expr_lncrna))) next

  expr_vals <- as.numeric(expr_lncrna[gene_id, rownames(surv_data)])

  # Optimal cutpoint
  cut_data <- data.frame(
    os_time = surv_data$os_time,
    os_status = surv_data$os_status,
    expression = expr_vals
  )

  tryCatch({
    cut_res <- surv_cutpoint(cut_data, time = "os_time", event = "os_status",
                              variables = "expression")
    cutpoint <- cut_res$cutpoint$cutpoint
  }, error = function(e) {
    cutpoint <<- median(expr_vals)
  })

  surv_data$group <- ifelse(expr_vals >= cutpoint, "High", "Low")
  surv_data$group <- factor(surv_data$group, levels = c("Low", "High"))

  # Kaplan-Meier + log-rank
  fit <- survfit(Surv(os_time, os_status) ~ group, data = surv_data)
  diff <- survdiff(Surv(os_time, os_status) ~ group, data = surv_data)
  pval <- 1 - pchisq(diff$chisq, df = 1)

  # Extract median OS from survfit summary (strata names include "group=")
  fit_tbl <- summary(fit)$table
  if (is.matrix(fit_tbl)) {
    med_high <- fit_tbl[grep("High", rownames(fit_tbl)), "median"]
    med_low  <- fit_tbl[grep("Low", rownames(fit_tbl)), "median"]
  } else {
    med_high <- NA
    med_low  <- NA
  }

  km_results[[gene_id]] <- data.frame(
    gene_id = gene_id,
    gene_name = gene_name,
    log2FC = sig_lncrna$log2FoldChange[i],
    de_padj = sig_lncrna$padj[i],
    km_pvalue = pval,
    cutpoint = cutpoint,
    n_high = sum(surv_data$group == "High"),
    n_low = sum(surv_data$group == "Low"),
    median_os_high = med_high,
    median_os_low = med_low,
    stringsAsFactors = FALSE
  )

  # Save KM plot for significant results
  if (pval < 0.05 || i <= 20) {
    p_km <- ggsurvplot(
      fit,
      data = surv_data,
      pval = TRUE,
      pval.method = TRUE,
      risk.table = TRUE,
      risk.table.col = "strata",
      palette = c("#2166AC", "#B2182B"),
      title = paste0(gene_name, " (", gene_id, ")"),
      xlab = "Time (days)",
      ylab = "Overall Survival Probability",
      legend.labs = c("Low", "High"),
      ggtheme = theme_classic()
    )

    pdf(file.path(km_dir, paste0("km_", gsub("[^A-Za-z0-9]", "_", gene_name), ".pdf")),
        width = 7, height = 6)
    print(p_km)
    dev.off()
  }
}

km_df <- rbindlist(km_results)
km_df <- km_df[order(km_df$km_pvalue), ]
fwrite(km_df, file.path(table_dir, "survival_km_results.tsv"), sep = "\t")

n_sig_km <- sum(km_df$km_pvalue < 0.05)
cat(sprintf("  DE-lncRNAs with significant OS association (p < 0.05): %d / %d\n",
            n_sig_km, nrow(km_df)))


# ============================================================================
# 3. Univariate Cox regression
# ============================================================================
cat("[3/6] Univariate Cox regression...\n")

cox_uni_results <- list()
for (i in seq_len(nrow(sig_lncrna))) {
  gene_id <- sig_lncrna$gene_id[i]
  gene_name <- sig_lncrna$gene_name[i]
  if (is.na(gene_name)) gene_name <- gene_id

  if (!(gene_id %in% rownames(expr_lncrna))) next

  surv_data$expr <- as.numeric(expr_lncrna[gene_id, rownames(surv_data)])

  tryCatch({
    cox_fit <- coxph(Surv(os_time, os_status) ~ expr, data = surv_data)
    s <- summary(cox_fit)

    cox_uni_results[[gene_id]] <- data.frame(
      gene_id = gene_id,
      gene_name = gene_name,
      HR = s$coefficients[, "exp(coef)"],
      HR_lower = s$conf.int[, "lower .95"],
      HR_upper = s$conf.int[, "upper .95"],
      cox_pvalue = s$coefficients[, "Pr(>|z|)"],
      concordance = s$concordance["C"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

cox_uni_df <- rbindlist(cox_uni_results)
cox_uni_df <- cox_uni_df[order(cox_uni_df$cox_pvalue), ]
fwrite(cox_uni_df, file.path(table_dir, "survival_cox_univariate.tsv"), sep = "\t")

cat(sprintf("  Significant in univariate Cox (p < 0.05): %d\n",
            sum(cox_uni_df$cox_pvalue < 0.05)))


# ============================================================================
# 4. Multivariate Cox regression (top candidates)
# ============================================================================
cat("[4/6] Multivariate Cox regression...\n")

# Select top lncRNAs significant in both KM and Cox
top_candidates <- intersect(
  km_df$gene_id[km_df$km_pvalue < 0.05],
  cox_uni_df$gene_id[cox_uni_df$cox_pvalue < 0.05]
)
cat(sprintf("  Candidates for multivariate Cox: %d\n", length(top_candidates)))

cox_multi_results <- list()
for (gene_id in top_candidates) {
  gene_name <- gene_annot$gene_name[gene_annot$gene_id == gene_id][1]

  surv_data$expr <- as.numeric(expr_lncrna[gene_id, rownames(surv_data)])

  # Build formula with available covariates
  covariates <- "expr"
  if ("age_years" %in% colnames(surv_data)) covariates <- c(covariates, "age_years")
  if ("gender" %in% colnames(surv_data) && !all(is.na(surv_data$gender)))
    covariates <- c(covariates, "gender")
  # Exclude MGMT_status from multivariate model to avoid losing 25 samples
  # (MGMT_status has 25 "Unknown" values which would become NA)

  formula_str <- paste0("Surv(os_time, os_status) ~ ", paste(covariates, collapse = " + "))

  tryCatch({
    cox_fit <- coxph(as.formula(formula_str), data = surv_data)
    s <- summary(cox_fit)

    # Extract lncRNA-specific results
    expr_row <- grep("^expr$", rownames(s$coefficients))
    cox_multi_results[[gene_id]] <- data.frame(
      gene_id = gene_id,
      gene_name = gene_name,
      HR = s$coefficients[expr_row, "exp(coef)"],
      HR_lower = s$conf.int[expr_row, "lower .95"],
      HR_upper = s$conf.int[expr_row, "upper .95"],
      multi_pvalue = s$coefficients[expr_row, "Pr(>|z|)"],
      covariates = paste(covariates[-1], collapse = ", "),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

if (length(cox_multi_results) > 0) {
  cox_multi_df <- rbindlist(cox_multi_results)
  cox_multi_df <- cox_multi_df[order(cox_multi_df$multi_pvalue), ]
  fwrite(cox_multi_df, file.path(table_dir, "survival_cox_multivariate.tsv"), sep = "\t")

  # Forest plot
  independent_lncrnas <- cox_multi_df$gene_id[cox_multi_df$multi_pvalue < 0.05]
  cat(sprintf("  Independent prognostic lncRNAs (multivariate p < 0.05): %d\n",
              length(independent_lncrnas)))
}


# ============================================================================
# 5. LASSO-Cox prognostic signature
# ============================================================================
cat("[5/6] Building LASSO-Cox prognostic signature...\n")

# Use all DE-lncRNAs as candidates
candidate_genes <- sig_lncrna$gene_id[sig_lncrna$gene_id %in% rownames(expr_lncrna)]
x_matrix <- t(expr_lncrna[candidate_genes, rownames(surv_data)])
y_surv <- Surv(surv_data$os_time, surv_data$os_status)

if (length(candidate_genes) >= 3 && nrow(surv_data) >= 20) {
  # LASSO-Cox with cross-validation
  cv_fit <- cv.glmnet(x_matrix, y_surv, family = "cox",
                       alpha = 1, nfolds = 10, maxit = 10000)

  # Coefficients at lambda.min
  coef_min <- coef(cv_fit, s = "lambda.min")
  selected_genes <- rownames(coef_min)[which(coef_min[, 1] != 0)]
  cat(sprintf("  LASSO selected genes (lambda.min): %d\n", length(selected_genes)))

  # Coefficients at lambda.1se (more parsimonious)
  coef_1se <- coef(cv_fit, s = "lambda.1se")
  selected_genes_1se <- rownames(coef_1se)[which(coef_1se[, 1] != 0)]
  cat(sprintf("  LASSO selected genes (lambda.1se): %d\n", length(selected_genes_1se)))

  # Calculate risk score
  risk_score <- predict(cv_fit, newx = x_matrix, s = "lambda.min", type = "link")
  surv_data$risk_score <- as.numeric(risk_score)
  surv_data$risk_group <- ifelse(risk_score > median(risk_score), "High-risk", "Low-risk")

  # KM curve for risk groups
  fit_risk <- survfit(Surv(os_time, os_status) ~ risk_group, data = surv_data)
  p_risk <- ggsurvplot(
    fit_risk,
    data = surv_data,
    pval = TRUE,
    risk.table = TRUE,
    palette = c("#B2182B", "#2166AC"),
    title = "LASSO-Cox Prognostic Signature",
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    ggtheme = theme_classic()
  )

  pdf(file.path(fig_dir, "lasso_cox_risk_km.pdf"), width = 8, height = 7)
  print(p_risk)
  dev.off()

  # LASSO CV plot
  pdf(file.path(fig_dir, "lasso_cox_cv_plot.pdf"), width = 7, height = 5)
  plot(cv_fit, main = "LASSO-Cox Cross-Validation")
  dev.off()

  # Save signature details
  sig_details <- data.frame(
    gene_id = selected_genes,
    coefficient = as.numeric(coef_min[selected_genes, 1])
  )
  sig_details <- merge(sig_details,
                        gene_annot[, c("gene_id", "gene_name")],
                        by = "gene_id", all.x = TRUE)
  fwrite(sig_details, file.path(table_dir, "lasso_signature_genes.tsv"), sep = "\t")

  # Time-dependent ROC
  tryCatch({
    roc_res <- timeROC(
      T = surv_data$os_time,
      delta = surv_data$os_status,
      marker = surv_data$risk_score,
      cause = 1,
      times = c(365, 365 * 2, 365 * 3),  # 1, 2, 3 years
      iid = TRUE
    )

    auc_df <- data.frame(
      Time = c("1 year", "2 years", "3 years"),
      AUC = roc_res$AUC,
      AUC_lower = roc_res$AUC - 1.96 * roc_res$inference$vect_sd_1,
      AUC_upper = roc_res$AUC + 1.96 * roc_res$inference$vect_sd_1
    )
    fwrite(auc_df, file.path(table_dir, "time_dependent_auc.tsv"), sep = "\t")
    cat("  Time-dependent AUC:\n")
    print(auc_df)

    pdf(file.path(fig_dir, "time_dependent_roc.pdf"), width = 7, height = 6)
    plot(roc_res, time = 365, title = "Time-dependent ROC Curves",
         col = "#B2182B", lwd = 2)
    plot(roc_res, time = 365 * 2, add = TRUE, col = "#2166AC", lwd = 2)
    plot(roc_res, time = 365 * 3, add = TRUE, col = "#4DAF4A", lwd = 2)
    legend("bottomright", legend = paste(c("1-year", "2-year", "3-year"),
           "AUC =", round(roc_res$AUC, 3)),
           col = c("#B2182B", "#2166AC", "#4DAF4A"), lwd = 2)
    dev.off()
  }, error = function(e) {
    cat("  Note: timeROC analysis skipped due to sample size.\n")
  })

} else {
  cat("  Insufficient candidates or samples for LASSO-Cox model.\n")
}


# ============================================================================
# 6. Save all results
# ============================================================================
cat("[6/6] Saving results...\n")

save(km_df, cox_uni_df,
     file = file.path(results_dir, "survival_objects.RData"))

cat("\n============================================\n")
cat(" Survival analysis complete.\n")
cat("============================================\n")

dir.create(file.path(results_dir, "reports"), showWarnings = FALSE, recursive = TRUE)
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, "reports/survival_session_info.txt"))
