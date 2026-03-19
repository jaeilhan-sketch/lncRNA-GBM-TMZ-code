#!/usr/bin/env Rscript
# ============================================================================
# Generate Manuscript Tables
#
# Table 1: Baseline Patient Characteristics (stratified by TMZ response)
# Table 2: Multivariate Cox + Bootstrap Validation (H19 & LINC01936)
# ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(writexl)
})

project_dir <- "/media/jaeil/ff387598-88d2-44f2-b68e-cf799d26fb68/lncRNA-GBM-TMZ"
results_dir <- file.path(project_dir, "results")
table_dir   <- file.path(project_dir, "manuscript", "tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("============================================\n")
cat(" Generating Manuscript Tables\n")
cat("============================================\n")

# ── Load data ──
load(file.path(results_dir, "de_analysis_objects.RData"))

# ============================================================================
# Table 1: Baseline Patient Characteristics
# ============================================================================
cat("[1/2] Generating Table 1: Baseline characteristics...\n")

cd <- as.data.frame(col_data)

# Helper functions
fmt_mean_sd <- function(x, digits = 1) {
  x <- as.numeric(x[!is.na(x)])
  sprintf("%.1f ± %.1f", mean(x), sd(x))
}

fmt_median_range <- function(x, digits = 1) {
  x <- as.numeric(x[!is.na(x)])
  sprintf("%.1f (%.1f–%.1f)", median(x), min(x), max(x))
}

fmt_n_pct <- function(x, total) {
  sprintf("%d (%.1f%%)", x, 100 * x / total)
}

# Split by response
resp   <- cd[cd$tmz_response == "Responder", ]
nresp  <- cd[cd$tmz_response == "NonResponder", ]
n_resp  <- nrow(resp)
n_nresp <- nrow(nresp)
n_total <- nrow(cd)

# ── Continuous: age ──
age_all   <- fmt_mean_sd(cd$age_years)
age_resp  <- fmt_mean_sd(resp$age_years)
age_nresp <- fmt_mean_sd(nresp$age_years)
age_p     <- format(t.test(resp$age_years, nresp$age_years)$p.value, digits = 3)

# ── Continuous: OS (days) ──
os_all    <- fmt_median_range(cd$os_days)
os_resp   <- fmt_median_range(resp$os_days)
os_nresp  <- fmt_median_range(nresp$os_days)
os_p_val  <- tryCatch({
  library(survival)
  sf <- survdiff(Surv(as.numeric(os_days), as.numeric(os_status)) ~ tmz_response,
                 data = cd)
  p <- pchisq(sf$chisq, df = 1, lower.tail = FALSE)
  if (p < 0.001) "<0.001" else format(round(p, 3), nsmall = 3)
}, error = function(e) "N/A")

# ── Continuous: PFS (months) ──
pfs_all   <- fmt_median_range(cd$pfs_months)
pfs_resp  <- fmt_median_range(resp$pfs_months)
pfs_nresp <- fmt_median_range(nresp$pfs_months)

# ── Categorical: Gender ──
male_all   <- sum(cd$gender == "male",   na.rm = TRUE)
male_resp  <- sum(resp$gender == "male",  na.rm = TRUE)
male_nresp <- sum(nresp$gender == "male", na.rm = TRUE)
gender_tab <- table(cd$gender, cd$tmz_response)
gender_p   <- format(round(chisq.test(gender_tab)$p.value, 3), nsmall = 3)

# ── Categorical: MGMT status ──
mgmt_lev <- c("METHYLATED", "UNMETHYLATED", "Unknown")
mgmt_counts <- function(df) {
  sapply(mgmt_lev, function(l) sum(df$MGMT_STATUS == l, na.rm = TRUE))
}
mgmt_all   <- mgmt_counts(cd)
mgmt_resp  <- mgmt_counts(resp)
mgmt_nresp <- mgmt_counts(nresp)
mgmt_tab   <- table(cd$MGMT_STATUS, cd$tmz_response)
mgmt_p     <- tryCatch(
  format(round(fisher.test(mgmt_tab)$p.value, 3), nsmall = 3),
  error = function(e) format(round(chisq.test(mgmt_tab)$p.value, 3), nsmall = 3)
)

# ── Categorical: Expression subtype ──
expr_lev <- c("Classical", "Mesenchymal", "Neural", "Proneural", "Unknown")
expr_counts <- function(df) {
  sapply(expr_lev, function(l) sum(df$EXPRESSION_SUBTYPE == l, na.rm = TRUE))
}
expr_all   <- expr_counts(cd)
expr_resp  <- expr_counts(resp)
expr_nresp <- expr_counts(nresp)
expr_tab   <- table(cd$EXPRESSION_SUBTYPE, cd$tmz_response)
expr_p     <- tryCatch(
  format(round(fisher.test(expr_tab, simulate.p.value = TRUE)$p.value, 3), nsmall = 3),
  error = function(e) format(round(chisq.test(expr_tab)$p.value, 3), nsmall = 3)
)

# ── OS events ──
event_all   <- sum(as.numeric(cd$os_status) == 1, na.rm = TRUE)
event_resp  <- sum(as.numeric(resp$os_status) == 1, na.rm = TRUE)
event_nresp <- sum(as.numeric(nresp$os_status) == 1, na.rm = TRUE)

# ── Build Table 1 ──
table1_rows <- data.frame(
  Characteristic = c(
    sprintf("N"),
    "Age at diagnosis (years), mean ± SD",
    "Gender, n (%)",
    "  Male",
    "  Female",
    "MGMT promoter status, n (%)",
    "  Methylated",
    "  Unmethylated",
    "  Unknown",
    "Expression subtype, n (%)",
    "  Classical",
    "  Mesenchymal",
    "  Neural",
    "  Proneural",
    "  Unknown",
    "PFS, median (range), months",
    "OS, median (range), days",
    "Deaths (OS events), n (%)"
  ),
  Total = c(
    as.character(n_total),
    age_all,
    "",
    fmt_n_pct(male_all, n_total),
    fmt_n_pct(n_total - male_all, n_total),
    "",
    fmt_n_pct(mgmt_all["METHYLATED"],   n_total),
    fmt_n_pct(mgmt_all["UNMETHYLATED"], n_total),
    fmt_n_pct(mgmt_all["Unknown"],      n_total),
    "",
    fmt_n_pct(expr_all["Classical"],   n_total),
    fmt_n_pct(expr_all["Mesenchymal"], n_total),
    fmt_n_pct(expr_all["Neural"],      n_total),
    fmt_n_pct(expr_all["Proneural"],   n_total),
    fmt_n_pct(expr_all["Unknown"],     n_total),
    pfs_all,
    os_all,
    fmt_n_pct(event_all, n_total)
  ),
  Responder = c(
    as.character(n_resp),
    age_resp,
    "",
    fmt_n_pct(male_resp, n_resp),
    fmt_n_pct(n_resp - male_resp, n_resp),
    "",
    fmt_n_pct(mgmt_resp["METHYLATED"],   n_resp),
    fmt_n_pct(mgmt_resp["UNMETHYLATED"], n_resp),
    fmt_n_pct(mgmt_resp["Unknown"],      n_resp),
    "",
    fmt_n_pct(expr_resp["Classical"],   n_resp),
    fmt_n_pct(expr_resp["Mesenchymal"], n_resp),
    fmt_n_pct(expr_resp["Neural"],      n_resp),
    fmt_n_pct(expr_resp["Proneural"],   n_resp),
    fmt_n_pct(expr_resp["Unknown"],     n_resp),
    pfs_resp,
    os_resp,
    fmt_n_pct(event_resp, n_resp)
  ),
  NonResponder = c(
    as.character(n_nresp),
    age_nresp,
    "",
    fmt_n_pct(male_nresp, n_nresp),
    fmt_n_pct(n_nresp - male_nresp, n_nresp),
    "",
    fmt_n_pct(mgmt_nresp["METHYLATED"],   n_nresp),
    fmt_n_pct(mgmt_nresp["UNMETHYLATED"], n_nresp),
    fmt_n_pct(mgmt_nresp["Unknown"],      n_nresp),
    "",
    fmt_n_pct(expr_nresp["Classical"],   n_nresp),
    fmt_n_pct(expr_nresp["Mesenchymal"], n_nresp),
    fmt_n_pct(expr_nresp["Neural"],      n_nresp),
    fmt_n_pct(expr_nresp["Proneural"],   n_nresp),
    fmt_n_pct(expr_nresp["Unknown"],     n_nresp),
    pfs_nresp,
    os_nresp,
    fmt_n_pct(event_nresp, n_nresp)
  ),
  p_value = c(
    "",          # N
    age_p,       # Age
    gender_p,    # Gender header
    "", "",      # Male, Female
    mgmt_p,      # MGMT header
    "", "", "",  # MGMT levels
    expr_p,      # Expression subtype header
    "", "", "", "", "",  # Subtype levels
    "",          # PFS (by definition differs)
    os_p_val,   # OS log-rank
    ""           # Deaths
  ),
  stringsAsFactors = FALSE
)

out1 <- file.path(table_dir, "Table1_baseline_characteristics.xlsx")
write_xlsx(table1_rows, out1)
cat(sprintf("  Saved: %s\n", out1))


# ============================================================================
# Table 2: Multivariate Cox Regression + Bootstrap Validation
# ============================================================================
cat("[2/2] Generating Table 2: Cox + Bootstrap results...\n")

cox  <- fread(file.path(results_dir, "tables", "survival_cox_multivariate.tsv"))
boot <- fread(file.path(results_dir, "validation", "bootstrap_validation.tsv"))

# Keep only H19 and LINC01936
genes_of_interest <- c("H19", "LINC01936")
cox  <- cox[gene_name %in% genes_of_interest]
boot <- boot[gene_name %in% genes_of_interest]

# Merge
tbl2 <- merge(cox, boot, by = c("gene_id", "gene_name"), all.x = TRUE)

# Format columns
tbl2[, HR_CI_cox     := sprintf("%.2f (%.2f–%.2f)", HR, HR_lower, HR_upper)]
tbl2[, p_value_cox   := ifelse(multi_pvalue < 0.001, "<0.001",
                                sprintf("%.3f", multi_pvalue))]
tbl2[, Boot_sig_rate := sprintf("%.1f%%", boot_sig_rate * 100)]
tbl2[, HR_CI_boot    := sprintf("%.2f (%.2f–%.2f)", median_HR, HR_95_lower, HR_95_upper)]
tbl2[, Covariates    := covariates]

table2 <- tbl2[, .(
  Gene           = gene_name,
  `HR (95% CI) — Cox Multivariate` = HR_CI_cox,
  `p-value`      = p_value_cox,
  Covariates,
  `Bootstrap Significance Rate (n=1000)` = Boot_sig_rate,
  `Median HR (95% CI) — Bootstrap`      = HR_CI_boot
)]

# Reorder: H19 first
table2 <- table2[order(match(Gene, genes_of_interest))]

out2 <- file.path(table_dir, "Table2_cox_bootstrap.xlsx")
write_xlsx(as.data.frame(table2), out2)
cat(sprintf("  Saved: %s\n", out2))

# ── Print preview ──
cat("\n── Table 1 preview ──\n")
print(table1_rows, row.names = FALSE)

cat("\n── Table 2 ──\n")
print(as.data.frame(table2), row.names = FALSE)

cat("\n============================================\n")
cat(" Tables saved to manuscript/tables/\n")
cat("============================================\n")
