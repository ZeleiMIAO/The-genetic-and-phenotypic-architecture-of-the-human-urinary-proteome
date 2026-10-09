#!/usr/bin/env Rscript

# ==============================================================================
# Radial Mendelian Randomization Analysis
#
# Description:
#   Performs harmonisation, RadialMR outlier detection and removal, and
#   two-sample MR across all exposure-outcome combinations.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(RadialMR)
  library(future)
  library(future.apply)
})

# ==============================================================================
# Settings
# ==============================================================================

options(future.globals.maxSize = 4 * 1024^3)
plan(multisession, workers = 64)

exp_dir <- "<PROJECT_DIR>/MR/trait/TSMR_input_10MB/"
out_dir <- "<PROJECT_DIR>/MR/outcome/"
res_dir <- "<PROJECT_DIR>/MR/TSresults_10MB_radial_each/"

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

exp_files <- list.files(exp_dir, pattern = "\\.clump\\.summary\\.txt$", full.names = TRUE)
out_files <- list.files(out_dir, pattern = "_gwas\\.fastGWA$", full.names = TRUE)


# ==============================================================================
# 1. Format exposure data
# ==============================================================================

format_exposure_std <- function(dt, exposure_name) {
  dt <- dt[, .(SNP, A1, A2, freq, b, se, p)]
  dt[, p := as.numeric(p)]
  dt <- dt[!is.na(p)]

  dt[, phenotype := exposure_name]
  dt[, exposure := exposure_name]
  dt[, id.exposure := exposure_name]

  format_data(
    as.data.frame(dt), type = "exposure",
    snp_col = "SNP", beta_col = "b", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "freq", pval_col = "p", phenotype_col = "phenotype"
  )
}


# ==============================================================================
# 2. Format outcome data
# ==============================================================================

read_outcome_fastgwa <- function(path, outcome_name) {
  dt <- fread(path, showProgress = FALSE)

  # Convert SNP ID from CHR:POS:A1:A2 to CHR:POS
  dt[, SNP_key := sub("^([^:]+:[^:]+):.*", "\\1", SNP)]

  dt2 <- dt[, .(SNP = SNP_key, beta = BETA, se = SE, A1, A2, p = P, eaf = AF1)]
  dt2[, p := as.numeric(p)]
  dt2 <- dt2[!is.na(p)]

  dt2[, phenotype := outcome_name]
  dt2[, outcome := outcome_name]
  dt2[, id.outcome := outcome_name]

  format_data(
    as.data.frame(dt2), type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "eaf", pval_col = "p", phenotype_col = "phenotype"
  )
}


# ==============================================================================
# 3. Run one exposure-outcome pair
# ==============================================================================

run_one_pair <- function(ef, of) {

  exp_name <- sub("\\.clump\\.summary\\.txt$", "", basename(ef))
  out_name <- sub("_gwas\\.fastGWA$", "", basename(of))
  out_file <- file.path(res_dir, paste0(exp_name, "__", out_name, "_MR.txt"))

  # Skip completed tasks
  if (file.exists(out_file)) return(NULL)

  # Read exposure
  exp_raw <- tryCatch(fread(ef, showProgress = FALSE), error = function(e) NULL)
  if (is.null(exp_raw)) return(NULL)

  exp_fmt <- tryCatch(format_exposure_std(exp_raw, exp_name), error = function(e) NULL)
  if (is.null(exp_fmt) || nrow(exp_fmt) < 1) return(NULL)

  # Read outcome
  out_fmt <- tryCatch(read_outcome_fastgwa(of, out_name), error = function(e) NULL)
  if (is.null(out_fmt) || nrow(out_fmt) < 1) return(NULL)

  # Keep overlapping SNPs
  snps <- intersect(exp_fmt$SNP, out_fmt$SNP)
  if (length(snps) < 1) return(NULL)

  # Harmonise alleles
  harm <- harmonise_data(
    exp_fmt[exp_fmt$SNP %in% snps, ],
    out_fmt[out_fmt$SNP %in% snps, ],
    action = 1
  )
  if (nrow(harm) < 1) return(NULL)

  # RadialMR outlier detection
  radial_res <- tryCatch(
    ivw_radial(harm, alpha = 0.05, weights = 1),
    error = function(e) NULL
  )

  outliers <- character(0)
  if (!is.null(radial_res$outliers) &&
      !identical(radial_res$outliers, "No significant outliers")) {
    outliers <- as.character(radial_res$outliers$SNP)
  }

  # Remove outliers
  harm2 <- harm
  if (length(outliers) > 0) harm2 <- harm[!harm$SNP %in% outliers, ]
  if (nrow(harm2) < 1) return(NULL)

  # Main MR analysis
  mr_dt <- as.data.table(mr(harm2))

  # Egger intercept test
  egger_int <- tryCatch(mr_pleiotropy_test(harm2), error = function(e) NULL)

  if (!is.null(egger_int) && nrow(egger_int) == 1) {
    egger_int$method <- "egger_intercept"
    egger_int$nsnp <- NA
    egger_int <- egger_int[, c(
      "id.exposure", "id.outcome", "outcome", "exposure",
      "method", "nsnp", "egger_intercept", "se", "pval"
    )]
    colnames(egger_int)[7] <- "b"
    mr_dt <- rbind(mr_dt, as.data.table(egger_int), fill = TRUE)
  }

  # Record outliers
  mr_dt[, radial_outliers_n := length(outliers)]
  mr_dt[, radial_outliers := ifelse(
    length(outliers) > 0, paste(outliers, collapse = ","), NA
  )]

  # Save results
  fwrite(mr_dt, out_file, sep = "\t")
  return(out_file)
}


# ==============================================================================
# 4. Main analysis
# ==============================================================================

tasks <- CJ(ef = exp_files, of = out_files)
message("Total tasks: ", nrow(tasks))

future_lapply(seq_len(nrow(tasks)), function(i) {
  run_one_pair(tasks$ef[i], tasks$of[i])
})

message("ALL TASKS DONE")