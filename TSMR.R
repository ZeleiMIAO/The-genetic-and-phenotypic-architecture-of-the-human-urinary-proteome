#!/usr/bin/env Rscript

# ==============================================================================
# Two-Sample Mendelian Randomization Analysis
#
# Description:
#   Performs two-sample MR using clumped exposure instruments and fastGWA
#   outcome summary statistics.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(future)
  library(future.apply)
})

options(future.globals.maxSize = 16 * 1024^3)

# ==============================================================================
# Settings
# ==============================================================================

exp_dir <- "<PROJECT_DIR>/MR/trait/3.MR_input/"
out_dir <- "<PROJECT_DIR>/MR/outcome/"
res_dir <- "<PROJECT_DIR>/MR/TSresults/"

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
plan(multisession, workers = 64)

exp_files <- list.files(exp_dir, pattern = "\\.clump\\.summary\\.txt$", full.names = TRUE)
out_files <- list.files(out_dir, pattern = "_gwas\\.fastGWA$", full.names = TRUE)


# ==============================================================================
# 1. Format exposure data
# ==============================================================================

format_exposure_std <- function(dt, name) {
  dt <- dt[, .(SNP, A1, A2, freq, b, se, p)]
  dt[, p := as.numeric(p)]
  dt <- dt[!is.na(p)]
  dt[, phenotype := name]

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

read_outcome_fastgwa <- function(path, name) {
  dt <- fread(path, showProgress = FALSE)
  dt[, SNP := sub("^([^:]+:[^:]+):.*", "\\1", SNP)]

  dt <- dt[, .(SNP, beta = BETA, se = SE, A1, A2, p = P, eaf = AF1)]
  dt[, phenotype := name]

  format_data(
    as.data.frame(dt), type = "outcome",
    snp_col = "SNP", beta_col = "beta", se_col = "se",
    effect_allele_col = "A1", other_allele_col = "A2",
    eaf_col = "eaf", pval_col = "p", phenotype_col = "phenotype"
  )
}


# ==============================================================================
# 3. Run MR
# ==============================================================================

run_one_exposure <- function(ef) {
  exp_name <- sub("\\.clump\\.summary\\.txt$", "", basename(ef))
  message("Exposure: ", exp_name)

  exp_fmt <- tryCatch(
    format_exposure_std(fread(ef), exp_name),
    error = function(e) {
      message("Exposure format failed: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(exp_fmt) || nrow(exp_fmt) == 0) return(NULL)

  res_list <- future_lapply(out_files, function(of) {
    out_name <- sub("_gwas\\.fastGWA$", "", basename(of))

    tryCatch({
      out_fmt <- read_outcome_fastgwa(of, out_name)

      snps <- intersect(exp_fmt$SNP, out_fmt$SNP)
      if (length(snps) == 0) return(NULL)

      exp_sub <- exp_fmt[exp_fmt$SNP %in% snps, ]
      out_sub <- out_fmt[out_fmt$SNP %in% snps, ]

      exp_sub$id.exposure <- exp_name
      out_sub$id.outcome <- out_name

      harm <- harmonise_data(exp_sub, out_sub, action = 1)
      if (nrow(harm) == 0) return(NULL)

      mr_res <- mr(harm)
      if (is.null(mr_res)) return(NULL)

      mr_res$exposure <- exp_name
      mr_res$outcome <- out_name
      mr_res

    }, error = function(e) {
      message("Failed: ", exp_name, " -> ", out_name, ": ", conditionMessage(e))
      NULL
    })
  })

  res_all <- rbindlist(res_list, fill = TRUE)
  if (nrow(res_all) == 0) return(NULL)

  fwrite(res_all, file.path(res_dir, paste0(exp_name, "_MR_all_results.txt")), sep = "\t")
  message("Saved: ", exp_name, " | ", nrow(res_all), " rows")

  invisible(res_all)
}


# ==============================================================================
# 4. Main analysis
# ==============================================================================

for (i in seq_along(exp_files)) {
  exp_name <- sub("\\.clump\\.summary\\.txt$", "", basename(exp_files[i]))
  message("[", i, "/", length(exp_files), "] ", exp_name)

  t0 <- Sys.time()

  tryCatch(
    run_one_exposure(exp_files[i]),
    error = function(e) message("Exposure failed: ", exp_name, ": ", conditionMessage(e))
  )

  message("Finished: ", exp_name, " | ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")

  gc()
}