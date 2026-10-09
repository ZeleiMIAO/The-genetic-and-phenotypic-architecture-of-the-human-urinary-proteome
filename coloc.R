#!/usr/bin/env Rscript
# ==============================================================================
# Protein–Trait Colocalization Analysis
#
# Description:
#   Bayesian colocalization of regional protein GWAS and quantitative or
#   case-control trait GWAS, with allele alignment and variant matching.
#
# Notes:
#   - Protein A1 is assumed to be the minor allele (AF1 = MAF).
#   - Case-control traits require total sample size and case count.
# ==============================================================================

rm(list = ls())
suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(dplyr)
  library(purrr)
  library(stringr)
})

# ==============================================================================
# Settings
# ==============================================================================
PROTEIN_INPUT_DIR <- "<PROJECT_DIR>/coloc/protein_input" #Pre-extracted regional fastGWA files.
RESULT_DIR <- "<PROJECT_DIR>/coloc/results"
N_PROTEIN <- <PROTEIN_SAMPLE_SIZE>
dir.create(RESULT_DIR, recursive = TRUE, showWarnings = FALSE)

# type: quant = quantitative; cc = case-control. n_case = NA for quantitative.
trait_config <- data.table(
  trait = c("eGFR", "BUN", "CKD"),
  input_file = c("<PROJECT_DIR>/traits/eGFR_full_standard.txt",
                 "<PROJECT_DIR>/traits/BUN_full_standard.txt",
                 "<PROJECT_DIR>/traits/CKD_full_standard.txt"),
  type = c("quant", "quant", "cc"),
  n_total = c(<EGFR_SAMPLE_SIZE>, <BUN_SAMPLE_SIZE>, <CKD_TOTAL_SAMPLE_SIZE>),
  n_case = c(NA, NA, <CKD_CASE_COUNT>)
)

# Optional CSV with exposure (protein ID) and outcome (trait); NULL = all pairs.
SELECTION_FILE <- NULL
# SELECTION_FILE <- "<PROJECT_DIR>/coloc/selected_protein_trait_pairs.csv"

# ==============================================================================
# 1. Read and prepare data
# ==============================================================================
complement_allele <- function(allele) chartr("ATCG", "TAGC", toupper(allele))

read_protein_data <- function(input_dir) {
  protein_files <- list.files(input_dir, pattern = "_coloc\\.txt$", full.names = TRUE)
  if (length(protein_files) == 0) stop("No regional protein GWAS files found in: ", input_dir)

  protein_data <- map_dfr(protein_files, function(path) {
    dt <- fread(path, showProgress = FALSE)
    protein_id <- sub("_chr.*_coloc\\.txt$", "", basename(path))
    dt[, protein_id := protein_id]
    dt
  })

  required_columns <- c("CHR", "POS", "SNP", "A1", "A2", "AF1", "BETA", "SE", "P")
  missing_columns <- setdiff(required_columns, names(protein_data))
  if (length(missing_columns) > 0) stop("Missing protein GWAS columns: ", paste(missing_columns, collapse = ", "))

  protein_data %>%
    transmute(
      protein_id = as.character(protein_id), chromosome = as.character(CHR),
      position = as.numeric(POS),
      snp = str_remove(as.character(SNP), ":[^:]+:[^:]+$"),
      ref_protein = toupper(as.character(A2)), alt_protein = toupper(as.character(A1)),
      beta_protein = as.numeric(BETA), varbeta_protein = as.numeric(SE)^2,
      pvalue_protein = as.numeric(P), maf = as.numeric(AF1)
    ) %>%
    filter(!is.na(protein_id), !is.na(snp), !is.na(beta_protein),
           !is.na(varbeta_protein), !is.na(maf), varbeta_protein > 0, maf > 0, maf < 1) %>%
    distinct()
}

read_trait_data <- function(path) {
  trait_data <- fread(path, showProgress = FALSE)
  required_columns <- c("SNP", "A1", "A2", "b", "se", "p")
  missing_columns <- setdiff(required_columns, names(trait_data))
  if (length(missing_columns) > 0) {
    stop("Missing trait GWAS columns in ", path, ": ", paste(missing_columns, collapse = ", "))
  }

  trait_data %>%
    transmute(
      snp = str_remove(as.character(SNP), ":[^:]+:[^:]+$"),
      ref_trait = toupper(as.character(A2)), alt_trait = toupper(as.character(A1)),
      beta_trait = as.numeric(b), varbeta_trait = as.numeric(se)^2,
      pvalue_trait = as.numeric(p)
    ) %>%
    filter(!is.na(snp), !is.na(beta_trait), !is.na(varbeta_trait), varbeta_trait > 0) %>%
    distinct()
}

load_selected_pairs <- function(selection_file) {
  if (is.null(selection_file) || length(selection_file) == 0 ||
      is.na(selection_file) || selection_file == "") return(NULL)
  if (!file.exists(selection_file)) stop("Selection file not found: ", selection_file)

  selected_pairs <- fread(selection_file, showProgress = FALSE)
  missing_columns <- setdiff(c("exposure", "outcome"), names(selected_pairs))
  if (length(missing_columns) > 0) {
    stop("Missing selection-file columns: ", paste(missing_columns, collapse = ", "))
  }
  unique(selected_pairs[, .(exposure = as.character(exposure), outcome = as.character(outcome))])
}

# ==============================================================================
# 2. Align trait alleles to protein alleles
# ==============================================================================
align_trait_to_protein <- function(protein_data, trait_data) {
  merged <- inner_join(protein_data, trait_data, by = "snp")
  if (nrow(merged) == 0) return(merged)

  merged %>%
    mutate(
      ref_trait_comp = complement_allele(ref_trait),
      alt_trait_comp = complement_allele(alt_trait),
      alignment_status = case_when(
        (ref_protein == ref_trait & alt_protein == alt_trait) |
          (ref_protein == ref_trait_comp & alt_protein == alt_trait_comp) ~ "Same",
        (ref_protein == alt_trait & alt_protein == ref_trait) |
          (ref_protein == alt_trait_comp & alt_protein == ref_trait_comp) ~ "Flipped",
        TRUE ~ "Mismatched"
      ),
      beta_trait_aligned = case_when(
        alignment_status == "Same" ~ beta_trait,
        alignment_status == "Flipped" ~ -beta_trait,
        TRUE ~ NA_real_
      ),
      variant_id = paste(snp, ref_protein, alt_protein, sep = ":")
    ) %>%
    filter(alignment_status != "Mismatched", !is.na(beta_trait_aligned)) %>%
    distinct(protein_id, variant_id, .keep_all = TRUE)
}

# ==============================================================================
# 3. Colocalization analysis
# ==============================================================================
run_coloc_for_trait <- function(aligned_data, trait_name, trait_type, n_total, n_case, n_protein) {
  protein_ids <- unique(aligned_data$protein_id)

  coloc_results <- lapply(protein_ids, function(protein_id_i) {
    data_i <- aligned_data %>% filter(protein_id == protein_id_i) %>%
      distinct(variant_id, .keep_all = TRUE)
    if (nrow(data_i) < 2) return(NULL)

    dataset_protein <- list(
      snp = data_i$variant_id, beta = data_i$beta_protein,
      varbeta = data_i$varbeta_protein, MAF = data_i$maf,
      type = "quant", N = n_protein
    )

    if (trait_type == "quant") {
      dataset_trait <- list(
        snp = data_i$variant_id, beta = data_i$beta_trait_aligned,
        varbeta = data_i$varbeta_trait, MAF = data_i$maf,
        type = "quant", N = n_total
      )
    } else if (trait_type == "cc") {
      if (is.na(n_case)) stop("Case count is required for case-control trait: ", trait_name)
      dataset_trait <- list(
        snp = data_i$variant_id, beta = data_i$beta_trait_aligned,
        varbeta = data_i$varbeta_trait, MAF = data_i$maf,
        type = "cc", N = n_total, s = n_case / n_total
      )
    } else {
      stop("Unsupported trait type for ", trait_name, ": ", trait_type)
    }

    coloc_result <- tryCatch(
      coloc.abf(dataset1 = dataset_protein, dataset2 = dataset_trait),
      error = function(e) {
        message("[WARN] Colocalization failed: ", protein_id_i, " vs ", trait_name,
                " | ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(coloc_result)) return(NULL)

    summary_result <- as.data.frame(t(coloc_result$summary))
    summary_result$protein_id <- protein_id_i
    summary_result$trait <- trait_name
    summary_result$n_variants <- nrow(data_i)
    summary_result
  })

  bind_rows(coloc_results)
}

# ==============================================================================
# 4. Main analysis
# ==============================================================================
protein_full <- read_protein_data(PROTEIN_INPUT_DIR)
selected_pairs <- load_selected_pairs(SELECTION_FILE)
all_results <- list()

for (i in seq_len(nrow(trait_config))) {
  trait_name <- trait_config$trait[i]
  trait_file <- trait_config$input_file[i]
  trait_type <- trait_config$type[i]
  n_total <- as.numeric(trait_config$n_total[i])
  n_case <- as.numeric(trait_config$n_case[i])

  message("[", i, "/", nrow(trait_config), "] Trait: ", trait_name)
  if (!file.exists(trait_file)) {
    message("[WARN] Trait file not found: ", trait_file)
    next
  }

  protein_data <- protein_full
  if (!is.null(selected_pairs)) {
    selected_proteins <- selected_pairs[outcome == trait_name, unique(exposure)]
    if (length(selected_proteins) == 0) {
      message("[WARN] No selected proteins for trait: ", trait_name)
      next
    }
    protein_data <- protein_data %>% filter(protein_id %in% selected_proteins)
  }
  if (nrow(protein_data) == 0) {
    message("[WARN] No protein data available for trait: ", trait_name)
    next
  }

  trait_data <- tryCatch(read_trait_data(trait_file), error = function(e) {
    message("[WARN] Failed to read trait ", trait_name, ": ", conditionMessage(e))
    NULL
  })
  if (is.null(trait_data)) next

  aligned_data <- align_trait_to_protein(protein_data = protein_data, trait_data = trait_data)
  if (nrow(aligned_data) == 0) {
    message("[WARN] No aligned variants for trait: ", trait_name)
    next
  }

  coloc_results <- run_coloc_for_trait(
    aligned_data = aligned_data, trait_name = trait_name, trait_type = trait_type,
    n_total = n_total, n_case = n_case, n_protein = N_PROTEIN
  )
  if (nrow(coloc_results) == 0) {
    message("[WARN] No colocalization results for trait: ", trait_name)
    next
  }

  trait_output <- file.path(RESULT_DIR, paste0(trait_name, "_coloc_results.csv"))
  fwrite(as.data.table(coloc_results), trait_output)
  all_results[[trait_name]] <- coloc_results
}

# ==============================================================================
# 5. Save combined results
# ==============================================================================
if (length(all_results) > 0) {
  combined_results <- bind_rows(all_results)
  combined_output <- file.path(RESULT_DIR, "all_colocalization_results.csv")
  fwrite(as.data.table(combined_results), combined_output)
  message("[DONE] Colocalization results saved to: ", RESULT_DIR)
} else {
  message("[WARN] No colocalization results were generated.")
}
