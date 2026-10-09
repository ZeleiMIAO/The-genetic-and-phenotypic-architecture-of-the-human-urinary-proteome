#!/usr/bin/env Rscript

# ==============================================================================
# Clinical Feature Variance Partitioning and UMAP Analysis
#
# Description:
#   This script estimates the contributions of clinical features to urinary
#   protein abundance using LMG decomposition and performs UMAP analysis
#   based on clinical, lifestyle, and genetic feature contributions.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(tibble)
  library(relaimpo)
  library(uwot)
})

# ==============================================================================
# 1. Input data
# ==============================================================================

feature <- read_xlsx("<PROJECT_DIR>/label_for_each_protein/all_three_cohorts_rerun/result_lasso_feature_summary.xlsx")
load("<PROJECT_DIR>/treated_data/discovery_all_conti_irnt.RData")

feature1 <- feature[!is.na(feature$Clinical_features), ]

feature_list <- feature1 %>%
  mutate(Feature = strsplit(Clinical_features, ";\\s*")) %>%
  unnest(Feature) %>%
  mutate(Feature = trimws(Feature)) %>%
  distinct(Protein, Feature)

df <- protein_micro_pheno

df$Gender <- as.factor(df$Gender)
df$antihypertensive_drugs <- as.factor(df$antihypertensive_drugs)
df$glucose_lowering_drugs <- as.factor(df$glucose_lowering_drugs)
df$lipid_lowering_drugs <- as.factor(df$lipid_lowering_drugs)


# ==============================================================================
# 2. Clinical feature variance partitioning
# ==============================================================================

r2_decomposition_list <- list()

protein_tasks <- feature_list %>%
  group_by(Protein) %>%
  summarise(Features = list(Feature))

for (i in 1:nrow(protein_tasks)) {

  target_prot <- protein_tasks$Protein[i]
  target_vars <- unlist(protein_tasks$Features[i])
  available_vars <- target_vars[target_vars %in% colnames(df)]

  if (length(available_vars) == 0 || !(target_prot %in% colnames(df))) next

  formula_str <- paste(target_prot, "~", paste(available_vars, collapse = " + "))

  try({
    fit <- lm(as.formula(formula_str), data = df, na.action = na.exclude)

    adj_r2_total <- summary(fit)$adj.r.squared
    if (adj_r2_total <= 0) next

    # LMG relative importance scaled by adjusted R²
    rel_imp <- calc.relimp(fit, type = "lmg", rela = TRUE)
    contribution <- rel_imp$lmg * adj_r2_total

    r2_decomposition_list[[target_prot]] <- data.frame(
      Protein = target_prot,
      Feature = names(contribution),
      Contribution_AdjR2 = as.numeric(contribution),
      Total_AdjR2 = adj_r2_total,
      stringsAsFactors = FALSE
    )
  }, silent = TRUE)

  if (i %% 100 == 0) message("Processed ", i, " proteins...")
}

final_r2_results <- bind_rows(r2_decomposition_list)


# ==============================================================================
# 3. Construct feature contribution matrix
# ==============================================================================

clinical_matrix <- final_r2_results %>%
  dplyr::select(Protein, Feature, Contribution_AdjR2) %>%
  pivot_wider(names_from = Feature, values_from = Contribution_AdjR2, values_fill = 0)

final_feature_matrix <- clinical_matrix %>%
  left_join(lifestyle_r2_table, by = "Protein") %>%
  left_join(genetics_r2_table, by = "Protein") %>%
  column_to_rownames("Protein") %>%
  as.matrix()


# ==============================================================================
# 4. UMAP analysis
# ==============================================================================

scaled_matrix <- scale(final_feature_matrix)

final_umap_res <- umap(
  scaled_matrix,
  n_neighbors = 80,
  min_dist = 0.5,
  metric = "euclidean"
)


# ==============================================================================
# 5. Save results
# ==============================================================================

result_dir <- "<PROJECT_DIR>/results"
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(final_r2_results, file.path(result_dir, "variancePartition_clinical.csv"), row.names = FALSE)
write.csv(clinical_matrix, file.path(result_dir, "clinical_allr2.csv"), row.names = FALSE)
saveRDS(final_umap_res, file.path(result_dir, "clinical_umap.rds"))