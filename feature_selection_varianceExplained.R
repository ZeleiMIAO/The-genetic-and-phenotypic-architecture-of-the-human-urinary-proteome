# ==============================================================================
# Protein Variance Analysis
#
# Description:
# This script estimates the proportion of variance in individual protein
# levels explained by clinical factors, lifestyle factors, and sentinel SNPs.
#
# ==============================================================================

## Set working directory
setwd("/path/to/project")

options(stringsAsFactors = FALSE)

rm(list = ls())

# Load required packages
# library(RNOmni)
library(caret)
library(glmnet)
library(doMC)
library(tidyr)
library(dplyr)
library(plyr)
library(openxlsx)
library(foreach)

# ==============================================================================
# Load phenotype and protein data
# ==============================================================================

pro_pheno_dis <- read.xlsx( "./data/dis_pro_pheno.xlsx" )

pro_pheno_internal <- read.xlsx( "./data/internal_pro_pheno.xlsx" )

pro_pheno_independent <- read.xlsx( "./data/independent_pro_pheno.xlsx" )

pro_pheno_external <- read.xlsx( "./data/external_pro_pheno.xlsx" )

# ==============================================================================
# Load and merge SNP data
# ==============================================================================

pro_snp_dis <- read.xlsx( "./data/dis_pro_snp.xlsx" )

pro_snp_dis <- merge( pro_pheno_dis[ , c("SampleID", "survey_year", "season", "Followup") ], pro_snp_dis,
  by = "SampleID", all = FALSE )

pro_snp_internal <- read.xlsx( "./data/internal_pro_snp.xlsx" )

pro_snp_internal <- merge( pro_pheno_internal[ , c("Id", "survey_year", "season", "Followup") ], pro_snp_internal,
  by = "Id", all = FALSE )

pro_snp_independent <- read.xlsx( "./data/independent_pro_snp.xlsx" )

pro_snp_independent <- merge( pro_pheno_independent[ , c("Id", "survey_year", "season", "Followup") ],
  pro_snp_independent, by = "Id", all = FALSE )

pro_snp_external <- read.xlsx( "./data/external_pro_snp.xlsx" )

pro_snp_external <- merge( pro_pheno_external[ , c("Id", "survey_year", "season", "Followup") ], pro_snp_external,
  by = "Id", all = FALSE )

# ==============================================================================
# Register parallel backend
# ==============================================================================

registerDoMC(32)

# ==============================================================================
# Define analysis datasets
# ==============================================================================

dat <- pro_pheno_dis
val <- pro_pheno_internal
val2 <- pro_pheno_independent
val_external <- pro_pheno_external

dat_snp <- pro_snp_dis
val_snp <- pro_snp_internal
val_snp2 <- pro_snp_independent
val_snp_external <- pro_snp_external

# ==============================================================================
# Define covariates
# ==============================================================================

# category <- c( "survey_year", "season", "Followup", "Gender", "antihypertensive_drugs", "glucose_lowering_drugs",
# "lipid_lowering_drugs", "Alcohol", "Smoke", "Tea", "MultiVit" )
#
# continuous <- c( "Age", "BMI", "WHR", "HDL", "LDL", "HbA1C", "TC", "TG", "Glu", "eGFR", "SBP", "DBP", "ALT", "AST",
# "ALP", "MET", "Sleep_hr", "Water_ml", "SaltPref", "Energy_kcal", "Carb_g", "Protein_g", "Fat_g",
# "RefinedGrain_serv", "WholeGrain_serv", "Fruit_serv", "Veg_serv", "Nuts_serv", "Legume_serv", "Potato_serv",
# "Meat_serv", "Egg_serv", "Dairy_serv", "FishSea_serv", "Beverage_serv", "Sweets_serv", "Oil_g" )
#
# technical <- c( "survey_year", "season", "Followup" )

clinical <- c( "Age", "Gender", "BMI", "WHR", "HDL", "LDL", "HbA1C", "TC", "TG", "Glu", "eGFR", "SBP", "DBP",
  "ALT", "AST", "ALP", "antihypertensive_drugs", "glucose_lowering_drugs", "lipid_lowering_drugs" )

lifestyle <- c( "Alcohol", "Smoke", "MET", "Sleep_hr", "Tea", "Water_ml", "MultiVit", "SaltPref", "Energy_kcal",
  "Carb_g", "Protein_g", "Fat_g", "RefinedGrain_serv", "WholeGrain_serv", "Fruit_serv", "Veg_serv", "Nuts_serv",
  "Legume_serv", "Potato_serv", "Meat_serv", "Egg_serv", "Dairy_serv", "FishSea_serv", "Beverage_serv",
  "Sweets_serv", "Oil_g" )

protein_cols <- colnames(dat)[11:2400]

snp_cols <- colnames(pro_snp_dis)[
  (grep("PHENOTYPE", colnames(pro_snp_dis)) + 1):
  ncol(pro_snp_dis)
]

snp_summary <- read.table( "./data/all_sentinel_associations.txt", header = TRUE )

snp_summary$TRAIT <- gsub( "_gwas", "", snp_summary$TRAIT )

# ==============================================================================
# Estimating variance of individual proteins
# ==============================================================================

lasso_variance <- function( tmp, data, data_rep, data_rep2, data_external, response ) {

  tmp <- data[
    ,
    c(response, tmp),
    drop = FALSE
  ]

  x <- as.matrix(tmp[, -1])
  y <- tmp[, 1]

  n_boot <- 100
  n_vars <- ncol(x)

  coef_matrix <- matrix( 0, nrow = n_vars + 1, ncol = n_boot )

  for (b in seq_len(n_boot)) {

    boot_idx <- sample( nrow(x), replace = TRUE )

    x_boot <- x[
      boot_idx,
      ,
      drop = FALSE
    ]

    y_boot <- y[boot_idx]

    cv <- glmnet::cv.glmnet( x_boot, y_boot, alpha = 1, nfolds = 10, type.measure = "mse" )

    coefs <- as.vector( coef( cv, s = "lambda.min" ) )

    coef_matrix[, b] <- coefs != 0
  }

  freq <- rowMeans(coef_matrix)

  selected_coefs <- freq >= 0.7
  selected_idx <- which(selected_coefs)

  if (sum(selected_idx) <= 1) {

    stats <- c( 0, 0, 1, 0, 1, 0, 1, 0, 1 )

    selected_features <- character(0)

  } else {

    selected_features <- names(tmp)[selected_idx][-1]

    tmp <- tmp[
      ,
      selected_coefs,
      drop = FALSE
    ]

    fit <- summary( lm( tmp[, 1] ~ ., data.frame(tmp[, -1]) ) )

    tmp_rep <- data_rep[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_rep <- summary( lm( tmp_rep[, 1] ~ ., data.frame(tmp_rep[, -1]) ) )

    tmp_rep2 <- data_rep2[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_rep2 <- summary( lm( tmp_rep2[, 1] ~ ., data.frame(tmp_rep2[, -1]) ) )

    tmp_external <- data_external[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_external <- summary( lm( tmp_external[, 1] ~ ., data.frame(tmp_external[, -1]) ) )

    stats <- c( ncol(tmp) - 1,
      fit$adj.r.squared, pf( fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3], lower.tail = FALSE ),
      fit_rep$adj.r.squared, pf( fit_rep$fstatistic[1], fit_rep$fstatistic[2], fit_rep$fstatistic[3], lower.tail = FALSE
      ),
      fit_rep2$adj.r.squared, pf( fit_rep2$fstatistic[1], fit_rep2$fstatistic[2], fit_rep2$fstatistic[3],
        lower.tail = FALSE ),
      fit_external$adj.r.squared, pf( fit_external$fstatistic[1], fit_external$fstatistic[2],
        fit_external$fstatistic[3], lower.tail = FALSE ) )
  }

  return( list( stats = stats, features = selected_features ) )
}

ElasticNet_variance <- function( tmp, data, data_rep, data_rep2, data_external, response ) {

  tmp <- data[
    ,
    c(response, tmp),
    drop = FALSE
  ]

  x <- as.matrix(tmp[, -1])
  y <- tmp[, 1]

  n_boot <- 100
  n_vars <- ncol(x)

  coef_matrix <- matrix( 0, nrow = n_vars + 1, ncol = n_boot )

  for (b in seq_len(n_boot)) {

    boot_idx <- sample( nrow(x), replace = TRUE )

    x_boot <- x[
      boot_idx,
      ,
      drop = FALSE
    ]

    y_boot <- y[boot_idx]

    cv <- glmnet::cv.glmnet( x_boot, y_boot, alpha = 0.5, nfolds = 10, type.measure = "mse" )

    coefs <- as.vector( coef( cv, s = "lambda.min" ) )

    coef_matrix[, b] <- coefs != 0
  }

  freq <- rowMeans(coef_matrix)

  selected_coefs <- freq >= 0.7
  selected_idx <- which(selected_coefs)

  if (sum(selected_idx) <= 1) {

    fit <- c( 0, 0, 1, 0, 1, 0, 1, 0, 1 )

    selected_features <- character(0)

  } else {

    selected_features <- names(tmp)[selected_idx][-1]

    tmp <- tmp[
      ,
      selected_coefs,
      drop = FALSE
    ]

    fit_lm <- summary( lm( tmp[, 1] ~ ., data.frame(tmp[, -1]) ) )

    tmp_rep <- data_rep[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_rep <- summary( lm( tmp_rep[, 1] ~ ., data.frame(tmp_rep[, -1]) ) )

    tmp_rep2 <- data_rep2[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_rep2 <- summary( lm( tmp_rep2[, 1] ~ ., data.frame(tmp_rep2[, -1]) ) )

    tmp_external <- data_external[
      ,
      colnames(tmp),
      drop = FALSE
    ]

    fit_external <- summary( lm( tmp_external[, 1] ~ ., data.frame(tmp_external[, -1]) ) )

    fit <- c( ncol(tmp) - 1,
      fit_lm$adj.r.squared, pf( fit_lm$fstatistic[1], fit_lm$fstatistic[2], fit_lm$fstatistic[3], lower.tail = FALSE ),
      fit_rep$adj.r.squared, pf( fit_rep$fstatistic[1], fit_rep$fstatistic[2], fit_rep$fstatistic[3], lower.tail = FALSE
      ),
      fit_rep2$adj.r.squared, pf( fit_rep2$fstatistic[1], fit_rep2$fstatistic[2], fit_rep2$fstatistic[3],
        lower.tail = FALSE ),
      fit_external$adj.r.squared, pf( fit_external$fstatistic[1], fit_external$fstatistic[2],
        fit_external$fstatistic[3], lower.tail = FALSE ) )
  }

  return( list( stats = fit, features = selected_features ) )
}

variance_snp <- function( tmp, data, data_rep, data_rep2, data_external, response ) {

  if (length(tmp) < 1) {

    fit <- c( 0, 0, 1, 0, 1, 0, 1, 0, 1 )

  } else {

    tmp_data <- data[
      ,
      c(response, tmp),
      drop = FALSE
    ]

    fit_lm <- summary( lm( tmp_data[, 1] ~ ., data.frame(tmp_data[, -1]) ) )

    tmp_rep <- data_rep[
      ,
      colnames(tmp_data),
      drop = FALSE
    ]

    fit_rep <- summary( lm( tmp_rep[, 1] ~ ., data.frame(tmp_rep[, -1]) ) )

    tmp_rep2 <- data_rep2[
      ,
      colnames(tmp_data),
      drop = FALSE
    ]

    fit_rep2 <- summary( lm( tmp_rep2[, 1] ~ ., data.frame(tmp_rep2[, -1]) ) )

    tmp_external <- data_external[
      ,
      colnames(tmp_data),
      drop = FALSE
    ]

    fit_external <- summary( lm( tmp_external[, 1] ~ ., data.frame(tmp_external[, -1]) ) )

    fit <- c( ncol(tmp_data) - 1,
      fit_lm$adj.r.squared, pf( fit_lm$fstatistic[1], fit_lm$fstatistic[2], fit_lm$fstatistic[3], lower.tail = FALSE ),
      fit_rep$adj.r.squared, pf( fit_rep$fstatistic[1], fit_rep$fstatistic[2], fit_rep$fstatistic[3], lower.tail = FALSE
      ),
      fit_rep2$adj.r.squared, pf( fit_rep2$fstatistic[1], fit_rep2$fstatistic[2], fit_rep2$fstatistic[3],
        lower.tail = FALSE ),
      fit_external$adj.r.squared, pf( fit_external$fstatistic[1], fit_external$fstatistic[2],
        fit_external$fstatistic[3], lower.tail = FALSE ) )
  }

  return(fit)
}

# ==============================================================================
# LASSO analysis
# ==============================================================================

protein <- protein_cols

results <- foreach(
  i = seq_along(protein),
  .combine = "list",
  .multicombine = TRUE,
  .export = c(
    "protein",
    "dat",
    "val",
    "val2",
    "val_external",
    "clinical",
    "lifestyle",
    "snp_summary",
    "dat_snp",
    "val_snp",
    "val_snp2",
    "val_snp_external",
    "lasso_variance",
    "variance_snp"
  ),
  .packages = c(
    "glmnet",
    "caret",
    "tidyr",
    "dplyr",
    "plyr"
  )
) %dopar% {

  set.seed(666)

  prot <- protein[i]

  cat( "Processing protein", i, "of", length(protein), ":", prot, "\n" )

  stats_list <- list()
  features_list <- list()

  # Clinical factors
  res_clinical <- lasso_variance( tmp = clinical, data = dat, data_rep = val, data_rep2 = val2,
    data_external = val_external, response = prot )

  stats_list$clinical <- res_clinical$stats
  features_list$clinical <- res_clinical$features

  # Lifestyle factors
  res_lifestyle <- lasso_variance( tmp = lifestyle, data = dat, data_rep = val, data_rep2 = val2,
    data_external = val_external, response = prot )

  stats_list$lifestyle <- res_lifestyle$stats
  features_list$lifestyle <- res_lifestyle$features

  # Sentinel SNPs
  snp_vec <- unique( as.character( snp_summary$SNP[ snp_summary$TRAIT == prot ] ) )

  res_snp <- variance_snp( tmp = snp_vec, data = dat_snp, data_rep = val_snp, data_rep2 = val_snp2,
    data_external = val_snp_external, response = prot )

  stats_list$snp <- res_snp

  features_list$snp <- if (res_snp[1] > 0) {
    snp_vec
  } else {
    character(0)
  }

  current_row <- data.frame( protein = prot,
    n_clinical = stats_list$clinical[1], r2.adj_clinical = stats_list$clinical[2],
    r2.adj_clinical_p = stats_list$clinical[3], r2.adj_clinical_internal = stats_list$clinical[4],
    r2.adj_clinical_internal_p = stats_list$clinical[5], r2.adj_clinical_independent = stats_list$clinical[6],
    r2.adj_clinical_independent_p = stats_list$clinical[7], r2.adj_clinical_external = stats_list$clinical[8],
    r2.adj_clinical_external_p = stats_list$clinical[9],
    n_lifestyle = stats_list$lifestyle[1], r2.adj_lifestyle = stats_list$lifestyle[2],
    r2.adj_lifestyle_p = stats_list$lifestyle[3], r2.adj_lifestyle_internal = stats_list$lifestyle[4],
    r2.adj_lifestyle_internal_p = stats_list$lifestyle[5], r2.adj_lifestyle_independent = stats_list$lifestyle[6],
    r2.adj_lifestyle_independent_p = stats_list$lifestyle[7], r2.adj_lifestyle_external = stats_list$lifestyle[8],
    r2.adj_lifestyle_external_p = stats_list$lifestyle[9],
    n_snp = stats_list$snp[1], r2.adj_snp = stats_list$snp[2], r2.adj_snp_p = stats_list$snp[3],
    r2.adj_snp_internal = stats_list$snp[4], r2.adj_snp_internal_p = stats_list$snp[5],
    r2.adj_snp_independent = stats_list$snp[6], r2.adj_snp_independent_p = stats_list$snp[7],
    r2.adj_snp_external = stats_list$snp[8], r2.adj_snp_external_p = stats_list$snp[9],
    stringsAsFactors = FALSE )

  list( prot = prot, row = current_row, features = features_list )
}

result_lasso <- do.call( rbind, lapply( results, function(x) x$row ) )

selected_features_all <- setNames( lapply( results, function(x) x$features ), sapply( results, function(x) x$prot )
)

rownames(result_lasso) <- result_lasso$protein

feature_summary <- do.call( rbind, lapply( names(selected_features_all), function(p) {
      lst <- selected_features_all[[p]]
      data.frame( Protein = p, Clinical_features = paste( lst$clinical, collapse = "; " ), Lifestyle_features = paste(
          lst$lifestyle, collapse = "; " ), SNP_features = paste( lst$snp, collapse = "; " ),
        stringsAsFactors = FALSE )
    } ) )

write.xlsx( result_lasso, "./result/result_lasso_for_each_protein.xlsx" )

write.xlsx( feature_summary, "./result/result_lasso_feature_summary.xlsx" )

# ==============================================================================
# Elastic Net analysis
# ==============================================================================

protein <- protein_cols

results <- foreach(
  i = seq_along(protein),
  .combine = "list",
  .multicombine = TRUE,
  .export = c(
    "protein",
    "dat",
    "val",
    "val2",
    "val_external",
    "clinical",
    "lifestyle",
    "snp_summary",
    "dat_snp",
    "val_snp",
    "val_snp2",
    "val_snp_external",
    "ElasticNet_variance",
    "variance_snp"
  ),
  .packages = c(
    "glmnet",
    "caret",
    "tidyr",
    "dplyr",
    "plyr"
  )
) %dopar% {

  set.seed(666)

  prot <- protein[i]

  cat( "Processing protein", i, "of", length(protein), ":", prot, "\n" )

  stats_list <- list()
  features_list <- list()

  # Clinical factors
  res_clinical <- ElasticNet_variance( tmp = clinical, data = dat, data_rep = val, data_rep2 = val2,
    data_external = val_external, response = prot )

  stats_list$clinical <- res_clinical$stats
  features_list$clinical <- res_clinical$features

  # Lifestyle factors
  res_lifestyle <- ElasticNet_variance( tmp = lifestyle, data = dat, data_rep = val, data_rep2 = val2,
    data_external = val_external, response = prot )

  stats_list$lifestyle <- res_lifestyle$stats
  features_list$lifestyle <- res_lifestyle$features

  # Sentinel SNPs
  snp_vec <- unique( as.character( snp_summary$SNP[ snp_summary$TRAIT == prot ] ) )

  res_snp <- variance_snp( tmp = snp_vec, data = dat_snp, data_rep = val_snp, data_rep2 = val_snp2,
    data_external = val_snp_external, response = prot )

  stats_list$snp <- res_snp

  features_list$snp <- if (res_snp[1] > 0) {
    snp_vec
  } else {
    character(0)
  }

  current_row <- data.frame( protein = prot,
    n_clinical = stats_list$clinical[1], r2.adj_clinical = stats_list$clinical[2],
    r2.adj_clinical_p = stats_list$clinical[3], r2.adj_clinical_internal = stats_list$clinical[4],
    r2.adj_clinical_internal_p = stats_list$clinical[5], r2.adj_clinical_independent = stats_list$clinical[6],
    r2.adj_clinical_independent_p = stats_list$clinical[7], r2.adj_clinical_external = stats_list$clinical[8],
    r2.adj_clinical_external_p = stats_list$clinical[9],
    n_lifestyle = stats_list$lifestyle[1], r2.adj_lifestyle = stats_list$lifestyle[2],
    r2.adj_lifestyle_p = stats_list$lifestyle[3], r2.adj_lifestyle_internal = stats_list$lifestyle[4],
    r2.adj_lifestyle_internal_p = stats_list$lifestyle[5], r2.adj_lifestyle_independent = stats_list$lifestyle[6],
    r2.adj_lifestyle_independent_p = stats_list$lifestyle[7], r2.adj_lifestyle_external = stats_list$lifestyle[8],
    r2.adj_lifestyle_external_p = stats_list$lifestyle[9],
    n_snp = stats_list$snp[1], r2.adj_snp = stats_list$snp[2], r2.adj_snp_p = stats_list$snp[3],
    r2.adj_snp_internal = stats_list$snp[4], r2.adj_snp_internal_p = stats_list$snp[5],
    r2.adj_snp_independent = stats_list$snp[6], r2.adj_snp_independent_p = stats_list$snp[7],
    r2.adj_snp_external = stats_list$snp[8], r2.adj_snp_external_p = stats_list$snp[9],
    stringsAsFactors = FALSE )

  list( prot = prot, row = current_row, features = features_list )
}

result_elastic <- do.call( rbind, lapply( results, function(x) x$row ) )

selected_features_all <- setNames( lapply( results, function(x) x$features ), sapply( results, function(x) x$prot )
)

rownames(result_elastic) <- result_elastic$protein

feature_summary <- do.call( rbind, lapply( names(selected_features_all), function(p) {
      lst <- selected_features_all[[p]]
      data.frame( Protein = p, Clinical_features = paste( lst$clinical, collapse = "; " ), Lifestyle_features = paste(
          lst$lifestyle, collapse = "; " ), SNP_features = paste( lst$snp, collapse = "; " ),
        stringsAsFactors = FALSE )
    } ) )

write.xlsx( result_elastic, "./result/result_elastic_for_each_protein.xlsx" )

write.xlsx( feature_summary, "./result/result_elastic_feature_summary.xlsx" )
