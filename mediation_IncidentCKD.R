#!/usr/bin/env Rscript

# ==============================================================================
# Urinary proteins and incident CKD
#
# Description:
#   This script performs the core analyses evaluating relationships between
#   modifiable-factor-driven urinary proteins and incident chronic kidney
#   disease (CKD).
#
#   Analyses include:
#
#   1. Determinant-domain enrichment:
#      Fisher's exact tests are used to evaluate whether urinary proteins
#      associated with incident CKD are enriched within specific modifiable
#      determinant domains.
#
#   2. Mediation analysis:
#      Tests whether blood pressure-related urinary proteins mediate
#      associations of systolic blood pressure (SBP), diastolic blood pressure
#      (DBP), or hypertension with incident CKD.
#
#   3. Hypertension-protein association:
#      Evaluates associations between hypertension and individual candidate
#      urinary proteins.
#
# ==============================================================================

rm(list = ls())

# ==============================================================================
# Packages
# ==============================================================================

suppressPackageStartupMessages({ library(data.table) library(dplyr) library(mediation) })

# ==============================================================================
# Analysis settings
# ==============================================================================

FDR_THRESHOLD <- 0.05

RANDOM_SEED <- 123

MEDIATION_BOOTSTRAP_SIMULATIONS <- 10000

modifiable <- c( "Anthropometric", "Blood pressure", "Glyaemic", "Lifestyle/Diet", "Lipids", "Renal", "Liver"
)

# ==============================================================================
# 1. Modifiable determinant-domain enrichment
# ==============================================================================

# Recalculate FDR among modifiable-factor-driven proteins using Model 1
# incident CKD association P values.

ckd_modifiable$q_new <- p.adjust( ckd_modifiable$p_value_model1, method = "BH" )

# ------------------------------------------------------------------------------
# Fisher's exact test
# ------------------------------------------------------------------------------

calculate_fisher <- function( data, domain_list ) {
  
  data_sub <- data %>%
    filter( Top_cat %in% domain_list, !is.na(q_new) ) %>%
    mutate( CKD_sig = q_new < FDR_THRESHOLD )
  
  results <- data.frame()
  
  for (d in domain_list) {
    
    # 2 x 2 table:
    #
    # Current determinant domain vs all other eligible determinant domains
    # by incident CKD association status.
    
    tab <- matrix( c( sum( data_sub$Top_cat == d & data_sub$CKD_sig ),
                      
                      sum( data_sub$Top_cat == d & !data_sub$CKD_sig ),
                      
                      sum( data_sub$Top_cat != d & data_sub$CKD_sig ),
                      
                      sum( data_sub$Top_cat != d & !data_sub$CKD_sig )
    ),
    
    nrow = 2,
    byrow = TRUE,
    
    dimnames = list( Group = c( d, "Other domains" ),
                     
                     CKD_association = c( "FDR < 0.05", "FDR >= 0.05" )
    )
    )
    
    ft <- fisher.test( tab )
    
    results <- rbind( results,
                      
                      data.frame( Domain = d,
                                  
                                  n_total = sum( data_sub$Top_cat == d ),
                                  
                                  n_CKD_associated = sum( data_sub$Top_cat == d & data_sub$CKD_sig ),
                                  
                                  proportion_CKD = mean( data_sub$CKD_sig[ data_sub$Top_cat == d ] ),
                                  
                                  OR = unname( ft$estimate ),
                                  
                                  CI_low = ft$conf.int[1],
                                  
                                  CI_high = ft$conf.int[2],
                                  
                                  P_value = ft$p.value
                      )
    )
  }
  
  results %>%
    mutate( FDR_fisher = p.adjust( P_value, method = "BH" ) ) %>%
    arrange( FDR_fisher )
}

# ------------------------------------------------------------------------------
# Seven-domain enrichment analysis
# ------------------------------------------------------------------------------

result_7 <- calculate_fisher( ckd_modifiable, modifiable )

# ------------------------------------------------------------------------------
# Six-domain non-renal enrichment analysis
#
# Renal proteins are excluded from the comparison background.
# ------------------------------------------------------------------------------

result_6 <- calculate_fisher( ckd_modifiable, setdiff( modifiable, "Renal" ) )

# ------------------------------------------------------------------------------
# Blood pressure-domain CKD-associated proteins used for mediation analyses
# ------------------------------------------------------------------------------

ckd_cardio <- ckd_modifiable[ ckd_modifiable$Top_cat == "Blood pressure" &
                                ckd_modifiable$q_new < FDR_THRESHOLD,
]

protein_ids <- ckd_cardio$protein_id

# ==============================================================================
# 2. SBP and DBP mediation analyses
# ==============================================================================

results <- data.frame()

for (i in seq_along(protein_ids)) {
  
  protein_name <- protein_ids[i]
  
  message( "[", i, "/", length(protein_ids), "] ", protein_name )
  
  temp_data <- gnhs1_pair_ck0_2[ , c( "SBP", "DBP", protein_name, "CKD", "age", "sex", "eGFR_new" ) ] %>%
    as.data.frame() %>%
    na.omit()
  
  # ---------------------------------------------------------------------------
  # SBP -> urinary protein -> incident CKD
  # ---------------------------------------------------------------------------
  
  tryCatch({
    
    med_fit <- lm( as.formula( paste( protein_name, "~ SBP + age + sex + eGFR_new" ) ), data = temp_data )
    
    y_fit <- glm( as.formula( paste( "CKD ~ SBP +", protein_name, "+ age + sex + eGFR_new" ) ),
                  data = temp_data,
                  family = poisson( link = "log" )
    )
    
    set.seed( RANDOM_SEED )
    
    mediate_result <- mediate( med_fit, y_fit, treat = "SBP", mediator = protein_name, boot = TRUE,
                               sims = MEDIATION_BOOTSTRAP_SIMULATIONS
    )
    
    temp_result <- data.frame( protein = protein_name, exposure = "SBP",
                               
                               ACME = mediate_result$d0,
                               ACME_ci_lower = mediate_result$d0.ci[1],
                               ACME_ci_upper = mediate_result$d0.ci[2],
                               ACME_p = mediate_result$d0.p,
                               
                               ADE = mediate_result$z0,
                               ADE_ci_lower = mediate_result$z0.ci[1],
                               ADE_ci_upper = mediate_result$z0.ci[2],
                               ADE_p = mediate_result$z0.p,
                               
                               Total_Effect = mediate_result$tau.coef,
                               Total_Effect_p = mediate_result$tau.p,
                               
                               Prop_Mediated = mediate_result$n0,
                               Prop_Mediated_ci_lower = mediate_result$n0.ci[1],
                               Prop_Mediated_ci_upper = mediate_result$n0.ci[2],
                               Prop_Mediated_p = mediate_result$n0.p,
                               
                               n = nrow( temp_data )
    )
    
    results <- bind_rows( results, temp_result )
    
  }, error = function(e) {
    
    message( "[WARN] SBP mediation failed for ", protein_name, ": ", conditionMessage(e) )
  })
  
  # ---------------------------------------------------------------------------
  # DBP -> urinary protein -> incident CKD
  # ---------------------------------------------------------------------------
  
  tryCatch({
    
    med_fit_dbp <- lm( as.formula( paste( protein_name, "~ DBP + age + sex + eGFR_new" ) ), data = temp_data )
    
    y_fit_dbp <- glm( as.formula( paste( "CKD ~ DBP +", protein_name, "+ age + sex + eGFR_new" ) ),
                      data = temp_data,
                      family = poisson( link = "log" )
    )
    
    set.seed( RANDOM_SEED )
    
    mediate_result_dbp <- mediate( med_fit_dbp, y_fit_dbp, treat = "DBP", mediator = protein_name,
                                   boot = TRUE,
                                   sims = MEDIATION_BOOTSTRAP_SIMULATIONS
    )
    
    temp_result_dbp <- data.frame( protein = protein_name, exposure = "DBP",
                                   
                                   ACME = mediate_result_dbp$d0,
                                   ACME_ci_lower = mediate_result_dbp$d0.ci[1],
                                   ACME_ci_upper = mediate_result_dbp$d0.ci[2],
                                   ACME_p = mediate_result_dbp$d0.p,
                                   
                                   ADE = mediate_result_dbp$z0,
                                   ADE_ci_lower = mediate_result_dbp$z0.ci[1],
                                   ADE_ci_upper = mediate_result_dbp$z0.ci[2],
                                   ADE_p = mediate_result_dbp$z0.p,
                                   
                                   Total_Effect = mediate_result_dbp$tau.coef,
                                   Total_Effect_p = mediate_result_dbp$tau.p,
                                   
                                   Prop_Mediated = mediate_result_dbp$n0,
                                   Prop_Mediated_ci_lower = mediate_result_dbp$n0.ci[1],
                                   Prop_Mediated_ci_upper = mediate_result_dbp$n0.ci[2],
                                   Prop_Mediated_p = mediate_result_dbp$n0.p,
                                   
                                   n = nrow( temp_data )
    )
    
    results <- bind_rows( results, temp_result_dbp )
    
  }, error = function(e) {
    
    message( "[WARN] DBP mediation failed for ", protein_name, ": ", conditionMessage(e) )
  })
}

# ==============================================================================
# 3. Hypertension mediation analysis
# ==============================================================================

results1 <- data.frame()

for (i in seq_along(protein_ids)) {
  
  protein_name <- protein_ids[i]
  
  message( "[", i, "/", length(protein_ids), "] ", protein_name )
  
  temp_data <- gnhs1_pair_ck0_3[ , c( "hyper", protein_name, "CKD", "age", "sex", "eGFR_new" ) ] %>%
    as.data.frame() %>%
    na.omit()
  
  tryCatch({
    
    # -------------------------------------------------------------------------
    # Hypertension -> urinary protein
    # -------------------------------------------------------------------------
    
    med_fit <- lm( as.formula( paste( protein_name, "~ hyper + age + sex + eGFR_new" ) ), data = temp_data )
    
    hyper_beta <- coef( med_fit )[["hyper"]]
    
    hyper_ci <- confint( med_fit, parm = "hyper", level = 0.95 )
    
    # -------------------------------------------------------------------------
    # Hypertension + urinary protein -> incident CKD
    # -------------------------------------------------------------------------
    
    y_fit <- glm( as.formula( paste( "CKD ~ hyper +", protein_name, "+ age + sex + eGFR_new" ) ),
                  data = temp_data,
                  family = poisson( link = "log" )
    )
    
    set.seed( RANDOM_SEED )
    
    mediate_result <- mediate( med_fit, y_fit, treat = "hyper", mediator = protein_name, boot = TRUE,
                               sims = MEDIATION_BOOTSTRAP_SIMULATIONS
    )
    
    temp_result <- data.frame( protein = protein_name, exposure = "hypertension",
                               
                               hyper_protein_beta = hyper_beta,
                               hyper_protein_ci_lower = hyper_ci[1],
                               hyper_protein_ci_upper = hyper_ci[2],
                               
                               ACME = mediate_result$d0,
                               ACME_ci_lower = mediate_result$d0.ci[1],
                               ACME_ci_upper = mediate_result$d0.ci[2],
                               ACME_p = mediate_result$d0.p,
                               
                               ADE = mediate_result$z0,
                               ADE_ci_lower = mediate_result$z0.ci[1],
                               ADE_ci_upper = mediate_result$z0.ci[2],
                               ADE_p = mediate_result$z0.p,
                               
                               Total_Effect = mediate_result$tau.coef,
                               Total_Effect_p = mediate_result$tau.p,
                               
                               Prop_Mediated = mediate_result$n0,
                               Prop_Mediated_ci_lower = mediate_result$n0.ci[1],
                               Prop_Mediated_ci_upper = mediate_result$n0.ci[2],
                               Prop_Mediated_p = mediate_result$n0.p,
                               
                               n = nrow( temp_data )
    )
    
    results1 <- bind_rows( results1, temp_result )
    
  }, error = function(e) {
    
    message( "[WARN] Hypertension mediation failed for ", protein_name, ": ", conditionMessage(e) )
  })
}

# ==============================================================================
# 3.1 Protein annotation and multiple-testing correction
# ==============================================================================

# ------------------------------------------------------------------------------
# SBP / DBP mediation results
# ------------------------------------------------------------------------------

results <- merge( gene, results, by.x = "protein_id", by.y = "protein" )

results$acme_fdr <- p.adjust( results$ACME_p, method = "BH" )

results$prop_fdr <- p.adjust( results$Prop_Mediated_p, method = "BH" )

# ------------------------------------------------------------------------------
# Hypertension mediation results
# ------------------------------------------------------------------------------

results1 <- merge( gene, results1, by.x = "protein_id", by.y = "protein" )

results1$acme_fdr <- p.adjust( results1$ACME_p, method = "BH" )

results1$prop_fdr <- p.adjust( results1$Prop_Mediated_p, method = "BH" )

# Proteins with significant mediated proportions.

results1_sig <- results1[ results1$prop_fdr < FDR_THRESHOLD, ]

# ==============================================================================
# 4. Supplementary hypertension-protein association
# ==============================================================================

results2 <- data.frame()

for (i in seq_along(protein_ids)) {
  
  protein_name <- protein_ids[i]
  
  temp_data <- gnhs1_pair_ck0_3[ , c( "hyper", protein_name, "age", "sex", "eGFR_new" ) ] %>%
    as.data.frame() %>%
    na.omit()
  
  tryCatch({
    
    fit_protein <- lm( as.formula( paste( protein_name, "~ hyper + age + sex + eGFR_new" ) ), data = temp_data
    )
    
    hyper_beta <- coef( fit_protein )[["hyper"]]
    
    hyper_ci <- confint( fit_protein, parm = "hyper", level = 0.95 )
    
    hyper_p <- coef( summary( fit_protein ) )[ "hyper", "Pr(>|t|)" ]
    
    temp_result <- data.frame( protein = protein_name, exposure = "hypertension",
                               
                               hyper_protein_beta = hyper_beta,
                               hyper_protein_ci_lower = hyper_ci[1],
                               hyper_protein_ci_upper = hyper_ci[2],
                               hyper_protein_p = hyper_p,
                               
                               n = nrow( temp_data )
    )
    
    results2 <- bind_rows( results2, temp_result )
    
  }, error = function(e) {
    
    message( "[WARN] Hypertension-protein model failed for ", protein_name, ": ", conditionMessage(e) )
  })
}

results2 <- merge( gene, results2, by.x = "protein_id", by.y = "protein" )


# ==============================================================================
# 5. Save results
# ==============================================================================

OUTPUT_DIR <- "<PROJECT_DIR>/results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(result_7), file.path(OUTPUT_DIR, "enrichment_7domains.csv"))
fwrite(as.data.table(result_6), file.path(OUTPUT_DIR, "enrichment_6domains.csv"))
fwrite(as.data.table(results), file.path(OUTPUT_DIR, "mediation_SBP_DBP.csv"))
fwrite(as.data.table(results1), file.path(OUTPUT_DIR, "mediation_hypertension.csv"))
fwrite(as.data.table(results1_sig), file.path(OUTPUT_DIR, "mediation_hypertension_significant.csv"))
fwrite(as.data.table(results2), file.path(OUTPUT_DIR, "hypertension_protein_association.csv"))
