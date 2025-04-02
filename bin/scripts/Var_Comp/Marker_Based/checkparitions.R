#!/usr/bin/env Rscript

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
path_SNP <- args[1]          # File containing SNP and BOA batch file paths for the group
path_BOA <- args[2]   
null_model_file <- args[3]     # Null model RData file
my_functions_R <- args[4] 

# path_SNP <- "/blue/mateescu/gzayas97/Gabe_Thesis/5.Breed_Specific_GWAS_V2/results/Thermo_Brangus/SWA/Marker_GWAS/snp_data.fst"
# path_BOA <- "/blue/mateescu/gzayas97/Gabe_Thesis/5.Breed_Specific_GWAS_V2/results/Thermo_Brangus/SWA/Marker_GWAS/boa_data.fst"
# null_model_file <- "/blue/mateescu/gzayas97/Gabe_Thesis/5.Breed_Specific_GWAS_V2/results/Thermo_Brangus/SWA/Marker_GWAS/Thermo_Brangus_null_model.Rdata"
# my_functions_R <-"/blue/mateescu/gzayas97/Gabe_Thesis/5.Breed_Specific_GWAS_V2/bin/scripts/my_functions.R"

i <- 1
source(my_functions_R)
# Measure time for performance monitoring
start_time <- proc.time()
# Load the G matrix and null model with combined data
cat("Loading G matrix, null model, and combined data...\n")
load(null_model_file)

combined_data$SNP <- read_fst(path_SNP ,columns=  as.character(i))
combined_data$BOA <- read_fst(path_BOA , columns=  as.character(i))
# Impute missing values in SNPs
combined_data$SNP[is.na(combined_data$SNP)] <- mean(combined_data$SNP, na.rm = TRUE)

# Ensure BOA is coded as a factor
combined_data$BOA <- as.factor(combined_data$BOA)


# Fit models with SNP and BOA effects
X_snp <- model.matrix(~ as.factor(CG) + SNP, data = combined_data)

# Fit models with SNP and BOA effects
X_snp <- model.matrix(~ as.factor(CG) + SNP, data = combined_data)
snp_model <- lmm.diago(Y = combined_data$phenotype, X = X_snp, eigenK = eiK , verbose = FALSE)
# Print the time taken for processing
X_boa <- model.matrix(~ as.factor(CG) + as.factor(BOA), data = combined_data)
boa_model <- lmm.diago(Y = combined_data$phenotype, X = X_boa, eigenK = eiK , verbose = FALSE)
# Compare models to null model
boa_vs_null_pval <- as.numeric(calculate_joint_wald_test(boa_model, colnames(X_boa ), "BOA*", "SNP*"))
snp_vs_null_pval  <- as.numeric(calculate_joint_wald_test(snp_model, colnames(X_snp ), "SNP*","BOA*"))
# Selected parameters for Wald tests
selected_parameters <- c("SNP","as.factor(BOA)1",
                              "as.factor(BOA)2","SNP:as.factor(BOA)1",
                              "as.factor(BOA)2","SNP:as.factor(BOA)2")
theta_starting <- c(snp_model$sigma2, snp_model$tau)
# Fit full model with SNP and BOA interaction
X_full <- model.matrix(~ as.factor(CG) + SNP * as.factor(BOA), data = combined_data)

# Check multicollinearity with the condition number
condition_number <- kappa(X_full)
cat("Condition number for full model: ", condition_number, "\n")

if (condition_number > 1000) {
  cat("High multicollinearity detected in full model. Removing interaction term...\n")
  
  # Fit a reduced model without interaction term
  X_reduced <- model.matrix(~ as.factor(CG) + SNP + as.factor(BOA), data = combined_data)
  condition_number_reduced <- kappa(X_reduced)
  cat("Condition number for reduced model: ", condition_number_reduced, "\n")
  
  if (condition_number_reduced > 1000) {
    cat("High multicollinearity still present. Printing null results for SNP: ", names(snp_partition)[1], "\n")
    print(as.character(i), SNP_P_value = NA, BOA_P_value = NA,
                     Interaction_P_value = NA, Combined_P_value = NA,
                     SNP_vs_Null_P_value = snp_vs_null_pval, BOA_vs_Null_P_value = boa_vs_null_pval)
  } else {
    # Fit the reduced model if multicollinearity is acceptable
    reduced_model <- lmm.diago(Y = combined_data$phenotype, X = X_reduced, eigenK = eiK , verbose = FALSE)
    
    # Extract p-values for SNP and BOA without interaction
    parameter_names <- colnames(X_reduced)
    snp_pval <- as.numeric(calculate_joint_wald_test(reduced_model, parameter_names, "SNP*", "BOA*"))
    boa_pval <- as.numeric(calculate_joint_wald_test(reduced_model, parameter_names, "BOA*", "SNP*"))
    combined_pval <- as.numeric(calculate_joint_wald_test(reduced_model, parameter_names , c("BOA*","SNP*")))
    
    print(data.table(SNP_ID = as.character(i) , SNP_P_value = snp_pval, BOA_P_value = boa_pval,
                     Interaction_P_value = NA, Combined_P_value = combined_pval,
                     SNP_vs_Null_P_value = snp_vs_null_pval, BOA_vs_Null_P_value = boa_vs_null_pval))
  }
} else {
  print("Doing Normal Model")
  # Fit the full model if multicollinearity is acceptable
  full_model <- lmm.diago(Y = combined_data$phenotype, X = X_full, eigenK = eiK , verbose = FALSE)
    # Extract p-values for SNP, BOA, and interaction
  parameter_names <- colnames(X_full)
  # Perform comparison between the null and full model to get the combined p-value
  snp_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "SNP*","BOA*"))
  boa_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "BOA*","SNP*"))
  interaction_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "SNP:as.factor(BOA)*"))
  combined_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names , c("BOA*","SNP*")))

  print(data.table(SNP_ID = as.character(i), SNP_P_value = snp_pval, BOA_P_value = boa_pval,
                   Interaction_P_value = interaction_pval, Combined_P_value = combined_pval,
                   SNP_vs_Null_P_value = snp_vs_null_pval, BOA_vs_Null_P_value = boa_vs_null_pval))
}

end_time <- proc.time() - start_time
cat("Processed SNP in " , end_time[3], "seconds\n")