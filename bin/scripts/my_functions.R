suppressPackageStartupMessages({
  # Load necessary libraries
  library(data.table)
  library(lme4)
  library(lattice)
  library(lme4qtl)
  library(AGHmatrix)
  library(gaston)
  # Load necessary libraries
  library(aod)
  library(ggplot2)
  library(dplyr)
  # Ensure patchwork is loaded for plot combination
  library(patchwork)
  library(arrow)
  library(fst)
  library(tidyverse)   # for data manipulation
  library(LDheatmap)   # for generating LD heatmaps
  library(genetics)    # for working with SNP data
})


# -----------------------------------
# 1) compare_models_gaston
# -----------------------------------
compare_models_gaston <- function(null_model, full_model) {
  
  # Extract log-likelihoods from both models
  loglik_null <- null_model$logL
  loglik_full <- full_model$logL
  
  # Check if the full model is worse than the null model
  if (loglik_full < loglik_null) {
    return(1)  # If full model is worse, return p-value of 1
  }
  
  # Calculate Likelihood Ratio Test statistic
  lrt_statistic <- 2 * (loglik_full - loglik_null)
  
  # Degrees of freedom difference
  df_diff <- length(full_model$BLUP_beta) - length(null_model$BLUP_beta)
  # Ensure degrees of freedom difference is valid
  if (df_diff <= 0) {
    stop("Full model must have more parameters than the null model")
  }
  # Calculate p-value for the Likelihood Ratio Test
  p_value_lrt <- pchisq(lrt_statistic, df = df_diff, lower.tail = FALSE)
  
  # Return the p-value
  return(p_value_lrt)
}


# -----------------------------------
# 2) calculate_wald_test_gaston
# -----------------------------------
calculate_wald_test_gaston <- function(full_model, parameter_names, selected_parameters) {
  
  # Extract fixed effects (coefficients) and their standard errors
  fixed_effects <- full_model$BLUP_beta
  se_fixed_effects <- sqrt(diag(full_model$varbeta))  
  # Standard errors = sqrt(diag(covariance_matrix))
  
  # Calculate Wald statistics and p-values
  wald_statistics <- fixed_effects / se_fixed_effects
  p_values_wald   <- 2 * (1 - pnorm(abs(wald_statistics)))  # Two-tailed p-values
  
  # Create a table of results
  wald_results <- data.frame(
    Parameter      = parameter_names,
    Estimate       = fixed_effects,
    Std_Error      = se_fixed_effects,
    Wald_Statistic = wald_statistics,
    P_value        = p_values_wald
  )
  
  # Select only the specified parameters
  selected_results <- wald_results[wald_results$Parameter %in% selected_parameters, ]
  rownames(selected_results) <- selected_results$Parameter
  
  # Return the selected results
  return(selected_results)
}


# -----------------------------------
# 3) calculate_joint_wald_test
# -----------------------------------
calculate_joint_wald_test <- function(full_model, parameter_names, parameter_patterns, exclude_pattern = NULL) {
  # Identify indices for the selected parameters based on multiple patterns
  selected_indices <- unlist(lapply(parameter_patterns, function(pattern) grep(pattern, parameter_names)))
  
  # If there is an exclusion pattern, remove those indices
  if (!is.null(exclude_pattern)) {
    exclude_indices <- grep(exclude_pattern, parameter_names)
    selected_indices <- setdiff(selected_indices, exclude_indices)
  }
  
  # Remove duplicate indices if any pattern overlap exists
  selected_indices <- unique(selected_indices)
  
  # Check if selected_indices is empty
  if (length(selected_indices) == 0) {
    warning("No parameters matched the given patterns.")
    return(NA)
  }
  
  # Extract relevant estimates and covariance matrix
  joint_estimates <- full_model$BLUP_beta[selected_indices]
  joint_cov       <- full_model$varbeta[selected_indices, selected_indices]
  
  # Check if joint_cov is valid
  if (is.null(joint_cov) || length(joint_cov) == 0 || any(is.na(joint_cov))) {
    stop("The covariance matrix for the selected parameters is invalid or empty.")
  }
  
  # Handle case when there is only one parameter selected
  if (length(joint_estimates) == 1) {
    wald_statistic <- (joint_estimates / sqrt(joint_cov))^2
    df <- 1
  } else {
    # Add ridge penalty to stabilize the covariance matrix
    epsilon <- 1e-7
    joint_cov_reg <- joint_cov + diag(epsilon, nrow(joint_cov))
    
    # Calculate Wald statistic for the joint parameters
    wald_statistic <- t(joint_estimates) %*% solve(joint_cov_reg) %*% joint_estimates
    df <- length(joint_estimates)
  }
  
  # Calculate the p-value from the chi-square distribution
  p_value_wald <- pchisq(wald_statistic, df = df, lower.tail = FALSE)
  
  # Return the p-value
  return(p_value_wald)
}


# -----------------------------------
# 4) create_boa_windows
# -----------------------------------
create_boa_windows <- function(snp_data, boa_data, max_window_size = 20, 
                               output_dir = "Partitions", prefix = "window") {
  
  snp_columns <- snp_data[, -(1:6), with = FALSE]  # Drop non-genotype columns
  boa_columns <- boa_data[, -(1:6), with = FALSE]  # Drop non-genotype columns
  snp_ids     <- colnames(snp_columns)             # Real SNP IDs
  
  num_markers     <- ncol(snp_columns)
  unique_windows  <- list()
  window_index    <- 1
  
  # Create output_dir if not exists
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  start_col <- 1
  while (start_col <= num_markers) {
    # Start with the maximum window size and reduce if necessary
    for (window_size in max_window_size:1) {
      end_col <- min(start_col + window_size - 1, num_markers)
      
      snp_window <- snp_columns[, start_col:end_col, with = FALSE]
      boa_window <- boa_columns[, start_col:end_col, with = FALSE]
      
      # Check if all BOA values within each row are the same across the markers in the window
      boa_unique <- all(apply(boa_window, 1, function(row) length(unique(row)) == 1))
      
      if (boa_unique) {
        # Save the SNP and BOA window
        snp_window_names <- paste0("BOA", window_index, "_SNP", seq_len(ncol(snp_window)))
        names(snp_window) <- snp_window_names
        
        boa_window_name <- paste0("BOA_window", window_index)
        
        fwrite(snp_window, file.path(output_dir, paste0("SNP_window_", window_index, ".csv")))
        fwrite(boa_window[, 1, with = FALSE], file.path(output_dir, 
               paste0(boa_window_name, ".csv")))  # Save just the first column
        
        # Create and save the SNP ID file for this window
        snp_id_file <- data.table(
          Real_SNP   = snp_ids[start_col:end_col],
          Custom_SNP = snp_window_names
        )
        fwrite(snp_id_file, file.path(output_dir, 
               paste0("SNP_IDs_window_", window_index, ".csv")))
        
        unique_windows[[window_index]] <- list(
          snp_file    = paste0("SNP_window_", window_index, ".csv"),
          boa_file    = paste0(boa_window_name, ".csv"),
          snp_id_file = paste0("SNP_IDs_window_", window_index, ".csv")
        )
        
        start_col <- end_col + 1
        break
      }
      
      # If no unique window found, reduce window size
      if (window_size == 1) {
        cat("Warning: Couldn't find a unique BOA window for markers", 
            start_col, "to", end_col, "\n")
        start_col <- end_col + 1
      }
    }
    
    window_index <- window_index + 1
  }
  
  return(unique_windows)
}


# -----------------------------------
# 5) partition_into_batches
# -----------------------------------
partition_into_batches <- function(unique_windows, num_partitions = 100, output_dir = "Partitions") {
  num_windows <- length(unique_windows)
  
  # Calculate base number of windows per partition and remainder
  windows_per_partition <- num_windows %/% num_partitions
  extra_windows         <- num_windows %% num_partitions
  
  # Create an empty list to hold partition data
  partition_data <- vector("list", num_partitions)
  for (i in 1:num_partitions) {
    partition_data[[i]] <- list(
      snp_files    = list(),
      boa_files    = list(),
      snp_id_files = list()
    )
  }
  
  # Assign windows to partitions in a round-robin fashion
  window_idx <- 1
  for (i in seq_along(unique_windows)) {
    partition_idx <- (window_idx - 1) %% num_partitions + 1
    
    if (length(partition_data[[partition_idx]]$snp_files) < windows_per_partition ||
        (partition_idx <= extra_windows && 
         length(partition_data[[partition_idx]]$snp_files) == windows_per_partition)) {
      
      partition_data[[partition_idx]]$snp_files    <- c(partition_data[[partition_idx]]$snp_files, 
                                                        unique_windows[[i]]$snp_file)
      partition_data[[partition_idx]]$boa_files    <- c(partition_data[[partition_idx]]$boa_files, 
                                                        unique_windows[[i]]$boa_file)
      partition_data[[partition_idx]]$snp_id_files <- c(partition_data[[partition_idx]]$snp_id_files, 
                                                        unique_windows[[i]]$snp_id_file)
      window_idx <- window_idx + 1
    }
  }
  
  # Save each partition's SNP and BOA files
  for (partition_idx in seq_along(partition_data)) {
    partition_dir <- file.path(output_dir, paste0("Partition_", partition_idx))
    if (!dir.exists(partition_dir)) dir.create(partition_dir)
    
    # Copy the relevant SNP, BOA, and SNP ID files to the partition directory
    for (snp_file in partition_data[[partition_idx]]$snp_files) {
      file.copy(file.path(output_dir, snp_file), file.path(partition_dir, basename(snp_file)))
    }
    for (boa_file in partition_data[[partition_idx]]$boa_files) {
      file.copy(file.path(output_dir, boa_file), file.path(partition_dir, basename(boa_file)))
    }
    for (snp_id_file in partition_data[[partition_idx]]$snp_id_files) {
      file.copy(file.path(output_dir, snp_id_file), 
                file.path(partition_dir, basename(snp_id_file)))
    }
    
    cat("Saved Partition:", partition_idx, "with",
        length(partition_data[[partition_idx]]$snp_files), "windows\n")
  }
}


# ===================================================================
# 6) Final "process" function kept: process_snp_boa.diago_Reduced.V2
# ===================================================================
process_snp_boa.diago_Reduced<- function(i, combined_data, snp_file, boa_file, null_model , eiK) {
  # Set SNP and BOA for the current SNP
  combined_data$SNP <- read_fst(snp_file, columns = as.character(i))
  combined_data$BOA <- read_fst(boa_file, columns = as.character(i))
  
  # Impute missing values in SNP with mean
  combined_data$SNP[is.na(combined_data$SNP)] <- mean(combined_data$SNP, na.rm = TRUE)
  
  # Ensure BOA is numeric (if that is indeed intended)
  combined_data$BOA <- as.numeric(combined_data$BOA)
  
  # Fit SNP-only model
  snp_model <- lmm.diago(
    Y       = combined_data$phenotype, 
    X       = model.matrix(~ as.factor(CG) + SNP, data = combined_data), 
    eigenK  = eiK, 
    verbose = FALSE
  )
  fixed_effects_1 <- snp_model$BLUP_beta
  fixed_var_1     <- snp_model$varbeta
  snp_effect_1    <- fixed_effects_1[length(fixed_effects_1)]
  snp_se_1        <- sqrt(fixed_var_1[length(fixed_effects_1), length(fixed_effects_1)])
  
  # Fit BOA-only model
  boa_model <- lmm.diago(
    Y       = combined_data$phenotype, 
    X       = model.matrix(~ as.factor(CG) + BOA, data = combined_data), 
    eigenK  = eiK, 
    verbose = FALSE
  )
  fixed_effects_2 <- boa_model$BLUP_beta
  fixed_var_2     <- boa_model$varbeta
  boa_effect_2    <- fixed_effects_2[length(fixed_effects_2)]
  boa_se_2        <- sqrt(fixed_var_2[length(fixed_effects_2), length(fixed_effects_2)])
  
  # Wald tests vs. null
  boa_vs_null_pval <- as.numeric(calculate_joint_wald_test(
    boa_model,
    colnames(model.matrix(~ as.factor(CG) + BOA, data = combined_data)),
    "BOA*", "SNP*"
  ))
  snp_vs_null_pval <- as.numeric(calculate_joint_wald_test(
    snp_model,
    colnames(model.matrix(~ as.factor(CG) + SNP, data = combined_data)),
    "SNP*", "BOA*"
  ))
  
  # Fit reduced model with SNP + BOA
  X_reduced    <- model.matrix(~ as.factor(CG) + SNP + BOA, data = combined_data)
  reduced_model <- lmm.diago(
    Y       = combined_data$phenotype, 
    X       = X_reduced, 
    eigenK  = eiK, 
    verbose = FALSE
  )
  fixed_effects_3 <- reduced_model$BLUP_beta
  fixed_var_3     <- reduced_model$varbeta
  
  # Indices for the last two effects (SNP & BOA) in the reduced model
  snp_effect_3 <- fixed_effects_3[length(fixed_effects_3) - 1]
  boa_effect_3 <- fixed_effects_3[length(fixed_effects_3)]
  snp_se_3     <- sqrt(fixed_var_3[length(fixed_effects_3) - 1, length(fixed_effects_3) - 1])
  boa_se_3     <- sqrt(fixed_var_3[length(fixed_effects_3), length(fixed_effects_3)])
  
  reduced_snp_pval <- as.numeric(calculate_joint_wald_test(
    reduced_model, colnames(X_reduced), "SNP*", "BOA*"
  ))
  reduced_boa_pval <- as.numeric(calculate_joint_wald_test(
    reduced_model, colnames(X_reduced), "BOA*", "SNP*"
  ))
  reduced_combined_pval <- as.numeric(calculate_joint_wald_test(
    reduced_model, colnames(X_reduced), c("BOA*", "SNP*")
  ))
  
  # Basic allele frequencies (assuming 0/1/2 coding for SNP)
  allele_freq <- mean(combined_data$SNP) / 2  
  maf         <- min(allele_freq, 1 - allele_freq)
  aaf         <- mean(combined_data$BOA) / 2
  
  # Final table
  out <- data.table(
    SNP_ID                  = i,
    SNP_vs_Null_P_value     = snp_vs_null_pval,
    BOA_vs_Null_P_value     = boa_vs_null_pval,
    Reduced_SNP_P_value     = reduced_snp_pval,
    Reduced_BOA_P_value     = reduced_boa_pval,
    Reduced_Combined_P_value= reduced_combined_pval,
    
    Model1_SNP_effect       = snp_effect_1,
    Model1_SNP_se           = snp_se_1,
    
    Model2_BOA_effect       = boa_effect_2,
    Model2_BOA_se           = boa_se_2,
    
    Model3_SNP_effect       = snp_effect_3,
    Model3_SNP_se           = snp_se_3,
    Model3_BOA_effect       = boa_effect_3,
    Model3_BOA_se           = boa_se_3,
    
    MAF                     = maf,
    Angus_af                = aaf
  )
  
  return(out)
}
