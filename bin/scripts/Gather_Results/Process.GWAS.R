#!/usr/bin/env Rscript
# Read input arguments
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
plot_dir <- args[2]
trait <- args[3]
bonferroni_threshold <- 0.1

# Load required library
library(dplyr)

process_gwas_hits <- function(input_csv, bonferroni_threshold, output_name) {
  
  # Read the input file
  data <- read.csv(input_csv)
  
  # Define the p-value columns to process
  pvalue_columns <- c("Reduced_SNP_P_value", "Reduced_BOA_P_value", 
                      "Reduced_Combined_P_value", "SNP_vs_Null_P_value", 
                      "BOA_vs_Null_P_value")
  
  # Create a mapping of models and effects
  model_mapping <- list(
    Reduced_SNP_P_value = list(Model = 3, Effect = "SNP"),
    Reduced_BOA_P_value = list(Model = 3, Effect = "BOA"),
    Reduced_Combined_P_value = list(Model = 3, Effect = "Joint"),
    SNP_vs_Null_P_value = list(Model = 1, Effect = "SNP"),
    BOA_vs_Null_P_value = list(Model = 2, Effect = "BOA")
  )
  
  # Initialize result dataframe
  results <- data.frame()
  
  # Process each p-value column independently
  for (col in pvalue_columns) {
    
    # Get unique p-values for the column
    unique_pvalues <- unique(data[[col]])

    # Identify Bonferroni-corrected threshold
    bonferroni_corrected_threshold <- bonferroni_threshold / length(unique_pvalues)
    bonferroni_suggestive_threshold <- 1/ length(unique_pvalues)
    # Annotate the data for this column
    temp_data <- data %>%
      mutate(
        pvalue = .data[[col]],
        Significance_Level = case_when(
          pvalue <= bonferroni_corrected_threshold ~ "Significant",
          pvalue <= bonferroni_suggestive_threshold ~ "Suggestive",
          TRUE ~ NA_character_
        )
      ) %>% 
      filter(!is.na(Significance_Level)) %>%
      mutate(
        SNP = SNP_ID,
        CHR = Chromosome,
        POS = Position,
        Model = model_mapping[[col]]$Model,
        Effect = model_mapping[[col]]$Effect
      ) %>%
      select(SNP, CHR, POS, pvalue, Model, Effect, Significance_Level)
    
    # Append to the results if there are any hits
    if (nrow(temp_data) > 0) {
      results <- bind_rows(results, temp_data)
    } else {
      message(paste("No significant or suggestive hits found for column:", col))
    }
  }
  
  # Write results to output directory
  write.csv(results, output_name, row.names = FALSE)
  message("Processing complete. Results saved to ", output_name)
}

output_file <- file.path(plot_dir, paste0(trait, "_Reduced_GWAS_hits.csv"))
process_gwas_hits(input_csv = input_file, bonferroni_threshold = bonferroni_threshold, output_name = output_file)