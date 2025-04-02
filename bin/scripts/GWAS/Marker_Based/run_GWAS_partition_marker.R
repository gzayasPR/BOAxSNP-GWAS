#!/usr/bin/env Rscript

# Input arguments
args <- commandArgs(trailingOnly = TRUE)
index_file <- args[1]  # File containing SNP and BOA batch file paths for the group
parition_num <- args[2]    # File containing SNP and BOA batch file paths for the group
G_matrix_path <- args[3]     # Null model RData file
output_file <- args[4]         # Output file for results
my_functions_R <- args[5] 
snp_file <- args[6] 
boa_file <- args[7]  

source(my_functions_R)

# Delete output file if it already exists
if (file.exists(output_file)) {
  file.remove(output_file)
  cat("Previous output file deleted: ", output_file, "\n")
}

# Measure time for performance monitoring
start_time <- proc.time()

index <- read.csv(index_file)
indexes <- index[index$Global_Partition == parition_num,]
chr <- unique(indexes$Chromosome)
# Load the G matrix and null model with combined data
cat("Loading G matrix, null model, and combined data...\n")
print(paste0(G_matrix_path,"/Gmatrix_chr",chr,".Rdata"))
load(paste0(G_matrix_path,"/Gmatrix_chr",chr,".Rdata"))
# Initialize results table
results <- c()

# Loop through the group of windows
for (i in seq(indexes$Start_Index,indexes$End_Index)) {
 result <- process_snp_boa.diago_Reduced(i, combined_data, snp_file,  boa_file , null_model ,eiK)
    # Ensure column names of result match those of results before appending
    # setnames(result, names(results))
    results <- rbind(results, result, fill = TRUE) }

# Save the results to a CSV file
cat("Saving results to", output_file, "\n")
fwrite(results, output_file)

# Print the time taken for processing
end_time <- proc.time() - start_time
cat("Processed ",nrow(results)," SNPs in " , end_time[3], "seconds\n")
