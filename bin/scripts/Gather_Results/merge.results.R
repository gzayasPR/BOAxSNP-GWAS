#!/usr/bin/env Rscript

# Load necessary libraries
library(data.table)
# Input arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
results_dir <- args[1]         # Directory containing results for all partitions
map_file <- args[2]            # Map file for SNP position information
output_file <- args[3]    
my_functions_R <- args[4] 
source(my_functions_R)
# Step 1: Merge results from all partitions
cat("Merging results from partitions...\n")
results_files <- list.files(results_dir, pattern = "partition.*\\.csv", full.names = TRUE)

# Initialize an empty data table to store results
all_results <- data.table()

# Read and merge all partition results
for (file in results_files) {
  partition_data <- fread(file)
  all_results <- rbind(all_results, partition_data, fill = TRUE)
}
head(  all_results)
cat("Merged all partition results.\n")

cat("Merging with SNP map file...\n")
map_data <- fread(map_file)
names(map_data) <- c("SNP_ID","Chromosome","Position")
map_data$Index <- 1:nrow(map_data) 
names(all_results)[1] <- "Index"
# Merge the map data (assuming map file has columns: Real_SNP, Chromosome, Position)
merged_data <- merge( all_results, map_data, by = "Index", all.x = TRUE)

cat("Merged with SNP map file.\n")

# Step 5: Save the merged results
cat("Saving merged results...\n")
fwrite(merged_data, file = output_file)
cat("Merged results saved to:", output_file, "\n")

# The data is now ready for further analysis, including plotting the results.
