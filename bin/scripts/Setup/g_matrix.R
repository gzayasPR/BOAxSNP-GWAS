library(AGHmatrix)
library(data.table)

# Step 1: Load the SNP data
# Read the SNP raw file generated from PLINK
snp_data <- fread(snp_file)

# Drop non-genotype columns (assuming first 6 columns are non-genetic)
genotype_matrix <- snp_data[, -(1:6), with = FALSE]

# Convert the genotype matrix to numeric (if necessary)
genotype_matrix <- as.matrix(genotype_matrix)

# Step 2: Create the G matrix using AGHmatrix
G_matrix <- Gmatrix(SNPmatrix = genotype_matrix, method = "VanRaden", ploidy = 2)