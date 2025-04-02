#!/usr/bin/env Rscript

# Input arguments from shell script
args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
phenotype_file <- args[2]
phenotype_column <- as.numeric(args[3])
cg_column <- as.numeric(args[4])
snp_file <- args[5]
boa_file <- args[6]
partitions <- as.numeric(args[7])
Window_size <- as.numeric(args[8])
pca_file <- args[9] # New argument for PCA results
map_file <- args[10] # New argument for PCA results
my_functions_R <- args[11] 


source(my_functions_R)
# Print the input arguments
cat("Name:", name, "\n")
cat("Phenotype file:", phenotype_file, "\n")
cat("Phenotype column:", phenotype_column, "\n")
cat("CG column:", cg_column, "\n")
cat("SNP file:", snp_file, "\n")
cat("BOA file:", boa_file, "\n")
cat("Partition number:", partitions , "\n")
cat("PCA file:", pca_file, "\n")
cat("map file:", map_file, "\n")



# Load SNP and BOA data
cat("Loading and sorting SNP and BOA data...\n")
snp_data <- fread(snp_file)
boa_data <- fread(boa_file)

# Load phenotype and PCA data
cat("Loading and sorting phenotype data...\n")
phenotypes <- fread(phenotype_file, header = TRUE)
pca_data <- fread(pca_file, header = FALSE)
colnames(pca_data) <- c("FID", "ID", paste0("PC", 1:10))

# Ensure IDs are character type to avoid merge issues
phenotypes$ID <- as.character(phenotypes[[1]])
phenotypes$CG <- as.character(phenotypes[[cg_column]])
phenotypes$phenotype <-  phenotypes[[phenotype_column ]]
pca_data$ID <- as.character(pca_data$ID)

# Merge PCA data with phenotype data
phenotypes <- merge(phenotypes, pca_data, by = "ID")
phenotypes <- phenotypes[order(phenotypes$ID)]

# Filter common IDs across SNP, BOA, phenotype, and G-matrix data
cat("Filtering common IDs across SNP, BOA, phenotype, and G-matrix data...\n")
common_ids <- Reduce(intersect, list(phenotypes$ID, snp_data[[2]], boa_data[[2]]))

# Subset each dataset to keep only rows with common IDs and order by ID
phenotypes <- phenotypes[phenotypes$ID %in% common_ids][order(phenotypes$ID)]
snp_data <- snp_data[snp_data[[2]] %in% common_ids][order(snp_data[[2]])]
boa_data <- boa_data[boa_data[[2]] %in% common_ids][order(boa_data[[2]])]

# Step 4: Create the G-matrix using AGHmatrix (VanRaden method) and ensure alignment
cat("Creating the G-matrix...\n")
genotype_matrix <- as.matrix(snp_data[, -(1:6), with = FALSE])  # Remove non-genotype columns (first 6 columns)
G_matrix <- Gmatrix(SNPmatrix = genotype_matrix, method = "VanRaden", ploidy = 2)
rownames(G_matrix) <- snp_data[[2]]  # Set row names to individual IDs
colnames(G_matrix) <- snp_data[[2]]

# Step 5: Check if the G matrix is singular using determinant
det_G <- det(G_matrix)
cat("Determinant of G matrix: ", det_G, "\n")


# Threshold for determining singularity (close to zero)
singularity_threshold <- 1e-7

# Calculate initial determinant
det_G <- det(G_matrix)
cat("Initial determinant of G matrix: ", det_G, "\n")

if (abs(det_G) <= singularity_threshold) {
  cat("G matrix is singular or near-singular. Applying manual blending...\n")
  
  lambda <- 0.01  # Start with a small value
  det_G_new <- det_G
  while (abs(det_G_new) <= singularity_threshold && lambda <= 1) {
    G_matrix <- G_matrix + lambda * diag(nrow(G_matrix))
    det_G_new <- det(G_matrix)
    lambda <- lambda + 0.01
  }
  
  if (abs(det_G_new) <= singularity_threshold) {
    stop("Blended G matrix is still singular even after increasing lambda. Please check your data.")
  } else {
    cat("Blended G matrix is now non-singular with lambda =", lambda - 0.01, "\n")
  }
} else {
  cat("G matrix is non-singular. No blending needed.\n")
}


# Step 5: Create the B-matrix using AGHmatrix (VanRaden method) and ensure alignment
cat("Creating the B-matrix...\n")
boa_matrix <- as.matrix(boa_data[, -(1:6), with = FALSE])  # Remove non-genotype columns (first 6 columns)
B_matrix<- Gmatrix(SNPmatrix = boa_matrix, method = "VanRaden", ploidy = 2)
rownames(B_matrix) <- boa_data[[2]]  # Set row names to individual IDs
colnames(B_matrix) <- boa_data[[2]]

# Step 5: Check if the det_B is singular using determinant
det_B <- det(B_matrix)
cat("Determinant of det_B: ", det_B, "\n")


# Threshold for determining singularity (close to zero)
singularity_threshold <- 1e-7

# Calculate initial determinant
det_B <- det(B_matrix)
cat("Initial determinant of det_B: ", det_B, "\n")

if (abs(det_B) <= singularity_threshold) {
  cat("det_B is singular or near-singular. Applying manual blending...\n")
  
  lambda <- 0.01  # Start with a small value
  det_B_new <- det_B
  while (abs(det_B_new) <= singularity_threshold && lambda <= 1) {
    B_matrix<- B_matrix+ lambda * diag(nrow(B_matrix))
    det_B_new <- det(B_matrix)
    lambda <- lambda + 0.01
  }
  
  if (abs(det_B_new) <= singularity_threshold) {
    stop("Blended det_B is still singular even after increasing lambda. Please check your data.")
  } else {
    cat("Blended det_B is now non-singular with lambda =", lambda - 0.01, "\n")
  }
} else {
  cat("det_B is non-singular. No blending needed.\n")
}

# Save the B-matrix for later use
cat("Saving the B-matrix...\n")

# Ensure G_matrix, phenotypes, snp_data, and boa_data have the same order of individuals
phenotypes <- phenotypes[order(phenotypes$ID)]
snp_data <- snp_data[match(phenotypes$ID, snp_data[[2]])]
boa_data <- boa_data[match(phenotypes$ID, boa_data[[2]])]
G_matrix <- G_matrix[phenotypes$ID, phenotypes$ID]  # Reorder G_matrix by ID
B_matrix <- B_matrix[phenotypes$ID, phenotypes$ID]  # Reorder B_matrix by ID

# Add SNP and BOA columns to combined data
combined_data <- data.table(ID = phenotypes$ID, 
                            phenotype = phenotypes$phenotype, 
                            CG = as.factor(phenotypes$CG))
cat("Fitting the null model...\n")

eiK <- eigen(G_matrix)
eiB <- eigen(B_matrix)
start_time_null <- proc.time()
# Fit the null model using gaston's lmm.aireml function
cat("Fitting the null model using gaston...\n")
null_model <- lmm.diago(Y = combined_data$phenotype,
                                X = model.matrix(~ as.factor(CG) + 0, data = combined_data),
                                eigenK = eiK ,
                                verbose = TRUE)

end_time_null <- proc.time() - start_time_null
cat("Null model fitting completed in", end_time_null[3], "seconds\n")
cat("Extracting variance components from gaston model...\n")
var_comp <- data.table(
  Component = c("Genetic", "Residual"),
  Variance = c(null_model$tau, null_model$sigma2),
  Proportion = c(null_model$tau /(null_model$tau + null_model$sigma2 ), null_model$sigma2 /(null_model$tau + null_model$sigma2 ))
)
cat("Variance Components (gaston):",  "\n")
print(var_comp)

# Save the null model
combined_data$BLUPs <- null_model$BLUP_omega
write.csv(combined_data,file = paste0(name, "_BLUPs.csv"))

cat("Fitting the null model using gaston...\n")
null_model2 <- lmm.aireml(Y = combined_data$phenotype,
                                X = model.matrix(~ as.factor(CG) + 0, data = combined_data),
                                K = list(G_matrix,B_matrix) ,get.P = TRUE,
                                verbose = TRUE)

cat("Null model G + B fitting completed in", end_time_null[3], "seconds\n")

cat("Extracting variance components from gaston model...\n")
var_comp2 <- data.table(
  Component = c("Genetic","BOA", "Residual"),
  Variance = c(null_model2$tau[1],null_model2$tau[2], null_model2$sigma2),
  Proportion = c(null_model2$tau[1] /(null_model2$tau[1] + null_model2$tau[2]  + null_model2$sigma2 ),
                 null_model2$tau[2] /(null_model2$tau[1] + null_model2$tau[2]  + null_model2$sigma2) , 
                 null_model2$sigma2 /(null_model2$tau[1] + null_model2$tau[2]  + null_model2$sigma2 )
  ))
cat("Variance Components (gaston):",  "\n")
print(var_comp2)



# Output variance components to a file
cat("Saving variance components...\n")
fwrite(var_comp, file = paste0(name, "_variance_components.csv"))

# Step 8: Calculate PCs based on genotypes
cat("Fitting the first SNP and BOA interaction...\n")
theta_starting <- c(null_model$sigma2,null_model$tau)
save(null_model,combined_data,G_matrix,theta_starting,eiK , file = paste0(name, "_null_model.Rdata"))
combined_data$SNP <- snp_data[[7]]
combined_data$BOA <- boa_data[[7]] # Ensure BOA is a factor
# Fit the full model with the first SNP * BOA interaction
start_time_update <- proc.time()
# Impute missing values in SNPs
combined_data$SNP[is.na(combined_data$SNP)] <- mean(combined_data$SNP, na.rm = TRUE)
# Fit the full model with SNPs using gaston
cat("Fitting the full model with SNP and BOA interaction using gaston...\n")
X_full <- model.matrix(~ as.factor(CG) + 0 + SNP * as.factor(BOA), data = combined_data)
full_model <- lmm.diago(Y = combined_data$phenotype,
                                X = X_full,
                                eigenK = eiK ,
                                verbose = TRUE)
end_time_update <- proc.time() - start_time_update
cat("New model fitting completed in", end_time_update[3], "seconds\n")
cat("New model summary:\n")
# Step 11: Perform ANOVA to compare null_model and full_model
# Extract the names of the columns from the design matrix X_full
parameter_names <- colnames(X_full)
snp_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "SNP*","BOA*"))
boa_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "BOA*","SNP*"))
interaction_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names, "SNP:as.factor(BOA)*"))
combined_pval <- as.numeric(calculate_joint_wald_test(full_model, parameter_names , c("BOA*","SNP*")))
# Perform ANOVA comparisons with the null model
results <- data.table(SNP_ID = colnames(snp_data)[7], SNP_P_value = snp_pval, BOA_P_value = boa_pval,
                    Interaction_P_value = interaction_pval, Combined_P_value = combined_pval)
print(results)

# ################################# Partition Data #################################
# # Step 8: Partition SNP and BOA data into smaller batches
# Calculate partition indexes and save to a CSV file
snp_data2 <- snp_data[, -(1:6), with = FALSE] 
boa_data2 <- boa_data[, -(1:6), with = FALSE] 
#Rename columns to be numeric based on their position
setnames(snp_data2, as.character(seq_len(ncol(snp_data2))))
setnames(boa_data2, as.character(seq_len(ncol(snp_data2))))

# Save the data in Parquet format
write_fst(snp_data2, "snp_data.fst")
write_fst(boa_data2, "boa_data.fst")
##############################
# ---- NEW SECTION: Create leave-one-chromosome out G_matrices ----
cat("Creating leave-one-chromosome out G_matrices...\n")

# Load the map file to obtain chromosome information for each SNP
map_data <- fread(map_file)
# Assume map_data has columns: CHR, SNP, CM, BP

# Extract the genotype (SNP) data from snp_data (columns 7 onward)
genotype_data <- snp_data[, -(1:6), with = FALSE]
names(genotype_data) <- sub("_[A-Za-z]+$", "", names(genotype_data))
# Reorder map_data to match the order of SNP columns in genotype_data based on SNP IDs
map_data_ordered <- map_data[match(names(genotype_data), map_data$SNP)]
if(any(is.na(map_data_ordered$CHR))) {
  stop("Some SNP names in the genotype data do not match the SNP identifiers in the map file.")
}
# Get a vector indicating the chromosome for each SNP (in the same order as genotype_data)
snp_chr <- map_data_ordered$CHR

# Create a directory for leave-one-out G matrices if it doesn't exist
if (!dir.exists("G_matrices")) {
  dir.create("G_matrices")
}

# For each unique chromosome, compute a G_matrix using all SNPs except those on that chromosome
unique_chrs <- sort(unique(snp_chr))
for(chr in unique_chrs) {
  cat("Processing leave-one-out G_matrix for chromosome:", chr, "\n")
  
  # Identify the columns (SNPs) NOT on the current chromosome
  cols_to_include <- which(snp_chr != chr)
  
  # Subset the genotype matrix accordingly
  genotype_subset <- as.matrix(genotype_data[, cols_to_include, with = FALSE])
  
  # Compute the G_matrix using the VanRaden method (from your AGHmatrix function)
  G_matrix_leave <- Gmatrix(SNPmatrix = genotype_subset, method = "VanRaden", ploidy = 2)
  rownames(G_matrix_leave) <- snp_data[[2]]  # set individual IDs as row names
  colnames(G_matrix_leave) <- snp_data[[2]]
  
  # Optionally check for singularity and blend if necessary (using the same threshold as before)
  det_G_leave <- det(G_matrix_leave)
  if (abs(det_G_leave) <= singularity_threshold) {
    cat("G_matrix for chromosome", chr, "is singular or near-singular. Applying blending...\n")
    lambda <- 0.01
    det_G_new <- det_G_leave
    while (abs(det_G_new) <= singularity_threshold && lambda <= 1) {
      G_matrix_leave <- G_matrix_leave + lambda * diag(nrow(G_matrix_leave))
      det_G_new <- det(G_matrix_leave)
      lambda <- lambda + 0.01
    }
    if (abs(det_G_new) <= singularity_threshold) {
      stop(paste("Blended G matrix for chromosome", chr, "is still singular. Please check your data."))
    } else {
      cat("Blended G matrix for chromosome", chr, "is now non-singular with lambda =", lambda - 0.01, "\n")
    }
  }
  eiK <- eigen(G_matrix_leave)
  # Save the leave-one-out G_matrix to the G_matrices directory
  save_file <- paste0("G_matrices/Gmatrix_chr", chr, ".Rdata")
  save(G_matrix_leave,null_model,combined_data,theta_starting,eiK, file = save_file)
  cat("Saved leave-one-out G_matrix for chromosome", chr, "to", save_file, "\n")
}

##############################
# ---- NEW SECTION: Chromosome-specific SNP Partitioning ----
cat("Creating chromosome-specific SNP partitions...\n")

# Total number of SNPs in the genotype data
total_snps <- ncol(genotype_data)

# The overall (global) partition size if SNPs were evenly distributed
global_partition_size <- total_snps / partitions

# Initialize a data frame to store partition indexes with an extra column for chromosome
partition_indexes_chr <- data.frame(Global_Partition = integer(),
                                    Chromosome = character(),
                                    Chr_Partition = integer(),
                                    Start_Index = integer(),
                                    End_Index = integer(),
                                    stringsAsFactors = FALSE)

global_counter <- 0

# Loop over each chromosome to partition SNPs within that chromosome
for(chr in unique_chrs) {
  # Get indices of SNPs that belong to the current chromosome
  chr_indices <- which(snp_chr == chr)
  n_chr <- length(chr_indices)
  
  # Determine the number of partitions for this chromosome proportionally
  # Here we use the proportion of SNPs on the chromosome relative to the global total.
  partitions_for_chr <- round(n_chr / total_snps * partitions)
  if (partitions_for_chr < 1) partitions_for_chr <- 1
  
  # Calculate the window (partition) size for this chromosome (number of SNPs per partition)
  partition_size_chr <- ceiling(n_chr / partitions_for_chr)
  
  cat("Chromosome", chr, "has", n_chr, "SNPs and will be partitioned into", partitions_for_chr, "windows.\n")
  
  # Create partitions for the current chromosome
  for (i in 1:partitions_for_chr) {
    global_counter <- global_counter + 1
    start_rel <- (i - 1) * partition_size_chr + 1
    end_rel <- min(i * partition_size_chr, n_chr)
    # Map the relative (within chromosome) indices to the global genotype_data indices
    start_index <- chr_indices[start_rel]
    end_index <- chr_indices[end_rel]
    
    partition_indexes_chr <- rbind(partition_indexes_chr,
                                   data.frame(Global_Partition = global_counter,
                                              Chromosome = chr,
                                              Chr_Partition = i,
                                              Start_Index = start_index,
                                              End_Index = end_index,
                                              stringsAsFactors = FALSE))
  }
}

# Save the chromosome-specific partition index file (with an added Chromosome column)
index_file_chr <- paste0(name, "_partition_indexes.csv")
write.csv(partition_indexes_chr, index_file_chr, row.names = FALSE)
cat("Chromosome-specific partition indexes saved to:", index_file_chr, "\n")
