library(data.table)
library(parallel)
print(detectCores())
# Read genotype and origin files, only select necessary columns
genotype_data <- fread("genotype.raw", header = TRUE, stringsAsFactors = FALSE, select = c(2,7:ncol(fread("genotype.raw", nrows = 0))))
origin_data <- fread("origin.raw", header = TRUE, stringsAsFactors = FALSE, select =  c(2,7:ncol(fread("origin.raw", nrows = 0))))
# Match origin data to genotype data by IID
origin_data <- origin_data[match(genotype_data$IID, origin_data$IID), ]
# Remove IID column and keep ID separately
genotype_ID <- genotype_data[, 1]
genotype_data <- genotype_data[, -1, with = FALSE]
origin_data <- origin_data[, -1, with = FALSE]
convert_snp <- function(genotype, origin) {
  # Initialize result matrix with default value (5, 5, 5)
  result <- matrix(5, nrow = length(genotype), ncol = 3)

  # Loop through each individual and apply the logic based on genotype and origin
  for (i in seq_along(genotype)) {
    # Skip if either genotype or origin is NA
    if (is.na(genotype[i]) || is.na(origin[i])) {
      next  # Leave the default value (5, 5, 5) for missing data
    }
    
    # Handle genotype == 0
    if (genotype[i] == 0) {
      result[i, ] <- c(0, 0, 0)
    
    # Handle genotype == 1 and various origin values
    } else if (genotype[i] == 1) {
      if (origin[i] == 2) {
        result[i, ] <- c(0, 0, 1)
      } else if (origin[i] == 1) {
        result[i, ] <- c(0, 1, 0)
      } else if (origin[i] == 0) {
        result[i, ] <- c(1, 0, 0)
      }
    
    # Handle genotype == 2 and various origin values
    } else if (genotype[i] == 2) {
      if (origin[i] == 2) {
        result[i, ] <- c(0, 0, 2)
      } else if (origin[i] == 1) {
        result[i, ] <- c(0, 2, 0)
      } else if (origin[i] == 0) {
        result[i, ] <- c(2, 0, 0)
      }
    }
  }
  
  return(result)
}


# Pre-allocate memory for the new data matrix
num_snps <- ncol(genotype_data)
new_data <- vector("list", num_snps)

# Use parallel lapply for efficiency
new_data <- mclapply(1:num_snps, function(i) {
  convert_snp(as.numeric(genotype_data[[i]]), as.numeric(origin_data[[i]]))
}, mc.cores = parallel::detectCores() - 6) # Use all but one core

print("made new_data")
# Combine the genotype data
genotype_matrix <- do.call(cbind, new_data)
print("Finished binding")
# Function to concatenate genotypes for each individual without spaces
concat_genotypes <- apply(genotype_matrix, 1, paste0, collapse = "")
print("Finished concat")
# Combine the IDs with the concatenated genotypes
new_geno <- data.frame(genotype_ID, concat_genotypes)

# Write to file without headers and with no spaces between genotype columns
fwrite(new_geno, "Geno.Xmatrix", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

print("Wrote Geno.Xmatrix")
# Process map file
SNP <- fread("Geno.SNP.map", header = TRUE)
setnames(SNP, c("SNPID", "CHR", "POS"))

# Vectorized approach to create the map file
map <- data.table(
  SNPID = rep(paste0(SNP$SNPID, c("_BB", "_AB", "_AA")), each = 3),
  CHR = rep(SNP$CHR, each = 3),
  POS = rep(SNP$POS, each = 3)
)

# Write map file efficiently
fwrite(map, "Geno.Xmatrix.map", quote = FALSE, row.names = FALSE, col.names = TRUE,sep= " ")