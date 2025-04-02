# Load necessary library
suppressPackageStartupMessages({
    library(data.table)})

# Read phenotype and genotype files
pheno <- fread("pheno.csv", header = TRUE)
id_snp <- fread("ID.SNP", header = FALSE)
id_bo <- fread("ID.BO", header = FALSE)

# Extract IID from the phenotype file (assuming it's in the first column of pheno)
iid_pheno <- pheno[[1]]

# Combine FID and IID from the SNP and BO genotype files
fid_iid_snp <- paste(id_snp[[1]], id_snp[[2]], sep = "_")
fid_iid_bo <- paste(id_bo[[1]], id_bo[[2]], sep = "_")

# Extract common IIDs between the phenotype and both genotype files
common_iids <- Reduce(intersect, list(iid_pheno, id_snp[[2]], id_bo[[2]]))

# Extract the full FID_IID for the common IIDs
intersect_ids_snp <- id_snp[id_snp[[2]] %in% common_iids, ]
intersect_ids_bo <- id_bo[id_bo[[2]] %in% common_iids, ]

# Write intersected FID_IID from SNP and BO files to intersect.ID
intersect_ids <- intersect(intersect_ids_snp, intersect_ids_bo)
fwrite(intersect_ids, "intersect.ID", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

# Filter phenotype data to only include rows with common IIDs
pheno_filtered <- pheno[pheno[[1]] %in% common_iids, ]

# Write the filtered phenotype data back to pheno.csv
fwrite(pheno_filtered, "pheno.csv", quote = FALSE, row.names = FALSE)

print("Common IDs written to intersect.ID and pheno.csv updated with common IIDs.")
