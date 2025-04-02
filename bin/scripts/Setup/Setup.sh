#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=8gb  # Increased memory per CPU
date

# Input arguments
proj_env=$1
NAME=$2
TRAIT=$3
PHENOTYPE_FILE=$4
PHENOTYPE_COLUMN=$5
CG_COLUMN=$6
SNPGeno=$7
BOGeno=$8

source ${proj_env}


# Navigate to the my_results directory for the specified trait
mkdir -p ${my_results}/${NAME}/${TRAIT}/Genotypes
setup_dir=${my_results}/${NAME}/${TRAIT}/Genotypes
cd ${setup_dir}

# Loading modules
ml R
ml plink
newname=${NAME}
# Define new names for intermediate files
newname_SNP=${newname}.SNP
newname_BO=${newname}.BOA
# Copy necessary files to the current directory
cp ${PHENOTYPE_FILE} .

# Get header line from the phenotype file
TRAITS_HEAD=$(echo "#head -1 $(basename ${PHENOTYPE_FILE})")

# Start process logging
echo "Step 1: Extracting IDs from files"

# Extract IDs from the phenotype file (assumes first column contains IDs)
cp ${PHENOTYPE_FILE} pheno.csv
echo "Checkpoint 1: Extracted IDs from phenotype file"

# Extract IDs from SNP and BOA PLINK .fam files
awk '{print $1,$2}' ${SNPGeno}.fam > ID.SNP
awk '{print $1,$2}' ${BOGeno}.fam > ID.BO
echo "Checkpoint 2: Extracted IDs from SNP and BOA files"

# Load R module and run the ID comparison script
Rscript ${my_bin}/scripts/Setup/IDs.Combined.R
echo "Checkpoint 3: Ran R script to find common IDs"

# Check if the intersecting ID file was generated
if [[ -f "intersect.ID" ]]; then
    echo "Checkpoint 4: Common IDs found and saved in intersect.ID"
else
    echo "Error: intersect.ID not generated"
    exit 1
fi

# Remove the temporary ID files to save space
echo "Checkpoint 5: Removed temporary ID files"

# Create binary PLINK files for SNP and BOA data with specific quality control filters
plink -cow -bfile ${SNPGeno} -chr 1-29 -keep intersect.ID -geno 0.05 -mind 0.1 -maf 0.01 -make-bed --keep-allele-order --out SNP.temp
plink -cow -bfile ${BOGeno} -chr 1-29 -keep intersect.ID -geno 0.05 -mind 0.1 -make-bed --keep-allele-order --out BO.temp
echo "Checkpoint 6: PLINK QC and binary files created"

awk '{print $1,$2}' SNP.temp.fam > ID.SNP
awk '{print $1,$2}' BO.temp.fam > ID.BO
Rscript ${my_bin}/scripts/Setup/IDs.Combined.R
# Create PED file using intersecting IDs
awk -NF " " '{print $2,"0","0"}' intersect.ID > ${newname}.PED

# Create reference file for Brahman genotypes
awk -F " " '{print $2, "A"}' BO.temp.bim > ref.Brahman.txt

# Extract common SNPs between SNP and BOA datasets
awk '{print $2}' SNP.temp.bim > SNP.snplist
awk '{print $2}' BO.temp.bim > BO.snplist
comm -12 <(sort SNP.snplist) <(sort BO.snplist) > common.snplist
echo "Checkpoint 7: Extracted common SNPs between SNP and BOA"

# Generate raw and map files for BOA and SNP genotypes using the common SNP list
plink -cow -bfile BO.temp -keep intersect.ID --extract common.snplist --ref-allele ref.Brahman.txt  --recode A -make-bed --out ${NAME}.BOA
plink -cow -bfile SNP.temp -keep intersect.ID --extract common.snplist --recode A -make-bed --out ${NAME}.SNP
echo "Checkpoint 8: Extracted and recoded genotypes for intersecting IDs with common SNPs"

# Generate the map files for BOA and SNP datasets
awk '{print $2,$1,$4}' ${NAME}.BOA.bim > Geno.BOA.map
awk '{print $2,$1,$4}' ${NAME}.SNP.bim > Geno.SNP.map
sed -i -e '1i\SNPID CHR POS' Geno.BOA.map
sed -i -e '1i\SNPID CHR POS' Geno.SNP.map
echo "Checkpoint 9: Map files generated for BOA and SNP"


echo "Checkpoint 8: Performing PCA on SNP data"

# Perform Principal Component Analysis (PCA) on SNP data using PLINK
plink -cow -bfile ${NAME}.SNP --pca 10 --out ${NAME}_PCA

# Check if PCA file was generated
if [[ -f "${NAME}_PCA.eigenvec" ]]; then
    echo "Checkpoint 9: PCA completed and results saved in ${NAME}_PCA.eigenvec"
else
    echo "Error: PCA results not generated"
    exit 1
fi

# Save PCA results for use in R script
echo "Saving PCA results for use in subsequent R analysis"
cp ${NAME}_PCA.eigenvec ${setup_dir}/${NAME}_PCA_results.txt
