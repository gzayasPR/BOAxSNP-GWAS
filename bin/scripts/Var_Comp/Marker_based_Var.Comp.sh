#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50gb
#SBATCH --time=72:00:00

# Marker-based GWAS analysis
# Author: Gabe Zayas, Date: 10/17/2024
# This script runs marker-based GWAS with partitioned data (SNP and BOA) and batch processing.
# Output is saved into separate folders by trait.

# Input arguments
proj_env=$1
NAME=$2
TRAIT=$3
PHENOTYPE_FILE=$4
PHENOTYPE_COLUMN=$5
CG_COLUMN=$6
partitions=$7    # Number of partitions to divide the data into
window_size=$8         # Number of batches within each partition


# Additional variables (adjustable based on your requirements)
id_col=1              # Column for individual IDs in the phenotype file
# Load environment variables
source ${proj_env}
source ${my_bin}/scripts/shell.functions.sh
my_functions_R=${my_bin}/scripts/my_functions.R
# Output directories
out_dir=${my_results}/${NAME}
gwas_dir=${out_dir}/${TRAIT}/Marker_GWAS
partition_dir="${gwas_dir}/Partitions/"

# Define files and directories
setup_dir=${out_dir}/${TRAIT}/Genotypes
geno_BOA="${setup_dir}/${NAME}.BOA.raw"
geno_SNP="${setup_dir}/${NAME}.SNP.raw"
pheno_file="${setup_dir}/pheno.csv"
ped_file="${setup_dir}/${NAME}.PED"
map_file="${setup_dir}/Geno.SNP.map"
pca_results="${setup_dir}/${NAME}_PCA_results.txt"
# Create output directory if it doesn't exist
mkdir -p $gwas_dir
cd $gwas_dir
ml R

# Step 1: Run the null model and partition the data
echo "Running null model and partitioning data into ${partitions} partitions with batch size ${window_size}..."

Rscript $my_bin/scripts/Var_Comp/Marker_Based/run_nullmodel_marker.R \
    ${NAME} ${pheno_file} ${PHENOTYPE_COLUMN} ${CG_COLUMN} ${setup_dir}/${NAME}.SNP.raw \
    ${setup_dir}/${NAME}.BOA.raw $partitions $window_size ${pca_results} ${map_file} ${my_functions_R}

path_SNP="${gwas_dir}/snp_data.fst"
path_BOA="${gwas_dir}/boa_data.fst"
Rscript $my_bin/scripts/Var_Comp/Marker_Based/checkparitions.R \
    ${path_SNP} ${path_BOA} $gwas_dir/${NAME}_null_model.Rdata ${my_functions_R}