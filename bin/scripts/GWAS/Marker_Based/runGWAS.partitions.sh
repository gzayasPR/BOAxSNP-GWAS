#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  # Set this depending on how many cores you want to allocate per partition

# Input arguments
proj_env=$1
NAME=$2
TRAIT=$3
PARTITION_NUM=$4  # This is the partition number to process
# Load environment variables
source ${proj_env}
source ${my_bin}/scripts/shell.functions.sh
my_functions_R=${my_bin}/scripts/my_functions.R
# Define directories and files
out_dir=${my_results}/${NAME}
gwas_dir=${out_dir}/${TRAIT}/Marker_GWAS
path_SNP="${gwas_dir}/snp_data.fst"
path_BOA="${gwas_dir}/boa_data.fst"
G_matrix_path="${gwas_dir}/G_matrices/"
output_dir="${gwas_dir}/Results_Partition/"
index_file="${gwas_dir}/${NAME}_partition_indexes.csv"
# Create output directory if it doesn't exist
mkdir -p ${output_dir}

# Load R module
ml R

Rscript ${my_bin}/scripts/GWAS/Marker_Based/run_GWAS_partition_marker.R "$index_file" "${PARTITION_NUM}" \
        "${G_matrix_path}" "${output_dir}/partition.${PARTITION_NUM}.csv" "${my_functions_R}" "${path_SNP}" "${path_BOA}"
    