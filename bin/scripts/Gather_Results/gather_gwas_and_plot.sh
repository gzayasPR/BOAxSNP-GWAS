#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20Gb
#SBATCH --time=72:00:00

## Single-step GWAS analysis using BLUPF90 ##
# First round: unweighted ssGWAS
# Second and third rounds: weighted ssGWAS (updating G matrix)
# Output is saved into separate folders by trait


# Input arguments
proj_env=$1
NAME=$2
TRAIT=$3 

# Load environment variables
source ${proj_env}
source ${my_bin}/scripts/shell.functions.sh
my_functions_R=${my_bin}/scripts/my_functions.R
# Output directories
out_dir=${my_results}/${NAME}
gwas_dir=${out_dir}/${TRAIT}/Marker_GWAS
results_dir="${gwas_dir}/Results_Partition/"
partition_dir="${gwas_dir}/Partitions/"

# Define files and directories
setup_dir=${out_dir}/${TRAIT}/Genotypes
map_file="${setup_dir}/Geno.SNP.map"

# Get current date in mm.dd.yyyy format
current_date=$(date +'%m.%d.%Y')
# Run R script to merge results after all jobs are done
merged_results="${gwas_dir}/${TRAIT}.${current_date}.marker.GWAS.results.csv"
ml R

# Step 4: Merge results from all partitions
echo "Merging GWAS results..."
Rscript $my_bin/scripts/Gather_Results/merge.results.R \
    ${results_dir} ${map_file} ${merged_results}  ${my_functions_R}

# Step 5: Generate Manhattan plot
echo "Generating Manhattan plot..."
Rscript $my_bin/scripts/Gather_Results/plot.GWAS.marker.R \
    ${merged_results} ${gwas_dir}/ ${TRAIT}  ${my_functions_R}

Rscript $my_bin/scripts/Gather_Results/Process.GWAS.R \
    ${merged_results} ${gwas_dir}/ ${TRAIT} 

# echo "Marker-based GWAS analysis completed."
