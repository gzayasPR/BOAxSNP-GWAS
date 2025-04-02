#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1gb
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
memory_limit=6Gb     # Memory limit for each partition job
cpu_limit=1           # CPU usage for each partition
max_time="72:00:00"   # Maximum time for each partition job

# Load environment variables
source ${proj_env}
source ${my_bin}/scripts/shell.functions.sh

# Output directories
out_dir=${my_results}/${NAME}
gwas_dir=${out_dir}/${TRAIT}/Marker_GWAS

# Define files and directories
setup_dir=${out_dir}/${TRAIT}/Genotypes
path_SNP="${gwas_dir}/snp_data.fst"
path_BOA="${gwas_dir}/boa_data.fst"
partition_file="${gwas_dir}/${NAME}_partition_indexes.csv"
# Create output directory if it doesn't exist
mkdir -p $gwas_dir/bash
cd $gwas_dir/bash
rm $gwas_dir/bash/*

    # Step 2: Submit GWAS jobs for each partition
    for i in $(seq 1 $partitions); do
        job_id=$(sbatch --account="mateescu" \
            --job-name="${i}_Marker_GWAS_${TRAIT}" \
            --cpus-per-task=${cpu_limit} \
            --mem-per-cpu=${memory_limit} \
            --time=${max_time} \
            --output="${i}_Marker_GWAS_${TRAIT}.out" \
            --error="${i}_Marker_GWAS_${TRAIT}.err" \
            $my_bin/scripts/GWAS/Marker_Based/runGWAS.partitions.sh \
            $proj_env $NAME $TRAIT $i | awk '{print $4}')
            job_ids+=($job_id)
    done

# echo $job_ids
result_files="${gwas_dir}/Results_Partition/"
mkdir -p $result_files
# Get the start time of the GWAS process
start_time_gwas=$(date +%s)
echo "Waiting for all GWAS jobs to complete..."
time_min=5
job_file=${my_bash}/${NAME}/${TRAIT}/GWAS/jobs.txt
echo "${job_ids[@]}" > ${job_file}
wait_job_check_progress_GWAS "${job_file}" "${time_min}" "${start_time_gwas}" "${result_files}"
cd $gwas_dir/bash
rm *.txt