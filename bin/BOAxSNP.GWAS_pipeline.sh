#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Run_Thermo_ALL
#SBATCH --output=Run_Thermo_ALL_%j.out
#SBATCH --error=Run_Thermo_ALL_%j.err
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=1gb

# Load project environment and functions
echo "Loading project environment and shell functions..."
source project_env.sh
source "${my_bin}/scripts/shell.functions.sh"

# Source the input file containing all parameters
source pipeline_input.sh

# Ensure directories exist
echo "Creating output directories..."
bash_out="${my_bash}"
mkdir -p "${bash_out}"
mkdir -p "${bash_out}/${NAME}"
out_dir="${my_results}/${NAME}"
job_ids=()

# Check if all step flags are set
if [[ -z $runSetup || -z $runVARCOMP || -z $runGWAS || -z $runGatherPlots ]]; then
  echo "Error: Missing arguments. Please provide TRUE/FALSE for each step (runSetup, runVARCOMP, runGWAS, runGatherPlots)."
  exit 1
fi

# Step 1: Setup
if [[ $runSetup == "TRUE" ]]; then
  echo "Starting BLUPF90 Setup for each trait..."
  for i in "${!TRAIT_NAMES[@]}"; do
    TRAIT="${TRAIT_NAMES[$i]}"
    echo "Processing trait: $TRAIT"
    mkdir -p "${bash_out}/${NAME}/${TRAIT}/Setup/"
    cd "${bash_out}/${NAME}/${TRAIT}/Setup/"

    PHENOTYPE_FILE="${phenotypes_files[$i]}"
    PHENOTYPE_COLUMN="${pheno_pos[$i]}"
    CG_COLUMN="${CG_pos[$i]}"

    echo "Submitting BLUPF90 setup job for ${TRAIT}..."
    job_id=$(sbatch --account="mateescu" \
             --job-name="Setup.${TRAIT}" \
             --output="Setup.${TRAIT}.out" \
             --error="Setup.${TRAIT}.err" \
             "${my_bin}/scripts/Setup/Setup.sh" "${proj_env}" "${NAME}" "${TRAIT}" "${PHENOTYPE_FILE}" "${PHENOTYPE_COLUMN}" "${CG_COLUMN}" "${SNPGeno}" "${BOGeno}" | awk '{print $4}')
    job_ids+=("$job_id")
    sleep 5
  done

  echo "Waiting for all Setup jobs to complete..."
  for job_id in "${job_ids[@]}"; do
    wait_for_job_completion "$job_id" 1
  done
fi

# Step 2: Var_Comp analysis
if [[ $runVARCOMP == "TRUE" ]]; then
  echo "Starting Var_Comp analysis for each trait..."
  for i in "${!TRAIT_NAMES[@]}"; do
    TRAIT="${TRAIT_NAMES[$i]}"
    echo "Processing trait: $TRAIT"
    mkdir -p "${bash_out}/${NAME}/${TRAIT}/Var_Comp/"
    cd "${bash_out}/${NAME}/${TRAIT}/Var_Comp/"

    PHENOTYPE_FILE="${phenotypes_files[$i]}"
    PHENOTYPE_COLUMN="${pheno_pos[$i]}"
    CG_COLUMN="${CG_pos[$i]}"

    echo "Submitting Var_Comp job for ${TRAIT}..."
    job_id=$(sbatch --account="mateescu" \
             --job-name="Var_Comp_${TRAIT}" \
             --output="Var_Comp_${TRAIT}.out" \
             --error="Var_Comp_${TRAIT}.err" \
             "${my_bin}/scripts/Var_Comp/Marker_based_Var.Comp.sh" "${proj_env}" "${NAME}" "${TRAIT}" "${PHENOTYPE_FILE}" "${PHENOTYPE_COLUMN}" "${CG_COLUMN}" "${partitions}" "${window_size}" | awk '{print $4}')
    job_ids+=("$job_id")
    sleep 5
  done

  echo "Waiting for all Var_Comp jobs to complete..."
  for job_id in "${job_ids[@]}"; do
    wait_for_job_completion "$job_id" 10
  done
fi

# Step 3: GWAS analysis
if [[ $runGWAS == "TRUE" ]]; then
  echo "Starting GWAS analysis for each trait..."
  for i in "${!TRAIT_NAMES[@]}"; do
    TRAIT="${TRAIT_NAMES[$i]}"
    echo "Processing trait: $TRAIT"
    mkdir -p "${bash_out}/${NAME}/${TRAIT}/GWAS/"
    cd "${bash_out}/${NAME}/${TRAIT}/GWAS/"

    PHENOTYPE_FILE="${phenotypes_files[$i]}"
    PHENOTYPE_COLUMN="${pheno_pos[$i]}"
    CG_COLUMN="${CG_pos[$i]}"

    echo "Submitting GWAS job for ${TRAIT}..."
    gwas_job_id=$(sbatch --account="mateescu" \
                      --job-name="GWAS_${TRAIT}" \
                      --output="GWAS_${TRAIT}.out" \
                      --error="GWAS_${TRAIT}.err" \
                      "${my_bin}/scripts/GWAS/Marker_based_GWAS.sh" "${proj_env}" "${NAME}" "${TRAIT}" "${PHENOTYPE_FILE}" "${PHENOTYPE_COLUMN}" "${CG_COLUMN}" "${partitions}" "${window_size}" | awk '{print $4}')
    echo "Waiting for GWAS job ${gwas_job_id} to complete..."
    wait_for_job_completion "$gwas_job_id" 5
  done
fi

# Step 4: Gather_Results and plotting
if [[ $runGatherPlots == "TRUE" ]]; then
  echo "Starting Gather_Results and plotting for each trait..."
  for i in "${!TRAIT_NAMES[@]}"; do
    TRAIT="${TRAIT_NAMES[$i]}"
    echo "Processing trait: $TRAIT"
    mkdir -p "${bash_out}/${NAME}/${TRAIT}/Gather_Results/"
    cd "${bash_out}/${NAME}/${TRAIT}/Gather_Results/"

    echo "Submitting Gather_Results job for ${TRAIT}..."
    job_id=$(sbatch --account="mateescu" \
                      --job-name="Gather_Results_${TRAIT}" \
                      --output="Gather_Results_${TRAIT}.out" \
                      --error="Gather_Results_${TRAIT}.err" \
                      "${my_bin}/scripts/Gather_Results/gather_gwas_and_plot.sh" "${proj_env}" "${NAME}" "${TRAIT}" "${PHENOTYPE_FILE}" "${PHENOTYPE_COLUMN}" "${CG_COLUMN}" "${partitions}" "${window_size}" | awk '{print $4}')
    job_ids+=("$job_id")
  done
  echo "Waiting for all Plotting jobs to complete..."
  for job_id in "${job_ids[@]}"; do
    wait_for_job_completion "$job_id" 5
  done
fi

echo "Pipeline completed."
date
