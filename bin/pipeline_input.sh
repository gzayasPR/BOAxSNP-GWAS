#!/bin/bash
# pipeline_input.sh: Configuration file for the pipeline.
# This file sets up all environment-specific variables, file paths, and parameters
# used by the pipeline run script. Adjust these settings as needed without modifying
# the execution logic.

# -----------------------------------------------------------
# Environment-specific Variables
# -----------------------------------------------------------
# Path to the project environment file containing additional settings.
proj_env="${my_bin}/project_env.sh"
# Name of the analysis/pipeline run.
NAME="Thermo_ALL"
# Path to the SNP genotype file used in the analysis.
SNPGeno="${my_data}/Genotypes/SNP/UFID_UF250K_Feb_2024_Illumina.ID_ARS"
# Path to the breed-of-origin genotype file.
BOGeno="/blue/mateescu/gzayas97/BOA_Estimation/results/BOA_UF.2024/UF.2024/Breed_of_Origin.files/UF.2024.BO"
# Number of partitions used to split the data for parallel processing.
partitions=100
# Window size parameter from BOA analsysis
window_size=15

# -----------------------------------------------------------
# Phenotype Files and Corresponding Column Positions
# -----------------------------------------------------------
# List of phenotype file paths for each trait.
phenotypes_files=(
  "${my_data}/Phenotypes/Thermo_ALL/Thermo_ALL_SCL.csv"   # Phenotype data for trait SCL
  "${my_data}/Phenotypes/Thermo_ALL/Thermo_ALL_SWA.csv"   # Phenotype data for trait SWA
  "${my_data}/Phenotypes/Thermo_ALL/Thermo_ALL_TSS2.csv"    # Phenotype data for trait TSS2
  "${my_data}/Phenotypes/Thermo_ALL/Thermo_ALL_LHL.csv"     # Phenotype data for trait LHL
)
# Define the column position in each phenotype file that contains the main phenotype data.
pheno_pos=(5 5 6 5)
# Define the column position for additional information (e.g., contemporary group as fixed effect) in each file.
CG_pos=(2 2 2 2)

# -----------------------------------------------------------
# Trait Names
# -----------------------------------------------------------
# Each trait name corresponds to a phenotype file above and is used for labeling outputs.
TRAIT_NAMES=( "SCL" "SWA" "TSS2" "LHL" )

# -----------------------------------------------------------
# Flags to Control Pipeline Steps
# -----------------------------------------------------------
# Set these flags to "TRUE" or "FALSE" to enable/disable specific analysis steps:
# runSetup: Run the initial setup (data preparation)
# runVARCOMP: Run the variance component analysis
# runGWAS: Run the genome-wide association study analysis
# runGatherPlots: Gather results and generate plots
runSetup="TRUE"
runVARCOMP="TRUE"
runGWAS="TRUE"
runGatherPlots="TRUE"
