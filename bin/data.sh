

#!/bin/bash
#SBATCH --account=mateescu
#SBATCH --job-name=Data
#SBATCH --output=Data_%j.out
#SBATCH --error=errorfile_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gzayas97@ufl.edu
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=7000mb
data=/blue/mateescu/gzayas97/Gabe_Thesis/5.Breed_Specific_GWAS_V2/data
SNP=/blue/mateescu/gzayas97/UF250K_Feb_2024
SNP_name=UFID_UF250K_Feb_2024_Illumina.ID_ARS
BO=/blue/mateescu/gzayas97/BOA_Estimation/results/Breed.of.Origin/UF.2024/Breed_of_Origin.files
BO_name=UF.2024.BO
cd $data
#Setting up phenotype files
cd $data/
mkdir -p Phenotypes
cd $data/Phenotypes
cd $data
mkdir -p Genotypes
cd $data/Genotypes/
mkdir -p SNP
cd $data/Genotypes/SNP
cp $SNP/${SNP_name}.* .


cd $data/Genotypes/
mkdir -p BOA
cd $data/Genotypes/BOA
cp $BO/${BO_name}.* .