# Thermotolerance GWAS Pipeline

This repository contains the configuration and scripts for running genome-wide association studies (GWAS) on thermotolerance traits in Angus-Brahman crossbred cattle using both SNP and Breed-of-Origin (BOA) genotypes.

The pipeline is controlled via a configurable input file: `pipeline_input.sh`, which allows easy customization of trait-specific phenotype files, genotype sources, and analysis steps without modifying core logic.

---

## 🔧 Configuration Overview: `pipeline_input.sh`

This file defines key variables and flags for running the pipeline.

### 🔬 Project Settings

- **Project Name:** `Thermo_ALL`
- **Project Environment File:** `project_env.sh`
- **SNP Genotypes:** `${my_data}/Genotypes/SNP/UFID_UF250K_Feb_2024_Illumina.ID_ARS`
- **BOA Genotypes:** `/blue/mateescu/gzayas97/BOA_Estimation/results/BOA_UF.2024/UF.2024/Breed_of_Origin.files/UF.2024.BO`
- **Window Size (BOA):** `15`
- **Partitions for Parallel Processing:** `100`

> **BOA assignment** is based on [BOA_Estimation](https://github.com/gzayasPR/BOA_Estimation), a custom pipeline using LAMP-LD for local ancestry inference and breed-of-origin encoding.

---

### 📂 Phenotypes

The pipeline supports four thermotolerance traits:

| Trait Name | File | Phenotype Column | CG Column |
|------------|------|------------------|-----------|
| `SCL`  | `Thermo_ALL_SCL.csv`  | Column 5 | Column 2 |
| `SWA`  | `Thermo_ALL_SWA.csv`  | Column 5 | Column 2 |
| `TSS2` | `Thermo_ALL_TSS2.csv` | Column 6 | Column 2 |
| `LHL`  | `Thermo_ALL_LHL.csv`  | Column 5 | Column 2 |

Each trait file is provided as a `.csv` and should follow a consistent structure with defined phenotype and fixed effect (contemporary group) columns.

---

### ⚙️ Pipeline Steps

The following flags enable/disable stages of the pipeline:

| Step            | Description                         | Flag Name         |
|-----------------|-------------------------------------|-------------------|
| Setup           | Data formatting and checks          | `runSetup`        |
| Variance Component Analysis | Fit null models for each trait | `runVARCOMP`     |
| GWAS            | Perform GWAS across partitions      | `runGWAS`         |
| Gather + Plot   | Merge results and generate figures  | `runGatherPlots`  |

All steps are **enabled by default** in the example configuration:
```bash
runSetup="TRUE"
runVARCOMP="TRUE"
runGWAS="TRUE"
runGatherPlots="TRUE"
```
## 📬 Contact

**Gabriel A. Zayas Santiago**  
gzayas97@ufl.edu  
