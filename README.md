# Genetic Architecture of Thermotolerance Traits in Angus-Brahman Crossbred Cattle

**Gabriel A. Zayas, et al.**  
Department of Animal Sciences, University of Florida  

This repository contains the code, data processing pipeline, and analysis scripts for the study:

> **"Novel Insights into the Genetic Architecture of Thermotolerance Traits in Angus-Brahman Crossbreeds Using SNP and BOA-Based Analyses"**

---

## 🧬 Overview

This study investigates the genetic basis of thermotolerance traits in crossbred Angus-Brahman cattle using both standard SNP markers and Breed of Origin of Allele (BOA) information. A total of 3,962 animals were phenotyped and genotyped, and genome-wide association studies (GWAS) were performed using:

- **Model 1**: SNP-only  
- **Model 2**: BOA-only  
- **Model 3**: Integrated SNP + BOA  
- All models implemented with a **Leave-One-Chromosome-Out (LOCO)** genomic relationship matrix to prevent proximal contamination during QTL detection.

---

## 🧪 Phenotypes

Four thermotolerance-related traits were analyzed:

- **SHL** – Short Hair Length  
- **LHL** – Long Hair Length  
- **SWA** – Sweat Gland Area  
- **TSS** – Thermal Stress Slope  
  *(calculated from body temperature response to daily thermal load using iButton loggers and THI data)*

---

## 🧬 Genotyping and BOA Assignment

- Genotyped with the **GGP F250 array (221K SNPs)**.  
- QC performed using **PLINK2**.  
- BOA assignment performed using **LAMP-LD**, with 123 Angus and 406 Brahman purebreds as the reference populations.  
- BOA genotypes encoded as:  
  `0` = Angus, `1` = Heterozygous, `2` = Brahman  

For details on the BOA assignment pipeline, please refer to the companion repository:  
👉 [`BOA_Estimation`](https://github.com/gzayasPR/BOA_Estimation)

---

## 🧠 GWAS Analysis

- GWAS performed using **linear mixed models** via the `lmm.diago()` function in the **gaston** R package.
- Models included fixed effects (SNP, BOA), random animal effects, and a **LOCO-based GRM** to control for confounding due to relatedness.
- Marker significance was tested via Wald statistics.  
- Bonferroni correction (α = 0.1) and suggestive thresholds (1/N markers) were applied.

---

## 📁 Repository Structure

```
BOAxSNP-GWAS/
├── data/                 # Input phenotypes, genotypes, BOA outputs
├── scripts/              # R scripts for QC, LOCO GRM generation, GWAS, plotting
├── results/              # GWAS outputs, summary tables
└── README.md             # This file
```

---

## 📚 Citation

If you use this pipeline or build on this work, please cite:  
Zayas et al., *"Novel Insights into the Genetic Architecture of Thermotolerance Traits in Angus-Brahman Crossbreeds Using SNP and BOA-Based Analyses"*, **Genetics Selection and Evolution**, (In Review).

---

## 📂 Data Availability

Phenotypes and genotypes used in this study are **available upon reasonable request**.  
Please contact the corresponding author (below) for data access.

---

## 🧰 Tools and Software

- `R v4.4`, [`gaston`](https://cran.r-project.org/package=gaston), `data.table`, `ggplot2`  
- [`PLINK2`](https://www.cog-genomics.org/plink/2.0/)  
- [`BOA_Estimation`](https://github.com/gzayasPR/BOA_Estimation) (BOA assignment with LAMP-LD)

---

## 📬 Contact

**Gabriel A. Zayas Santiago**  
gzayas97@ufl.edu  
