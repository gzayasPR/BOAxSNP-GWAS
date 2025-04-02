# BOAxSNP-GWAS

## Genetic Architecture of Thermotolerance Traits in Angus-Brahman Crossbred Cattle

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
- BOA assignment performed using **LAMP-LD** with 123 Angus and 406 Brahman purebreds.  
- BOA genotypes encoded as:  
  `0` = Angus, `1` = Heterozygous, `2` = Brahman  

---

## 🧠 GWAS Analysis

- GWAS performed with **linear mixed models** using the `lmm.diago()` function in the **gaston** R package.
- Models include fixed effects (SNP, BOA), random animal effects, and a **chromosome-specific GRM** excluding the chromosome of interest (**LOCO**).
- Marker significance assessed via Wald tests.  
- Correction for multiple testing with Bonferroni (α = 0.1) and a suggestive threshold of 1/N markers.

---

## 📁 Repository Structure

```
BOA_Thermotolerance/
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

## 🧰 Tools and Software

- `R v4.4`, [`gaston`](https://cran.r-project.org/package=gaston), `data.table`, `ggplot2`  
- [`PLINK2`](https://www.cog-genomics.org/plink/2.0/)  
- [`LAMP-LD`](https://github.com/gzayasPR/BOA_Estimation) (for BOA estimation)

---

## 📬 Contact

**Gabriel A. Zayas Santiago**  
gzayas97@ufl.edu  

