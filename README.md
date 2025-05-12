# Genome assembly (RAD-Seq), GEA Analyses & Evolutionary Rescue

This repo contains the core R and shell scripts used in:

> **Azevedo J. A. R. et al.**  
> *Deforestation limits evolutionary rescue under climate change in Amazonian lizards.*

> **Status:**  
> - These scripts are tailored to the original datasets and analyses.  
> - A cleaned, generalized version will be released soon or upon request.

---

## 📂 Script overview

### 1. `Variable selection for GEA clean.R`  
Performs genotype–environment association (GEA) analyses to identify candidate SNPs under selection.  
- **Methods implemented:**  
  - Redundancy Analysis (RDA)  
  - Latent Factor Mixed Models (LFMM)  
- **Outputs:**  
  - Lists of climate-associated SNPs  
  - Genotype subsets for downstream analyses  

### 2. `Modelling_conectivity_Ken_calc_v3.R`  
Models landscape connectivity and projects species distributions to simulate evolutionary rescue.  
- **Workflow:**  
  1. Load individuals classified by climate-associated genotypes  
  2. Fit Species Distribution Models (SDMs) under present and future climate scenarios  
  3. Run diffusion simulations across connectivity surfaces  
  4. Calculate evolutionary rescue metrics and generate spatial rescue maps  

---

## 🔧 Dependencies

- **R** ≥ 4.0  
- CRAN & Bioconductor packages (non-exhaustive):  
  - `vegan` (RDA)  
  - `lfmm` (GEA)  
  - `adegenet` / `vcfR` (genotype handling)  
  - `raster`, `sp`, `sf` (spatial data)  
  - `sdm` / `maxnet` (SDMs)  
  - `gdistance` (connectivity/diffusion)  

---

## 🚀 How to get a cleaner version

If you need a more generalized workflow or stripped-down scripts for your own data, please open an issue or contact me:

---

## 📄 License

All code in this repository is released under a **CC0 Public Domain Waiver**.  
Feel free to reuse, adapt, and cite as:

> Azevedo J. A. R. et al. (2025) Deforestation limits evolutionary rescue under climate change in Amazonian lizards. *Manuscript in review.*

