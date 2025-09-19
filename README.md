
# GAM Trajectories (mgcv)

**Author:** Dr. Adyasha Tejaswi Khuntia  
**Date:** 19 Sep 2025  
**R version:** ≥ 4.1  
**Key packages:** `mgcv`, `ggplot2`, `dplyr`, `tidyr`, `readxl`, `patchwork`, `writexl`, `readr`

---

## 📌 Overview

This project fits **four Generalized Additive Models (GAMs)** using `mgcv::bam()` and generates **population-level predicted trajectories** (with 95% CIs) for cognitive and functional outcomes over time, by risk class. Plots are formatted to match the provided client reference figures.

### Responses
- **mPACC** → `zscore` (Gaussian)
- **MMSE** → `mlm` (Binomial logit, scaled to 0–30)
- **MEM** → `m_score_std` (Standardized version of `m_score`)
- **EXEC** → `func_score` (Gaussian)

---

## 📂 Data Requirements

Input file: **`cdf.xlsx`** with the following columns (exact names expected):

- **Outcomes**: `zscore`, `mlm`, `m_score`, `func_score`  
- **Predictors**: `class`, `gender`, `age_at_base`, `ap4`, `yr_educ`  
- **Time**: `yrs_from_base` (numeric or coercible to numeric)  
- **Random effect**: `patient` (subject ID)

---

## ⚙️ Model Specification

For each response:
- **Random effects:** patient-specific smooth.  
- **Population-level predictions:** random effects excluded.  
- **MMSE:** fit as `cbind(mlm, 30 - mlm)` with logit link; predictions inverse-linked and scaled to 0–30.

### 🔎 About `invlink`
- **Gaussian** → identity (η = y).  
- **Binomial logit** → η = log(p/(1-p)); p = 1/(1+exp(-η)); then scaled as `p*30`.  
- This restores predictions to the **original MMSE scale**.

### 📊 MEM Standardization
- `m_score` is standardized relative to the **baseline distribution** (per-patient minimum time).  
- This yields `m_score_std`, so baseline ~0 and decline mirrors client figures.

---

## 🛠️ Installation

Install dependencies:

```r
install.packages(c(
  "mgcv", "ggplot2", "dplyr", "tidyr", "readxl", 
  "patchwork", "writexl", "readr"
), repos = "https://cloud.r-project.org")
