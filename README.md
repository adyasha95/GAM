# Clinical Biostatistics Modeling  
Survival Analysis • Cox Models • Generalized Additive Models (GAM)

This repository contains modular R scripts for core biostatistical workflows used in clinical and biomedical research.  
The code demonstrates reproducible methods for **time-to-event analysis**, **hazard modeling**, and **nonlinear effect estimation** using **Generalized Additive Models (GAM)**.

> **Important:**  
> These scripts were originally developed as part of projects using **sensitive biomedical datasets** (clinical imaging, psychiatric cohort data, or patient-level outcomes).  
> **Due to GDPR and data-sharing restrictions, the original datasets cannot be made public.**  
>  
> The versions provided here include **generic, synthetic, or user-replaceable data-loading placeholders**, ensuring the methods are reproducible without exposing confidential data.

---

## 📦 Repository Contents

### 1. `survival_analysis_v2.R`
A complete time-to-event modeling pipeline using R.

**Includes:**
- Kaplan–Meier estimator  
- Log-rank test  
- Cox proportional hazards model (univariate + multivariate)  
- Hazard ratio extraction and visualization  
- Time-to-event prediction workflow  
- Clean separation of loading, modeling, and plotting functions  

**Key Packages:**  
`survival`, `survminer`, `dplyr`, `ggplot2`

---

### 2. `code_GAM_v2.r`
Script for fitting **Generalized Additive Models (GAM)** to clinical/biomedical variables.

**Includes:**
- Model specification with nonlinear smooth terms  
- Covariate adjustment  
- Partial effect plots (smooth splines)  
- Handling continuous and categorical predictors  
- Clear model summary + interpretation outputs  

**Key Packages:**  
`mgcv`, `ggplot2`, `dplyr`

---

## 🧠 Why These Methods Matter in Biomedical Research

### **Survival Analysis**
Used for:
- Clinical progression  
- Hospitalization risk  
- Transition prediction (e.g., CHR → psychosis)  
- Treatment response duration  
- Risk factor modeling  

Provides:
- Hazard ratios  
- Time-dependent probabilities  
- Survival curves  
- Clinically interpretable risk profiles  

---

### **Generalized Additive Models (GAM)**
Used for:
- Nonlinear biomarker–outcome relationships  
- Brain/clinical covariate effects  
- Age effects with nonlinearity  
- Modeling without unrealistic linear assumptions  

GAMs are especially helpful in:
- Neuroimaging  
- Clinical psychiatry  
- Biomarker discovery  
- Developmental and lifespan modeling  

---

## 🔐 Data Privacy & Compliance Notice

This repository includes **code only**, not data.  
The original analyses were performed on **GDPR-protected biomedical datasets**, including patient-level clinical, imaging, or questionnaire data.

To comply with:
- **GDPR** (EU)  
- **Institutional Ethics approvals**  
- **Data sharing agreements**  

➡️ **No raw data or identifiable variables can be shared.**

All current examples are built to work with:
- synthetic data  
- user-provided datasets  
- de-identified tables  

You may adapt the scripts to your own dataset structure.

---

## ▶️ How to Use

### **1. Install dependencies**
```r
install.packages(c("survival", "survminer", "mgcv", "ggplot2", "dplyr"))
---
