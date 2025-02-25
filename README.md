# ILC-Spatial-Transcriptomics

This repository contains original functions and scripts for analyzing the spatial transcriptomics (ST) data presented in the publication:

**Serra M. et al., [Publication Title] ([Link])**

The repository is organized into the following folders:

---

## 📁 Folder Structure

### **Python**  
Contains Python scripts for:
- Transferring morphological annotations from the annotated H&E slide to the ST spot level.

### **R**  
Contains R scripts for:
- Pre-processing and normalization of ST data.
- Co-occurrence analysis using morphological data (from annotation) at the ST spot level.
- Identification of ILC tumor microenvironment (TME)-based subtypes from RNA-seq data.

### **input_data**  
Contains:
- TME-based subtype-related gene signatures used for RNA-seq data classification.

---

## 📊 Data Availability
All other data required to replicate the analyses and use these scripts can be accessed on **Zenodo**:  
https://zenodo.org/records/14924871 

---

## 🛠️ Software Requirements
The scripts in this repository have been tested using the following software versions:
- **R**: v4.4.0  
- **Python**: v3.9.1  

---

## 🔑 ILC TME-Based Subtypes Dictionary

| Subtype | Description                 | Hex Color Code |
|---------|-----------------------------|----------------|
| **NSE** | Normal/Stroma-Enriched      | `#F57C00`      |
| **P**   | Proliferative               | `#2E7D32`      |
| **ARE** | Androgen Receptor-Enriched  | `#1E88E5`      |
| **MIE** | Metabolic/Immune-Enriched   | `#8E24AA`      |

---

## 📢 Citation
If you use this classification and these scripts in your research, please cite our work:

> **Serra M., et al., [Title of the Paper]**  
> **[Publication Link](#)**

---

## 📘 Further Information
This analysis and the identification of ILC TME-based subtypes were conducted by **Matteo Serra** at the **Breast Cancer Translational Research Laboratory (BCTL), Institut Jules Bordet, ULB**, under the supervision of **Christos Sotiriou (MD, PhD)**.

For inquiries, feel free to contact matteo95serra@gmail.com.

