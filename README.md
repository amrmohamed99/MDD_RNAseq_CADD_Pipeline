# MDD_RNAseq_CADD_Pipeline
Integrated RNA-Seq analysis, Machine learning and virtual screening pipeline for Major Depressive Disorder using DESeq2, ClusterProfiler, RDKit, and AutoDock Vina 

This repository contains a complete pipeline integrating **RNA-Seq analysis**, **Machine Learning (ML)**-based classification, and **Computer-Aided Drug Design (CADD)** targeting **BDNF** and **NOTCH1** for **Major Depressive Disorder (MDD)**.

---

## 🔬 Overview

1. **RNA-Seq Analysis** (R)
   - Identification of Differentially Expressed Genes (DEGs) using DESeq2.
   - Functional enrichment analysis using ClusterProfiler.

2. **Machine Learning for MDD Prediction** (Python)
   - Applied supervised ML classifiers (e.g., Random Forest, SVM, Logistic Regression) on gene expression data.
   - Goal: Classify patients as **MDD vs Control** based on their transcriptomic profile.

3. **CADD Virtual Screening Pipeline** (Python/Shell)
   - Energy minimization of FDA-approved drugs using RDKit.
   - Conversion of drug structures from PDB to PDBQT using OpenBabel.
   - Virtual screening with AutoDock Vina against BDNF and NOTCH1 proteins.

---

## 📁 Folder Structure

RNAseq_Analysis/
│ ├── DESeq2_analysis.R
│ ├── enrichment_analysis.R
│ └── ML_classification_MDD.ipynb

CADD_VirtualScreening/
│ ├── energy_minimization_rdkit.py
│ ├── pdb_to_pdbqt_conversion.sh
│ └── virtual_screening_vina.sh



---

## 🚀 Usage

### 1. RNA-Seq Analysis
- Run `DESeq2_analysis.R` to identify differentially expressed genes.
- Perform enrichment analysis using `enrichment_analysis.R`.

### 2. Machine Learning Classification
- Run `ML_classification_MDD.ipynb` to train and evaluate classifiers on expression data.

### 3. CADD Virtual Screening
- Minimize FDA drug structures using `energy_minimization_rdkit.py`.
- Convert structures to PDBQT format using `pdb_to_pdbqt_conversion.sh`.
- Run docking with `virtual_screening_vina.sh`.

---

## 📦 Dependencies

### R Packages
- DESeq2
- ClusterProfiler
- org.Hs.eg.db

### Python Packages (see `requirements.txt`)
- rdkit
- openbabel
- scikit-learn
- pandas
- numpy
- matplotlib / seaborn

---

## 📊 Machine Learning Methods

- Random Forest
- Support Vector Machine (SVM)
- Logistic Regression
- ROC Curve, Accuracy, Precision, and Recall for performance evaluation

---
