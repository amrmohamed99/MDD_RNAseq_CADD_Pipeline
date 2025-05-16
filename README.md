# MDD_RNAseq_CADD_Pipeline
Integrated RNA-Seq analysis, Machine learning and virtual screening pipeline for Major Depressive Disorder using DESeq2, ClusterProfiler, RDKit, and AutoDock Vina 

This repository provides an end-to-end **bioinformatics and cheminformatics pipeline** to analyze RNA-Seq data from **Major Depressive Disorder (MDD)** patients, apply **machine learning classifiers**, and perform **in silico drug repurposing** using FDA-approved compounds on target proteins **BDNF** and **NOTCH1**.

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

MDD_RNAseq_CADD_Pipeline/
├── RNAseq_Analysis/
│   ├── DESeq2.R                  # Differential expression analysis
│   ├── GSEA.R                    # Enrichment analysis (GO, KEGG)
│   └── ML_research.ipynb        # ML classification of MDD vs. control
├── CADD_VirtualScreening/
│   ├── E_rdkit.sh               # Energy minimization using RDKit
│   ├── PDB2PDBQT                # OpenBabel-based file conversion
│   └── virtual_screening_super.sh # AutoDock Vina batch screening
├── .gitignore
├── requirements.txt
└── README.md


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

💊 Drug Repurposing Workflow
Step	Tool	Description
Energy Minimization	RDKit	Minimizes FDA drugs
File Conversion to PDBQT	OpenBabel	Prepares for docking
Virtual Screening	AutoDock Vina	Screens against BDNF, NOTCH1
Output	Log files + best docking poses
---

📦 Installation

conda create -n mdd_pipeline python=3.10
conda activate mdd_pipeline
pip install -r requirements.txt

For R dependencies:

install.packages(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "ggplot2"))

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
