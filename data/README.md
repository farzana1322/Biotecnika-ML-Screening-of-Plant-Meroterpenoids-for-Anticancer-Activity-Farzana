# 📁 Data Folder

This folder contains curated SMILES, descriptor matrices, IC50 values, and raw compound data.  
These files are essential for reproducibility, modeling, and recruiter-grade clarity.

## 🔒 Do Not Delete
- These files are restored from backups and protected by version control.
- Every molecule matters. Every descriptor is sponsor-grade.

## 🧬 Contents
- `curated_smiles.csv` – Canonicalized molecular structures
- `descriptor_matrix.csv` – Morgan fingerprints (radius=2, size=1024)
- `ic50_values.csv` – Bioactivity values for regression
- `ic50_values_resimulated.csv` – Refined scaffolded IC50 values

This folder is locked by ritual. Treat it with respect.

`ic50_scaffolded.csv` contains merged descriptors and IC50 values used for regression modeling (Step 9).  
Generated via feature selection and Gradient Boosting. See main README for details.
