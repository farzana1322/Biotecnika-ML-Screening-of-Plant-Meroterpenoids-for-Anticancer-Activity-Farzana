# 🧬 Biotecnika ML Screening of Plant Meroterpenoids for Anticancer Activity  
A reproducible machine learning pipeline for in silico screening of plant-derived meroterpenoids using RDKit, Pandas, and Python.

**Author:** Mohammed Farzana Begum  
**Start Date:** August 2025  
**Location:** Guntur, Andhra Pradesh, India (Remote Development)

---

## 🎯 Objective  
To build a sponsor-grade ML pipeline for anticancer screening of meroterpenoids through literature mining, molecular standardization, descriptor generation, and predictive modeling.

---

## 📁 Project Structure

Biotecnika-ML-Screening/  
├── data/  
│   ├── raw_meroterpenoids.csv  
│   ├── curated_smiles.csv  
│   ├── descriptor_matrix.csv  
│   └── ic50_scaffolded.csv  
├── models/  
│   └── random_forest_model.joblib  
├── notebooks/  
│   ├── step1_literature_mining.ipynb  
│   ├── step2_smiles_cleaning.ipynb  
│   ├── step3_descriptor_setup.ipynb  
│   ├── step5_virtual_screening.ipynb  
│   └── step9_ic50_modeling.ipynb  
├── results/  
│   └── screening_results.csv  
├── backups/  
│   └── Biotecnika_Release_0609.zip  

---

## 🧪 Pipeline Summary

### Step 1 — Literature Mining  
- Extracted SMILES and compound IDs from published sources  
- Saved to `raw_meroterpenoids.csv`

### Step 2 — SMILES Cleaning  
- Removed duplicates and invalid entries  
- Standardized naming and formatting

### Step 3 — Canonical SMILES Standardization  
- Applied RDKit SaltRemover and canonicalization  
- Saved to `curated_smiles.csv`

### Step 3B — Descriptor Matrix Generation  
- Generated 1024-bit ECFP fingerprints using RDKit  
- Saved to `descriptor_matrix.csv`

### Step 4 — Model Training  
- Created dummy labels for supervised learning simulation  
- Trained Random Forest model (200 trees, max depth=10)  
- Saved to `models/random_forest_model.joblib`

### Step 5 — Virtual Screening  
- Predicted anticancer probabilities  
- Ranked results saved to `results/screening_results.csv`

---

## 📊 Step 9 — IC50 Regression Modeling

Evaluates the predictive relationship between descriptors and IC50 bioactivity.

### 🔬 Data Sources
- `descriptor_matrix.csv` – Morgan fingerprints  
- `ic50_values.csv` – Literature-curated IC50 values  
- `ic50_scaffolded.csv` – Merged descriptor + IC50 matrix

### 🧠 Modeling Strategy
- Feature selection via Pearson correlation  
- Top 20 descriptors selected  
- Model: GradientBoostingRegressor (n_estimators=100, learning_rate=0.1)

### 📈 Results
- **MSE**: 0.1895  
- **R² Score**: 0.0436 (positive, weak signal)

### 🧪 Interpretation
- Signal present but weak; descriptors show low IC50 correlation  
- Future improvement via MolWt, LogP, TPSA, and ensemble models

### 🔐 Reproducibility
- All files backed up in `Biotecnika_Release_0609.zip`  
- Pipeline reproducible via `step9_ic50_modeling.ipynb`

---

## 📦 Input Files

| File | Description |
|------|-------------|
| `raw_meroterpenoids.csv` | Raw SMILES from literature |
| `curated_smiles.csv` | Canonicalized SMILES |
| `descriptor_matrix.csv` | Fingerprint-based descriptors |
| `ic50_values.csv` | Bioactivity values |
| `ic50_scaffolded.csv` | Merged descriptors + IC50 |

---

## 🔧 How to Run

1. Clone the repo  
2. Open notebooks in Jupyter or VS Code  
3. Run each step sequentially  
4. Final predictions saved to `results/screening_results.csv`

---

## 💼 Highlights

- Hands-on RDKit, Pandas, and SMILES standardization  
- Fully documented pipeline with Markdown and comments  
- Clean GitHub structure for reproducibility  
- Ready for descriptor generation, ML modeling, and virtual screening  
- In silico predictions prepared for ADMET filtering and PAINS removal

---

## 🚀 Next Steps

- Apply ADMET filtering via SwissADME or pkCSM  
- Remove PAINS using RDKit or FAF-Drugs4  
- Validate top hits with literature bioactivity  
- Prepare RA documentation for sponsor-facing review  
- Extend pipeline to include IC50 classification and ensemble modeling

---

## 💚 Author Note

Mohammed Farzana Begum – Clinical AI Developer, ethics-certified researcher, and mentor.    
Every molecule matters. Every line of code is a step toward sponsor-grade clarity.
