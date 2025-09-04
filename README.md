# 🧬 Biotecnika ML Screening of Plant Meroterpenoids for Anticancer Activity  
A machine learning pipeline for anticancer screening of plant-derived meroterpenoids using RDKit and Python.

**Author:** Mohammed Farzana Begum  
**Start Date:** August 2025  
**Location:** Guntur, Andhra Pradesh, India (remote project development)

---

## 🎯 Objective  
To build a reproducible ML pipeline for screening plant-derived meroterpenoids for anticancer activity using literature mining, molecular standardization, and predictive modeling.

---

## 📁 Project Structure

Biotecnika-ML-Screening/  
├── data/  
│   ├── raw_meroterpenoids.csv  
│   ├── curated_smiles.csv  
│   └── descriptor_matrix.csv  
├── models/  
│   └── random_forest_model.joblib  
├── notebooks/  
│   ├── step1_literature_mining.ipynb  
│   ├── step2_smiles_cleaning.ipynb  
│   ├── step3_descriptor_setup.ipynb  
│   └── step5_virtual_screening.ipynb  
├── results/  
│   └── screening_results.csv  

---

## 🧪 Pipeline Summary

### Step 1 — Literature Mining  
- Extracted SMILES and compound IDs from published sources  
- Saved to `raw_meroterpenoids.csv`

### Step 2 — SMILES Cleaning  
- Removed duplicates and invalid entries  
- Standardized naming and formatting

### Step 3 — Canonical SMILES Standardization  
- Standardized SMILES using RDKit SaltRemover  
- Generated canonical SMILES for valid compounds  
- Saved to `curated_smiles.csv`

### Step 3B — Descriptor Matrix Generation  
- Computed physicochemical descriptors using RDKit  
- Generated 1024-bit ECFP fingerprints for each compound  
- Saved ML-ready feature matrix to `descriptor_matrix.csv`

  ### Step 4 — Model Training  
- Created balanced dummy labels for initial testing  
- Split data into training and testing sets using `train_test_split`  
- Trained Random Forest model with 200 trees and max depth of 10  
- Saved trained model to `models/random_forest_model.joblib`  
- This step simulates supervised learning and prepares the model for virtual screening


### Step 5 — Virtual Screening  
- Loaded descriptor matrix from `data/descriptor_matrix.csv`  
- Trained Random Forest model using balanced dummy labels  
- Saved model to `models/random_forest_model.joblib`  
- Predicted anticancer probabilities for each compound  
- Saved ranked results to `results/screening_results.csv`  
- This step simulates in silico screening and prepares compounds for ADMET filtering and PAINS removal

---

## 📦 Input Files

| File | Description |
|------|-------------|
| `raw_meroterpenoids.csv` | Raw SMILES extracted from literature |
| `curated_smiles.csv` | Canonical SMILES after RDKit standardization |
| `descriptor_matrix.csv` | ML-ready feature matrix with descriptors and fingerprints |

---

## 🔧 How to Run

1. Clone the repo  
2. Open `notebooks/step5_virtual_screening.ipynb` in Jupyter or VS Code  
3. Run each cell stepwise  
4. Final predictions will be saved to `results/screening_results.csv`

---

## 💼 Highlights

- Hands-on RDKit, Pandas, and SMILES standardization  
- Fully documented pipeline with Markdown and comments  
- Clean GitHub structure for reproducibility  
- Ready for descriptor generation, ML modeling, and virtual screening  
- In silico predictions saved for downstream filtering and validation

---

## 🚀 Next Steps

- Apply ADMET filtering using SwissADME or pkCSM  
- Remove PAINS using RDKit or FAF-Drugs4  
- Validate top hits with literature bioactivity data  
- Prepare RA documentation for sponsor-facing review  
- Extend pipeline to include regression models for IC50 prediction

---

