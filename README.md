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
│ ├── raw_meroterpenoids.csv 
│ └── curated_smiles.csv 
├── notebooks/ 
│ ├── step1_literature_mining.ipynb 
│ ├── step2_smiles_cleaning.ipynb 
│ └── step3_descriptor_setup.ipynb


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

---

## 📦 Input Files

| File | Description |
|------|-------------|
| `raw_meroterpenoids.csv` | Raw SMILES extracted from literature |
| `curated_smiles.csv` | Canonical SMILES after RDKit standardization |

---

## 💼 Highlights

- Hands-on RDKit, Pandas, and SMILES standardization  
- Fully documented pipeline with Markdown and comments  
- Clean GitHub structure for reproducibility  
- Ready for descriptor generation and ML modeling in next steps
