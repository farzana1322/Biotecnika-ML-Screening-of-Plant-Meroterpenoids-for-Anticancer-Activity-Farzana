# 🧬 Step 3 — Descriptor Calculation  
**📅 Date:** 21 August 2025  
**👩‍🔬 Author:** Mohammed Farzana Begum  
**📁 Notebook Path:** notebooks/step3_descriptor_calculation.md  

---

## 🎯 Objective

Generate molecular descriptors for curated meroterpenoids using RDKit and Mordred. These descriptors will serve as input features for ML modeling in Step 4.

---

## 🧪 Environment Setup

```bash
conda activate meroterp-ml
pip install mordred

from rdkit import Chem
from mordred import Calculator, descriptors
import pandas as pd

# Initialize Mordred calculator
calc = Calculator(descriptors, ignore_3D=True)

# Define descriptor function
def compute_descriptors(smiles_list):
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    df = calc.pandas(mols)
    return df

# Load curated SMILES
df_smiles = pd.read_csv("data/curated_smiles.csv")  # Adjust path if needed
smiles_list = df_smiles["SMILES"].dropna().tolist()

# Compute descriptors
descriptor_df = compute_descriptors(smiles_list)

# Save to CSV
descriptor_df.to_csv("data/descriptors.csv", index=False)

# Inspect the descriptor matrix
print(descriptor_df.head())
print(descriptor_df.shape)

---

## 🧹 Descriptor Filtering Strategy

- **Missing Value Handling:**  
  - Drop descriptors with >30% missing values  
  - Impute remaining missing values using median strategy

- **Low Variance Removal:**  
  - Remove descriptors with near-zero variance (threshold = 1e-4)

- **Scaling:**  
  - Apply StandardScaler to normalize descriptor ranges

- **Output:**  
  - `filtered_descriptors.csv` — cleaned and scaled descriptor matrix  
  - `descriptor_metadata.md` — notes on removed features and preprocessing steps
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold

# Step 1: Drop descriptors with >30% missing values
missing_thresh = 0.3
valid_cols = descriptor_df.columns[descriptor_df.isnull().mean() < missing_thresh]
filtered_df = descriptor_df[valid_cols]

# Step 2: Impute missing values (median)
imputer = SimpleImputer(strategy='median')
imputed_array = imputer.fit_transform(filtered_df)

# Step 3: Remove low-variance descriptors
vt = VarianceThreshold(threshold=1e-4)
reduced_array = vt.fit_transform(imputed_array)

# Step 4: Scale descriptors
scaler = StandardScaler()
scaled_array = scaler.fit_transform(reduced_array)

# Step 5: Save cleaned matrix
filtered_df_final = pd.DataFrame(scaled_array, columns=[f'F{i}' for i in range(scaled_array.shape[1])])
filtered_df_final.to_csv("data/filtered_descriptors.csv", index=False)

---

## ✅ Status

Step 3 completed successfully. Descriptor matrix cleaned, scaled, and saved as `filtered_descriptors.csv`. Ready to proceed to **Step 4: ML Model Development**.

