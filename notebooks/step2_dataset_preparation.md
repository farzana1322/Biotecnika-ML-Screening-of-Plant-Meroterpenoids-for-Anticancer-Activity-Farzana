# 🧪 Step 2 — Dataset Preparation & Annotation  
**📅 Date:** 20 August 2025  
**👩‍🔬 Author:** Mohammed Farzana Begum  
**📁 Notebook Path:** notebooks/step2_dataset_preparation.md  

---

## 📚 Data Sources

| Source     | Type                      | Access Date | Notes |
|------------|---------------------------|-------------|-------|
| PubChem    | Chemical + Literature     | 18 Aug 2025 | SMILES, CID, IC₅₀ values |
| NPASS      | Bioactivity Database      | 18 Aug 2025 | IC₅₀, EC₅₀, assay metadata |
| COCONUT    | Natural Products Database | 19 Aug 2025 | Meroterpenoid scaffolds |
| ChEMBL     | Bioactivity               | 19 Aug 2025 | Target-specific IC₅₀ values |
| PubMed     | Literature                | 17–20 Aug 2025 | Mechanistic insights |
| Scopus     | Literature Index          | 17–20 Aug 2025 | Compound-target mapping |

---

## 🧬 Environment Setup

```bash
conda create -n meroterp-ml python=3.10 -c conda-forge rdkit pandas scikit-learn
conda activate meroterp-ml

from rdkit import Chem
from rdkit.Chem import rdMolStandardize

def standardize_smiles(smiles):
    m = Chem.MolFromSmiles(smiles)
    normalizer = rdMolStandardize.LargestFragmentChooser()
    m = normalizer.choose(m)
    Chem.SanitizeMol(m)
    uncharger = rdMolStandardize.Uncharger()
    m = uncharger.uncharge(m)
    can_smiles = Chem.MolToSmiles(m, canonical=True)
    inchi = Chem.MolToInchi(m)
    inchikey = Chem.InchiToInchiKey(inchi)
    return can_smiles, inchikey

---

## 🎯 Activity Labeling Strategy

- **Binary Classification:**  
  - Active if IC₅₀ ≤ 10 µM  
  - Inactive if IC₅₀ > 10 µM or not reported

- **Regression Target:**  
  - Convert IC₅₀ to pIC₅₀ using:  
    $$ \text{pIC}_{50} = -\log_{10}(\text{IC}_{50} \text{ in M}) $$  
    > **Note:** This transformation is standard in pharmacology and cheminformatics. It is derived from the internship PDF’s labeling strategy, which recommends converting IC₅₀ to pIC₅₀ using log₁₀ scale. Example: IC₅₀ = 1 µM = \(1 \times 10^{-6}\) M → pIC₅₀ = 6.

---

## 📁 Output Files

| Filename               | Description |
|------------------------|-------------|
| `curated_smiles.csv`   | Cleaned and standardized SMILES with InChIKeys |
| `curated_dataset.csv`  | Final annotated dataset with activity labels |
| `metadata.csv`         | Assay details, source organisms, and references |

---

## ✅ Status

Step 2 completed successfully using curated compounds from Step 1. Ready to proceed to **Step 3: Descriptor Calculation**.

