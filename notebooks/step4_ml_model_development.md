# 🤖 Step 4 — ML Model Development  
**📅 Date:** 21 August 2025  
**👩‍🔬 Author:** Mohammed Farzana Begum  
**📁 Notebook Path:** notebooks/step4_ml_model_development.md  

---

## 🎯 Objective

Train supervised ML models to predict anticancer activity of meroterpenoids using filtered molecular descriptors. Evaluate performance using classification and regression metrics.

---

## 🧪 Environment Setup

```bash
conda activate meroterp-ml
pip install scikit-learn xgboost

import pandas as pd

# Load descriptors and labels
X = pd.read_csv("data/filtered_descriptors.csv")
y_df = pd.read_csv("data/curated_dataset.csv")

# Extract targets
y_class = y_df["Activity_Label"]  # Binary: 1 = Active, 0 = Inactive
y_reg = y_df["pIC50"]             # Continuous: potency

# Check alignment
print(X.shape, y_class.shape, y_reg.shape)

---

## 🧠 Model Training — Random Forest (Classification)

Train a Random Forest classifier to predict binary anticancer activity using filtered descriptors.

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report

# Split data
X_train, X_test, y_train, y_test = train_test_split(X, y_class, test_size=0.2, stratify=y_class, random_state=42)

# Train model
rf = RandomForestClassifier(n_estimators=500, random_state=42)
rf.fit(X_train, y_train)

# Predict and evaluate
y_pred = rf.predict(X_test)
y_proba = rf.predict_proba(X_test)[:, 1]

print("Accuracy:", accuracy_score(y_test, y_pred))
print("ROC-AUC:", roc_auc_score(y_test, y_proba))
print(classification_report(y_test, y_pred))
