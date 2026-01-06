"""
Verify Supervised Learning Claims for Thesis Manuscript
Runs Random Forest classification on cluster labels and outputs metrics.
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
from pathlib import Path

# Load data
data_path = Path(__file__).parent / "data" / "processed" / "clustered_dataset.csv"
df = pd.read_csv(data_path)

# Define feature columns (22 encoded antibiotics)
# Get all columns with '_encoded' suffix
encoded_cols = [c for c in df.columns if '_encoded' in c]
print(f"Using {len(encoded_cols)} encoded features: {encoded_cols}")

# Antibiotic name mapping for display
antibiotic_names = [c.replace('_encoded', '') for c in encoded_cols]

# Prepare X and y
X = df[encoded_cols].values
y = df['CLUSTER'].values  # Target is cluster assignment

# Random Forest with 5-fold CV
rf = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

# Cross-validated predictions
y_pred = cross_val_predict(rf, X, y, cv=cv)

# Compute metrics
accuracy = accuracy_score(y, y_pred)
precision_macro = precision_score(y, y_pred, average='macro')
recall_macro = recall_score(y, y_pred, average='macro')
f1_macro = f1_score(y, y_pred, average='macro')
precision_weighted = precision_score(y, y_pred, average='weighted')
recall_weighted = recall_score(y, y_pred, average='weighted')

print("\n" + "="*60)
print("RANDOM FOREST CLUSTER VALIDATION RESULTS")
print("="*60)
print(f"Accuracy:           {accuracy:.3f} ({accuracy*100:.1f}%)")
print(f"Macro F1-Score:     {f1_macro:.3f}")
print(f"Weighted Precision: {precision_weighted:.2f}")
print(f"Weighted Recall:    {recall_weighted:.2f}")
print("="*60)

# Train full model for feature importance
rf_full = RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=-1)
rf_full.fit(X, y)

# Feature importance
importances = rf_full.feature_importances_
importance_df = pd.DataFrame({
    'Antibiotic': antibiotic_names,
    'Importance': importances
}).sort_values('Importance', ascending=False)

print("\nTOP 10 FEATURE IMPORTANCES:")
print("-"*40)
for i, row in importance_df.head(10).iterrows():
    print(f"  {row['Antibiotic']:4s}: {row['Importance']:.4f}")
print("-"*40)

# Save results
output_path = Path(__file__).parent / "rf_cluster_validation_results.txt"
with open(output_path, 'w') as f:
    f.write("RANDOM FOREST CLUSTER VALIDATION RESULTS\n")
    f.write("="*50 + "\n")
    f.write(f"Accuracy:           {accuracy:.4f} ({accuracy*100:.1f}%)\n")
    f.write(f"Macro F1-Score:     {f1_macro:.4f}\n")
    f.write(f"Macro Precision:    {precision_macro:.4f}\n")
    f.write(f"Macro Recall:       {recall_macro:.4f}\n")
    f.write(f"Weighted Precision: {precision_weighted:.4f}\n")
    f.write(f"Weighted Recall:    {recall_weighted:.4f}\n")
    f.write("\nTOP 10 FEATURE IMPORTANCES:\n")
    f.write("-"*40 + "\n")
    for i, row in importance_df.head(10).iterrows():
        f.write(f"  {row['Antibiotic']:4s}: {row['Importance']:.4f}\n")

print(f"\nResults saved to: {output_path}")
