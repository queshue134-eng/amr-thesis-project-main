"""
Split Ratio Sensitivity Analysis for Cluster Discrimination
This script runs actual simulations to compare different train-test splits
and cross-validation configurations.
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import f1_score, accuracy_score
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Load data
project_root = Path(__file__).parent.parent
data_path = project_root / "data" / "processed" / "clustered_dataset.csv"

print("=" * 70)
print("SPLIT RATIO SENSITIVITY ANALYSIS")
print("Cluster Discrimination Task")
print("=" * 70)

df = pd.read_csv(data_path)
feature_cols = [c for c in df.columns if c.endswith('_encoded')]
X = df[feature_cols].values
y = df['CLUSTER'].values

print(f"\nDataset: {len(df)} samples, {len(feature_cols)} features")
print(f"Clusters: {np.unique(y)}")
print(f"Cluster distribution: {dict(zip(*np.unique(y, return_counts=True)))}")

# Models
models = {
    'Logistic Regression': LogisticRegression(max_iter=1000, random_state=42),
    'Random Forest': RandomForestClassifier(n_estimators=100, random_state=42),
    'KNN': KNeighborsClassifier(n_neighbors=5)
}

# Test different split ratios
split_ratios = [0.30, 0.20, 0.10]  # test sizes
seeds = [42, 123, 456, 789, 1011]

print("\n" + "=" * 70)
print("PART 1: SPLIT RATIO COMPARISON")
print("=" * 70)

split_results = []

for test_size in split_ratios:
    train_pct = int((1 - test_size) * 100)
    test_pct = int(test_size * 100)
    split_name = f"{train_pct}/{test_pct}"
    
    print(f"\n--- {split_name} Split ---")
    
    for model_name, model_template in models.items():
        f1_scores = []
        acc_scores = []
        
        for seed in seeds:
            # Split data
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=test_size, random_state=seed, stratify=y
            )
            
            # Impute and Scale
            imputer = SimpleImputer(strategy='median')
            X_train_imputed = imputer.fit_transform(X_train)
            X_test_imputed = imputer.transform(X_test)
            
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train_imputed)
            X_test_scaled = scaler.transform(X_test_imputed)
            
            # Train model (fresh instance)
            model = model_template.__class__(**model_template.get_params())
            model.fit(X_train_scaled, y_train)
            
            # Evaluate
            y_pred = model.predict(X_test_scaled)
            f1 = f1_score(y_test, y_pred, average='macro', zero_division=0)
            acc = accuracy_score(y_test, y_pred)
            
            f1_scores.append(f1)
            acc_scores.append(acc)
        
        f1_mean = np.mean(f1_scores)
        f1_std = np.std(f1_scores)
        acc_mean = np.mean(acc_scores)
        
        split_results.append({
            'Split': split_name,
            'Model': model_name,
            'F1_Mean': f1_mean,
            'F1_Std': f1_std,
            'Accuracy': acc_mean,
            'Test_Samples': int(len(df) * test_size)
        })
        
        print(f"  {model_name}: F1={f1_mean:.3f} ± {f1_std:.3f}, Acc={acc_mean:.3f}")

# Cross-validation comparison
print("\n" + "=" * 70)
print("PART 2: CROSS-VALIDATION COMPARISON")
print("=" * 70)

cv_results = []
cv_configs = [5, 10]

# Impute and scale entire dataset for CV
imputer = SimpleImputer(strategy='median')
X_imputed = imputer.fit_transform(X)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_imputed)

for n_folds in cv_configs:
    print(f"\n--- {n_folds}-Fold Cross-Validation ---")
    
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=42)
    
    for model_name, model_template in models.items():
        model = model_template.__class__(**model_template.get_params())
        
        # Run cross-validation
        scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='f1_macro')
        
        # Also get accuracy
        acc_scores = cross_val_score(model, X_scaled, y, cv=cv, scoring='accuracy')
        
        cv_results.append({
            'CV_Folds': f"{n_folds}-fold",
            'Model': model_name,
            'F1_Mean': np.mean(scores),
            'F1_Std': np.std(scores),
            'Accuracy': np.mean(acc_scores)
        })
        
        print(f"  {model_name}: F1={np.mean(scores):.3f} ± {np.std(scores):.3f}, Acc={np.mean(acc_scores):.3f}")

# Print formatted tables for manuscript
print("\n" + "=" * 70)
print("FORMATTED OUTPUT FOR MANUSCRIPT")
print("=" * 70)

print("\n--- Split Ratio Table ---")
print("Split | Model | F1 Score | Accuracy | Stability (std)")
for r in split_results:
    print(f"[{r['Split']}], [{r['Model']}], [{r['F1_Mean']:.3f} ± {r['F1_Std']:.3f}], [{r['Accuracy']:.3f}], [{r['F1_Std']:.3f}],")

print("\n--- Cross-Validation Table ---")
print("CV Folds | Model | F1 Score | Accuracy | Stability (std)")
for r in cv_results:
    print(f"[{r['CV_Folds']}], [{r['Model']}], [{r['F1_Mean']:.3f} ± {r['F1_Std']:.3f}], [{r['Accuracy']:.3f}], [{r['F1_Std']:.3f}],")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
