"""
Manuscript Verification Script
Compares claimed values in manuscript against actual code outputs.
"""
import pandas as pd
import pickle
from pathlib import Path

print("="*80)
print("MANUSCRIPT VS CODEBASE VERIFICATION REPORT")
print("="*80)

# Load actual data
df = pd.read_csv('data/processed/clustered_dataset.csv')

discrepancies = []
matches = []

# ============================================================================
# 1. TOTAL ISOLATES
# ============================================================================
print("\n1. TOTAL ISOLATES")
actual_n = len(df)
claimed_n = 491
if actual_n == claimed_n:
    matches.append(f"Total isolates: {actual_n} ✓")
    print(f"   MATCH: Claimed={claimed_n}, Actual={actual_n} ✓")
else:
    discrepancies.append(f"Total isolates: Claimed={claimed_n}, Actual={actual_n}")
    print(f"   DISCREPANCY: Claimed={claimed_n}, Actual={actual_n} ✗")

# ============================================================================
# 2. CLUSTER SIZES
# ============================================================================
print("\n2. CLUSTER SIZES")
claimed_sizes = {1: 23, 2: 93, 3: 123, 4: 252}
actual_sizes = df['CLUSTER'].value_counts().sort_index().to_dict()
for c in [1, 2, 3, 4]:
    claimed = claimed_sizes.get(c, 0)
    actual = actual_sizes.get(c, 0)
    if claimed == actual:
        matches.append(f"C{c} size: {actual} ✓")
        print(f"   C{c}: Claimed={claimed}, Actual={actual} ✓")
    else:
        discrepancies.append(f"C{c} size: Claimed={claimed}, Actual={actual}")
        print(f"   C{c}: DISCREPANCY Claimed={claimed}, Actual={actual} ✗")

# ============================================================================
# 3. MDR RATES PER CLUSTER
# ============================================================================
print("\n3. MDR RATES PER CLUSTER")
claimed_mdr = {1: 4.3, 2: 2.2, 3: 53.7, 4: 0.4}
for c in [1, 2, 3, 4]:
    cluster_df = df[df['CLUSTER'] == c]
    mdr_count = cluster_df['MDR_FLAG'].sum()
    actual_rate = (mdr_count / len(cluster_df)) * 100
    claimed_rate = claimed_mdr.get(c, 0)
    diff = abs(actual_rate - claimed_rate)
    if diff < 0.2:  # Allow 0.2% tolerance for rounding
        matches.append(f"C{c} MDR rate: {actual_rate:.1f}% ✓")
        print(f"   C{c}: Claimed={claimed_rate}%, Actual={actual_rate:.1f}% ({mdr_count}/{len(cluster_df)}) ✓")
    else:
        discrepancies.append(f"C{c} MDR rate: Claimed={claimed_rate}%, Actual={actual_rate:.1f}%")
        print(f"   C{c}: DISCREPANCY Claimed={claimed_rate}%, Actual={actual_rate:.1f}% ✗")

# ============================================================================
# 4. TOTAL MDR COUNT
# ============================================================================
print("\n4. TOTAL MDR COUNT")
actual_mdr = int(df['MDR_FLAG'].sum())
claimed_mdr_total = 70  # From manuscript: "all 70 MDR isolates"
if actual_mdr == claimed_mdr_total:
    matches.append(f"Total MDR: {actual_mdr} ✓")
    print(f"   MATCH: Claimed={claimed_mdr_total}, Actual={actual_mdr} ✓")
else:
    discrepancies.append(f"Total MDR: Claimed={claimed_mdr_total}, Actual={actual_mdr}")
    print(f"   DISCREPANCY: Claimed={claimed_mdr_total}, Actual={actual_mdr} ✗")

# ============================================================================
# 5. NUMBER OF ANTIBIOTICS
# ============================================================================
print("\n5. NUMBER OF ANTIBIOTICS")
feature_cols = [c for c in df.columns if c.endswith('_encoded')]
actual_antibiotics = len(feature_cols)
claimed_antibiotics = 22
if actual_antibiotics == claimed_antibiotics:
    matches.append(f"Antibiotics: {actual_antibiotics} ✓")
    print(f"   MATCH: Claimed={claimed_antibiotics}, Actual={actual_antibiotics} ✓")
else:
    discrepancies.append(f"Antibiotics: Claimed={claimed_antibiotics}, Actual={actual_antibiotics}")
    print(f"   DISCREPANCY: Claimed={claimed_antibiotics}, Actual={actual_antibiotics} ✗")

# ============================================================================
# 6. SPECIES DISTRIBUTION (C1 = 100% Salmonella)
# ============================================================================
print("\n6. SPECIES IN CLUSTER 1 (Claimed: 100% Salmonella)")
c1_df = df[df['CLUSTER'] == 1]
c1_species = c1_df['ISOLATE_ID'].value_counts()
salmonella_pct = 0
for sp in c1_species.index:
    if 'salmonella' in sp.lower():
        salmonella_pct = (c1_species[sp] / len(c1_df)) * 100
        break
if salmonella_pct >= 99:  # Allow small tolerance
    matches.append(f"C1 Salmonella: {salmonella_pct:.1f}% ✓")
    print(f"   MATCH: Claimed=100%, Actual={salmonella_pct:.1f}% ✓")
else:
    discrepancies.append(f"C1 Salmonella: Claimed=100%, Actual={salmonella_pct:.1f}%")
    print(f"   DISCREPANCY: Claimed=100%, Actual={salmonella_pct:.1f}% ✗")
print(f"   Actual species: {c1_species.to_dict()}")

# ============================================================================
# 7. C2 SPECIES (Claimed: Enterobacter cloacae 71.0%)
# ============================================================================
print("\n7. SPECIES IN CLUSTER 2 (Claimed: E. cloacae 71.0%)")
c2_df = df[df['CLUSTER'] == 2]
c2_species = c2_df['ISOLATE_ID'].value_counts()
ec_count = 0
for sp in c2_species.index:
    if 'enterobacter cloacae' in sp.lower():
        ec_count = c2_species[sp]
        break
actual_ec_pct = (ec_count / len(c2_df)) * 100
claimed_ec_pct = 71.0
if abs(actual_ec_pct - claimed_ec_pct) < 1:
    matches.append(f"C2 E. cloacae: {actual_ec_pct:.1f}% ✓")
    print(f"   MATCH: Claimed={claimed_ec_pct}%, Actual={actual_ec_pct:.1f}% ✓")
else:
    discrepancies.append(f"C2 E. cloacae: Claimed={claimed_ec_pct}%, Actual={actual_ec_pct:.1f}%")
    print(f"   DISCREPANCY: Claimed={claimed_ec_pct}%, Actual={actual_ec_pct:.1f}% ✗")
print(f"   Actual species: {c2_species.head(3).to_dict()}")

# ============================================================================
# 8. C3 SPECIES (Claimed: E. coli 77.2%, K. pneumoniae 22.0%)
# ============================================================================
print("\n8. SPECIES IN CLUSTER 3 (Claimed: E. coli 77.2%, K. pneumoniae 22.0%)")
c3_df = df[df['CLUSTER'] == 3]
c3_species = c3_df['ISOLATE_ID'].value_counts()
ecoli_count = kp_count = 0
for sp in c3_species.index:
    if 'escherichia coli' in sp.lower():
        ecoli_count = c3_species[sp]
    if 'klebsiella pneumoniae' in sp.lower():
        kp_count = c3_species[sp]
actual_ecoli_pct = (ecoli_count / len(c3_df)) * 100
actual_kp_pct = (kp_count / len(c3_df)) * 100
print(f"   E. coli: Claimed=77.2%, Actual={actual_ecoli_pct:.1f}%")
print(f"   K. pneumoniae: Claimed=22.0%, Actual={actual_kp_pct:.1f}%")
if abs(actual_ecoli_pct - 77.2) < 1 and abs(actual_kp_pct - 22.0) < 1:
    matches.append(f"C3 E. coli: {actual_ecoli_pct:.1f}% ✓")
    matches.append(f"C3 K. pneumoniae: {actual_kp_pct:.1f}% ✓")
else:
    if abs(actual_ecoli_pct - 77.2) >= 1:
        discrepancies.append(f"C3 E. coli: Claimed=77.2%, Actual={actual_ecoli_pct:.1f}%")
    if abs(actual_kp_pct - 22.0) >= 1:
        discrepancies.append(f"C3 K. pneumoniae: Claimed=22.0%, Actual={actual_kp_pct:.1f}%")

# ============================================================================
# 9. C4 SPECIES (Claimed: E. coli 51.2%, K. pneumoniae 47.2%)
# ============================================================================
print("\n9. SPECIES IN CLUSTER 4 (Claimed: E. coli 51.2%, K. pneumoniae 47.2%)")
c4_df = df[df['CLUSTER'] == 4]
c4_species = c4_df['ISOLATE_ID'].value_counts()
ecoli_count = kp_count = 0
for sp in c4_species.index:
    if 'escherichia coli' in sp.lower():
        ecoli_count = c4_species[sp]
    if 'klebsiella pneumoniae' in sp.lower():
        kp_count = c4_species[sp]
actual_ecoli_pct = (ecoli_count / len(c4_df)) * 100
actual_kp_pct = (kp_count / len(c4_df)) * 100
print(f"   E. coli: Claimed=51.2%, Actual={actual_ecoli_pct:.1f}%")
print(f"   K. pneumoniae: Claimed=47.2%, Actual={actual_kp_pct:.1f}%")

# ============================================================================
# 10. REGIONAL DISTRIBUTIONS
# ============================================================================
print("\n10. REGIONAL DISTRIBUTION")
region_counts = df['REGION'].value_counts()
print(f"   Regions: {region_counts.to_dict()}")

# C3 regional claim: 53.7% from BARMM
c3_barmm = len(c3_df[c3_df['REGION'].str.contains('BARMM', na=False)])
c3_barmm_pct = (c3_barmm / len(c3_df)) * 100
print(f"   C3 from BARMM: Claimed=53.7%, Actual={c3_barmm_pct:.1f}%")
if abs(c3_barmm_pct - 53.7) < 1:
    matches.append(f"C3 BARMM: {c3_barmm_pct:.1f}% ✓")
else:
    discrepancies.append(f"C3 BARMM: Claimed=53.7%, Actual={c3_barmm_pct:.1f}%")

# ============================================================================
# 11. CLUSTERING METRICS (Silhouette, etc.)
# ============================================================================
print("\n11. CLUSTERING METRICS")
try:
    with open('data/processed/clustering_artifacts/clustering_info.pkl', 'rb') as f:
        info = pickle.load(f)
    print(f"   n_clusters: {info.get('n_clusters', 'N/A')}")
    print(f"   method: {info.get('method', 'N/A')}")
    print(f"   metric: {info.get('metric', 'N/A')}")
    
    if 'robustness_check' in info and info['robustness_check']:
        ari = info['robustness_check'].get('agreement_score', {}).get('adjusted_rand_index', 'N/A')
        print(f"   Robustness ARI: {ari}")
except Exception as e:
    print(f"   Could not load clustering_info: {e}")

# Load validation metrics
try:
    val_df = pd.read_csv('data/processed/figures/cluster_validation_metrics.csv')
    k4_row = val_df[val_df['k'] == 4]
    if not k4_row.empty:
        actual_sil = k4_row['silhouette_score'].values[0]
        claimed_sil = 0.466
        print(f"   Silhouette (k=4): Claimed={claimed_sil}, Actual={actual_sil:.3f}")
        if abs(actual_sil - claimed_sil) < 0.001:
            matches.append(f"Silhouette k=4: {actual_sil:.3f} ✓")
        else:
            discrepancies.append(f"Silhouette k=4: Claimed={claimed_sil}, Actual={actual_sil:.3f}")
except Exception as e:
    print(f"   Could not load validation metrics: {e}")

# ============================================================================
# 12. SUPERVISED LEARNING METRICS
# ============================================================================
print("\n12. SUPERVISED LEARNING METRICS")
try:
    with open('data/processed/models/cluster_metrics.txt', 'r') as f:
        metrics_str = f.read()
    print(f"   (Raw metrics file found)")
    # Check if RF accuracy ~ 99%
    if 'Random Forest' in metrics_str:
        # Parse if possible
        pass
except Exception as e:
    print(f"   Could not load cluster_metrics: {e}")

# Try to load the importance file
try:
    import json
    with open('data/processed/models/cluster_antibiotic_importance.json', 'r') as f:
        importance = json.load(f)
    print(f"   Feature importance file found with {len(importance)} entries")
    # Check top feature
    if importance:
        sorted_imp = sorted(importance.items(), key=lambda x: x[1], reverse=True)
        print(f"   Top feature: {sorted_imp[0][0]} = {sorted_imp[0][1]:.3f}")
        # Claimed: Tetracycline = 0.241
        if 'TE' in sorted_imp[0][0] or 'Tetracycline' in sorted_imp[0][0]:
            if abs(sorted_imp[0][1] - 0.241) < 0.01:
                matches.append(f"Top feature TE: {sorted_imp[0][1]:.3f} ✓")
            else:
                print(f"   Note: TE importance {sorted_imp[0][1]:.3f} vs claimed 0.241")
except Exception as e:
    print(f"   Could not load importance: {e}")

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "="*80)
print("VERIFICATION SUMMARY")
print("="*80)
print(f"\n✓ MATCHES ({len(matches)}):")
for m in matches:
    print(f"   {m}")

print(f"\n✗ DISCREPANCIES ({len(discrepancies)}):")
if discrepancies:
    for d in discrepancies:
        print(f"   {d}")
else:
    print("   None found!")

print(f"\n{'='*80}")
print(f"FINAL STATUS: {len(matches)} matches, {len(discrepancies)} discrepancies")
print("="*80)
