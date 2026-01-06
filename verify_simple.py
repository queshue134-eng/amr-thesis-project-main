import pandas as pd
import pickle

df = pd.read_csv('data/processed/clustered_dataset.csv')
print("VERIFICATION RESULTS")
print("="*60)

# 1. Total isolates
n = len(df)
print("1. Total: Claimed=491, Actual=" + str(n) + " [" + ("MATCH" if n==491 else "MISMATCH") + "]")

# 2. Cluster sizes
sizes = df['CLUSTER'].value_counts().sort_index().to_dict()
claimed = {1: 23, 2: 93, 3: 123, 4: 252}
for c in [1,2,3,4]:
    actual = sizes.get(c, 0)
    status = "MATCH" if actual == claimed[c] else "MISMATCH"
    print("   C" + str(c) + ": Claimed=" + str(claimed[c]) + ", Actual=" + str(actual) + " [" + status + "]")

# 3. MDR rates per cluster
print("3. MDR RATES:")
mdr_claimed = {1: 4.3, 2: 2.2, 3: 53.7, 4: 0.4}
for c in [1,2,3,4]:
    cdf = df[df['CLUSTER']==c]
    mdr = int(cdf['MDR_FLAG'].sum())
    rate = mdr/len(cdf)*100
    diff = abs(rate - mdr_claimed[c])
    status = "MATCH" if diff < 0.5 else "MISMATCH"
    print("   C" + str(c) + ": " + str(mdr) + "/" + str(len(cdf)) + " = " + str(round(rate,1)) + "% (Claimed=" + str(mdr_claimed[c]) + "%) [" + status + "]")

# 4. Total MDR
total_mdr = int(df['MDR_FLAG'].sum())
print("4. Total MDR: " + str(total_mdr) + " (Claimed=70) [" + ("MATCH" if total_mdr==70 else "MISMATCH") + "]")

# 5. Antibiotics
n_ab = len([c for c in df.columns if c.endswith('_encoded')])
print("5. Antibiotics: " + str(n_ab) + " (Claimed=22) [" + ("MATCH" if n_ab==22 else "MISMATCH") + "]")

# 6. Species distribution
print("6. Species:")
for sp, n in df['ISOLATE_ID'].value_counts().items():
    pct = round(n/len(df)*100, 1)
    print("   " + sp + ": " + str(n) + " (" + str(pct) + "%)")

# 7. Regions
print("7. Regions:")
for r, n in df['REGION'].value_counts().items():
    pct = round(n/len(df)*100, 1)
    print("   " + r + ": " + str(n) + " (" + str(pct) + "%)")

# 8. Silhouette
try:
    val = pd.read_csv('data/processed/figures/cluster_validation_metrics.csv')
    k4 = val[val['k']==4]['silhouette_score'].values[0]
    diff = abs(k4 - 0.466)
    status = "MATCH" if diff < 0.001 else "MISMATCH"
    print("8. Silhouette k=4: " + str(round(k4,4)) + " (Claimed=0.466) [" + status + "]")
except Exception as e:
    print("8. Could not load silhouette: " + str(e))

# 9. C1 species (100% Salmonella)
c1 = df[df['CLUSTER']==1]
c1_species = c1['ISOLATE_ID'].value_counts()
print("9. C1 Species (Claimed: 100% Salmonella):")
for sp, n in c1_species.items():
    print("   " + sp + ": " + str(n))

# 10. C2 species (E. cloacae 71%)
c2 = df[df['CLUSTER']==2]
c2_species = c2['ISOLATE_ID'].value_counts()
print("10. C2 Species (Claimed: E. cloacae 71.0%):")
for sp, n in c2_species.items():
    pct = round(n/len(c2)*100, 1)
    print("    " + sp + ": " + str(n) + " (" + str(pct) + "%)")

# 11. C3 species (E. coli 77.2%, K. pneumoniae 22.0%)
c3 = df[df['CLUSTER']==3]
c3_species = c3['ISOLATE_ID'].value_counts()
print("11. C3 Species (Claimed: E.coli 77.2%, K.pneumoniae 22.0%):")
for sp, n in c3_species.items():
    pct = round(n/len(c3)*100, 1)
    print("    " + sp + ": " + str(n) + " (" + str(pct) + "%)")

# 12. C4 species (E. coli 51.2%, K. pneumoniae 47.2%)
c4 = df[df['CLUSTER']==4]
c4_species = c4['ISOLATE_ID'].value_counts()
print("12. C4 Species (Claimed: E.coli 51.2%, K.pneumoniae 47.2%):")
for sp, n in c4_species.items():
    pct = round(n/len(c4)*100, 1)
    print("    " + sp + ": " + str(n) + " (" + str(pct) + "%)")

print("="*60)
print("VERIFICATION COMPLETE")
