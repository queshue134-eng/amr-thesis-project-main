
import pickle
import numpy as np
from scipy.cluster.hierarchy import fcluster
import sys

# Load linkage matrix
try:
    with open(r'c:\Users\quesh\Downloads\amr-thesis-project-main-v3.3\amr_thesis_project_code\data\processed\clustering_artifacts\linkage_matrix.pkl', 'rb') as f:
        Z = pickle.load(f)
    print("Linkage Matrix Loaded Successfully.")
    
    # The distance column is index 2
    distances = Z[:, 2]
    
    # Sort distances to understand the merge levels
    sorted_distances = np.sort(distances)
    
    # For k clusters, the cut is made between the (k-1)th last merge and the kth last merge
    # The dendrogram height (distance) at which the merge happens
    
    # Last merge (k=2 -> k=1)
    # 2nd last merge (k=3 -> k=2)
    # 3rd last merge (k=4 -> k=3)
    # 4th last merge (k=5 -> k=4)
    
    # To get exactly k=4, the threshold t must be:
    # distance(5->4) < t <= distance(4->3)
    
    merge_dist_4_to_3 = sorted_distances[-3] # Distance of the merge that creates 2 clusters from 3? No.
    # The last element [-1] is the distance of merging 2 clusters into 1.
    # [-2] is merging 3 into 2.
    # [-3] is merging 4 into 3.
    # [-4] is merging 5 into 4.
    
    print(f"\n--- Analysis for k=4 ---")
    print(f"Distance for merging into 3 clusters (dist[-3]): {sorted_distances[-3]:.4f}")
    print(f"Distance for merging into 4 clusters (dist[-4]): {sorted_distances[-4]:.4f}")
    
    # The threshold usually cited is just above the forming merge, or the midpoint.
    # Often standard implementations use values just below the next merge.
    
    # Let's verify with fcluster
    t = sorted_distances[-3]
    labels_at_t = fcluster(Z, t, criterion='distance')
    print(f"At t={t:.4f}, n_clusters={len(np.unique(labels_at_t))}")
    
    # Usually we want the smallest distance that still yields k=4, or the "Height" of the cut.
    # The 'height' of the cut is often described as the distance of the merge that forms the system.
    
    print(f"\nComputed Euclidean Value (Cut Height for k=4): {sorted_distances[-3]:.4f}")
    
except Exception as e:
    print(f"Error: {e}")
