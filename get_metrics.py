import pandas as pd
import sys
from pathlib import Path
import os
import builtins

sys.path.append(os.getcwd())

from src.supervised.supervised_learning import run_supervised_pipeline

def main():
    p = Path('data/processed/clustered_dataset.csv')
    df = pd.read_csv(p)
    cols = [c for c in df.columns if c.endswith('_encoded')]
    
    # Quiet run
    original_print = builtins.print
    builtins.print = lambda *args, **kwargs: None
    
    try:
        res = run_supervised_pipeline(df, cols, 'CLUSTER', random_state=42, task_type='cluster')
    finally:
        builtins.print = original_print

    m = res['model_results']['Random Forest']
    
    # Print line by line
    print(f"ACC={m['accuracy']:.4f}")
    print(f"F1M={m['f1_score_macro']:.4f}")
    print(f"PREW={m.get('precision_weighted', 'N/A')}")
    print(f"RECW={m.get('recall_weighted', 'N/A')}")

if __name__ == '__main__':
    main()
