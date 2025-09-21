import sys
import os

sys.path.append('../../')

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
from src.data.data_loader import inspect_embeddings, load_esm_embeddings
import pandas as pd

def main():
    # First, inspect your embedding structure
    print("=== Inspecting Embedding Structure ===")
    sample_data = inspect_embeddings('../../input_data/hack_seqs-esm_embeddings/')
    
    # Load all embeddings
    print("\n=== Loading All Embeddings ===")
    embeddings, sequence_ids = load_esm_embeddings('../../input_data/hack_seqs-esm_embeddings/')
    print(f"Loaded {len(embeddings)} sequences")
    print(f"Embedding dimension: {embeddings.shape[1]}")
    print(f"Sample sequence IDs: {sequence_ids[:5]}")
    
    # If you have activity data, load it
    # activity_data = pd.read_csv('../../data/raw/activity_data.csv', index_col=0)
    # embeddings, activities, matched_ids = load_esm_embeddings(
    #     '../../input_data/hack_seqs-esm_embeddings/', 
    #     activity_data
    # )
    # print(f"Matched {len(matched_ids)} sequences with activity data")

if __name__ == "__main__":
    main()