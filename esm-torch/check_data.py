from pathlib import Path
import torch
import numpy as np
import pandas as pd

def load_and_inspect_data():
    """Load and inspect the complete dataset"""
    # Load activity data
    activity_df = pd.read_csv('src/data/raw/seo_activity_data.csv')
    activity_df = activity_df.set_index('Accession code')
    
    embeddings_dir = Path('input_data/seo_act_embeddings/')
    embeddings = []
    activities = []
    sequence_ids = []
    
    matched_count = 0
    total_count = 0
    
    for pt_file in embeddings_dir.glob('*.pt'):
        total_count += 1
        try:
            data = torch.load(pt_file, map_location='cpu')
            seq_id = data['label']
            
            # Check if this sequence has activity data
            if seq_id in activity_df.index:
                # Extract layer 33 mean embedding
                embedding = data['mean_representations'][33].numpy()
                
                embeddings.append(embedding)
                activities.append(activity_df.loc[seq_id, 'Activity'])
                sequence_ids.append(seq_id)
                matched_count += 1
                
        except Exception as e:
            print(f"Error loading {pt_file.name}: {e}")
            continue
    
    print(f"Matched {matched_count}/{total_count} sequences with activity data")
    
    embeddings = np.array(embeddings)
    activities = np.array(activities)
    
    print(f"\n=== Dataset Summary ===")
    print(f"Number of sequences: {len(embeddings)}")
    print(f"Embedding shape: {embeddings.shape}")
    print(f"Activity statistics:")
    print(f"  Min: {activities.min()}")
    print(f"  Max: {activities.max()}")
    print(f"  Mean: {activities.mean():.4f}")
    print(f"  Std: {activities.std():.4f}")
    print(f"  Unique values: {len(np.unique(activities))}")
    print(f"  Value counts:")
    unique, counts = np.unique(activities, return_counts=True)
    for val, count in zip(unique, counts):
        print(f"    {val}: {count} sequences")
    
    # Show some sample data
    print(f"\nSample sequence IDs: {sequence_ids[:5]}")
    print(f"Sample activities: {activities[:10]}")
    
    return embeddings, activities, sequence_ids

if __name__ == "__main__":
    embeddings, activities, seq_ids = load_and_inspect_data()