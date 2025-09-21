from pathlib import Path
import torch
import numpy as np
import pandas as pd

def load_activity_data(filepath='src/data/raw/seo_activity_data.csv'):
    """Load activity data from CSV"""
    df = pd.read_csv(filepath)
    # Set 'Accession code' as index for easy lookup
    df_indexed = df.set_index('Accession code')
    return df_indexed

def load_esm_embeddings_with_activity(embeddings_dir='input_data/seo_act_embeddings/'):
    """Load embeddings and match with activity data"""
    # Load activity data
    activity_df = load_activity_data()
    
    embeddings_dir = Path(embeddings_dir)
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
    
    if matched_count == 0:
        raise ValueError("No sequences matched between embeddings and activity data")
    
    embeddings = np.array(embeddings)
    activities = np.array(activities)
    
    return embeddings, activities, sequence_ids

def inspect_full_dataset():
    """Inspect the complete matched dataset"""
    embeddings, activities, seq_ids = load_esm_embeddings_with_activity()
    
    print(f"=== Complete Dataset ===")
    print(f"Number of sequences: {len(embeddings)}")
    print(f"Embedding shape: {embeddings.shape}")
    print(f"Activity statistics:")
    print(f"  Min: {activities.min()}")
    print(f"  Max: {activities.max()}")
    print(f"  Mean: {activities.mean():.4f}")
    print(f"  Std: {activities.std():.4f}")
    print(f"  Unique values: {len(np.unique(activities))}")
    print(f"  Value distribution:")
    unique, counts = np.unique(activities, return_counts=True)
    for val, count in zip(unique, counts):
        print(f"    {val}: {count} sequences")
    
    return embeddings, activities, seq_ids