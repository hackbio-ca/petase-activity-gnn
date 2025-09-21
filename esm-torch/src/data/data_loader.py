from pathlib import Path
import torch
import numpy as np
import pandas as pd

def load_activity_data(filepath='src/data/raw/seo_activity_data.csv'):
    """Load activity data from CSV"""
    df = pd.read_csv(filepath)
    
    # Try to identify the sequence ID and activity columns
    id_col='Accession code'
    activity_col='Activity'
    
    for col in df.columns:
        if any(keyword in col.lower() for keyword in ['id', 'sequence', 'protein', 'name']):
            id_col = col
        if any(keyword in col.lower() for keyword in ['activity', 'value', 'measurement', 'target']):
            activity_col = col
    
    if id_col is None or activity_col is None:
        print("Available columns:", list(df.columns))
        raise ValueError("Could not automatically identify ID and activity columns")
    
    # Set sequence ID as index
    df_indexed = df.set_index(id_col)
    return df_indexed, activity_col

def load_esm_embeddings_with_activity(embeddings_dir='input_data/seo_act_embeddings/'):
    """Load embeddings and match with activity data"""
    # Load activity data
    activity_df, activity_col = load_activity_data()
    
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
                if isinstance(data['mean_representations'], dict):
                    embedding = data['mean_representations'][33].numpy()
                else:
                    embedding = data['mean_representations'].numpy()
                
                embeddings.append(embedding)
                activities.append(activity_df.loc[seq_id, activity_col])
                sequence_ids.append(seq_id)
                matched_count += 1
                
        except Exception as e:
            print(f"Error loading {pt_file.name}: {e}")
            continue
    
    print(f"Matched {matched_count}/{total_count} sequences with activity data")
    
    if matched_count == 0:
        raise ValueError("No sequences matched between embeddings and activity data")
    
    return np.array(embeddings), np.array(activities), sequence_ids

def inspect_data_match():
    """Quick function to check data matching"""
    try:
        embeddings, activities, seq_ids = load_esm_embeddings_with_activity()
        print(f"Successfully loaded {len(embeddings)} matched sequences")
        print(f"Embedding shape: {embeddings.shape}")
        print(f"Activity range: {activities.min():.3f} to {activities.max():.3f}")
        print(f"Sample sequence IDs: {seq_ids[:5]}")
        return True
    except Exception as e:
        print(f"Data loading failed: {e}")
        return False