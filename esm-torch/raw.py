# Extract ESM embeddings for your activity dataset
# Use layer 33 for optimal functional prediction (already done)
# esm-extract esm2_t33_650M_UR50D activity_sequences.fasta embeddings_dir \
#   --repr_layers 33 --include mean

# Load embeddings and activity labels
import torch
import numpy as np
import pandas as pd

def load_esm_embeddings(embeddings_dir, activity_data):
    embeddings = []
    activities = []
    sequence_ids = []
    
    for pt_file in Path(embeddings_dir).glob('*.pt'):
        data = torch.load(pt_file)
        seq_id = data['label']
        if seq_id in activity_data.index:
            embeddings.append(data['mean_representations'][33].numpy())
            activities.append(activity_data.loc[seq_id, 'activity'])
            sequence_ids.append(seq_id)
    
    return np.array(embeddings), np.array(activities), sequence_ids