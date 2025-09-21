import sys
sys.path.append('.')

import torch
import esm
import numpy as np
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from torch.utils.data import DataLoader
import torch.optim as optim
import torch.nn as nn

from src.data.data_loader import load_esm_embeddings_with_activity
from src.data.dataset import ESMActivityDataset
from src.models.baseline import ESMActivityPredictor

def extract_small_esm_embeddings(embeddings_dir='input_data/seo_act_embeddings/'):
    """Extract embeddings using the smallest ESM model"""
    print("Loading ESM2 6M model...")
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    
    # Load activity data to get sequence list
    from src.data.data_loader import load_activity_data
    activity_df, activity_col = load_activity_data()
    
    embeddings = []
    activities = []
    sequence_ids = []
    
    # Read sequences from your embedding files to get the actual sequences
    # This is a bit hacky - ideally you'd have the original FASTA
    embeddings_dir = Path(embeddings_dir)
    
    for pt_file in embeddings_dir.glob('*.pt'):
        data = torch.load(pt_file, map_location='cpu')
        seq_id = data['label']
        
        if seq_id in activity_df.index:
            # You'll need the actual sequence here
            # For now, we'll skip this and assume you have a FASTA file
            print(f"Would process {seq_id} if sequence available")
    
    print("Note: You need to provide the original FASTA file to re-extract with small ESM")
    print("Using existing embeddings for now...")
    
    # Fall back to loading existing embeddings (but with warning)
    return load_esm_embeddings_with_activity()

def train_small_esm_model():
    """Train model with small ESM embeddings"""
    # Load data (you'll need to modify this to use small ESM)
    embeddings, activities, sequence_ids = extract_small_esm_embeddings()
    
    print(f"Training with embedding dimension: {embeddings.shape[1]}")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        embeddings, activities, test_size=0.2, random_state=42
    )
    
    # Create datasets
    train_dataset = ESMActivityDataset(X_train, y_train)
    test_dataset = ESMActivityDataset(X_test, y_test)
    
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    
    # Model with correct input dimension for small ESM
    input_dim = embeddings.shape[1]  # Should be 320 for 6M model
    model = ESMActivityPredictor(
        input_dim=input_dim,
        hidden_dims=[256, 128, 64],  # Smaller network for smaller embeddings
        dropout=0.3
    )
    
    optimizer = optim.Adam(model.parameters(), lr=0.001)
    criterion = nn.MSELoss()
    
    print("Training small ESM model...")
    for epoch in range(100):
        model.train()
        for batch_x, batch_y in train_loader:
            optimizer.zero_grad()
            predictions = model(batch_x)
            loss = criterion(predictions.squeeze(), batch_y)
            loss.backward()
            optimizer.step()
        
        if epoch % 20 == 0:
            model.eval()
            test_preds = []
            test_targets = []
            with torch.no_grad():
                for batch_x, batch_y in test_loader:
                    pred = model(batch_x)
                    test_preds.extend(pred.squeeze().numpy())
                    test_targets.extend(batch_y.numpy())
            
            r2 = r2_score(test_targets, test_preds)
            print(f"Epoch {epoch}: Test R² = {r2:.4f}")
    
    # Save model
    results_dir = Path('experiments/small_esm/results/')
    results_dir.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), results_dir / 'small_esm_model.pth')
    
    print(f"Model saved to {results_dir / 'small_esm_model.pth'}")
    return model

if __name__ == "__main__":
    train_small_esm_model()