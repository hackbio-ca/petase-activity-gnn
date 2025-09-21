import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import json
import os

class HybridLoss(nn.Module):
    def __init__(self, alpha=0.5):
        super().__init__()
        self.alpha = alpha
        self.mse_loss = nn.MSELoss()
        self.bce_loss = nn.BCEWithLogitsLoss()
    
    def forward(self, pred_abs, true_abs, pred_diff, true_diff):
        abs_loss = self.mse_loss(pred_abs, true_abs)
        diff_loss = self.bce_loss(pred_diff, true_diff)
        total_loss = self.alpha * abs_loss + (1 - self.alpha) * diff_loss
        return total_loss, abs_loss.item(), diff_loss.item()

class PairwiseTrainer:
    def __init__(self, model, optimizer, device='cpu'):
        self.model = model.to(device)
        self.optimizer = optimizer
        self.device = device
        self.loss_fn = HybridLoss()
        self.history = {'total_loss': [], 'abs_loss': [], 'diff_loss': []}
    
    def train_epoch(self, individual_loader, pairwise_loader):
        self.model.train()
        total_loss = total_abs = total_diff = batch_count = 0
        
        for (ind_x, ind_y), (pair_x1, pair_x2, pair_y) in zip(individual_loader, pairwise_loader):
            # Move to device
            ind_x, ind_y = ind_x.to(self.device), ind_y.to(self.device)
            pair_x1, pair_x2, pair_y = pair_x1.to(self.device), pair_x2.to(self.device), pair_y.to(self.device)
            
            self.optimizer.zero_grad()
            
            # Forward pass
            pred_abs = self.model(ind_x)
            pred_diff = self.model(pair_x1, pair_x2)
            
            # Calculate loss
            loss, abs_loss, diff_loss = self.loss_fn(
                pred_abs.squeeze(), ind_y,
                pred_diff.squeeze(), pair_y.squeeze()
            )
            
            loss.backward()
            self.optimizer.step()
            
            total_loss += loss.item()
            total_abs += abs_loss
            total_diff += diff_loss
            batch_count += 1
        
        # Update history
        avg_total = total_loss / batch_count
        avg_abs = total_abs / batch_count
        avg_diff = total_diff / batch_count
        
        self.history['total_loss'].append(avg_total)
        self.history['abs_loss'].append(avg_abs)
        self.history['diff_loss'].append(avg_diff)
        
        return {'total_loss': avg_total, 'abs_loss': avg_abs, 'diff_loss': avg_diff}
    
    def evaluate(self, test_loader, individual=True):
        self.model.eval()
        y_true, y_pred = [], []
        
        with torch.no_grad():
            for batch in test_loader:
                if individual:
                    x, y = batch
                    x, y = x.to(self.device), y.to(self.device)
                    pred = self.model(x)
                else:  # pairwise
                    x1, x2, y = batch
                    x1, x2, y = x1.to(self.device), x2.to(self.device), y.to(self.device)
                    pred = torch.sigmoid(self.model(x1, x2))
                
                y_true.extend(y.cpu().numpy())
                y_pred.extend(pred.squeeze().cpu().numpy())
        
        return {
            'mse': mean_squared_error(y_true, y_pred),
            'mae': mean_absolute_error(y_true, y_pred),
            'r2': r2_score(y_true, y_pred),
            'predictions': (y_true, y_pred)
        }
    
    def save_model(self, path):
        torch.save(self.model.state_dict(), path)
    
    def plot_training(self, save_path=None):
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        epochs = range(1, len(self.history['total_loss']) + 1)
        
        axes[0].plot(epochs, self.history['total_loss'])
        axes[0].set_title('Total Loss')
        axes[0].grid(True)
        
        axes[1].plot(epochs, self.history['abs_loss'])
        axes[1].set_title('Absolute Loss')
        axes[1].grid(True)
        
        axes[2].plot(epochs, self.history['diff_loss'])
        axes[2].set_title('Difference Loss')
        axes[2].grid(True)
        
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
        plt.show()

def train_model(model, train_loaders, test_loaders, epochs=50, lr=0.001, device='cpu'):
    """Main training function"""
    individual_train, pairwise_train = train_loaders
    individual_test, pairwise_test = test_loaders
    
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    trainer = PairwiseTrainer(model, optimizer, device)
    
    print(f"Training for {epochs} epochs on {device}")
    
    for epoch in range(epochs):
        # Train
        metrics = trainer.train_epoch(individual_train, pairwise_train)
        
        # Print progress
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch + 1}: Loss={metrics['total_loss']:.4f}")
    
    # Evaluate
    individual_results = trainer.evaluate(individual_test, individual=True)
    pairwise_results = trainer.evaluate(pairwise_test, individual=False)
    
    print(f"\nResults:")
    print(f"Individual R²: {individual_results['r2']:.4f}")
    print(f"Pairwise R²: {pairwise_results['r2']:.4f}")
    
    return trainer, individual_results, pairwise_results

# Example usage - replace with your actual imports and data
if __name__ == "__main__":
    from src.models.residuals import HybridResidualESMPredictor
    from src.data.data_loader import load_esm_embeddings_with_activity
    from src.data.dataset import ESMActivityDataset, PairwiseESMDataset
    from torch.utils.data import DataLoader
    from sklearn.model_selection import train_test_split
    
    # Load and prepare data
    print("Loading data...")
    embeddings, activities, sequence_ids = load_esm_embeddings_with_activity()
    
    # Split data
    X_train, X_test, y_train, y_test, _, _ = train_test_split(
        embeddings, activities, sequence_ids, test_size=0.2, random_state=42
    )
    
    # Create datasets and loaders
    individual_train = DataLoader(ESMActivityDataset(X_train, y_train), batch_size=32, shuffle=True)
    individual_test = DataLoader(ESMActivityDataset(X_test, y_test), batch_size=32)
    pairwise_train = DataLoader(PairwiseESMDataset(X_train, y_train), batch_size=32, shuffle=True)
    pairwise_test = DataLoader(PairwiseESMDataset(X_test, y_test), batch_size=32)
    
    # Initialize model
    model = HybridResidualESMPredictor(input_dim=embeddings.shape[1])
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    # Train
    trainer, individual_results, pairwise_results = train_model(
        model=model,
        train_loaders=(individual_train, pairwise_train),
        test_loaders=(individual_test, pairwise_test),
        epochs=50,
        device=device
    )
    
    # Plot results
    trainer.plot_training()
    
    # Save model if needed
    trainer.save_model('model.pth')