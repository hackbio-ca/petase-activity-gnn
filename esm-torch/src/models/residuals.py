import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import numpy as np

class ResidualESMPredictor(nn.Module):
    """
    A residual neural network for protein embedding prediction using ESM embeddings.
    
    Args:
        input_dim (int): Dimension of input features (default: 1280 for ESM-2 embeddings)
    
    The network consists of:
    - An input projection layer
    - Three residual blocks with decreasing dimensions
    - A final output layer predicting a single value
    """
    def __init__(self, input_dim=1280):
        super().__init__()
        self.input_proj = nn.Linear(input_dim, 512)
        
        self.block1 = self._make_block(512, 512)
        self.block2 = self._make_block(512, 256)
        self.block3 = self._make_block(256, 128)
        
        self.output = nn.Linear(128, 1)
    
    def _make_block(self, in_dim, out_dim):
        """Creates a residual block with the specified input and output dimensions."""
        return nn.Sequential(
            nn.Linear(in_dim, out_dim),
            nn.ReLU(),
            nn.BatchNorm1d(out_dim),
            nn.Dropout(0.3),
            nn.Linear(out_dim, out_dim),
            nn.ReLU(),
            nn.BatchNorm1d(out_dim)
        )
    
    def forward(self, x):
        """
        Forward pass of the model.
        
        Args:
            x (torch.Tensor): Input tensor of shape (batch_size, input_dim)
            
        Returns:
            torch.Tensor: Predicted values of shape (batch_size, 1)
        """
        x = self.input_proj(x)
        x = torch.relu(x)  # Add activation after projection
        
        # Residual blocks with proper skip connections
        residual = x
        x = self.block1(x)
        if x.shape == residual.shape:  # Only add residual if dimensions match
            x = x + residual
        
        x = self.block2(x)
        x = self.block3(x)
        return self.output(x)

class HybridResidualESMPredictor(nn.Module):
    """
    Hybrid model that can handle both individual predictions and pairwise comparisons.
    """
    def __init__(self, input_dim=1280):
        super().__init__()
        self.base_model = ResidualESMPredictor(input_dim)
        
    def forward(self, x1, x2=None):
        """
        Forward pass that handles both individual and pairwise inputs.
        
        Args:
            x1 (torch.Tensor): First input (or only input for individual prediction)
            x2 (torch.Tensor, optional): Second input for pairwise comparison
            
        Returns:
            torch.Tensor: Individual prediction if x2 is None, 
                         pairwise comparison probability if x2 is provided
        """
        if x2 is None:
            # Individual prediction
            return self.base_model(x1)
        else:
            # Pairwise comparison
            pred1 = self.base_model(x1)
            pred2 = self.base_model(x2)
            # Return the difference (will be passed through sigmoid in loss function)
            return pred1 - pred2

def plot_predictions(model, test_loader, save_path=None):
    """
    Plot predicted vs true values and calculate R².
    
    Args:
        model (nn.Module): Trained model
        test_loader (DataLoader): Test data loader
        save_path (str, optional): Path to save the plot
    
    Returns:
        float: R² score
    """
    model.eval()
    y_true = []
    y_pred = []
    
    with torch.no_grad():
        for x, labels in test_loader:
            predictions = model(x)
            y_true.extend(labels.numpy())
            y_pred.extend(predictions.squeeze().numpy())
    
    # Calculate R²
    r2 = r2_score(y_true, y_pred)
    
    # Create plot
    plt.figure(figsize=(8, 8))
    plt.scatter(y_true, y_pred, alpha=0.5)
    
    # Add perfect prediction line
    min_val = min(min(y_true), min(y_pred))
    max_val = max(max(y_true), max(y_pred))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='Perfect prediction')
    
    plt.xlabel('True Values')
    plt.ylabel('Predicted Values')
    plt.title(f'Predicted vs True Values (R² = {r2:.3f})')
    plt.legend()
    
    if save_path:
        plt.savefig(save_path)
    plt.show()
    
    return r2

if __name__ == "__main__":
    print("HybridResidualESMPredictor module loaded successfully")