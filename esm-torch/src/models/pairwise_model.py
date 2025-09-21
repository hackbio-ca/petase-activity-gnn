import torch
import torch.nn as nn
from .baseline import ESMActivityPredictor

class PairwiseESMPredictor(nn.Module):
    def __init__(self, input_dim=1280, hidden_dims=[512, 256, 128], dropout=0.3):
        super().__init__()
        # Shared network for both proteins
        self.shared_network = ESMActivityPredictor(input_dim, hidden_dims, dropout)
    
    def forward(self, x1, x2=None):
        if x2 is None:
            # Single protein prediction
            return self.shared_network(x1)
        else:
            # Pairwise prediction
            pred1 = self.shared_network(x1)
            pred2 = self.shared_network(x2)
            return pred1 - pred2
    
    def predict_single(self, x):
        """Convenience method for single protein prediction"""
        return self.shared_network(x)

class HybridLoss(nn.Module):
    def __init__(self, alpha=0.7):
        super().__init__()
        self.alpha = alpha  # Weight for absolute loss
        self.mse = nn.MSELoss()
    
    def forward(self, pred_abs, true_abs, pred_diff, true_diff):
        abs_loss = self.mse(pred_abs, true_abs)
        diff_loss = self.mse(pred_diff, true_diff)
        total_loss = self.alpha * abs_loss + (1 - self.alpha) * diff_loss
        return total_loss, abs_loss.item(), diff_loss.item()