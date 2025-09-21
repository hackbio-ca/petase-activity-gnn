import torch
import torch.nn as nn

# Deeper and more Gradual Reduction in Layer Sizes
class ImprovedESMPredictor(nn.Module):
    def __init__(self, input_dim=1280, dropout=0.4):
        super().__init__()
        self.network = nn.Sequential(
            nn.Linear(input_dim, 1024),
            nn.ReLU(),
            nn.BatchNorm1d(1024),
            nn.Dropout(dropout),
            
            nn.Linear(1024, 768),
            nn.ReLU(),
            nn.BatchNorm1d(768),
            nn.Dropout(dropout),
            
            nn.Linear(768, 512),
            nn.ReLU(),
            nn.BatchNorm1d(512),
            nn.Dropout(dropout),
            
            nn.Linear(512, 256),
            nn.ReLU(),
            nn.BatchNorm1d(256),
            nn.Dropout(dropout * 0.7),  # Reduce dropout in later layers
            
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.BatchNorm1d(128),
            nn.Dropout(dropout * 0.5),
            
            nn.Linear(128, 1)
        )

    def forward(self, x):
        return self.network(x)

if __name__ == "__main__":
    print("Testing ImprovedESMPredictor...")
    model = ImprovedESMPredictor()
    dummy_input = torch.randn(4, 1280)
    output = model(dummy_input)
    print("Output shape:", output.shape)