import torch.nn as nn

class EnhancedPetaseModel(nn.Module):
    def __init__(self, input_dim=1280):
        super().__init__()
        
        # Feature extractor (based on actual model dimensions)
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, 1024),     # 1280 -> 1024
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(1024, 1024),          # 1024 -> 1024  
            nn.BatchNorm1d(1024),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(1024, 512),           # 1024 -> 512
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(0.3)
        )
        
        # Classifier branch
        self.classifier = nn.Sequential(
            nn.Linear(512, 128),            # 512 -> 128
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),             # 128 -> 64
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)                # 64 -> 1
        )
        
        # Regressor branch 
        self.regressor = nn.Sequential(
            nn.Linear(512, 256),            # 512 -> 256
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 128),            # 256 -> 128
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),             # 128 -> 64
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(64, 1)                # 64 -> 1 (this is layer 9, not 8)
        )
    
    def forward(self, x, task='regression'):
        features = self.feature_extractor(x)
        
        if task == 'classification':
            return self.classifier(features)
        else:  # regression (default for activity prediction)
            return self.regressor(features)