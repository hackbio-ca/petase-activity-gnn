import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import seaborn as sns

class TwoStageESMPredictor(nn.Module):
    """
    Two-stage model: 
    1. Classifier to identify active vs inactive PETases
    2. Regressor to predict log-activity for active ones
    """
    def __init__(self, input_dim, hidden_dim=512, dropout=0.3):
        super().__init__()
        
        # Shared feature extractor
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Classification head (active/inactive)
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim // 2, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1)  # Binary classification
        )
        
        # Regression head (log-activity for active sequences)
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim // 2, 64),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(64, 1)  # Log-activity prediction
        )
    
    def forward(self, x, stage='both'):
        features = self.feature_extractor(x)
        
        if stage == 'classify':
            return torch.sigmoid(self.classifier(features))
        elif stage == 'regress':
            return self.regressor(features)
        else:  # both
            classification = torch.sigmoid(self.classifier(features))
            regression = self.regressor(features)
            return classification, regression

class PETaseDataset(Dataset):
    """Dataset with log-transformed activities and binary labels"""
    def __init__(self, embeddings, activities, activity_threshold=0.01):
        self.embeddings = torch.FloatTensor(embeddings)
        
        # Create binary labels (active/inactive)
        self.is_active = (activities > activity_threshold).astype(int)
        
        # Log-transform activities (add small constant to avoid log(0))
        log_activities = np.log(activities + 1e-8)
        
        # Only keep log-activities for active sequences
        self.log_activities = np.where(self.is_active, log_activities, 0)
        
        self.activities = torch.FloatTensor(activities)
        self.is_active = torch.LongTensor(self.is_active)
        self.log_activities = torch.FloatTensor(self.log_activities)
        
    def __len__(self):
        return len(self.embeddings)
    
    def __getitem__(self, idx):
        return (self.embeddings[idx], 
                self.is_active[idx], 
                self.log_activities[idx],
                self.activities[idx])

class TwoStageTrainer:
    def __init__(self, model, device='cpu', activity_threshold=0.01):
        self.model = model.to(device)
        self.device = device
        self.activity_threshold = activity_threshold
        
        # Separate optimizers for each stage
        self.classifier_optimizer = torch.optim.Adam(
            list(model.feature_extractor.parameters()) + list(model.classifier.parameters()),
            lr=0.001
        )
        self.regressor_optimizer = torch.optim.Adam(
            list(model.feature_extractor.parameters()) + list(model.regressor.parameters()),
            lr=0.001
        )
        
        # Loss functions
        self.classification_loss = nn.BCELoss()
        self.regression_loss = nn.MSELoss()
        
        # Training history
        self.history = {
            'classification_loss': [],
            'regression_loss': [],
            'classification_acc': [],
            'regression_r2': []
        }
    
    def train_classification_epoch(self, dataloader):
        """Train only the classification head"""
        self.model.train()
        total_loss = 0
        correct = 0
        total = 0
        
        for embeddings, is_active, _, _ in dataloader:
            embeddings = embeddings.to(self.device)
            is_active = is_active.to(self.device).float()
            
            self.classifier_optimizer.zero_grad()
            
            # Forward pass (classification only)
            pred_active = self.model(embeddings, stage='classify').squeeze()
            
            # Classification loss
            loss = self.classification_loss(pred_active, is_active)
            loss.backward()
            self.classifier_optimizer.step()
            
            total_loss += loss.item()
            
            # Calculate accuracy
            predicted = (pred_active > 0.5).float()
            correct += (predicted == is_active).sum().item()
            total += is_active.size(0)
        
        accuracy = correct / total
        avg_loss = total_loss / len(dataloader)
        
        return avg_loss, accuracy
    
    def train_regression_epoch(self, dataloader):
        """Train only the regression head on active sequences"""
        self.model.train()
        total_loss = 0
        batch_count = 0
        y_true_all = []
        y_pred_all = []
        
        for embeddings, is_active, log_activities, _ in dataloader:
            embeddings = embeddings.to(self.device)
            is_active = is_active.to(self.device)
            log_activities = log_activities.to(self.device)
            
            # Only train on active sequences
            active_mask = is_active.bool()
            if active_mask.sum() == 0:
                continue
                
            active_embeddings = embeddings[active_mask]
            active_log_activities = log_activities[active_mask]
            
            self.regressor_optimizer.zero_grad()
            
            # Forward pass (regression only)
            pred_log_activity = self.model(active_embeddings, stage='regress').squeeze()
            
            # Regression loss
            loss = self.regression_loss(pred_log_activity, active_log_activities)
            loss.backward()
            self.regressor_optimizer.step()
            
            total_loss += loss.item()
            batch_count += 1
            
            # Store for R2 calculation
            y_true_all.extend(active_log_activities.cpu().numpy())
            y_pred_all.extend(pred_log_activity.detach().cpu().numpy())
        
        avg_loss = total_loss / batch_count if batch_count > 0 else 0
        r2 = r2_score(y_true_all, y_pred_all) if len(y_true_all) > 0 else 0
        
        return avg_loss, r2
    
    def evaluate(self, dataloader):
        """Evaluate both classification and regression"""
        self.model.eval()
        
        # Classification metrics
        class_correct = 0
        class_total = 0
        all_class_true = []
        all_class_pred = []
        
        # Regression metrics (only for active sequences)
        reg_true = []
        reg_pred = []
        
        with torch.no_grad():
            for embeddings, is_active, log_activities, original_activities in dataloader:
                embeddings = embeddings.to(self.device)
                is_active = is_active.to(self.device)
                log_activities = log_activities.to(self.device)
                
                # Get predictions
                pred_active, pred_log_activity = self.model(embeddings, stage='both')
                pred_active = pred_active.squeeze()
                pred_log_activity = pred_log_activity.squeeze()
                
                # Classification evaluation
                predicted_active = (pred_active > 0.5).float()
                class_correct += (predicted_active == is_active.float()).sum().item()
                class_total += is_active.size(0)
                
                all_class_true.extend(is_active.cpu().numpy())
                all_class_pred.extend(predicted_active.cpu().numpy())
                
                # Regression evaluation (only for truly active sequences)
                active_mask = is_active.bool()
                if active_mask.sum() > 0:
                    reg_true.extend(log_activities[active_mask].cpu().numpy())
                    reg_pred.extend(pred_log_activity[active_mask].cpu().numpy())
        
        # Calculate metrics
        classification_acc = class_correct / class_total
        regression_r2 = r2_score(reg_true, reg_pred) if len(reg_true) > 0 else 0
        regression_mse = mean_squared_error(reg_true, reg_pred) if len(reg_true) > 0 else 0
        
        return {
            'classification_accuracy': classification_acc,
            'classification_true': all_class_true,
            'classification_pred': all_class_pred,
            'regression_r2': regression_r2,
            'regression_mse': regression_mse,
            'regression_true': reg_true,
            'regression_pred': reg_pred
        }
    
    def train(self, train_loader, val_loader, epochs=100, classification_epochs=50):
        """
        Two-stage training:
        1. Train classification for first N epochs
        2. Alternate between classification and regression
        """
        print("Stage 1: Training classifier...")
        
        # Stage 1: Focus on classification
        for epoch in range(classification_epochs):
            class_loss, class_acc = self.train_classification_epoch(train_loader)
            
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch + 1}: Classification Loss = {class_loss:.4f}, Accuracy = {class_acc:.4f}")
        
        print("\nStage 2: Joint training...")
        
        # Stage 2: Joint training
        for epoch in range(classification_epochs, epochs):
            # Train both stages
            class_loss, class_acc = self.train_classification_epoch(train_loader)
            reg_loss, reg_r2 = self.train_regression_epoch(train_loader)
            
            # Update history
            self.history['classification_loss'].append(class_loss)
            self.history['regression_loss'].append(reg_loss)
            self.history['classification_acc'].append(class_acc)
            self.history['regression_r2'].append(reg_r2)
            
            if (epoch + 1) % 10 == 0:
                print(f"Epoch {epoch + 1}:")
                print(f"  Classification: Loss = {class_loss:.4f}, Accuracy = {class_acc:.4f}")
                print(f"  Regression: Loss = {reg_loss:.4f}, R² = {reg_r2:.4f}")
        
        # Final evaluation
        print("\nFinal evaluation...")
        results = self.evaluate(val_loader)
        print(f"Validation Classification Accuracy: {results['classification_accuracy']:.4f}")
        print(f"Validation Regression R²: {results['regression_r2']:.4f}")
        
        return results

def plot_results(results, save_path=None):
    """Plot classification and regression results"""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Classification confusion matrix
    cm = confusion_matrix(results['classification_true'], results['classification_pred'])
    sns.heatmap(cm, annot=True, fmt='d', ax=axes[0,0], cmap='Blues')
    axes[0,0].set_title('Classification Confusion Matrix')
    axes[0,0].set_xlabel('Predicted')
    axes[0,0].set_ylabel('True')
    
    # Classification ROC-like plot (predicted probabilities would be better)
    axes[0,1].hist([results['classification_true'], results['classification_pred']], 
                   label=['True', 'Predicted'], alpha=0.7)
    axes[0,1].set_title('Classification Distribution')
    axes[0,1].legend()
    
    # Regression scatter plot
    if len(results['regression_true']) > 0:
        axes[1,0].scatter(results['regression_true'], results['regression_pred'], alpha=0.6)
        axes[1,0].plot([min(results['regression_true']), max(results['regression_true'])], 
                       [min(results['regression_true']), max(results['regression_true'])], 'r--')
        axes[1,0].set_xlabel('True Log-Activity')
        axes[1,0].set_ylabel('Predicted Log-Activity')
        axes[1,0].set_title(f'Regression Results (R² = {results["regression_r2"]:.3f})')
    
    # Regression residuals
    if len(results['regression_true']) > 0:
        residuals = np.array(results['regression_pred']) - np.array(results['regression_true'])
        axes[1,1].scatter(results['regression_pred'], residuals, alpha=0.6)
        axes[1,1].axhline(y=0, color='r', linestyle='--')
        axes[1,1].set_xlabel('Predicted Log-Activity')
        axes[1,1].set_ylabel('Residuals')
        axes[1,1].set_title('Regression Residuals')
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
    plt.show()

# Example usage
if __name__ == "__main__":
    # Load your data (replace with actual data loading)
    # embeddings, activities = load_your_data()
    
    # Example with dummy data
    embeddings = np.random.randn(1000, 1280)  # ESM embeddings
    activities = np.concatenate([
        np.zeros(500),  # Inactive PETases
        np.random.exponential(0.1, 500)  # Active PETases with exponential distribution
    ])
    
    # Create dataset
    dataset = PETaseDataset(embeddings, activities, activity_threshold=0.01)
    
    # Split data
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    
    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False)
    
    # Initialize model and trainer
    model = TwoStageESMPredictor(input_dim=1280)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    trainer = TwoStageTrainer(model, device=device)
    
    # Train model
    results = trainer.train(train_loader, val_loader, epochs=100, classification_epochs=30)
    
    # Plot results
    plot_results(results)