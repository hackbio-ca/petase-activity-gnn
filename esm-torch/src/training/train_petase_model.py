#!/usr/bin/env python3

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset, DataLoader
import seaborn as sns
import os
import sys
import json

# Add the current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

class TwoStageESMPredictor(nn.Module):
    """Two-stage model: Classifier + Regressor for PETase activity"""
    def __init__(self, input_dim, hidden_dim=128, dropout=0.5):  # Much smaller, more dropout
        super().__init__()
        
        # Simpler shared feature extractor
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Simpler classification head
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim // 2, 1)
        )
        
        # Simpler regression head
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim // 2, 1)
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
            lr=0.0005, weight_decay=0.01  # Lower learning rate, add regularization
        )
        self.regressor_optimizer = torch.optim.Adam(
            list(model.feature_extractor.parameters()) + list(model.regressor.parameters()),
            lr=0.0005, weight_decay=0.01
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
            
            # Handle single sample case
            if pred_active.dim() == 0:
                pred_active = pred_active.unsqueeze(0)
            
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
            
            # Handle single sample case
            if pred_log_activity.dim() == 0:
                pred_log_activity = pred_log_activity.unsqueeze(0)
            if active_log_activities.dim() == 0:
                active_log_activities = active_log_activities.unsqueeze(0)
            
            # Regression loss
            loss = self.regression_loss(pred_log_activity, active_log_activities)
            loss.backward()
            self.regressor_optimizer.step()
            
            total_loss += loss.item()
            batch_count += 1
            
            # Store for R2 calculation - ensure we can iterate
            y_true_all.extend(active_log_activities.cpu().numpy().flatten())
            y_pred_all.extend(pred_log_activity.detach().cpu().numpy().flatten())
        
        avg_loss = total_loss / batch_count if batch_count > 0 else 0
        r2 = r2_score(y_true_all, y_pred_all) if len(y_true_all) > 1 else 0
        
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
                
                # Handle single sample case
                if pred_active.dim() == 0:
                    pred_active = pred_active.unsqueeze(0)
                if pred_log_activity.dim() == 0:
                    pred_log_activity = pred_log_activity.unsqueeze(0)
                
                # Classification evaluation
                predicted_active = (pred_active > 0.5).float()
                class_correct += (predicted_active == is_active.float()).sum().item()
                class_total += is_active.size(0)
                
                all_class_true.extend(is_active.cpu().numpy().flatten())
                all_class_pred.extend(predicted_active.cpu().numpy().flatten())
                
                # Regression evaluation (only for truly active sequences)
                active_mask = is_active.bool()
                if active_mask.sum() > 0:
                    reg_true.extend(log_activities[active_mask].cpu().numpy().flatten())
                    reg_pred.extend(pred_log_activity[active_mask].cpu().numpy().flatten())
        
        # Calculate metrics
        classification_acc = class_correct / class_total
        regression_r2 = r2_score(reg_true, reg_pred) if len(reg_true) > 1 else 0
        regression_mse = mean_squared_error(reg_true, reg_pred) if len(reg_true) > 1 else 0
        
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
        """Two-stage training"""
        print("Training classifier...")
        
        # Stage 1: Focus on classification
        for epoch in range(classification_epochs):
            class_loss, class_acc = self.train_classification_epoch(train_loader)
            
            if (epoch + 1) % 20 == 0:
                print(f"Epoch {epoch + 1}: Classification Acc = {class_acc:.3f}")
        
        print("Joint training...")
        
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
            
            if (epoch + 1) % 20 == 0:
                print(f"Epoch {epoch + 1}: Class Acc = {class_acc:.3f}, Reg R² = {reg_r2:.3f}")
        
        # Final evaluation
        results = self.evaluate(val_loader)
        return results

def load_data(embeddings_folder, activities_path):
    """Load embeddings and activities data"""
    # Load activities first to get the accession codes
    activities_df = pd.read_csv(activities_path)
    accession_codes = activities_df['Accession code'].values
    activities = activities_df['Activity'].values
    
    # Load embeddings from .pt files
    embeddings_list = []
    missing_embeddings = []
    
    for accession in accession_codes:
        embedding_file = os.path.join(embeddings_folder, f"{accession}.pt")
        
        if os.path.exists(embedding_file):
            try:
                embedding = torch.load(embedding_file, map_location='cpu')
                # Handle different possible formats
                if isinstance(embedding, dict):
                    # If it's a dictionary, try common keys
                    if 'representations' in embedding:
                        embedding = embedding['representations'][33]  # Last layer for ESM
                    elif 'embedding' in embedding:
                        embedding = embedding['embedding']
                    elif 'mean_representations' in embedding:
                        embedding = embedding['mean_representations'][33]
                    else:
                        # Take the first tensor value if it's a dict
                        embedding = list(embedding.values())[0]
                
                # Ensure it's a 1D tensor and convert to numpy
                if embedding.dim() > 1:
                    embedding = embedding.mean(dim=0)  # Average over sequence length
                
                embeddings_list.append(embedding.numpy())
                
            except Exception as e:
                missing_embeddings.append(accession)
        else:
            missing_embeddings.append(accession)
    
    # Filter to only sequences with embeddings
    valid_indices = [i for i, acc in enumerate(accession_codes) if acc not in missing_embeddings]
    
    if len(valid_indices) == 0:
        raise ValueError("No valid embeddings found!")
    
    embeddings = np.array([embeddings_list[i] for i in range(len(embeddings_list))])
    activities = activities[valid_indices]
    valid_accessions = accession_codes[valid_indices]
    
    return embeddings, activities, valid_accessions

def plot_results(results, history, save_dir='results'):
    """Plot and save results"""
    os.makedirs(save_dir, exist_ok=True)
    
    # Main results figure
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Training curves
    if history['classification_loss']:
        epochs = range(1, len(history['classification_loss']) + 1)
        
        axes[0,0].plot(epochs, history['classification_loss'], 'b-')
        axes[0,0].set_title('Classification Loss')
        axes[0,0].set_xlabel('Epoch')
        axes[0,0].grid(True)
        
        axes[0,1].plot(epochs, history['classification_acc'], 'g-')
        axes[0,1].set_title('Classification Accuracy')
        axes[0,1].set_xlabel('Epoch')
        axes[0,1].grid(True)
        
        axes[0,2].plot(epochs, history['regression_r2'], 'purple')
        axes[0,2].set_title('Regression R²')
        axes[0,2].set_xlabel('Epoch')
        axes[0,2].grid(True)
    
    # Confusion matrix
    cm = confusion_matrix(results['classification_true'], results['classification_pred'])
    sns.heatmap(cm, annot=True, fmt='d', ax=axes[1,0], cmap='Blues',
                xticklabels=['Inactive', 'Active'], yticklabels=['Inactive', 'Active'])
    axes[1,0].set_title('Classification Confusion Matrix')
    
    # Regression results
    if len(results['regression_true']) > 0:
        axes[1,1].scatter(results['regression_true'], results['regression_pred'], alpha=0.6)
        min_val = min(min(results['regression_true']), min(results['regression_pred']))
        max_val = max(max(results['regression_true']), max(results['regression_pred']))
        axes[1,1].plot([min_val, max_val], [min_val, max_val], 'r--')
        axes[1,1].set_xlabel('True Log-Activity')
        axes[1,1].set_ylabel('Predicted Log-Activity')
        axes[1,1].set_title(f'Regression (R² = {results["regression_r2"]:.3f})')
        axes[1,1].grid(True, alpha=0.3)
        
        # Residuals
        residuals = np.array(results['regression_pred']) - np.array(results['regression_true'])
        axes[1,2].scatter(results['regression_pred'], residuals, alpha=0.6)
        axes[1,2].axhline(y=0, color='r', linestyle='--')
        axes[1,2].set_xlabel('Predicted Log-Activity')
        axes[1,2].set_ylabel('Residuals')
        axes[1,2].set_title('Residuals')
        axes[1,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'training_results.png'), dpi=300, bbox_inches='tight')
    plt.close()

def save_results_data(results, history, accession_codes, save_dir='results'):
    """Save results data"""
    os.makedirs(save_dir, exist_ok=True)
    
    # Ensure arrays have same length
    min_length = min(len(accession_codes), len(results['classification_true']), len(results['classification_pred']))
    
    # Save predictions
    predictions_df = pd.DataFrame({
        'accession_code': accession_codes[:min_length],
        'true_active': results['classification_true'][:min_length],
        'pred_active': results['classification_pred'][:min_length]
    })
    
    # Add regression results for active sequences
    if len(results['regression_true']) > 0:
        active_indices = [i for i, active in enumerate(results['classification_true']) if active == 1]
        
        predictions_df['true_log_activity'] = np.nan
        predictions_df['pred_log_activity'] = np.nan
        
        for i, (true_log, pred_log) in enumerate(zip(results['regression_true'], results['regression_pred'])):
            if i < len(active_indices):
                idx = active_indices[i]
                predictions_df.loc[idx, 'true_log_activity'] = true_log
                predictions_df.loc[idx, 'pred_log_activity'] = pred_log
    
    predictions_df.to_csv(os.path.join(save_dir, 'predictions.csv'), index=False)
    
    # Save training history
    if history['classification_loss']:
        history_df = pd.DataFrame(history)
        history_df.to_csv(os.path.join(save_dir, 'training_history.csv'), index=False)
    
    # Save summary metrics
    summary_metrics = {
        'classification_accuracy': float(results['classification_accuracy']),
        'regression_r2': float(results['regression_r2']),
        'regression_mse': float(results['regression_mse']),
        'total_sequences': len(results['classification_true']),
        'true_active_count': int(sum(results['classification_true'])),
        'pred_active_count': int(sum(results['classification_pred'])),
        'regression_samples': len(results['regression_true'])
    }
    
    with open(os.path.join(save_dir, 'summary_metrics.json'), 'w') as f:
        json.dump(summary_metrics, f, indent=2)

def main():
    # Configuration - more conservative settings
    EMBEDDINGS_FOLDER = 'input_data/seo_act_embeddings'
    ACTIVITIES_PATH = 'src/data/raw/seo_activity_data.csv'
    ACTIVITY_THRESHOLD = 0.01
    BATCH_SIZE = 16  # Smaller batch size
    EPOCHS = 50      # Fewer epochs to prevent overfitting
    CLASSIFICATION_EPOCHS = 20
    
    print("Loading data...")
    try:
        embeddings, activities, accession_codes = load_data(EMBEDDINGS_FOLDER, ACTIVITIES_PATH)
        print(f"Loaded {len(embeddings)} sequences with {embeddings.shape[1]} features")
        print(f"Active sequences: {sum(activities > ACTIVITY_THRESHOLD)}")
        print(f"Inactive sequences: {sum(activities <= ACTIVITY_THRESHOLD)}")
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Create dataset
    dataset = PETaseDataset(embeddings, activities, activity_threshold=ACTIVITY_THRESHOLD)
    
    # Split data
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    
    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)
    
    # Initialize model and trainer
    input_dim = embeddings.shape[1]
    model = TwoStageESMPredictor(input_dim=input_dim)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")
    
    trainer = TwoStageTrainer(model, device=device, activity_threshold=ACTIVITY_THRESHOLD)
    
    # Train model
    results = trainer.train(train_loader, val_loader, 
                           epochs=EPOCHS, classification_epochs=CLASSIFICATION_EPOCHS)
    
    # Debug: Check array lengths
    print(f"Debug - Array lengths:")
    print(f"  accession_codes: {len(accession_codes)}")
    print(f"  classification_true: {len(results['classification_true'])}")
    print(f"  classification_pred: {len(results['classification_pred'])}")
    
    # Create results directory and save everything
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    
    # Save model
    torch.save(model.state_dict(), os.path.join(results_dir, 'petase_model.pth'))
    
    # Save results
    plot_results(results, trainer.history, save_dir=results_dir)
    save_results_data(results, trainer.history, accession_codes, save_dir=results_dir)
    
    # Print final results
    print(f"\nFinal Results:")
    print(f"Classification Accuracy: {results['classification_accuracy']:.3f}")
    print(f"Regression R²: {results['regression_r2']:.3f}")
    print(f"Results saved to: {os.path.abspath(results_dir)}")

if __name__ == "__main__":
    main()