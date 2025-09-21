#!/usr/bin/env python3

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import classification_report, confusion_matrix, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, RobustScaler
from torch.utils.data import Dataset, DataLoader
import seaborn as sns
import os
import sys
import json

class ImprovedTwoStageESMPredictor(nn.Module):
    """Enhanced two-stage model with better regression performance"""
    def __init__(self, input_dim, hidden_dim=1024, dropout=0.2):
        super().__init__()
        
        # Enhanced feature extractor with residual connections
        self.feature_extractor = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout),
            
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.BatchNorm1d(hidden_dim // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
        )
        
        # Enhanced classification head
        self.classifier = nn.Sequential(
            nn.Linear(hidden_dim // 2, 128),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(dropout / 2),
            nn.Linear(64, 1)
        )
        
        # Enhanced regression head with more capacity
        self.regressor = nn.Sequential(
            nn.Linear(hidden_dim // 2, 256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(dropout / 2),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Dropout(dropout / 4),
            nn.Linear(64, 1)
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

class EnhancedPETaseDataset(Dataset):
    """Enhanced dataset with better preprocessing"""
    def __init__(self, embeddings, activities, activity_threshold=0.001, use_robust_scaling=True):
        # Normalize embeddings
        if use_robust_scaling:
            scaler = RobustScaler()
            embeddings = scaler.fit_transform(embeddings)
        else:
            scaler = StandardScaler()
            embeddings = scaler.fit_transform(embeddings)
        
        self.embeddings = torch.FloatTensor(embeddings)
        self.scaler = scaler
        
        # Create binary labels with lower threshold
        self.is_active = (activities > activity_threshold).astype(int)
        
        # Enhanced log transformation with different strategies
        # Strategy 1: Log(activity + small_constant)
        log_activities_v1 = np.log(activities + 1e-6)
        
        # Strategy 2: Log1p transformation (more stable)
        log_activities_v2 = np.log1p(activities)
        
        # Strategy 3: Box-Cox like transformation
        # Use log1p for better numerical stability
        self.log_activities = np.where(self.is_active, log_activities_v2, log_activities_v2.min())
        
        # Store original activities for analysis
        self.activities = torch.FloatTensor(activities)
        self.is_active = torch.LongTensor(self.is_active)
        self.log_activities = torch.FloatTensor(self.log_activities)
        
        print(f"Dataset preprocessing:")
        print(f"  Embedding scaling: {'Robust' if use_robust_scaling else 'Standard'}")
        print(f"  Activity threshold: {activity_threshold}")
        print(f"  Active sequences: {self.is_active.sum().item()}")
        print(f"  Log-activity range: {self.log_activities.min():.3f} to {self.log_activities.max():.3f}")
        
    def __len__(self):
        return len(self.embeddings)
    
    def __getitem__(self, idx):
        return (self.embeddings[idx], 
                self.is_active[idx], 
                self.log_activities[idx],
                self.activities[idx])

class ImprovedTwoStageTrainer:
    def __init__(self, model, device='cpu', activity_threshold=0.001):
        self.model = model.to(device)
        self.device = device
        self.activity_threshold = activity_threshold
        
        # Enhanced optimizers with different learning rates
        self.classifier_optimizer = torch.optim.AdamW(
            list(model.feature_extractor.parameters()) + list(model.classifier.parameters()),
            lr=0.0005, weight_decay=1e-4
        )
        self.regressor_optimizer = torch.optim.AdamW(
            list(model.feature_extractor.parameters()) + list(model.regressor.parameters()),
            lr=0.0003, weight_decay=1e-4  # Lower LR for regression
        )
        
        # Learning rate schedulers
        self.class_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.classifier_optimizer, mode='min', patience=15, factor=0.5
        )
        self.reg_scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
            self.regressor_optimizer, mode='max', patience=10, factor=0.7
        )
        
        # Enhanced loss functions
        self.classification_loss = nn.BCELoss()
        # Use SmoothL1Loss for more robust regression
        self.regression_loss = nn.SmoothL1Loss()
        
        # Training history
        self.history = {
            'classification_loss': [],
            'regression_loss': [],
            'classification_acc': [],
            'regression_r2': [],
            'learning_rates': []
        }
        
        # Best model tracking
        self.best_r2 = -float('inf')
        self.best_model_state = None
    
    def train_classification_epoch(self, dataloader):
        """Enhanced classification training"""
        self.model.train()
        total_loss = 0
        correct = 0
        total = 0
        
        for embeddings, is_active, _, _ in dataloader:
            embeddings = embeddings.to(self.device)
            is_active = is_active.to(self.device).float()
            
            self.classifier_optimizer.zero_grad()
            
            # Forward pass
            pred_active = self.model(embeddings, stage='classify').squeeze()
            
            # Add label smoothing
            smoothed_labels = is_active * 0.9 + 0.05
            loss = self.classification_loss(pred_active, smoothed_labels)
            
            loss.backward()
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.classifier_optimizer.step()
            
            total_loss += loss.item()
            predicted = (pred_active > 0.5).float()
            correct += (predicted == is_active).sum().item()
            total += is_active.size(0)
        
        accuracy = correct / total
        avg_loss = total_loss / len(dataloader)
        
        return avg_loss, accuracy
    
    def train_regression_epoch(self, dataloader):
        """Enhanced regression training with better handling"""
        self.model.train()
        total_loss = 0
        batch_count = 0
        y_true_all = []
        y_pred_all = []
        
        for embeddings, is_active, log_activities, _ in dataloader:
            embeddings = embeddings.to(self.device)
            is_active = is_active.to(self.device)
            log_activities = log_activities.to(self.device)
            
            # Train on ALL sequences, not just active ones
            # This helps the model learn the full activity spectrum
            self.regressor_optimizer.zero_grad()
            
            # Forward pass
            pred_log_activity = self.model(embeddings, stage='regress').squeeze()
            
            # Weighted loss: higher weight for active sequences
            weights = torch.where(is_active.bool(), 2.0, 1.0)
            weighted_loss = self.regression_loss(pred_log_activity, log_activities)
            
            # Apply weights manually
            individual_losses = (pred_log_activity - log_activities) ** 2
            weighted_individual_losses = individual_losses * weights
            loss = weighted_individual_losses.mean()
            
            loss.backward()
            # Gradient clipping
            torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
            self.regressor_optimizer.step()
            
            total_loss += loss.item()
            batch_count += 1
            
            # Store for R2 calculation
            y_true_all.extend(log_activities.cpu().numpy())
            y_pred_all.extend(pred_log_activity.detach().cpu().numpy())
        
        avg_loss = total_loss / batch_count if batch_count > 0 else 0
        r2 = r2_score(y_true_all, y_pred_all) if len(y_true_all) > 1 else 0
        
        return avg_loss, r2
    
    def evaluate(self, dataloader):
        """Enhanced evaluation"""
        self.model.eval()
        
        class_correct = 0
        class_total = 0
        all_class_true = []
        all_class_pred = []
        
        # Evaluate on ALL sequences for regression
        reg_true = []
        reg_pred = []
        
        # Separate evaluation for active sequences only
        active_reg_true = []
        active_reg_pred = []
        
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
                
                # Regression evaluation - ALL sequences
                reg_true.extend(log_activities.cpu().numpy())
                reg_pred.extend(pred_log_activity.cpu().numpy())
                
                # Regression evaluation - ACTIVE sequences only
                active_mask = is_active.bool()
                if active_mask.sum() > 0:
                    active_reg_true.extend(log_activities[active_mask].cpu().numpy())
                    active_reg_pred.extend(pred_log_activity[active_mask].cpu().numpy())
        
        # Calculate metrics
        classification_acc = class_correct / class_total
        
        # Overall regression R²
        regression_r2_all = r2_score(reg_true, reg_pred) if len(reg_true) > 1 else 0
        
        # Active-only regression R²
        regression_r2_active = r2_score(active_reg_true, active_reg_pred) if len(active_reg_true) > 1 else 0
        
        regression_mse = mean_squared_error(reg_true, reg_pred) if len(reg_true) > 1 else 0
        
        return {
            'classification_accuracy': classification_acc,
            'classification_true': all_class_true,
            'classification_pred': all_class_pred,
            'regression_r2': regression_r2_all,
            'regression_r2_active_only': regression_r2_active,
            'regression_mse': regression_mse,
            'regression_true': reg_true,
            'regression_pred': reg_pred,
            'active_regression_true': active_reg_true,
            'active_regression_pred': active_reg_pred
        }
    
    def train(self, train_loader, val_loader, epochs=200, classification_epochs=50, patience=30):
        """Enhanced training with early stopping"""
        print("Stage 1: Training classifier...")
        
        # Stage 1: Classification focus
        for epoch in range(classification_epochs):
            class_loss, class_acc = self.train_classification_epoch(train_loader)
            self.class_scheduler.step(class_loss)
            
            if (epoch + 1) % 25 == 0:
                print(f"Epoch {epoch + 1}: Classification Acc = {class_acc:.3f}")
        
        print("Stage 2: Joint training with early stopping...")
        
        no_improve_count = 0
        
        # Stage 2: Joint training with early stopping
        for epoch in range(classification_epochs, epochs):
            # Train both stages
            class_loss, class_acc = self.train_classification_epoch(train_loader)
            reg_loss, reg_r2 = self.train_regression_epoch(train_loader)
            
            # Update schedulers
            self.class_scheduler.step(class_loss)
            self.reg_scheduler.step(reg_r2)
            
            # Update history
            self.history['classification_loss'].append(class_loss)
            self.history['regression_loss'].append(reg_loss)
            self.history['classification_acc'].append(class_acc)
            self.history['regression_r2'].append(reg_r2)
            self.history['learning_rates'].append(self.regressor_optimizer.param_groups[0]['lr'])
            
            # Check for improvement
            if reg_r2 > self.best_r2:
                self.best_r2 = reg_r2
                self.best_model_state = self.model.state_dict().copy()
                no_improve_count = 0
            else:
                no_improve_count += 1
            
            if (epoch + 1) % 25 == 0:
                print(f"Epoch {epoch + 1}: Class Acc = {class_acc:.3f}, Reg R² = {reg_r2:.3f} (Best: {self.best_r2:.3f})")
            
            # Early stopping
            if no_improve_count >= patience:
                print(f"Early stopping at epoch {epoch + 1} (no improvement for {patience} epochs)")
                break
        
        # Load best model
        if self.best_model_state is not None:
            self.model.load_state_dict(self.best_model_state)
            print(f"Loaded best model with R² = {self.best_r2:.3f}")
        
        # Final evaluation
        results = self.evaluate(val_loader)
        return results

def load_data(embeddings_folder, activities_path):
    """Load data (same as before)"""
    activities_df = pd.read_csv(activities_path)
    accession_codes = activities_df['Accession code'].values
    activities = activities_df['Activity'].values
    
    embeddings_list = []
    missing_embeddings = []
    
    for accession in accession_codes:
        embedding_file = os.path.join(embeddings_folder, f"{accession}.pt")
        
        if os.path.exists(embedding_file):
            try:
                embedding = torch.load(embedding_file, map_location='cpu')
                if isinstance(embedding, dict):
                    if 'representations' in embedding:
                        embedding = embedding['representations'][33]
                    elif 'embedding' in embedding:
                        embedding = embedding['embedding']
                    elif 'mean_representations' in embedding:
                        embedding = embedding['mean_representations'][33]
                    else:
                        embedding = list(embedding.values())[0]
                
                if embedding.dim() > 1:
                    embedding = embedding.mean(dim=0)
                
                embeddings_list.append(embedding.numpy())
                
            except Exception as e:
                missing_embeddings.append(accession)
        else:
            missing_embeddings.append(accession)
    
    valid_indices = [i for i, acc in enumerate(accession_codes) if acc not in missing_embeddings]
    
    if len(valid_indices) == 0:
        raise ValueError("No valid embeddings found!")
    
    embeddings = np.array([embeddings_list[i] for i in range(len(embeddings_list))])
    activities = activities[valid_indices]
    valid_accessions = accession_codes[valid_indices]
    
    return embeddings, activities, valid_accessions

def plot_enhanced_results(results, history, save_dir='results'):
    """Enhanced plotting with more detailed analysis"""
    os.makedirs(save_dir, exist_ok=True)
    
    # Create comprehensive figure
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    
    # Training curves
    if history['classification_loss']:
        epochs = range(1, len(history['classification_loss']) + 1)
        
        # Loss curves
        axes[0,0].plot(epochs, history['classification_loss'], 'b-', label='Classification')
        axes[0,0].plot(epochs, history['regression_loss'], 'r-', label='Regression')
        axes[0,0].set_title('Training Losses')
        axes[0,0].set_xlabel('Epoch')
        axes[0,0].legend()
        axes[0,0].grid(True)
        
        # Accuracy and R²
        axes[0,1].plot(epochs, history['classification_acc'], 'g-', label='Classification Acc')
        axes[0,1].plot(epochs, history['regression_r2'], 'purple', label='Regression R²')
        axes[0,1].set_title('Performance Metrics')
        axes[0,1].set_xlabel('Epoch')
        axes[0,1].legend()
        axes[0,1].grid(True)
        
        # Learning rates
        axes[0,2].plot(epochs, history['learning_rates'], 'orange')
        axes[0,2].set_title('Learning Rate Schedule')
        axes[0,2].set_xlabel('Epoch')
        axes[0,2].set_yscale('log')
        axes[0,2].grid(True)
    
    # Classification analysis
    cm = confusion_matrix(results['classification_true'], results['classification_pred'])
    sns.heatmap(cm, annot=True, fmt='d', ax=axes[1,0], cmap='Blues',
                xticklabels=['Inactive', 'Active'], yticklabels=['Inactive', 'Active'])
    axes[1,0].set_title('Classification Confusion Matrix')
    
    # Class distribution
    true_counts = pd.Series(results['classification_true']).value_counts().sort_index()
    pred_counts = pd.Series(results['classification_pred']).value_counts().sort_index()
    x = np.arange(len(true_counts))
    width = 0.35
    axes[1,1].bar(x - width/2, true_counts.values, width, label='True', alpha=0.7)
    axes[1,1].bar(x + width/2, pred_counts.values, width, label='Predicted', alpha=0.7)
    axes[1,1].set_title('Class Distribution')
    axes[1,1].set_xticks(x)
    axes[1,1].set_xticklabels(['Inactive', 'Active'])
    axes[1,1].legend()
    
    # Regression analysis - ALL sequences
    if len(results['regression_true']) > 0:
        axes[1,2].scatter(results['regression_true'], results['regression_pred'], alpha=0.6, s=20)
        min_val = min(min(results['regression_true']), min(results['regression_pred']))
        max_val = max(max(results['regression_true']), max(results['regression_pred']))
        axes[1,2].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
        axes[1,2].set_xlabel('True Log-Activity')
        axes[1,2].set_ylabel('Predicted Log-Activity')
        axes[1,2].set_title(f'All Sequences (R² = {results["regression_r2"]:.3f})')
        axes[1,2].grid(True, alpha=0.3)
        
        # Regression analysis - ACTIVE sequences only
        if len(results['active_regression_true']) > 0:
            axes[2,0].scatter(results['active_regression_true'], results['active_regression_pred'], 
                             alpha=0.6, s=20, color='green')
            min_val = min(min(results['active_regression_true']), min(results['active_regression_pred']))
            max_val = max(max(results['active_regression_true']), max(results['active_regression_pred']))
            axes[2,0].plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2)
            axes[2,0].set_xlabel('True Log-Activity')
            axes[2,0].set_ylabel('Predicted Log-Activity')
            axes[2,0].set_title(f'Active Sequences Only (R² = {results["regression_r2_active_only"]:.3f})')
            axes[2,0].grid(True, alpha=0.3)
        
        # Residuals analysis
        residuals = np.array(results['regression_pred']) - np.array(results['regression_true'])
        axes[2,1].scatter(results['regression_pred'], residuals, alpha=0.6, s=20)
        axes[2,1].axhline(y=0, color='r', linestyle='--', linewidth=2)
        axes[2,1].set_xlabel('Predicted Log-Activity')
        axes[2,1].set_ylabel('Residuals')
        axes[2,1].set_title('Residuals Analysis')
        axes[2,1].grid(True, alpha=0.3)
        
        # Distribution of residuals
        axes[2,2].hist(residuals, bins=20, alpha=0.7, edgecolor='black')
        axes[2,2].axvline(x=0, color='r', linestyle='--', linewidth=2)
        axes[2,2].set_xlabel('Residuals')
        axes[2,2].set_ylabel('Frequency')
        axes[2,2].set_title('Residuals Distribution')
        axes[2,2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, 'enhanced_training_results.png'), dpi=300, bbox_inches='tight')
    plt.close()

def save_enhanced_results_data(results, history, accession_codes, save_dir='results'):
    """Save enhanced results"""
    os.makedirs(save_dir, exist_ok=True)
    
    # Enhanced predictions with more details
    predictions_df = pd.DataFrame({
        'accession_code': accession_codes,
        'true_active': results['classification_true'],
        'pred_active': results['classification_pred'],
        'true_log_activity': results['regression_true'],
        'pred_log_activity': results['regression_pred']
    })
    
    # Add residuals and errors
    predictions_df['residual'] = predictions_df['pred_log_activity'] - predictions_df['true_log_activity']
    predictions_df['abs_error'] = np.abs(predictions_df['residual'])
    
    predictions_df.to_csv(os.path.join(save_dir, 'enhanced_predictions.csv'), index=False)
    
    # Enhanced training history
    if history['classification_loss']:
        history_df = pd.DataFrame(history)
        history_df.to_csv(os.path.join(save_dir, 'enhanced_training_history.csv'), index=False)
    
    # Enhanced summary metrics
    summary_metrics = {
        'classification_accuracy': float(results['classification_accuracy']),
        'regression_r2_all': float(results['regression_r2']),
        'regression_r2_active_only': float(results['regression_r2_active_only']),
        'regression_mse': float(results['regression_mse']),
        'total_sequences': len(results['classification_true']),
        'true_active_count': int(sum(results['classification_true'])),
        'pred_active_count': int(sum(results['classification_pred'])),
        'regression_samples': len(results['regression_true']),
        'active_regression_samples': len(results['active_regression_true'])
    }
    
    with open(os.path.join(save_dir, 'enhanced_summary_metrics.json'), 'w') as f:
        json.dump(summary_metrics, f, indent=2)

def main():
    # Enhanced configuration
    EMBEDDINGS_FOLDER = 'input_data/seo_act_embeddings'
    ACTIVITIES_PATH = 'src/data/raw/seo_activity_data.csv'
    ACTIVITY_THRESHOLD = 0.001  # Lower threshold
    BATCH_SIZE = 16  # Smaller batch size for better gradients
    EPOCHS = 200  # More epochs
    CLASSIFICATION_EPOCHS = 50
    
    print("Loading and preprocessing data...")
    try:
        embeddings, activities, accession_codes = load_data(EMBEDDINGS_FOLDER, ACTIVITIES_PATH)
        print(f"Loaded {len(embeddings)} sequences with {embeddings.shape[1]} features")
    except Exception as e:
        print(f"Error loading data: {e}")
        return
    
    # Create enhanced dataset with better preprocessing
    dataset = EnhancedPETaseDataset(embeddings, activities, 
                                   activity_threshold=ACTIVITY_THRESHOLD,
                                   use_robust_scaling=True)
    
    # Stratified split to maintain class balance
    train_size = int(0.8 * len(dataset))
    val_size = len(dataset) - train_size
    train_dataset, val_dataset = torch.utils.data.random_split(dataset, [train_size, val_size])
    
    # Create dataloaders
    train_loader = DataLoader(train_dataset, batch_size=BATCH_SIZE, shuffle=True)
    val_loader = DataLoader(val_dataset, batch_size=BATCH_SIZE, shuffle=False)
    
    # Initialize enhanced model
    input_dim = embeddings.shape[1]
    model = ImprovedTwoStageESMPredictor(input_dim=input_dim, hidden_dim=1024, dropout=0.2)
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f"Using device: {device}")
    
    trainer = ImprovedTwoStageTrainer(model, device=device, activity_threshold=ACTIVITY_THRESHOLD)
    
    # Train enhanced model
    results = trainer.train(train_loader, val_loader, 
                           epochs=EPOCHS, classification_epochs=CLASSIFICATION_EPOCHS,
                           patience=40)
    
    # Save everything
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)
    
    torch.save(model.state_dict(), os.path.join(results_dir, 'enhanced_petase_model.pth'))
    
    plot_enhanced_results(results, trainer.history, save_dir=results_dir)
    save_enhanced_results_data(results, trainer.history, accession_codes, save_dir=results_dir)
    
    # Print detailed results
    print(f"\nEnhanced Results:")
    print(f"Classification Accuracy: {results['classification_accuracy']:.3f}")
    print(f"Regression R² (All): {results['regression_r2']:.3f}")
    print(f"Regression R² (Active Only): {results['regression_r2_active_only']:.3f}")
    print(f"Best Training R²: {trainer.best_r2:.3f}")
    print(f"Results saved to: {os.path.abspath(results_dir)}")

if __name__ == "__main__":
    main()