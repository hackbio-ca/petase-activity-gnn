import sys
import os
import matplotlib
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error, roc_auc_score
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import importlib

from src.data.data_loader import load_esm_embeddings_with_activity
from src.data.dataset import ESMActivityDataset, PairwiseESMDataset
from src.models.deep_gradual import ImprovedESMPredictor

print("train_baseline.py started")

def get_model_class(model_name):
    """
    Dynamically import and return the model class based on model name.
    
    Args:
        model_name (str): Name of the model to import (e.g., 'deep_gradual.ImprovedESMPredictor')
    
    Returns:
        class: The requested model class
    """
    module_path, class_name = model_name.rsplit('.', 1)
    module = importlib.import_module(f'src.models.{module_path}')
    return getattr(module, class_name)

# Deeper and more Gradual Reduction in Layer Sizes
def train_baseline_model(model_class_name):
    """
    Train a protein activity prediction model.
    
    Args:
        model_class_name (str): Name of the model class to use
    """
    # Create results directory
    results_dir = Path('results')
    results_dir.mkdir(exist_ok=True)

    # Get model class and set model name
    model_class = get_model_class(model_class_name)
    model_name = model_class.__name__.lower()

    # Load data
    print("Loading embeddings and activity data...")
    embeddings, activities, sequence_ids = load_esm_embeddings_with_activity()
    
    print(f"Dataset size: {len(embeddings)} sequences")
    print(f"Embedding dimension: {embeddings.shape[1]}")
    print(f"Activity range: {activities.min():.3f} to {activities.max():.3f}")
    
    # Split data
    X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(
        embeddings, activities, sequence_ids, 
        test_size=0.2, random_state=42
    )
    
    X_train, X_val, y_train, y_val, ids_train, ids_val = train_test_split(
        X_train, y_train, ids_train,
        test_size=0.2, random_state=42
    )
    
    print(f"Train: {len(X_train)}, Val: {len(X_val)}, Test: {len(X_test)}")
    
    # Create pairwise datasets
    train_dataset = PairwiseESMDataset(X_train, y_train)
    val_dataset = PairwiseESMDataset(X_val, y_val)
    test_dataset = PairwiseESMDataset(X_test, y_test)
    
    # Update loss function for pairwise learning
    criterion = nn.BCELoss()  # Binary Cross Entropy for pairwise comparison
    
    # Modified training loop
    for epoch in range(100):
        model.train()
        epoch_train_loss = 0
        for (batch_x1, batch_x2), batch_labels in train_loader:
            optimizer.zero_grad()
            predictions = model((batch_x1, batch_x2))
            loss = criterion(predictions, (batch_labels + 1) / 2)  # Convert -1/1 to 0/1
            loss.backward()
            optimizer.step()
            epoch_train_loss += loss.item()
            
        # Validation loop
        model.eval()
        val_correct = 0
        val_total = 0
        with torch.no_grad():
            for (batch_x1, batch_x2), batch_labels in val_loader:
                predictions = model((batch_x1, batch_x2))
                val_correct += ((predictions > 0.5).float() == (batch_labels + 1) / 2).sum().item()
                val_total += len(batch_labels)
        
        val_accuracy = val_correct / val_total
        print(f"Epoch {epoch}: Train Loss {avg_train_loss:.4f}, Val Accuracy {val_accuracy:.4f}")
    
    # Final evaluation on test set
    model.eval()
    test_predictions = []
    test_targets = []
    
    with torch.no_grad():
        for (batch_x1, batch_x2), batch_y in test_loader:
            predictions = model((batch_x1, batch_x2))
            test_predictions.extend(predictions.squeeze().numpy())
            test_targets.extend(batch_y.numpy())
    
    # Calculate metrics
    test_mse = mean_squared_error(test_targets, test_predictions)
    test_rmse = np.sqrt(test_mse)
    test_mae = mean_absolute_error(test_targets, test_predictions)
    test_r2 = r2_score(test_targets, test_predictions)
    
    print(f"\n=== Final Test Results ===")
    print(f"MSE: {test_mse:.4f}")
    print(f"RMSE: {test_rmse:.4f}")
    print(f"MAE: {test_mae:.4f}")
    print(f"R²: {test_r2:.4f}")
    
    # Save results
    torch.save(model.state_dict(), results_dir / f'{model_name}_model.pth')

    # Plot training curves
    plt.figure(figsize=(12, 4))

    plt.subplot(1, 2, 1)
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Val Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.legend()
    plt.title('Training Curves')

    plt.subplot(1, 2, 2)
    plt.scatter(test_targets, test_predictions, alpha=0.6)
    plt.plot([min(test_targets), max(test_targets)], [min(test_targets), max(test_targets)], 'r--')
    plt.xlabel('True Activity')
    plt.ylabel('Predicted Activity')
    plt.title(f'Test Set Predictions (R² = {test_r2:.3f})')

    plt.tight_layout()
    plt.savefig(results_dir / f'{model_name}_results.png', dpi=300, bbox_inches='tight')
    print(f"Results saved to {results_dir} as {model_name}_results.png")

    # Save detailed results
    results = {
        'test_mse': test_mse,
        'test_rmse': test_rmse,
        'test_mae': test_mae,
        'test_r2': test_r2,
        'train_losses': train_losses,
        'val_losses': val_losses,
        'val_r2_scores': val_r2_scores,
        'test_predictions': test_predictions,
        'test_targets': test_targets,
        'test_sequence_ids': ids_test
    }

    torch.save(results, results_dir / f'{model_name}_detailed_results.pth')

    return results

def evaluate_pairwise_metrics(model, test_loader):
    """Evaluate pairwise prediction metrics"""
    model.eval()
    correct = 0
    total = 0
    all_probs = []
    all_labels = []
    
    with torch.no_grad():
        for (batch_x1, batch_x2), batch_labels in test_loader:
            predictions = model((batch_x1, batch_x2))
            correct += ((predictions > 0.5).float() == (batch_labels + 1) / 2).sum().item()
            total += len(batch_labels)
            all_probs.extend(predictions.numpy())
            all_labels.extend(((batch_labels + 1) / 2).numpy())
    
    accuracy = correct / total
    auc = roc_auc_score(all_labels, all_probs)
    
    return {
        'accuracy': accuracy,
        'auc': auc
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train a protein activity prediction model')
    parser.add_argument('--model', type=str, default='deep_gradual.ImprovedESMPredictor',
                      help='Model to use (e.g., deep_gradual.ImprovedESMPredictor, residuals.ResidualESMPredictor)')
    
    args = parser.parse_args()
    print(f"Training with model: {args.model}")
    
    try:
        results = train_baseline_model(args.model)
    except Exception as e:
        import traceback
        print("Error during training:", e)
        traceback.print_exc()