import torch
from pathlib import Path

def inspect_model(model_path):
    """Inspect the saved model to understand its architecture"""
    print(f"Inspecting model: {model_path}")
    
    # Load model state dict
    checkpoint = torch.load(model_path, map_location='cpu')
    
    print("Model state dict keys:")
    for key, tensor in checkpoint.items():
        print(f"  {key}: {tensor.shape}")
    
    # Try to determine input dimension from first layer
    first_layer_key = list(checkpoint.keys())[0]
    if 'weight' in first_layer_key:
        input_dim = checkpoint[first_layer_key].shape[1]
        output_dim = checkpoint[first_layer_key].shape[0]
        print(f"\nDetected input dimension: {input_dim}")
        print(f"First layer output dimension: {output_dim}")
    
    # Check if it's a pairwise model (has shared_network prefix)
    has_shared_network = any('shared_network' in key for key in checkpoint.keys())
    print(f"Is pairwise model: {has_shared_network}")
    
    return checkpoint, input_dim if 'input_dim' in locals() else None, has_shared_network

if __name__ == "__main__":
    model_path = "/home/ec2-user/petase-activity-gnn/esm-torch/results/two_stage_p2/enhanced_petase_model.pth"
    inspect_model(model_path)