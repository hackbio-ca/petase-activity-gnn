import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../../')))
print("sys.path:", sys.path)
print("Script started...")

try:
    import torch
    import torch.nn as nn
    print("PyTorch imported")
except Exception as e:
    print(f"PyTorch import error: {e}")
    exit()

try:
    from src.data.data_loader import load_esm_embeddings_with_activity
    print("Data loader imported")
except Exception as e:
    print(f"Data loader import error: {e}")
    exit()

try:
    embeddings, activities, sequence_ids = load_esm_embeddings_with_activity()
    print(f"Data loaded: {len(embeddings)} sequences")
except Exception as e:
    print(f"Data loading error: {e}")
    exit()

print("All basic components working. Issue likely in missing model/dataset classes.")