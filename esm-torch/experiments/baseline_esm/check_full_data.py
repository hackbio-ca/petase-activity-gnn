import sys
sys.path.append('../../')

from src.data.data_loader import inspect_full_dataset

if __name__ == "__main__":
    embeddings, activities, seq_ids = inspect_full_dataset()