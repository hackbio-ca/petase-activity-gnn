import torch
from torch.utils.data import Dataset
import numpy as np

class ESMActivityDataset(Dataset):
    """Dataset for individual protein embeddings and their activities."""
    def __init__(self, embeddings, activities):
        self.embeddings = torch.FloatTensor(embeddings)
        self.activities = torch.FloatTensor(activities)
    
    def __len__(self):
        return len(self.embeddings)
    
    def __getitem__(self, idx):
        return self.embeddings[idx], self.activities[idx]

class PairwiseESMDataset(Dataset):
    """Dataset for pairwise protein comparisons."""
    def __init__(self, embeddings, activities):
        self.embeddings = torch.FloatTensor(embeddings)
        self.activities = torch.FloatTensor(activities)
        self.pairs = self._generate_pairs()
    
    def _generate_pairs(self):
        """Generate all possible pairs and their relative activity labels"""
        n = len(self.embeddings)
        pairs = []
        for i in range(n):
            for j in range(i + 1, n):
                if abs(self.activities[i] - self.activities[j]) > 1e-6:  # Only include pairs with different activities
                    # Label: 1 if protein i is more active than protein j, 0 otherwise
                    label = 1.0 if self.activities[i] > self.activities[j] else 0.0
                    pairs.append((i, j, label))
        return pairs
    
    def __len__(self):
        return len(self.pairs)
    
    def __getitem__(self, idx):
        i, j, label = self.pairs[idx]
        return self.embeddings[i], self.embeddings[j], torch.FloatTensor([label])