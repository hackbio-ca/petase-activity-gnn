import torch
from torch.utils.data import Dataset
import numpy as np
from itertools import combinations

class PairwiseDataset(Dataset):
    def __init__(self, embeddings, activities, sequence_ids, max_pairs=None):
        """
        Create pairwise dataset from embeddings and activities
        max_pairs: limit number of pairs to prevent memory issues
        """
        self.pairs = []
        self.targets = []
        self.pair_ids = []
        
        # Generate all possible pairs
        indices = list(range(len(embeddings)))
        pair_indices = list(combinations(indices, 2))
        
        # Limit pairs if specified
        if max_pairs and len(pair_indices) > max_pairs:
            np.random.seed(42)
            pair_indices = np.random.choice(len(pair_indices), max_pairs, replace=False)
            pair_indices = [pair_indices[i] for i in range(max_pairs)]
        
        for i, j in pair_indices:
            emb1, emb2 = embeddings[i], embeddings[j]
            activity_diff = activities[i] - activities[j]
            
            self.pairs.append((emb1, emb2))
            self.targets.append(activity_diff)
            self.pair_ids.append((sequence_ids[i], sequence_ids[j]))
    
    def __len__(self):
        return len(self.pairs)
    
    def __getitem__(self, idx):
        emb1, emb2 = self.pairs[idx]
        return torch.FloatTensor(emb1), torch.FloatTensor(emb2), torch.FloatTensor([self.targets[idx]])

class HybridDataset(Dataset):
    """Dataset that provides both individual and pairwise data"""
    def __init__(self, embeddings, activities, pairwise_dataset):
        self.individual_data = torch.FloatTensor(embeddings)
        self.individual_targets = torch.FloatTensor(activities)
        self.pairwise_data = pairwise_dataset
        
    def __len__(self):
        return len(self.individual_data)
    
    def get_individual(self, idx):
        return self.individual_data[idx], self.individual_targets[idx]
    
    def get_random_pair(self):
        idx = np.random.randint(len(self.pairwise_data))
        return self.pairwise_data[idx]