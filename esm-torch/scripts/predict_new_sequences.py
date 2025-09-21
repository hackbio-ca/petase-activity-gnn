import sys
sys.path.append('.')

import torch
import esm
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from Bio import SeqIO

from src.models.baseline import ESMActivityPredictor
from src.models.pairwise_model import PairwiseESMPredictor

class ActivityPredictor:
    def __init__(self, model_path, model_type='baseline'):
        """
        Initialize predictor with trained model
        Args:
            model_path: Path to saved model weights
            model_type: 'baseline' or 'pairwise'
        """
        # Load smallest ESM model
        print("Loading ESM2 6M model...")
        self.esm_model, self.alphabet = esm.pretrained.esm2_t6_8M_UR50D()
        self.esm_model.eval()
        self.batch_converter = self.alphabet.get_batch_converter()
        
        # Load trained prediction model
        self.model_type = model_type
        if model_type == 'baseline':
            self.predictor = ESMActivityPredictor(input_dim=320)  # 6M model has 320 dim
        else:
            self.predictor = PairwiseESMPredictor(input_dim=320)
        
        self.predictor.load_state_dict(torch.load(model_path, map_location='cpu'))
        self.predictor.eval()
        
        print(f"Loaded {model_type} model from {model_path}")
    
    def extract_embeddings(self, sequences):
        """Extract ESM embeddings for list of sequences"""
        embeddings = []
        
        with torch.no_grad():
            for label, sequence in sequences:
                # Tokenize sequence
                batch_labels, batch_strs, batch_tokens = self.batch_converter([(label, sequence)])
                
                # Get embeddings from layer 6 (last layer of 6M model)
                results = self.esm_model(batch_tokens, repr_layers=[6])
                
                # Extract per-residue embeddings and average
                token_repr = results["representations"][6]
                seq_repr = token_repr[0, 1:len(sequence)+1].mean(0)  # Remove BOS/EOS, average
                
                embeddings.append(seq_repr.numpy())
                print(f"Processed {label}: {len(sequence)} residues")
        
        return np.array(embeddings)
    
    def predict_activities(self, sequences):
        """
        Predict activities for new sequences
        Args:
            sequences: List of (label, sequence) tuples
        Returns:
            DataFrame with predictions
        """
        # Extract embeddings
        embeddings = self.extract_embeddings(sequences)
        
        # Make predictions
        with torch.no_grad():
            predictions = self.predictor(torch.FloatTensor(embeddings))
            predictions = predictions.squeeze().numpy()
        
        # Format results
        labels = [label for label, _ in sequences]
        results_df = pd.DataFrame({
            'sequence_id': labels,
            'predicted_activity': predictions,
            'sequence_length': [len(seq) for _, seq in sequences]
        })
        
        return results_df

def load_sequences_from_fasta(fasta_path):
    """Load sequences from FASTA file"""
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append((record.id, str(record.seq)))
    return sequences

def main():
    parser = argparse.ArgumentParser(description='Predict protein activities from sequences')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--model', required=True, help='Path to trained model')
    parser.add_argument('--model_type', choices=['baseline', 'pairwise'], default='baseline')
    parser.add_argument('--output', default='predictions.csv', help='Output CSV file')
    
    args = parser.parse_args()
    
    # Load sequences
    print(f"Loading sequences from {args.fasta}")
    sequences = load_sequences_from_fasta(args.fasta)
    print(f"Loaded {len(sequences)} sequences")
    
    # Initialize predictor
    predictor = ActivityPredictor(args.model, args.model_type)
    
    # Make predictions
    print("Making predictions...")
    results = predictor.predict_activities(sequences)
    
    # Save results
    results.to_csv(args.output, index=False)
    print(f"Predictions saved to {args.output}")
    
    # Print summary
    print(f"\nPrediction Summary:")
    print(f"Mean predicted activity: {results['predicted_activity'].mean():.2f}")
    print(f"Activity range: {results['predicted_activity'].min():.2f} - {results['predicted_activity'].max():.2f}")
    print(f"Top 5 most active sequences:")
    print(results.nlargest(5, 'predicted_activity')[['sequence_id', 'predicted_activity']])

if __name__ == "__main__":
    main()