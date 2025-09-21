import sys
sys.path.append('.')

import torch
import esm
import numpy as np
from pathlib import Path
from Bio import SeqIO
import argparse

def extract_embeddings_small_esm(fasta_file, output_dir):
    """Extract embeddings using ESM2 6M model"""
    # Load smallest ESM model
    print("Loading ESM2 6M model...")
    model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Process FASTA file
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append((record.id, str(record.seq)))
    
    print(f"Processing {len(sequences)} sequences...")
    
    for label, sequence in sequences:
        # Tokenize
        batch_labels, batch_strs, batch_tokens = batch_converter([(label, sequence)])
        
        # Extract embeddings
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[6])  # Last layer of 6M model
        
        # Get representations
        token_representations = results["representations"][6]
        
        # Mean pooled sequence representation
        seq_repr = token_representations[0, 1:len(sequence)+1].mean(0)
        
        # Per-token representations
        per_tok_repr = token_representations[0, 1:len(sequence)+1]
        
        # Save in same format as original
        output_data = {
            'label': label,
            'mean_representations': {6: seq_repr},
            'representations': {6: per_tok_repr}
        }
        
        output_file = output_dir / f"{label.replace('/', '_')}.pt"
        torch.save(output_data, output_file)
        
        print(f"Processed {label}: {len(sequence)} residues, saved to {output_file.name}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--output', required=True, help='Output directory for embeddings')
    
    args = parser.parse_args()
    extract_embeddings_small_esm(args.fasta, args.output)

if __name__ == "__main__":
    main()