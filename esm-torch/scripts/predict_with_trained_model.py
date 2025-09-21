import sys
sys.path.append('.')

import torch
import esm
import numpy as np
import pandas as pd
from pathlib import Path
import argparse
from Bio import SeqIO
import csv
import os
import time
from src.utils.live_writer import LiveResultsWriter


class LiveResultsWriter:
    def __init__(self, output_path, buffer_size=100):
        self.output_path = output_path
        self.buffer_size = buffer_size
        self.buffer = []
        self.total_written = 0
        
        # Create output file with headers if it doesn't exist
        if not os.path.exists(output_path):
            with open(output_path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['sequence_id', 'predicted_activity', 'sequence_length', 'timestamp'])
        
        # Create progress file
        self.progress_file = output_path.replace('.csv', '_progress.txt')
        
    def add_result(self, sequence_id, activity, seq_length):
        result = {
            'sequence_id': sequence_id,
            'predicted_activity': activity,
            'sequence_length': seq_length,
            'timestamp': time.time()
        }
        self.buffer.append(result)
        
        if len(self.buffer) >= self.buffer_size:
            self.flush()
    
    def flush(self):
        if not self.buffer:
            return
            
        # Append to CSV
        with open(self.output_path, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['sequence_id', 'predicted_activity', 'sequence_length', 'timestamp'])
            writer.writerows(self.buffer)
        
        self.total_written += len(self.buffer)
        self.buffer = []
        
        # Update progress
        self.update_progress()
    
    def update_progress(self):
        with open(self.progress_file, 'w') as f:
            f.write(f"Processed: {self.total_written}\n")
            f.write(f"Last update: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            if hasattr(self, 'start_time') and self.total_written > 0:
                elapsed = time.time() - self.start_time
                rate = self.total_written / elapsed
                f.write(f"Rate: {rate:.2f} sequences/second\n")
                f.write(f"Elapsed: {elapsed/3600:.2f} hours\n")
    
    def start_timing(self):
        self.start_time = time.time()
    
    def finish(self):
        self.flush()
        print(f"Batch completed! Total processed in this batch: {self.total_written}")

class FlexibleActivityPredictor:
    def __init__(self, model_path, device=None):
        """
        Initialize predictor with your trained model
        Automatically detects model type and architecture
        """
        if device is None:
            self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        else:
            self.device = device
            
        print(f"Using device: {self.device}")
        print(f"Loading trained model from {model_path}")
        
        # Load model checkpoint
        self.checkpoint = torch.load(model_path, map_location='cpu')
        
        # Determine model architecture
        self.is_pairwise = any('shared_network' in key for key in self.checkpoint.keys())
        
        # Determine input dimension from first layer
        first_layer_key = list(self.checkpoint.keys())[0]
        if 'weight' in first_layer_key:
            self.input_dim = self.checkpoint[first_layer_key].shape[1]
        else:
            # Fallback - assume it matches your ESM model
            self.input_dim = 1280  # Default for ESM2 33 layer model
        
        print(f"Detected input dimension: {self.input_dim}")
        print(f"Model type: {'Pairwise' if self.is_pairwise else 'Baseline'}")
        
        # Load appropriate ESM model based on input dimension
        if self.input_dim == 320:
            print("Loading ESM2 6M model (320 dim)")
            self.esm_model, self.alphabet = esm.pretrained.esm2_t6_8M_UR50D()
            self.esm_layer = 6
        elif self.input_dim == 480:
            print("Loading ESM2 35M model (480 dim)")
            self.esm_model, self.alphabet = esm.pretrained.esm2_t12_35M_UR50D()
            self.esm_layer = 12
        elif self.input_dim == 640:
            print("Loading ESM2 150M model (640 dim)")
            self.esm_model, self.alphabet = esm.pretrained.esm2_t30_150M_UR50D()
            self.esm_layer = 30
        elif self.input_dim == 1280:
            print("Loading ESM2 650M model (1280 dim)")
            self.esm_model, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
            self.esm_layer = 33
        else:
            raise ValueError(f"Unknown input dimension: {self.input_dim}")
        
        self.esm_model = self.esm_model.to(self.device)
        self.esm_model.eval()
        self.batch_converter = self.alphabet.get_batch_converter()
        
        # Load prediction model architecture
        self._load_prediction_model()
        
    def _load_prediction_model(self):
        """Load the trained prediction model"""
        # Check model architecture from state dict keys
        state_dict_keys = list(self.checkpoint.keys())
        
        if any('feature_extractor' in key for key in state_dict_keys):
            # Enhanced multi-task model
            from src.models.enhanced_model import EnhancedPetaseModel
            self.predictor = EnhancedPetaseModel(input_dim=self.input_dim)
            print("Loading Enhanced PETase Model")
        elif self.is_pairwise:
            from src.models.pairwise_model import PairwiseESMPredictor
            self.predictor = PairwiseESMPredictor(input_dim=self.input_dim)
        else:
            from src.models.baseline import ESMActivityPredictor
            self.predictor = ESMActivityPredictor(input_dim=self.input_dim)
        
        # Load weights
        self.predictor.load_state_dict(self.checkpoint)
        self.predictor = self.predictor.to(self.device)
        self.predictor.eval()
        print("Prediction model loaded successfully")
    
    def predict_activities_live(self, sequences, writer):
        """Predict activities with live output"""
        for i, (label, sequence) in enumerate(sequences):
            try:
                # Extract embedding for single sequence
                batch_labels, batch_strs, batch_tokens = self.batch_converter([(label, sequence)])
                batch_tokens = batch_tokens.to(self.device)
                
                with torch.no_grad():
                    results = self.esm_model(batch_tokens, repr_layers=[self.esm_layer])
                    seq_repr = results["representations"][self.esm_layer][0, 1:len(sequence)+1].mean(0)
                    
                    # Predict activity
                    if hasattr(self.predictor, 'regressor'):
                        prediction = self.predictor(seq_repr.unsqueeze(0).to(self.device), task='regression')
                    else:
                        prediction = self.predictor(seq_repr.unsqueeze(0).to(self.device))
                    
                    activity = prediction.item()
                    
                    # Write result immediately
                    writer.add_result(label, activity, len(sequence))
                    
                    if i % 50 == 0 and i > 0:
                        print(f"  Processed {i+1}/{len(sequences)}: {label} -> {activity:.3f}")
                        
            except Exception as e:
                print(f"Error processing {label}: {e}")
                writer.add_result(label, -999, len(sequence))  # Error marker

    def process_large_fasta_live(self, fasta_path, output_path, chunk_size=500):
        """Process large FASTA with live updates and resumption capability"""
        
        # Check for existing progress
        processed_ids = set()
        if os.path.exists(output_path):
            print("Found existing output file, checking progress...")
            with open(output_path, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    processed_ids.add(row['sequence_id'])
            print(f"Resuming from {len(processed_ids)} already processed sequences")
        
        chunk_sequences = []
        total_count = 0
        skipped_count = 0
        
        print("Starting live processing...")
        
        for record in SeqIO.parse(fasta_path, "fasta"):
            # Skip if already processed
            if record.id in processed_ids:
                skipped_count += 1
                continue
                
            chunk_sequences.append((record.id, str(record.seq)))
            total_count += 1
            
            if len(chunk_sequences) >= chunk_size:
                print(f"\nProcessing chunk starting at sequence {total_count-chunk_size+1}")
                
                writer = LiveResultsWriter(output_path)
                writer.start_timing()
                
                self.predict_activities_live(chunk_sequences, writer)
                writer.finish()
                
                chunk_sequences = []
                
                # Memory cleanup
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
        
        # Process remaining sequences
        if chunk_sequences:
            print(f"\nProcessing final chunk: {len(chunk_sequences)} sequences")
            writer = LiveResultsWriter(output_path)
            writer.start_timing()
            self.predict_activities_live(chunk_sequences, writer)
            writer.finish()
        
        print(f"\nCompleted! Total new sequences: {total_count}, Skipped: {skipped_count}")

    def predict_activities(self, sequences):
        """Original batch prediction method for compatibility"""
        embeddings = []
        
        with torch.no_grad():
            for label, sequence in sequences:
                # Tokenize sequence
                batch_labels, batch_strs, batch_tokens = self.batch_converter([(label, sequence)])
                batch_tokens = batch_tokens.to(self.device)
                
                # Get embeddings from appropriate layer
                results = self.esm_model(batch_tokens, repr_layers=[self.esm_layer])
                
                # Extract mean-pooled sequence representation
                token_repr = results["representations"][self.esm_layer]
                seq_repr = token_repr[0, 1:len(sequence)+1].mean(0)
                
                embeddings.append(seq_repr.cpu().numpy())
        
        embeddings = np.array(embeddings)
        
        # Make predictions
        with torch.no_grad():
            if hasattr(self.predictor, 'regressor'):
                # Enhanced model - use regression task
                predictions = self.predictor(torch.FloatTensor(embeddings).to(self.device), task='regression')
            elif self.is_pairwise:
                # For pairwise model, use single protein prediction
                predictions = self.predictor.predict_single(torch.FloatTensor(embeddings).to(self.device))
            else:
                predictions = self.predictor(torch.FloatTensor(embeddings).to(self.device))
            
            predictions = predictions.squeeze().cpu().numpy()
        
        # Format results
        labels = [label for label, _ in sequences]
        results_df = pd.DataFrame({
            'sequence_id': labels,
            'predicted_activity': predictions,
            'sequence_length': [len(seq) for _, seq in sequences]
        })
        
        return results_df

def load_sequences_from_fasta(fasta_path, max_sequences=None):
    """Load sequences from FASTA file with optional limit"""
    sequences = []
    for i, record in enumerate(SeqIO.parse(fasta_path, "fasta")):
        sequences.append((record.id, str(record.seq)))
        if max_sequences and i + 1 >= max_sequences:
            break
    return sequences

def main():
    parser = argparse.ArgumentParser(description='Predict protein activities using trained model')
    parser.add_argument('--fasta', required=True, help='Input FASTA file')
    parser.add_argument('--model', 
                       default='/home/ec2-user/petase-activity-gnn/esm-torch/results/two_stage_p2/enhanced_petase_model.pth',
                       help='Path to trained model')
    parser.add_argument('--output', default='predictions.csv', help='Output CSV file')
    parser.add_argument('--max_sequences', type=int, help='Limit number of sequences (for testing)')
    parser.add_argument('--chunk_size', type=int, default=500, help='Chunk size for live processing')
    parser.add_argument('--live', action='store_true', help='Use live processing with resume capability')
    
    args = parser.parse_args()
    
    # Initialize predictor
    predictor = FlexibleActivityPredictor(args.model)
    
    if args.live:
        # Use live processing
        print(f"Starting live processing of {args.fasta}")
        predictor.process_large_fasta_live(args.fasta, args.output, args.chunk_size)
    else:
        # Original batch processing
        print(f"Loading sequences from {args.fasta}")
        sequences = load_sequences_from_fasta(args.fasta, args.max_sequences)
        print(f"Loaded {len(sequences)} sequences")
        
        # Process in batches
        all_results = []
        batch_size = args.chunk_size
        
        for i in range(0, len(sequences), batch_size):
            batch_sequences = sequences[i:i+batch_size]
            print(f"Processing batch {i//batch_size + 1}/{(len(sequences)-1)//batch_size + 1}")
            
            batch_results = predictor.predict_activities(batch_sequences)
            all_results.append(batch_results)
        
        # Combine all results
        final_results = pd.concat(all_results, ignore_index=True)
        
        # Save results
        final_results.to_csv(args.output, index=False)
        print(f"Predictions saved to {args.output}")
        
        # Print summary
        print(f"\nPrediction Summary:")
        print(f"Total sequences processed: {len(final_results)}")
        print(f"Mean predicted activity: {final_results['predicted_activity'].mean():.2f}")
        print(f"Activity range: {final_results['predicted_activity'].min():.2f} - {final_results['predicted_activity'].max():.2f}")
        print(f"Top 5 most active sequences:")
        print(final_results.nlargest(5, 'predicted_activity')[['sequence_id', 'predicted_activity']])

if __name__ == "__main__":
    main()