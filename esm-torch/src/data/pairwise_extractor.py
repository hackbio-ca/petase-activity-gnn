#!/usr/bin/env python3

import pandas as pd
import numpy as np
from itertools import combinations
import argparse
import os

def extract_pairwise_data(input_csv, output_file=None, sample_size=None, random_state=42):
    """
    Extract pairwise comparisons from activity CSV and create accession1, accession2, y_true format
    
    Parameters:
    - input_csv: path to input CSV with columns ['Accession code', 'Activity']
    - output_file: path for output file (default: 'pairwise_data.tsv')
    - sample_size: number of pairs to sample (default: all possible pairs)
    - random_state: random seed for reproducibility
    """
    
    # Read the input CSV
    df = pd.read_csv(input_csv)
    print(f"Loaded {len(df)} sequences from {input_csv}")
    
    # Clean the data - strip whitespace and handle any encoding issues
    df['Accession code'] = df['Accession code'].astype(str).str.strip()
    df['Activity'] = pd.to_numeric(df['Activity'], errors='coerce')
    
    # Remove any rows with NaN activities
    df = df.dropna(subset=['Activity'])
    print(f"After cleaning: {len(df)} valid sequences")
    
    # Validate columns
    required_cols = ['Accession code', 'Activity']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Generate all possible pairs
    pairs = list(combinations(range(len(df)), 2))
    print(f"Generated {len(pairs)} possible pairs")
    
    # Sample pairs if requested
    if sample_size and sample_size < len(pairs):
        np.random.seed(random_state)
        selected_indices = np.random.choice(len(pairs), size=sample_size, replace=False)
        pairs = [pairs[i] for i in selected_indices]
        print(f"Sampled {len(pairs)} pairs")
    
    # Create pairwise data
    pairwise_data = []
    
    # Debug: Print first few activities to see what we're working with
    print(f"\nFirst 10 activities:")
    for idx in range(min(10, len(df))):
        print(f"  {df.iloc[idx]['Accession code']}: {df.iloc[idx]['Activity']}")
    
    for i, j in pairs:
        accession_1 = str(df.iloc[i]['Accession code']).strip()
        activity_1 = float(df.iloc[i]['Activity'])
        accession_2 = str(df.iloc[j]['Accession code']).strip()
        activity_2 = float(df.iloc[j]['Activity'])
        
        # Binary comparison: 1 if accession_1 > accession_2, 0 otherwise
        y_true = int(1 if activity_1 > activity_2 else 0)
        
        # Debug: Print first few comparisons
        if len(pairwise_data) < 5:
            print(f"  Compare: {accession_1}({activity_1}) vs {accession_2}({activity_2}) → y_true={y_true}")
        
        pairwise_data.append({
            'accession_1': accession_1,
            'accession_2': accession_2,
            'y_true': y_true
        })
    
    # Create DataFrame with explicit data types
    pairwise_df = pd.DataFrame(pairwise_data)
    pairwise_df['accession_1'] = pairwise_df['accession_1'].astype(str)
    pairwise_df['accession_2'] = pairwise_df['accession_2'].astype(str)
    pairwise_df['y_true'] = pairwise_df['y_true'].astype(int)
    
    # Set output filename
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(input_csv))[0]
        output_file = f'{base_name}_pairwise.csv'
    
    # Save to CSV
    pairwise_df.to_csv(output_file, index=False, encoding='utf-8')
    
    # Debug: Print first few rows to verify format
    print(f"\nFirst 5 rows of output:")
    print(pairwise_df.head().to_string(index=False))
    
    # Print summary statistics
    print(f"\nPairwise data saved to: {output_file}")
    print(f"Total pairs: {len(pairwise_df)}")
    print(f"y_true distribution:")
    print(f"  0 (accession_1 <= accession_2): {sum(pairwise_df['y_true'] == 0)}")
    print(f"  1 (accession_1 > accession_2): {sum(pairwise_df['y_true'] == 1)}")
    
    # Show activity statistics
    activities = df['Activity'].values
    print(f"\nActivity statistics:")
    print(f"  Min: {activities.min():.6f}")
    print(f"  Max: {activities.max():.6f}")
    print(f"  Mean: {activities.mean():.6f}")
    print(f"  Unique values: {len(np.unique(activities))}")
    
    # If all activities are the same, warn user
    if len(np.unique(activities)) == 1:
        print(f"\nWARNING: All activities have the same value ({activities[0]})")
        print("This means all y_true values will be 0 (no meaningful comparisons)")
    
    return pairwise_df

def main():
    parser = argparse.ArgumentParser(description='Extract pairwise comparison data from activity CSV')
    parser.add_argument('input_csv', help='Input CSV file with activity data')
    parser.add_argument('--output', '-o', help='Output TSV file (default: auto-generated)')
    parser.add_argument('--sample_size', '-s', type=int, help='Number of pairs to sample')
    parser.add_argument('--random_state', '-r', type=int, default=42, help='Random seed')
    
    args = parser.parse_args()
    
    try:
        extract_pairwise_data(
            args.input_csv,
            args.output,
            args.sample_size,
            args.random_state
        )
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())

# Example usage:
# python extract_pairwise.py src/data/raw/seo_activity_data.csv
# python extract_pairwise.py src/data/raw/seo_activity_data.csv --output my_pairs.tsv --sample_size 5000