
from Bio import SeqIO
import csv

# --- Step 1: Map residues in template to alignment columns ---
def map_residue_to_alignment(msa_file, template_id, residues_of_interest):
    records = list(SeqIO.parse(msa_file, "fasta"))
    
    # Find template sequence
    template = None
    for record in records:
        if template_id in record.id:
            template = record
            break
    if template is None:
        raise ValueError(f"Template {template_id} not found in alignment.")
    
    alignment_seq = str(template.seq)
    ungapped_pos = 0
    pos_map = {}
    for aln_index, aa in enumerate(alignment_seq, start=1):
        if aa != "-":
            ungapped_pos += 1
            pos_map[ungapped_pos] = (aln_index, aa)
    
    # Map residues of interest to alignment columns
    results = {}
    for res in residues_of_interest:
        aa = res[0]
        pos = int(res[1:])
        if pos not in pos_map:
            results[res] = None
        else:
            aln_index, found_aa = pos_map[pos]
            results[res] = aln_index  # alignment column
    return results


# --- Step 2: For each protein, compute ungapped position up to that column ---
def compute_ungapped_positions(msa_file, residue_alignment_map, output_csv):
    records = list(SeqIO.parse(msa_file, "fasta"))
    results = []

    for record in records:
        seq_id = record.id
        seq = str(record.seq)

        for residue, aln_col in residue_alignment_map.items():
            if aln_col is None:
                results.append([seq_id, residue, None])
                continue

            # Substring up to the alignment column
            substring = seq[:aln_col]  # up to and including aln_col
            ungapped_length = len(substring.replace("-", ""))  # remove gaps
            results.append([seq_id, residue, ungapped_length])

    # Save CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ProteinID", "Residue", "Position_in_ungapped_sequence"])
        writer.writerows(results)


# === Example usage ===
msa_file = "C:/cygwin64/home/jared/jared/PETaseProject/aligned_structures/PETase.afa"
template_id = "WP_054022242_A"
residues_of_interest = ["S160", "D206", "H237"]
output_csv = "C:/cygwin64/home/jared/jared/PETaseProject/residue_ungapped_positions.csv"

# Map residues in template to alignment columns
residue_alignment_map = map_residue_to_alignment(msa_file, template_id, residues_of_interest)

# Compute ungapped positions for all proteins
compute_ungapped_positions(msa_file, residue_alignment_map, output_csv)

print(f"Results saved to {output_csv}")
