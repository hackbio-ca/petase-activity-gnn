import pandas as pd
import itertools
import numpy as np

# Load input CSV
df = pd.read_csv("atom_coords.csv")

# Prepare output
output_rows = []

# Group by Protein_ID
for protein_id, group in df.groupby("Protein_ID"):
    atoms = group[["Residue_ID", "Atom_ID", "x", "y", "z"]].values
    
    # Generate all unique pairs of atoms
    for (res1, atom1, x1, y1, z1), (res2, atom2, x2, y2, z2) in itertools.combinations(atoms, 2):
        distance = np.sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
        
        output_rows.append({
            "Protein_ID": protein_id,
            "First_Atom": f"{res1}_{atom1}",
            "Second_Atom": f"{res2}_{atom2}",
            "x1": x1, "y1": y1, "z1": z1,
            "x2": x2, "y2": y2, "z2": z2,
            "distance": distance
        })

# Save to CSV
output_df = pd.DataFrame(output_rows)
output_df.to_csv("atom_pair_distances.csv", index=False)

print("✅ Output saved to atom_pair_distances.csv")
