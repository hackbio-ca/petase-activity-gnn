import os
import csv

conf_dir = "./conf"
out_csv = "atom_coords.csv"

# Atoms to extract
target_atoms = {4, 5, 9, 20, 21, 25}

rows = []

for pdbqt_file in os.listdir(conf_dir):
    if pdbqt_file.endswith(".pdbqt"):
        protein = os.path.splitext(pdbqt_file)[0]
        file_path = os.path.join(conf_dir, pdbqt_file)

        with open(file_path, "r") as f:
            for line in f:
                if line.startswith("ATOM"):
                    parts = line.split()
                    atom_num = int(parts[1])
                    if atom_num in target_atoms:
                        atom_name = parts[2]
                        x = float(parts[6])
                        y = float(parts[7])
                        z = float(parts[8])
                        rows.append([protein, atom_num, atom_name, x, y, z])

# write CSV
with open(out_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["protein", "atom_num", "atom", "x", "y", "z"])
    writer.writerows(rows)

print(f"✅ CSV written to {out_csv} with {len(rows)} entries")
