from Bio.PDB import PDBParser
import os
import csv


##UNFINISHED some protein xyz coordinates were not recorded########
##Prorgam read all the proteins from the afa file but is missing the xyz coordinates even though they are
#in the afa file

def extract_catalytic_atoms(pdb_dir, residue_map_csv, output_csv):
    parser = PDBParser(QUIET=True)

    # Define catalytic atoms of interest
    catalytic_atoms = {
        "S160": ["OG"],          # Ser hydroxyl
        "H237": ["NE2"],         # His imidazole nitrogen
        "D206": ["OD1"]          #Asp hydroxyl 
    }

    # Load residue mapping from CSV
    residue_map = {}  # {ProteinID: {Residue: UngappedPosition}}
    with open(residue_map_csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pid = row["ProteinID"]
            res = row["Residue"]
            pos = int(row["Position_in_ungapped_sequence"])
            residue_map.setdefault(pid, {})[res] = pos

    results = []

    # Iterate through PDB files
    for pdb_file in os.listdir(pdb_dir):
        if not pdb_file.endswith(".pdb"):
            continue

        pdb_id = os.path.splitext(pdb_file)[0]

        # Find corresponding ID in residue_map (AFA IDs end with "_A")
        matching_pid = None
        for pid in residue_map:
            if pid.startswith(pdb_id):
                matching_pid = pid
                break
        if not matching_pid:
            continue

        pdb_path = os.path.join(pdb_dir, pdb_file)
        structure = parser.get_structure(pdb_id, pdb_path)

        for res, pos in residue_map[matching_pid].items():
            atoms_to_find = catalytic_atoms.get(res, [])
            found = False
            for model in structure:
                for chain in model:
                    for residue in chain:
                        het, resseq, icode = residue.id
                        if het == " " and resseq == pos:  # match residue number
                            for atom_name in atoms_to_find:
                                if atom_name in residue:
                                    atom = residue[atom_name]
                                    x, y, z = atom.coord
                                    results.append([pdb_id, res, atom_name, x, y, z])
                                    found = True
            if not found:
                # Residue or atom missing
                for atom_name in atoms_to_find:
                    results.append([pdb_id, res, atom_name, None, None, None])

    # Write output CSV
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ProteinID", "Residue", "Atom", "X", "Y", "Z"])
        writer.writerows(results)


# === Example usage ===
pdb_dir = "C:/cygwin64/home/jared/jared/PETaseProject/aligned_structures/"  # directory with pdb files
residue_map_csv = "C:/cygwin64/home/jared/jared/PETaseProject/residue_ungapped_positions.csv"
output_csv = "C:/cygwin64/home/jared/jared/PETaseProject/catalytic_atom_coordinates6csv"

extract_catalytic_atoms(pdb_dir, residue_map_csv, output_csv)

print(f"XYZ coordinates saved to {output_csv}")
