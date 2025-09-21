import os
import re

dlg_dir = "./dlg"
conf_dir = "./conf"
os.makedirs(conf_dir, exist_ok=True)

def extract_best_conformation(dlg_file, out_file):
    with open(dlg_file, "r") as f:
        lines = f.readlines()

    conf_blocks = []
    energies = []
    block = []
    recording = False
    current_energy = None

    for line in lines:
        # start of a model
        if line.startswith("MODEL"):
            recording = True
            block = [line]
            current_energy = None
            continue

        if recording:
            block.append(line)

            if "Estimated Free Energy of Binding" in line:
                match = re.search(r"=\s+(-?\d+\.\d+)", line)
                if match:
                    current_energy = float(match.group(1))

            if line.strip() == "ENDMDL":
                recording = False
                if current_energy is not None:
                    conf_blocks.append(block)
                    energies.append(current_energy)
                else:
                    print(f"⚠️ Skipping MODEL without energy in {dlg_file}")

    if not conf_blocks:
        print(f"⚠️ No valid conformations with energy found in {dlg_file}")
        return

    # select the block with lowest energy
    best_idx = min(range(len(energies)), key=lambda i: energies[i])
    best_conf = conf_blocks[best_idx]

    with open(out_file, "w") as out:
        out.writelines(best_conf)

    print(f"✅ {os.path.basename(dlg_file)} → {os.path.basename(out_file)} "
          f"(best binding energy {energies[best_idx]} kcal/mol)")

# process all dlg files
for dlg_file in os.listdir(dlg_dir):
    if dlg_file.endswith(".dlg"):
        in_path = os.path.join(dlg_dir, dlg_file)
        out_path = os.path.join(conf_dir, dlg_file.replace(".dlg", ".pdbqt"))
        extract_best_conformation(in_path, out_path)
