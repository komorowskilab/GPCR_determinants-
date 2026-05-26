from Bio.PDB import PDBParser, PDBIO

# -------- paths --------
pdb_path = "human_model.pdb"
values_path = "normalized_distances.txt"
output_pdb = "human_model_REPLACED.pdb"

# -------- read residue-level values --------
with open(values_path) as f:
    line = f.readline().strip()
    parts = line.split()

# skip first column (model name)
residue_values = [float(x) for x in parts[1:]]

# -------- parse PDB --------
parser = PDBParser(QUIET=True)
structure = parser.get_structure("model", pdb_path)

res_idx = 0

# -------- replace B-factors per residue --------
for model in structure:
    for chain in model:
        for residue in chain:
            # ATOM records only
            if residue.id[0] != " ":
                continue

            if res_idx >= len(residue_values):
                break

            new_b = residue_values[res_idx]

            # replace for ALL atoms in this residue
            for atom in residue:
                atom.set_bfactor(new_b)

            res_idx += 1

# -------- write new PDB --------
io = PDBIO()
io.set_structure(structure)
io.save(output_pdb)

print(f"Replaced B-factors for {res_idx} residues")
print("Saved:", output_pdb)
