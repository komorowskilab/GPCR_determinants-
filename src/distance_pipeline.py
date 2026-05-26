import os
import numpy as np
from collections import defaultdict

def compute_distances(tsv_file):
    """Calculate distances from COG of residues to COG ligand."""
    ligand_coords = []
    residue_atoms = defaultdict(list)

    with open(tsv_file) as f:
        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip().split("\t")

            if line.startswith("HETATM"):
                try:
                    ligand_coords.append([float(fields[5]),
                                           float(fields[6]),
                                           float(fields[7])])
                except:
                    pass

            elif line.startswith("ATOM"):
                try:
                    resnum = fields[5]
                    residue_atoms[resnum].append([
                        float(fields[6]),
                        float(fields[7]),
                        float(fields[8])
                    ])
                except:
                    pass

    if not ligand_coords:
        return None

    ligand_center = np.mean(np.array(ligand_coords), axis=0)

    distances = {}
    for resnum, coords in residue_atoms.items():
        residue_center = np.mean(np.array(coords), axis=0)
        distances[resnum] = np.linalg.norm(residue_center - ligand_center)

    return distances



# Main program
print("Looking for TSV files...")

all_results = {}

# Process each TSV file
for filename in os.listdir("."):
    if filename.endswith(".tsv"):
        print(f"Processing: {filename}")
        
        distances = compute_distances(filename)
        
        if distances is None:
            print("  No ligand found - skipping")
            continue
        
        sample_name = filename.replace(".tsv", "")
        all_results[sample_name] = distances
        print(f"  Found {len(distances)} residues")

# Check if we have any results
if len(all_results) == 0:
    print("No valid files found!")
    exit()

# Get all unique residue numbers and sort them
all_residues = set()
for distances in all_results.values():
    all_residues.update(distances.keys())

all_residues = sorted(all_residues, key=lambda x: float(x))

# Write output file
output_file = "distance_table_v5.tsv"
with open(output_file, "w") as f:
    # Write header
    f.write("Sample\t" + "\t".join(all_residues) + "\n")
    
    # Write data for each sample
    for sample_name in sorted(all_results.keys()):
        distances = all_results[sample_name]
        
        row = [sample_name]
        for resnum in all_residues:
            if resnum in distances:
                row.append(f"{distances[resnum]:.3f}")
            else:
                row.append("")
        
        f.write("\t".join(row) + "\n")

print(f"\nDone! Created {output_file}")
print(f"Samples: {len(all_results)}")
print(f"Residues: {len(all_residues)}")
