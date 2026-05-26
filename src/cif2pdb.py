import pymol
pymol.finish_launching()

from pymol import cmd
import os

# Current working directory
folder = "."  
files = [f for f in os.listdir(folder) if f.endswith(".cif")]

for i, f in enumerate(files):
    temp_name = f"tmp_{i}"
    obj_name = os.path.splitext(f)[0]

    print(f"\n=== Processing {f} ===")

    # Load the CIF file
    cmd.load(f, temp_name)

    # Get extents (min, max)
    minc, maxc = cmd.get_extent(temp_name)

    # Compute bounding-box center
    center_x = (minc[0] + maxc[0]) / 2.0
    center_y = (minc[1] + maxc[1]) / 2.0
    center_z = (minc[2] + maxc[2]) / 2.0

    # Translation vector to move center → (0,0,0)
    # Molecule is moved to origin to have decent number in its coordinates
    tx = -center_x
    ty = -center_y
    tz = -center_z

    print(f"Translation: [{tx:.3f}, {ty:.3f}, {tz:.3f}]")

    # Apply the translation to coordinates
    cmd.translate([tx, ty, tz], temp_name, camera=0)

    # Sort atoms
    cmd.sort(temp_name)

    # Save as PDB in the same directory
    pdb_path = obj_name + ".pdb"
    cmd.save(pdb_path, temp_name, format="pdb")


    print(f"Saved PDB → {pdb_path}")

    # Clean up
    cmd.delete(temp_name)
