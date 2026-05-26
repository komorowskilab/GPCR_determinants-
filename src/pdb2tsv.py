
import glob
import re

for pdb in glob.glob("*.pdb"):
    out = pdb.replace(".pdb", ".tsv")

    with open(pdb) as fin, open(out, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "HETATM")):
                tsv_line = re.sub(r'\s+', '\t', line.strip())
                fout.write(tsv_line + "\n")

print("Done! Converted all PDB files to TSV.")
