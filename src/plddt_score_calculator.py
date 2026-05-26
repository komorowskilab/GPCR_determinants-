import os
import numpy as np

def compute_plddt_score(tsv_file):
    
    all_plddt_scores = []
    
    # Read the file and collect coordinates
    with open(tsv_file, "r") as f:
        for line in f:
            if not line.strip():
                continue
            
            
            # Collect all pLDDT scores from all atoms 
            if line.startswith("ATOM"):
                fields = line.split("\t")
                try:
                    plddt = float(fields[10])
                    all_plddt_scores.append(plddt)
                except (IndexError, ValueError):
                    continue
    
    
 
    # Calculate average pLDDT score
    if len(all_plddt_scores) > 0:
        return np.mean(all_plddt_scores)
    else:
        return None




# Main program
print("Looking for TSV files...")

all_results = {}

# Process each TSV file
for filename in os.listdir("."):
    if filename.endswith(".tsv"):
        print(f"Processing: {filename}")
        
        score = compute_plddt_score(filename)
        
        if score is None:
            print("  No valid ATOM lines found - skipping")
            continue
        
        sample_name = filename.replace(".tsv", "")
        all_results[sample_name] = score
        print(f"  Average pLDDT score: {score:.3f}")

# Check if we have any results
if len(all_results) == 0:
    print("No valid files found!")
    exit()


# Write output file
output_file = "plddt_score_table.tsv"
with open(output_file, "w") as f:
    # Write header
    f.write("Sample\tAverage_pLDDT_Score\n")
    
    # Write data for each sample
    for sample_name in sorted(all_results.keys()):
        score = all_results[sample_name]
        f.write(f"{sample_name}\t{score:.3f}\n")

print(f"\nDone! Created {output_file}")


