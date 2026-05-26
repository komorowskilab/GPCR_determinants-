
# Distance table and aligned fasta file should be adjusted accordingly, in the main usage part

import pandas as pd
import numpy as np

def parse_aligned_sequences(tsv_file):
    """
    Parse aligned sequences from TSV or FASTA file.
    Extracts only the core accession ID (up to the first slash or space)
    and lowercases it to match distance table IDs.
    """
    sequences = {}

    if tsv_file.endswith(('.fasta', '.fa', '.faa')):
        current_id = None
        current_seq = []

        with open(tsv_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id:
                        sequences[current_id] = ''.join(current_seq)

                    # Extract only the accession part before '/' or space, lowercase
                    header = line[1:].split()[0]       # e.g., 'XP_076580365/26-155'
                    core_id = header.split('/')[0].lower()  # e.g., 'xp_076580365'
                    current_id = core_id
                    current_seq = []
                elif line:
                    current_seq.append(line)

            if current_id:
                sequences[current_id] = ''.join(current_seq)

    else:
        # TSV format: ID\tSequence
        with open(tsv_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        seq_id = parts[0].split('_')[0].lower()  # e.g., 'xp_071635654'
                        seq = parts[1]
                        sequences[seq_id] = seq

    return sequences

def create_position_mapping(aligned_seq):
    """
    Create mapping from original (ungapped) position to aligned position.
    
    Returns: list of tuples (original_pos, aligned_pos, amino_acid)
    """
    mapping = []
    original_pos = 0
    
    for aligned_pos, aa in enumerate(aligned_seq):
        if aa not in ['.', '-', 'X']:  # Skip gaps
            mapping.append((original_pos, aligned_pos, aa))
            original_pos += 1
    
    return mapping

def align_distances_to_positions(distance_df, aligned_sequences, aligned_length=135):
    aligned_data = []

    for idx, row in distance_df.iterrows():
        distance_id = row['Clean_ID'].lower()  # e.g., xp_076580365_temnothorax

        # Find a matching FASTA ID that is a prefix of distance_id
        matched_fasta_id = next((fid for fid in aligned_sequences if distance_id.startswith(fid)), None)

        if not matched_fasta_id:
            print(f"Warning: {distance_id} not found in aligned sequences, skipping...")
            continue

        aligned_seq = aligned_sequences[matched_fasta_id]
        original_distances = row.drop('Clean_ID').values
        mapping = create_position_mapping(aligned_seq)
        aligned_distances = np.full(aligned_length, np.nan)

        for i, (orig_pos, aligned_pos, aa) in enumerate(mapping):
            if i < len(original_distances):
                aligned_distances[aligned_pos] = original_distances[i]

        aligned_row = {'Clean_ID': row['Clean_ID']}
        for i in range(aligned_length):
            aligned_row[f'P{i+1}'] = aligned_distances[i]

        aligned_data.append(aligned_row)

    return pd.DataFrame(aligned_data)


def create_detailed_mapping_all_species(distance_df, aligned_sequences):
    """
    Create a unified detailed mapping table for all species:
    Columns:
        Clean_ID, FASTA_ID, Original_Pos, Aligned_Pos, AA, Distance
    """
    all_mappings = []

    for idx, row in distance_df.iterrows():
        seq_id = row['Clean_ID'].lower()

        # Match FASTA ID the same way align_distances_to_positions() does
        matched_fasta_id = next(
            (fid for fid in aligned_sequences if seq_id.startswith(fid)),
            None
        )

        if not matched_fasta_id:
            print(f"Warning: {seq_id} missing in aligned sequences — skipped")
            continue

        aligned_seq = aligned_sequences[matched_fasta_id]
        original_distances = row.drop('Clean_ID').values
        mapping = create_position_mapping(aligned_seq)

        for i, (orig_pos, aligned_pos, aa) in enumerate(mapping):
            dist = original_distances[i] if i < len(original_distances) else np.nan

            all_mappings.append({
                'Clean_ID': row['Clean_ID'],
                'FASTA_ID': matched_fasta_id,
                'Original_Pos': orig_pos + 1,
                'Aligned_Pos': aligned_pos + 1,
                'AA': aa,
                'Distance': dist
            })

    return pd.DataFrame(all_mappings)

# ============================================================================
# MAIN USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    
    # Load distance data (TSV format)
    # Format: ID in first column (no header), then positions
    distance_df = pd.read_csv('distance_table.tsv', sep='\t', index_col=0) # Adjust the distance table 
    
    # Reset index to make the ID a regular column called 'Clean_ID'
    distance_df = distance_df.reset_index()
    distance_df = distance_df.rename(columns={distance_df.columns[0]: 'Clean_ID'})
    
    print(f"Loaded {len(distance_df)} sequences with {len(distance_df.columns)-1} positions")
    print(f"First sequence ID: {distance_df.iloc[0]['Clean_ID']}")
    
    # Load aligned sequences from FASTA
    aligned_sequences = parse_aligned_sequences('aligned.fa') # Adjust the aligned fasta file
    print(f"Loaded {len(aligned_sequences)} aligned sequences")
    

    # Create aligned distance table
    # Automatically compute aligned length from FASTA
    aligned_length = max(len(seq) for seq in aligned_sequences.values())

    aligned_df = align_distances_to_positions(
        distance_df, 
        aligned_sequences, 
        aligned_length=aligned_length
    )

    
    # Save aligned distances as TSV
    aligned_df.to_csv('aligned_distances.tsv', sep='\t', index=False)
    print(f"\nAligned distance table saved: {aligned_df.shape}")
    print(f"\nFirst few rows:\n{aligned_df.head()}")
    

    # Create detailed mapping for inspection
    mapping_df = create_detailed_mapping_all_species(distance_df, aligned_sequences)
    mapping_df.to_csv('position_mapping_all_species.tsv', sep='\t', index=False)
    print(f"\nDetailed mapping saved: {mapping_df.shape}")
    
    # Summary statistics
    print("\n=== Summary ===")
    print(f"Total sequences processed: {len(aligned_df)}")
    print(f"Aligned length: {len([c for c in aligned_df.columns if c.startswith('P')])}")
    print(f"Average non-NA positions per sequence: {aligned_df.iloc[:, 1:].notna().sum(axis=1).mean():.1f}")
    
    # Show example of one sequence
    if len(aligned_df) > 0:
        print(f"\n=== Example: {aligned_df.iloc[0]['Clean_ID']} ===")
        seq_id = aligned_df.iloc[0]['Clean_ID']
        if seq_id in aligned_sequences:
            print(f"Aligned sequence: {aligned_sequences[seq_id][:60]}...")
            non_na_count = aligned_df.iloc[0, 1:].notna().sum()
            print(f"Non-NA positions: {non_na_count}")




def print_aligned_table_example(aligned_df, aligned_sequences, n_sequences=5):
    """
    Print aligned sequences along with their mapped distances
    for the first few sequences in the aligned_df.
    """
    for idx in range(min(n_sequences, len(aligned_df))):
        row = aligned_df.iloc[idx]
        distance_id = row['Clean_ID'].lower()

        # Find matched FASTA ID
        matched_fasta_id = next((fid for fid in aligned_sequences if distance_id.startswith(fid)), None)
        if not matched_fasta_id:
            print(f"Warning: {distance_id} not found in aligned sequences, skipping...")
            continue

        aligned_seq = aligned_sequences[matched_fasta_id]
        distances = row[1:].values  # all columns

        # Print header
        print(f"\n=== {row['Clean_ID']} ===")
        print(f"Distance table ID: {distance_id}")
        print(f"Matched FASTA ID: {matched_fasta_id}")

        # Print aligned sequence (first 60 aa for readability)
        seq_display = aligned_seq[:60] + '...' if len(aligned_seq) > 60 else aligned_seq
        print(f"Aligned sequence (first 60 aa):\n{seq_display}")

        # Print distances as aligned positions
        aligned_str = ''.join(['X' if not np.isnan(d) else '-' for d in distances])
        aligned_display = aligned_str[:60] + '...' if len(aligned_str) > 60 else aligned_str
        print(f"Distance mapping (X=has value, -=gap):\n{aligned_display}")

print_aligned_table_example(aligned_df, aligned_sequences, n_sequences=5)
