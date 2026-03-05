## Script written by Claude 3.7 Sonnet

## Takes a fasta alignment and a txt file with samples grouped, outputs an alignment with a consensus sequence of the given groups. Ungrouped samples are copied straight to the output.

# Example txt file: Flindersia_brassii = "Flindersia_brassii_Costion_C_1763","Flindersia_brassii_MJB1889" per line
# To run: python condense_alignments.py Flindersia_confident_monophyletics.txt 4932_1.ortho.selected_stripped.aln.trimmed.fasta --output flindersia_consensus.fasta


from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import re
import os

def parse_groups_file(file_path):
    """Parse the groups file to extract sample groupings."""
    groups = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or "=" not in line:
                continue
                
            parts = line.split("=")
            if len(parts) != 2:
                continue
                
            group_name = parts[0].strip()
            # Extract sample names from the quoted strings
            samples = re.findall(r'"([^"]+)"', parts[1])
            groups[group_name] = samples
    
    return groups

def get_sequences_for_group(alignment, sample_names):
    """Extract sequences for specified sample names from the alignment."""
    group_sequences = []
    matched_ids = []
    
    for record in alignment:
        record_id = record.id
        
        # Check if this record matches any sample in the group
        for sample in sample_names:
            if sample in record_id:
                group_sequences.append(record)
                matched_ids.append(record_id)
                break
    
    return group_sequences, matched_ids

def calculate_consensus(sequences):
    """Calculate consensus sequence from a list of sequence records."""
    if not sequences:
        return None
    
    # Create a temporary alignment from the sequences
    from Bio.Align import MultipleSeqAlignment
    temp_alignment = MultipleSeqAlignment(sequences)
    
    # Calculate the summary info
    summary_info = AlignInfo.SummaryInfo(temp_alignment)
    
    # Get the consensus with a threshold of 0.7 (customizable)
    consensus = summary_info.dumb_consensus(threshold=0.7, ambiguous='-')
    
    return consensus

def main():
    parser = argparse.ArgumentParser(description='Generate consensus sequences for groups of samples.')
    parser.add_argument('groups_file', help='File containing sample groupings')
    parser.add_argument('alignment_file', help='Alignment file in supported format (FASTA, Clustal, etc.)')
    parser.add_argument('--format', default='fasta', help='Alignment file format (default: fasta)')
    parser.add_argument('--output', default='consensus_sequences.fasta', help='Output file for consensus sequences')
    parser.add_argument('--threshold', type=float, default=0.7, help='Consensus threshold (0.0-1.0)')
    
    args = parser.parse_args()
    
    # Parse the groups file
    groups = parse_groups_file(args.groups_file)
    print(f"Found {len(groups)} groups in the file")
    
    # Read the alignment
    try:
        alignment = AlignIO.read(args.alignment_file, args.format)
        print(f"Loaded alignment with {len(alignment)} sequences, each of length {alignment.get_alignment_length()}")
    except Exception as e:
        print(f"Error reading alignment: {e}")
        return
    
    # Calculate consensus for each group
    output_records = []
    all_grouped_ids = set()  # Track all sequence IDs that belong to groups
    
    for group_name, sample_names in groups.items():
        print(f"Processing group: {group_name} with {len(sample_names)} samples")
        
        # Get sequences for this group and their IDs
        group_sequences, matched_ids = get_sequences_for_group(alignment, sample_names)
        print(f"  Found {len(group_sequences)} matching sequences in alignment")
        
        # Add matched IDs to our tracking set
        all_grouped_ids.update(matched_ids)
        
        if group_sequences:
            # Calculate consensus
            consensus = calculate_consensus(group_sequences)
            
            # Create a SeqRecord for the consensus
            consensus_record = SeqRecord(
                seq=consensus,
                id=f"consensus_{group_name}",
                description=f"Consensus sequence for {group_name} (from {len(group_sequences)} sequences)"
            )
            
            output_records.append(consensus_record)
    
    # Find and add sequences not in any group
    ungrouped_count = 0
    
    for record in alignment:
        if record.id not in all_grouped_ids:
            output_records.append(record)
            ungrouped_count += 1
    
    print(f"Found {ungrouped_count} sequences not assigned to any group")
    
    # Write all sequences to output file
    with open(args.output, 'w') as output_handle:
        SeqIO.write(output_records, output_handle, 'fasta')
    
    print(f"Wrote {len(output_records)} sequences to {args.output} ({len(output_records) - ungrouped_count} consensus sequences and {ungrouped_count} ungrouped sequences)")

if __name__ == "__main__":
    main()