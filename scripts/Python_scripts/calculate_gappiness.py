#!/usr/bin/env python3
### CHATGPT

from Bio import AlignIO
from pathlib import Path
import sys

def calculate_gap_percentage(alignment):
    total_bases = sum(len(seq) for seq in alignment)
    gap_bases = sum(seq.count('-') for seq in alignment)
    return (gap_bases / total_bases) * 100

def process_alignment_files(folder_path):
    with open("calculate_gappiness_results.txt", 'w') as output:
        folder = Path(folder_path)
        
        # Iterate over each FASTA file in the folder
        for file_path in folder.glob('*.fasta'):
            print(f"Processing file: {file_path.name}")
            
            # Read the alignment from the file
            alignment = AlignIO.read(file_path, 'fasta')
            
            # Calculate the percentage of gap sites
            gap_percentage = calculate_gap_percentage(alignment)
            
            # Print the result
            print(f"Percentage of gap sites: {gap_percentage:.2f}%\n")

            # Write the results to the file
            output.write(f"{file_path.name}\t{gap_percentage:.2f}\n")

if __name__ == "__main__":
    folder_path = str(sys.argv[1])
    
    # Process the alignment files in the specified folder
    process_alignment_files(folder_path)
