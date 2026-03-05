#!/usr/bin/env python3

# Script generated using Claude 3.7 Sonnet

import re
import csv
import statistics
import os
import sys

def parse_iqtree_file(file_path):
    """
    Parse an IQ-TREE output file and extract alignment statistics.
    Returns a list of dictionaries containing the parsed data.
    """
    # Initialize variables
    records = []
    capture = False
    sequence_section_found = False
    lines_after_sequence = []
    header_found = False
    
    print(f"Attempting to open file: {file_path}")
    if not os.path.exists(file_path):
        print(f"Error: File {file_path} does not exist")
        return []
    
    # Open and read the file
    with open(file_path, 'r') as file:
        line_num = 0
        for line in file:
            line_num += 1
            
            # Check if we've reached the sequence alignment section
            if "SEQUENCE ALIGNMENT" in line:
                print(f"Found SEQUENCE ALIGNMENT section at line {line_num}")
                capture = True
                sequence_section_found = True
                continue
            
            # Collect lines after sequence alignment for debugging
            if sequence_section_found and len(lines_after_sequence) < 10:
                lines_after_sequence.append((line_num, line.strip()))
            
            # Look for the header line with "ID", "Type", "Seq", etc.
            if capture and "ID" in line and "Type" in line and "Seq" in line:
                header_found = True
                print(f"Found header line at line {line_num}: {line.strip()}")
                continue
            
            # Skip separator line after header (usually dashes)
            if capture and header_found and line.strip().startswith("-"):
                continue
                
            # Stop capturing if we reach a blank line after the table
            if capture and header_found and line.strip() == "":
                break
                
            # Capture and parse the table rows
            if capture and header_found and line.strip():
                # Parse the line using regex to handle multiple spaces
                parts = re.split(r'\s+', line.strip())
                
                # Skip if the line doesn't match our expected format
                if len(parts) < 9:
                    continue
                
                # Try to parse the ID part to check if it's a valid row
                try:
                    id_part = int(parts[0])
                    
                    # Create a record from the parsed data
                    record = {
                        'ID': id_part,
                        'Type': parts[1],
                        'Seq': int(parts[2]),
                        'Site': int(parts[3]),
                        'Unique': int(parts[4]),
                        'Infor': int(parts[5]),
                        'Invar': int(parts[6]),
                        'Const': int(parts[7]),
                        'Name': parts[8]
                    }
                    records.append(record)
                    
                except ValueError:
                    continue
    
    print("\nFirst 10 lines after SEQUENCE ALIGNMENT section:")
    for idx, (line_num, line) in enumerate(lines_after_sequence):
        print(f"Line {line_num}: '{line}'")
    
    if not sequence_section_found:
        print("Warning: No 'SEQUENCE ALIGNMENT' section found in the file")
    elif not header_found:
        print("Warning: Found 'SEQUENCE ALIGNMENT' section but couldn't find header line with ID, Type, Seq, etc.")
    elif not records:
        print("Warning: Found header but couldn't parse any records")
        
    return records

def write_to_csv(records, output_file):
    """
    Write the parsed records to a CSV file.
    """
    fieldnames = ['ID', 'Type', 'Seq', 'Site', 'Unique', 'Infor', 'Invar', 'Const', 'Name']
    
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

def calculate_statistics(records, output_file):
    """
    Calculate statistics for specific columns and write to a text file.
    """
    # Extract the columns of interest
    seq_values = [record['Seq'] for record in records]
    site_values = [record['Site'] for record in records]
    infor_values = [record['Infor'] for record in records]
    
    # Calculate means
    seq_mean = statistics.mean(seq_values)
    site_mean = statistics.mean(site_values)
    infor_mean = statistics.mean(infor_values)
    
    # Write statistics to file
    with open(output_file, 'w') as f:
        f.write("Statistics Summary for IQ-TREE Alignment\n")
        f.write("=======================================\n\n")
        f.write(f"Mean Seq: {seq_mean:.2f}\n")
        f.write(f"Mean Site: {site_mean:.2f}\n")
        f.write(f"Mean Infor: {infor_mean:.2f}\n")

def main():
    # You can adjust the input file name here
    input_file = str(sys.argv[1])  # Input IQ-TREE file
    csv_output = "alignment_data.csv"  # Output CSV file
    stats_output = "alignment_statistics.txt"  # Output statistics file
    
    # Display current working directory to help with file path issues
    print(f"Current working directory: {os.getcwd()}")
    
    try:
        # Parse the IQ-TREE file
        print(f"Parsing file: {input_file}")
        records = parse_iqtree_file(input_file)
        
        if not records:
            print("No alignment data found in the input file.")
            
            # Additional debug option: dump the entire file content
            print("\nWould you like to see the entire content of the file? (y/n)")
            response = input().strip().lower()
            if response == 'y':
                print("\nFile content:")
                with open(input_file, 'r') as f:
                    print(f.read())
            return
        
        print(f"Found {len(records)} records")
        
        # Write data to CSV
        write_to_csv(records, csv_output)
        print(f"Successfully wrote alignment data to {csv_output}")
        
        # Calculate and write statistics
        calculate_statistics(records, stats_output)
        print(f"Successfully wrote statistics to {stats_output}")
        
    except Exception as e:
        print(f"An error occurred: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()