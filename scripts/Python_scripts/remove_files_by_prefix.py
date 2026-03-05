#!/usr/bin/env python3

# Author: Harvey Orel
# Date 21 Feb 2023

# Removes files from a specified directory that don't have the prefixes specified in a prefixes.txt file

import os
import sys

def read_prefixes(file_path):
    with open(file_path, 'r') as f:
        prefixes = set([line.strip() for line in f])
    return prefixes

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <directory_path> <prefixes_file>")
        sys.exit(1)

    directory_path = sys.argv[1]
    prefixes_file = sys.argv[2]

    # Read prefixes to keep from file
    prefixes_to_keep = read_prefixes(prefixes_file)

    # Iterate over files in directory
    for file_name in os.listdir(directory_path):
        # Check if file name starts with a prefix not in the list to keep
        if file_name.startswith(tuple(prefixes_to_keep)):
            continue

        # Remove file if prefix not in the list to keep
        file_path = os.path.join(directory_path, file_name)
        os.remove(file_path)
        print(f"Removed {file_path}")
