#!/usr/bin/env python3

import os

# Specify the path to each directory
processed_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/Kyle2008_DynamicCorrection/'
dc_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/Kyle2008_Processed/'

# Define an empty list to place the files that need to be removed
unmatched_files = []

# Loop through each of the raw files 
for processed_file in os.listdir(processed_path):
    unmatched = True
    # Separate the raw file name and its extension
    processed = os.path.splitext(processed_file)[0]
    # Loop through each of the .frd files 
    for dc_file in os.listdir(dc_path):
        # Separate the .frd file name and its extension
        dc = os.path.splitext(dc_file)[0]
        # If there is a match, do not remove
        if processed == dc:
            unmatched = False
            break
    # If there is no match, add to the list of files to be removed
    if unmatched:
        unmatched_files.append(dc_file)

# Print the unmatched files
if (len(unmatched_files) > 0):
    for file_name in unmatched_files:
        print(file_name)
else:
    print("No unmatched files found.")