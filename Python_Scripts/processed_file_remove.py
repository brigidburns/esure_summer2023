#!/usr/bin/env python3

'''
This Python script is used to remove the processed frd files without matching raw/frd files. 
The script runs through each processed .frd file, comparing it to each of the available raw files. 
If the sonde numbers match, the file is left. If there is no matching file, the processed .frd file is placed 
in a list of files to remove. After all files have been compared, all of the files that need to be removed
are removed. 
Use this script when the number of processed .frd files does not match the number of .frd and raw files. 
'''
import os

path_base = '../All_Storms/'

# Prompt the user for the hurricane
hurricane = input("Enter the hurricane name and year, in the form NameYear (ex: Katrina2005): ")

# Specify the path to each directory
raw_path = path_base + hurricane + '/' + hurricane + '_Raw/'
processed_path = path_base + hurricane + '/' + hurricane + '_Processed/'

if not os.path.exists(raw_path) or not os.path.exists(processed_path): 
    print("Error: Invalid path.")
    exit(1)

# Define an empty list to place the files that need to be removed
remove_files = []

# Loop through each of the raw files 
for processed_file in os.listdir(processed_path):
    remove = True
    # Separate the raw file name and its extension
    processed = os.path.splitext(processed_file)[0][:-2]
    # Loop through each of the .frd files 
    for raw_file in os.listdir(raw_path):
        # Separate the .frd file name and its extension
        raw = os.path.splitext(raw_file)[0]
        # If there is a match, do not remove
        if raw == processed:
            remove = False
            break
    # If there is no match, add to the list of files to be removed
    if remove:
        remove_files.append(processed_file)

if (len(remove_files) == 0):
    print("No files to remove.")
    exit(0)

# For each file in the list of files to remove
for file_name in remove_files:
    # Specify the path of the file to remove
    file_to_remove = processed_path + file_name
    # If the file exists
    if os.path.exists(file_to_remove):
        print(f'Removing {file_name}...')
        # Remove it
        os.remove(file_to_remove)