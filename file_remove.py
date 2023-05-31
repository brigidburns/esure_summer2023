'''
This Python script is used to remove the raw files without matching .frd files. 
The script runs through each raw file, comparing it to each of the available .frd files. 
If the sonde numbers match, the file is left. If there is no matching file, the raw file is placed 
in a list of files to remove. After all files have been compared, all of the files that need to be removed
are removed. 
'''
import os

# Specify the path to each directory
raw_path = '/escnfs/home/bburns4/richter-lab/kyle2008/test/Kyle2008_AVAPS/'
frd_path = '/escnfs/home/bburns4/richter-lab/kyle2008/test/Kyle2008/'

# Define an empty list to place the files that need to be removed
remove_files = []

# Loop through each of the raw files 
for raw_file in os.listdir(raw_path):
    remove = True
    # Separate the raw file name and its extension
    raw = os.path.splitext(raw_file)[0]
    # Loop through each of the .frd files 
    for frd_file in os.listdir(frd_path):
        # Separate the .frd file name and its extension
        frd = os.path.splitext(frd_file)[0]
        # If there is a match, do not remove
        if raw == frd:
            remove = False
            break
    # If there is no match, add to the list of files to be removed
    if remove:
        remove_files.append(raw_file)

# For each file in the list of files to remove
for file_name in remove_files:
    # Specify the path of the file to remove
    file_to_remove = '/escnfs/home/bburns4/richter-lab/kyle2008/test/Kyle2008_AVAPS/'+file_name
    # If the file exists
    if os.path.exists(file_to_remove):
        print(f'Removing {file_name}...')
        # Remove it
        os.remove(file_to_remove)