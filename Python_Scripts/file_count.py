#!/usr/bin/env python3

import os
import subprocess
import shlex

path_base = '../All_Storms/'

# Prompt the user to enter the hurricane
hurricane = input('Enter the hurricane name and year, in the form NameYear (ex: Katrina2005): ')

# Specify the path to each directory
raw_path = path_base + hurricane + '/' + hurricane + '_Raw/'
frd_path = path_base + hurricane + '/' + hurricane + '_FRD/'
processed_path = path_base + hurricane + '/' + hurricane + '_Processed/'
dc_path = path_base + hurricane + '/' + hurricane + '_DynamicCorrection/'

if not os.path.exists(raw_path) or not os.path.exists(frd_path):
    print('Error: Invalid path.')
    exit(1)

# Count the number of raw and .frd files
num_raw_files = sum(len(files) for _, _, files in os.walk(raw_path))
num_frd_files = sum(len(files) for _, _, files in os.walk(frd_path))

# Print the number of files
print(f'Number of raw files: {num_raw_files}')
print(f'Number of FRD files: {num_frd_files}')

if os.path.exists(processed_path):
    num_processed_files = sum(len(files) for _, _, files in os.walk(processed_path))
    print(f'Number of processed files: {num_processed_files}')

if os.path.exists(dc_path):
    num_dc_files = sum(len(files) for _, _, files in os.walk(dc_path))
    print(f'Number of DC files: {num_dc_files}')