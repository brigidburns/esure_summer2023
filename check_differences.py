#!/usr/bin/env python3

'''
This Python script is used to compare the original .frd files to the files processed using the ASPEN
software directly. 
'''

import subprocess
import tempfile
import os

# Specify the file path
processed_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/Kyle2008_Processed/'
dc_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/Kyle2008_DynamicCorrection/'

for processed_file in os.listdir(processed_path):
    for dc_file in os.listdir(dc_path):
        if processed_file == dc_file:
            # Define the command
            command = ['diff', '-w', '-q', processed_path+processed_file, dc_path+dc_file]

            # Execute the command in a subprocess
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # Wait for the command to finish, then capture the output
            output, error = process.communicate()

            # Print the output
            if output:
                print("Output:")
                print(output.decode())

            # Print errors, if any
            if error:
                print("Error:")
                print(error.decode())

'''
for processed_file in os.listdir(processed_path):
    processed = os.path.splitext(processed_file)[0]
    for frd_file in os.listdir(frd_path):
        frd = os.path.splitext(frd_file)[0]
        if processed == frd:
            print(processed_file)
            print(frd_file)
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_processed_file, tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_frd_file:
                with open(processed_file, 'r') as orig_processed_file, open(frd_file, 'r') as orig_frd_file:
                    # Skip the first 3 lines and write the rest to temporary files
                    temp_processed_file.writelines(orig_processed_file.readlines()[3:])
                    temp_frd_file.writelines(orig_frd_file.readlines()[3:])

                    temp_processed_file.flush()
                    temp_frd_file.flush()

                    # Specify the command with the temporary file paths
                    command = ['diff', '-w', temp_processed_file.name, temp_frd_file.name]
                    print(command)
                    # Execute the command in a subprocess
                    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
                    # Wait for the command to finish and capture the output
                    output, error = process.communicate()
                
                    # Print the output
                    if output:
                        print("Output:")
                        print(output.decode())

                    # Print any error messages
                    if error:
                        print("Error:")
                        print(error.decode())

                    break

processed_file = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/processed/test.txt'
frd_file = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/frd/test.txt'

with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_processed_file, tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_frd_file:
    with open(processed_file, 'r') as orig_processed_file, open(frd_file, 'r') as orig_frd_file:
        # Skip the first 3 lines and write the rest to temporary files
        temp_processed_file.writelines(orig_processed_file.readlines()[3:])
        temp_frd_file.writelines(orig_frd_file.readlines()[3:])

        temp_processed_file.flush()
        temp_frd_file.flush()

        # Specify the command with the temporary file paths
        command = ['diff', '-w', temp_processed_file.name, temp_frd_file.name]
        print(command)
        # Execute the command in a subprocess
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
        # Wait for the command to finish and capture the output
        output, error = process.communicate()
                
        # Print the output
        if output:
            print("Output:")
            print(output.decode())

        # Print any error messages
        if error:
            print("Error:")
            print(error.decode())

'''