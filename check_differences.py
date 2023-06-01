'''
This Python script is used to compare the original .frd files to the files processed using the ASPEN
software directly. 
'''

import subprocess
import tempfile
import os

# Specify the files to be compared
file_a = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/test.txt'
file_b = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/test2.txt'

# Designate the file paths
# Note to self: need to probably add the data to the python branch 
# but can I access it if I'm still on this branch? Idk something to look into 
processed_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023'
frd_path = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023'

for processed_file in os.listdir(processed_path):
    processed = os.splitext(processed_file)[0]
    for frd_file in os.listdir(frd_path):
        frd = os.splitext(frd_file)[0]
        if processed == frd:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_processed_file, tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_frd_file:
                with open(processed_file, 'r') as orig_processed_file, open(frd_file, 'r') as orig_frd_file:
                # Skip the first 3 lines and write the rest to temporary files
                    temp_processed_file.writelines(orig_processed_file.readlines()[20:])
                    temp_frd_file.writelines(orig_frd_file.readlines()[20:])

                temp_file_a.flush()
                temp_file_b.flush()

                # Specify the command with the temporary file paths
                command = ['diff', '-w', temp_file_a.name, temp_file_b.name]

                # Execute the command in a subprocess
                process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

                # Wait for the command to finish and capture the output
                output, error = process.communicate()

                # Print the output
                if output:
                    print(output.decode())

                # Print any error messages
                if error:
                    print(error.decode())