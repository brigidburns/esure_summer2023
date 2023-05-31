'''
This Python script is used to compare the original .frd files to the files processed using the ASPEN
software directly. 
'''

import subprocess
import tempfile

# Specify the files to be compared
file_a = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/test.txt'
file_b = '/escnfs/home/bburns4/richter-lab/kyle2008/esure_summer2023/test2.txt'

# Create temporary files
with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file_a, tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file_b:
    with open(file_a, 'r') as orig_file_a, open(file_b, 'r') as orig_file_b:
        # Skip the first 3 lines and write the rest to temporary files
        temp_file_a.writelines(orig_file_a.readlines()[3:])
        temp_file_b.writelines(orig_file_b.readlines()[3:])

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