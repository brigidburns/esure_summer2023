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

if not os.path.exists(raw_path) or not os.path.exists(frd_path):
    print('Error: Invalid path.')
    exit(1)

# Define the commands to execute
raw_ls_command = 'ls {}'.format(shlex.quote(raw_path))
raw_wc_command = 'wc -l'
frd_ls_command = 'ls {}'.format(shlex.quote(frd_path))
frd_wc_command = 'wc -l'

# Execute the commands
raw_ls_process = subprocess.Popen(raw_ls_command, shell=True, stdout=subprocess.PIPE)
raw_wc_process = subprocess.Popen(raw_wc_command, shell=True, stdin=raw_ls_process.stdout, stdout=subprocess.PIPE)
frd_ls_process = subprocess.Popen(frd_ls_command, shell=True, stdout=subprocess.PIPE)
frd_wc_process = subprocess.Popen(frd_wc_command, shell=True, stdin=frd_ls_process.stdout, stdout=subprocess.PIPE)

# Capture the output of the command
raw_output, raw_error = raw_wc_process.communicate()
frd_output, frd_error = frd_wc_process.communicate()

# Print the output of the command 
print('Number of raw files:')
print(raw_output.decode().strip())
print('Number of FRD files:')
print(frd_output.decode().strip())