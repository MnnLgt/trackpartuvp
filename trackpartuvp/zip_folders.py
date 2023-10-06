# -*- coding: utf-8 -*-
"""
Created on Thu Sep  7 08:18:48 2023

@author: Manon Laget
"""


import os
from glob import glob
import shutil


directory = ''

# Move to the user specified directory
os.chdir(directory)

# Make sure that the "raw" directory exists
dir_raw = os.path.join(directory, 'raw')
if not os.path.exists(dir_raw):
    print('Raw directory does not exist')

os.chdir(dir_raw)
subdirs = sorted(glob('*'))

print('--------------------')

# Creating and saving dataframe with particles for each sequence

for subdir in subdirs:
    
    print('Raw sequence: ', subdir)
    
    # Check if there is a zip file
    path_zip = os.path.join(dir_raw, subdir, 'images.zip')
    if not os.path.exists(path_zip):
        shutil.make_archive(path_zip, 'zip', os.path.join(dir_raw, subdir))
    
