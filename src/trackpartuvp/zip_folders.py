# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Zip images for each sequence of a project
2022
@author: Manon Laget
------------------------------
"""


import os
from glob import glob
import shutil


def zip_image_folders(directory):
    """
    Zip images for each sequence of a project.

    Parameters
    ----------
    directory: str
        Path to UVP6 project.

    Returns
    -------
    None.

    """

    # Move to the user specified directory
    os.chdir(directory)

    # Make sure that the "raw" directory exists
    dir_raw = os.path.join(directory, 'raw')
    if not os.path.exists(dir_raw):
        print('Raw directory does not exist')
        return

    os.chdir(dir_raw)
    subdirs = sorted(glob('*'))

    print('--------------------')

    for subdir in subdirs:

        print('Raw sequence: ', subdir)

        # Check if there is a zip file
        path_zip = os.path.join(dir_raw, subdir, 'images.zip')
        if not os.path.exists(path_zip):
            shutil.make_archive(path_zip, 'zip', os.path.join(dir_raw, subdir))
            print('Zipped')

    return
