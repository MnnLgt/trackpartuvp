# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Detect particles for a given project
2022
@author: Manon Laget
------------------------------
"""


import os
from glob import glob
from pathlib import Path
import pandas as pd
import zipfile
from .utils._utils import create_dir
from .funcparticles._thresholding import Thresholding
import shutil
import tarfile


def get_particles(source_directory, dest_directory, A, B, threshold,
                  platform_speed_path=None):
    """
    Extracts the particles of a given sequence of images. Save the dataframe
    containing all particles.

    Parameters
    ----------
    source_directory: str
        Path to UVP6 project. It should contain a subdirectory 'raw' with
        sequences of images.
    dest_directory: str
        Path to the folder where you want to save the tracking results.
    A, B: float
        UVP6 calibration constants.
    threshold: int
        Threshold for image segmentation (depends on UVP6).
    platform_speed_path: str
        Path to the file containing speed data of the platform (float or
        sediment trap).

    Returns
    -------
    None.

    """

    # 1. Setup

    # Move to the user-specified source directory, containing raw data
    os.chdir(source_directory)

    # Make sure that the "raw" directory exists
    dir_raw = os.path.join(source_directory, 'raw')
    if not os.path.exists(dir_raw):
        print('Raw directory does not exist')
        return

    # Move to destination directory
    os.chdir(dest_directory)

    # Create a "results" directory
    create_dir('results')

    # Create folder to store the temporary copy of images.zip
    create_dir('copy_raw')
    dir_copy_raw = os.path.join(dest_directory, 'copy_raw')

    # Create a directory for particles data
    dir_results = os.path.join(dest_directory, 'results')
    os.chdir(dir_results)
    create_dir('particles_data')

    # Create subdirectories in the results directory for each sequence
    os.chdir(dir_raw)
    subdirs = sorted(glob('*'))
    print(len(subdirs), "sequences")

    # Prepare a dataframe to store particle size distributions
    # df_size_spectra = pd.DataFrame(
    #     columns=['sequence', 'area', 'imgcount', 'nbr'])
    # path_csv_psd = os.path.join(dir_results, 'size_spectra.csv')

    # 2. Creating and saving dataframe with particles for each sequence

    print('--------------------')

    for subdir in subdirs:

        print('Raw sequence: ', subdir)

        # Prepare name of result file
        filename = "".join(['particles_', subdir, '.csv'])
        path_csv1 = os.path.join(dir_results, 'particles_data', filename)

        # Check if the particles dataframe already exists
        if (Path(path_csv1).is_file()):  # and Path(path_csv_psd).is_file()
            print('Particle dataframe already exists')
            continue

        # Check if there is a zip file
        path_zip = os.path.join(dir_raw, subdir, 'images.zip')
        path_tar = os.path.join(dir_raw, subdir)

        if os.path.exists(path_zip):
            archive = zipfile.ZipFile(path_zip)
            names = [i for i in archive.namelist()
                     if not i.startswith('images')]
            print('extract images from server to local')
            for file in names:
                archive.extract(file, os.path.join(dir_copy_raw, subdir))

            # Uncomment if you want to delete zip folder (to free space) DANGER
            # os.remove(path_zip)

        elif os.path.exists(path_tar):
            print('tar file found')
            with tarfile.open(path_tar) as tar:
                tar.extractall(path=os.path.join(dir_copy_raw, subdir))

        else:
            print('No zipped file found, moving to next sequence')
            continue

        # Creating particles dataframe

        print('Getting image list...')
        img_names = sorted(
            filter(os.path.isfile, glob(os.path.join(
                dir_copy_raw, subdir, '**', '**.png'), recursive=True)))

        # Run threshold on current subdir and images
        print('Looking for particles...')
        thrs = Thresholding(A=A, B=B, threshold=threshold)

        if len(img_names) > 5:
            print('Enough images to initiate tracking')

            df_particles = thrs.run(
                subdir, img_names, platform_speed_path)

            # df_particles, df_areas = thrs.run(
            #    subdir, img_names, platform_speed_path)

            print('Saving dataframe particles')
            path_csv1 = path_csv1.replace('.tar', '')
            df_particles.to_csv(path_csv1, index=False)

            # Method append disappeared with pandas 2.0
            # df_size_spectra = df_size_spectra.append(
            #    df_areas, ignore_index = True, sort = True)
            # Instead use the method concat
            # df_size_spectra = pd.concat([df_size_spectra, df_areas])

            print('--------------------')

        else:
            print('NOT enough images to initiate tracking')
            print('--------------------')

        # Delete the temporary copy raw folder
        print('Delete images from local')
        shutil.rmtree(os.path.join(dir_copy_raw, subdir))

    # Saving dataframe with particle size spectra
    # Only one for the whole project to avoid generating too many files
    # if not Path(path_csv_psd).is_file():
    #    os.chdir(dir_results)
    #    df_size_spectra.to_csv(path_csv_psd, index=False)

    return
