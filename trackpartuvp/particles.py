# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Detect particles for a given project
23-08-2022
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


def get_particles(directory, A, B, threshold):
    """
    Extracts the particles of a given set of images. Save the dataframe
    containing all particles.

    Parameters
    ----------
    directory : str
        Path to UVP6 project. It should contain a subdirectory 'raw' with  
        sequences of images.
    A, B : float
        UVP6 calibration constants.
    threshold : int
        Threshold for image segmentation (depends on UVP6).

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
    
    # Create a "results" directory
    create_dir('results')
    
    # Create a directory for particles data
    dir_results = os.path.join(directory, 'results')
    os.chdir(dir_results)
    create_dir('particles_data')
    
    # Create subdirectories in the results directory for each sequence
    os.chdir(dir_raw)
    subdirs = sorted(glob('*'))[1:-1]

    df_size_spectra = pd.DataFrame(
        columns = ['sequence', 'area', 'imgcount', 'nbr'])
    path_csv_all = os.path.join(dir_results, 'size_spectra.csv')
    
    print('--------------------')
    
    # Creating and saving dataframe with particles for each sequence
    
    for subdir in subdirs:
        
        print('Raw sequence: ', subdir)
        
        # Check if there is a zip file
        path_zip = os.path.join(dir_raw, subdir, 'images.zip')
        if os.path.exists(path_zip):
            archive = zipfile.ZipFile(path_zip)
            names = [i for i in archive.namelist() 
                     if not i.startswith('images')]
            for file in names:
                archive.extract(file, os.path.join(dir_raw, subdir))
        
            # Uncomment if you want to delete zip folder (to free space)
            #os.remove(path_zip)
            
        # Creating particles dataframe
        
        filename = "".join(['particles_' , subdir, '.csv'])
        path_csv1 = os.path.join(dir_results, 'particles_data', filename)
        
        if (Path(path_csv1).is_file() and Path(path_csv_all).is_file()):
            continue
            
        print('Getting image list...')
        img_names = sorted(
            filter(os.path.isfile, glob(os.path.join(
                dir_raw, subdir, '**', '**.png'), recursive = True)))
        
        # Run threshold on current subdir and images
        print('Looking for particles...')
        thrs = Thresholding(A = A, B = B, threshold = threshold)
        df_particles, df_areas = thrs.run(subdir, img_names)
        
        print('Saving dataframe particles')
        df_particles.to_csv(path_csv1, index = False)
        
        # Method append disappeared with pandas 2.0
        #df_size_spectra = df_size_spectra.append(
        #    df_areas, ignore_index = True, sort = True)
        # Instead use the method concat
        df_size_spectra = pd.concat([df_size_spectra, df_areas])
        
        print('--------------------')
        
    # Saving dataframe with particle size spectra
    # Only one for the whole project to avoid generating too much files
    if not Path(path_csv_all).is_file():
        os.chdir(dir_results)
        df_size_spectra.to_csv(path_csv_all, index = False)
    
    return

