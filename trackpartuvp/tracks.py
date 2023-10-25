# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Pipeline to track detected particles and save tracks
23-08-2022
@author: Manon Laget
------------------------------
"""


import os
from glob import glob
from pathlib import Path
from shutil import copy2, make_archive
import pandas as pd
from trackpartuvp.utils._utils import parse_datetime, create_dir
from trackpartuvp.functracks._inclino import add_inclination
from trackpartuvp.functracks._tracking import Tracking
from trackpartuvp.functracks._get_tracks import GetTracks
from trackpartuvp.functracks._plots import Plots
from trackpartuvp.functracks._plot_ecotaxa import plot_ecotaxa
from trackpartuvp.functracks._plots_global import PlotsGlobal


def get_tracks(directory, deployment, depth, path_inclino = '', axis = None,
               sign = None, save_zip = True, angle = 100):
    """
    Pipeline to track particles for a given project.

    Parameters
    ----------
    directory : str
        Path to UVP6 project. It should contain a subdirectory 'raw' with  
        sequences of images.
    deployment : str
        Deployment name for this project.
    depth : str
        Deployment depth.
    path_inclino : str, optional
        Path to inclinometer data file, to correct the tracks. The default is
        None.
    save_zip : bool, optional
        Whether to save the results for each sequence as a .zip file, to
        facilitate export. The default is True.
    angle : int, optional
        Set angle for particle tracking.
    
    Returns
    -------
    None.

    """
    
    os.chdir(directory)
    
    # Check if deployment and depth are strings
    
    if not type(deployment) is str:
        raise TypeError('Deployment name should be of string type')
    
    if not type(depth) is str:
        raise TypeError('Depth should be of string type')
    
    # Make sure that the "particles_data" directory exists in "results"
    dir_results = os.path.join(directory, 'results')
    if os.path.exists(dir_results):
        os.chdir(dir_results)
        if not os.path.exists('particles_data'):
            print('particles_data directory does not exist')
            return
    
    dir_particles = os.path.join(dir_results, 'particles_data')
    
    # Get subdirectories names (=sequence names)
    os.chdir(os.path.join(directory, 'raw'))
    subdirs = sorted(glob('*'))[1:-1]
    
    # Create directory to store tracks
    os.chdir(dir_results)
    create_dir('detected_tracks')
    dir_tracks = os.path.join(dir_results, 'detected_tracks')
    
    # Create directory to store plots for EcoTaxa
    os.chdir(dir_tracks)
    create_dir('plots_ecotaxa')
    dir_ecotaxa = os.path.join(dir_tracks, 'plots_ecotaxa')
    
    # To sort them by orientation
    os.chdir(dir_ecotaxa)
    for orientation in ['asc', 'desc', 'mix']:
        create_dir(orientation)
    
    # Create directory to store the 1st vignette of a track (to be displayed
    # in further analyses)
    os.chdir(dir_tracks)
    create_dir('img_to_show')
    
    # Create a sub directory for each sequecne
    for subdir in subdirs:
        create_dir(subdir)
    
    # For plots
    lengths, angles, speeds, orientations = ([] for i in range(4))
    
    # For Ecotaxa file
    list_ecotaxa = [{'img_file_name': '[t]',
                     'img_rank': '[t]',
                     'object_id': '[t]',
                     'object_esd': '[f]',
                     'object_orientation': '[t]',
                     'object_sequence': '[t]',
                     'sample_id': '[t]',
                     'sample_deployment': '[t]',
                     'sample_depth': '[t]'}]

    print('--------------------')
    
    # Tracking particles and saving tracks for each sequence
    for subdir in subdirs:
        
        # Dataframe containing particles
        print('Sequence: ', subdir)
        filename = "".join(('particles_', subdir, '.csv'))
        path_part = os.path.join(dir_particles, filename)
        
        # If the file doesn't exist, let's continue
        if not Path(path_part).is_file():
            continue
            
        # Read paticle file
        df_particles = pd.read_csv(path_part)
        
        # If no image (no particle), let's continue
        if not 'img_name' in df_particles.columns:
            continue
        
        #Make sure that datetime has the good format
        df_particles['datetime'] = [
            parse_datetime(dt) for dt in df_particles['img_name'].tolist()]
        df_particles.sort_values(by = 'datetime')
        
        # Get correction for angle for this sequence
        df_particles['roll'] = [None]*len(df_particles)
        if path_inclino:
            df_particles = add_inclination(
                path_inclino, df_particles, axis, sign)
        
        res_subdir = os.path.join(dir_tracks, subdir)
        os.chdir(res_subdir)
        
        # Check if a file already exists for tracks
        if glob(os.path.join(res_subdir, '**.csv'), recursive = True):
            print('csv files found: moving to the next sequence')
            continue

        # Looking for tracks
        print('Tracking particles...')
        tracker = Tracking(df_particles = df_particles, angle = angle)
        list_tracks = tracker.run()
        
        if not list_tracks:
            print('List of tracks is empty')
            os.chdir(dir_tracks)
            continue
        
        # Saving tracks
        print('Saving tracks...')
        
        # Appending to list of lengths, angles, orientations
        for track in list_tracks:
            lengths.append(track.length)
            angles.append(track.mean_angle)
            speeds.append(track.vertical_speed)
            orientations.append(track.orientation)
        
        GetTracks(list_tracks = list_tracks, subdir = subdir, 
                  deployment = deployment, depth = depth).run()
        
        # Creating subdirectories for the tracks and saving plots
        Plots(list_tracks, df_particles, res_subdir).run()
        
        # Saving plots_ecotaxa
        os.chdir(dir_ecotaxa)
        dics = plot_ecotaxa(list_tracks = list_tracks, subdir = subdir, 
                            deployment = deployment, depth = depth)
        list_ecotaxa.extend(dics)
        
        # Save the 1st vignette in img_to_show
        for track in list_tracks:
            # Filename in its original folder
            filename_src = "".join(
                (str(track.img_names[0]), '-', str(track.id), '.png'))
            src = os.path.join(
                dir_tracks, subdir, str(track.id), filename_src)
            # Filename in destination folder (img_to_show)
            filename_dst = "".join((
                deployment, depth, '_', str(track.img_names[0]), '-', 
                str(track.id), '.png'))
            dst = os.path.join(dir_tracks, 'img_to_show', filename_dst)
            # Paste in destination folder
            copy2(src, dst)
        
        # Make a zip file with results
        if save_zip:
            make_archive(res_subdir, 'zip', res_subdir)

        print('--------------------')

        os.chdir(dir_tracks)
                
    # Save global plots
    os.chdir(dir_tracks)
    PlotsGlobal(lengths = lengths, angles = angles, speeds = speeds, 
                orientations = orientations, deployment = deployment, 
                depth = depth).run()
    
    # Save tsv
    os.chdir(dir_ecotaxa)
    tsv_name = "".join(('ecotaxa_', deployment, depth, '.tsv'))
    pd.DataFrame(list_ecotaxa).to_csv(tsv_name, sep = "\t", index = False)

    # Create a global file
    bool_file = False
    ct_subdir = 0 
    
    # Read all summary files and append them to global file
    while ct_subdir<len(subdirs):
        subdir = subdirs[ct_subdir]
        res_subdir = os.path.join(dir_tracks, subdir)
        filename = os.path.join(
            res_subdir, "".join(
                ('tracks_summary_', deployment, '_', depth, '_', subdir, 
                 '.tsv')))
        ct_subdir += 1
        if Path(filename).is_file():
            if not bool_file:
                df_global = pd.read_csv(filename, sep = "\t")
                bool_file = True
            else:
                df_summary = pd.read_csv(filename, sep = "\t")
                # Method append disappeared with pandas 2.0
                #df_global = df_global.append(df_summary, ignore_index = True)
                # Instead use the method concat
                df_global = pd.concat([df_global, df_summary])
                
    # Save global file
    filename = "".join(('tracks_summary_', deployment, '_', depth, '.tsv'))
    filename = os.path.join(dir_tracks, filename)
    df_global.to_csv(filename, sep = "\t", index = False)

    return


