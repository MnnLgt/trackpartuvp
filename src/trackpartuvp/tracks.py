# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Pipeline to track detected particles and save tracks
2022
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
# from trackpartuvp.functracks._plots import Plots
# from trackpartuvp.functracks._plot_ecotaxa import plot_ecotaxa
from trackpartuvp.functracks._plot_ecotaxa_from_df import plot_ecotaxa_from_dataframe
# from trackpartuvp.functracks._plots_global import PlotsGlobal


def get_tracks(source_directory, dest_directory, deployment, depth,
               platform_speed_path='', angle=90, inclino_params=None,
               save_zip=False):
    """
    Pipeline to track particles for a given project.

    Parameters
    ----------
    source_directory: str
        Path to UVP6 project. It should contain a subdirectory 'raw' with
        sequences of images.
    dest_directory: str
        Path to the folder where you want to save the tracking results.
    deployment: str
        Deployment name for this project.
    depth: str
        Deployment depth.
    platform_speed_path: str
        Path to file containing speeds of the platform (float, sediment trap).
    angle: int, optional
        Set angle for particle tracking.
    inclino_params: str, optional
        Path to dataframe containing parameters to use to prepare inclino data.
    save_zip: bool, optional
        Whether to save the results for each sequence as a .zip file, to
        facilitate export. The default is True.

    Returns
    -------
    None.

    """

    # 1. Setup

    os.chdir(source_directory)

    # Check if deployment and depth are strings

    if not type(deployment) is str:
        raise TypeError('Deployment name should be of string type')

    if not type(depth) is str:
        raise TypeError('Depth should be of string type')

    # Make sure that the "particles_data" directory exists in "results"
    dir_results = os.path.join(dest_directory, 'results')
    if os.path.exists(dir_results):
        os.chdir(dir_results)
        if not os.path.exists('particles_data'):
            print('particles_data directory does not exist')
            return

    # Path to particle data
    dir_particles = os.path.join(dir_results, 'particles_data')

    # Get subdirectories names (=sequence names)
    os.chdir(os.path.join(source_directory, 'raw'))
    subdirs = sorted(glob('*'))

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

    # Create a sub directory for each sequence
    for subdir in subdirs:
        subdir = subdir.replace('.tar', '')
        create_dir(subdir)

    # For plots
    lengths, angles, speeds, orientations = ([] for i in range(4))

    # For header of Ecotaxa file
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

    # 2.  Tracking particles and saving tracks for each sequence

    for subdir in subdirs:

        subdir = subdir.replace('.tar', '')

        # Dataframe containing particles
        print('Sequence: ', subdir)
        filename = "".join(('particles_', subdir, '.csv'))
        path_part = os.path.join(dir_particles, filename)

        # If the file doesn't exist, let's continue
        if not Path(path_part).is_file():
            continue

        # Read paticle file
        try:
            df_particles = pd.read_csv(path_part)

        except:
            print(f'df_particles {subdir} is empty')
            continue

        # If no image (no particle), let's continue
        # if not 'img_name' in df_particles.columns:
        if len(df_particles) == 0:
            continue

        # Make sure that datetime has the good format
        df_particles['datetime'] = [
            parse_datetime(dt) for dt in df_particles['img_name'].tolist()]
        df_particles.sort_values(by='datetime')

        # Get correction for angle for this sequence
        df_particles['roll'] = [None]*len(df_particles)
        if inclino_params is None:
            print('')
        else:
            df_particles = add_inclination(inclino_params, df_particles)

        res_subdir = os.path.join(dir_tracks, subdir)
        os.chdir(res_subdir)

        # Check if a file already exists for tracks
        if glob(os.path.join(res_subdir, '**.csv'), recursive=True):
            print('csv files found: moving to the next sequence')
            continue

        # Looking for tracks
        print('Tracking particles...')
        tracker = Tracking(df_particles=df_particles, angle=angle,
                           platform_speed_path=platform_speed_path)  # AA
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

        GetTracks(
            list_tracks=list_tracks, subdir=subdir, deployment=deployment,
            depth=depth).run()

        # Creating subdirectories for the tracks and saving plots
        # Plots(list_tracks, df_particles, res_subdir).run()

        # Saving plots_ecotaxa
        # os.chdir(dir_ecotaxa)

        # dics = plot_ecotaxa(
        #     list_tracks=list_tracks, subdir=subdir, deployment=deployment,
        #     depth=depth, save_plots=False)
        # list_ecotaxa.extend(dics)

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
            try:
                copy2(src, dst)
            except:
                print('no such directory')

        # Make a zip file with results
        if save_zip:
            make_archive(res_subdir, 'zip', res_subdir)

        print('--------------------')

        os.chdir(dir_tracks)

    # Save global plots
    os.chdir(dir_tracks)
    # PlotsGlobal(lengths=lengths, angles=angles, speeds=speeds,
    #             orientations=orientations, deployment=deployment,
    #             depth=depth).run()

    # Save tsv
    os.chdir(dir_ecotaxa)
    tsv_name = "".join(('ecotaxa_', deployment, depth, '.tsv'))
    pd.DataFrame(list_ecotaxa).to_csv(tsv_name, sep="\t", index=False)

    # Create a global file
    bool_file_global = False
    bool_file_all = False  # AA
    ct_subdir = 0

    # Initialize global and all dataframes (AA)
    df_global = pd.DataFrame()
    df_all = pd.DataFrame()

    # Read all summary files and append them to global file
    while ct_subdir < len(subdirs):

        subdir = subdirs[ct_subdir]
        subdir = subdir.replace('.tar', '')
        res_subdir = os.path.join(dir_tracks, subdir)
        res_subdir = res_subdir.replace('.tar', '')

        # For tracks_summary files
        filename_summary = os.path.join(
            res_subdir, "".join(
                ('tracks_summary_', deployment, '_', depth, '_', subdir,
                 '.tsv')))

        # For tracks_all files (AA)
        filename_all = os.path.join(
            res_subdir, "".join(
                ('tracks_all_', deployment, '_', depth, '_', subdir, '.tsv')))

        ct_subdir += 1

        # Read an append to the global Dataframe
        if Path(filename_summary).is_file():
            if not bool_file_global:
                df_global = pd.read_csv(filename_summary, sep="\t")
                bool_file_global = True
            else:
                df_summary = pd.read_csv(filename_summary, sep="\t")
                df_global = pd.concat([df_global, df_summary])

        # Read an append to the all Dataframe (AA)
        if Path(filename_all).is_file():
            if not bool_file_all:
                df_all = pd.read_csv(filename_all, sep="\t")
                bool_file_all = True
            else:
                df_all_summary = pd.read_csv(filename_all, sep="\t")
                df_all = pd.concat([df_all, df_all_summary])

    # Save the all dataframe (AA)
    filename_all = "".join(('tracks_all_', deployment, '_', depth, '.tsv'))
    filename_all = os.path.join(dir_tracks, filename_all)
    df_all.to_csv(filename_all, sep="\t", index=False)

    # Save the summary dataframe
    filename_global = "".join(
        ('tracks_summary_', deployment, '_', depth, '.tsv'))
    filename_global = os.path.join(dir_tracks, filename_global)
    df_global.to_csv(filename_global, sep="\t", index=False)

    # Save Ecotaxa files and plots
    # plot_ecotaxa_from_dataframe(filename_all, subdir, deployment, depth)

    # plot the distribution of euclidean distance (AA)
    # sns.histplot(df_all['euclidean_distance'], kde=True, color="blue")
    # Add labels and title
    # plt.xlabel("Euclidean distance")
    # plt.ylabel("Density")
    # Save the plot as a PNG file
    # path = "".join(('hist_dist_tracks_all_', deployment, '_', depth, '.png'))
    # path = os.path.join(dir_tracks, path)
    # plt.savefig(path)

    return
