# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Estimate vertical oscillations of the platform from UVP6 data files
2023
@author: Manon Laget
------------------------------
"""


import os
from glob import glob
import pandas as pd
import tarfile
import shutil
from pathlib import Path
from .utils._utils import parse_datetime, create_dir


def read_data_files(source_directory, dest_directory, deployment, depth_depl):
    """
    Read UVP6 data files (there should be one per raw sequence)

    Parameters
    ----------
    source_directory: str
        Path to UVP6 project. It should contain a subdirectory 'raw' with
        sequences of images.
    dest_directory: str
        Path to directory where to store data.
    deployment: str
        Deployment name for this project.
    depth_depl: str
        Deployment depth. Do not mistake it with actual platform depth.

    Returns
    -------
    df_all : pandas dataframe
        A dataframe with all depths recorded for each sequence.

    """

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

    # Create subdirectories in the results directory for each sequence
    os.chdir(dir_raw)
    subdirs = sorted(glob('*'))
    print(len(subdirs), "sequences")

    # Create global dataframe
    df_all = pd.DataFrame(
        columns=['sequence', 'imgname', 'datetime', 'averaged_depth'])

    for subdir in subdirs:

        # Check if there is a tar file
        path_tar = os.path.join(dir_raw, subdir)

        if os.path.exists(path_tar):
            print('tar file found')
            with tarfile.open(path_tar) as tar:
                tar.extractall(path=os.path.join(dir_copy_raw, subdir))

        # To store a list of dictionaries with file lines
        list_dics = []

        # Path to data file
        filename = "".join((subdir, '_data.txt'))
        filename = filename.replace('.tar', '')
        if os.path.exists(path_tar):
            subdir2 = subdir.replace('.tar', '')
            path_to_file = os.path.join(dir_copy_raw, subdir, subdir2, filename)
        else:
            path_to_file = os.path.join(dir_copy_raw, subdir, filename)

        # Read data file
        with open(path_to_file, 'r') as file:

            lines = file.readlines()
            for line in lines[4:]:

                # Skip blank lines
                if not line.strip():
                    continue

                # Get 1st and 2nd columns
                datetime = line.split(',')[0]
                depth = float(line.split(',')[1])
                dic = {'sequence': subdir,
                       'depth': depth,
                       'imgname': datetime,
                       'datetime': parse_datetime(datetime)
                       }
                list_dics.append(dic)

        # Make a dataframe from all dictionaries
        df_seq = pd.DataFrame(list_dics)

        # Moving average for the depth
        df_seq['averaged_depth'] = df_seq['depth'].rolling(
            3, center=True).mean()
        df_seq = df_seq.drop(['depth'], axis=1)
        # Method append disappeared with pandas 2.0
        # df_all = df_all.append(df_seq, ignore_index = True)
        # Instead use the method concat
        df_all = pd.concat([df_all, df_seq])

        if os.path.exists(path_tar):
            try:
                shutil.rmtree(os.path.join(dir_copy_raw, subdir))
            except:
                continue

    # Calculate difference with mean deployment depth
    mean_depth = df_all['averaged_depth'].mean()
    df_all['depth_anomaly'] = df_all['averaged_depth'] - mean_depth

    # Save csv file
    filename = "".join(
        ('platform_depth_', deployment, '_', depth_depl, '.csv'))
    df_all.to_csv(os.path.join(dest_directory, 'results', filename),
                  index=False)

    return df_all


def calc_pltf_speed(source_directory, dest_directory, deployment, depth):
    """
    Calculate platform speed for the duration of each track.

    Parameters
    ----------
    source_directory: str
        Path to UVP6 project. It should contain a subdirectory 'raw' with
        sequences of images.
    dest_directory: str
        Path to directory where to store data.
    deployment : str
        Deployment name for this project.
    depth : str
        Deployment depth. Do not mistake it with actual trap depth.

    Raises
    ------
    Exception
        If cannot read track files (tracking must be performed first).

    Returns
    -------
    None.

    """

    # Result directory to read trap depth file
    dir_results = os.path.join(dest_directory, 'results')
    filename = os.path.join(
        dir_results, "".join(
            ('platform_depth_', deployment, '_', depth, '.csv')))

    # If file doesn't exist, run the function 'read_data_files'
    if Path(filename).is_file():
        df_depths = pd.read_csv(filename)
    else:
        df_depths = read_data_files(source_directory, dest_directory,
                                    deployment, depth)

    # Parse datetime
    list_datetimes = df_depths['imgname'].tolist()
    df_depths['datetime'] = [parse_datetime(dt) for dt in list_datetimes]
    df_depths.sort_values(by='datetime')

    # Track filename
    filename = os.path.join(
        dir_results, 'detected_tracks', "".join(
            ('tracks_summary_', deployment, '_', depth, '.tsv')))

    # Try to read track file. If not present, raise exception
    try:
        df_tracks = pd.read_csv(filename, sep="\t")
    except:
        raise Exception('Could not read track file')

    list_dic = []

    #
    for i in range(len(df_tracks)):

        # Get time range of the track
        datetime_ini = parse_datetime(df_tracks['img_ini'].iloc[i])
        datetime_fin = parse_datetime(df_tracks['img_fin'].iloc[i])

        # Locate depth of the trap
        depth_trap = df_depths.loc[
            (df_depths['datetime'] >= datetime_ini) &
            (df_depths['datetime'] <= datetime_fin)]['averaged_depth'].tolist()

        if not depth_trap:
            continue

        # Calculate trap speed
        dist = 0
        for j in range(1, len(depth_trap)):
            dist += (depth_trap[j] - depth_trap[j-1])

        dt = (datetime_fin - datetime_ini).total_seconds()

        # Convert it to cm/s
        trap_speed = (dist/dt)*100

        track_speed = df_tracks['vertical_speed'].iloc[i]
        if df_tracks['angle_mean'].iloc[i] < 180:
            track_speed = -track_speed

        dic = {'sequence': df_tracks['seq'].iloc[i],
               'datetime_ini': datetime_ini,
               'datetime_fin': datetime_fin,
               'track_id': df_tracks['track_id'].iloc[i],
               'track_speed': track_speed,
               'trap_speed': trap_speed
               }

        list_dic.append(dic)

    df_speeds = pd.DataFrame(list_dic)

    filename = os.path.join(
        dir_results, "".join(
            ('platform_speeds_', deployment, '_', depth, '.csv')))
    df_speeds.to_csv(filename, index=False)

    return
