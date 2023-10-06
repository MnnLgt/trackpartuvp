# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Estimate vertical oscillations of the trap
19-01-2023
@author: Manon Laget
------------------------------
"""


import os
from glob import glob
import pandas as pd
from pathlib import Path
from .utils._utils import parse_datetime


def read_data_files(directory, deployment, depth_depl):
    """
    Read UVP6 data files (there should be one per raw sequence)

    Parameters
    ----------
    directory : str
        Path to UVP6 project. It should contain a subdirectory 'raw' with  
        sequences of images.
    deployment : str
        Deployment name for this project.
    depth_depl : str
        Deployment depth. Do not mistake it with actual trap depth.

    Returns
    -------
    df_all : pandas dataframe
        A dataframe with all depths recorded for each sequence.

    """
    
    # Make sure that the "raw" directory exists
    dir_raw = os.path.join(directory, 'raw')
    if not os.path.exists(dir_raw):
        print('Raw directory does not exist')
        return
    
    # Get subdirectories names (=sequence names)
    os.chdir(dir_raw)
    subdirs = sorted(glob('*'))[1:-1]
    
    # Create global dataframe
    df_all = pd.DataFrame(
        columns = ['sequence', 'imgname', 'datetime', 'averaged_depth'])
    
    for subdir in subdirs:
        
        # To store a list of dictionaries with file lines
        list_dics = []
        
        # Path to data file
        filename = "".join((subdir, '_data.txt'))
        path_to_file = os.path.join(dir_raw, subdir, filename)
    
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
            3, center = True).mean()
        df_seq = df_seq.drop(['depth'], axis=1)
        df_all = df_all.append(df_seq, ignore_index = True)
    
    # Calculate difference with mean deployment depth
    mean_depth = df_all['averaged_depth'].mean()
    df_all['depth_anomaly'] = df_all['averaged_depth'] - mean_depth
    
    # Save csv file
    filename = "".join(('trap_depth_', deployment, '_', depth_depl, '.csv'))
    df_all.to_csv(os.path.join(directory, 'results', filename), index = False)
    
    return(df_all)


def calc_trap_speed(directory, deployment, depth):
    """
    Calculate trap speed for the duration of each track.

    Parameters
    ----------
    directory : str
        Path to UVP6 project. It should contain a subdirectory 'raw' with  
        sequences of images.
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
    dir_results = os.path.join(directory, 'results')
    filename = os.path.join(
        dir_results, "".join((deployment, '_', depth, '_trap_depth.csv')))
    
    # If file doesn't exist, run the function 'read_data_files'
    if Path(filename).is_file():
        df_depths = pd.read_csv(filename)
    else:
        df_depths = read_data_files(directory, deployment, depth)
    
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
        df_tracks = pd.read_csv(filename, sep = "\t")
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
            (df_depths['datetime']>=datetime_ini) &
            (df_depths['datetime']<=datetime_fin)]['averaged_depth'].tolist()
        
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
        if df_tracks['angle_mean'].iloc[i]<180:
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
            ('trap_speeds_', deployment, '_', depth, '.csv')))
    df_speeds.to_csv(filename, index = False)
    
    return


