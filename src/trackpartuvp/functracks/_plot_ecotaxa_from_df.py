# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Save plots for EcoTaxa
2024
@author: Manon Laget, modified by Alexandre Accardo
------------------------------
"""


import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.io import imread
from trackpartuvp.utils._utils import invert_grayscale, resize_image, convert_to_gray
from trackpartuvp.functracks._plots import add_points, make_title


# Define a Track class to represent a track and hold its information
class Track_from_df:

    def __init__(self, track_data):

        self.id = track_data['track_id'].values[0]
        self.length = len(track_data)
        self.coords = list(zip(track_data['coord_x'], track_data['coord_y']))
        self.vertical_speed = track_data['vertical_speed'].values[0]
        self.datetimes = pd.to_datetime(track_data['datetime'])
        self.img_names = track_data['img_name'].tolist()
        self.mean_angle = track_data['angle'].mean()
        self.orientation = track_data['orientation'].values[0]
        self.filenames = track_data['filename'].tolist()
        self.particles = [
            {'esd_px': esd_px,
             'esd_um': esd_um,
             'datetime': dt,
             'filename': fn,
             'img_name': img_name,
             'above_threshold': above_threshold}
            for esd_px, esd_um, dt, fn, img_name, above_threshold in zip(
                    track_data['esd_px'],
                    track_data['esd_um'],
                    self.datetimes,
                    self.filenames,
                    self.img_names,
                    track_data['above_threshold'])]
        self.mean_esd = track_data['esd_px'].mean()


def plot_ecotaxa_from_dataframe(df, subdir, deployment, depth,
                                save_plots=True):
    """
    For each track in the DataFrame, plot vignettes of its particles
    alongside its plots. In the meantime, create a .tsv file to import
    data on the server. You can mute plots if you just want the .tsv file.

    Parameters
    ----------
    df : pandas DataFrame
        DataFrame containing track information.
    subdir : str
        Name of the sequence.
    deployment : str
        Name of the deployment of the project.
    depth : str
        Depth of the deployment.
    save_plots: bool, optional
        Whether to save individual plots or not.

    Returns
    -------
    list_dic_tsv : list of dictionaries
    """

    list_dic_tsv = []

    # Iterate over each unique track_id
    for track_id in df['track_id'].unique():

        # Get track data from dataframe
        track_data = df[df['track_id'] == track_id]

        # Create a Track object from the DataFrame rows
        track = Track_from_df(track_data)

        if save_plots:

            # General layout

            ncols = 4 + math.ceil(track.length/4)
            nrows = 4

            if ncols >= 6:
                fig_width = 16
            elif ncols == 5:
                fig_width = 14
            elif ncols == 4:
                fig_width = 12
            elif ncols == 3:
                fig_width = 10
            elif ncols == 2:
                fig_width = 8
            elif ncols == 1:
                fig_width = 7

            fig = plt.figure(figsize=(fig_width, 6), constrained_layout=True)
            gs = fig.add_gridspec(nrows, ncols)

            # Plotting the track

            ax1 = fig.add_subplot(gs[0:5, :4])
            x = int(round(2464 * 73 * 1e-3, 0))
            y = int(round(2056 * 73 * 1e-3, 0))
            img = np.ones(shape=[y, x, 3], dtype=np.uint8)*255
            plt.imshow(img, cmap='gray')

            add_points(track.coords, 'b.')

            # Adding title
            ttl = make_title(track)
            plt.title(ttl)
            dt1 = track.datetimes[0]
            dt2 = track.datetimes[-1]
            plt.suptitle('From {dt1} to {dt2}'.format(dt1=dt1, dt2=dt2), y=1)
            ax1.set_title(ttl)
            ax1.set_xlabel('x (mm)', size=12)
            ax1.set_ylabel('y (mm)', size=12)

            # Plotting vignettes

            ct_col, ct_row = 4, 4

            for i in range(track.length):

                ncol = math.trunc(ct_col)
                ct_col += 0.25
                nrow = ct_row % 4
                ct_row += 1

                # Reading image
                img = imread(track.filenames[i])

                # Defining the window
                esd_px = np.mean([particle['esd_px']
                                  for particle in track.particles])
                window = esd_px + 4
                right = min(int(track.coords[i][0] + window), 2464)
                left = max(0, int(track.coords[i][0] - window))
                top = max(0, int(track.coords[i][1] - window))
                bott = min(int(track.coords[i][1] + window), 2056)

                # Cropping the original image
                crop = img[top:bott, left:right]

                # Dimensions un µm
                width, height = crop.shape[1] * 73, crop.shape[0] * 73

                # Inverting gray scale and resizing image
                crop = invert_grayscale(resize_image(crop, 750))

                ax = fig.add_subplot(gs[nrow, ncol])
                ttl = str(i) + ' | (µm)'
                ax.set_title(ttl, fontsize=6)
                ax.tick_params(axis="x", labelsize=6)
                ax.tick_params(axis="y", labelsize=6)
                plt.imshow(crop, cmap='gray', extent=[0, width, height, 0],
                           vmin=0, vmax=255)

            fig.tight_layout(h_pad=0.2)

        filename = "".join(
            (deployment, depth, '_',  str(track.img_names[0]), '-',
             str(track.id), '.png'))
        filename = os.path.join(track.orientation, filename)

        if save_plots:

            try:
                plt.savefig(filename, dpi=80)
                plt.close()
            except:
                print('Plot for Ecotaxa: could not save figure')

            convert_to_gray(filename, filename)

        list_dic_tsv.append({
            'img_file_name': filename,
            'img_rank': '0',
            'object_id': "".join(
                (deployment, depth, '_',  str(track.img_names[0]), '-',
                 str(track.id))),
            'object_esd': track.mean_esd,
            'object_orientation': track.orientation,
            'object_sequence': subdir,
            'sample_id': "_".join((deployment, depth)),
            'sample_deployment': deployment,
            'sample_depth': depth
            })

    return list_dic_tsv


# Application
# df = pd.read_csv('/home/aaccardo/uvp6_sn000146lp_2023_apero_float_cycle5/results/detected_tracks/tracks_all_cycle_5_.tsv', sep='\t')
# df['above_threshold'] = df['euclidean_distance'] > df['95_percentile']
# len(df)
# true_df = df[df['above_threshold']== True]
# len(true_df)
# true_track_id = true_df['track_id'].unique().tolist()
# true_track_id
# len(true_track_id)
# track_above = df[df['track_id'].isin(true_track_id)]
# track_above.head(10)
# track_above.shape
# result_list = plot_ecotaxa_from_dataframe(track_above, '', '', '')









def plot_ecotaxa_95_percentile(df_all_path, saving_path):
    """
    For each track of a sequence, plots thumbnails of its particles
    alongside its plots. In the meantime, create a .tsv file to import
    data on the server.

    Parameters
    ----------
    df_all_path : str
    saving_path : str
        Name of the sequence.

    Returns
    -------
    list_dic_tsv : list of dictionaries

    """

    list_dic_tsv = []

    df_all = pd.read_csv(df_all_path, sep='\t')

    # Group the DataFrame by 'track_id'
    df_all = df_all.groupby('track_id')

    # Iterate over track_id
    #for track in list_tracks:
    for track_id, group in df_all:

        # General layout

        ncols = 4 + math.ceil(group['track_length'].mean() / 4)
        nrows = 4
        if ncols >= 6:
            fig_width = 16
        elif ncols == 5:
            fig_width = 14
        elif ncols == 4:
            fig_width = 12
        elif ncols == 3:
            fig_width = 10
        elif ncols == 2:
            fig_width = 8
        elif ncols == 1:
            fig_width = 7

        fig = plt.figure(figsize=(fig_width, 6), constrained_layout=True)
        gs = fig.add_gridspec(nrows, ncols)

        # Plotting the track

        ax1 = fig.add_subplot(gs[0:5, :4])
        x = int(round(2464 * 73 * 1e-3, 0))
        y = int(round(2056 * 73 * 1e-3, 0))
        img = np.ones(shape=[y, x, 3], dtype=np.uint8)*255
        plt.imshow(img, cmap='gray')

        track_coord = list(zip(group['coord_x'], group['coord_y']))
        add_points(track_coord, 'b.')

        # Calculating mean speed and ESD to make title
        ss = int(round(group['vertical_speed'].mean(), 0))
        ss = -ss if group['angle'].mean() < 180 else ss

        esd = int(round(group['esd_px'].mean(), 0))
        esd_std = int(round(group['esd_px'].std(), 0))

        # Creating title
        title = 'Vertical speed = {ss} m/day - ESD = {esd} +/- {std} µm'
        title = title.format(ss=ss, esd=esd, std=esd_std)

        plt.title(title)
        # dt1 = track.datetimes[0]
        # dt2 = track.datetimes[-1]
        dt1 = group['datetime'].min()
        dt2 = group['datetime'].max()
        plt.suptitle('From {dt1} to {dt2}'.format(dt1=dt1, dt2=dt2), y=1)
        ax1.set_title(title)
        ax1.set_xlabel('x (mm)', size=12)
        ax1.set_ylabel('y (mm)', size=12)

        # Plotting vignettes

        ct_col, ct_row = 4, 4

        group = group.reset_index(drop=True)

        for i in range(int(group['track_length'].mean())):

            ncol = math.trunc(ct_col)
            ct_col += 0.25
            nrow = ct_row % 4
            ct_row += 1

            # Reading image
            img = imread(group['filename'][i])

            # Defining the window
            # esd_px = np.mean([particle['esd_px']
            #                 for particle in track.particles])
            esd_px = group['esd_px'].mean()
            window = esd_px + 4
            right = min(int(track_coord[i][0] + window), 2464)
            left = max(0, int(track_coord[i][0] - window))
            top = max(0, int(track_coord[i][1] - window))
            bott = min(int(track_coord[i][1] + window), 2056)

            # Cropping the original image
            crop = img[top:bott, left:right]

            # Dimensions un µm
            width, height = crop.shape[1] * 73, crop.shape[0] * 73

            # Inverting gray scale and resizing image
            crop = invert_grayscale(resize_image(crop, 750))

            ax = fig.add_subplot(gs[nrow, ncol])
            ttl = str(i) + ' | (µm)'
            ax.set_title(ttl, fontsize=6)
            ax.tick_params(axis="x", labelsize=6)
            ax.tick_params(axis="y", labelsize=6)
            plt.imshow(crop, cmap='gray', extent=[0, width, height, 0],
                       vmin=0, vmax=255)

        fig.tight_layout(h_pad=0.2)

        filename = "".join(
            ('_',  str(group['img_name'].iloc[0]), '-',
             str(track_id), '.png'))
        filename = os.path.join(str(group['orientation'].unique()[0]),
                                filename)
        filename = os.path.join(saving_path, filename)

        try:
            plt.savefig(filename, dpi=80)
            plt.close()
        except:
            print('Plot for Ecotaxa: could not save figure')

        convert_to_gray(filename, filename)

    return


'''
import pandas as pd
df_all = pd.read_csv('/home/aaccardo/uvp6_sn000146lp_2023_apero_float_cycle5/results/detected_tracks/tracks_all_cycle_5_.tsv', sep='\t')
df_all.head(5)

list_tracks = df_all['track_id'].unique().tolist()
list_tracks

subdir = df_all['sequence'].unique().tolist()
subdir[:5]

deployment = ''
depth = ''

coordinates_list = list(zip(df_all['coord_x'], df_all['coord_y']))
coordinates_list

df_all = df_all.groupby('track_id')

    # Iterate over track_id
    #for track in list_tracks:
for track_id, group in df_all:
    esd_px = df_all['esd_px'].mean()


esd_px
'''

# plot_ecotaxa(df_all_path = '/home/aaccardo/uvp6_sn000146lp_2023_apero_float_cycle5/results/detected_tracks/tracks_all_cycle_5_.tsv', saving_path = '/home/aaccardo/uvp6_sn000146lp_2023_apero_float_cycle5/results/detected_tracks/plot_ecotaxa_95_percentile/')
