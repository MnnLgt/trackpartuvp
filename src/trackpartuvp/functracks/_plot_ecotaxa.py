# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Save plots for EcoTaxa
2022
@author: Manon Laget
------------------------------
"""


import math
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.io import imread
from trackpartuvp.utils._utils import invert_grayscale, resize_image, convert_to_gray
from ._plots import add_points, make_title


def plot_ecotaxa(list_tracks, subdir, deployment, depth, save_plots=True):
    """
    For each track of a sequence, plot vignettes of its particles
    alongside its plots. In the meantime, create a .tsv file to import
    data on the server. You can mute plots if you just want the .tsv file.

    Parameters
    ----------
    list_tracks: list of Track objects
        Contains the tracks of a sequence.
    subdir: str
        Name of the sequence.
    deployment: str
        Name of the deployment of the project.
    depth: str
        Depth of the deployment.
    save_plots: bool, optional
        Whether to save individual plots or not.

    Returns
    -------
    list_dic_tsv: list of dictionaries

    """

    list_dic_tsv = []

    for track in list_tracks:

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
