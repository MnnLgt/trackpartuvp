# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Plot each track
2022
@author: Manon Laget
------------------------------
"""

import os
import numpy as np
import pandas as pd
from skimage.io import imread, imsave
from cv2 import circle
import matplotlib.pyplot as plt
from trackpartuvp.utils._utils import invert_grayscale, resize_image, create_window


def add_points(coords, color='b.', size=8):
    """
    Adds points of a track to a plot and annotates them.

    Parameters
    ----------
    coords: list of tuples
        Coordinates of each point of the track.
    color : str, optional
        Color and type of the points. Default is 'b.'.
    size : int, optional
        Size of the annotations. Default is 8.

    """

    ct = -1

    for j in range(len(coords)):

        pos_x, pos_y = coords[j][0] * 73 * 1e-3, coords[j][1] * 73 * 1e-3
        plt.plot(pos_x, pos_y, color)
        ct += 1

        if ct % 4 == 0 or ct == (len(coords)-1):
            xtext = pos_x - 2 if pos_x > (179.8/2) else pos_x + 2
            ytext = pos_y - 2
            plt.text(xtext, ytext, 't' + str(j), size=size)


def make_title(track):
    """ Create a title with the size and speed of a particle.

    Parameters
    ----------
    track : Track object
        Track of a particle.

    Returns
    -------
    title : str
        The formatted title.

    """

    # Calculating mean speed and ESD to make title
    ss = int(round(track.vertical_speed, 0))
    ss = -ss if track.mean_angle < 180 else ss

    esd = int(round(np.mean([particle['esd_um']
                             for particle in track.particles]), 0))
    esd_std = int(round(np.std([particle['esd_um']
                                for particle in track.particles]), 0))

    # Creating title
    title = 'Vertical speed = {ss} m/day - ESD = {esd} +/- {std} µm'
    title = title.format(ss=ss, esd=esd, std=esd_std)

    return title


class Plots:
    """
    Various plots to visualize the detected tracks.

    Attributes
    ----------
    list_tracks : list of Track objects
        The tracks to be plotted.
    df_particles : pandas DataFrame
        Contains all the particles detected in a sequence.
    directory : str
        Path of the directory where are saved the plots.

    Methods
    -------
    plot
    extract_particles
    save_raw_images
    run

    """

    def __init__(self, list_tracks, df_particles, directory):

        self.list_tracks = list_tracks
        self.df_particles = df_particles
        self.directory = directory
        self.current_track_dir = None
        self.track = None

    def plot(self):
        """
        Plots the positions of a track on a single graph.

        Returns
        -------
        None.

        """

        track = self.track

        # Getting particle coordinates
        coords = track.coords

        # Creating the window
        create_window()

        # Plotting the track and the timesteps
        add_points(coords, 'b.')

        # Adding title
        dt1, dt2 = track.datetimes[0], track.datetimes[-1]
        plt.title(make_title(track))
        plt.suptitle('From {dt1} to {dt2}'.format(dt1=dt1, dt2=dt2), y=1)

        # Saving plot
        filename = "".join((str(track.id), '.png'))
        filename = os.path.join(self.current_track_dir, filename)
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def extract_particles(self):
        """
        Extracts the particles of a track from raw images and create
        a stroboscope.

        Returns
        -------
        None.

        """

        track = self.track

        # Setting the image to stack particles
        stroboscope = np.ones(shape=[2056, 2464], dtype=np.uint8)*0

        for i in range(track.length):

            # Reading image
            img_name = track.filenames[i]
            img = imread(img_name)

            # Defining the window
            esd_px = np.mean([particle['esd_px']
                              for particle in track.particles])
            window = esd_px*1.5
            right = min(int(track.coords[i][0] + window), 2464)
            left = max(0, int(track.coords[i][0] - window))
            top = max(0, int(track.coords[i][1] - window))
            bott = min(int(track.coords[i][1] + window), 2056)

            # Stacking image
            to_stack = img
            to_stack[0:top, :] = 0
            to_stack[:, 0:left] = 0
            to_stack[bott:2056, :] = 0
            to_stack[:, right:2464] = 0
            stroboscope = stroboscope + to_stack

            # Cropping the original image
            crop = img[top:bott, left:right]

            # Inverting gray scale and resizing image
            crop = resize_image(crop, 300)
            crop = invert_grayscale(crop)

            # Dimensions un µm
            width = (crop.shape[1] * 73) / 5
            height = (crop.shape[0] * 73) / 5

            plt.imshow(crop, cmap='gray', extent=[0, width, height, 0],
                       vmin=0, vmax=255)
            plt.xlabel('x (µm)')
            plt.ylabel('y (µm)')

            # Saving the vignette
            filename = "".join(
                (str(track.img_names[i]), '-', str(track.id), '.png'))
            filename = os.path.join(self.current_track_dir, filename)
            imsave(filename, crop)

        # Stroboscope
        stroboscope = invert_grayscale(stroboscope)

        # Cropping the stroboscope

        # Zooming in
        x = [track.coords[i][0] for i in range(track.length)]
        y = [track.coords[i][1] for i in range(track.length)]
        left = int(min(x)) - 100
        right = int(max(x)) + 100
        top = int(min(y)) - 100
        bott = int(max(y)) + 100

        right = min(right, 2464)
        left = max(0, left)
        top = max(0, top)
        bott = min(bott, 2056)

        crop_strobo = stroboscope[top:bott, left:right]
        crop_strobo = resize_image(crop_strobo, 300)

        # Adding axes to the cropped stroboscope
        im_crop = np.array(crop_strobo)
        ext = np.array([left, right, bott, top])*73*1e-3
        plt.imshow(im_crop, cmap='gray', extent=ext, vmin=0, vmax=255)
        plt.xlabel('x (mm)')
        plt.ylabel('y (mm)')

        # Saving cropped stroboscope
        filename = "".join(('strobo-', str(track.id), '-crop.png'))
        filename = os.path.join(self.current_track_dir, filename)
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def save_raw_images(self):
        """
        For a set of tracks, saves the raw images of each with
        contoured particles and saves the 2 frames before and after.

        Returns
        -------
        None.

        """

        track = self.track
        df_particles = self.df_particles

        # Saving frames of the track with contoured particles

        for i in range(track.length):

            # Reading image and circling particle
            img = imread(track.filenames[i])
            x, y = int(track.coords[i][0]), int(track.coords[i][1])
            esd_px = np.mean([particle['esd_px']
                              for particle in track.particles])
            rad = int(esd_px*4)
            img = circle(img, center=(x, y), radius=rad, color=(255, 0, 0),
                         thickness=8)
            img = invert_grayscale(img)

            # Saving the image
            filename = "".join(
                (str(track.img_names[i]), '_circled-parts.jpeg'))
            filename = os.path.join(self.current_track_dir, filename)
            imsave(filename, img, quality=70)

        # Saving the two frames before and after

        n_beginning = track.particles[0]['n_img']
        n_end = track.particles[-1]['n_img']

        list_img = df_particles['n_img'].unique()

        for i in range(1, 3, 1):

            n_img = n_beginning - i

            if n_img in list_img:

                df_sub = df_particles.loc[df_particles['n_img'] == n_img]

                im_name = df_sub.iloc[0]['filename']
                frame = os.path.split(im_name)
                frame = frame[1].split('.')[0]
                img = invert_grayscale(imread(im_name))

                filename = "".join((str(frame), '-frame', str(n_img), '.jpeg'))
                filename = os.path.join(self.current_track_dir, filename)
                imsave(filename, img, quality=70)

            n_img = n_end + i

            if n_img in list_img:

                df_sub = df_particles.loc[df_particles['n_img'] == n_img]

                im_name = df_sub.iloc[0]['filename']
                frame = os.path.split(im_name)
                frame = frame[1].split('.')[0]
                img = invert_grayscale(imread(im_name))

                filename = "".join((str(frame), '.jpeg'))
                filename = os.path.join(self.current_track_dir, filename)
                imsave(filename, img, quality=70)

        return

    def run(self):
        """ Saves all the plots. """

        for track in self.list_tracks:

            self.current_track_dir = os.path.join(
                self.directory, str(track.id))

            self.track = track

            # Create it folder
            if not os.path.exists(self.current_track_dir):
                os.mkdir(self.current_track_dir)

            # Save plots
            # self.plot() # to keep
            # self.extract_particles() # to keep
            # self.save_raw_images()

        return
