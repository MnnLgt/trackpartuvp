#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Various functions
2022
@author: Manon Laget
------------------------------
"""


import os
import numpy as np
import cv2
from datetime import datetime
import matplotlib.pyplot as plt
# from skimage import img_as_ubyte


def create_dir(dir_name):
    """
    Checks if a directory exists; if not, creates it.

    Parameters
    ----------
    dir_name : str
        Name of the directory

    """

    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

    return


def resize_image(img, scale_percent):
    """
    Resizes an image to a given percentage of the original.

    Parameters
    ----------
    img : numpy darray
        Image to resize
    scale_percent : float
        Scaling percent for resizing

    Returns
    -------
    resized : numpy darray
        The resized image

    """

    width = int(img.shape[1] * scale_percent / 100)
    height = int(img.shape[0] * scale_percent / 100)

    if scale_percent <= 100:
        resized = cv2.resize(img, (width, height),
                             interpolation=cv2.INTER_AREA)
    else:
        resized = cv2.resize(
            img, (width, height), interpolation=cv2.INTER_LINEAR)

    return resized


def invert_grayscale(img):
    """ Inverts the grayscale of an image. """

    inverted_img = 255 - img

    return inverted_img


def convert_to_gray(input_image_path, output_image_path):
    """ Converts an image to grayscale and saves it. """

    img = cv2.imread(input_image_path)
    grayscale = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)

    try:
        os.remove(input_image_path)
    except input_image_path.DoesNotExist:
        print('')

    cv2.imwrite(output_image_path, grayscale)


def parse_datetime(filename):
    """ Converts UVP6 image name to standard Python datetime format.

    Parameters
    ----------
    filename: str
        the filename of an image

    Returns
    -------
    datetime object
        the formated datetime
    """

    date = filename.split('-')[0]
    time = filename.split('-')[1]

    y, mo, d = int(date[0:4]), int(date[4:6]), int(date[6:])
    h, mi, s = int(time[0:2]), int(time[2:4]), int(time[4:6])

    if (len(filename.split('-')) > 2):

        nb = filename.split('-')[2]
        ms = 0 if int(nb) == 1 else 5*(10**5)

        return datetime(y, mo, d, h, mi, s, ms)

    else:
        ms = 0
        return datetime(y, mo, d, h, mi, s)


def create_window():
    """ Creates a window 2464 x 2056 px window to plot tracks. """

    x = int(round(2464*73*1e-3, 0))
    y = int(round(2056*73*1e-3, 0))
    img = np.ones(shape=[y, x, 3], dtype=np.uint8)*255
    plt.imshow(img, cmap='gray')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
