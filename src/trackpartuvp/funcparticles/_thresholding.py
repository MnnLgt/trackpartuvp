# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Functions to detect particles on raw images
2022
@author: Manon Laget, modified by Alexandre Accardo
------------------------------
"""


import os
from math import ceil
import numpy as np
from skimage.io import imread
import pandas as pd
from skimage.measure import label, regionprops
from skimage.restoration import denoise_bilateral
from trackpartuvp.utils._utils import parse_datetime
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re
from datetime import datetime


class Thresholding():
    """
    Detects particles on a set of UVP6 images by applying a threshold.

    Attributes
    ----------
    A, B: float
        Calibration factor (constant depending of the UVP) to convert surface
        rea in px to surface area in µm.
    threshold: int
        Threshold to binarize images.

    Methods
    -------
    pca_morpho
        Perform pca on morphological features.
    extract_datetime_from_filename
        Transform filename into datetime object.
    run
        Run thresholding on a sequence of images

    """

    def __init__(self, A, B, threshold):

        # Set UVP6 parameters
        self.A = A
        self.B = B
        self.threshold = threshold

        # Set a minimum ESD to keep morphological data
        self.min_esd = 100
        self.min_area_px = 3
        # Originally set up to 500 microns

        # Set image counter to keep image order
        self.img_counter = 0

        # Set morphological properties we want to get
        self.properties = [
            'area', 'area_convex', 'axis_major_length', 'axis_minor_length',
            'centroid', 'coords', 'eccentricity', 'extent', 'label',
            'orientation', 'perimeter', 'solidity']

    def pca_morpho(self, df_particles):
        """
        Function to perform pca on morphological features. The pca is build on
        all the particles detected on a sequence.
        Coordinates of each particle on the principal components are stored in
        the particle_df csv file.
        """

        # List morphological features that will be used for PCA
        # (I removed the correlated ones)
        morpho_features = ['area_um', 'esd_um', 'major_px', 'minor_px',
                           'orientation_particle', 'perimmajor', 'meangrey',
                           'eccentricity', 'extent', 'solidity']
        # 'skewgrey','kurtgrey' 'elongation', 'circularity'

        # morpho_features = ['coord_x', 'coord_y', 'esd_um', 'meangrey']

        # Select the desired columns
        df_morpho = df_particles[morpho_features]

        # Normalize features
        df_morpho_norm = StandardScaler().fit_transform(df_morpho)
        df_morpho_norm = pd.DataFrame(df_morpho_norm, columns=morpho_features)

        # Perform pca
        pca = PCA(n_components=5)
        PrincComp = pca.fit_transform(df_morpho_norm)
        pca_df = pd.DataFrame(data=PrincComp,
                              columns=['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])

        # join to df_particles
        df_particles = pd.concat([df_particles, pca_df], axis=1)

        return df_particles

    def extract_datetime_from_filename(filename):

        # Define the regex pattern to match the filename structure
        pattern = r"livecamera(\d{8})_(\d{6})_\d{3}\.png"

        # Use regex to extract date and time components
        match = re.match(pattern, filename)
        if not match:
            raise ValueError("Invalid filename format")

        # Extract date and time values from the regex match
        date_str = match.group(1)  # e.g., '20240424'
        time_str = match.group(2)  # e.g., '154943'

        # Convert extracted strings into datetime object
        date_time_str = f"{date_str}_{time_str}"
        date_time = datetime.strptime(date_time_str, "%Y%m%d_%H%M%S")

        # Format datetime object into desired string format
        formatted_datetime = date_time.strftime("%Y-%m-%d %H:%M:%S")

        return formatted_datetime

    def run(self, subdir, img_names, platform_speed_path=None):
        """
        Runs thresholding on a sequence of images.

        Parameters
        ----------
        subdir : str
            Name of the sequence.
        img_names : list of str
            A sequence of image names.

        Returns
        -------
        df_particles : pandas dataframe
            Particles morphological properties (particle size >min_esd)
        df_areas : pandas dataframe
            Areas of particles >100 µm
        png plot with the centroid positions of object > 5 pixels during the
        entire sequence

        """

        # Create lists to store particles and areas
        all_particles, all_areas = ([] for i in range(2))
        centroid_positions = []

        # Initialize an empty list to store particle data for PCA (AA)
        # pca_data = []

        # Load platform vertical speed for each sequence if file exists (AA)
        if platform_speed_path:
            csv_file_path = platform_speed_path  # path to the csv
            csv_data = pd.read_csv(csv_file_path)  # import the data
            seq_speed_mapping = dict(
                zip(csv_data['seq'],
                    csv_data['float_vertical_speed'])
                )  # create a dictionary mapping
            # Get the 'mean_speed_sequence' value from the dictionary
            mean_speed_sequence = seq_speed_mapping.get(subdir, "N/A")

        img_count = 0
        total_img = len(img_names)

        for img_name in img_names:

            img_count += 1
            print(f'Processing image {img_count}/{total_img}: {img_name}')

            try:
                self.img_counter += 1
                img = imread(img_name)

                # Binarize image
                img_binary = img > self.threshold

                # reduce noise
                img_binary = denoise_bilateral(img_binary, sigma_spatial=1)

                # Get blobs properties
                props = regionprops(label(img_binary))

                # list_area = []

                # Datetime
                datetime = os.path.split(img_name)
                datetime = datetime[1].split('.')[0]

                for i in range(len(props)):

                    # Coordinates
                    center = props[i]['centroid']

                    # Area (Area of the region i.e. number of pixels of the
                    # region scaled by pixel-area)
                    raw_area = props[i]['area']

                    if not raw_area >= 2:
                        continue

                    # Append area for size spectrum
                    # list_area.append(ceil(raw_area))

                    # Correcting area with A and B
                    area = self.A * (raw_area**self.B)

                    # if area >= 5:
                    #   centroid_positions.append(center)
                    # select only the centroid positions of object >5 µm^2

                    # Append center to the list
                    centroid_positions.append(center)

                    # Equivalent Spherical Diameter
                    esd = np.sqrt(4*area/np.pi)

                    if not esd >= self.min_esd:
                        # if not raw_area >= self.min_area_px:
                        continue

                    # Grey values
                    grey_values = []
                    for coord in props[i]['coords']:
                        grey_values.append(img[coord[0], coord[1]])

                    elong = None
                    circ = None
                    if raw_area >= 5:
                        try:
                            elong = props[i][
                                'major_axis_length']/props[i]['minor_axis_length']
                        except:
                            print(" ")
                        try:
                            circ = (4*np.pi*props[i][
                                'area'])/(props[i]['perimeter']**2)
                        except:
                            print(" ")
                        # it has to
                        # be a value between 0 and 1; 0 = elongated object and
                        # 1 = perfect circle but I have values higher than 1
                        # and I have an inversed relationship

                    # Creating a dictonary for the particle
                    p_dict = {
                        'n_img': self.img_counter,
                        'img_name': datetime,
                        'datetime': parse_datetime(datetime),
                        # 'datetime': extract_datetime_from_filename(datetime),
                        'coord_x': center[1],
                        'coord_y': center[0],
                        'area_px': raw_area,
                        'esd_px': np.sqrt(4*props[i]['area']/np.pi),
                        # Perimeter of object which approximates the contour as
                        # a line through the centers of border pixels using a
                        # 4-connectivity
                        'perim_px': props[i]['perimeter'],
                        'area_um': area,
                        'esd_um': esd,
                        'major_px': props[i]['major_axis_length'], # The length of the major axis of the ellipse that has the same normalized second central moments as the region
                        'minor_px': props[i]['minor_axis_length'], # The length of the minor axis of the ellipse that has the same normalized second central moments as the region
                        'orientation_particle': props[i]['orientation'], # Angle between the 0th axis (rows) and the major axis of the ellipse that has the same second moments as the region, ranging from -pi/2 to pi/2 counter-clockwise
                        'perimmajor': props[i]['perimeter'
                                               ]/props[i]['major_axis_length'],
                        'elongation': elong,
                        'circularity': circ,
                        # Mean grey: 0 = white ; 255 = black (in this case
                        # because the images are binarized)
                        'meangrey': np.mean(grey_values),
                        'stdgrey': np.std(grey_values),
                        'cvgrey': 100*np.std(grey_values)/np.mean(grey_values),
                        'intgrey': np.sum(grey_values),
                        'mediangrey': np.median(grey_values),
                        'mingrey': np.min(grey_values),
                        'maxgrey': np.max(grey_values),
                        'rangegrey': np.max(grey_values) - np.min(grey_values),
                        # 'skewgrey': stats.skew(grey_values),
                        # 'kurtgrey': stats.kurtosis(grey_values),
                        'area_convex': props[i]['convex_area'],
                        'eccentricity': props[i]['eccentricity'],
                        'extent': props[i]['extent'],
                        'solidity': props[i]['solidity'],
                        'filename': img_name
                        }

                    all_particles.append(p_dict)

            except:
                continue

        # all_areas = list()
        # for val in sorted(np.unique(list_area)):
        #    all_areas.append({'area': val,
        #                      'nbr': list_area.count(val)})

        # Display a plot with the stacked position of all centroids for objects > 5 pixels only # AA

        # centroid_positions = np.array(centroid_positions)
        # plt.figure(figsize=(8, 8))
        # plt.scatter(centroid_positions[:, 1], centroid_positions[:, 0], c='black', marker='.', s=1)
        # plt.gca().invert_yaxis()
        # plt.title(f'Centroid Positions for Sequence: {subdir} ; Float speed : {mean_speed_sequence} m/d')
        # plt.xlabel('X Position')
        # plt.ylabel('Y Position')
        # plt.xlim(0, 2500)
        # plt.ylim(2100, 0)

        # Generate a unique name # AA
        # plot_filename = f"{subdir}_stacked_position.png"

        # Specify the output directory # AA
        # output_directory = '/home/aaccardo/uvp6_sn000146lp_2023_apero_float_cycle1/results/object_detection_results' # have to change at each new run on new deployment (I have to adapt it !)

        # Create the full path for the output file # AA
        # output_path = os.path.join(output_directory, plot_filename)

        # Save the plot as a PNG file # AA
        # plt.savefig(output_path)

        # print(subdir + ' plot successfully saved') # AA

        # Remove black images from total image count
        imgcount = len(img_names) - (ceil(len(img_names)/20))

        # Convert list of dictionaries to dataframe
        df_particles = pd.DataFrame(all_particles)
        # df_areas = pd.DataFrame(all_areas)

        # Add unique id to each particle
        df_particles.insert(
            0, 'part_id', range(1, (len(df_particles) + 1)), True)

        # Add subdir name to the dataframe
        df_particles.insert(0, 'sequence', subdir, True)
        # df_areas.insert(0, 'imgcount', imgcount, True)
        # df_areas.insert(0, 'sequence', subdir, True)
        # print('df_particles: ', df_particles)
        df_particles = self.pca_morpho(df_particles)

        # return (df_particles, df_areas)
        return df_particles
