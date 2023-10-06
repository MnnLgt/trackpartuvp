# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Functions to detect particles on raw images
23-08-2022
@author: Manon Laget
------------------------------
"""


import os
from math import ceil
import numpy as np
from skimage.io import imread
import pandas as pd
from scipy import stats
from skimage.measure import label, regionprops
from trackpartuvp.utils._utils import parse_datetime


class Thresholding():
    """
    Detects particles on a set of UVP6 images by applying a threshold.
    
    Attributes
    ----------
    A, B : float
        Calibration factor (constant depending of the UVP) to convert 
        surface area in px to surface area in µm
    threshold : integer
        Threshold to binarize image
    
    Methods
    -------
    run
        Run thresholding on a sequence of images
    
    """


    def __init__(self, A, B, threshold):
        
        # Set UVP6 parameters
        self.A = A
        self.B = B
        self.threshold = threshold
        
        # Set a minimum ESD to keep morphological data
        self.min_esd = 500
        
        # Set image counter to keep image order
        self.img_counter = 0
            
        # Set morphological properties we want to get
        self.properties = [
            'area', 'area_convex', 'axis_major_length', 'axis_minor_length',
            'centroid', 'coords', 'eccentricity', 'extent', 'label', 
            'orientation', 'perimeter', 'solidity']
        
    
    def run(self, subdir, img_names):
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
        
        """
        
        all_particles, all_areas = ([] for i in range(2))
        
        for img_name in img_names:
            
            try:
                self.img_counter += 1
                img = imread(img_name)
            
                # Binarize image
                img_binary = img > self.threshold
                
                # Get blobs properties
                props = regionprops(label(img_binary))
                
                list_area = []
    
                # Datetime
                datetime = os.path.split(img_name)
                datetime = datetime[1].split('.')[0]
    
                for i in range(len(props)):
    
                    # Coordinates
                    center = props[i]['centroid']
                    
                    # Area
                    raw_area = props[i]['area']
                    
                    if not raw_area>=2:
                        continue
                    
                    # Append area for size spectrum
                    list_area.append(ceil(raw_area))
                    
                    # Correcting area with A and B
                    area = self.A * (raw_area**self.B)

                    # Equivalent Spherical Diameter
                    esd = np.sqrt(4*area/np.pi)
                    
                    if not esd>=self.min_esd:
                        continue
    
                    # Grey values
                    grey_values = []
                    for coord in props[i]['coords']:
                        grey_values.append(img[coord[0], coord[1]])
                    
                    # Creating a dictonary for the particle
                    p_dict = {
                        'n_img': self.img_counter,
                        'img_name': datetime,
                        'datetime':parse_datetime(datetime),
                        'coord_x': center[1],
                        'coord_y': center[0],
                        'area_px': raw_area,
                        'esd_px': np.sqrt(4*props[i]['area']/np.pi),
                        'perim_px': props[i]['perimeter'],
                        'area_um': area,
                        'esd_um': esd,
                        'major_px': props[i]['major_axis_length'],
                        'minor_px': props[i]['minor_axis_length'],
                        'orientation_particle': props[i]['orientation'],
                        'perimmajor':props[i][
                            'perimeter']/props[i]['major_axis_length'],
                        'elongation': 
                            props[i]['major_axis_length']/props[i][
                                'minor_axis_length'],
                        'circularity': (4*np.pi*props[i][
                            'area'])/(props[i]['perimeter']**2),
                        'meangrey': np.mean(grey_values),
                        'stdgrey': np.std(grey_values),
                        'cvgrey': 100*np.std(grey_values)/np.mean(grey_values),
                        'intgrey': np.sum(grey_values),
                        'mediangrey': np.median(grey_values),
                        'mingrey': np.min(grey_values),
                        'maxgrey': np.max(grey_values),
                        'rangegrey': np.max(grey_values) - np.min(grey_values),
                        'skewgrey': stats.skew(grey_values),
                        'kurtgrey': stats.kurtosis(grey_values),
                        'area_convex': props[i]['convex_area'],
                        'eccentricity': props[i]['eccentricity'], 
                        'extent': props[i]['extent'],
                        'solidity': props[i]['solidity'],
                        'filename': img_name
                        }
                    
                    all_particles.append(p_dict)

            except:
                continue
    
        all_areas = list()
        for val in sorted(np.unique(list_area)):
            all_areas.append({'area': val,
                              'nbr': list_area.count(val)})
        
        # Remove black images from total image count
        imgcount = len(img_names) - (ceil(len(img_names)/20))
        
        # Convert list of dictionaries to dataframe
        df_particles = pd.DataFrame(all_particles)
        df_areas = pd.DataFrame(all_areas)
        
        # Add unique id to each particle
        df_particles.insert(
            0, 'part_id', range(1, (len(df_particles) + 1)), True)
        
        # Add subdir name to the dataframe
        df_particles.insert(0, 'sequence', subdir, True)
        df_areas.insert(0, 'imgcount', imgcount, True)
        df_areas.insert(0, 'sequence', subdir, True)
        
        return(df_particles, df_areas)
    
    