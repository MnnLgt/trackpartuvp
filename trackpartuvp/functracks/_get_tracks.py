# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Save tracks
23-08-2022
@author: Manon Laget
------------------------------
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trackpartuvp.utils._utils import create_window


class GetTracks:
    """
    Gets tracks of marine particles.
    
    Attributes
    ----------
    list_tracks : list of Track objects.
        Tracks of a project.
    subdir : str
        Name of the sequence.
    deployment : str
        Name of the deployment for this project.
    depth : str
        Depth of the deployment.
    
    Methods
    -------
    save_tracks
    plot_all_tracks
    run
    
    """
    
    
    def __init__(self, list_tracks, subdir, deployment, depth):
        
        # Check if deployment and depth are strings
        
        if not type(deployment) is str:
            raise TypeError('Deployment name should be of string type')
        
        if not type(depth) is str:
            raise TypeError('Depth should be of string type')
        
        self.list_tracks = list_tracks
        self.subdir = subdir
        
        # Name and depth of the deployment
        self.deployment = deployment
        self.depth = depth


    def save_df_summary(self):
        """ Save dataframe with summary of each track. """
        
        list_keys = [
            'area_px', 'esd_px', 'perim_px','area_um', 'esd_um', 'major_px', 
            'minor_px', 'perimmajor', 'elongation', 'circularity', 'meangrey',
            'stdgrey', 'cvgrey', 'intgrey', 'mediangrey', 'mingrey', 'maxgrey',
            'rangegrey', 'skewgrey', 'kurtgrey', 'area_convex', 'eccentricity',
            'extent', 'solidity']
        
        list_dicts = []
        
        for track in self.list_tracks:
            
            track_id = "".join(
                (self. deployment, self.depth, '_',  str(track.img_names[0]), 
                 '-', str(track.id)))
            
            list_rad = [particle['orientation_particle'] 
                        for particle in track.particles]
            n = len(list_rad)
            sin_alpha = np.sum(np.sin(list_rad))
            cos_alpha = np.sum(np.cos(list_rad))
            sin_alpha, cos_alpha = sin_alpha/n, cos_alpha/n
            alpha = np.arctan2(sin_alpha, cos_alpha)
            alpha = 360 - ((alpha*180/np.pi)%360)
            
            dict_track = {
                'seq': self.subdir,
                'track_id': track_id,
                'length': track.length,
                'img_ini': track.img_names[0],
                'img_fin': track.img_names[-1],
                'datetime_ini': track.datetimes[0],
                'datetime_fin': track.datetimes[-1],
                'angle_mean': track.mean_angle,
                'angle_std': track.std_angle,
                'mean_orientation_particle': alpha,
                'speed': track.speed,
                'vertical_speed': track.vertical_speed,
                'vx': track.vx,
                'orientation': track.orientation,
                'sinusoity_index': track.sinusoity_index,
                'step_length_mean': track.mean_step,
                'step_length_std': track.std_step,
                'longest_distance': track.longest_dist,
                'rmsd': track.rmsd,
                #'rayleigh_statistics': track.rayleigh,
                'roll_correction': track.roll_correction,
                'vig_name': "".join((
                    self.deployment, self.depth, '_', str(track.img_names[0]),
                    '-', str(track.id), '.png'))
                }
            
            for key in list_keys:
                
                med = np.median([particle[key] for particle in track.particles])                
                dict_track[key] = med

            list_dicts.append(dict_track)
        
        # Saving dataframe summarizing the tracks
        df_summary = pd.DataFrame(list_dicts)
        df_summary.insert(0, 'Depth', self.depth)
        df_summary.insert(0, 'Cycle', self.deployment)
        filename = "".join(('tracks_summary_', self.deployment, '_',
                                self.depth, '_', self.subdir, '.tsv'))
        df_summary.to_csv(filename, sep = "\t", index=False)

        return


    def save_df_all(self):
        """ Convert a track object to a list of dictionary containing its
        particles properties. """
    
        list_dicts = []
        
        for track in self.list_tracks:
            
            track_id = "".join(
                (self. deployment, self.depth, '_',  str(track.img_names[0]), 
                 '-', str(track.id)))
            
            for particle in track.particles:
                dictpart = {
                    'track_id': track_id,
                    'track_length': track.length,
                    'speed': track.speed,
                    'angle': track.mean_angle,
                    'vertical_speed': track.vertical_speed,
                    'orientation': track.orientation
                    }
                list_dicts.append({**dictpart, **particle})

        # Saving dataframe of the tracks
        filename = "".join(('tracks_all_', self.deployment, '_', 
                                self.depth, '_', self.subdir, '.tsv'))
        pd.DataFrame(list_dicts).to_csv(filename, sep = "\t", index = False)

        return


    def plot_all_tracks(self):
        """
        Plots on a single figure all the tracks of a sequence, for each
        type of orientation.

        """
        
        orientations = ['desc', 'asc', 'mix']
        
        for orient in orientations:
            
            list_tr = [track for track in self.list_tracks 
                       if track.orientation==orient]

            if not list_tr:
                continue

            create_window()
            
            for track in list_tr:
                pos_x = [x[0] for x in track.coords] 
                pos_x = [x * 73 * 1e-3 for x in pos_x]
                pos_y = [x[1] for x in track.coords]
                pos_y = [y * 73 * 1e-3 for y in pos_y]
                plt.plot(pos_x, pos_y)

            # Adding title
            dt1, dt2 = list_tr[0].img_names[0], list_tr[-1].img_names[-1]
            plt.title('From {dt1} to {dt2}'.format(dt1 = dt1, dt2 = dt2))
            nb = len(list_tr)
            supttl = '{nb} detected {orient} tracks'
            plt.suptitle(supttl.format(nb = nb, orient = orient), y = 1)
    
            filename = "".join(('tracks_', orient, '.png'))
            plt.savefig(filename, facecolor = 'white',) 
            plt.close()
                
        return


    def run(self):
        """ Saving dataframes and plots. """
        
        # Saving dataframes tracks
        self.save_df_summary()
        self.save_df_all()

        # Plotting all tracks
        self.plot_all_tracks()
        
        return

