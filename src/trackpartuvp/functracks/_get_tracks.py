# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Save tracks
2022
@author: Manon Laget, modified by Alexandre Accardo
------------------------------
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from trackpartuvp.utils._utils import create_window
from scipy.spatial import distance  # AA
import seaborn as sns


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

    def save_df_summary(self, tracks_to_keep):
        """ Save dataframe with summary of each track. """

        list_keys = [
            'area_px', 'esd_px', 'perim_px', 'area_um', 'esd_um', 'major_px',
            'minor_px', 'perimmajor', 'meangrey',
            'stdgrey', 'cvgrey', 'intgrey', 'mediangrey', 'mingrey', 'maxgrey',
            'rangegrey', 'area_convex', 'eccentricity',
            'extent', 'solidity']
        # 'skewgrey', 'kurtgrey', 'elongation', 'circularity'

        list_dicts = []

        # Loop over tracks
        for track in self.list_tracks:

            # Get unique ID for each track
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
            alpha = 360 - ((alpha*180/np.pi) % 360)

            # Define variables to put in the dataframe
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
                # 'rayleigh_statistics': track.rayleigh,
                'roll_correction': track.roll_correction,
                'vig_name': "".join((
                    self.deployment, self.depth, '_', str(track.img_names[0]),
                    '-', str(track.id), '.png'))
                }

            # For each morphological variable, get the median value
            for key in list_keys:
                med = np.median([particle[key]
                                 for particle in track.particles])
                dict_track[key] = med

            # Append track dictionary to the list
            list_dicts.append(dict_track)

        # Turn dict list to dataframe
        df_summary = pd.DataFrame(list_dicts)

        # Add columns with deployment name and depth
        df_summary.insert(0, 'Depth', self.depth)
        df_summary.insert(0, 'Cycle', self.deployment)

        summary = df_summary[df_summary['track_id'].isin(tracks_to_keep)]
        df_summary = summary.drop_duplicates(subset=['track_id'])

        # Saving dataframe summarizing the tracks
        filename = "".join(('tracks_summary_', self.deployment, '_',
                            self.depth, '_', self.subdir, '.tsv'))
        df_summary.to_csv(filename, sep="\t", index=False)

        # Load the dataframe saved by save_df_all
        # filename_all = "".join(('tracks_all_', self.deployment, '_',
        #                       self.depth, '_', self.subdir, '.tsv'))
        # df_all = pd.read_csv(filename_all, sep="\t")

        # Calculate mean and standard deviation of 'euclidean_distance' for
        # each unique track_id
        # euclidean_stats = df_all.groupby(
        # 'track_id')['euclidean_distance'].agg(['mean', 'std'])

        # Merge summary dataframe with the mean and standard deviation values
        # df_summary = pd.merge(
        # df_summary, euclidean_stats, left_on='track_id', right_index=True,
        # how='left')

        # Rename the new columns
        # df_summary.rename(
        # columns={'mean': 'euclidean_mean', 'std': 'euclidean_std'},
        # inplace=True)

        # save it
        # df_summary.to_csv(filename, sep="\t", index=False)

        return

    def euclidean_hist(self, track_particle, hist_path):  # AA
        """
        Hisogram of euclidean distances.
        TO DO -> only display 10^-3 distance

        Parameters
        ----------
        track_particle: str
            dataframe that contains the morphological features of particle
            belonging to a track (i.e. same track_id).
        hist_path: str
            Path to save the histogram.

        Returns
        -------
        None.

        """

        g = sns.FacetGrid(track_particle, col="track_id", col_wrap=3, height=4)
        g.map(sns.histplot, "euclidean_distance", bins=10, kde=True,
              color="black")
        g.set_axis_labels("Euclidean Distance", "Density")
        # Save the plot as a PNG file
        plt.savefig(hist_path)

        return

    def save_df_all(self):  # AA
        """
        Convert a track object to a list of dictionary containing its
        particles properties

        Compute the euclidean distance between particle of a same track

        Create a density plot of euclidean distance distribution for each track
        """

        list_dicts = []
        list_keys = [
            'area_px', 'esd_px', 'perim_px', 'area_um', 'esd_um', 'major_px',
            'minor_px', 'perimmajor', 'meangrey',
            'stdgrey', 'cvgrey', 'intgrey', 'mediangrey', 'mingrey', 'maxgrey',
            'rangegrey', 'area_convex', 'eccentricity',
            'extent', 'solidity']
        # 'skewgrey', 'kurtgrey', 'elongation', 'circularity'

        for track in self.list_tracks:

            track_id = "".join(
                (self. deployment, self.depth, '_',  str(track.img_names[0]),
                 '-', str(track.id)))

            morpho_features_track = None

            for i, particle in enumerate(track.particles):
                dictpart = {
                    'track_id': track_id,
                    'track_length': track.length,
                    'speed': track.speed,
                    'angle': track.mean_angle,
                    'vertical_speed': track.vertical_speed,
                    'orientation': track.orientation
                    }

                # Add particle properties to the dictionary
                dictpart.update(particle)
                list_dicts.append(dictpart)

                # Compute the euclidean distance for particle in the same track
                if i > 0:
                    if morpho_features_track is not None:
                        morpho_features_particle = np.array(
                            [particle[key] for key in list_keys])
                        euclidean_distance = distance.euclidean(
                            morpho_features_track, morpho_features_particle)
                        dictpart['euclidean_distance'] = euclidean_distance

                # update morpho features track for the next iteration
                morpho_features_track = np.array(
                    [particle[key] for key in list_keys])

        # Saving dataframe of the tracks
        filename = "".join(('tracks_all_', self.deployment, '_',
                            self.depth, '_', self.subdir, '.tsv'))
        list_dicts = (pd.DataFrame(list_dicts)).fillna(0)

        tracks = list_dicts
        # Group by 'track_id' and find the minimum 'datetime' for each group
        min_date_per_track = list_dicts.groupby(
            ['sequence', 'track_id'])['datetime'].min().reset_index()

        # Merge the result back to the original DataFrame
        list_dicts = pd.merge(
            list_dicts, min_date_per_track, on='track_id', how='left',
            suffixes=('', '_ini'))

        list_dicts = list_dicts.rename(columns={'datetime_ini': 'date_ini'})

        # Sort tracks in chronological order
        list_dicts = list_dicts.sort_values(
            by=['date_ini', 'track_id', 'part_id'])
        print('list_dicts: ', list_dicts)

        # Identify tracks that have common particles
        common_particles = list_dicts.groupby(
            ['sequence', 'part_id'])['track_id'].apply(set).reset_index()
        common_particles = common_particles[
            common_particles['track_id'].apply(len) > 1]

        print('common particles step 1: ', common_particles)

        # remove duplicates only if there are some duplicates
        if len(common_particles) != 0:

            # Keep the track with the highest length

            def get_smaller_track(row):
                tracks = row['track_id']
                track_lengths = list_dicts[
                    list_dicts['track_id'].isin(tracks)]['track_length']
                return track_lengths.idxmin()

            common_particles['track_to_remove'] = common_particles.apply(
                get_smaller_track, axis=1)
            common_particles['track_to_remove'] = common_particles[
                'track_to_remove'].map(list_dicts['track_id'])

            # Remove duplicate tracks from original dataset
            list_dicts_without_dupl = list_dicts[
                ~list_dicts['track_id'].isin(
                    common_particles['track_to_remove'])]

        else:  # if there are no duplicates
            list_dicts_without_dupl = list_dicts

        # extract the tracks to keep
        tracks_to_keep = list_dicts_without_dupl['track_id'].unique()

        list_dicts_without_dupl.to_csv(filename, sep="\t", index=False)

        # Plot the density distribution of euclidean distance for each track
        # hist_path = "".join(('tracks_density', self.subdir, '.png'))
        # self.euclidean_hist(track_particle = pd.DataFrame(list_dicts),
        # hist_path = hist_path)

        # TO DO -> find a distance threshold to stop the trakc because at this
        # time I select the particle with the lowest euclidean distance within
        # the research ellipse

        return list_dicts_without_dupl, tracks_to_keep

    def plot_all_tracks(self, df_all):
        """
        Plots on a single figure all the tracks of a sequence, for each
        type of orientation.

        """

        orientations = ['desc', 'asc', 'mix']

        for orient in orientations:
            '''
            list_tr = [track for track in self.list_tracks
                       if track.orientation==orient]

            if not list_tr:
                continue
            '''

            list_tr = df_all[
                df_all['orientation'] == orient]['track_id'].unique()

            colors = plt.cm.rainbow(np.linspace(0, 1, len(list_tr)))

            create_window()

            for i, track_id in enumerate(list_tr):
                track_data = df_all[(
                    df_all['track_id'] == track_id) & (
                        df_all['orientation'] == orient)]
                plt.plot(track_data['coord_x'], track_data['coord_y'],
                         color=colors[i])

            # Adding title
            dt1, dt2 = df_all.iloc[0]['img_name'], df_all.iloc[-1]['img_name']
            plt.title('From {dt1} to {dt2}'.format(dt1=dt1, dt2=dt2))
            nb = len(list_tr)
            supttl = '{nb} detected {orient} tracks'
            plt.suptitle(supttl.format(nb=nb, orient=orient), y=1)

            filename = "".join(('tracks_', orient, '.png'))
            plt.savefig(filename, facecolor='white',)
            plt.close()

        return

    def run(self):
        """ Saving dataframes and plots. """

        # Saving dataframes tracks
        df_all, tracks_to_keep = self.save_df_all()
        self.save_df_summary(tracks_to_keep)

        # Plotting all tracks
        self.plot_all_tracks(df_all)

        return
