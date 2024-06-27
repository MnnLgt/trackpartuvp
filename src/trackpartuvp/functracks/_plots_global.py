# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Save various plots for the whole project
2022
@author: Manon Laget
------------------------------
"""


import matplotlib.pyplot as plt


class PlotsGlobal:
    """
    Various plots for the whole deployment.

    Attributes
    ----------
    lengths
    angles
    speeds
    orientations
    deployment
    depth

    Methods
    -------
    hist_lengths
    polar_speeds
    hist_orientations

    """

    def __init__(self, lengths, angles, speeds, orientations,
                 deployment, depth):

        self.lengths = lengths
        self.angles = angles
        self.speeds = speeds
        self.orientations = orientations
        self.deployment = deployment
        self.depth = depth

    def hist_lengths(self):
        """ Histogram of the lengths of the tracks. """

        # Histogram of track lengths
        plt.hist(self.lengths)
        plt.xlabel('Length')
        plt.ylabel('Counts')
        plt.title('Number of tracks per length')
        filename = "".join(
            (self.deployment, self.depth, '_hist_lengths.png'))
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def polar_speeds(self):
        """ Polar plot of speeds as a function of angles. """

        # Polar plot of track angles
        plt.polar(self.angles, self.speeds, '.',
                  markersize=3)
        plt.title('Speed as a function of angle')
        filename = "".join(
            (self.deployment, self.depth, '_polar.png'))
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def hist_orientations(self):
        """ Histogram of orientation types. """

        # Proportion of asc, desc and mix tracks
        plt.hist(self.orientations)
        plt.xlabel('Orientation')
        plt.ylabel('Counts')
        plt.title('Number of tracks per orientation type')
        filename = "".join(
            (self.deployment, self.depth, '_hist_orientations.png'))
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def hist_speeds(self):
        """ Histogram of speeds. """

        # Proportion of asc, desc and mix tracks
        plt.hist(self.speeds)
        plt.xlabel('Speeds')
        plt.ylabel('Counts')
        plt.title('Histogram of speeds')
        filename = "".join(
            (self.deployment, self.depth, '_hist_speeds.png'))
        plt.savefig(filename, facecolor='white')
        plt.close()

        return

    def run(self):
        """ Runs all the plots. """

        plt.close()
        self.hist_lengths()
        self.polar_speeds()
        self.hist_orientations()
        self.hist_speeds()
        plt.close()

        return
