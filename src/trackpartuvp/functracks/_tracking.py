# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Particle tracking functions
2022
@author: Manon Laget, Alexandre Accardo
------------------------------
"""


import numpy as np
from scipy.stats import circstd
# from astropy.stats.circstats import rayleightest
import pandas as pd
from scipy.spatial import distance


def is_angle_between(n, a, b):
    """
    Returns true if an angle (in degrees) is between two angles (on the
    trigonometric circle), else false.

    Parameters
    ----------
    n: float
        Angle in degrees to test
    a: float
        Lower limit of the circle section (in degrees)
    b: float
        Upper limit of the circle section (in degrees)

    Returns
    -------
    bool
        True if angle is in range, else false.

    """

    n = (360 + (n % 360)) % 360
    a = (3600000 + a) % 360
    b = (3600000 + b) % 360

    if (a < b):
        return (a < n and n < b)

    return (a < n or n < b)


class Track:
    """
    Contains a track's properties.

    Attributes
    ----------
    particle: a dictionary contaning a particle's properties
        Particle that initiates the track.

    Methods
    -------
    calc_angle
        Calculates the mean and std of the track direction (angle) in degrees.
    calc_orientation
        Calculates the sign on the y-axis between points of a track and
        determines if a track goes descending, ascending or a mix.
    calc_speed
        Calculates the mean speed between points of a track.
    calc_vertical_speed
        Calculates the mean speed on the vertical (y) axis of a track.
    calc_vx
        Calculates the mean speed on the horizontal (x) axis of a track.
    calc_sin_index
        Calculates the sinusoity index of a track, defined as the ratio
        between effective track length and the length of the segment
        beginning-end.
    calc_rmsd
        Calculates root mean squared displacement of the track.
    calc_longest_distance
        Calculate the longest distance between two points.
    update
        Add a new particle to a track.
    end
        End a track.

    """

    def __init__(self, particle):

        # ID of the track given by the ID of its first particle
        self.id = particle['part_id']

        # Initial length
        self.length = 1

        # List of coordinates
        self.coords = [(particle['coord_x'], particle['coord_y'])]

        # List of datetimes
        self.datetimes = [particle['datetime']]

        # List of filenames
        self.filenames = [particle['filename']]

        # List of image names
        self.img_names = [particle['img_name']]

        # List of particle dictionaries
        self.particles = [particle]

        # Various properties of the track set to None
        self.mean_esd = particle['esd_um']
        self.track_angle = None
        self.steps = []
        self.angles_deg = []
        self.angles_rad = []
        self.mean_step = None
        self.std_step = None
        self.mean_angle = None
        self.std_angle = None
        self.speed = None
        self.vertical_speed = None
        self.vx = None
        self.orientation = None
        self.sinusoity_index = None
        self.longest_dist = None
        self.rmsd = None
        self.roll_correction = None
        self.last_particle = []

        # Set pixel length
        self.px_to_um = 73

    def calc_angle(self):
        """
        Calculates the mean and std of the track direction (angle) in degrees.

        Returns
        -------
        alpha: float
            Mean angle of the particle track, in degrees.
        alpha_std: float
            Std of the angle, in degrees.

        """

        list_rad = self.angles_rad

        if not list_rad:
            return (None, None)

        # Mean angle
        n = len(list_rad)
        sin_alpha = np.sum(np.sin(list_rad))
        cos_alpha = np.sum(np.cos(list_rad))
        sin_alpha, cos_alpha = sin_alpha/n, cos_alpha/n
        alpha = np.arctan2(sin_alpha, cos_alpha)
        alpha = 360 - ((alpha*180/np.pi) % 360)
        alpha_std = circstd(list_rad)*180/np.pi

        return (alpha, alpha_std)

    # def calc_track_angle(self):
    #    """
    #    Calculates the mean angle between consecutive particles in degrees.
    #
    #    Returns
    #    -------
    #    track_angle : float
    #        Mean angle between consecutive particles in degrees.
    #    """
    #    if len(self.coords) < 2:
    #        return None
    #
    #    # Calculate angle between consecutive particles
    #    angles = [np.arctan2(self.coords[i + 1][1] - self.coords[i][1],
    #                         self.coords[i + 1][0] - self.coords[i][0])
    #              for i in range(len(self.coords) - 1)]
    #
    #    # Convert angles to degrees
    #    angles_deg = [360 - ((angle * 180 / np.pi) % 360) for angle in angles]
    #
    #    # Calculate mean angle
    #    # /!\ CANNOT CALCULATE MEAN OF CIRCULAR DATA LIKE THAT
    #    mean_angle = np.mean(angles_deg)
    #
    #    return mean_angle

    def calc_orientation(self):
        """
        Calculates the sign on the y-axis between points of a track and
        determines if a track goes descending, ascending or a mix.

        Returns
        -------
        orientation: str
            string indicating the orientation
            'a' for ascending
            'd' for descending
            'm' for mix

        """

        crds = self.coords

        # Calculate dy between points
        # /!\ The signs are inverted because of matrix coordinates,
        #    - means ascending and + descending
        list_dy = [crds[i+1][1] - crds[i][1] for i in range(len(crds)-1)]

        # Sign of dy between points
        test_sign = (np.sign(list_dy) == 1)

        if all(i for i in test_sign) is True:
            return 'desc'
        elif any(i for i in test_sign) is True:
            return 'mix'
        else:
            return 'asc'

    def calc_speed(self):
        """
        Calculates the mean speed between points of a track.

        Returns
        -------
        speed of the particle over all the track, in m/day

        """

        crds = self.coords

        # Calculate dt
        dt = (self.datetimes[-1] - self.datetimes[0]).total_seconds()

        # Calculate speed
        distance = np.sum(
            [((crds[i+1][0] - crds[i][0]) ** 2 + (
                crds[i+1][1] - crds[i][1]) ** 2) ** 0.5
                for i in range(len(crds)-1)]) * self.px_to_um * (1e-6)

        return (86400*distance/dt)  # m/day

    def calc_vertical_speed(self):
        """
        Calculates the mean speed on the vertical (y) axis of a track.

        Returns
        -------
        speed of the particle on the y-axis, in m/day
        """

        # Calculate absolute distance
        dy = (self.coords[-1][1] - self.coords[0][1])
        dx = (self.coords[-1][0] - self.coords[0][0])
        dist_abs = np.sqrt(dx**2+dy**2) * self.px_to_um * 1e-6

        # Calculate vertical distance
        if self.mean_angle < 90:
            angle = 90 - self.mean_angle
        elif self.mean_angle < 180:
            angle = self.mean_angle - 90
        elif self.mean_angle < 270:
            angle = 270 - self.mean_angle
        elif self.mean_angle < 360:
            angle = self.mean_angle - 270

        angle = angle*np.pi/180
        dist_vert = dist_abs * np.cos(angle)

        # Calculate dt
        dt = (self.datetimes[-1] - self.datetimes[0]).total_seconds()

        return (86400*dist_vert/dt)  # m/day

    def calc_vx(self):
        """
        Calculates the mean speed on the horizontal (x) axis of a track.

        Returns
        -------
        speed in m/day

        """

        crds = self.coords
        distance = np.sum(
            [((crds[i+1][0] - crds[i][0])**2)**0.5
             for i in range(len(crds)-1)]) * self.px_to_um * (1e-6)

        # Calculate dt
        dt = (self.datetimes[-1] - self.datetimes[0]).total_seconds()

        return (86400*distance/dt)  # m/day

    def calc_sin_index(self):
        """ Calculates the sinusoity index of a track, defined as the ratio
        between effective track length and the length of the segment
        beginning-end. """

        crds = self.coords

        distance = np.sum(
            [((crds[i+1][0] - crds[i][0])**2 +
              (crds[i+1][1] - crds[i][1])**2)**0.5
             for i in range(len(crds)-1)])

        # Calculate length of beginning-end segment
        distance2 = ((crds[-1][0] - crds[1][0])**2
                     + (crds[-1][1] - crds[1][1])**2)**0.5

        return (distance/distance2)

    def calc_rmsd(self):
        """ Calculates root mean squared displacement of the track. """

        crds = self.coords

        xdata = np.asarray([crds[i][0] for i in range(len(crds))])
        ydata = np.asarray([crds[i][1] for i in range(len(crds))])

        r = np.sqrt(xdata**2 + ydata**2)
        diff = np.diff(r)  # this calculates r(t + dt) - r(t)
        diff_sq = diff**2
        rmsd = np.sqrt(np.mean(diff_sq))

        return rmsd

    def calc_longest_distance(self):
        """
        Calculate the longest distance between two points.

        """

        crds = self.coords

        p_ini = np.asarray(crds[0])
        p_fin = np.asarray(crds[-1])

        longest_dist = 0
        for crd in crds[1:-1]:
            pi = np.asarray(crd)
            d = np.cross(
                p_fin - p_ini, pi-p_ini) / np.linalg.norm(
                    p_fin - p_ini)
            if d > longest_dist:
                longest_dist = d

        return longest_dist

    def update(self, particle):
        """
        Add a new particle to a track.

        Parameters
        ----------
        particle: dictionary
            Particle to be added to the track.

        Returns
        -------
        The updated track.

        """

        # Increment length
        self.length += 1

        # Adding to list of coords, datetimes and particles
        self.coords.append((particle['coord_x'], particle['coord_y']))
        self.datetimes.append(particle['datetime'])
        self.particles.append(particle)
        self.filenames.append(particle['filename'])
        self.img_names.append(particle['img_name'])
        self.last_particle = particle['part_id']

        # Mean esd
        self.mean_esd = np.mean([p['esd_um'] for p in self.particles])

        # Calculate dx and dy
        dy = (self.coords[-1][1] - self.coords[-2][1])
        dx = (self.coords[-1][0] - self.coords[-2][0])

        # Calculate the angle
        alpha = np.arctan2(dy, dx)
        self.angles_rad.append(alpha)
        self.angles_deg.append(360 - ((alpha*180/np.pi) % 360))
        self.mean_angle, self.std_angle = self.calc_angle()
        self.track_angle = self.mean_angle

        # Updating step list
        self.steps.append(np.sqrt(dx**2+dy**2) * self.px_to_um)

        # Updating mean_step
        self.mean_step = np.mean(self.steps)
        self.std_step = np.std(self.steps)

        return self

    def end(self):
        """
        End a track.

        Returns
        -------
        The finished track.

        """

        rolls = [p['roll'] for p in self.particles]
        rolls = [roll for roll in rolls if roll is not None]
        self.roll_correction = np.mean(rolls) if rolls else 0

        # Correct angle
        self.mean_angle = (self.mean_angle + self.roll_correction) % 360

        # Calculate various values
        self.speed = self.calc_speed()
        self.vertical_speed = self.calc_vertical_speed()
        self.vx = self.calc_vx()
        self.orientation = self.calc_orientation()
        self.sinusoity_index = self.calc_sin_index()
        self.rmsd = self.calc_rmsd()
        self.longest_dist = self.calc_longest_distance()
        # self.rayleigh = rayleightest(np.array(self.angles_rad))

        return self


class Tracking:
    """
    Runs tracking on a set of particles.

    Attributes
    ----------
    df_particles: panda dataframe
        Dataframe containing all particles detected on a sequence.
    angle: float
        Angle to set research area.
    platform_speed_path: str, optional
        Path to dataframe with speeds of the platform.
    px_to_um: float, optional
        Conversion factor from px to µm.
    min_length: integer, optional
        Minimum length to keep a track.

    Methods
    -------
    load_platform_speed_data
    is_in_ellipse
    is_in_rotated_ellipse
    is_same_particle
    is_same_particle_in_next_frame
    linking
    run

    """

    def __init__(self, df_particles, angle, platform_speed_path, px_to_um=73):

        self.df_particles = df_particles
        self.px_to_um = px_to_um
        self.angle = angle
        self.platform_speed_path = platform_speed_path  # AA
        self.df_particles['potential_match'] = False

        # Define lists to store running and finished tracks
        self.all_tracks = []
        # self.platform_speed_dict = self.load_platform_speed_data() # AA

    def load_platform_speed_data(self):  # AA
        """
        Create a dictionary with the the float speed associated with the
        sequence (used in is_same_particle).

        """

        # Read the CSV file containing sequence-specific float speed values
        speed_data = pd.read_csv(self.platform_speed_path)

        return dict(zip(speed_data['seq'],
                        speed_data['float_vertical_speed']))

    def is_in_ellipse(self, x, y, x2, y2, minor_axis, major_axis):  # AA
        """
        Test if coordinates of a new particle is inside the research ellipse.
        The research ellipse doesn't have orientation (initialisation).

        Parameters
        ----------
        x, y: float
            Detected particle coordinates at t(n).
        x2, y2: float
            Potential particle coordinates at t(n+1).
        minor_axis: float
            Minor axis length of the ellipse.
        major_axis: float
            Major axis length of the ellipse.

        Returns
        -------
        Whether the particle is in ellipse or no.

        """

        if ((np.square((x2 - x) * self.px_to_um) / np.square(major_axis)) +
            (np.square((y2 - y) * self.px_to_um) / np.square(minor_axis))) > 1:
            return ('Not in the ellipse')
        else:
            return ('In the ellipse')

    def is_in_rotated_ellipse(self, x, y, x2, y2, major_axis, minor_axis,
                              track_angle_degrees):  # AA
        """
        Test if coordinates of a new particle is inside the research ellipse.
        Only for the post-initialisation phase, when a track is already
        composed of minimun two particles. Ellipse with an orientation
        depending on the angle of the track considered.

        Parameters
        ----------
        x, y: float
            Detected particle coordinates at t(n).
        x2, y2: float
            Potential particle coordinates at t(n+1).
        minor_axis: float
            Minor axis length of the ellipse.
        major_axis: float
            Major axis length of the ellipse.
        track_angle_degrees: float
            Mean orientation of the track (in degrees).

        Returns
        -------
        Whether the particle is in ellipse or no.

        """

        if track_angle_degrees is None:
            track_angle_degrees = 0
        else:
            track_angle_degrees = track_angle_degrees

        # print('track angle: ', track_angle_degrees)
        angle_rad = np.radians(track_angle_degrees)
        # print('ellipse angle in deg: ', track_angle_degrees)
        # print('ellispe angle in rad: ', angle_rad)

        term1_numerator = np.square(
            ((x2 - x) * self.px_to_um) * np.cos(2*np.pi - angle_rad) -
            ((y2 - y) * self.px_to_um) * np.sin(2*np.pi - angle_rad))

        term1_denominator = np.square(minor_axis)

        term2_numerator = np.square(
            ((x2 - x) * self.px_to_um) * np.sin(2*np.pi - angle_rad) +
            ((y2 - y) * self.px_to_um) * np.cos(2*np.pi - angle_rad))

        term2_denominator = np.square(major_axis)

        if (term1_numerator / term1_denominator) + (
                term2_numerator / term2_denominator) > 1:
            return False
        else:
            return True

    def is_same_particle(self, track, particle):
        """
        Check if particle is within search area.

        Parameters
        ----------
        track : track object
        particle : pandas dataframe
            Dataframe row containing a particle's properties

        Returns
        -------
        bool
            Whether a particle is within the search range of a track.

        """

        # if particle['esd_um']>600:
        #    if particle['esd_um']<(0.8*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)
        #    elif particle['esd_um']>(1.2*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)

        # else:
        #    if particle['esd_um']<(0.1*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)
        #    elif particle['esd_um']>(1.9*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)

        # Get the float speed for the current sequence
        current_seq = particle['sequence']
        # float_speed = self.float_speed_dict.get(current_seq, 'NA')
        # Replace 'NA' with 0
        # float_speed = float(float_speed) if float_speed != 'NA' else 0
        # float_speed = abs(float_speed)

        # Set up the research zone as an ellipse

        # major axis (parallel to the x-axis, i.e. horizontal component of
        # the ellipse)
        major_axis = track.mean_esd * 10  # dependant on the size of the
        # particle to track (Manon put the threshold to 25*esd
        # and Alexandre to 10*esd)

        # set a max distance
        if major_axis > 30000:
            major_axis = 30000

        # minor axis (parallel to the y-axis, i.e. vertical component of the
        # ellipse)

        # if float_speed == 0:
        #    minor_axis = 100000 # create a circle rather than an ellispe when
        #                        # the float is stable
        # else:
        #    minor_axis = 3*np.sqrt(float_speed)*major_axis # weight the minor
        #    # axis by the the float speed (maybe add a factor before to
        #    # intensify the influence of the float speed)
        #
        # set a max distance
        # if minor_axis > 100000:
        #    minor_axis = 100000

        minor_axis = major_axis  # create a circle if minor and major equal

        # print('major_axis: ', major_axis)
        # print('minor_axis: ', minor_axis)

        # Get coordinates of the last point in the track and the particle
        # center of the research ellipse
        x1, y1 = track.coords[-1][0], track.coords[-1][1]
        # coordinates to test if there are in the ellipse
        x2, y2 = particle['coord_x'], particle['coord_y']

        # Calculate the change in x and y, scaled by 'px_to_um'
        dx, dy = (x2-x1)*self.px_to_um, (y2-y1)*self.px_to_um

        # Angle in radians
        angle = np.arctan2(dy, dx)
        # convert radian to degree (*180/pi)
        angle = 360 - ((angle*180/np.pi) % 360)

        dx, dy = dx**2, dy**2

        # If we have at least one step

        # if track.angles_deg:
        # convert angle (into radians)
        #   track_angle_degrees = track.angles_deg[0]
        # else:
        #   track_angle_degrees = 0

        if not self.is_in_rotated_ellipse(
                x=x1, y=y1, x2=x2, y2=y2,
                major_axis=major_axis, minor_axis=minor_axis,
                track_angle_degrees=track.track_angle):
            # print('Not in the ellipse (is_same_particle)')
            return False

        # Check for angle
        if track.track_angle and track.mean_step:

            lim_a = ((track.track_angle % 360) - self.angle) % 360
            lim_b = ((track.track_angle % 360) + self.angle) % 360
            # print('lim_a: ', lim_a)
            # print('lim_b: ', lim_b)

            if not is_angle_between(angle, lim_a, lim_b):
                # print('not in angle (is_same_particle)')
                return False

            elif np.sqrt(dx + dy) > (4*track.mean_step):
                # print('too far from particle (is_same_particle)')
                return False

        return True

    def is_same_particle_in_next_frame(self, track, particle):
        """
        Check if particle is within search area.

        Parameters
        ----------
        track : track object
        particle : pandas dataframe
            Dataframe row containing a particle's properties

        Returns
        -------
        bool
            Whether a particle is within the search range of a track.

        """

        # if particle['esd_um']>600:
        #    if particle['esd_um']<(0.8*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)
        #    elif particle['esd_um']>(1.2*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)

        # else:
        #    if particle['esd_um']<(0.1*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)
        #    elif particle['esd_um']>(1.9*track.mean_esd):
        #        print('particle size too different (is_same_particle)')
        #        return(False)

        # Get the float speed for the current sequence

        current_seq = particle['sequence']
        # float_speed = self.float_speed_dict.get(current_seq, 'NA')
        # Replace 'NA' with 0
        # float_speed = float(float_speed) if float_speed != 'NA' else 0
        # float_speed = abs(float_speed)

        # Set up the research zone as an ellipse

        # major axis (parallel to the x-axis, i.e. horizontal component of
        # the ellipse)
        major_axis = track.mean_esd * 25  # dependant on the size of the
        # particle to track

        # set a max distance
        if major_axis > 30000:
            major_axis = 30000

        # minor axis (parallel to the y-axis, i.e. vertical component of the
        # ellipse)

        # if float_speed == 0:
        #    minor_axis = 100000 # create a circle rather than an ellispe when
        #                        # the float is stable
        # else:
        #    minor_axis = 3*np.sqrt(float_speed)*major_axis # weight the minor
        #    # axis by the the float speed (maybe add a factor before to
        #    # intensify the influence of the float speed)
        #
        # set a max distance
        # if minor_axis > 100000:
        #    minor_axis = 100000

        minor_axis = major_axis
        # print('major_axis: ', major_axis)
        # print('minor_axis: ', minor_axis)

        # Get coordinates of the last point in the track and the particle
        # center of the research ellipse
        x1, y1 = track.coords[-1][0], track.coords[-1][1]
        # coordinates to test if there are in the ellipse
        x2, y2 = particle['coord_x'], particle['coord_y']

        # Calculate the change in x and y, scaled by 'px_to_um'
        dx, dy = (x2-x1)*self.px_to_um, (y2-y1)*self.px_to_um

        # Angle in radians
        angle = np.arctan2(dy, dx)
        # convert radian to degree (*180/pi)
        angle = 360 - ((angle*180/np.pi) % 360)

        dx, dy = dx**2, dy**2

        # If we have at least one step

        # if track.angles_deg:
        # convert angle (into radians)
        #   track_angle_degrees = track.angles_deg[0]
        # else:
        #   track_angle_degrees = 0

        if not self.is_in_rotated_ellipse(
                x=x1, y=y1, x2=x2, y2=y2,
                major_axis=major_axis, minor_axis=minor_axis,
                track_angle_degrees=track.track_angle):
            # print('Not in the ellipse (is_same_particle)')
            return False

        # Check for angle
        if track.track_angle and track.mean_step:

            lim_a = ((track.track_angle % 360) - self.angle) % 360
            lim_b = ((track.track_angle % 360) + self.angle) % 360
            # print('lim_a: ', lim_a)
            # print('lim_b: ', lim_b)

            if not is_angle_between(angle, lim_a, lim_b):
                # print('not in angle (is_same_particle_next)')
                return False

            elif np.sqrt(dx + dy) > (4*track.mean_step):
                # print('too far from particle (is_same_particle_next)')
                return False

        return True

    def linking(self,
                # list of dict containing particles detected on the image n
                list_particles,
                list_particles_n2,
                # list of track objects, containing currrent unfinished tracks
                running_tracks
                ):
        """
        Associates particles on an image to a set of tracks.

        Parameters
        ----------
        list_particles : list of dictionaries
            Contains particles detected on an image.
        list_particles_n2 : list of dictionaries
            Contains particles detected on the image before.
        running_tracks : list of track objects
            Contains current unfinished tracks

        Returns
        -------
        new_running_tracks : list of track objects
            Contains new unfinished tracks

        """

        # create empty list
        assigned_particles, new_running_tracks = ([] for i in range(2))
        # assigned_particles -> track of particles that have been assigned to an existing track
        # new_running_tracks -> track of particles that have been assigned to a new track

        for track in running_tracks:

            good_particle = None
            euclidean_dist_max = 1e10  # set a high score to start

            # Create a vector with the mean values of each morphological
            # feature defined by the PCA
            track_PC1 = np.mean([p['PC1'] for p in track.particles])
            track_PC2 = np.mean([p['PC2'] for p in track.particles])
            track_PC3 = np.mean([p['PC3'] for p in track.particles])
            track_PC4 = np.mean([p['PC4'] for p in track.particles])
            track_PC5 = np.mean([p['PC5'] for p in track.particles])

            morpho_features_track = np.array(
                [track_PC1, track_PC2, track_PC3, track_PC4, track_PC5])

            track_particle = [p['part_id'] for p in track.particles]
            # print('particles in current track: ', track_particle)

            # Check if the particle is a potential match to the track
            for particle in list_particles:

                # print('particles in the current frame: ', particle)
                # rint('potential match in current frame? ',
                #      self.is_same_particle(track, particle), 'for ',
                #      particle['part_id'])

                # Apply the euclidean distance
                if particle['part_id'] not in track_particle and self.is_same_particle(track, particle):

                    PC1_part = particle['PC1']
                    PC2_part = particle['PC2']
                    PC3_part = particle['PC3']
                    PC4_part = particle['PC4']
                    PC5_part = particle['PC5']

                    morpho_features_particle = np.array(
                        [PC1_part, PC2_part, PC3_part, PC4_part, PC5_part])

                    euclidean_distance = distance.euclidean(
                        morpho_features_track, morpho_features_particle)

                    # score = np.abs(mean_esd-particle['esd_um'])
                    # if self.is_same_particle() is true, then it computes a
                    # score for the particle based on its size
                    # Compute the euclidean distance between the morpho
                    # features of this particle and the mean morpho of the
                    # current track

                    if euclidean_distance < euclidean_dist_max:
                        # TO DO
                        # add a minimum threshold here (95th percentile)

                        # update euclidean_dist_max and good_particle if
                        # euclidean_distance < euclidean_dist_max
                        euclidean_dist_max, good_particle = euclidean_distance, particle

            if good_particle:

                new_running_tracks.append(track.update(good_particle))
                assigned_particles.append(good_particle)

            elif len(list_particles_n2) != 0:

                # print('no potential match in current frame, check in img n+1')

                good_particle = None
                euclidean_dist_max = 1e10

                # checks if the particle is a potential match to the track
                for particle in list_particles_n2:

                    # print('particles in the next frame: ', particle)
                    # print('potential match in the next frame: ',
                    #      self.is_same_particle_in_next_frame(track, particle),
                    #      'for ', particle['part_id'])

                    # Apply the euclidean distance
                    if particle['part_id'] not in track_particle and self.is_same_particle_in_next_frame(track, particle):

                        PC1_part = particle['PC1']
                        PC2_part = particle['PC2']
                        PC3_part = particle['PC3']
                        PC4_part = particle['PC4']
                        PC5_part = particle['PC5']

                        morpho_features_particle = np.array(
                            [PC1_part, PC2_part, PC3_part, PC4_part, PC5_part])

                        euclidean_distance = distance.euclidean(
                            morpho_features_track, morpho_features_particle)

                        # score = np.abs(mean_esd-particle['esd_um'])
                        # if self.is_same_particle() is true, then it computes a
                        # score for the particle based on its size
                        # Compute the euclidean distance between the morpho
                        # features of this particle and the mean morpho of the
                        # current track

                        if euclidean_distance < euclidean_dist_max:
                            # TO DO
                            # add a minimum threshold here (95th percentile)

                            # update euclidean_dist_max and good_particle if
                            # euclidean_distance < euclidean_dist_max
                            euclidean_dist_max, good_particle = euclidean_distance, particle

                if good_particle:
                    # print('good_particle found in the next frame')
                    new_running_tracks.append(track.update(good_particle))
                    assigned_particles.append(good_particle)

                # elif track.length >= 4: #4
                elif (track.length >= 4 and track.mean_esd >= 100):  # 200
                    # print('no potential match in next frame and track long enough to be saved')
                    self.all_tracks.append(track.end())

                # If track is not long enough the last particle becomes a new
                # track
                elif track.length > 1:
                    # print('no potential match in next frame and track too short to be saved, so start a new one with particle')
                    # Append to running_track so it'll join the loop
                    running_tracks.append(Track(track.particles[-1]))

        # Create new tracks with unassigned particles
        # Append them to list of running_tracks

        # print('running track kk:', running_tracks)
        list_ids = [part['part_id'] for part in assigned_particles]

        for particle in list_particles:
            if particle['part_id'] not in list_ids:
                new_running_tracks.append(Track(particle))

        return new_running_tracks

    def run(self):
        """ Runs tracking on a set of particles.

        Returns
        -------
        all_tracks : list of Track objects
            List of detected tracks.

        """

        # 1. Initialization
        # Get image ids
        ids_img = sorted(self.df_particles['n_img'].unique())

        # Get the particles on the first image and convert to dict
        particles = self.df_particles.loc[
            self.df_particles['n_img'] == ids_img[0]
            ].to_dict('records')  # image number ONE

        running_tracks = [Track(particle) for particle in particles]
        # print('particles in the first image: ', particles)

        # print('last particle in the initial track: ',
        # (running_tracks.last_particle for running_tracks in running_tracks))

        # 2. Loop over the images
        for i in range(1, len(ids_img)):

            # print('----start----')
            particles_n1 = self.df_particles.loc[
                self.df_particles['n_img'] == ids_img[i]].to_dict('records')

            if i+1 in range(1, len(ids_img)-1):
                particles_n2 = self.df_particles.loc[
                    self.df_particles['n_img'] == ids_img[i+1]
                    ].to_dict('records')

            else:
                particles_n2 = []

            running_tracks = self.linking(
                particles_n1, particles_n2, running_tracks)
            # print('----end----')
            # print('  ')

        # 3. End
        self.all_tracks.extend(
            [track.end() for track in running_tracks
             if (track.length >= 4 and track.mean_esd >= 100)])  # 200

        return self.all_tracks
