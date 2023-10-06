# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Track particles.
23-08-2022
@author: Manon Laget
------------------------------
"""


import numpy as np
from scipy.stats import circstd
#from astropy.stats.circstats import rayleightest


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
        return(a < n and n < b)
    
    return(a < n or n < b)


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
    calc_orientation
    calc_speed
    calc_vertical_speed
    calc_vx
    calc_sin_index
    calc_rmsd
    calc_longest_distance
    update
    end
    
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
        

    def calc_angle(self):
        """
        Calculates the mean and std of the direction (angle) in 
        degrees.
        
        Returns
        -------
        alpha : float
            Mean angle of the particle, in degrees.
        alpha_std : float
            Std of the angle, in degrees.
        
        """
        
        list_rad = self.angles_rad
        
        if not list_rad:
            return(None, None)

        # Mean angle
        n = len(list_rad)
        sin_alpha = np.sum(np.sin(list_rad))
        cos_alpha = np.sum(np.cos(list_rad))
        sin_alpha, cos_alpha = sin_alpha/n, cos_alpha/n
        alpha = np.arctan2(sin_alpha, cos_alpha)
        alpha = 360 - ((alpha*180/np.pi)%360)
        alpha_std = circstd(list_rad)*180/np.pi

        return(alpha, alpha_std)


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
        #/!\ The signs are inverted because of matrix coordinates,
        #    - means ascending and + descending
        list_dy = [crds[i+1][1] - crds[i][1] for i in range(len(crds)-1)]

        # Sign of dy between points
        test_sign = (np.sign(list_dy) == 1)

        if all(i for i in test_sign) is True:
            return('desc')
        elif any(i for i in test_sign) is True:
            return('mix')
        else:
            return('asc')


    def calc_speed(self):
        """ 
        Calculates the mean speed between points of a track.
        
        Returns
        -------
        speed of the particle over all the track, in m/day
        
        """
        
        crds = self.coords
        px_to_um = 73

        # Calculate dt
        dt = (self.datetimes[-1] - self.datetimes[0]).total_seconds()
         
        # Calculate speed
        distance = np.sum([((crds[i+1][0] - crds[i][0]) ** 2 + 
                            (crds[i+1][1] - crds[i][1]) ** 2) ** 0.5
                           for i in range(len(crds)-1)]) * px_to_um * (1e-6)

        return(86400*distance/dt) # m/day


    def calc_vertical_speed(self):
        """
        Calculates the mean speed on the vertical (y) axis of a track.
        
        Returns
        -------
        speed of the particle on the y-axis, in m/day
        """
        
        px_to_um = 73
        
        # Calculate absolute distance
        dy = (self.coords[-1][1] - self.coords[0][1])
        dx = (self.coords[-1][0] - self.coords[0][0])
        dist_abs = np.sqrt(dx**2+dy**2) * px_to_um * 1e-6
        
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
        
        return(86400*dist_vert/dt) # m/day
    
    
    def calc_vx(self):
        """
        Calculates the mean speed on the horizontal (x) axis of a track.
        
        Returns
        -------
        speed in m/day
        
        """
        
        crds = self.coords
        px_to_um = 73
        distance = np.sum([((crds[i+1][0] - crds[i][0])**2)**0.5
                           for i in range(len(crds)-1)]) * px_to_um * (1e-6)
        
        # Calculate dt
        dt = (self.datetimes[-1] - self.datetimes[0]).total_seconds()

        return(86400*distance/dt) # m/day
    
    
    def calc_sin_index(self):
        """ Calculates the sinusoity index of a track, defined as the ratio
        between effective track length and the length of the segment
        beginning-end. """
        
        crds = self.coords
        
        distance = np.sum([((crds[i+1][0] - crds[i][0])**2 + 
                             (crds[i+1][1] - crds[i][1])**2)**0.5
                            for i in range(len(crds)-1)])
        
        # Calculate length of beginning-end segment
        distance2 = ((crds[-1][0] - crds[1][0])**2 
                     + (crds[-1][1] - crds[1][1])**2)**0.5
        
        return(distance/distance2)
    
    
    def calc_rmsd(self):
        """ Calculates root mean squared displacement of the track. """
        
        crds = self.coords
        
        xdata = np.asarray([crds[i][0] for i in range(len(crds))])
        ydata = np.asarray([crds[i][1] for i in range(len(crds))])
        
        r = np.sqrt(xdata**2 + ydata**2)
        diff = np.diff(r) #this calculates r(t + dt) - r(t)
        diff_sq = diff**2
        rmsd = np.sqrt(np.mean(diff_sq))
        
        return(rmsd)
    
    
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
            d = np.cross(p_fin-p_ini,pi-p_ini)/np.linalg.norm(p_fin-p_ini)
            if d>longest_dist:
                longest_dist = d
        
        return(longest_dist)
    
    
    def update(self, particle):
        """
        Adds a new particle to a track.
        
        Parameters
        ----------
        particle : dictionary
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
        
        # Mean esd 
        self.mean_esd = np.mean([p['esd_um'] for p in self.particles])
        
        # Calculate dx and dy
        dy = (self.coords[-1][1] - self.coords[-2][1])
        dx = (self.coords[-1][0] - self.coords[-2][0])
        
        # Calculate the angle
        alpha = np.arctan2(dy, dx)
        self.angles_rad.append(alpha)
        self.angles_deg.append(360 - ((alpha*180/np.pi)%360))
        self.mean_angle, self.std_angle = self.calc_angle()
        
        # Updating step list
        self.steps.append(np.sqrt(dx**2+dy**2) * 73)
        # Updating mean_step
        self.mean_step = np.mean(self.steps)
        self.std_step = np.std(self.steps)
        
        return(self)


    def end(self):
        """
        Ends a track.
        
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
        #self.rayleigh = rayleightest(np.array(self.angles_rad))

        return(self)


class Tracking:
    """
    Runs tracking on a set of particles.
    
    Attributes
    ----------
    px_to_um : float
        Conversion factor from px to µm.
    min_length : integer, optional
        Minimum length to keep a track.
        
    Methods
    -------
    is_same_particle
    linking
    run
    
    """


    def __init__(self, df_particles, angle):
        
        self.df_particles = df_particles
        self.px_to_um = 73
        self.angle = angle
        
        # Define lists to store running and finished tracks
        self.all_tracks = []
    

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
        
        # Check for the size
        if particle['esd_um']<(0.8*track.mean_esd):
            return(False)
        elif particle['esd_um']>(1.2*track.mean_esd):
            return(False)

        # Set a max distance
        max_dist = track.mean_esd * 25
        if max_dist>30000:
            max_dist = 30000
        
        x1, y1 = track.coords[-1][0], track.coords[-1][1]
        x2, y2 = particle['coord_x'], particle['coord_y']

        dx, dy = (x2-x1)*self.px_to_um, (y2-y1)*self.px_to_um
        
        # Angle in radians
        angle = np.arctan2(dy, dx)
        angle = 360 - ((angle*180/np.pi)%360)
        
        dx, dy = dx**2, dy**2

        # Checking if the step length is >max_dist
        if np.sqrt(dx + dy) > max_dist:
            return(False)

        # If we have at least one step
        
        if track.mean_step and track.mean_angle:
            
            # Check for angle
            
            lim_a = (track.angles_deg[-1] % 360 - self.angle) % 360
            lim_b = (track.angles_deg[-1] % 360 + self.angle) % 360
            
            if not is_angle_between(angle, lim_a, lim_b):
                return(False)

            elif np.sqrt(dx + dy)>(2*track.mean_step):
                return(False)

        return(True)


    def linking(self, list_particles, running_tracks):
        """
        Associates particles on an image to a set of tracks.
        
        Parameters
        ----------
        list_particles : list of dictionaries
            Contains particles detected on an image.
        running_tracks : list of track objects
            Contains current unfinished tracks
        
        Returns
        -------
        new_running_tracks : list of track objects
            Contains new unfinished tracks
        
        """
        
        assigned_particles, new_running_tracks = ([] for i in range(2))
        
        for track in running_tracks:
            
            good_particle = None
            score_max = 1e6
            mean_esd = np.mean([p['esd_um'] for p in track.particles])
    
            for particle in list_particles:
                
                if self.is_same_particle(track, particle):
                    # Now calculate the size difference
                    score = np.abs(mean_esd-particle['esd_um'])
                    if score < score_max:
                        score_max, good_particle = score, particle
                        
            # Assign the most similar to the track
            if good_particle:
                new_running_tracks.append(track.update(good_particle))
                assigned_particles.append(good_particle)
                
            #If good_particle is still none then we end the track
            elif (track.length >= 4 and track.mean_esd >= 600):
                self.all_tracks.append(track.end())
            #If track is not long enough the last particle becomes a new track
            elif track.length>1: 
                # Append to running_track so it'll join the loop
                running_tracks.append(Track(track.particles[-1]))
            
        # Create new tracks with unassigned particles
        # Append them to list of running_tracks
        list_ids = [part['part_id'] for part in assigned_particles]
        
        for particle in list_particles:
            if particle['part_id'] not in list_ids:
                new_running_tracks.append(Track(particle))
        
        return(new_running_tracks)
        

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
            self.df_particles['n_img'] == ids_img[0]].to_dict('records')
        running_tracks = [Track(particle) for particle in particles]
        
        # 2. Loop over the images
        for i in range(1, len(ids_img)):
            particles = self.df_particles.loc[
                self.df_particles['n_img'] == ids_img[i]].to_dict('records')
            running_tracks = self.linking(particles, running_tracks)
        
        # 3. End
        self.all_tracks.extend(
            [track.end() for track in running_tracks 
             if (track.length >= 4 and track.mean_esd >= 600)])

        return(self.all_tracks)
    