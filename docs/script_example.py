# -*- coding: utf-8 -*-
"""
------------------------------
UVP6 Particle tracking
------------------------------
Pipeline to track detected particles and save tracks
2022
@author: Manon Laget, Alexandre Accardo
------------------------------
"""






"""


Before initiate tracking, make sure that the data are well structured.
The algorithm needs sequences of png file. For each sequence, png files have to be stacked into a zipfile called 'images.zip'.
Each sequence needs to be organized in individual folder. They need to be named with the date-time of the sequence acquisiton (for example : '20230614-141034' for a sequence acquired on 2023/06/14 14:10:34).
Each sequence folder have to be in a global folder called 'raw'.
Overall your folder architecture needs to be like this :

↳ raw
  ↳ 20230614-141034 (sequence folder)
    ↳ images.zip
  ↳ 20230614-151034 (sequence folder)
    ↳ images.zip
  ↳ ...

You also need to know the orientation of the UVP. Because you have to know where are the top and the bottom of the frame (in order to know if particles goes up or down !)

"""

# Libraries
import sys
import os

# Path to modules directory
# Where you store the 'trackpartuvp' directory
path_modules = "these_alex/script/trackpartuvp-main/trackpartuvp"
sys.path.append(path_modules)

print(sys.path)


# Current working directory
# Where you store your UVP6 projects
path = "/home/aaccardo/"

# Import modules
from trackpartuvp.particles import get_particles
from trackpartuvp.tracks import get_tracks
from trackpartuvp.trap_depth import calc_trap_speed

'''
C1_path = 'these_alex/Analysis/APERO/TZEX/cycle_1'
C2_path = 'these_alex/Analysis/APERO/TZEX/cycle_2'
C3_path = 'these_alex/Analysis/APERO/TZEX/cycle_3'
C4_path = 'these_alex/Analysis/APERO/TZEX/cycle_4'
C5_path = 'these_alex/Analysis/APERO/TZEX/cycle_5'
Vlfr_path = 'these_alex/Analysis/APERO/TZEX/vlfr_test'
UVP_calibr = 'these_alex/Analysis/UVP_calibration'
sequence_selected = 'these_alex/Analysis/APERO/TZEX/all_cycle_analysis/selected_sequences'
'''

DeepTrap_UVP = 'plankton/uvp6_missions/uvp6_sn000110lp/uvp6_sn000110lp_202206_atlanteco_pps5_mooring_track'

# Path to UVP6 project
source_directory = path + f'{DeepTrap_UVP}' #mettre chemin jusqu'au projet uvp (avant dossier raw)
os.chdir(source_directory)

# Path to the destination folder
dest_directory = path + '/these_alex/Analysis/APERO/DeepTrap_UVP'


float_speed_path = ''

# Modify UVP6 parameters (A, B and threshold)
# They can be found in 'data' files of UVP projects

#get_particles(directory = directory, A = 2300, B = 1.136, threshold = 21, float_speed_path = '') # for visutrap tristan-manon
get_particles(source_directory = source_directory, dest_directory = dest_directory, A = 2300, B = 1.136, threshold = 22, float_speed_path = float_speed_path) # for TZEX

# Add deployment name and depth
get_tracks(source_directory = source_directory, dest_directory = dest_directory, deployment = 'DeepUVP', depth = '', float_speed_path = float_speed_path) # add the path to the float vertical speed


#calc_trap_speed(directory, 'vertical', '') #calcule la vitesse d'oscillation du piège (à adapter pour le fotteur)


