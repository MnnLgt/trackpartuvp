# trackpartuvp: Tracking particles from UVP6 images

This code is desgined to measure the *in situ* settling velocity of marine
particles from UVP6 data, whether the UVP6 is mounted on a sediment trap or a
float.

As the associated paper is in the review process, please cite my PhD dissertation ([accessible here](https://theses.hal.science/tel-04552165/)) if you use part or all of this code.

## Overview

The main goal of this code is to:

1. Detect particles on a sequence of UVP6 full images with Python *skimage* module
2. Link particles to track them and calculate their speed
3. Visualize detected tracks

## Running the code

### Input data

When running the function ```get_particles```, you have to specify:

- ```source_directory```: path to UVP6 project containing raw images. Make sure it includes a ```raw``` folder.
- ```dest_directory```: path to result directory.
- ```A```, ```B``` and ```threshold```: these parameters depend on the calibration of your UVP6. Please check in the config file associated to your project.

Optional arguments:

- ```platform_speed_path```: path to platform speed data.

When running the function ```get_tracks```, you have to specify:

- ```source_directory```: path to UVP6 project containing raw images. Make sure it includes a ```raw``` folder.
- ```dest_directory```: path to result directory.
- The name of your ```deployment```
- The deployment ```depth```

Optional arguments:

- ```platform_speed_path```: path to platform speed data.
- ```angle```: angle used for tracking particles.
- ```inclino_params```: Path to dataframe containing path to inclino data, axis and sign to use.
- ```save_zip```: whether to save a zip file of track results.

### Output data

A ```results/``` folder is created. It is organized as follows:

- ```particles_data/``` contains ```.csv``` files with particles data, one per sequence
- ```detected_tracks/``` contais a subdirectory with the tracks for each sequence
  - ```<sequence_name>/```

### Repository organization

- ```funcparticles/```
- ```functracks/```
- ```utils/```
- ```particles.py```
- ```tracks.py```
- ```platform_depth.py```
- ```plot_ecotaxa_from_df.py```
- ```zip_folders.py```

## Running the code




