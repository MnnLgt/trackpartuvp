# trackpartuvp: Tracking particles from UVP6 images

*In situ* settling velocity of marine particles from UVP6 data

## Overview

The main goal of this code is to:

1. Detect particles on a sequence of UVP6 full images with Python *skimage* module
2. Link the particles to track them and calculate their speed
3. Visualize detected tracks

## Structure

### Input data

When running the code, you have to specify:

- ```directory```
- ```A```, ```B``` and ```threshold```
- ```deployment```
- ```depth```

You can set up optional arguments:

- ```path_inclino```, ```axis``` and ```sign```
- ```save_zip```
- ```angle```

### Output data

- ```results/``` folder is created in the global directory. It is organized as follows:

- ```binned_size_spectra.csv``` file containing the size spectra of particles, binned per sequence
- ```particles_data/``` contains ```.csv``` files with particles data, one per sequence
- ```detected_tracks/``` contais a subdirectory with the tracks for each sequence
  - ```<sequence_name>/```

### Repository organization

- ```func_particles/```
- ```func_tracks/```
- ```utils/``` 
- ```particles.py```
- ```tracks.py```
- ```trap_depth.py```

## Running the code

## Examples


