# REPERIO SONDA V1.0 - Technical Documentation
Reperio B.V. 
## Content:
- Introduction
- System requirements
- Installation 
- Data acquisition
- Data analysis
- Misc

## Introduction
This document contains the *technical* documentation of the REPERIO SONDA V1.0 MATLAB toolbox. This toolbox allows to collect gaze data with video eye-tracking using the SONDA continuous tracking task, and then provides several oculomotor metrics (Grillini et. al 2018, Grillini et al. 2020). This document is intended for researchers and developers who wants to use (and potentially expand) this toolbox for research and non-profit purposes only. Regular users can refer to the User Manual. 

## System requirements 
### Eye-tracker models
SONDA can be interfaced with desktop-mounted eye-trackers from two manufacturer: SR-Research (Ottawa, Ontario, Canada) and Tobii (Stockholm, Sweden). The validation has been performed on EyeLink 1000 (SR-Research), EyeLink Portable Duo (SR-Research) and Tobii T60XL (Tobii).

### Operative System
MacOS (latest stable version tested on MacOS 10.12 Sierra) is the preferred OS for SR-Research eye-trackers. 
Windows (7 or later) 64bit is the preferred OS for Tobii eye-trackers. 

### Display
It is recommended to use displays that can run at a stable refresh rate of 120 Hz or above, with a minimum resolution of 1920x1080.

## Installation
Run the sonda_install.m file. This script will automatically: 

   - install Psychtoolbox, necessary for displaying the visual stimuli and ensure proper time-sync between the display and the eye-tracker.
   - create and store the stimulus trajectory paths, for smooth and saccadic pursuit modalities
   - determine which operative system is running
   - create a init_config.log file for future reference
   
## Data acquisition
acquire\_gaze\_data.m 


