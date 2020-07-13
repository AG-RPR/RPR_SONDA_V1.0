# REPERIO SONDA V1.0 - Technical Documentation
Reperio B.V. (www.reperio-medtech.com) provides time-efficient, patient-friendly diagnostic solutions for neurological and ophthalmological impairments. The proprietary technology is a combination of high-speed oculography and Deep Learning that allows achieving extremely fast quantifications of visual functions. It has applications in, amongst others, Glaucoma, different forms of Retinopathies, Parkinson's Disease, Multiple Sclerosis and Acquired Brain Injury.
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
The version in this repository can be interfaced exclusively with SR-Research eyetrackers.

### Operative System
MacOS (latest stable version tested on MacOS 10.12 Sierra) is the preferred OS for SR-Research eye-trackers. 
Windows (7 or later) 64bit is the preferred OS for Tobii eye-trackers. (version not present in this repo)

### Display
It is recommended to use displays that can run at a stable refresh rate of 120 Hz or above, with a minimum resolution of 1920x1080.

## Installation
Run the sonda_install.m file. This script will automatically: 

   - install Psychtoolbox, necessary for displaying the visual stimuli and ensure proper time-sync between the display and the eye-tracker.
   - create and store the stimulus trajectory paths, for smooth and saccadic pursuit modalities
   - determine which operative system is running
   - create a init_config.log file for future reference
   
## Data acquisition
- acquire\_gaze\_data.m 

- rnd\_list.m: create a .txt file containing the list of conditions of the gaze acquisition session

- rnd\_fixpath.m: pre-compute the smooth trajectories of the stimulus at different velocity levels, linearly scaled with each other. It is a rather inefficent way to generate random walks that respect the boundaries of the screen, but it's 100% guaranteed that the resulting trajectories do NOT contain any periodic auto-correlations (which are detrimental for the cross-correlogram analysis).

 