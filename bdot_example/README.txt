README

R. Datta 
6/25/2024

This code processes the raw inductive probe signals from MARZ.

Use this code to:

(1) Load and plot raw voltage signals. (Section 1)
(2) Combine opposite polarity signals to calculate electrostatic and inductive components, and magnetic field. (Section 2)
(3) Calculate velocity from time-of-flight. (Section 3).

The main scripts are:

bdot_analysis_shot4.m - Analysis of MARZ 4
bdot_analysis_shot3.m - Analysis of MARZ 3
bdot_analysis_shot2.m - Analysis of MARZ 2
bdot_analysis_shot1.m - Analysis of MARZ 1

All other functions are helper functions.

Please see bdot_analysis_shot3.m or bdot_analysis_shot4.m for detailed instructions and methodology.

Note: The scripts and helper functions differ shot-to-shot. This is because of different naming schemes, incorrect labeling, or different data storage between the shots.

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.