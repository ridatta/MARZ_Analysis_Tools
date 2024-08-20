README

R. Datta 
6/25/2024

This code processes the UXI data from MARZ.

Use this code to:

(1) Load and plot raw UXI images.
(2) Add a scale to the UXI images.
(3) Estimate velocity of hotspots from UXI images.


The main scripts are:

UXI_load.mlx
VelocityUXI1.mlx
VelocityUXI2.mlx

All other functions are helper functions.

Note: The scripts and helper functions differ shot-to-shot. This is because of different naming schemes, incorrect labeling, or different data storage between the shots.

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.