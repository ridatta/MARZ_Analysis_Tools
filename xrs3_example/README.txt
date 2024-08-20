README

R. Datta 
7/30/2024

This code visualizes and plots the raw XRS3 X-Ray spectroscopy Data, and analyses the experimental spectra through comparison with SCRAM + Radiation Transport Synthetic Spectra.

Use XRS3_marz<sht num>.mlx code to:

(1) Load and plot raw XRS3 data.
(2) Perform energy calibration.
(3) Characterize variations in line position along axial position.
(4) Export experimental spectra for analysis.

Please see scripts [e.g. XRS3_marz3.mlx] for detailed instructions and methodology.

ANALYSIS (in .\SCRAM_V2_Analysis\):

(1) You will need access to SCRAM tables. Use loadSCRAM.ipynb to load SCRAM tables and convert units. 
(2) Use KNR.ipynb to interpolate the SCRAM tables, and to create functions that output opacity and emissivity given any temp. and density.
(3) Use NearestNeighborModelValidation.ipynb to validate the interpolation model.
(4) Use marz<shot_num>_z=<position>_analysis.ipynb to compare synthetic spectra to the experimental spectra. Calculate bounds on density and temperature based on a radiation transport calculations, based on a single hotspot model. Generate corner plots showing bounds on hotspot density, size, and temperature. 

Library of functions required for this analysis in contained in scram.py.

NOTE: You should NOT have to repeat steps (1-3), unless the SCRAM tables change. Proceed with Step 4 directly in this case. 

Please see these scripts for more detail. 

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.