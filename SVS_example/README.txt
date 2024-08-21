README

R. Datta 
6/25/2024

This code processes the raw SVS streak data from MARZ to export corrected intensity signals for analysis.

Use this code to:

(1) Load and display raw streak data from SVS. 
(2) Perform corrections to the raw image.
(3) Wavelength and time calibration. 
(4) Find instrument broadening.
(5) Correct for instrument response.
(6) Transit time correction.
(7) Take lineouts of corrected and calibrated spectra, and export them for radiation transport / CR analysis.

The analysis code is stored in the library SVSAnalysis.m.

A good example analysis with detailed instructions and methodology is SVS5_marz1.mlx.

Other codes in this folder are helper functions.

Analysis scripts are named as follows SVS<system no.>_marz_<shot no.>.mlx e.g. SVS5_marz1.mlx is SVS5 from MARZ 1.

Note: You may require Sandia's SMASH toolbox for some analysis https://github.com/SMASHtoolbox/release.git 



Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.
