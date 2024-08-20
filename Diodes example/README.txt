README

R. Datta 
7/30/2024

This code processes the X-Ray Diode (TadPole and LOS 170 Diode Signals).

Use DiodeSignals.mlx to:

(1) Load and plot voltage data from X ray diodes.
(2) Convert the diode voltage to measured X ray power.

Use DiodeAnalysis.ipynb to:
(1) Using Spk, generate synthetic diode signals for unfiltered and filtered diodes.
(2) Infer bounds on density and temperature from the diode signals.
(3) The synthetic signals plot emissivity and opacity, as well as intensity, accounting for radiation transport
(4) Radiation transport can be solved in 1-D, assuming different geometries - sphere or infinite planar slab
(5) Estimate volumetric radiation power loss rate [Wm-3] for given plasma parameters.

Please see these scripts for detailed instructions and methodology. Additional code in this folder are helper functions. 

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.
