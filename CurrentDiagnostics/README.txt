README

R. Datta 
7/30/2024


Use this code to:

(1) Load and plot the IDTL, PDV, and Current B-dot signals.

No processing is required, as the signals are pre-processed by Sandia staff.

Please see these scripts for detailed instructions and methodology.

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.