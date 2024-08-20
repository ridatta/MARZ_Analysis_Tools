
README

R. Datta 
7/30/2024


Use this code to perform SVS analysis. Note that the experimental SVS data must be pre-processed and exported as CSV files for comparison before analysis.

(1) Generate PrismSPECT Emissivity and Opacity Tables.

- Use PrismSPECT (or your code of choice) to run CR simulations to generate emissivity and opacity data. 
- Previously run simulations are stored in Prism. You can determine the density and temperature ranges of these simulations by opening the transpwr_<density idx>_<temp. idx>.dat files, and the possible atomic transitions by opening the corresponding .pop files using a text editor.

(2) Use data-processing.ipynb to:
- Load and visualize the PrismSPECT emissivity and opacity data. 
- Create an interpolation model to generate  spectral emissivity and opacity value for arbitrary ni and Te.
- In the current analysis, Prism simulation data for density and temp. in the range 1e16 < ni < 1e19 per cc and 0.5 < Te < 15 eV are used. 
- Interpolator models are stored in ./models/ and called directly in step (3). Please download these models from https://zenodo.org/records/13351147 and add them to the working directory.
- These models are stored in radTran.py and called using the getEmi and getOpa functions.
These functions can be replaced by your own functions that return emissivity and opacity given density and temperature.

(3) Use Simulator.ipynb to:
- Generate synthetic intensity spectra by solving radiation transport using the PrismSPECT emissivity and opacity data.
- Library of functions required for this analysis in contained in radTran.py
- Radiation transport can be solved in 1D geometry using this code for arbitrary density and temp. profiles; for our analysis, we solve for Gaussian or uniform density and constant temp. 

(4) Use ExpComparison_<shot num>_svs<svs num>.ipynb to determine density and temp. though curve fitting.
- Load the experimental spectrum, and perform continuum subtraction.
- Isolate the lines of interest, and get density by fitting to a well-separated line.
- Get temp. by fitting to interstage lines. 
- Compare the overall fit to the experimental spectrum.

Preliminary analysis for MARZ4 data has been performed by R. Datta in ExpComparison_z3978_svs5.ipynb.

NOTE: You should NOT have to repeat steps (1-2), unless you need new emissivity and opacity data. Proceed with Step 3 onwards directly.

Please see these scripts for more detail. 

Important considerations:

(1) Ensure input directory path containing raw data and/or metadata points to the correct location on your machine.
    
The checkDir() function is NOT required.

(2) Ensure save directory path where processed data is stored is accurate on your machine.

The checkDir() function is NOT required.

(3) Additional libraries may be required for some analysis. 

These can be downloaded from https://github.com/ridatta/PlasmaFormulary.git. 

Use addpath(/path/to/downloaded/code/) to add the downloaded libraries to the MATLAB workspace.
