# MARZ_Analysis_Tools

R. Datta, Aug 15, 2024

This repo contains code for post-processing and analyzing data collected on the MARZ platform. 

Please download and unzip each sub-folder to access the code and instructions specific to a given diagnostic.

Experimental details are in:

[1] R. Datta et al. (2024) PRL. https://doi.org/10.1103/physrevlett.132.155102 and
[2] R. Datta, et al. (2024). Phys. Plasmas. https://doi.org/10.1063/5.0201683

The diagnostics include:

# (1) Inductive / B-dot Probes

Data from each probe is stored as .CSV files.
We perfom common mode rejection to find the inductive signal, and integrate the signal(s) using known calibration factors to determine the magnetic field.

# (2) Visible Spectroscopy (SVS)

Data is stored in HDF files. For pre-processing, we need shot data, as well as calibration data, which includes pre-shot laser images, LDLS fast and slow images, and Tungsten lamp images. Details of this pre-processing are provided in M. Schaeuble et al. Phys. Plasmas, 28(6):062902, 2021. https://doi.org/10.1063/5.0047931

The pre-processed spectra are then compared to synthetic spectra to determine density and temperature. This is done using a custom library that solves the 1-D radiation transport equation.

# (3) X Ray Spectroscopy (XRS3)

Data is stored in TIFF files. Measured spectra are calibrated using the energies of known lines.

The spectra are then compared to synthetic spectra to determine density and temperature. This is done using a custom library that solves the 1-D radiation transport equation.

# (4) Ultra Fast Imaging (UXI) Pinhole Cameras



# (5) Self-Emission Gated Optical Imager (SEGOI)

# (6) X Ray Diodes (TADPoles and LOS 170 Silicon Diodes)

# (7) Current Diagnostics (IDTLs, PDV, and Machine B-dots)




