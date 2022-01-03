# sFCS_BTS

Simulation and analysis code for scanning fluorescence correlation spectroscopy (sFCS) data. 
This repository will accompany a manuscript describing the use of sFCS and brightness transit statistics (BTS) histograms to quantify biomolecular organisation.

## The folder "Example_data" contains:

- "ExampleFittingParameter_SLB.xlsx"
This is an exemplary file from a series of sFCS measurements of an SLB labelled with EGFP-His. The Excel sheet contains the sFCS fitting parameters (such as transit time, brightness or number of molecules), which are used for down-stream statistical analysis. The fitting and exporting of the data was performed using the [FoCuS_scan software package](https://github.com/dwaithe/FCS_scanning_correlator) written by [Dr. Dominic Waithe](https://github.com/dwaithe). 

- "sFCS_2069Hz_3p94us.lsm"
Experimental sFCS acquisition on an SLB labelled with AF488-His. For analysis, the file should be correlated using a scanning frequency of 2069 Hz and a pixel dwell time of 0.00394 ms. The measurement was taken on a Zeiss780 and should allow to load the acquisition settings (pixel size, acquisition time, laser power, MBS etc.) using the "Reuse" button. 

## The folder "Scripts" contains the following analysis or simulations scripts written in Python:

- "sFCS_freeDiff_simulation_parallel.py"
Script to simulate sFCS data for free diffusion. Molecules will be randomly distributed in a simulation box of variable size. Brownian motion will be simulated (diffusion coefficient, track length etc can be adjusted). A sFCS measurement will be taken on those diffusing particles (line length, scanning frequency, PSF size etc can be adjusted). To allow simulation of a large number of molecules in the simulation box, the computational load is divided using parallel computing ([iPyParallel](https://pypi.org/project/ipyparallel/)). 

- "sFCS_freeDiff_simulation_oligomerisation_parallel.py"
Script to simulate sFCS for free diffusion in the presence of oligomers/clusters. Similar to the script on simple diffusion here fraction, size and diffusion coefficient of oligomers can be defined and adjusted. 

- "raw2fittingpara.py"
End-to-end script for sFCS data analysis. Intensity fluctuation data can be loaded as .tif file (from simulation or converted experimental data) or as lsm5 files (experimental data). The script will generate intensity traces for every file in a specified folder, correlate the intensity carpets and fit the data to a simple 2D diffusion model. The output is a .csv file with the fitting parameters for further downstream processing, ie. statistical analysis. The analysis represents an automated workflow but is not meant to replace the analysis in the [FoCuS_scan software package](https://github.com/dwaithe/FCS_scanning_correlator) which also allows for cropping or photobleaching correction. Rather this script provides additional functionality such as correction for detector dark counts or the calculation of quality control parameters such as signal-to-noise-ratio (SNR) or normalised root mean squared residuals (nRMSD) values for data evaluation. 

- "BTS_analysis.py"
Script to load and process sFCS fitting parameters (as excel sheet or .csv file). Transit time data should be given as a column named 'txy1' and brightness as a column 'cpm (kHz)'. Multiple conditions can be loaded for comparison and plot parameters can be user-adjusted. 

- "Analyse_sFCS_data_MLE.py"
Script to statistically analyse transit time data and reveal the underlying Brownian or non-Brownian diffusion dynamics. The script performs a maximum likelihood estimation on the data using a normal distribution, a lognormal, and a double lognormal distribution. This allows to calculate the Bayesian Information Criterion for each model from which it calculates a relative likelihood (RL) value for each model. This RL value indicates which model represents the data best. Free diffusion is best described by a lognormal distribution. Hindered diffusion (ie. trapping with transient nano-scale interactions in the diffusion path) results in a two-component / double lognormal distribution. The concept of the analysis is described [here](https://pubs.acs.org/doi/10.1021/acsnano.8b04080) and in details in the [SI](https://pubs.acs.org/doi/suppl/10.1021/acsnano.8b04080/suppl_file/nn8b04080_si_001.pdf).

All scripts contain further information on functionality, input and output parameters. 
Some parts of the simulations and analysis were further developed or re-implemented in Python based on our previously published works (e.g. 10.1016/j.ymeth.2017.09.010, 10.1021/acsnano.8b04080, 10.1021/acs.nanolett.8b01190).   

Questions and comments welcome. 
@faldalf
