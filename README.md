# sFCS_BTS

Simulation and analysis code for scanning fluorescence correlation spectroscopy (sFCS) data. 
This repository will accompany a manuscript describing the use of sFCS and brightness transit statistics (BTS) histograms to quantify biomolecular organisation. Collaboration with the [Fritzsche Lab](https://www.bpi-oxford.com/). 

## The folder "Example_data" contains:

- "ExampleFittingParameter_EGFP-His_lowConc_780" and "ExampleFittingParameter-EGFP-His_HighConc_980.xlsx"

These are examples for fitting results from a series of sFCS measurements of an SLB labelled with EGFP-His (at high or low EGFP-His concentration). The Excel sheet contains the sFCS fitting parameters (such as transit time, brightness or number of molecules), which are used for down-stream statistical analysis. The fitting and exporting of the data was performed using the [FoCuS_scan software package](https://github.com/dwaithe/FCS_scanning_correlator) written by [Dr. Dominic Waithe](https://github.com/dwaithe). The low concentration data were taken on a Zeiss LSM 780 versus the high concentration data were taken on a Zeiss LSM 980. 

- "sFCS_780_2069Hz_3p94us.lsm" and "sFCS_980_3415_0p98.czi"

sFCS raw data from a Zeiss 780 (AF488 His on a bilayer) and from a Zeiss 980 (EGFP-His on a bilayer). For analysis, the files should be correlated using a either scanning frequency of 2069 Hz and a pixel dwell time of 0.00394 ms or a scanning frequencies of 3415 Hz and a pixel dwell time of 0.00098 ms. The measurement should allow to load the acquisition settings on the respective microscopes (pixel size, acquisition time, laser power, MBS etc.) using the "Reuse" button. 

## The folder "Scripts" contains the following analysis or simulations scripts written in Python:

- "sFCS_freeDiff_simulation_parallel.py"

Script to simulate sFCS data for free diffusion. Molecules will be randomly distributed in a simulation box of variable size. Brownian motion will be simulated (diffusion coefficient, track length etc can be adjusted). A sFCS measurement will be taken on those diffusing particles (line length, scanning frequency, PSF size etc can be adjusted). To allow simulation of a large number of molecules in the simulation box, the computational load is divided using parallel computing ([iPyParallel](https://pypi.org/project/ipyparallel/)). Before running the simulation, the enginges need to be initialised. Run "ipcluster start -n XX" from the anaconda prompt. Replace XX with the number of CPU cores available or to be used. Then the IDE and the code can be run (it should print out the number of initialised engines/cores).

- "sFCS_freeDiff_simulation_singleCPU.py"
Script to simulate sFCS data for free diffusion as above but not using parallel computing. 


- "sFCS_freeDiff_simulation_oligomerisation_parallel.py"

Script to simulate sFCS for free diffusion in the presence of oligomers/clusters. Similar to the script on simple diffusion here fraction, size and diffusion coefficient of oligomers can be defined and adjusted. 

- "raw2fittingpara.py"

End-to-end script for sFCS data analysis. Intensity fluctuation data can be loaded as .tif file (from simulation or converted experimental data) or as lsm5 files (experimental data). The script will generate intensity traces for every file in a specified folder, correlate the intensity carpets and fit the data to a simple 2D diffusion model. The output is a .csv file with the fitting parameters for further downstream processing, ie. statistical analysis. The analysis represents an automated workflow but is not meant to replace the analysis in the [FoCuS_scan software package](https://github.com/dwaithe/FCS_scanning_correlator) which also allows for cropping or photobleaching correction. Rather this script provides additional functionality such as correction for detector dark counts or the calculation of quality control parameters such as signal-to-noise-ratio (SNR) or normalised root mean squared residuals (nRMSD) values for data evaluation. 

- "BTS_analysis.py"

Script to load and process sFCS fitting parameters (as excel sheet or .csv file). Transit time data should be given as a column named 'txy1' and brightness as a column 'cpm (kHz)'. Multiple conditions can be loaded for comparison and plot parameters can be user-adjusted. 

- "Analyse_sFCS_data_MLE.py"

Script to statistically analyse transit time data and reveal the underlying Brownian or non-Brownian diffusion dynamics. The script performs a maximum likelihood estimation on the data using a normal distribution, a lognormal, and a double lognormal distribution. This allows to calculate the Bayesian Information Criterion for each model from which it calculates a relative likelihood (RL) value for each model. This RL value indicates which model represents the data best. Free diffusion is best described by a lognormal distribution. Hindered diffusion (ie. trapping with transient nano-scale interactions in the diffusion path) results in a two-component / double lognormal distribution. The concept of the analysis is described [here](https://pubs.acs.org/doi/10.1021/acsnano.8b04080) and in details in the [SI](https://pubs.acs.org/doi/suppl/10.1021/acsnano.8b04080/suppl_file/nn8b04080_si_001.pdf).

- "Compare2BTS_datasets.py"

Script to statistically compare two BTS datasets, e.g., control versus treated. The files should be loaded as above as excel sheet or .csv file. Transit time data should be given as a column named 'txy1' and brightness as a column 'cpm (kHz)'. Plot parameters can be user-adjusted. The script calculates the Battacharyya distance as metric for the difference between the two conditions. The data is then concatenated, permutated, and subsampled numerous times (bootstrapping approach) to create a null distribution of Battacharyya distances. Comparing the Battacharyya distance between the conditions and the null distribution allows to calculate a p-value to indicate if the two conditions were likely from the same underlying distribution. A p-value < 0.01 indicates a significant difference between the two input BTS histograms. 
The analysis was inspired by [this approach](https://thenode.biologists.com/user-friendly-p-values/research/). 


All scripts contain further information on functionality, input and output parameters. 
Some parts of the simulations and analysis were further developed or re-implemented in Python based on our previously published works (e.g. 10.1016/j.ymeth.2017.09.010, 10.1021/acsnano.8b04080, 10.1021/acs.nanolett.8b01190).   

Questions and comments welcome. 
@faldalf



## Technical details and requirements

All scripts were written and executed under Python 3.7.7 in Spyder (4.1.5, conda version 4.10.3) on a Windows 10 (64-bit) laptop with an Intel(R) Core(TM) i7-6600U CPU @ 2.60GHz and 16 GB RAM. 

**Python dependencies**
The following Python libraries are required to run the above scripts. The indicated versions were the ones used during initial realease of the repository. 

- numpy 1.18.5
- matplotlib 3.2.1 
- scipy 1.5.2
- pandas 1.1.3
- tkinter (or rather tk) 8.6.10
- tifffile 2020.6.3
- multipletau 0.3.3
- lmfit 1.0.1
- ipyparallel 6.2.4
