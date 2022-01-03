# BTS_analysis.py
# Script to generate/analyse BTS histograms from sFCS data.
# This is a cleaned-up version to accompany the BTS histogram manuscript.
#
# sFCS fitting parameters (as for example generated by FoCuS_scan) are required to be .xlsx or .csv files as input. 
# The columns should be named 'txy1' and 'cpm (kHz)'
#
# User needs to define the sample names/conditions in the list "samples"
# User can adapt plotting parameters as appropriate. 
# Plots need to be saved to .png or .pdf manually
# Falk Schneider @faldalf 02/01/2022
#
# Questions and feedback are welcome 
#

#%% Imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from tkinter import filedialog
from tkinter import *


#%% Functions

def load_sFCS_data (samples): 
# Small helper function to read in sFCS fittingparameters from FoCuS_scan or simulations    

    root = Tk()
    root.filename =  filedialog.askopenfilename(title = ["Load file to analyse: "+ samples],filetypes = (("FCS Fittingparameters",".csv .xlsx"),("all files","*.*")))
    filename = root.filename
    filename = filename.replace("/","\\\\")
    root.withdraw()

    if 'xlsx' in filename:
        data = pd.read_excel (filename)
        data = data.drop(data.index[-1]) # removes the last enrty ("end" or "NaN" from FoCuS_scan)
        print (os.path.basename (filename), ' has been loaded')
        print ('Warning: Last row has been removed')
    elif 'csv' in filename:
        data = pd.read_csv (filename)
    else:
        print ('Error: Can\'t read file format. Please provide .xlsx or .csv')
        
    return data
    
#%% Load data

samples = ['Condition 1', 'Condition 2'] # Sample names (conditions) in a list. Length of this defines how many plots are generated! 

data = {} # Predefine data dict

simulation = False # Indicate if data are simulated (dealing with differences in naming of variables). For analysing experimental data keep False. 

for i in range (0, len (samples)):
    data [i] =  load_sFCS_data (samples[i])
    if simulation: 
        data [i] = data[i].drop(data[i].index[-1]) 
        data[i]['txy1'] = data[i]['txy'].astype(float)
        data[i]['cpm (kHz)'] = data[i]['cpm(kHz)']
        


#%% Plotting

# Define layout of BTS histogram
x_min = 0
cpm_max = 30
txy_max = np.log (100)
binning = 40 # Most of the paper plots are in 40
vmax = 0.50 # colorbar max. 

fig1, axes = plt.subplots(nrows=1, ncols=len(samples),  figsize=(12, 4)) 

for i, ax in enumerate(axes.flatten()):
    ax.set_xlabel ('Ln (Transit time (ms)/(ms))', fontsize = 14)
    ax.set_ylabel ('Brightness (kHz)', fontsize = 14)
    ax.tick_params (axis='x', labelsize = 14)
    ax.tick_params (axis='y', labelsize = 14)
    

for data_set in data.keys ():
    histo = axes [data_set].hist2d (np.log(data [data_set] ['txy1']), data [data_set] ['cpm (kHz)'], bins = binning, range = [[x_min, txy_max],[0, cpm_max]], density = True, vmax = vmax) 
    
    axes [data_set].set_title (samples[data_set])
    

plt.tight_layout()

# Generating the LUT colorbar in a different figure. 
plt.figure ()
cbar1 = plt.colorbar(histo[3], format = "%.2f")
cbar1.set_label (label = 'Probability (a.u.)', size = 14)
cbar1.ax.tick_params (labelsize=14)


#%% Generate transit time histogram. 

fig1, axes = plt.subplots(nrows=1, ncols=1,  figsize=(6, 6))

for data_set in data.keys ():
    histo = axes.hist (np.log(data [data_set] ['txy1']), bins = binning, range = [x_min, txy_max], density = True, edgecolor = 'k', label = [samples[data_set]], alpha = 0.7) 
    
axes.set_title ('Transit time histogram')
axes.set_xlabel ('Ln (Transit time (ms)/(ms))', fontsize = 14)
axes.set_ylabel ('Probability (a.u.)', fontsize = 14)
axes.tick_params (axis='x', labelsize = 14)
axes.tick_params (axis='y', labelsize = 14)
axes.legend ()
