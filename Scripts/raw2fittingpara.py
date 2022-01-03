# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 15:29:49 2021

@author: Falk Schneider
"""
#
# raw2fittingpara.py
#
# Script for automatically analysing sFCS data from raw fluctuation files (tiff = simulations, lsm5 = measurements).
# This is a cleaned-up version to accompany the BTS histogram manuscript.
# This script is not meant to replace analysis with FoCuS_scan (https://github.com/dwaithe/FCS_scanning_correlator) but is based on it and can be used in addition to calculate nRMSD, SNR or apply detector dark count corrrections for N&B type analysis.
# The script reads all files to analyse in a folder, correlates the data and fits them. Fittingparameters are exported as .csv file. 
# 
# Falk Schneider @faldalf 02/01/2022
#
# Questions and feedback are welcome 


#%% Imports


from tkinter import filedialog
from tkinter import *
import os
import tifffile
import numpy as np
from multipletau import autocorrelate
from lmfit import minimize, Parameters,report_fit,report_errors, fit_report
import matplotlib.pyplot as plt
import pandas as pd
import copy

import fnmatch

def residual(param, x, data):
    A = equation_(param, x)
    residuals = np.array(data)-np.array(A)
    return np.array(residuals).astype(np.float64)

# This is our standard equation we fit. 
def equation_(param, tc):
    """This is equation for fitting"""
    

    #A1 is relative component of fluorescent species
    #tc is tau.
    #txy1 is xy difusion   for fluorescent species one.

    offset =param['offset'].value; 
    GN0 =param['GN0'].value; 
    A1 = param['A1'].value; txy1 = param['txy1'].value; alpha1 = param['alpha1'].value;
    #For one diffusing species
    GDiff = (A1*(((1+((tc/txy1)**alpha1))**-1)))
    return offset + (GN0*GDiff)

#%% The script

# Important global input parameters
scanning_frequency = 2000 #2069 # Hz
dwell_time = 4 #3.94 #4 # pixel dwell time in us; please remember it's 4 us for the simulations!!!
lagtime_max = 1000 # maximum lag time for fitting in ms
num_of_channels = 1 # Number of channels in the files

std_limit = 10 # number of points -1 used to calculate the SNR e.g. first 10: First nine correlation points

line_dur = ((1.0/scanning_frequency)*1000.0) # in ms

#Global parameters for fitting
param = Parameters()
param.add('offset', value=0.01, min=-0.5, max=2.5, vary=True);
param.add('GN0', value=1.000, min=-0.0001, max=10000.0, vary=True);
param.add('A1', value=1.000, min=0.0001, max=1.0000, vary=False);
param.add('txy1',value=20, min=0.001, max=5000.0, vary=True);
param.add('alpha1',value=1.0, min=0.600, max=2.0, vary=False);
options = {'Dimen':1,'Diff_eq':1,'Triplet_eq':1,'Diff_species':1}

# Correction factors for moment estimation of the brightness
# These values need to be measured for every detector and ideally be checked for every experiment
var_0 = 0# -0.016# Readout noise, default is 0 (Needs to be measured experimentally for every detector)
S = 1.04# Correction factor for variance(intensity) e.g. 1.06#1.04 #default is 1 (Needs to be measured experimentally for every detector)
detector_offset =  0.072  # Dark counts of the detector in kHzn#default is 0 (Needs to be measured experimentally for every detector)


# Easy selection of file path with a small GUI window
root = Tk()
root.withdraw()
path = filedialog.askdirectory()
path = path.replace ("/","\\\\")
path = path + '\\\\'

print ('Starting analysis')

# Here we loop over the channels in the file(s)
# For now I only accept either one or two channel files. 

ch = 0
for ch in range (0, num_of_channels):
    print ('\n Channel ' + str (ch))

    # Create placeholder variables to store all the parameters
    corr_carpet = {}
    corr_carpet ['curves']={}
    
    corr_carpet_temp = {}
    corr_carpet_temp ['curves']={}
        
    results = {}
        
    txy_array = np.array([])
    
    G0_array = np.array([])
    int_array = np.array ([])
    cpm_array = np.zeros([len (G0_array)])
    nRMSD_array = np.array ([])
    SNR_array = np.array ([])
    
    cpm_mom_array = np.array ([]) # assuming ideal photon counting detector
    cpm_mom_corr_array = np.array ([]) # corrected for the detector
    
    curve_counter = 0
    fit_counter = 0
      
    
    
    for filename in os.listdir(path):
        if filename.endswith(".lsm") or filename.endswith(".tif"):
          print (filename, ' is being correlated')
          files = (path + filename)
          int_data_in = tifffile.imread(files,key=0) # For reading lsm files. This is very important!
          if filename.endswith(".lsm"):
              if num_of_channels > 1:
                  int_data = int_data_in [ch,:,:]
              else:
                  int_data = int_data_in
          elif filename.endswith(".tif"): 
              int_data = int_data_in
          # Let's save the int_data
          if ch == 0:
              time_ms = np.cumsum (np.ones(len (int_data)))*line_dur
              int_data_avg = np.sum (int_data,1)
              fig_int, axes = plt.subplots (nrows = 1, ncols = 1, figsize = (8,4))
              axes.plot (time_ms, int_data_avg, label = 'Channel: 0' )
              axes.set_ylabel ('Intensity (photons)')
              axes.set_xlabel ('Time (ms)')
              axes.set_title ('Intensity time trace')
              if num_of_channels == 2 and ch == 0:
                  int_data_avg = np.sum (int_data_in [1,:,:],1)
                  axes.plot (time_ms, int_data_avg, label = 'Channel: 1')
              axes.legend ()
              plt.tight_layout ()
              
              #Save Intensity trace for file
              fig_int_path = files[:-4]+'CH_'+str(ch)+'_Int_trace.png'   
              fig_int.savefig(fig_int_path, bbox_inches='tight')
          
          
          #Correlate file
          for curve in range (0,int_data.shape[1]):
              corr_carpet_temp ['curves'][curve] = autocorrelate(int_data[:,curve] ,normalize=True,deltat=line_dur, copy=True, dtype=None)
              corr_carpet ['curves'] [curve_counter]= copy.deepcopy (corr_carpet_temp ['curves'][curve])
              avg_int = (((np.average (int_data [:,curve])) / dwell_time) * 10**3) # go from count per bin to counts in kHz
              int_array = np.append (int_array,avg_int)
              
              #Calculate SNR by dividing the whole trace into three chunks and extracting mean and standard deviation.
              chunk1 = autocorrelate(int_data[0:np.int(np.floor(len(int_data)/3)),curve] ,normalize=True,deltat=line_dur)
              chunk2 = autocorrelate(int_data[np.int(np.floor(len(int_data)/3)):2*np.int(np.floor(len(int_data)/3)),curve] ,normalize=True,deltat=line_dur)
              chunk3 = autocorrelate(int_data[2*np.int(np.floor(len(int_data)/3)):np.int(np.floor(len(int_data)/3))*3,curve] ,normalize=True,deltat=line_dur)
              corr_stack = np.stack ((chunk1[:,1],chunk2[:,1],chunk3[:,1]),1)
              std = np.std (corr_stack[1:std_limit,:],axis=1)
              SNR_list = corr_carpet_temp ['curves'][curve] [1:std_limit,1] / std
              SNR_value = np.average(SNR_list)
              SNR_array = np.append (SNR_array,SNR_value)
              
              # Calculate the brightness from int_data using the moment analysis
              cpm_mom = (np.var (int_data[:,curve])- np.average (int_data[:,curve])) / (np.average (int_data[:,curve])) / (dwell_time*10**(-3)) # in kHz!!!!
              cpm_mom_array = np.append (cpm_mom_array, cpm_mom)  
              
              # Let's do the same for the corrected moment analysis 
              cpm_mom_corr =(( (np.var (int_data[:,curve]) - var_0) / (S* (np.average (int_data[:,curve])-(detector_offset*10**(3)*dwell_time*10**(-6)))) ) - 1) / (dwell_time*10**(-3)) # in kHz  :)
              cpm_mom_corr_array = np.append (cpm_mom_corr_array, cpm_mom_corr)
                    
              
              curve_counter += 1
          
                    
          #Fit file
          lagtime = corr_carpet_temp ['curves'][0][:,0]
          lagtime_max_index = np.min(np.where (lagtime > lagtime_max))
          for iii in range (0,len (corr_carpet_temp ['curves'])):
              res = minimize(residual, param, args=(lagtime[1:lagtime_max_index],corr_carpet_temp ['curves'][iii][1:lagtime_max_index,1])) # ignoring 0 lag time for fitting

              results [fit_counter] = res
              txy_array = np.append (txy_array, res.params ['txy1']._val)
              G0_array = np.append(G0_array, res.params ['GN0']._val)
              fit_counter += 1
              
    
      
    cpm_array = G0_array*int_array 
    
    

    
    
    curve_counter = 0

    # Calculating the nRMSD values based on the fitting residuals
    for curve_counter in range (0, len(results)):
        residuals = (equation_(results [curve_counter].params,corr_carpet ['curves'][curve_counter][:,0]) - corr_carpet ['curves'][curve_counter][:,1])
        residuals = residuals [1:lagtime_max_index] # Now going up to the maximum lag time, excluding 0 lag time. This can of course be set if needed. 
        nRMSD = (np.sqrt (np.sum(residuals**2)/(len(residuals))))/(results [curve_counter].params ['GN0'].value)
        nRMSD_array = np.append (nRMSD_array, nRMSD)


    
    
    #% Save everything to .csv

    data = np.stack ((txy_array,G0_array,cpm_array,int_array,nRMSD_array,SNR_array,cpm_mom_array,cpm_mom_corr_array),1)
    file_name = 'Channel_'+str (ch) + 'FittingParameters_auto.csv'
    out_path = path
    header = ('txy1,GN0,cpm (kHz),Int (kHz),nRMSD,SNR,bri (kHz),bri_corr (kHz)')
    np.savetxt (out_path+file_name, (data), delimiter=',', comments='', header= header)
    
    #Plot and save figures

    fig1 = plt.figure()
##Plot Data-model

    frame3=fig1.add_axes((.1,.3,.8,.6))
    for plot_counter in range (0,len (results)):
        plt.semilogx(corr_carpet ['curves'][plot_counter][1:lagtime_max_index,0], corr_carpet ['curves'][plot_counter][1:lagtime_max_index,1],'ob', alpha = 0.2)
    for plot_counter in range (0,len (results)):    
        plt.semilogx(corr_carpet ['curves'][plot_counter][1:lagtime_max_index,0],equation_(results [plot_counter].params,corr_carpet ['curves'][plot_counter][1:lagtime_max_index,0]),'--r')
        plt.xlabel("tau (ms)")
        plt.ylabel("G(tau)")
        plt.title ('CH_'+ str(ch)+' FCS fit with standard 2D model')
      
        
        
    frame4=fig1.add_axes((.1,.1,.8,.2)) 
    for plot_counter in range (0,len (results)):
        presiduals_standard = (equation_(results [plot_counter].params,corr_carpet ['curves'][plot_counter][1:lagtime_max_index,0]) - corr_carpet ['curves'][plot_counter][1:lagtime_max_index,1])
        plt.semilogx(corr_carpet ['curves'][plot_counter][1:lagtime_max_index,0],presiduals_standard)
        plt.xlabel("tau (ms)")
        plt.ylabel("Residual")
        plt.title ('')
        plt.yticks ([-0.5, 0, 0.5 ])

    

    plt.show()
    
    fig1_path = path+'CH_'+str(ch)+'_Fits.png'   
    fig1=fig1.savefig(fig1_path, bbox_inches='tight')        
    
    #
    # Curve quality 
    fig5, axes = plt.subplots (nrows=1, ncols=2,  figsize=(8, 8))
    ax0, ax1 = axes.flatten ()          
    
    ax0.plot (nRMSD_array,label = 'nRMSD')
    ax1.plot (SNR_array, label = 'SNR')
    ax0.set_xlabel('# of curve')
    ax0.set_ylabel('nRMSD value')
    ax0.set_title (str(ch)+'  nRMSD Analysis')
    ax0.set_ylim ([0, 0.2])
    
    ax1.set_xlabel('# of curve')
    ax1.set_ylabel('SNR value')
    ax1.set_title (str(ch)+'  SNR analysis')
    
    plt.tight_layout()
    
    fig5_path = path+'CH_'+str(ch)+'_CurveQuality.png'   
    fig5.savefig(fig5_path, bbox_inches='tight') 
    
    # Plot BTS histogram (quick and dirty)
    x_min = 0
    x_max = np.log(np.round (max (txy_array)+10))
    binning = 50
    
    fig_bi, ax = plt.subplots(nrows=1, ncols=2,  figsize=(12, 6))
    ax0,ax1 = ax.flatten()
    
    histo = ax0.hist2d (np.log(txy_array), cpm_array, bins = binning, range = [[x_min, x_max],[0, np.round (max (cpm_mom_array))+5]], density = True) 
    ax0.set_xlabel ('Ln (Transit time (ms)/(ms))')
    ax0.set_ylabel ('cpm (kHz)')
    ax0.set_title ('Ch'+ str(ch)+' BTS Histogram cpm')
    plt.colorbar(histo[3], ax = ax0, label = 'Probability (a.u.)')
    
    hist1 = ax1.hist2d (np.log(txy_array), cpm_mom_corr_array, bins = binning, range = [[x_min, x_max],[0, np.round (max (cpm_mom_corr_array))+5]], density = True) 
    ax1.set_xlabel ('Ln (Transit time (ms)/(ms))')
    ax1.set_ylabel ('Bri_corr (kHz)')
    ax1.set_title ('Ch'+str(ch)+' BTS Histogram Brightness')
    plt.colorbar(histo[3], ax = ax1, label = 'Probability (a.u.)')
    
    plt.tight_layout()
              
    fig_bi_path = path+'CH_'+str(ch)+'_Bi-var-hist.png'   
    fig_bi.savefig(fig_bi_path, bbox_inches='tight')   
              

print ('\n All done')


          


              
          
          
          
          
      
      
 