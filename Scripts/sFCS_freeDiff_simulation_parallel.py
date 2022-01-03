# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:37:50 2020

@author: Falk Schneider
"""

# sFCS_freeDiff_simulation_parallel
#
# Python script to simulate free diffusion and analyse it with scanning FCS (sFCS).
# This is a cleaned-up version to accompany the BTS histogram manuscript. 
#
# The code is based on "sFCS_free_parallel_v03_forServer"(for internal reference only).
# The simulations were developed further from the ones used in 10.1016/j.ymeth.2017.09.010, in 10.1021/acsnano.8b04080, and in 10.1021/acs.nanolett.8b01190. A lot of credit and thanks go to Dominic Waithe who has done a lot of the initial work on the simulations for the previous papers. 
# This script displays an attempt to make all of this open-source available. 
#
#
# You can change most of the simulation parameters (like diffusion coefficient or simulation time or brightness at the top of the script (after the functions)).
# The scrip requires the use of iPyParallel to distribute computational load to multiple CPUs!
# Please set an output path for the resulting .tif files containing the intensity fluctuations
# The sFCS data can be analysed using the FoCuS_scan software package (10.1016/j.ymeth.2017.09.010 or https://github.com/dwaithe/FCS_scanning_correlator).
# Alternatively, the data can be analysed with raw2fittingpara.py (which is based on FoCuS_scan).
#
# Falk Schneider 02/01/2022 @faldalf
#
# Questions and feedback are welcome

#%% Imports

import time
from scipy.stats import norm
import tifffile
import numpy as np 
import ipyparallel as ipp

#%% Functions
def calculate_psf_cut (psfs,distance,cut=4):
    # calculate excitation psf
    # Allows "cutting" of the PSF (setting it to 0) at some multiple of the FWHMs "cut"
    #
    # psfs: array of FWHMs of the PSFs to compute
    # distance: maximum_x (related to size of the simulationn box)
    # cut is multiplier of FWHM from which on every intensity is 0. default is 4
    #
    psf = {}
    psf['FWHMs'] = psfs # np.array with the FWHMs of the Gaussian PSFs to generate 
    psf['pixel_size'] = 1;
    psf['ri'] = np.meshgrid(np.arange(0,distance,psf['pixel_size']))[0]; # generate x-values

    psf['number_FWHMs'] = psf['FWHMs'].__len__()
    psf['V'] = {}
    for ki in range(0, psf['number_FWHMs']):
        psf['V'][ki] = 2.0**(- psf['ri']**2 / (psf['FWHMs'][ki]/2.0)**2); #Gaussian distribution
        psf['V'][ki][np.int(psf['FWHMs'][ki]*cut):] = 0 # Cutting off.
    return psf
    

    
def save_to_tiff(out_point,fwhm,file_name):
    # This function converts the data_store dict to a tiff file
    #
    # out_point is the data_store dict.
    # fwhm is index of the intensity trace meaning which FWHM on psf['FWHMs']    
    matrix_im = np.zeros((out_point.__len__(),out_point[0]['trace'][fwhm].shape[0]-1))
    
    for b in out_point:
        
        matrix_im[b,:] = out_point[b]['trace'][fwhm][:out_point[0]['trace'][fwhm].shape[0]-1]
    
    with tifffile.TiffWriter(file_name, bigtiff=True) as tif:
        tif.save(matrix_im.T.astype(np.float64))
    return 

def tracks_to_int (total_sim_time,time_step,num_of_mol,D, width, height, pos_y, pos_x, psf,maxphotonflux,steps_to_jump):
   #
   # Actual simulation function to generate random walks and perform the sFCS measurement
   # Uses parallel computing
   # 
   # We need to reimport the modules on the engines for parallel computing. 
    import numpy as np 
    from scipy.stats import norm
    
    num_of_steps = int(round(float(total_sim_time)/float(time_step),0))

    # Calculates length scales
    scale_in = np.sqrt(2.0 * (float(D)*1e3) * float(time_step)) #[nm/time_step]
    dT = time_step/1000 # convert to seconds. 
    
    
    # Create data store
    cc = 0 # Pixel 0 to end of line
    data_store = {}
    temp = len(np.round(np.arange(cc,num_of_steps,steps_to_jump),0).astype(np.int32))
    for psy, psx in zip(pos_y,pos_x):
        data_store[cc] = {}
        data_store[cc]['trace'] = {}
        data_store[cc]['psy'] = psy
        data_store[cc]['psx'] = psx
        for ki in range(0, psf['number_FWHMs']):
            data_store[cc]['trace'][ki]= np.zeros ([temp])
        cc += 1
        
   # Now we walk molecucle by molecule and get the intensity added up, molecule by molecule
    for b in range (0,num_of_mol):
        # Initialise molecule
        start_coord_x = (np.random.uniform(0.0, 1.0, 1))*width
        start_coord_y = (np.random.uniform(0.0, 1.0, 1))*height
    
        track_arr = np.zeros((2,num_of_steps))
        track_arr[0,0] = start_coord_y
        track_arr[1,0] = start_coord_x
        rand_in  = norm.rvs(size=[2,num_of_steps])*scale_in
        track_arr[:,1:] += rand_in[:,1:]
        track_arr = np.cumsum(track_arr,1)
        mod = np.zeros((track_arr.shape)) # for periodic boundary.
        mod[0,:] = np.floor(track_arr[0,:].astype(np.float64)/height)
        mod[1,:] = np.floor(track_arr[1,:].astype(np.float64)/width)
        track_arr = np.array(track_arr-([mod[0,:]*height,mod[1,:]*width]))
            
        # Here comes the sFCS part :) 
        # Just integrating the PSFs
        cc = 0 # cc ensures that the pixels are sampled one after each other and not at the same time.
        for psy, psx in zip(pos_y,pos_x): #This function returns a list of tuples, where the i-th tuple contains the i-th element from each of the argument sequences or iterables. 
            psf['trace'] ={}
            for ki in range(0, psf['number_FWHMs']):
                psf['trace'][ki] = 0                
                to_sample = np.round(np.arange(cc,track_arr.shape[1],steps_to_jump),0).astype(np.int32)
                to_sample[-1] = to_sample [-1]-1
                track_x = track_arr[1][to_sample] # That generate a "down-sampled" track array depending on the scanning frequency (and time steps). Say only every to_sampleth time step will be used
                track_y = track_arr[0][to_sample]
                intensity = psf['V'][ki][np.round(np.sqrt((track_x-psx)**2+(track_y-psy)**2),0).astype(np.int32)]*dT*maxphotonflux # Intensities across space based on PSF
                photon_trace = np.random.poisson (intensity,len(intensity)) # Discretize intensity to obtain photon counting data
                data_store[cc]['trace'][ki] += photon_trace
                    
            
            cc +=1 
        
        
  
            
    #That's it.          
    

    return data_store # intensity data over space



#%% The simulation parameters and the loop.
# Change parameters as appropriate

#Output directory - where the out data goes:
out_path_tiff = '\\Users\\Falk Schneider\\Desktop\\Temp\\test\\'

# Parallel computing: 
#https://ipyparallel.readthedocs.io/en/latest/direct.html#starting-the-ipython-controller-and-engines
#http://people.duke.edu/~ccc14/sta-663-2018/notebooks/S14D_IPyParallel.html
#ipcluster start -n 4 # don't forget to start clusters before running parallel code. 
clients = ipp.Client () # Start the engines
dview = clients[:] # control over all running engines
print ('I see .... we are trying to do parallel computing ... Alright then... ')
print ('You are using ',np.str(len(clients)),' engines at the moment')

#Time simulation
time_start = time.time()

# All the input parameters for our simulation
# These parameters do not change over the long loop or number of simulations. 

width = 8000.0 #The width of the simulation in nm.
height = 5000.0 #The height of the simulation in nm. 
total_sim_time = 10e03 # Total time in ms
num_of_mol = 40 #Number of molecules in simulation.
time_step = 0.004 #Time step in ms = pixel dwell time!!! (time increment between to pixels)
D = [0.5] #diffusion coefficient [um^2/s]; SHOULD BE A LIST

line_time = 2000 #Hz - scanning frequency in Hz = overall NOT per channel.
line_dur = ((1.0/line_time)*1000.0) # line time in ms
scan_length = 5000 #nm 
num_of_pixels = 50 # along the line
pixel_size = scan_length/num_of_pixels # 
steps_to_jump = line_dur/(time_step) # down_sampling of the tracks by frequency ... only look every ith time
px_gap = (scan_length)/num_of_pixels 
maxphotonflux = [25e03]#, 50e03, 100e03] # List of brightness -> 1/2 of measured cpm

foci_to_generate = np.array([240]) # Array of PSFs to generate (indicate FWHM of every PSF).

size = max(width+1000,height+1000)

# cut PSF cut = 2 means cutting at FWHM*2. Remember though FWHM is beam diameter. 
# to cut actually at FWHM set cut = 0.5
psf = calculate_psf_cut (foci_to_generate,size,cut=3)

pos_x = np.arange((np.max(width)-scan_length)//2,scan_length+(np.max(width)-scan_length)//2,px_gap) # gives the position of the psf along the scan line
pos_y = [height//2]*pos_x.shape[0] # This defines the scanning line. constant y
    
# Definitions for parallel computing
dview.push ({'tracks_to_int': tracks_to_int})       
# Migrate the name space
dview['total_sim_time'] = total_sim_time
dview['time_step'] = time_step
dview['width'] = width
dview['height'] = height
dview['pos_y'] = pos_y
dview['pos_x'] = pos_x
dview['psf'] = psf
dview['steps_to_jump'] = steps_to_jump
num_of_mol_per_engine = list ((np.ones (len(clients))*num_of_mol//len(clients)).astype(int)) # Dividing the computational load by splitting up the total number of molecules
dview['num_of_mol_per_engine'] = num_of_mol_per_engine


#Looooping
numofsimu = 1
jc = 1

# Generate the simulations: 
    
# Loop over repetitions
for N in range (0,numofsimu):
    print ('Starting simulation round number:',jc)
    
    # Loop over maxphotonflux
    for cpm in maxphotonflux:
        dview['cpm'] = cpm
        
        # Loop over D
        for diff in D:
            dview['diff'] = diff
            print ('Starting simulation with cpm = ',cpm, ' Hz', ' and with D = ',D, ' um^2/s')
            parallel_data_stores = dview.map_sync (lambda n: tracks_to_int (total_sim_time,time_step,n,diff, width, height, pos_y, pos_x, psf,cpm,steps_to_jump), num_of_mol_per_engine)
               
            # Now we put all together. 
            # Create empty data_store
            num_of_steps = int(round(float(total_sim_time)/float(time_step),0)) 
            cc = 0 # Pixel 0 to end of line
            data_store = {}
            temp = len(np.round(np.arange(cc,num_of_steps,steps_to_jump),0).astype(np.int32))
            for psy, psx in zip(pos_y,pos_x):
                data_store[cc] = {}
                data_store[cc]['trace'] = {}
                data_store[cc]['psy'] = psy
                data_store[cc]['psx'] = psx
                for ki in range(0, psf['number_FWHMs']):
                    data_store[cc]['trace'][ki]= np.zeros ([temp])
                cc += 1
            # Put all the parallel data_stores additive into the final data_store
            cc = 0
            for i in range (0,len(parallel_data_stores)):
                cc = 0
                for psy, psx in zip(pos_y,pos_x):
                    for ki in range(0, psf['number_FWHMs']):
                        data_store [cc]['trace'][ki] += parallel_data_stores [i] [cc]['trace'][ki]
                    cc += 1
            
            # Save it all to tiff. 
            
            
            file_name_240 = 'MaxPhoto'+str(cpm//1000).replace('.','p')+'_kHz'+'_D_'+str(diff).replace('.','p')+'_acquTime_'+str(total_sim_time//1000)+'s'+'_numofmol_is_'+str(num_of_mol).replace('.','p')+'_at_'+str (foci_to_generate[0]).replace('.','p')+'_nm_repetition_'+str(N)+'.tif'
            save_to_tiff(data_store,0,out_path_tiff+file_name_240)
            
            # Save it all to tiff. 
            # down-sampling to half the frequency
            # frequency for correlation will be half of simulation frequency!
            for index in range (0,len(data_store)):
                for ki in range(0, psf['number_FWHMs']):
                    data_store [index] ['trace'] [ki] = data_store [index]['trace'][ki] [ki:len(data_store [index]['trace'][ki]):2]
            
            # save the down-sampled file
            file_name_240 = 'MaxPhoto'+str(cpm//1000).replace('.','p')+'_kHz'+'_D_'+str(diff).replace('.','p')+'_acquTime_'+str(total_sim_time//1000)+'s'+'_numofmol_is_'+str(num_of_mol).replace('.','p')+'_at_'+str (foci_to_generate[0]).replace('.','p')+'_nm_repetition_'+str(N)+'half_frequency'+'.tif'
            save_to_tiff(data_store,0,out_path_tiff+file_name_240)
                
        
            
    time_elapsed = np.round ((time.time() - time_start),1)
    print ('Elapsed time (Python_Parallel): ', np.str(time_elapsed),' seconds for one repetition (all cpm and D).')
         
          
    ccc = jc/numofsimu*100
    print (ccc,'% of simulations done')
    jc = jc + 1
            
    



print ('All done!')
    
