# -*- coding: utf-8 -*-

# =============================================================================
# Description
# =============================================================================

 # Code for the Step Detection Algorithm (SDA) and its variants

# =============================================================================


# =============================================================================
# Imports
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import copy 

import sys
sys.path.append('../../')

# my toolkits
from APIs.SDACore.CostFunctions import CostFunctions
from APIs.SDACore.SDAHistogram import SDAHistogram

# =============================================================================

#%%


# =============================================================================
# Class containing the code for the SDA and variants
# =============================================================================
class SDA():
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):
        a = 1
   
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------        
    # SDA standard
    #--------------------------------------------------------------------------        
    def SDA_standard(self, noise_std, y, resolution = 0.1, iterations = 5, diagnostics = 0, make_plots = 0):
                
        y_range = np.ptp(y) # Range of the dataset
        #resol_N = 100 # No. of points for reslution
        #resol = abs(R/(resol_N-1)) # Resolution to be used
        
        resol = resolution;
        resol_N = np.int(np.ceil(abs(y_range/resol)))
        seq_len = np.size(y,1) # Signal length
        
        no_bins = resol_N

        est_array = []
        hist_array = []
                
        SDAHistObj = SDAHistogram() 
        #diff_y = np.diff(y)
        #bin_boundaries = SDAHistObj.get_boundaries(-np.max(np.abs(diff_y)), np.max(np.abs(diff_y)), no_bins)
        bin_boundaries = SDAHistObj.get_boundaries(-np.max(np.abs(y)), np.max(np.abs(y)), no_bins)
        
        for iter_no in range(0,iterations,1):
            
            if(iter_no == 0):
                step_hist = []
                 
        
            cost_array = np.matrix(np.zeros((resol_N,seq_len))) # Array to store cost functions
            min_x_array = np.matrix(np.zeros((resol_N,seq_len))) # Array to store x_candidates that produces least cost for that step
            min_x_loc_array = np.matrix(np.zeros((resol_N,seq_len))) 
            
            # Computing cost for different path possibilities
            costObj = CostFunctions()
            est_noise_std = noise_std # Estimate noise std
            x_candidates = np.matrix(np.linspace(np.min(y),np.max(y),resol_N)).T
            cost_prev = np.matrix(np.zeros((resol_N,1)))
            for ii in range(0,seq_len,1):    
                
                x_measured = y[0,ii]
                
                if(iter_no == 0):
                    [local_best_cost, local_best_estimate, local_best_locs] = costObj.step_cost_gen(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist)
                else:
                    [local_best_cost, local_best_estimate, local_best_locs] = costObj.step_cost_gen(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist)
                    
                    
                cost_array[:,ii] =  local_best_cost
                min_x_array[:,ii] = local_best_estimate
                min_x_loc_array[:,ii] = local_best_locs
                cost_prev = local_best_cost    
            
            
            # Selecting the optimal path
            x_estimate = np.matrix(np.zeros((seq_len,1)))
            
            min_x_loc = np.argmin(cost_array[:,seq_len-1])
            x_estimate[seq_len-1,0] = x_candidates[min_x_loc,0]
            min_x_loc = min_x_loc_array[min_x_loc,seq_len-1]
            min_x_loc = np.int(min_x_loc)
            
            for ii in range(seq_len-2,-1,-1):
                x_estimate[ii,0] = x_candidates[min_x_loc,0]
                min_x_loc = min_x_loc_array[min_x_loc,ii]
                min_x_loc = np.int(min_x_loc)
                        
            
            step_estimate = np.diff(x_estimate.T) 
            #step_estimate_no_0 = step_estimate[step_estimate!=0]
            
            # Unidirectional steps
            #step_hist = SDAHistObj.build_USC_hist(-np.abs(step_estimate),bin_boundaries, density = True)
            
            # Bidirectional steps
            #step_hist = SDAHistObj.build_USC_hist(-step_estimate, bin_boundaries, density = True)
            step_hist = SDAHistObj.build_USC_hist(step_estimate, bin_boundaries, density = True)

            #-- Make diagnostic plots if required --#
            if(make_plots == 1):
                # Measurements and Estimates 
                plt.figure()
                plt.plot(y.T)
                plt.plot(x_estimate)
                                
                SDAHistObj.plot_hist(step_hist[0], step_hist[1])
            
            #-- Collect diagnostic data if required --#
            if(diagnostics == 1):
                est_array.append(x_estimate)
                hist_array.append(step_hist)
                         
        return(x_estimate, est_array, hist_array)
    #--------------------------------------------------------------------------    
    
    
    #--------------------------------------------------------------------------        
    # SDA with sensor dynamics (optional, if dynamics are 1, reverts to standard SDA)
    #--------------------------------------------------------------------------        
    def SDA_dynamics(self, noise_std, y, resolution = 0.1, hist_tol = 0.001, iterations = 5, diagnostics = 0, make_plots = 0, sysd = None):
                
        y_range = np.ptp(y) # Range of the dataset
        #resol_N = 100 # No. of points for reslution
        #resol = abs(R/(resol_N-1)) # Resolution to be used
        
        resol = resolution;
        resol_N = np.int(np.ceil(abs(y_range/resol)))
        seq_len = np.size(y,1) # Signal length
        
        no_bins = resol_N

        est_array = []
        hist_array = []
        hist_ks_dist = []
                
        SDAHistObj = SDAHistogram() 
        #diff_y = np.diff(y)
        bin_boundaries = SDAHistObj.get_boundaries(-np.max(np.abs(y)), np.max(np.abs(y)), no_bins)
        
        for iter_no in range(0,iterations,1):
            
            if(iter_no == 0):
                step_hist = []
                 
        
            cost_array = np.matrix(np.zeros((resol_N,seq_len))) # Array to store cost functions
            min_x_array = np.matrix(np.zeros((resol_N,seq_len))) # Array to store x_candidates that produces least cost for that step
            min_x_loc_array = np.matrix(np.zeros((resol_N,seq_len))) 
            
            # Computing cost for different path possibilities
            costObj = CostFunctions()
            est_noise_std = noise_std # Estimate noise std
            x_candidates = np.matrix(np.linspace(np.min(y),np.max(y),resol_N)).T
            cost_prev = np.matrix(np.zeros((resol_N,1)))
            
            
            # If sensor dynamics are given
            if(sysd != None):
                z_init = np.zeros((np.shape(sysd.A)[0],1))          
                z_k_best_arr = []       
                for kk in range(0,np.size(x_candidates),1):                                
                    z_k_best_arr.append(z_init)                    
                                 
            
            for ii in range(0,seq_len,1):    
                
                x_measured = y[0,ii]
                
                #if(iter_no == 0):
                #    [local_best_cost, local_best_estimate, local_best_locs] = costObj.step_cost_gen(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist)
                #else:
                #    [local_best_cost, local_best_estimate, local_best_locs] = costObj.step_cost_gen(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist)
                
                if(sysd != None):
                    [local_best_cost, local_best_estimate, local_best_locs, z_k_best_arr] = costObj.step_cost_dyn_fast(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist, sysd, z_k_best_arr)
                else:
                    [local_best_cost, local_best_estimate, local_best_locs] = costObj.step_cost_gen(x_candidates, x_measured, est_noise_std, cost_prev, iter_no, step_hist)
                    
                    
                cost_array[:,ii] =  local_best_cost
                min_x_array[:,ii] = local_best_estimate
                min_x_loc_array[:,ii] = local_best_locs
                cost_prev = local_best_cost    
            
            
            # Selecting the optimal path
            x_estimate = np.matrix(np.zeros((seq_len,1)))
            
            min_x_loc = np.argmin(cost_array[:,seq_len-1])
            x_estimate[seq_len-1,0] = x_candidates[min_x_loc,0]
            min_x_loc = min_x_loc_array[min_x_loc,seq_len-1]
            min_x_loc = np.int(min_x_loc)
            
            for ii in range(seq_len-2,-1,-1):
                x_estimate[ii,0] = x_candidates[min_x_loc,0]
                min_x_loc = min_x_loc_array[min_x_loc,ii]
                min_x_loc = np.int(min_x_loc)
                   
            # Compute steps from estimate                 
            step_estimate = np.diff(x_estimate.T) 

            # Computing the histogram from the step estimate
            step_hist = SDAHistObj.build_USC_hist(step_estimate, bin_boundaries, density = True)

            # Compute the KS distance from the previous Histogram
            cumul_dist = np.cumsum(step_hist[0])
            if(iter_no == 0):
                cumul_dist_prev = np.zeros(np.shape(cumul_dist))

            ks_dist = np.max(np.abs(cumul_dist - cumul_dist_prev))
            hist_ks_dist.append(ks_dist)          
            cumul_dist_prev = copy.deepcopy(cumul_dist)
                                        
            
            #-- Make diagnostic plots if required --#
            if(make_plots == 1):
                # Measurements and Estimates 
                plt.figure()
                plt.plot(y.T)
                plt.plot(x_estimate)
                                
                SDAHistObj.plot_hist(step_hist[0], step_hist[1])
            
            #-- Collect diagnostic data if required --#
            if(diagnostics == 1):
                est_array.append(x_estimate)
                hist_array.append(step_hist)
            
            print('Iteration count = ' + str(iter_no))
            print('Hist Tol = ' + str(ks_dist)+'\n')
            
            if(ks_dist <= hist_tol):
                print('Hist tolerance reached! \n')
                print('Hist Tol = ' + str(ks_dist)+'\n')
                print('Iteration count = ' + str(iter_no))
                
                if(iter_no >=5 ):
                    break
                         
        return(x_estimate, est_array, hist_array)
    #--------------------------------------------------------------------------        
    
    
    #--------------------------------------------------------------------------        
    # Plot the SDA stages
    #--------------------------------------------------------------------------        
    def plot_SDA_stages(self, x_actual, est_array, hist_array, no_zeros=0):    
        
        SDAHistObj = SDAHistogram()
        
        arr_len = np.size(est_array,0)
        
        for ii in range(0,arr_len,1):
            plt.figure()
            plt.plot(x_actual)
            plt.plot(est_array[ii])
            
        
            SDAHistObj.plot_hist(hist_array[ii][0], hist_array[ii][1], no_zeros= no_zeros)        
        return()        
    #--------------------------------------------------------------------------        

# =============================================================================

#%%
'''

# =============================================================================
# Generating noisy stepping data 
# =============================================================================
d = 2 # No. of steps
dt = 0.1 # Sampling interval
N = 200 # total number of samples 
L = 1 # Step size

# Noise standard deviation
sigma = 0.5 #0.3162 = SNR of 5 for signal of size 1 and a duty ratio of 0.5 


sigObj = SampleSignals()

x = sigObj.get_step_sample(step_size = L, no_steps = d, no_samples = N)
noise = sigObj.get_gaussian_noise(std = sigma, no_samples = N)
#noise = sigObj.get_gamma_noise(shape = 5, scale = 1, no_samples = N)

y = x + noise  # Noisy measurements



# =============================================================================



SDAObj = SDA()

[step_estimate, est_array, hist_array]= SDAObj.SDA_standard(sigma, y, make_plots = 0, diagnostics = 1)


plt.figure()
plt.plot(step_estimate)
plt.plot(x.T, 'r--')

'''