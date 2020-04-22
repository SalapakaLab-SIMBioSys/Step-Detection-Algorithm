# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:48:36 2020

@author: siva_
"""

# =============================================================================
# Description
# =============================================================================

# Code for generating example plots for the paper

# =============================================================================


# =============================================================================
# Imports
# =============================================================================
import random
import numpy as np
import matplotlib.pyplot as plt
import control as ctrl


import os

import sys
sys.path.append('../')

# my toolkits
from APIs.SDACore.CostFunctions import CostFunctions
from APIs.SDACore.SDAHistogram import SDAHistogram
from APIs.SDACore.SDA import SDA
from APIs.Validation.ErrorMetrics import ErrorMetrics
from APIs.Validation.SampleSignals import SampleSignals
from APIs.SDACore.SDADynamics import SDADynamics
from APIs.AlternateTools.AlternateMethods import AlternateMethods
# =============================================================================



altMethodObj = AlternateMethods()
SDynObj = SDADynamics()
sigObj = SampleSignals()
SDAObj = SDA()
ERRObj = ErrorMetrics()
# =============================================================================
# Default signal settings 
# =============================================================================
d = 10 # No. of steps
N_d_ratio = 50 # No. of samples per step
t_avg = 1 # s (Average duration per step in seconds)

dt = t_avg/N_d_ratio # Sampling interval
N = d*N_d_ratio # total number of samples 
L = [1,0.5] # Step size

L_prob = [0.5,0.5] # Step Size probabilities

# Noise standard deviation
sigma = 0.2 #0.3162 # 0.3162 = SNR of 5 for signal of size 1 and a duty ratio of 0.5 

include_dynamics = 0
num = [1]
den = [0.2, 1]
sysd = SDynObj.get_DT_sys(num, den, dt)


# Packaging the standard stepping parameters
step_params = {}
step_params['no_of_steps'] = d
step_params['tot_no_samples'] = N
step_params['step_size'] = L
step_params['sigma'] = sigma
step_params['samples_per_step'] = N_d_ratio
step_params['average_step_duration'] = t_avg
step_params['sensor_dynamics_ON'] = 0
step_params['sensor_num'] = [1]
step_params['sensor_den'] = [0.2, 1]
step_params['sampling_time'] = dt 


# Error tolerance
tol = 0.00
# =============================================================================



# =============================================================================
# Generating noisy stepping data 
# =============================================================================

seed = 9
x = sigObj.get_m_step_sample(step_sizes = L, no_steps = d, no_samples = N, bidirectional = 1,step_size_prob = L_prob, seed = seed)
noise = sigObj.get_gaussian_noise(std = sigma, no_samples = N, seed = seed)

if(include_dynamics == 1):
    [T, x_sensor] = sigObj.sensor_dynamics(num, den, dt, x)
    x_sensor = np.expand_dims(x_sensor, axis = 0)
else:
    x_sensor = x
    sysd = None

y = x_sensor + noise  # Noisy measurements

plt.figure()
plt.plot(y.T)
plt.plot(x.T)
plt.plot(x_sensor.T)

# =============================================================================

#%%

resolution = 0.05
hist_tol = 1e-7 
iterations = 10
[step_estimate, est_array, hist_array]= SDAObj.SDA_dynamics(sigma, y, resolution = resolution,  iterations = iterations, hist_tol = hist_tol , make_plots = 0, diagnostics = 1, sysd = sysd)

#SDAObj.plot_SDA_stages(x.T, est_array, hist_array, no_zeros = 1)


#%%



err_func = ERRObj.get_RMSE    
normalize = 'mean_estimate'
err_array = ERRObj.get_err_array(est_array, x.T, err_func, tol = 0, normalize = 'None')   

save_fig = 1
font = 12
xlim = (0.5,1.5)
ylim = (0, 0.022)

SDAHistObj = SDAHistogram()

arr_len = np.size(est_array,0)

# Plot clean Histograms
for ii in range(0,arr_len,1):    
    [bin_centers, bin_heights, bin_width, bin_boundaries] = SDAHistObj.plot_hist(hist_array[ii][0], hist_array[ii][1], no_zeros= 1, suppress_plots = 1)     
    
    

    plt.figure()        
    plt.bar(bin_centers, bin_heights, width = bin_width, color=(0.05, 0.05, 0.05, 0.2),\
            edgecolor='blue', align = 'center', lw = 2)
    
    
    txt_str = 'Iteration = ' + str(ii) +'\n'+\
                'NRMSE = ' + str(np.around(err_array[ii], 4 ))
    
    plt.yticks(fontsize = font)
    plt.xticks(bin_boundaries, fontsize = font, rotation = 'horizontal')
    plt.ylim(ylim)
    plt.xlim(xlim)
    plt.xlabel('Step Size', fontsize = font)
    plt.ylabel('Probability', fontsize = font)
    ax = plt.gca()
    ax.yaxis.grid()
    
    ax.text(xlim[1]*0.8, ylim[1]*0.8, txt_str, fontsize = font)
    
    plt.tight_layout()
    plt.show() 
    
    
    if(save_fig == 1):        
        plt.savefig('Stages/StepHistogram_Stage_'+str(ii)+'.png')
        plt.savefig('Stages/StepHistogram_Stage_'+str(ii)+'.svg')
        
     

#%%
# Plot clean sample plots
for ii in range(0,arr_len,1):
    plt.figure()
    plt.plot(y.T, color = (0, 0.5, 0, 0.5) )
    plt.plot(x.T, lw = 2, color = 'r')
    plt.plot(est_array[ii], ls='--',lw = 2, color = 'b')
    

    plt.yticks(fontsize = font)
    plt.xticks(fontsize = font, rotation = 'horizontal')
    plt.xlabel('Time', fontsize = font)
    plt.ylabel('Magnitude', fontsize = font)
    plt.show() 
    
    plt.grid()
    
    if(save_fig == 1):        
        plt.savefig('Stages/StepFit_Stage_'+str(ii)+'.png')
        plt.savefig('Stages/StepFit_Stage_'+str(ii)+'.svg')

#%%

#step_sizes = [1,5]
#step_size_prob = [0.5, 0.5]
#bidirectional = 1
#no_steps = 5
#no_samples = 1000
#seed = 9
#bidirec_prob = 0.5
#
## Initialize both randomizations with seed if a valid value is given 
#if(seed is not None):
#    np.random.seed(seed)
#
#seed_arr = np.random.randint(2**31, size = 2)
#
#L = step_sizes
#d = no_steps
#N = no_samples
#
#np.random.seed(seed_arr[0])
#s = np.floor(np.random.rand(no_steps)*N)
#s = s.astype(int)
#
#x = np.matrix(np.zeros((1,N)))
#
#
#np.random.seed(seed_arr[1])
#for ii in range(d):
#    
#    L_local = np.random.choice(L, p = step_size_prob)
#    if(bidirectional == 1):        
#        x[0,s[ii]] = L_local*np.sign(bidirec_prob-np.random.rand(1))
#    else:
#        x[0,s[ii]] = L_local
#
#x = np.cumsum(x)   # Underlying signal
#
#step_signal = x        
#
#
#
#
#plt.figure()
#plt.plot(x.T)
