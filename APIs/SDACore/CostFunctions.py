# -*- coding: utf-8 -*-

# =============================================================================
# Description
# =============================================================================

 # Code for SDA cost functions and variants

# =============================================================================
 
 
# =============================================================================
# Imports
# ============================================================================= 
import numpy as np

import sys
sys.path.append('../../')

from APIs.SDACore.SDAHistogram import SDAHistogram
from APIs.SDACore.SDADynamics import SDADynamics
# =============================================================================

# =============================================================================
# Class containing cost functions 
# =============================================================================
class CostFunctions():
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):        
        self.SDAHistObj = SDAHistogram() 
        self.SDynObj = SDADynamics()
    #--------------------------------------------------------------------------    

    #--------------------------------------------------------------------------
    # Cost function for a single iteration
    #--------------------------------------------------------------------------
    def step_cost(self, x_candidates, x_measured, noise_std, cost_prev_step):
                        
        energy_cost = np.power((x_candidates - x_measured),2) 
        
        penalty_weight = (9*(noise_std**2))    
        x_prev_step = x_candidates.T
        penalty_possibilities = np.multiply((x_prev_step != x_candidates), penalty_weight)
        
        cost_possibilities = energy_cost + penalty_possibilities
        
        local_best_cost = np.min(cost_possibilities,0).T
        local_best_estimate_locs = np.argmin(cost_possibilities,0).T
        local_best_estimate = x_candidates[local_best_estimate_locs]
        local_best_estimate = local_best_estimate.reshape((np.size(local_best_estimate,0),1))
        
        
        return(local_best_cost, local_best_estimate, local_best_estimate_locs)
    #--------------------------------------------------------------------------
    
    
    
    
    #--------------------------------------------------------------------------
    # Cost function for multiple iterations
    # Based on histogram weigths, without interpolation, assumes 0 step size not 
    # included while computing histogram. Step size of 0 given a penalty of 0.    
    #--------------------------------------------------------------------------
    def step_cost_gen(self, x_candidates, x_measured, noise_std, cost_prev_step, iter_no, hist_data):
                        
        energy_cost = np.power((x_candidates - x_measured),2) 
        x_prev_step = x_candidates.T
        
        if(iter_no == 0):
            penalty_weight = 9*(noise_std**2)      
            penalty_possibilities = np.multiply((x_prev_step != x_candidates), penalty_weight)
                                       
        else:                                              
            step_sizes = x_candidates-x_prev_step
            step_prob = self.SDAHistObj.SDA_hist(step_sizes, hist_data[0], hist_data[1])
            
            #-- Fixing the issue of finite bin around the zero --#
            
            boundary_diffs = np.diff(hist_data[1])
            resol_act = np.min(np.abs(boundary_diffs))            
            zero_dist_chk = (np.abs(step_sizes) <= resol_act) & (step_sizes !=0) 
            step_prob[zero_dist_chk] = 0
            
            #----------------------------------------------------#
            
            penalty_weight = (-2*(noise_std**2))*np.log(step_prob)
            penalty_possibilities = penalty_weight
            
        
        penalty_possibilities[np.isnan(penalty_possibilities)] = 0  
         
        #cost_possibilities = energy_cost + penalty_possibilities + cost_prev_step        
        
        #local_best_cost = np.min(cost_possibilities,0).T
        #local_best_estimate_locs = np.argmin(cost_possibilities,0).T

        cost_possibilities = energy_cost + penalty_possibilities + cost_prev_step.T
        
        local_best_cost = np.min(cost_possibilities,1)
        local_best_estimate_locs = np.argmin(cost_possibilities,1)       
        
        local_best_estimate = x_candidates[local_best_estimate_locs]
        local_best_estimate = local_best_estimate.reshape((np.size(local_best_estimate,0),1))
  
        
        return(local_best_cost, local_best_estimate, local_best_estimate_locs)
    #--------------------------------------------------------------------------    



    #--------------------------------------------------------------------------
    # Cost function for multiple iterations
    # Based on histogram weigths, without interpolation, assumes 0 step size not 
    # included while computing histogram. Step size of 0 given a penalty of 0.    
    #--------------------------------------------------------------------------
    def step_cost_dyn(self, x_candidates, x_measured, noise_std, cost_prev_step, iter_no, hist_data, sysd, z_k_best_arr):
                        
        z_measured = x_measured

        z_k_tot_arr = []
        z_candidates_tot_arr = []
        u_k_arr = []
        for jj in range(0, np.size(x_candidates),1):
                    
            u_k_arr.append(x_candidates[jj])
            
            z_k_arr = []   
            z_candidates_arr = []
            for kk in range(0,np.size(x_candidates),1):                    
            
                [z_k, z_candidates] = self.SDynObj.get_1step_update(sysd, z_k_best_arr[kk], u_k_arr[jj])
                
                z_candidates_arr.append(z_candidates)
                z_k_arr.append(z_k)
            
            z_candidates_tot_arr.append(z_candidates_arr)
            z_k_tot_arr.append(z_k_arr)
            
        #z_k_1_tot_arr = z_k_tot_arr     
        
        #energy_cost = np.power((z_candidates - z_measured),2)
        
        z_can_tot = np.matrix(np.array(z_candidates_tot_arr))
        
        energy_cost = np.power((z_can_tot-z_measured),2)
        
        
           
        #energy_cost = np.power((x_candidates - x_measured),2) 
        x_prev_step = x_candidates.T
        
        if(iter_no == 0):
            penalty_weight = 9*(noise_std**2)      
            penalty_possibilities = np.multiply((x_prev_step != x_candidates), penalty_weight)
                                       
        else:                                              
            step_sizes = x_candidates-x_prev_step
            step_prob = self.SDAHistObj.SDA_hist(step_sizes, hist_data[0], hist_data[1])
            
            #-- Fixing the issue of finite bin around the zero --#
            
            boundary_diffs = np.diff(hist_data[1])
            resol_act = np.min(np.abs(boundary_diffs))            
            zero_dist_chk = (np.abs(step_sizes) <= resol_act) & (step_sizes !=0) 
            step_prob[zero_dist_chk] = 0
            
            #----------------------------------------------------#
            
            penalty_weight = (-2*(noise_std**2))*np.log(step_prob)
            penalty_possibilities = penalty_weight
            
      
        penalty_possibilities[np.isnan(penalty_possibilities)] = 0  
        cost_possibilities = energy_cost + penalty_possibilities + cost_prev_step.T
        
        
        #local_best_cost = np.min(cost_possibilities,0).T
        #local_best_estimate_locs = np.argmin(cost_possibilities,0).T
        
        local_best_cost = np.min(cost_possibilities,1)
        local_best_estimate_locs = np.argmin(cost_possibilities,1)
                    
        local_best_estimate = x_candidates[local_best_estimate_locs]
        local_best_estimate = local_best_estimate.reshape((np.size(local_best_estimate,0),1))
      
        #local_best_locs = local_best_estimate_locs
        
        z_k_best_arr = []
        for kk in range(0,np.size(x_candidates),1):                
            z_k_best_arr.append(z_k_tot_arr[kk][local_best_estimate_locs[kk][0,0]])
                
        
        return(local_best_cost, local_best_estimate, local_best_estimate_locs, z_k_best_arr)

    #--------------------------------------------------------------------------



    #--------------------------------------------------------------------------
    # Cost function for multiple iterations
    # Based on histogram weigths, without interpolation, assumes 0 step size not 
    # included while computing histogram. Step size of 0 given a penalty of 0.    
    #--------------------------------------------------------------------------
    def step_cost_dyn_fast(self, x_candidates, x_measured, noise_std, cost_prev_step, iter_no, hist_data, sysd, z_k_best_arr):                        

        z_measured = x_measured

        #--------------- Vectorized Dynamics Update ------------------#
        
        dim_z_k = np.shape(sysd.A)[0]
        dim_x_can = np.size(x_candidates)

        #-- Initialize --#      
        if(dim_z_k == 1):
            zz_k_best_arr = np.reshape(np.array(z_k_best_arr),(dim_z_k, dim_x_can))
        else:
            zz_k_best_arr = np.matrix(np.array(z_k_best_arr)).T        
        
        
        #-- Next state --#
        zzx = sysd.A * zz_k_best_arr 
        zzu = sysd.B * x_candidates.T
        
        zz_k_arr = []
        for jj in range(0, dim_z_k, 1):                
            zz_k_arr.append(zzx[jj] + zzu[jj].T)
        zz_k_arr = np.array(zz_k_arr)
        

        #-- Output from state --#
        zz_canx = []
        for jj in range(0,dim_z_k,1):                
            zz_canx.append(np.multiply(sysd.C[0,jj], zz_k_arr[jj,:,:]))            

        zz_canx = np.array(zz_canx)
        zz_canx = np.sum(zz_canx, axis =0)

        zz_canu = sysd.D * x_candidates.T
        z_can_tot = zz_canx + zz_canu.T
        #--------------------------------------------------------#
        
        
        energy_cost = np.power((z_can_tot-z_measured),2)
        
        x_prev_step = x_candidates.T
        
        if(iter_no == 0):
            penalty_weight = 9*(noise_std**2)      
            penalty_possibilities = np.multiply((x_prev_step != x_candidates), penalty_weight)
                                       
        else:                                              
            step_sizes = x_candidates-x_prev_step
            step_prob = self.SDAHistObj.SDA_hist(step_sizes, hist_data[0], hist_data[1])
            
            #-- Fixing the issue of finite bin around the zero --#
            
            boundary_diffs = np.diff(hist_data[1])
            resol_act = np.min(np.abs(boundary_diffs))            
            zero_dist_chk = (np.abs(step_sizes) <= resol_act) & (step_sizes !=0) 
            step_prob[zero_dist_chk] = 0
            
            #----------------------------------------------------#
            
            penalty_weight = (-2*(noise_std**2))*np.log(step_prob)
            penalty_possibilities = penalty_weight
            
      
        penalty_possibilities[np.isnan(penalty_possibilities)] = 0  
        cost_possibilities = energy_cost + penalty_possibilities + cost_prev_step.T
                
        local_best_cost = np.min(cost_possibilities,1)
        local_best_estimate_locs = np.argmin(cost_possibilities,1)
                    
        local_best_estimate = x_candidates[local_best_estimate_locs]
        local_best_estimate = local_best_estimate.reshape((np.size(local_best_estimate,0),1))
      
        z_k_best_arr = []
        for kk in range(0,np.size(x_candidates),1):
            z_k_best_arr.append(np.expand_dims(zz_k_arr[:,kk,local_best_estimate_locs[kk][0,0]],1))
                   
               
        return(local_best_cost, local_best_estimate, local_best_estimate_locs, z_k_best_arr)
    #--------------------------------------------------------------------------






    #--------------------------------------------------------------------------
    # Modified cost function for multiple iterations - 9 sigma^2 noise weight always present
    # Based on histogram weigths, without interpolation, assumes 0 step size not 
    # included while computing histogram. Step size of 0 given a penalty of 0.    
    #--------------------------------------------------------------------------
    def step_cost_parsimonious(self, x_candidates, x_measured, noise_std, cost_prev_step, iter_no, hist_data):
                        
        energy_cost = np.power((x_candidates - x_measured),2) 
        x_prev_step = x_candidates.T
        
        if(iter_no == 0):
            penalty_weight = 9*(noise_std**2)                     
            
        else:
            step_sizes = x_candidates-x_prev_step
            step_prob = self.SDAHistObj.SDA_hist(step_sizes, hist_data[0], hist_data[1])
            penalty_weight = (-2*(noise_std**2))*np.log(step_prob) + 9*(noise_std**2)  
            
            
        penalty_possibilities = np.multiply((x_prev_step != x_candidates),penalty_weight)  
        penalty_possibilities[np.isnan(penalty_possibilities)] = 0         
        cost_possibilities = energy_cost + penalty_possibilities + cost_prev_step
        
        local_best_cost = np.min(cost_possibilities,0).T
        local_best_estimate_locs = np.argmin(cost_possibilities,0).T
        local_best_estimate = x_candidates[local_best_estimate_locs]
        local_best_estimate = local_best_estimate.reshape((np.size(local_best_estimate,0),1))
        
        
        return(local_best_cost, local_best_estimate, local_best_estimate_locs)
    #-------------------------------------------------------------------------- 
    
    
# =============================================================================
