# -*- coding: utf-8 -*-

# =============================================================================
# Description
# =============================================================================

 # Code for fit evaluation and fit errors

# =============================================================================
 
 
 
# =============================================================================
# Imports
# =============================================================================
import numpy as np
from scipy.stats import ks_2samp
import copy

import sys
sys.path.append('../../')
# =============================================================================


# =============================================================================
# Class containing error metric functions
# =============================================================================
class ErrorMetrics():

    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------        
    def __init__(self):
        a = 1
   
    #--------------------------------------------------------------------------    


    #--------------------------------------------------------------------------
    # RMSE: Root Mean Square Error (use column vectors)
    #--------------------------------------------------------------------------
    def get_RMSE(self, x_estimate, x_actual, tol = 0, normalize = 'None'):
                       
        err = x_estimate - x_actual
        
        # absolute values of err less than tol considered to be 0
        err_tol_chk = ( np.abs(err) <= tol )
        
        err_tol = copy.deepcopy(err)
        err_tol[err_tol_chk] = 0.0
        
        n = np.size(err,0)
        rmse = np.linalg.norm(err_tol,2)
        #rmse = np.sqrt(rmse/n) # Error in code
        rmse = rmse/np.sqrt(n)
        
        if(normalize == 'None'):
            rmse = rmse
        elif(normalize == 'range'):
            x_ptp = np.ptp(x_actual)
            rmse = rmse/np.abs(x_ptp)
            
        elif(normalize == 'mean_estimate'):
            mu = np.mean(x_actual)
            err_mu = x_actual - mu
            m_rmse = (np.linalg.norm(err_mu, 2))/np.sqrt(n)
            rmse = rmse/m_rmse
            
        
        return(rmse)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # AEE: Average Euclidean Error (use column vectors)
    #--------------------------------------------------------------------------
    def get_AEE(self, x_estimate, x_actual, tol = 0, normalize = 'None'):
                       
        err = x_estimate - x_actual
        
        # absolute values of err less than tol considered to be 0
        err_tol_chk = ( np.abs(err) <= tol )
        
        err_tol = copy.deepcopy(err)
        err_tol[err_tol_chk] = 0.0        
        
        n = np.size(err,0)
        aee = np.absolute(err_tol)
        aee = np.sum(aee)/n
        
        return(aee)
    #--------------------------------------------------------------------------
    
    
    
    #--------------------------------------------------------------------------
    # Error function over an array
    #--------------------------------------------------------------------------
    def get_err_array(self, est_array, x_actual, err_func, tol = 0, normalize = 'None'):
                       
        arr_len = np.size(est_array,0)
        err_array = []
        for ii in range(0,arr_len,1):
            err = err_func(est_array[ii],x_actual, tol = tol, normalize = normalize)
            err_array.append(err)
        
        return(err_array)
    #--------------------------------------------------------------------------    


    #--------------------------------------------------------------------------
    # Generate 2 sample ks metric
    #--------------------------------------------------------------------------
    def get_ks_2samp(self, est_array, x_actual):
                       
        a = np.array(est_array)
        a = np.squeeze(a)
        a = np.diff(a)
        a = a[a!=0]
        
        b = np.array(x_actual)
        b = np.squeeze(b)
        b = np.diff(b)
        b = b[b!=0]
        
        [ks_stat, ks_pval] = ks_2samp(a,b)
        
        return([ks_stat, ks_pval])
    #-------------------------------------------------------------------------- 


# =============================================================================