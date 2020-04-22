# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 14:21:28 2019

@author: siva_
"""

# =============================================================================
# Description
# =============================================================================
# Code for performance evaluation functions
# =============================================================================



# =============================================================================
# Imports
# ============================================================================= 
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sc
from scipy import stats


import sys
sys.path.append('../../')


# my toolkits
from APIs.SDACore.CostFunctions import CostFunctions
from APIs.SDACore.SDAHistogram import SDAHistogram
from APIs.SDACore.SDA import SDA
from APIs.SDACore.SDADynamics import SDADynamics
from APIs.Validation.ErrorMetrics import ErrorMetrics
from APIs.Validation.SampleSignals import SampleSignals
from APIs.AlternateTools.AlternateMethods import AlternateMethods
from APIs.Utils.Utils import Utils


# =============================================================================
#%%
# =============================================================================
# Class containing utility functions
# =============================================================================
class PerformanceEvaluation():
    
    #--------------------------------------------------------------------------        
    # Init
    #--------------------------------------------------------------------------  
    def __init__(self):         
        self.sigObj = SampleSignals()
        self.errMetObj = ErrorMetrics()
        self.UtilObj = Utils()
        self.SDynObj = SDADynamics()
        
        self.SDAObj = SDA()
        self.altMethodObj = AlternateMethods()        
    #--------------------------------------------------------------------------        


    #--------------------------------------------------------------------------        
    # Plot metric comparisons
    #--------------------------------------------------------------------------  
    def plot_metric_comparisons(self, error_stat_df, method_list, metric_name, param_interest, fmt_list = None, labels = None, save_fig = 1, ebars = 0): 
        
        if(labels == None):
            labels = {}
            labels['title'] = ''
            labels['xlabel'] = ''
            labels['ylabel'] = ''
        
        # Plot error means
        plt.figure()
        
        for ii in range(0,len(method_list),1):
            
            method_str = method_list[ii]
            
            if(fmt_list == None):
                fmt = ''
            else:
                fmt = fmt_list[ii]
            
            if(ebars == 0):
                plt.plot(error_stat_df.loc[:,param_interest], error_stat_df.loc[:,'mean_'+method_str], fmt = fmt)
            
            elif(ebars == 1):
                ebar_low = error_stat_df.loc[:,'CI_low_'+method_str]
                ebar_high = error_stat_df.loc[:,'CI_high_'+method_str]
                #plt.fill_between(error_stat_df.loc[:,'sigma'], ebar_low, ebar_high)                
                plt.errorbar(error_stat_df.loc[:,param_interest],error_stat_df.loc[:,'mean_'+method_str], yerr = (ebar_high-ebar_low)/2, fmt = fmt)
        
        
        
        # Adding a custom dictionary for the legend
        legend_dict = {}
        legend_dict['SDA'] = 'SDA'
        legend_dict['SDA_OR'] = 'STEP FIT'
        legend_dict['HMM'] = 'VS-HMM'
        legend_dict['KER'] = '$\chi ^2$'
        legend_dict['TVR'] = 'TVR'
        
        legend_str = []        
        for jj in range(len(method_list)):
            legend_str.append(legend_dict[method_list[jj]])
            
            
        #legend_str = method_list
        
        plt.grid()
        plt.title(labels['title'])
        plt.xlabel(labels['xlabel'])
        plt.ylabel(labels['ylabel'])
        plt.legend(legend_str)
        
        
        if(save_fig == 1):
            
            plt.savefig(metric_name +'_'+ labels['title'] +'_' + '_Comparisons.png')
            plt.savefig(metric_name +'_'+ labels['title'] +'_' + '_Comparisons.svg')
    
    
        return()
    #--------------------------------------------------------------------------        


    #--------------------------------------------------------------------------        
    # Extract Performance Statistics
    #-------------------------------------------------------------------------- 
    def get_perf_stats(self, error_df, method_list, param_interest, conf_val = 0.95):
        
#        method_list = ['SDA', 'SDA_OR', 'TVR', 'HMM']
        
        
#        error_df = AEE_array
        
        error_stat_df = pd.DataFrame([])
        error_stat_df[param_interest] = []
        
        for ii in range(0,len(method_list),1):
            
            method_str = method_list[ii]
            
            error_stat_df['mean_'+method_str] = []
            error_stat_df['CI_low_'+method_str] = []
            error_stat_df['CI_high_'+method_str] = []
            
            
        # Extract error statistics
        
        sigma_arr = error_df.loc[:,param_interest]
        sigma_unique = np.unique(sigma_arr)
            
        for ii in range(0,np.size(sigma_unique),1):
            
            row_sel = error_df.loc[:, param_interest] == sigma_unique[ii]    
            error_stat_df.loc[ii, param_interest] = sigma_unique[ii]
            
            for jj in range(0,len(method_list),1):
                
                method_str = method_list[jj]
                
                error_stat_df.loc[ii,'mean_'+method_str] = np.mean(error_df.loc[row_sel, method_str])
                conf_int = self.get_conf_interval(error_df.loc[row_sel, method_str],confidence = conf_val)
                error_stat_df.loc[ii,'CI_low_'+method_str] = conf_int[0]
                error_stat_df.loc[ii,'CI_high_'+method_str] = conf_int[1]            
    
        return(error_stat_df)
    #-------------------------------------------------------------------------- 



    #--------------------------------------------------------------------------
    # Get Student's t confidence interval
    #--------------------------------------------------------------------------
    def get_conf_interval(self, data, confidence=0.95):
        data = np.array(data)
        n = len(data)
        data_mean = np.mean(data)
        data_sem = sc.stats.sem(data)
        #h = data_sem * sc.stats.t.ppf((1 + confidence)/2., n-1)
        conf_int = sc.stats.t.interval(confidence, n-1, loc = data_mean, scale = data_sem)
        return(conf_int)
    #--------------------------------------------------------------------------
    
    
    #--------------------------------------------------------------------------
    # Get shmoo data
    #--------------------------------------------------------------------------
    def get_shmoo_data(self, step_params, no_runs, var_to_change, var_range, method_list, master_seed = 9, tol = 0, plot_diagnostics = 1, save_diagnostics = 1, save_rawdata = 1):
                
        d = step_params['no_of_steps']
        N = step_params['tot_no_samples']
        L = step_params['step_size']
        sigma = step_params['sigma']
        N_d_ratio = step_params['samples_per_step']
        t_avg = step_params['average_step_duration'] 
        
        sensor_dynamics = step_params['sensor_dynamics_ON'] 
        num_sensor = step_params['sensor_num']
        den_sensor = step_params['sensor_den']
        #Ts = step_params['sampling_time']
        
        np.random.seed(master_seed)
        seed_arr = np.random.randint(2**31, size = no_runs)
                     
        var_arr = np.arange(var_range[0],var_range[1],var_range[2])
                
        path_diagnostics = '../../../Plots/Diagnostics/'
        path_metrics = '../../../Metrics/'
        
        val_cols = [var_to_change,'Run','Seed']
        AEE_array = pd.DataFrame(columns = val_cols)
        RMSE_array = pd.DataFrame(columns = val_cols)
        KS_array = pd.DataFrame(columns = val_cols)
                
        #-- Data Collection DataFrame --#
        data_df = pd.DataFrame(columns = val_cols)
        data_df['x_true'] = []
        data_df['x_true'] = data_df['x_true'].astype(object)
        
        data_df['x_measured'] = []
        data_df['x_measured'] = data_df['x_measured'].astype(object)
                
        data_df['SDA'] = []
        data_df['SDA'] = data_df['SDA'].astype(object)
        
        data_df['SDA_OR'] = []
        data_df['SDA_OR'] = data_df['SDA_OR'].astype(object)
        
        data_df['TVR'] = []
        data_df['TVR'] = data_df['TVR'].astype(object)
        
        data_df['HMM'] = []
        data_df['HMM'] = data_df['HMM'].astype(object)
        
        data_df['KER'] = []
        data_df['KER'] = data_df['KER'].astype(object)
        #-------------------------------#        
        
        exp_cnt = 0
        for var in var_arr:     
            run = 0
            
            for seed in seed_arr:
                
                AEE_array = AEE_array.append({var_to_change: var, 'Run': run, 'Seed': seed},ignore_index = True)
                RMSE_array = RMSE_array.append({var_to_change: var, 'Run': run, 'Seed': seed},ignore_index = True)
                KS_array = KS_array.append({var_to_change: var, 'Run': run, 'Seed': seed},ignore_index = True)
                
                data_df.loc[exp_cnt, var_to_change] =  var
                data_df.loc[exp_cnt, 'Run'] = run
                data_df.loc[exp_cnt, 'Seed'] = seed
                
                if(var_to_change == 'sigma'):
                    sigma = var*L
                    
                elif(var_to_change == 'step_size'):
                    L = var
                    
                elif(var_to_change == 'tot_no_samples'):
                    N = var
                    
                elif(var_to_change == 'no_of_steps'):
                    # Keeping the number of samples per step constant
                    d = var
                    N = N_d_ratio*d
                    
                elif(var_to_change == 'N_d_ratio'):
                    N_d_ratio = var
                    N = N_d_ratio*d
                    
                elif(var_to_change == 'Sensor_Time_Constant'):
                    num_sensor = [1]
                    den_sensor = [var, 1]
                    
                    sensor_dynamics = 1

                # Computing the correct sampling time for each specific example
                dt = float(t_avg/N_d_ratio)

                # Generate step signal                    
                x = self.sigObj.get_step_sample(step_size = L, no_steps = d, no_samples = N, seed = seed)
                noise = self.sigObj.get_gaussian_noise(std = sigma, no_samples = N, seed = seed)     
                # Optional sensor dynamics
                if(sensor_dynamics==1):                    
                    [T, x_sensor] = self.sigObj.sensor_dynamics(num_sensor, den_sensor, dt, x)
                    #[T, x_sensor] = self.sigObj.sensor_dynamics(num_sensor, den_sensor, dt, x+noise)
                    x_sensor = np.expand_dims(x_sensor, axis = 0)
                else:
                    T = np.arange(0,N)*dt
                    x_sensor = x
                               
                y = x_sensor + noise  # Noisy measurements
                #y = x_sensor   
                
                if(save_rawdata == 1): 
                    data_df.loc[exp_cnt,'x_true'] = [x]
                    data_df.loc[exp_cnt,'x_measured'] = [y]                 
                
                for method_cnt in range(0,len(method_list),1):
                    
                    method_str = method_list[method_cnt]
                    
                    if(method_str == 'SDA'):
                        
                        iterations = 100 # No. of SDA iterations
                        resolution = 0.05 # Resolution of the SDA   
                        hist_tol = 1e-7 # Tolerance of histogram differences (KS distance)
                        
                        #[step_estimate, est_array, hist_array]= self.SDAObj.SDA_standard(sigma, y, resolution = resolution,  iterations = iterations, make_plots = 0, diagnostics = 1)
                        
                        if(sensor_dynamics == 0):
                            sensor_sysd = None
                        else:
                            sensor_sysd = self.SDynObj.get_DT_sys(num_sensor, den_sensor, dt)
                            
                        [step_estimate, est_array, hist_array]= self.SDAObj.SDA_dynamics(sigma, y, resolution = resolution, hist_tol = hist_tol,  iterations = iterations, make_plots = 0, diagnostics = 1, sysd = sensor_sysd)
                        
                    elif(method_str == 'SDA_OR'):
                        
                        if(sensor_dynamics == 0):
                            # No sensor dynamics
                            step_estimate = self.altMethodObj.SDA_original(y)
                        else:
                            # Incorporating sensor dynamics (assumes first order dyanmics )
                            step_estimate = self.altMethodObj.SDA_original(y, Ts = dt, tau = den_sensor[1])
                            
                        step_estimate = step_estimate.T
                        
                    elif(method_str == 'TVR'):
                        
                        step_estimate = self.altMethodObj.total_variation_regularization(sigma, y)
                        step_estimate = np.array(step_estimate)
                        
                    elif(method_str == 'HMM'):
                        step_estimate = self.altMethodObj.HMMsForMotors(y)
                        
                        est_ptp = np.ptp(step_estimate)
                        x_ptp = np.ptp(x)
                        scaling = est_ptp/x_ptp                
                        
                        step_estimate = step_estimate/scaling
                        
                    elif(method_str == 'KER'):
                        step_estimate = self.altMethodObj.Kerssemakers(y)
                        step_estimate = np.array(step_estimate)
                     
                    # Save error metrics
                    AEE = self.errMetObj.get_err_array(step_estimate.T, x, self.errMetObj.get_AEE, tol = tol)
                    RMSE = self.errMetObj.get_err_array(step_estimate.T, x, self.errMetObj.get_RMSE, tol = tol)        
                    [ks_stat, ks_pval] = self.errMetObj.get_ks_2samp(step_estimate, x)
                                        
                    AEE_array.loc[exp_cnt, method_str] = AEE[0]
                    RMSE_array.loc[exp_cnt, method_str] = RMSE[0]                
                    KS_array.loc[exp_cnt, method_str] = ks_stat
                           
                    if(save_rawdata == 1):                                                     
                        data_df.loc[exp_cnt, method_str] = [step_estimate] 
                    
                    # Save diagnostic plots
                    if(plot_diagnostics == 1):
                        plt.figure()
                        plt.plot(T,y.T)
                        plt.plot(T,x.T)
                        plt.plot(T,x_sensor.T)
                        plt.plot(T,step_estimate)
                        plt.title(method_str)
                        
                        if(save_diagnostics == 1):
                            fig_path = path_diagnostics + method_str +'/'                                
                            fig_name = fig_path+'StepFit_'+var_to_change+'_'+str(var)+'_Run_' + str(run) 
                            self.UtilObj.folder_check(fig_path)
                            plt.savefig(fig_name + '.png')
                            plt.savefig(fig_name + '.svg')
                            
                        plt.close()
                
                run += 1
                exp_cnt += 1
                
        # Pickle (save) the metric arrays
        metric_path = path_metrics + var_to_change + '/'
        self.UtilObj.folder_check(metric_path)
        AEE_array.to_pickle(metric_path + 'AEE_array.pkl')
        RMSE_array.to_pickle(metric_path + 'RMSE_array.pkl')
        KS_array.to_pickle(metric_path + 'KS_array.pkl')     
        
        data_df.to_pickle(metric_path+'Data.pkl')               

        return(AEE_array, RMSE_array, KS_array, data_df) #, [x, x_sensor], step_estimate)
    #--------------------------------------------------------------------------    
    
    
    #--------------------------------------------------------------------------
    # Generate error metrics from collected & saved simulated data
    #--------------------------------------------------------------------------
    def get_error_df(self, data_df, method_list, error_type = 'RMSE', normalize = 'None', tol = 0):
    
        if(error_type == 'RMSE'):
            err_func = self.errMetObj.get_RMSE
        elif(error_type == 'AEE'):
            err_func = self.errMetObj.get_AEE
        elif(error_type == 'KS_Distance'):
            err_func = self.errMetObj.get_ks_2samp
        elif(error_type == 'Classification_metric'):
            err_func = self.get_rho
        elif(error_type == 'F1_score'):
            err_func = self.get_f1
            
        param_interest = data_df.columns[0]
        error_df = pd.DataFrame(columns = [param_interest,'Run', 'Seed']+method_list)
        nb_rows = data_df.shape[0]
        
        for ii in range(0,nb_rows,1):
#            print(ii)
            error_df.loc[ii, param_interest] = data_df.loc[ii, param_interest] 
            error_df.loc[ii, 'Run'] = data_df.loc[ii, 'Run']
            error_df.loc[ii, 'Seed'] = data_df.loc[ii, 'Seed']
            
            x = data_df.loc[ii, 'x_true'][0]
                
            for method in method_list:
                
#                print(method)
                step_estimate = data_df.loc[ii, method][0]   
                if(error_type == 'Classification_metric'):
                    err_val = self.get_rho(step_estimate.T, x, tol=tol, normalize = normalize)
                elif(error_type == 'F1_score'):
                    err_val = self.get_f1(step_estimate.T, x, tol=tol, normalize = normalize)
                else:
                    err_val = self.errMetObj.get_err_array(step_estimate.T, x, err_func, tol = tol,\
                                                  normalize = normalize)        
                error_df.loc[ii, method] = err_val[0]
        
        
        
        if(param_interest == 'sigma'):
            labels = {}
            labels['title'] = 'Performance_vs_Noise_Mangitude'
            labels['xlabel'] = 'Standard Deviation of Noise'
            
        elif(param_interest == 'no_of_steps'):
            labels = {}
            labels['title'] = 'Performance_vs_No_of_Steps'
            labels['xlabel'] = 'No. of steps'
            
        elif(param_interest == 'N_d_ratio'):
            labels = {}
            labels['title'] = 'Performance_vs_No_of_Samples_per_Step'
            labels['xlabel'] = 'No. of samples per step'
        
        elif(param_interest == 'Sensor_Time_Constant'):
            labels = {}
            labels['title'] = 'Performance_vs_Sensor_Time_Constant'
            labels['xlabel'] = 'Sensor Time Constant'
        

        return(error_df, labels, param_interest)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # Obtain classification metrics by approximation to a tolerance
    #--------------------------------------------------------------------------
    def get_classification_metrics(self, x_est, x_true, tol_time, tol_mag):
        # Get true step locations
        x_t_diff = np.diff(x_true)
        x_t_locs = (x_t_diff != 0)
        
        # Total positive and total negative
        P = np.sum(x_t_locs)
        N = np.size(x_t_locs) - P
        # Generate a tolerance band around the true steps
        #tol_time = 5 # Number of samples of tolerance on either side
        #tol_mag = 3 # Percentage magnitude tolerance in either direction
        
        
        # Determine the true location of steps and the tolerance band
        x_t_locs_tol = np.zeros(np.shape(x_t_locs))
        
        len_x = np.size(x_t_locs)
        for ii in range(0,len_x,1):
        
            min_t = np.max([0, ii-tol_time])
            max_t = np.min([len_x, ii+tol_time])
        
            if(np.sum(x_t_locs[0,min_t:max_t+1])>0):
                
                x_t_locs_tol[0,ii] = 1
        
        
        #plt.plot(x_t_locs_tol.T)
        #plt.plot(x_t_locs.T)
        
        # Determine the estimated location of steps
        x_e_diff = np.diff(x_est.T)
        x_e_locs = (x_e_diff != 0)
        #plt.plot(x_e_locs.T)
        
        x_t_locs_array = []
        for ii in range(0,len_x,1):    
            if(x_t_locs[0,ii] == 1):        
                x_t_locs_array.append(ii)
        
        x_e_locs_array = []
        for ii in range(0,len_x,1):   
            if(x_e_locs[0,ii] == 1):        
                x_e_locs_array.append(ii)
                
        TP = 0
        TN = 0
        FP = 0
        FN = 0
        
        for ii in range(0,np.size(x_t_locs_array),1):    
            loc_t = x_t_locs_array[ii]
            step_est_exists = 0  
            step_est_is_correct = 0
            for jj in range(0,np.size(x_e_locs_array),1):        
                loc_e = x_e_locs_array[jj]        
                if((loc_e < loc_t + tol_time) or (loc_e > loc_t - tol_time)):            
                    step_est_exists = 1    
                    step_est_is_correct = 0
                    # Step exists, now check if its magnitude is correct 
                    if(np.abs(100*(x_e_diff[0,loc_e]-x_t_diff[0,loc_t])/x_t_diff[0,loc_t])<tol_mag):
                        step_est_is_correct = 1
                        x_e_locs_array.pop(jj)
                        break  
                        
            if(step_est_is_correct == 1):
                TP = TP + 1
            else:
                FN = FN + 1
        
        FP = np.size(x_e_locs_array)
        TN = len_x - TP - FP - FN
        return(TP, FP, TN, FN, P, N)
    #--------------------------------------------------------------------------


    #--------------------------------------------------------------------------
    # Obtain TPR, FPR, and ROC metric by approximation to a tolerance
    #--------------------------------------------------------------------------
    def get_cls_performance(self, TP, FP, TN, FN, P, N):
            
        # True Positive Rate | Sensitivity | Recall 
        TPR = (TP/P)*100
        
        # True Negative Rate | Specificity | Selectivity
        TNR = (TN/N)*100
        
        # False Positive Rate
        FPR = (FP/N)*100
                        
        # Get performance metric rho
        rho = (np.sqrt(2)/100)*(np.sin(np.arctan(TPR/FPR) - (np.pi/4)) )*(np.sqrt(TPR**2+FPR**2))
            
        return(rho)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Obtain ROC metric by approximation to a tolerance
    #--------------------------------------------------------------------------
    def get_rho(self, x_est, x_true, tol = [5,5], normalize = 'None'):
    
        [tol_time, tol_mag] = tol
        
        [TP, FP, TN, FN, P, N] = self.get_classification_metrics(x_est.T, x_true, tol_time, tol_mag)
        
        rho = self.get_cls_performance(TP, FP, TN, FN, P, N)
        
        return([rho])
        
    #--------------------------------------------------------------------------

        
    
    #--------------------------------------------------------------------------
    # Find F1 score
    #--------------------------------------------------------------------------
    def get_f1(self, x_est, x_true, tol = [5,5], normalize = 'None'):
    
        [tol_time, tol_mag] = tol
        
        [TP, FP, TN, FN, P, N] = self.get_classification_metrics(x_est.T, x_true, tol_time, tol_mag)
        
        f1 = (2*TP)/(2*TP + FP + FN) 
                
        return([f1])
    #--------------------------------------------------------------------------
    
    
# =============================================================================

    
    
    
    
    


