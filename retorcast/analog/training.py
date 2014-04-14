#!/usr/bin/env python
# Modules for severe_reforecast.py
# Written by F.Alvarez, 1/2014

"""
 Python module for extracting data for verification.

"""

import numpy as np
from netCDF4 import Dataset,date2index

def get_training_data(analog_dates,leadtime,train_fname,predictor_name,reforecast_dir):
    """
    Function that collects training data for statistical
    post-processing.
    
    analog_dates - list of datetime objects of potential analogous dates
    leadtime - forecast leadtime
    reforecast_dir - path to verification data
    """
    
    # --- First, get analog dates
    nc_fname = '{}refcst2_{}_day{}.nc'.format(reforecast_dir,train_fname,leadtime)
    trainData = Dataset(nc_fname,'r')
    train_idxs = date2index(analog_dates,trainData.variables['time'])
    trainData = np.asfortranarray(trainData.variables[predictor_name][train_idxs,:,:])

    
    return trainData    

def get_pct(analog_dates,leadtime,train_fname,predictor_name,reforecast_dir):
    """
    Function that collects training data for statistical
    post-processing.
    
    analog_dates - list of datetime objects of potential analogous dates
    leadtime - forecast leadtime
    reforecast_dir - path to verification data
    """
    
    # --- First, get analog dates
    
    nc_fname = '{}rfcst2_{}_pct_day{}.nc'.format(reforecast_dir,train_fname,leadtime)
    print nc_fname
    trainData = Dataset(nc_fname,'r')
    train_idxs = date2index(analog_dates,trainData.variables['time'])
    trainData = np.asfortranarray(trainData.variables['percentiles'][train_idxs,:,:])
    
    return trainData    