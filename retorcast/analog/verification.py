#!/usr/bin/env python
# Modules for severe_reforecast.py
# Written by F.Alvarez, 1/2014

"""
 Python module for extracting data for verification.

"""

import numpy as np
from netCDF4 import Dataset,date2index

def get_verification_data(analog_dates,leadtime,tordata_dir,radius):
    """
    Function that collects verification data for statistica
    post-processing.
    
    analog_dates - list of datetime objects of potential analogous dates
    leadtime - forecast leadtime
    tordata_dir - path to verification data
    radius - some ROI.
    """
    
    # --- First, get analog dates
    tor_fname = "{}/tor_day_data_1985-2011_day{}.nc".format(tordata_dir,leadtime)
    tor_nc = Dataset("{}".format(tor_fname),'r')
    tor_idxs = date2index(analog_dates,tor_nc.variables['time'])
    roi_idx = np.where(tor_nc.variables['roi'][:] == radius)[0] 
    tornado_array = np.asfortranarray(tor_nc.variables['Tornado_reports'][tor_idxs,roi_idx,:,:])
    tor_mask = (tornado_array >= 1)
    tornado_array[tor_mask] = 1
    tor_mask = (tornado_array < 1)
    tornado_array[tor_mask] = 0
    
    return tornado_array


def get_climo_data(analog_dates,leadtime,tordata_dir,radius):
    """
    Function that collects verification data for statistica
    post-processing.
    
    analog_dates - list of datetime objects of potential analogous dates
    leadtime - forecast leadtime
    tordata_dir - path to verification data
    radius - some ROI.
    """
    
    # --- First, get analog dates
    tor_fname = "{}/tor_day_data_1985-2011_day{}.nc".format(tordata_dir,leadtime)
    tor_nc = Dataset("{}".format(tor_fname),'r')
    tor_idxs = date2index(analog_dates,tor_nc.variables['time'])
    roi_idx = np.where(tor_nc.variables['roi'][:] == radius)[0] 
    tornado_array = np.asfortranarray(tor_nc.variables['Tornado_reports'][tor_idxs,roi_idx,:,:])
    tor_mask = (tornado_array >= 1)
    tornado_array[tor_mask] = 1
    tor_mask = (tornado_array < 1)
    tornado_array[tor_mask] = 0
    
    
    return np.mean(tornado_array,axis=0)*100.

def remap_verification_data(leadtime,tordata_dir,radius):
    """
    Function that collects verification data for statistical
    post-processing.
    
    leadtime - forecast leadtime
    tordata_dir - path to verification data
    radius - some ROI.
    """
    
    # --- First, get analog dates
    tor_fname = "{}/tor_day_data_1985-2011_day{}.nc".format(tordata_dir,leadtime)
    tor_nc = Dataset("{}".format(tor_fname),'r')
    roi_idx = np.where(tor_nc.variables['roi'][:] == radius)[0] 
    tornado_array = np.asfortranarray(tor_nc.variables['Tornado_reports'][:,roi_idx,:,:])
    tor_mask = (tornado_array >= 1)
    tornado_array[tor_mask] = 1
    tor_mask = (tornado_array < 1)
    tornado_array[tor_mask] = 0
    
    return tornado_array
   