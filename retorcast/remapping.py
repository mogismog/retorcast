#!/usr/bin/env python
# Modules for severe_reforecast.py
# Written by F.Alvarez, 1/2014

"""
 Python module for remapping probabilities from larger -> smaller ROI.

"""

import numpy as np
from netCDF4 import Dataset,date2index
from fortran_routines import cdf_loop,pct_loop
from utils import find_nearest_idx
from verification import remap_verification_data

def remap_probs(probabilities,leadtime,old_forecasts_dir,smaller_radius,**kwargs):
    """
    remap_probs(probabilities,leadtime,reforecast_dir,tordata_dir,smaller_radius)
    
    Function that remaps forecast probabilities from larger -> smaller ROI
    
    Required Arguments:
        probabilities - some 2-D numpy array of forecast probabilities [0,100]
        leadtime - forecast leadtime
        reforecast_dir - path to remapped probs data
        tordata_dir - path to verification data
        smaller_radius - a lower radius to which to remap
    Optional Arguments:
        train_fname - basically the var combo we used to generate older fcsts (default = "cape06shearcin")
        larger_radius - Larger ROI (default = 240)
    """
    train_fname = kwargs.get('train_fname',"1mocdf_cape06shearcin")
    larger_radius = kwargs.get('larger_radius',240)

    
    # --- First, get remapped information...
    nc_fname = '{}rfcst2_{}_remapped_probs_1985-2011_{}km_to_{}km_day{}.nc'.format(old_forecasts_dir,train_fname,larger_radius,smaller_radius,leadtime)
    trainData = Dataset(nc_fname,'r')
    old_probs = trainData.variables['remapped_probabilities'][:]
    prob_levs = trainData.variables['probability_levels'][:]

    # --- Now, remap the data...
    new_probabilities = np.rint(probabilities)
    
    for i in prob_levs[:]:
        new_probabilities[((new_probabilities < i+1)&(new_probabilities>=i-1))] = old_probs[i]
    
    return new_probabilities    

def prob_members(forecastDate,forecastData,probabilities,leadtime,reforecast_dir,cdf_dir,**kwargs):
    """
    prob_members(forecastDate,forecastData,probabilities,leadtime,reforecast_dir,cdf_dir,**kwargs)
    
    Function that remaps forecast probabilities from larger -> smaller ROI
    
    Required Arguments:
        forecastDate - some datetime object, the date of the forecast initialization
        probabilities - some 2-D numpy array of forecast probabilities [0,100]
        leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)
        reforecast_dir - path to remapped probs data
        cdf_dir - path to CDF data
    Optional Arguments:
        train_fname - basically the var combo we used to generate older fcsts (default = "cape06shearcin")
        mems_to_keep - members to keep based on a list of percentiles (default = lots...), probably shouldn't change...
        members - numpy array which has all number of members to keep for each individual forecast (default = lots...)
    """
    train_fname = kwargs.get('train_fname',"1mocdf_corr_cape06shearcin_avg")
    mems_to_keep = kwargs.get('mems_to_keep',np.array([200,200,175,150,150,150,125,125,100,75,30]))
    #mems_to_keep = kwargs.get('mems_to_keep',np.array([200,200,200,175,175,175,175,150,100,75,30]))
    members = kwargs.get('members',np.array([10,15,20,25,30,35,40,45,50,75,100,125,150,175,200,500,1000]))
    pct_min_list = [0,.2,.3,.4,.5,.6,.7,.8,.9,.99]
    pct_max_list = [.19,.29,.39,.49,.59,.69,.79,.89,.98,1.]
    
    # --- Get CDF percentiles
    infile = '{}/rfcst2_{}_{}_day{}.nc'.format(cdf_dir,forecastDate.strftime("%b"),train_fname,leadtime)
    nc = Dataset(infile)
    quantiles = nc.variables['stp_quantiles'][:,:,:]
    pct = nc.variables['pct'][:]      

    idx2 = find_nearest_idx(quantiles[:,:,:],forecastData[:,:])
    per = pct_loop(idx2,pct,pct.shape[0],idx2.shape[0],idx2.shape[1])
        
    nc.close()
    
    new_probabilities = np.zeros((probabilities.shape[-2:]))
    for i in xrange(len(pct_min_list)):
        per_mask = np.array((per[:,:] >= pct_min_list[i]) & (per[:,:] <= pct_max_list[i]))
        mem_idx = np.where(members == mems_to_keep[i])[0]
        new_probabilities[per_mask] = probabilities[mem_idx,per_mask] 
    
    return new_probabilities, per
    
def fcst_quants(forecastDate,forecastData,leadtime,reforecast_dir,cdf_dir,**kwargs):
    """
    prob_members(forecastDate,forecastData,probabilities,leadtime,reforecast_dir,cdf_dir,**kwargs)
    
    Function that remaps forecast probabilities from larger -> smaller ROI
    
    Required Arguments:
        forecastDate - some datetime object, the date of the forecast initialization
        probabilities - some 2-D numpy array of forecast probabilities [0,100]
        leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)
        reforecast_dir - path to remapped probs data
        cdf_dir - path to CDF data
    Optional Arguments:
        train_fname - basically the var combo we used to generate older fcsts (default = "cape06shearcin")
        mems_to_keep - members to keep based on a list of percentiles (default = lots...), probably shouldn't change...
        members - numpy array which has all number of members to keep for each individual forecast (default = lots...)
    """
    train_fname = kwargs.get('train_fname',"1mocdf_corr_cape06shearcin_avg")
    
    # --- Get CDF percentiles
    infile = '{}/rfcst2_{}_{}_day{}.nc'.format(cdf_dir,forecastDate.strftime("%b"),train_fname,leadtime)
    nc = Dataset(infile)
    quantiles = nc.variables['stp_quantiles'][:,:,:]
    pct = nc.variables['pct'][:]      

    idx2 = find_nearest_idx(quantiles[:,:,:],forecastData[:,:])
    per = pct_loop(idx2,pct,pct.shape[0],idx2.shape[0],idx2.shape[1])
    print per.max(),per.min()
    nc.close()
    
    return per    