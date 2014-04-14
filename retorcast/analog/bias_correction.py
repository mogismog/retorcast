#!/usr/bin/env python
# Modules for severe_reforecast.py
# Written by F.Alvarez, 1/2014

"""
 Python module for extracting data for verification.

"""

import numpy as np
from netCDF4 import Dataset,date2index
from retorcast.utils import find_nearest_idx
from fortran_routines import cdf_loop

def correct_bias(forecastDate,field_name,forecastField,cdf_dir,cdf_name,begin_fcst_hour,end_fcst_hour,maxlat,minlat,maxlon,minlon):
    """
        Function that corrects bias in forecast fields based on past reforecast data.
    """    
    
    # --- Checking the forecast hour...
    if end_fcst_hour <= 180:
        
        # --- Get domain lat/lon indices...
        lts = np.arange(minlat,maxlat+1,1)
        lns = np.arange(minlon,maxlon+1,1)
        
        forecast_cdf_nc = Dataset('{0}/rfcst2_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        forecast_fhrs = forecast_cdf_nc.variables['fhour'][:]
        forecast_quantiles = forecast_cdf_nc.variables["{}_quantiles".format(field_name)][:,:,:,:]        
        
        analysis_cdf_nc = Dataset('{0}/cfsr_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        analysis_fhrs = analysis_cdf_nc.variables['fhour'][:]
        analysis_quantiles = analysis_cdf_nc.variables["{0}_quantiles".format(field_name)][:,:,:,:]        
        
        # --- array of forecast hours
        hrs = np.arange(begin_fcst_hour,end_fcst_hour+1,6)
        
        cfsr_counter = np.where(analysis_fhrs[:] == 12)[0]
        field = np.ones(forecastField.shape[:])*-9999.9
        
        for fhr_idx in xrange(0,hrs.shape[0]):
            if cfsr_counter == (np.where(analysis_fhrs[:] == 12)[0])+2: # --- Past 24 hour mark
                cfsr_counter = np.where(analysis_fhrs[:] == 0)[0]
            rfcst2_fhr_idx = np.where(forecast_fhrs[:] == hrs[fhr_idx])[0] 
            dummy_rfcst_quants = forecast_quantiles[rfcst2_fhr_idx[0],:,:,:]
            dummy_rfcst = np.absolute(forecastField[fhr_idx,:,:])
            idx = find_nearest_idx(dummy_rfcst_quants,dummy_rfcst)
            field[fhr_idx,:,:] = cdf_loop(analysis_quantiles[cfsr_counter,:,:,:],idx,analysis_quantiles.shape[1],lts.shape[0],lns.shape[0])
            cfsr_counter += 1

        # --- Close up shop
        forecast_cdf_nc.close();analysis_cdf_nc.close();
        
    elif begin_fcst_hour > 180:
        
        # --- Get domain lat/lon indices...
        lts = np.arange(minlat,maxlat+1,1)
        lns = np.arange(minlon,maxlon+1,1)
        latindex = np.where(lts[:]==minlat)
        lonindex = np.where(lns[:]==minlon)    
        
        forecast_cdf_nc = Dataset('{0}/rfcst2_{1}_{2}_6hr_t190.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        forecast_fhrs = forecast_cdf_nc.variables['fhour'][:]
        forecast_quantiles = forecast_cdf_nc.variables["{}_quantiles".format(field_name)][:,:,:,:]        
        
        analysis_cdf_nc = Dataset('{0}/cfsr_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        analysis_fhrs = analysis_cdf_nc.variables['fhour'][:]
        analysis_quantiles = analysis_cdf_nc.variables["{0}_quantiles".format(field_name)][:,:,:,:]           
        
        
        # --- array of forecast hours
        hrs = np.arange(begin_fcst_hour,end_fcst_hour+1,6)
        
        cfsr_counter = np.where(analysis_fhrs[:] == 12)[0]
        field = np.ones(forecastField.shape[:])*-9999.9
        
        for fhr_idx in xrange(0,hrs.shape[0]):
            if cfsr_counter == (np.where(analysis_fhrs[:] == 12)[0])+2: # --- Past 24 hour mark
                cfsr_counter = np.where(analysis_fhrs[:] == 0)[0]
            rfcst2_fhr_idx = np.where(forecast_fhrs[:] == hrs[fhr_idx])[0] 
            dummy_rfcst_quants = forecast_quantiles[rfcst2_fhr_idx[0],:,:,:]
            dummy_rfcst = np.absolute(forecastField[fhr_idx,:,:])
            idx = find_nearest_idx(dummy_rfcst_quants,dummy_rfcst)
            field[fhr_idx,:,:] = cdf_loop(analysis_quantiles[cfsr_counter,:,:,:],idx,analysis_quantiles.shape[1],lts.shape[0],lns.shape[0])
            cfsr_counter += 1

        # --- Close up shop
        forecast_cdf_nc.close();analysis_cdf_nc.close();



    elif begin_fcst_hour == 180:
        
        # --- Get domain lat/lon indices...
        lts = np.arange(minlat,maxlat+1,1)
        lns = np.arange(minlon,maxlon+1,1)
        latindex = np.where(lts[:]==minlat)
        lonindex = np.where(lns[:]==minlon)    
        
        forecast_cdf_nc = Dataset('{0}/rfcst2_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        forecast_fhrs = forecast_cdf_nc.variables['fhour'][:]
        forecast_quantiles = forecast_cdf_nc.variables["{}_quantiles".format(field_name)][:,:,:,:]        
        
        analysis_cdf_nc = Dataset('{0}/cfsr_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        analysis_fhrs = analysis_cdf_nc.variables['fhour'][:]
        analysis_quantiles = analysis_cdf_nc.variables["{0}_quantiles".format(field_name)][:,:,:,:]        
        
        
        # --- array of forecast hours
        hrs = np.arange(begin_fcst_hour,end_fcst_hour+1,6)
        
        cfsr_counter = np.where(analysis_fhrs[:] == 12)[0]
        field = np.ones(forecastField.shape[:])*-9999.9
        
        for fhr_idx in xrange(0,3):
            if cfsr_counter == (np.where(analysis_fhrs[:] == 12)[0])+2: # --- Past 24 hour mark
                cfsr_counter = np.where(analysis_fhrs[:] == 0)[0]
            rfcst2_fhr_idx = np.where(forecast_fhrs[:] == hrs[fhr_idx])[0] 
            dummy_rfcst_quants = forecast_quantiles[rfcst2_fhr_idx[0],:,:,:]
            dummy_rfcst = np.absolute(forecastField[fhr_idx,:,:])
            idx = find_nearest_idx(dummy_rfcst_quants,dummy_rfcst)
            field[fhr_idx,:,:] = cdf_loop(analysis_quantiles[cfsr_counter,:,:,:],idx,analysis_quantiles.shape[1],lts.shape[0],lns.shape[0])
            cfsr_counter += 1

        # --- Close up shop
        forecast_cdf_nc.close();analysis_cdf_nc.close(); 

        forecast_cdf_nc = Dataset('{0}/rfcst2_{1}_{2}_6hr_t190.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        forecast_fhrs = forecast_cdf_nc.variables['fhour'][:]
        forecast_quantiles = forecast_cdf_nc.variables["{}_quantiles".format(field_name)][:,:,:,:]        
        
        analysis_cdf_nc = Dataset('{0}/cfsr_{1}_{2}_6hr.nc'.format(cdf_dir,forecastDate.strftime("%b"),cdf_name),'r')
        analysis_fhrs = analysis_cdf_nc.variables['fhour'][:]
        analysis_quantiles = analysis_cdf_nc.variables["{0}_quantiles".format(field_name)][:,:,:,:]           
        
        
        # --- array of forecast hours
        hrs = np.arange(begin_fcst_hour,end_fcst_hour+1,6)
        
        for fhr_idx in xrange(3,hrs.shape[0]):
            if cfsr_counter == (np.where(analysis_fhrs[:] == 12)[0])+2: # --- Past 24 hour mark
                cfsr_counter = np.where(analysis_fhrs[:] == 0)[0]
            rfcst2_fhr_idx = np.where(forecast_fhrs[:] == hrs[fhr_idx])[0] 
            dummy_rfcst_quants = forecast_quantiles[rfcst2_fhr_idx[0],:,:,:]
            dummy_rfcst = np.absolute(forecastField[fhr_idx,:,:])
            idx = find_nearest_idx(dummy_rfcst_quants,dummy_rfcst)
            field[fhr_idx,:,:] = cdf_loop(analysis_quantiles[cfsr_counter,:,:,:],idx,analysis_quantiles.shape[1],lts.shape[0],lns.shape[0])
            cfsr_counter += 1

        # --- Close up shop
        forecast_cdf_nc.close();analysis_cdf_nc.close();
    
    return field
