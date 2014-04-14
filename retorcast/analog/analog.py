#!/usr/bin/env python
# The module that does most of the heavy lifting
# Written by F.Alvarez, 1/2014

from datetime import timedelta

import numpy as np

from retorcast.forecast import get_bias_corrected_forecast_sigtor
from training import get_training_data,get_pct
from verification import get_verification_data,get_climo_data
from remapping import remap_probs,prob_members,fcst_quants
from retorcast.plot import plot_probs,plot_probs_timespan
from retorcast.utils import config_setup,get_analog_dates
from fortran_routines import rank_analog,rank_analog_pct


def generate_probabilities(forecastDate,leadtime,**kwargs):

    """
        generate_probabilities(forecastDate,leadtime,**kwargs)
        
        Function used to generate probabilities of tornado occurrence based on a rank-analog method
        and plot them out into pngs.
    
        Required arguments:
            forecastDate - some datetime object, the date of the forecast initialization
            leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)

        Optional arguments (probably shouldn't change these, but whatevs...):
            byear - earliest year for which to look for analogous dates (default = 1985)
            eyear - latest year for which to look for analogous dates (default = 2011)
            train_fname - variable name found in netCDF4 files with the training data (default = "cape06shearcin")
            predictor_name - needed to extract training data from netCDF4 file (default = "Significant_tornado_parameter")
            radius - some radius (km) with which to consider a tornado event to be yes/no from a grid point (default = 240)
            smaller_radius - some radius (km) with which to remap forecasts (default = 80)
            gridpt_window - number of grid points n/s/e/w of forecast grid point to match pattern (Default = 3)
            members - numpy array which has all number of members to keep for each individual forecast (default = lots...)
            return_probs - Boolean, returns numpy array of probabilities if true
            
        Returns:
            Well, nothing really, but it plots out the maps, so we have that going for us... which is nice...
            
    """
    byear = kwargs.get('byear',1985)
    eyear = kwargs.get('eyear',2011)
    train_fname = kwargs.get('train_fname',"1mocdf_cape06shearcin_avg")
    predictor_name = kwargs.get('predictor_name',"Significant_tornado_parameter")
    radius = kwargs.get('radius',240)
    smaller_radius = kwargs.get('radius',80)
    window = kwargs.get('window',3)
    members = kwargs.get('members',np.array([10,15,20,25,30,35,40,45,50,75,100,125,150,175,200,500,1000]))
    return_probs = kwargs.get('return_probs',False)
    use_pct = kwargs.get('use_pct',False)

    td1 = timedelta(days=(leadtime-1),hours=12.)
    td2 = timedelta(days=(leadtime-1),hours=36.)
    print "Initial date: {}".format(forecastDate)
    print "Forecasting from {} --> {}".format(forecastDate+td1,forecastDate+td2)
    
    # --- First, set up all the vars from the config file
    forecast_dir,reforecast_dir,old_forecasts_dir,tordata_dir,cdf_dir,image_dir,\
        maxlat,minlat,maxlon,minlon,allLats,allLons,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon = config_setup()

    # --- Now, get list of analog dates (all datetime objs)
    training_dates = get_analog_dates(forecastDate.year,forecastDate.month,forecastDate.day,byear,eyear,bias_corr=True)

    # --- Get the training data
    trainData = get_training_data(training_dates,leadtime,train_fname,predictor_name,reforecast_dir)

    # --- Get the verification data
    verifData = get_verification_data(training_dates,leadtime,tordata_dir,radius)

    # --- Get the verification data
    climo_dates = get_analog_dates(forecastDate.year,forecastDate.month,forecastDate.day,byear,eyear,bias_corr=False)
    climoData = get_climo_data(climo_dates,leadtime,tordata_dir,smaller_radius)
    
    # --- Get the forecast data
    first_fhour = (12+((leadtime-1)*24)) # --- First forecast hour in 12z-12z format
    last_fhour = first_fhour + 24 # --- Last forecast hour in 12z-12z format
    
    forecastData = get_bias_corrected_forecast_sigtor(forecastDate,forecast_dir,cdf_dir,first_fhour,last_fhour,\
        maxlat,minlat,maxlon,minlon)

    # --- We need to flip the data along the latitude axis.
    #forecastData = np.flipud(forecastData)
    
    # --- Now, let's get us some probabilities!
    if use_pct:

        print "pct!"
        quantField = get_pct(training_dates,leadtime,train_fname,predictor_name,reforecast_dir)
        quantFcst = fcst_quants(forecastDate,forecastData,leadtime,reforecast_dir,cdf_dir)

        temporary_probs = rank_analog_pct(trainData,forecastData,verifData,quantField,quantFcst,trainData.shape[0],members,members.shape[0],\
            allLats.shape[0],allLons.shape[0],fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon,\
            allLats,allLons,window)
            
    if not use_pct:
        temporary_probs = rank_analog(trainData,forecastData,verifData,trainData.shape[0],members,members.shape[0],\
            allLats.shape[0],allLons.shape[0],fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon,\
            allLats,allLons,window)    
    
    # --- Here, we keep n number of forecasts based on the rarity of the SigTor/CAPE*Shear/whatever forecast
    forecastProbs,quantiles = prob_members(forecastDate,forecastData,temporary_probs,leadtime,reforecast_dir,cdf_dir)
    
    print forecastProbs.max()
    # --- Let's remap the probs from a larger -> smaller ROI before we start plotting them      
    smaller_roi_forecastProbs = remap_probs(forecastProbs,leadtime,old_forecasts_dir,smaller_radius)
    
        
    if return_probs == True:
        return smaller_roi_forecastProbs
        
    # --- Now, we generate three forecast plots: One with just the forecast probabilities,
    # --- another with the remapped probabilities, and another with the original
    # --- probs, the remapped probs and two other maps showing atmospheric conditions
    print "Plotting out day {} forecast, initialized at {}".format(leadtime,forecastDate)
    plot_probs(forecastDate,leadtime,forecastProbs,image_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon,radius=radius)    
    plot_probs(forecastDate,leadtime,smaller_roi_forecastProbs,image_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon,radius=smaller_radius)    
    #quadplot_probs_quantiles(forecastDate,leadtime,quantiles,smaller_roi_forecastProbs,\
    #    cdf_dir,image_dir,forecast_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon)
    #quadplot_probs_climo(forecastDate,leadtime,climoData,smaller_roi_forecastProbs,\
    #    cdf_dir,image_dir,forecast_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon)
    #quadplot_probs(forecastDate,leadtime,forecastProbs,smaller_roi_forecastProbs,\
    #    image_dir,forecast_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon)

    print "All done!"
    

def multi_plot_analog(forecastDate,leadtime,):
    
    # --- First, set up all the vars from the config file
    forecast_dir,reforecast_dir,old_forecasts_dir,tordata_dir,cdf_dir,image_dir,\
        maxlat,minlat,maxlon,minlon,allLats,allLons,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon = config_setup()
        
    all_probs = np.ones((10,allLats.shape[0],allLons.shape[0]))*-9999.9
    
    for idx in reversed(xrange(1,11,1)):
        if forecastDate-timedelta(days=(idx-leadtime)) > forecastDate:
            continue
        all_probs[10-idx,...] = generate_probabilities((forecastDate-timedelta(days=(idx-leadtime))),idx,return_probs=True)
    
    plot_probs_timespan(forecastDate,leadtime,all_probs,image_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon,radius=80)    
    