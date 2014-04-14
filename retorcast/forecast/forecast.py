#!/usr/bin/env python
# Modules for severe_reforecast.py
# Written by F.Alvarez, 1/2013

"""
 Python module for extracting/calculating data from grib2 files, along
 with some other fun stuff...
     
 CHANGE LOG:
     11/13 - Added 0-1 km SRH, Significant Tornado Parameter (sans LCL), new
             function to collect potential analog dates.
     06/13 - Added 0-3 km storm relative helicity, cleaned up some stuff
     01/13 - Initial version, F.Alvarez

"""

import pygrib as pgb
import numpy as np
from analog.bias_correction import correct_bias

def get_fcst_06shear(date,dir_name,beginning_forecast_hour,ending_forecast_hour\
                        ,maxlat,minlat,maxlon,minlon):
    """

    Used to get U/V wind data from grib files and get the
    vector magnitude vertical shear from a forecast dataset.... that is:

    shear = sqrt((u_level1 - u_level2)**2 +
                      (v_level1 - v_level2)**2)


    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    shear06 - A 2-d (lat/lon) NumPy array of the mean 0-6 km shear values over the forecast hour...

    """
    if ending_forecast_hour <= 180:

        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour > 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour == 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon) 
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)


        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
    
        #Make an array to fill in data values (faster than vstack)
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
            cont_idx+=1                                    
            
        u_lower_level.close();u_upper_level.close()
        v_lower_level.close();v_upper_level.close()

        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
           
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr+cont_idx,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        

    return shear06
                

def get_fcst_plot_06shear(date,dir_name,beginning_forecast_hour,ending_forecast_hour\
                        ,maxlat,minlat,maxlon,minlon):
    """

    Used to get U/V wind data from grib files and get the
    vector magnitude vertical shear from a forecast dataset.... that is:

    shear = sqrt((u_level1 - u_level2)**2 +
                      (v_level1 - v_level2)**2)


    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    shear06 - A 2-d (lat/lon) NumPy array of the mean 0-6 km shear values over the forecast hour...

    """
    if ending_forecast_hour <= 180:

        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = np.mean(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))],axis=0)
        
    elif beginning_forecast_hour > 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
        shear06 = np.mean(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))],axis=0)
        
    elif beginning_forecast_hour == 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon) 
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)


        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
    
        #Make an array to fill in data values (faster than vstack)
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
            cont_idx+=1                                    
            
        u_lower_level.close();u_upper_level.close()
        v_lower_level.close();v_upper_level.close()

        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
           
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr+cont_idx,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = np.mean(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))],axis=0)
        

    return shear06

def get_fcst_plot_cdf_06shear(date,dir_name,cdf_dir,var_name,cdf_name,beginning_forecast_hour,ending_forecast_hour\
                        ,maxlat,minlat,maxlon,minlon):
    """
    Used to get U/V wind data from grib files and get the
    vector magnitude vertical shear from a forecast dataset.... that is:

    shear = sqrt((u_level1 - u_level2)**2 +
                      (v_level1 - v_level2)**2)


    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    shear06 - A 2-d (lat/lon) NumPy array of the mean 0-6 km shear values over the forecast hour...

    """
    if ending_forecast_hour <= 180:

        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]


        shear06 = np.mean(correct_bias(date,var_name,shear06,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                maxlat,minlat,maxlon,minlon),axis=0)  
        
    elif beginning_forecast_hour > 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        #Getting data for all levels

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)

    
        #Make an array to fill in data values (faster than vstack)
        val = np.zeros([len(u_lower_values),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]


        shear06 = np.mean(correct_bias(date,var_name,shear06,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                maxlat,minlat,maxlon,minlon),axis=0)         
    elif beginning_forecast_hour == 180:
    
        #Opening data files...
        
        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
    
    
    
        #Get lat/lon indices...
        grbs = u_lower_level.message(1)
        #Get domain lat/lon indices...
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon) 
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)


        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
    
        #Make an array to fill in data values (faster than vstack)
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
            cont_idx+=1                                    
            
        u_lower_level.close();u_upper_level.close()
        v_lower_level.close();v_upper_level.close()

        u_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_lower_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_10m_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        u_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/ugrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))
        v_upper_level = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/vgrd_pres_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d')))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        u_lower_values = u_lower_level(level=10,forecastTime=hrs)
        v_lower_values = v_lower_level(level=10,forecastTime=hrs)
        u_upper_values = u_upper_level(level=500,forecastTime=hrs)
        v_upper_values = v_upper_level(level=500,forecastTime=hrs)
           
        for hr in xrange(len(u_lower_values)):
            dummy_u_low = np.flipud(u_lower_values[hr].values)
            dummy_v_low = np.flipud(v_lower_values[hr].values)
            dummy_u_high = np.flipud(u_upper_values[hr].values)
            dummy_v_high = np.flipud(v_upper_values[hr].values)
            
            val[hr+cont_idx,:,:] = np.sqrt(((dummy_u_low\
                                    -dummy_u_high)**2)\
                                    +((dummy_v_low\
                                    -dummy_v_high)**2))
                                    
        shear06 = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]


        shear06 = np.mean(correct_bias(date,var_name,shear06,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                maxlat,minlat,maxlon,minlon),axis=0)    

    return shear06



def get_fcst_data(date,dir_name,var_name,level,beginning_forecast_hour,ending_forecast_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that extracts data from some grib file.
    
    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    var_name - name of the variable, specifically how it is found in the grib file name (e.g. cape, apcp)
    level - where this var is found, specifically how it looks in the file name (e.g. sfc, pres)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    var_data - A 2-d (lat/lon) NumPy array of the mean variable data over the forecast hours...

    """

    #Checking the forecast hour...
    if ending_forecast_hour <= 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        #extracted_data = main_grib_data(forecastTime=lambda l: (l <= ending_forecast_hour and l >= beginning_forecast_hour and np.mod(l,6) == 0))
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour > 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        #extracted_data = main_grib_data(forecastTime=lambda l: (l <= ending_forecast_hour and l >= beginning_forecast_hour and np.mod(l,6) == 0))
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)

    
        var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]

    elif beginning_forecast_hour == 180:
        #Getting the forecast data
        #Opening data files...
        #Getting fhour=180 data
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)

        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
            cont_idx+=1
        main_grib_data.close()

        #now getting the rest of the data...
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        hrs = np.arange(198,ending_forecast_hour+1,6)
        extracted_data = main_grib_data(forecastTime=hrs)  
        for hr in xrange(len(extracted_data)):
            val[hr+cont_idx,:,:] = np.flipud(extracted_data[hr].values)
         
        var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        

    return var_data

def get_fcst_plot_data(date,dir_name,var_name,level,beginning_forecast_hour,ending_forecast_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that extracts data from some grib file.
    
    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    var_name - name of the variable, specifically how it is found in the grib file name (e.g. cape, apcp)
    level - where this var is found, specifically how it looks in the file name (e.g. sfc, pres)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    var_data - A 2-d (lat/lon) NumPy array of the mean variable data over the forecast hours...

    """

    #Checking the forecast hour...
    if ending_forecast_hour <= 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        #extracted_data = main_grib_data(forecastTime=lambda l: (l <= ending_forecast_hour and l >= beginning_forecast_hour and np.mod(l,6) == 0))
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
        #val = np.zeros((latshape,lonshape))

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour > 180:
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)    
        
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]

    elif beginning_forecast_hour == 180:

        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(extracted_data)):
            dummy = extracted_data[hr].values
            val[hr,:,:] = np.flipud(dummy)  
            #val[hr,:,:] = np.flipud(dummy)  
            cont_idx+=1
        main_grib_data.close()

        #now getting the rest of the data...
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        extracted_data = main_grib_data(forecastTime=hrs)  
        for hr in xrange(len(extracted_data)):
            val[hr+cont_idx,:,:] = np.flipud(extracted_data[hr].values)
         
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        

    return var_data

def get_fcst_plot_cdf_data(date,dir_name,cdf_dir,var_name,var_short_name,cdf_name,level,beginning_forecast_hour,ending_forecast_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that extracts data from some grib file.
    
    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    var_name - name of the variable, specifically how it is found in the grib file name (e.g. cape, apcp)
    level - where this var is found, specifically how it looks in the file name (e.g. sfc, pres)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    var_data - A 2-d (lat/lon) NumPy array of the mean variable data over the forecast hours...

    """

    #Checking the forecast hour...
    if ending_forecast_hour <= 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        #extracted_data = main_grib_data(forecastTime=lambda l: (l <= ending_forecast_hour and l >= beginning_forecast_hour and np.mod(l,6) == 0))
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
        #val = np.zeros((latshape,lonshape))

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        #var_data = np.mean(val,axis=0)

        if var_short_name == 'cin':        
            var_data = np.absolute(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))])
    
            var_data = np.mean((correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon))*-1.,axis=0)
        else:
            var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
    
            var_data = np.mean(correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon),axis=0)  
                    
                    
    elif beginning_forecast_hour > 180:
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)    
        
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        if var_short_name == 'cin':        
            var_data = np.absolute(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))])
    
            var_data = np.mean((correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon))*-1.,axis=0)
        else:
            var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]

            var_data = np.mean(correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon),axis=0)   

    elif beginning_forecast_hour == 180:

        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(extracted_data)):
            dummy = extracted_data[hr].values
            val[hr,:,:] = np.flipud(dummy)  
            #val[hr,:,:] = np.flipud(dummy)  
            cont_idx+=1
        main_grib_data.close()

        #now getting the rest of the data...
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        extracted_data = main_grib_data(forecastTime=hrs)  
        for hr in xrange(len(extracted_data)):
            val[hr+cont_idx,:,:] = np.flipud(extracted_data[hr].values)
         
        if var_short_name == 'cin':        
            var_data = np.absolute(val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))])
    
            var_data = np.mean((correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon))*-1.,axis=0)
        else:
            var_data = val[:,latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
    
            var_data = np.mean(correct_bias(date,var_short_name,var_data,cdf_dir,cdf_name,beginning_forecast_hour,ending_forecast_hour,\
                    maxlat,minlat,maxlon,minlon),axis=0)  
        

    return var_data



def get_hgt_plot_data(date,dir_name,var_name,level,beginning_forecast_hour,ending_forecast_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that extracts data from some grib file.
    
    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    var_name - name of the variable, specifically how it is found in the grib file name (e.g. cape, apcp)
    level - where this var is found, specifically how it looks in the file name (e.g. sfc, pres)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    var_data - A 2-d (lat/lon) NumPy array of the mean variable data over the forecast hours...

    """

    #Checking the forecast hour...
    if ending_forecast_hour <= 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        #extracted_data = main_grib_data(forecastTime=lambda l: (l <= ending_forecast_hour and l >= beginning_forecast_hour and np.mod(l,6) == 0))
        extracted_data = main_grib_data(level=level,forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
        #val = np.zeros((latshape,lonshape))

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour > 180:
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)    
        
        extracted_data = main_grib_data(level=level,forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
    
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]

    elif beginning_forecast_hour == 180:

        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,6)
    
        extracted_data = main_grib_data(level=level,forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)  
            cont_idx+=1
        main_grib_data.close()

        #now getting the rest of the data...
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        extracted_data = main_grib_data(level=level,forecastTime=hrs)  
        for hr in xrange(len(extracted_data)):
            val[hr+cont_idx,:,:] = np.flipud(extracted_data[hr].values)
         
        var_data = np.mean(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        

    return var_data


def get_fcst_plot_apcp(date,dir_name,var_name,level,beginning_forecast_hour,ending_forecast_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that extracts data from some grib file.
    
    In:

    date - a datetime object, specifically the forecast date
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    var_name - name of the variable, specifically how it is found in the grib file name (e.g. cape, apcp)
    level - where this var is found, specifically how it looks in the file name (e.g. sfc, pres)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain

    Returns:

    var_data - A 2-d (lat/lon) NumPy array of the mean variable data over the forecast hours...

    """

    #Checking the forecast hour...
    
    if ending_forecast_hour <= 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name))
        
        #Getting data for all levels
        hrs = np.arange(beginning_forecast_hour,ending_forecast_hour+1,3)
    
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)
    
        var_data = np.sum(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        
    elif beginning_forecast_hour > 180:
        #Getting the forecast data
        #Opening data files...
        
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        #Getting data for all levels
    
        hrs = np.arange(beginning_forecast_hour+6,ending_forecast_hour+1,6)
    
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)   


        val = np.zeros([len(extracted_data),lts.shape[0],lns.shape[1]])
        #val = np.zeros((latshape,lonshape))

        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)

    
        var_data = np.sum(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]

    elif beginning_forecast_hour == 180:
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))
        
        hrs = np.arange(beginning_forecast_hour+6,ending_forecast_hour+1,6)
    
        extracted_data = main_grib_data(forecastTime=hrs)
        
        #Get domain lat/lon indices...
        grbs = main_grib_data.message(1)
        lts,lns = grbs.latlons()
        lts = np.flipud(lts)
        latindex = np.where(lts[:,0]==minlat)
        lonindex = np.where(lns[0,:]==minlon)  
    
        data_time_shape = (ending_forecast_hour - beginning_forecast_hour)/6
        val = np.zeros([(data_time_shape)+1,lts.shape[0],lns.shape[1]])
        cont_idx=0
        for hr in xrange(len(extracted_data)):
            val[hr,:,:] = np.flipud(extracted_data[hr].values)  
            cont_idx+=1
        main_grib_data.close()

        #now getting the rest of the data...
        main_grib_data = pgb.open('{0}/{1}/{2}/{3}00/mean/latlon/{4}_{3}00_mean_t190.grib2'.format(dir_name,date.year,date.strftime('%Y%m'),date.strftime('%Y%m%d'),var_name,level))

        hrs = np.arange(198,ending_forecast_hour+1,6)

        extracted_data = main_grib_data(forecastTime=hrs)  
        for hr in xrange(len(extracted_data)):
            val[hr+cont_idx,:,:] = np.flipud(extracted_data[hr].values)
         
        var_data = np.sum(val,axis=0)
        var_data = var_data[latindex[0]:int(latindex[0]+(maxlat-minlat+1)),lonindex[0]:int(lonindex[0]+(maxlon-minlon+1))]
        

    return var_data

def get_forecast_sigtor(dattim,dir_name,begin_fcst_hour,end_fcst_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that calculates SigTor.
    
    In:

    dattim - datetime object
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain
    sigtor - A 3-d (time,lat,lon) NumPy array of the mean variable data over the forecast hours...

    Returns:

    sigtor - A 3-d (time,lat,lon) NumPy array of the mean variable data over the forecast hours...

    """
    
    # --- Define some arrays to make this faster...
    #allLats = np.arange(minlat,maxlat+1,1)
    #allLons = np.arange(minlon,maxlon+1,1)  
    #dumCape = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    #dumShear = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    #dumCin = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    
    
    # --- CAPE
    dumCape = get_fcst_data(dattim,dir_name,'cape_sfc',0,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon)/1500.
    # --- 0-6 km Shear                           
    dumShear = get_fcst_06shear(dattim,dir_name,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon)/20.
                                               
    # --- Some 0-6 km Shear restrictions so it doesn't influence SigTor too much                                           
    msk = (dumShear > 1.5)
    dumShear[msk] = 1.5
    msk = (dumShear < 0.625)
    dumShear[msk] = 0.0
    
    # --- CIN
    dumCin = (250. + (get_fcst_data(dattim,dir_name,'cin_sfc',0,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon)))/200.

    # --- Now, calculate SigTor and make sure there aren't any negative values
    sigtor = dumCape*dumShear*dumCin
    sigtor[sigtor < 0.] = 0.
    sigtor = np.mean(sigtor,axis=0)
    return sigtor

def get_bias_corrected_forecast_sigtor(dattim,dir_name,cdf_dir,begin_fcst_hour,end_fcst_hour,maxlat,minlat,maxlon,minlon):
    """
    Function that calculates SigTor.
    
    In:

    dattim - datetime object
    dir_name - top directory of where grib files are found (e.g. /mnt_emc/Reforecast2)
    cdf_dir - top directory of where cdf files are found (e.g. ./data/cdfs/)
    beginning_forecast_hour - First hour of the forecast (e.g. 12)
    ending_forecast_hour - Last hour of the forecast (e.g. 36)
    maxlat/minlat/maxlon/minlon: corners of domain
    sigtor - A 3-d (time,lat,lon) NumPy array of the mean variable data over the forecast hours...

    Returns:

    sigtor - A 3-d (time,lat,lon) NumPy array of the mean variable data over the forecast hours...

    """
    
    # --- Define some arrays to make this faster...
    #allLats = np.arange(minlat,maxlat+1,1)
    #allLons = np.arange(minlon,maxlon+1,1)  
    #dumCape = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    #dumShear = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    #dumCin = np.zeros((allLats.shape[0],allLons.shape[0])) # --- Temp array for SigTor Calc
    
    
    # --- CAPE
    temp_dumCape = get_fcst_data(dattim,dir_name,'cape_sfc',0,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon)
                                               
    dumCape = (correct_bias(dattim,'cape',temp_dumCape,cdf_dir,'1mo_cdf_cape',begin_fcst_hour,\
        end_fcst_hour,maxlat,minlat,maxlon,minlon))/1500.
                                               
                                               
    # --- 0-6 km Shear                           
    temp_dumShear = get_fcst_06shear(dattim,dir_name,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon)

    dumShear = (correct_bias(dattim,'06shear',temp_dumShear,cdf_dir,'1mo_cdf_06shear',begin_fcst_hour,\
        end_fcst_hour,maxlat,minlat,maxlon,minlon))/20.
                                               
    # --- Some 0-6 km Shear restrictions so it doesn't influence SigTor too much                                           
    msk = (dumShear > 1.5)
    dumShear[msk] = 1.5
    msk = (dumShear < 0.625)
    dumShear[msk] = 0.0
    
    # --- CIN
    temp_dumCin = np.absolute(get_fcst_data(dattim,dir_name,'cin_sfc',0,begin_fcst_hour,end_fcst_hour,\
                                               maxlat,minlat,maxlon,minlon))

    temp_dumCin = correct_bias(dattim,'cin',temp_dumCin,cdf_dir,'1mo_cdf_cin',begin_fcst_hour,\
        end_fcst_hour,maxlat,minlat,maxlon,minlon)

    dumCin = (250. + (temp_dumCin*-1))/200.

    # --- Now, calculate SigTor and make sure there aren't any negative values
    sigtor = dumCape*dumShear*dumCin
    sigtor[sigtor < 0.] = 0.
    sigtor = np.mean(sigtor,axis=0)
    return sigtor
        
        