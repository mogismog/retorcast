#!/usr/bin/env python

import numpy as np
from datetime import datetime,timedelta
import ConfigParser
import itertools

# --- Used for percentiles
def find_nearest_idx(array,value):
    idx = (np.abs(array-value)).argmin(axis=0)
    return idx

def ConfigSectionMap(Config,section):
    """
    Function used to parse a configuration file and put
    variables found in different sections of a config file into
    a Python dictionary for easy use.
    """
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

def get_1mo_dates(inyr,inmo,indate,byear,eyear):
    """
 Used with netCDF4 files and py-netCDF. 

 From 1985-2011, a number of days range to search
 for dates centered around the given date and forecast time,
 returns the applicable dates in a list of daetime objects. 
 In other words, returns a list of file name dates for us to 
 search for analogs/use in logistic regression/whatever.

 indate - Initial date (1,31)
 inyr - Initial year, YYYY (1985 - )
 inmo - Initial month, (1,12)
 window - range of dates in past years to search, e.g. 45 will find dates 45 days before/after indate

 Returns:
 outdates - List of dates meeting the criteria
    """
    
    fnlist = []
    
    #print inmo,indate
    try:
        xdate = datetime(byear,inmo,indate)
    except ValueError:
        xdate = datetime(byear,inmo,indate-1)
    else:
        xdate = datetime(byear,inmo,indate)
    while xdate < datetime(eyear+1,1,1):
        #print xdate
        if xdate.year == inyr:
            try:
                xdate = datetime((xdate.year + 1),inmo,indate)
            except ValueError:
                xdate = datetime((xdate.year + 1),inmo,indate-1)
            continue
        for datechange in xrange(0,35):
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                #print analogdate,xdate
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate.month != xdate.month:
                    continue
                fnlist.append(analogdate)                
        try:
            xdate = datetime((xdate.year + 1),inmo,indate)
        except ValueError:
            xdate = datetime(xdate.year+1,inmo,indate-1)
          
    return fnlist

def get_analog_dates(inyr,inmo,indate,byear,eyear,**kwargs):
    """
 Used with netCDF4 files and py-netCDF. 

 From 1985-2011, a number of days range to search
 for dates centered around the given date and forecast time,
 returns the applicable dates in a list of daetime objects. 
 In other words, returns a list of file name dates for us to 
 search for analogs/use in logistic regression/whatever.

 indate - Initial date (1,31)
 inyr - Initial year, YYYY (1985 - )
 inmo - Initial month, (1,12)
 window - range of dates in past years to search, e.g. 45 will find dates 45 days before/after indate
 byear - earliest year for potential dates (usually 1985)
 eyear - latest year for potential dates (usually 2011)

 Returns:
 outdates - List of dates meeting the criteria
    """

    bias_corr = kwargs.get('bias_corr',False)

    
    fnlist = []
    date_list = []
    #print inmo,indate
    try:
        xdate = datetime(byear,inmo,indate)
    except ValueError:
        xdate = datetime(byear,inmo,indate-1)
    else:
        xdate = datetime(byear,inmo,indate)
    while xdate < datetime(eyear+1,1,1):
        #print xdate
        if xdate.year == inyr:
            try:
                xdate = datetime((xdate.year + 1),inmo,indate)
            except ValueError:
                xdate = datetime((xdate.year + 1),inmo,indate-1)
            continue
        for datechange in reversed(xrange(0,100)):
            if xdate.month > 1:
                tdelta = timedelta(days=datechange)
                analogdate = xdate - tdelta
                if analogdate.year < byear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate < datetime(xdate.year,xdate.month-1,1):
                    continue
                fnlist.append(analogdate)
            elif xdate.month == 1:
                tdelta = timedelta(days=datechange)
                analogdate = xdate - tdelta
                if analogdate.year < byear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate < datetime(xdate.year-1,12,1):
                    continue
                fnlist.append(analogdate)                
        for datechange in xrange(1,101):
            if xdate.month < 12:
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                try:
                    datetime(xdate.year,xdate.month+2,1)
                except ValueError: # --- xdate.month == 11
                    if analogdate >= datetime(xdate.year+1,1,1):
                        continue
                else:
                    if analogdate >= datetime(xdate.year,xdate.month+2,1):
                        continue
                fnlist.append(analogdate)
            elif xdate.month == 12:
                tdelta = timedelta(days=datechange)
                analogdate = xdate + tdelta
                if analogdate.year > eyear:
                    continue
                if analogdate.year == inyr:
                    continue
                if analogdate >= datetime(xdate.year+1,2,1):
                    continue
                fnlist.append(analogdate)                
        try:
            xdate = datetime((xdate.year + 1),inmo,indate)
        except ValueError:
            xdate = datetime(xdate.year+1,inmo,indate-1)


    # --- Here, since we are now using bias-corrected data, we can get additional potetial analog dates!
    if bias_corr:

        date_list.append(fnlist)    
        
        #for n_mo in xrange(1,13,1):
        #    if (n_mo >= inmo-1) and (n_mo <= inmo+1):
        #        continue
        #    else:
        #        date_list.append(get_1mo_dates(int(inyr),n_mo,1,byear,eyear))
        if (inmo < 2) or (inmo > 9):
           date_list.append(get_1mo_dates(int(inyr),3,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),4,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))
        if (inmo == 2):
           date_list.append(get_1mo_dates(int(inyr),4,1,byear,eyear))
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))          
        if (inmo == 3):
           date_list.append(get_1mo_dates(int(inyr),5,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))
        if (inmo == 4):
           date_list.append(get_1mo_dates(int(inyr),9,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),10,1,byear,eyear))          
           date_list.append(get_1mo_dates(int(inyr),11,1,byear,eyear))

        # --- Now flatten and return the list
        date_list = list(itertools.chain.from_iterable(date_list)) 
        return date_list         
    else:
        return fnlist

def winddir(u,v):
    """
    Returns the wind direction in degrees, given a u and v value.

    in:
       u - zonal wind component.
       v - meridional wind component.
    out:
       wind_direction - direction of wind in degrees.
    """
    wind_direction = (270-(np.arctan2(v,u)*(180/np.pi)))%360
    return wind_direction

def windspeed(u,v):
    """
    Returns the wind speed in m/s, given u/v component vectors.

    in:
       u - zonal wind component.
       v - meridional wind component.
    out:
       wind_speed - speed of wind in m/s.
    """
    wind_speed=np.sqrt((u**2)+(v**2))
    return wind_speed

def ucomp(wind_speed,wind_direction):
    """
    Returns the u component magnitude in m/s.

    in:
       wind_speed - wind_speed in m/s.
       wind_direction - wind direction in degrees.
    out:
       u_component - u-component of the wind.
    """
    wind_direction = wind_direction*(np.pi/180)
    u_component = wind_speed * np.sin(wind_direction)*-1
    return u_component

def vcomp(wind_speed,wind_direction):
    """
    Returns the v component magnitude in m/s.

    in:
       wind_speed - wind_speed in m/s.
       wind_direction - wind direction in degrees.
    out:
       v_component - v-component of the wind.
    """
    wind_direction = wind_direction*(np.pi/180)
    v_component = wind_speed * np.cos(wind_direction)*-1
    return v_component

def storm_motion_spddir(u,v):
    """
    Returns storm motion speed and direction given u,v of
    some level defined as mean wind. Uses 75/30 method.

    in:
       u - zonal wind component.
       v - meridional wind component.
    out:
       sm_speed - storm-motion wind speed.
       sm_direction - storm-motion wind direction.
    """
    speed = windspeed(u,v)
    direction = winddir(u,v)
    sm_speed = (75*speed)/100
    sm_direction = (direction+30)%360
    return sm_speed,sm_direction


def srwind_comps(u,v,sm_speed,sm_direction):
    """
    Returns the storm-relative wind component vectors.

    in:
       u - zonal wind component.
       v - meridional wind component.
       sm_speed - storm-motion wind speed.
       sm_direction - storm-motion wind direction.
    out:
       sru - storm-relative zonal wind component.
       srv - storm-relative meridional wind component.
    """
    smu = ucomp(sm_speed,sm_direction)
    smv = vcomp(sm_speed,sm_direction)
    sru = u-smu
    srv = v-smv
    return sru,srv

def srh(u,v,mean_u,mean_v):
    """
    Returns the storm-relative helicity (SRH) given a vector of
    u/v components at different levels and some mean-wind u/v at some
    level (e.g. 0-6km, 700 mb, etc.) from which we want to get
    the storm motion.

    in:
       u - zonal wind component.
       v - meridional wind component.
       mean_u - mean zonal wind component.
       mean_v - mean meridional wind component.
    out:
       storm_relative_helicity - in m2 s-2.
    """
    if u.shape != v.shape:
        raise ValueError, "u & v should be the same shape!"
    
    #Mean wind and wind difference between layers
    ndim = len(u.shape)
    shp = [(u.shape[0]-1)]
    for i in xrange(1,ndim,1):
        shp.append(u.shape[i])
    mean_u_layers = np.zeros((shp))
    mean_v_layers =  np.zeros((shp))
    delta_u_layers =  np.zeros((shp))
    delta_v_layers =  np.zeros((shp))
    u = np.asarray(u)
    v = np.asarray(v)
    for layer in xrange(shp[0]):
        mean_u_layers[layer,:] = (u[layer,:]+u[layer+1,:])/2
        mean_v_layers[layer,:] = (v[layer,:]+v[layer+1,:])/2
        delta_u_layers[layer,:] = (u[layer+1,:]-u[layer,:])
        delta_v_layers[layer,:] = (v[layer+1,:]-v[layer,:])
    #Storm-motion speed and direction
    sm_speed,sm_direction = storm_motion_spddir(mean_u,mean_v)

    #Storm-relative u,v
    sru,srv = srwind_comps(mean_u_layers,mean_v_layers,sm_speed,sm_direction)

    #Storm_relative_helicity
    storm_relative_helicity = np.sum(((srv*delta_u_layers) - (sru*delta_v_layers)),axis=0)
    return storm_relative_helicity


def config_setup():
    # --- First, parse the config.ini file...

    Config = ConfigParser.ConfigParser()
    Config.read('./setup.cfg')
    
    # --- Now, extract variables from config file
    
    forecast_dir = ConfigSectionMap(Config,"path_names")['forecast_dir']
    reforecast_dir = ConfigSectionMap(Config,"path_names")['reforecast_dir']
    old_forecasts_dir = ConfigSectionMap(Config,"path_names")['old_forecasts_dir']
    tordata_dir = ConfigSectionMap(Config,"path_names")['tordata_dir']
    cdf_dir = ConfigSectionMap(Config,"path_names")['cdf_dir']
    image_dir = ConfigSectionMap(Config,"path_names")['image_dir']
    
    # --- Basically, defining a domain around the CONUS
    maxlat = float(ConfigSectionMap(Config,"total_domain")['ur_latitude'])
    minlat = float(ConfigSectionMap(Config,"total_domain")['ll_latitude'])
    maxlon = float(ConfigSectionMap(Config,"total_domain")['ur_longitude'])
    minlon = float(ConfigSectionMap(Config,"total_domain")['ll_longitude'])
    allLats = np.arange(minlat,maxlat+1,1)
    allLons = np.arange(minlon,maxlon+1,1)
    
    # --- Now, let's define the forecast domain
    fcst_minlat = float(ConfigSectionMap(Config,"forecast_domain")['ll_fcst_latitude'])
    fcst_maxlat = float(ConfigSectionMap(Config,"forecast_domain")['ur_fcst_latitude']) 
    fcst_minlon = float(ConfigSectionMap(Config,"forecast_domain")['ll_fcst_longitude'])
    fcst_maxlon = float(ConfigSectionMap(Config,"forecast_domain")['ur_fcst_longitude'])
    
    return forecast_dir,reforecast_dir,old_forecasts_dir,tordata_dir,cdf_dir,image_dir,\
        maxlat,minlat,maxlon,minlon,allLats,allLons,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon

