#!/usr/bin/env python

import numpy as np
import ConfigParser

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
    print dict1
    return dict1

def AdditiveConfigSectionMap(Config,section):
    """
    Function used to parse a configuration file and put
    variables found in different sections of a config file into
    a Python dictionary for easy use.
    """
    vars = []
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None

    return dict1


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

    n_vars = int(ConfigSectionMap(Config,"forecast_domain")['ll_fcst_latitude'])

    return forecast_dir,reforecast_dir,old_forecasts_dir,tordata_dir,cdf_dir,image_dir,\
        maxlat,minlat,maxlon,minlon,allLats,allLons,fcst_minlat,fcst_maxlat,fcst_minlon,fcst_maxlon

