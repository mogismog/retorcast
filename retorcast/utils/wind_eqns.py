import numpy as np

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
