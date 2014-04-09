#!/usr/bin/env python

import numpy as np
from datetime import datetime,timedelta

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from forecast import get_fcst_plot_06shear,get_fcst_plot_data,get_hgt_plot_data,get_fcst_plot_apcp,get_fcst_plot_cdf_06shear,get_fcst_plot_cdf_data
matplotlib.rcParams['lines.linewidth'] = .25

def plot_probs(forecastDate,leadtime,probs,image_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon,**kwargs):
    """
        plot_probs(forecastDate,leadtime,probs,image_dir,**kwargs)
  
      A function that plots out the forecast probabilities trained on a larger ROI on a single panel.
        
        Required arguments:
            forecastDate - some datetime object, the date of the forecast initialization
            leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)        
            probs - some 2-D numpy array of forecast probabilities [0,100]
            image_dir - path to which the images are saved
            minlat,maxlat,minlon,maxlon - forecast domain
        Optional Arguments:
            radius - some radius for the tornado probabilities (default = 240km)  
    """
    radius = kwargs.get('radius',240)    


    lns = np.arange(minlon,maxlon+1)
    lts = np.arange(minlat,maxlat+1)
    lns,lts = np.meshgrid(lns,lts)

    # --- Get forecast/analysis times in an ok looking format
    start_fdattim_str = (forecastDate+timedelta(days=int(leadtime)-1,hours=12)).strftime("%Y-%m-%d %H UTC")
    end_fdattim_str = (forecastDate+timedelta(days=int(leadtime),hours=12)).strftime("%Y-%m-%d %H UTC")
    current_fdattim_str = forecastDate.strftime("%Y-%m-%d %H UTC")

    
    fig1 = plt.figure(figsize=(4.8,4.4))
    #title = 'Tornado probabilities (EF1+), {0} km ROI.\n {1} to {2}\n'.format(radius,start_fdattim_str,end_fdattim_str)+\
    #    'Initialization time = {0} UTC'.format(current_fdattim_str)
    colorst = ['BurlyWood','LightGray','LightSkyBlue','LightGreen','Green','Gold',\
               'Orange','Red','Plum','Orchid']
    if radius == 80:
        clevs=[1,2,3,4,5,7,10,12,15,20,25]
    elif radius == 160:
        clevs=[2,5,7,10,12,15,20,25,30,35,40]
    else:
        clevs=[5,10,15,20,25,30,35,40,45,50,60]
    
    # ---- make first plot, the probability forecast
    
    ax = fig1.add_axes([0.06,.13,0.88,.76])
    #ax.set_title(title,fontsize=9)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    #m.fillcontinents(color='beige')
    x, y = m(lns, lts)
    CS1 = m.contour(x,y,probs,levels=clevs,linewidth=0.3,colors='k',cmap=None)
    CS2 = m.contourf(x,y,probs,levels=clevs,colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)   

    #ax.text(1.0,-0.025,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)

### fig1 = plt.figure(figsize=(7.5,8.))
### ax = fig1.add_axes([0.12,0.12,0.76,0.76])
### plt.figtext(0.88, 0.16, psdtitle, ha='right', color='black', size='small')

    cax = fig1.add_axes([0.1,0.08,0.8,0.03])
    cbar = fig1.colorbar(CS2,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
    cbar.ax.tick_params(labelsize=8)
    cax.set_xlabel('Forecast tornado probability (%)',fontsize=9)
    
    plt.savefig('{0}tornado_fcst_{2}km_{1}_day{3}.png'.format(image_dir,forecastDate.strftime("%Y%m%d"),radius,leadtime),dpi=300)
    plt.clf()
    plt.close()


def quadplot_probs_quantiles(forecastDate,leadtime,quantiles,smaller_roi_probs,cdf_dir,image_dir,forecast_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon,**kwargs):
    """
        plot_probs(forecastDate,leadtime,probs,image_dir,**kwargs)
  
      A function that plots out the forecast probabilities trained on a larger ROI on a single panel.
        
        Required arguments:
            forecastDate - some datetime object, the date of the forecast initialization
            leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)        
            larger_roi_probs - some 2-D numpy array of forecast probabilities [0,100]
            smaller_roi_probs - some 2-D numpy array of forecast probabilities [0,100]
            image_dir - path to which the images are saved
            minlat,maxlat,minlon,maxlon - forecast domain
        Optional Arguments:
            larger_radius - some radius for the tornado probabilities (default = 240km)  
            smaller_radius - some radius for the tornado probabilities (default = 80km)  
    """
    larger_radius = kwargs.get('larger_radius',240)    
    smaller_radius = kwargs.get('smaller_radius',80)    

    lns = np.arange(minlon,maxlon+1)
    lts = np.arange(minlat,maxlat+1)
    lns,lts = np.meshgrid(lns,lts)
    
    # --- Get forecast/analysis times in an ok looking format
    start_fdattim_str = (forecastDate+timedelta(days=int(leadtime)-1,hours=12)).strftime("%Y-%m-%d %H UTC")
    end_fdattim_str = (forecastDate+timedelta(days=int(leadtime),hours=12)).strftime("%Y-%m-%d %H UTC")
    current_fdattim_str = forecastDate.strftime("%Y-%m-%d %H UTC")

    # --- Get the forecast data
    first_fhour = (12+((leadtime-1)*24)) # --- First forecast hour in 12z-12z format
    last_fhour = first_fhour + 24 # --- Last forecast hour in 12z-12z format
    
    # --- CAPE
    cape = get_fcst_plot_cdf_data(forecastDate,forecast_dir,cdf_dir,'cape_sfc','cape','1mo_cdf_cape',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- 0-6 km Shear                           
    shear = get_fcst_plot_cdf_06shear(forecastDate,forecast_dir,cdf_dir,'06shear','1mo_cdf_06shear',first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- CIN
    cin = get_fcst_plot_cdf_data(forecastDate,forecast_dir,cdf_dir,'cin_sfc','cin','1mo_cdf_cin',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- 500 mb Heights
    hgts_500mb = get_hgt_plot_data(forecastDate,forecast_dir,'hgt_pres',500,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)
    # --- mslp
    mslp = get_fcst_plot_data(forecastDate,forecast_dir,'pres_msl',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- precip
    precip = get_fcst_plot_apcp(forecastDate,forecast_dir,'apcp_sfc',0,first_fhour+6,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)
                                               
    fig1 = plt.figure(figsize=(7.2,7.5))
    title = '(a) Tornado Probabilities (F1+), {0} km ROI\n{1} to {2} \n'.format(smaller_radius,start_fdattim_str,end_fdattim_str)+\
        'Initialization time = {0}'.format(current_fdattim_str)
    colorst = ['BurlyWood','LightGray','LightSkyBlue','LightGreen','Green','Gold',\
               'Orange','Red','Plum','Orchid']
    clevs=[1,2,3,4,5,7,10,12,15,20,25]  
    
    # ---- make first plot, the probability forecast
    
    ax = fig1.add_axes([0.02,.585,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    #m.fillcontinents(color='beige')
    x, y = m(lns,lts)
    CS1 = m.contour(x,y,smaller_roi_probs,levels=clevs,linewidth=0.3,colors='k',cmap=None)
    CS2 = m.contourf(x,y,smaller_roi_probs,levels=clevs,colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    ax.text(1.0,-0.042,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)
    
    cax = fig1.add_axes([0.1,0.55,0.3,0.015])
    cbar = fig1.colorbar(CS2,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Forecast tornado probability (%)',fontsize=7)
    
    # ---- make second plot, the quantiles of the forecast
    
    title = '(c) Quantiles of sbCAPE*0-6 km Shear*CIN\n{0} to {1}\n'.format(start_fdattim_str,end_fdattim_str)+\
        'Initialization time = {0}'.format(current_fdattim_str)
    clevs = [80,85,90,95, 96, 97, 98, 99,99.5, 99.9]
    ax = fig1.add_axes([0.02,.08,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    #m.fillcontinents(color='beige')
    x,y = m(lns,lts)

    CS1 = m.contour(x,y,quantiles*100.,levels=clevs,linewidth=0.1,colors='k',cmap=None)
    CS2 = m.contourf(x,y,quantiles*100.,levels=clevs,colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)

    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)
    
    cax = fig1.add_axes([0.1,0.045,0.3,0.015])
    cbar = fig1.colorbar(CS2,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=[80,85,90,95, 96, 97, 98, 99,99.5, 99.9],format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Quantile',fontsize=7)
    
    # ---- third plot, MSLP, precip, and Z500
    
    title = '(b) Daily Avg MSLP, 500 hPa Heights, Total Precip\n{0} to {1}\n'.format(start_fdattim_str,end_fdattim_str)+\
            'Initialization time = {0}'.format(current_fdattim_str)
    ax = fig1.add_axes([0.52,.585,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    x,y = m(lns,lts)
    matplotlib.rcParams['lines.linewidth'] = 0.7
    CS1 = m.contour(x,y,mslp/100.,levels=range(940,1061,4),colors='Black',cmap=None,label='Sea-level pressure')
    ax.clabel(CS1,fmt='%4d',fontsize=7)
    CS2 = m.contour(x,y,hgts_500mb/10.,levels=range(480,600,6),linewidth=0.7,colors='Gray',cmap=None,label='500 hPa height')
    ax.clabel(CS2,fmt='%4d',fontsize=7)
    #xo, yo = m(lns, lts)
    matplotlib.rcParams['lines.linewidth'] = .25
    
    CS3 = m.contour(x,y,precip,levels=[1,2,5,10,25,50],linewidth=0.1,colors='Black',cmap=None)
    CS4 = m.contourf(x,y,precip,levels=[1,2,5,10,25,50],colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    
    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)

    cax = fig1.add_axes([0.6,0.55,0.3,0.015])
    cbar = fig1.colorbar(CS4,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=[1,2,5,10,25,50],format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Ensemble-mean precipitation amount (mm)',fontsize=7)
    
    # ---- fourth plot, CAPE, CIN, and 0-1 km shear
    
    
    title = '(d) Daily Avg sbCAPE, CIN & 0-6 km Bulk Shear\n{0} to {1}\n'.format(start_fdattim_str,end_fdattim_str)+\
            'Initialization time = {0}'.format(current_fdattim_str)

    ax = fig1.add_axes([0.52,.08,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    matplotlib.rcParams['lines.linewidth'] = .7
    CS1 = m.contour(x,y,cape,levels=[50,125,250,500,1000,1500,2000,2500,3000,4000],linewidth=.7,colors='Black',cmap=None,label='CAPE')
    CS4 = m.contour(x,y,shear,levels=range(10,80,5),linewidth=.7,colors='Black',cmap=None,label='Shear')
    ax.clabel(CS4,fmt='%3d',fontsize=7)
    CS2 = m.contourf(x,y,cape,levels=[50,125,250,500,1000,1500,2000,2500,3000,4000],colors=colorst,cmap=None,label='CAPE',extend="max")
    #ax.clabel(CS2,fmt='%2d',fontsize=9)
    CS3 = m.contour(x,y,cin*-1,levels=[50,100,150,250,500],\
                            linewidth=0.5,colors='Gray',cmap=None,label='CIN')
    ax.clabel(CS3,fmt='%3d',fontsize=7)
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    
    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)

    cax = fig1.add_axes([0.57,0.045,0.36,0.015])
    cbar = fig1.colorbar(CS2,extend='max', \
       orientation='horizontal',cax=cax,drawedges=True,\
       ticks=[50,125,250,500,1000,1500,2000,2500,3000,4000],format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('sbCAPE',fontsize=7)

    plt.savefig('{0}tornado_fcst_4panel_{1}_day{2}.png'.format(image_dir,forecastDate.strftime("%Y%m%d"),leadtime),dpi=200)
    plt.clf()
    plt.close()
  
def quadplot_probs_climo(forecastDate,leadtime,climo,smaller_roi_probs,cdf_dir,image_dir,forecast_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon,**kwargs):
    """
        plot_probs(forecastDate,leadtime,probs,image_dir,**kwargs)
  
      A function that plots out the forecast probabilities trained on a larger ROI on a single panel.
        
        Required arguments:
            forecastDate - some datetime object, the date of the forecast initialization
            leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)        
            larger_roi_probs - some 2-D numpy array of forecast probabilities [0,100]
            smaller_roi_probs - some 2-D numpy array of forecast probabilities [0,100]
            image_dir - path to which the images are saved
            minlat,maxlat,minlon,maxlon - forecast domain
        Optional Arguments:
            larger_radius - some radius for the tornado probabilities (default = 240km)  
            smaller_radius - some radius for the tornado probabilities (default = 80km)  
    """
    smaller_radius = kwargs.get('smaller_radius',80)    

    lns = np.arange(minlon,maxlon+1)
    lts = np.arange(minlat,maxlat+1)
    lns,lts = np.meshgrid(lns,lts)
    
    # --- Get forecast/analysis times in an ok looking format
    start_fdattim_str = (forecastDate+timedelta(days=int(leadtime)-1,hours=12)).strftime("%Y-%m-%d %H UTC")
    end_fdattim_str = (forecastDate+timedelta(days=int(leadtime),hours=12)).strftime("%Y-%m-%d %H UTC")
    current_fdattim_str = forecastDate.strftime("%Y-%m-%d %H UTC")
    
    if forecastDate.month == 1:
        start_month_str = (datetime(forecastDate.year,12,1).strftime("%b"))
        end_month_str = (datetime(forecastDate.year,forecastDate.month+1,1).strftime("%b"))
    elif forecastDate.month == 12:
        start_month_str = (datetime(forecastDate.year,forecastDate.month-1,1).strftime("%b"))
        end_month_str = (datetime(forecastDate.year,1,1).strftime("%b"))        
    else:
        start_month_str = (datetime(forecastDate.year,forecastDate.month-1,1).strftime("%b"))
        end_month_str = (datetime(forecastDate.year,forecastDate.month+1,1).strftime("%b"))
        
        
    # --- Get the forecast data
    first_fhour = (12+((leadtime-1)*24)) # --- First forecast hour in 12z-12z format
    last_fhour = first_fhour + 24 # --- Last forecast hour in 12z-12z format
    
    # --- CAPE
    cape = get_fcst_plot_cdf_data(forecastDate,forecast_dir,cdf_dir,'cape_sfc','cape','1mo_cdf_cape',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- 0-6 km Shear                           
    shear = get_fcst_plot_cdf_06shear(forecastDate,forecast_dir,cdf_dir,'06shear','1mo_cdf_06shear',first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- CIN
    cin = get_fcst_plot_cdf_data(forecastDate,forecast_dir,cdf_dir,'cin_sfc','cin','1mo_cdf_cin',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- 500 mb Heights
    hgts_500mb = get_hgt_plot_data(forecastDate,forecast_dir,'hgt_pres',500,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)
    # --- mslp
    mslp = get_fcst_plot_data(forecastDate,forecast_dir,'pres_msl',0,first_fhour,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)

    # --- precip
    precip = get_fcst_plot_apcp(forecastDate,forecast_dir,'apcp_sfc',0,first_fhour+6,last_fhour,\
                                               maxlat,minlat,maxlon,minlon)
                                               
    fig1 = plt.figure(figsize=(7.2,7.5))
    title = '(a) Tornado Probabilities (F1+), {0} km ROI\n{1} to {2} \n'.format(smaller_radius,start_fdattim_str,end_fdattim_str)+\
        'Initialization time = {0}'.format(current_fdattim_str)
    colorst = ['BurlyWood','LightGray','LightSkyBlue','LightGreen','Green','Gold',\
               'Orange','Red','Plum','Orchid']
    clevs=[1,2,3,4,5,7,10,12,15,20,25]  
    
    # ---- make first plot, the probability forecast
    
    ax = fig1.add_axes([0.02,.585,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    #m.fillcontinents(color='beige')
    x, y = m(lns,lts)
    CS1 = m.contour(x,y,smaller_roi_probs,levels=clevs,linewidth=0.3,colors='k',cmap=None)
    CS2 = m.contourf(x,y,smaller_roi_probs,levels=clevs,colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    ax.text(1.0,-0.042,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)
    
    cax = fig1.add_axes([0.1,0.55,0.3,0.015])
    cbar = fig1.colorbar(CS2,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Forecast tornado probability (%)',fontsize=7)
    
    # ---- make second plot, the quantiles of the forecast
    
    title = '(c) F1+ Tornado Raw Climatological Probability\n{0}-{1}-{2}, 1985-2011'.format(start_month_str,forecastDate.strftime('%b'),end_month_str)
    clevs=[1,2,3,4,5,7,10,12,15,20,25]  
    #clevs=[.5,1,1.5,2,3,4,5,7,10,12,15,]  
    ax = fig1.add_axes([0.02,.08,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    #m.fillcontinents(color='beige')
    x,y = m(lns,lts)

    CS1 = m.contour(x,y,climo,levels=clevs,linewidth=0.1,colors='k',cmap=None)
    CS2 = m.contourf(x,y,climo,levels=clevs,colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)

    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)
    
    cax = fig1.add_axes([0.1,0.045,0.3,0.015])
    cbar = fig1.colorbar(CS2,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Climatological Probability',fontsize=7)
    
    # ---- third plot, MSLP, precip, and Z500
    
    title = '(b) Daily Avg MSLP, 500 hPa Heights, Total Precip\n{0} to {1}\n'.format(start_fdattim_str,end_fdattim_str)+\
            'Initialization time = {0}'.format(current_fdattim_str)
    ax = fig1.add_axes([0.52,.585,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    x,y = m(lns,lts)
    matplotlib.rcParams['lines.linewidth'] = 0.7
    CS1 = m.contour(x,y,mslp/100.,levels=range(940,1061,4),colors='Black',cmap=None,label='Sea-level pressure')
    ax.clabel(CS1,fmt='%4d',fontsize=7)
    CS2 = m.contour(x,y,hgts_500mb/10.,levels=range(480,600,6),linewidth=0.7,colors='Gray',cmap=None,label='500 hPa height')
    ax.clabel(CS2,fmt='%4d',fontsize=7)
    #xo, yo = m(lns, lts)
    matplotlib.rcParams['lines.linewidth'] = .25
    
    CS3 = m.contour(x,y,precip,levels=[1,2,5,10,25,50],linewidth=0.1,colors='Black',cmap=None)
    CS4 = m.contourf(x,y,precip,levels=[1,2,5,10,25,50],colors=colorst,cmap=None,extend='max')
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    
    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)

    cax = fig1.add_axes([0.6,0.55,0.3,0.015])
    cbar = fig1.colorbar(CS4,extend='neither', \
       orientation='horizontal',cax=cax,drawedges=True,ticks=[1,2,5,10,25,50],format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('Ensemble-mean precipitation amount (mm)',fontsize=7)
    
    # ---- fourth plot, CAPE, CIN, and 0-1 km shear
    
    
    title = '(d) Daily Avg sbCAPE, CIN & 0-6 km Bulk Shear\n{0} to {1}\n'.format(start_fdattim_str,end_fdattim_str)+\
            'Initialization time = {0}'.format(current_fdattim_str)

    ax = fig1.add_axes([0.52,.08,0.46,0.33])
    ax.set_title(title,fontsize=10)
    m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
    matplotlib.rcParams['lines.linewidth'] = .7
    CS1 = m.contour(x,y,cape,levels=[50,125,250,500,1000,1500,2000,2500,3000,4000],linewidth=.7,colors='Black',cmap=None,label='CAPE')
    CS4 = m.contour(x,y,shear,levels=range(10,80,5),linewidth=.7,colors='Black',cmap=None,label='Shear')
    ax.clabel(CS4,fmt='%3d',fontsize=7)
    CS2 = m.contourf(x,y,cape,levels=[50,125,250,500,1000,1500,2000,2500,3000,4000],colors=colorst,cmap=None,label='CAPE',extend="max")
    #ax.clabel(CS2,fmt='%2d',fontsize=9)
    CS3 = m.contour(x,y,cin*-1,levels=[50,100,150,250,500],\
                            linewidth=0.5,colors='Gray',cmap=None,label='CIN')
    ax.clabel(CS3,fmt='%3d',fontsize=7)
    m.drawcoastlines(linewidth=.5)
    m.drawstates(linewidth=.5)
    m.drawcountries(linewidth=.5)
    
    ax.text(1.0,-0.040,'NOAA/ESRL/PSD and F. Alvarez (St. Louis U.)',ha='right',transform=ax.transAxes,fontsize=7)

    cax = fig1.add_axes([0.57,0.045,0.36,0.015])
    cbar = fig1.colorbar(CS2,extend='max', \
       orientation='horizontal',cax=cax,drawedges=True,\
       ticks=[50,125,250,500,1000,1500,2000,2500,3000,4000],format='%g')
    cbar.ax.tick_params(labelsize=7)
    cax.set_xlabel('sbCAPE',fontsize=7)

    plt.savefig('{0}tornado_fcst_4panel_{1}_day{2}.png'.format(image_dir,forecastDate.strftime("%Y%m%d"),leadtime),dpi=200)
    plt.clf()
    plt.close()
  
def plot_probs_timespan(forecastDate,leadtime,probs,image_dir,maxlat,minlat,maxlon,minlon,fcst_minlat,\
        fcst_maxlat,fcst_minlon,fcst_maxlon,**kwargs):
    """
        plot_probs(forecastDate,leadtime,probs,image_dir,**kwargs)
  
      A function that plots out the forecast probabilities trained on a larger ROI on a single panel.
        
        Required arguments:
            forecastDate - some datetime object, the date of the forecast initialization
            leadtime - lead time for the forecast, from 12z-12z (so leadtime=1 is 12z-36z)        
            probs - some 2-D numpy array of forecast probabilities [0,100]
            image_dir - path to which the images are saved
            minlat,maxlat,minlon,maxlon - forecast domain
        Optional Arguments:
            radius - some radius for the tornado probabilities (default = 240km)  
    """
    radius = kwargs.get('radius',240)     


    lns = np.arange(minlon,maxlon+1)
    lts = np.arange(minlat,maxlat+1)
    lns,lts = np.meshgrid(lns,lts)
    
    fdattim_str = (forecastDate+timedelta(days=int(leadtime)-1,hours=12)).strftime("%Y-%m-%d 12 UTC")
    fdattim2_str = (forecastDate+timedelta(days=int(leadtime),hours=12)).strftime("%Y-%m-%d 12 UTC")

    fig1 = plt.figure(figsize=(14.5,6.5))
    title = 'Tornado probabilities (EF1+), {0} km ROI. Valid: {1} - {2}'.format(radius,fdattim_str,fdattim2_str)
        #'Initialization time = {0} UTC'.format(forecastDate)
    fig1.suptitle(title,fontsize = 14)
    colorst = ['BurlyWood','LightGray','LightSkyBlue','LightGreen','Green','Gold',\
               'Orange','Red','Plum','Orchid']
    
    if radius == 80:
        clevs=[1,2,3,4,5,7,10,12,15,20,25]
    elif radius == 160:
        clevs=[2,5,7,10,12,15,20,25,30,35,40]
    else:
        clevs=[5,10,15,20,25,30,35,40,45,50,60]
    
    # ---- make first plot, the probability forecast
    
    for idx in reversed(xrange(1,11,1)):
        fdattim_str = (forecastDate-timedelta(days=int(idx - leadtime))).strftime("%Y-%m-%d")
        ax_title = 'Lead-time Day +{0}'.format(idx)
        ax = fig1.add_subplot(2,5,(11-idx))
        ax.set_title(ax_title,fontsize=9)
        m = Basemap(llcrnrlon=fcst_minlon,llcrnrlat=fcst_minlat,urcrnrlon=fcst_maxlon,urcrnrlat=fcst_maxlat,projection='mill',resolution='l')
        #m.fillcontinents(color='beige')
        x, y = m(lns, lts)
        if (probs[10-idx,...].reshape(-1).max() < 0.):
            CS3 = m.contourf(x,y,probs[10-idx,...],levels=[-1,0],colors=['DarkGray','DarkGray'],cmap=None,extend='min',alpha=.75)
        CS1 = m.contour(x,y,probs[10-idx,...],levels=clevs,linewidth=0.3,colors='k',cmap=None)
        CS2 = m.contourf(x,y,probs[10-idx,...],levels=clevs,colors=colorst,cmap=None,extend='max')
        m.drawcoastlines(linewidth=.5)
        m.drawstates(linewidth=.5)
        m.drawcountries(linewidth=.5)   
        
    cax = fig1.add_axes([0.25,0.06,0.5,0.03])
    cbar = fig1.colorbar(CS2,extend='neither', \
    orientation='horizontal',cax=cax,drawedges=True,ticks=clevs,format='%g')
    cbar.ax.tick_params(labelsize=8)
    cax.set_xlabel('Forecast tornado probability (%)',fontsize=9)

    plt.tight_layout()
    plt.savefig('{0}tornado_10day_panel_{1}_day{2}.png'.format(image_dir,forecastDate.strftime("%Y%m%d"),leadtime),dpi=200)
    plt.clf()
    plt.close()
  
