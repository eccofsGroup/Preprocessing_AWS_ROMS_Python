import numpy as np
import pandas as pd
import datetime 
from datetime import date,timedelta,datetime
import cartopy.crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import xarray as xr
import dask
import netCDF4 as nc
import time,os
import requests


baseurl='https://www.star.nesdis.noaa.gov/thredds/dodsC/gridG19ABINRTL3CWW00/'
suffix='-NOAA-L3C_GHRSST-SSTsubskin-ABI_G19-ACSPO_V3.00-v02.1-fv01.0.nc'
nday=3
outdirec='/home/om/cron/ECCOFS_OBS/PO.DAAC/data/raw/GOES19_L3/'
filepart='_NOAA-L3C_GHRSST-SSTsubskin-ABI_G19-ACSPO_V3.00'



def main(fconfig):
    print('getting GOES19 SST')
    encoding={}
    baseurl=fconfig['obs']['sst']['gbaseurl']
    suffix=fconfig['obs']['sst']['gsuffix']
    nday=fconfig['obs']['sst']['gnday']
    outdirec=fconfig['obs']['sst']['goutdirec']
    filepart=fconfig['obs']['sst']['gfilepart']
    dhour=range(0, 24)    
    start_date = date.today()-timedelta(days=1)
    date_times = [start_date  - timedelta(days=i) for i in range(nday)] 
    for day in date_times:
        yd=day.timetuple().tm_yday
        url=f'{baseurl}{day.year}/{yd:03d}'
        print(day)
        print(f"Today is day {yd} of the year.")
        sdate=day.strftime('%Y%m%d')
        ofile=outdirec+sdate+filepart
        urllist=[]
        for h in dhour:
            dh=timedelta(hours=h)
            ctime= datetime.combine(day, datetime.min.time())+dh
            tstamp=ctime.strftime('%Y%m%d%H%M%S')
            url=f'{baseurl}{day.year}/{yd:03d}/{tstamp}{suffix}'
            urllist.append(url)
        
        #print(urllist)
        print('Opening Datset')
        ds=xr.open_mfdataset(urllist,decode_timedelta=False)
        ds=ds.drop_vars('sst_dtime')
        print('subsetting')
        ds=ds.sel(lat=slice(52,0),lon=slice(-100,-36))
    
    
        compression_settings = {
            var: {"zlib": True, "complevel": 2}  # Enable zlib compression with level 5
            for var in ds.data_vars
        }
        
        encoding.update(compression_settings)
    
        print('saving dataset')
        ds.to_netcdf(ofile,encoding=encoding)
    #    ds.to_netcdf('test.nc')
            
