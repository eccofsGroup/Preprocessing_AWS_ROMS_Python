#!/usr/bin/env python
#
import os,sys,fnmatch,glob
from shutil import copyfile
from datetime import datetime,timedelta,date
import subprocess,shutil,os,glob,cartopy
import xarray as xr
import copernicusmarine as copernicus_marine
import pandas as pd
import matplotlib.pyplot as plt
outdir='/home/om/cron/ECCOFS_OBS/MERCATOR/data/raw/'

tmpdir='/home/om/cron/ECCOFS_OBS/MERCATOR/tmp/'
days=14

gfile='/home/om/roms/eccofs/grid_eccofs_6km_01.nc'

proj = cartopy.crs.Mercator(central_longitude=-74)
pc = cartopy.crs.PlateCarree()


latmin=0.0
latmax=52.0
lonmin=-100.0
lonmax=-36.0
zmin=-10.0;
zmax=6000.0
variables=['so' ,'thetao','uo','vo','zos']
tmpfiles=['TMP_so.nc','TMP_thetao.nc','TMP_uv.nc','TMP_zos.nc']



#Choose the appropriate dataset. 
ids='cmems_mod_glo_phy_my_0.083deg_P1D-m' #REANALYIS PRODUCT through 2021/6/30
ids='cmems_mod_glo_phy_myint_0.083deg_P1D-m' #interim REANALYIS PRODUCT 2021/7/1 -2024/10


#Set the time range for download

tstart=date(2021,7,1)
tend=date(2023,1,1)

days=pd.date_range(tstart,tend)

#mss : mss dtu15

def main(argv):
    print('MAIN')
    
    #ds=xr.open_dataset(gfile2)
    #plot_2D_LL(ds['lon_rho'],ds['lat_rho'],ds['h'])
    
    
    for day in days:
        print(day)


        day2=day+timedelta(days=0.5)
        ofilename = day.strftime(outdir+'mercator_doppio_%Y_%m_%d.nc')    
        
        copernicus_marine.subset(
            dataset_id = ids,
            variables = variables,
            minimum_longitude = lonmin,
            maximum_longitude = lonmax,
            minimum_latitude = latmin,
            maximum_latitude = latmax,
            minimum_depth = zmin,
            maximum_depth = zmax,
            start_datetime = day.strftime('%Y-%m-%d %H:%M:%S'),
            end_datetime = day2.strftime('%Y-%m-%d %H:%M:%S'),
            output_filename = ofilename,
            force_download=True,
            overwrite_output_data =True,
    #        dataset_version="202406",
            output_directory = tmpdir,
            )
    

            
def plot_2D_LL(lon,lat,var):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1,projection=proj)
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines(resolution='50m')

    ax.set_extent([lonmin, lonmax, latmin,latmax])
    
    
    pcid=ax.pcolormesh(lon,lat,var,transform=pc)
    colorbar = plt.colorbar(pcid)
    plt.show()
    
def plot_points_LL(lon,lat,var):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1,projection=proj)
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines(resolution='50m')

    ax.set_extent([lonmin, lonmax, latmin,latmax])
        
        
    pcid=ax.scatter(lon,lat,c=var,s=20,transform=pc,marker='o',vmin=0,vmax=30)
    colorbar = plt.colorbar(pcid)
    plt.show()
if __name__ == "__main__":

    
    print('RUNNING')
    main(sys.argv)
    # for f in glob.glob(tmpdir+'*.nc'):
    #     os.remove(f)
    # for f in glob.glob(tmpdir+'*.tmp'):
    #     os.remove(f)
