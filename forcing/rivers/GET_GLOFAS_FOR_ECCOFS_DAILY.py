#!/usr/bin/env python3
"""
GET_GLOFAS_FOR_ECCOFS_DAILY.py

Download GLOFAS forecasts for the last 3 days and create a ROMS rivers file 
for each.

Author Elias Hunter hunter@marine.rutgers.edu created 2024/05/02
"""
import os,glob
import numpy as np
import xarray as xr
import cartopy 
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import datetime,time
from datetime import timedelta,date
import romsutil as rutil
from scipy.io import loadmat
import netCDF4 as nc
import cdsapi

"""
Index file created by John Wilkin. The I,J variables from this file are 
    locations on the GLOFAS frid for dhsiarge near the coastline. It also     
    Xpostion and Eposition for the ROMS focing file. 
Some notes from Johns m-file:
    % Updating the GLOFAS rivers
% Go to: https://cds.climate.copernicus.eu/cdsapp#!/dataset/cems-glofas-historical?tab=form
%
% THIS MIGHT HAVE CHANGED TO
% https://ewds.climate.copernicus.eu/datasets/efas-historical?tab=download
%
% Select: LISFLOOD, Consolidated
% Year  - select
% Month - select all
% Day   - select all
% Subregion - W E S N 92.5 145.0 -20.0 25.0    <<<<<<<<<<<<< MINTIE
% Subregion - W E S N -105.0 -35.0 -5.0 55.0   <<<<<<<<<<<<< ECCOFS
% Select grib
% Login jwilkin@rutgers.edu nf0yawmd
% Submit request and wait.
% Download the netcdf file. It is zipped and will unpack into the same directory with the
% name data.nc (so it will overwrite other years if you have not renamed them)
% Go to /Users/wilkin/Dropbox/_roms-db/eccofs/rivers
% See script add_another_glofas_rivers_year.m

"""
rIJfile='/home/om/cron/ECCOFS_OBS/RIVERS/work/glofasV4_river_points_for_eccofs.mat'
rdata=loadmat(rIJfile)
lat_indices = xr.DataArray(rdata['J'][:][0], dims="points")-2    
lon_indices = xr.DataArray(rdata['I'][:][0], dims="points")-2

#Only for plotting
proj = cartopy.crs.Mercator(central_longitude=-74)
pc = cartopy.crs.PlateCarree()
latmin=-3.0
latmax=53.0
lonmin=-100.0
lonmax=-38.0

#NetCDF parameters
rpre='/home/om/cron/ECCOFS_OBS/RIVERS/data/processed/daily/'
fname='ECCOFS_rivers_file'
reftime=datetime.datetime(2011,1,1)
tunits=reftime.strftime('days since %Y-%m-%d')
rtime=np.datetime64(reftime.date())
reftime=datetime.datetime(2011,1,1)
rlist=["Naming the first several major rivers only: \n",
                        "(1) Manicouagen \n",
                        "(2) Saguenay \n",
                        "(3) St Lawrence \n",
                        "(4) St John ME \n",
                        "(5) Susquehanna \n",
                        "(6) Mobile \n",
                        "(7) Mississippi \n",
                        "(8) Papaloapan \n",
                        "(9) Usumacinta \n",
                        "(10) Colorado CR \n",
                        "(11) Atrato \n",
                        "(12) Magdalena 1 \n",
                        "(13) Orinoco \n",
                        "(14) Nétagamiou \n",
                        "(15) Connecticut \n",
                        "(16) Hudson \n",
                        "(17) Delaware \n",
                        "(18) St Johns FL \n",
                        "(19) Apalachicola \n",
                        "(20) Rio Grande \n",
                        "(21) Pánuco \n",
                        "(22) Coatzacoalcos \n",
                        "(23) Rio Ulúa \n",
                        "(24) Rio Coco \n",
                        "(25) Rio Kurinwás \n",
                        "(26) Escondido \n",
                        "(27) Rio Sinu \n",
                        "(28) Magdalena 2 \n",
                        "(29) Catatumbou \n",
                        "(30) Rio Barima "]
                        
                        
                        
routdata = {
    "nriver": len(lon_indices),
    "nz": 50,
    "type": 'ROMS FORCING file',
    "title":'ECCOFS River Forcing from GloFAS',
    "history":['Created by create_roms_rivers_file  on '+ datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')+" (based on J. Wilkin matlab scripts)"],
    "gridfile": '/home/om/cron/ECCOFS_OBS/MERCATOR/work/grid_eccofs_3km_08_b7.nc',
    "source": ["Daily discharge from GloFAS Global Flood Awareness Program globalfloods.eu accessed from Copernicus"+
                        "  https://ewds.climate.copernicus.eu/datasets/cems-glofas-historical?tab=overview"+
                        "  Source point locations identified using driver_glofas_rivers.m John Wilkin 2023-05-17"],
    "citation":"Harrigan, S., Zsoter, E., Alfieri, L., Prudhomme, C., Salamon, P., Wetterhall, F., Barnard, C., Cloke, H. and Pappenberger, F., 2020. GloFAS-ERA5 operational global river discharge reanalysis 1979present. Earth System Science Data, 12(3), pp.2043-2060., https://doi.org/10.5194/essd-12-2043-2020",   
    "rivers": rlist,
    "river": 1,
    "river_Xposition": 1,
    "river_Eposition": 1,
    "river_Direction": 1,
    "river_lat": 1,
    "river_lon": 1,
    "river_Vshape": 1,
    "river_salt": 1,
    "river_sign": 1,
    "river_temp": 1,
    "river_time": 1,
    "tunits":tunits,
    "river_transport": 1,
    "river_Vshape": 1,
}




#Choose file format
ncflag=True
if ncflag:
    indir='/home/om/cron/ECCOFS_OBS/RIVERS/data/raw/daily/'
    inname='GLOFAS_GRIDDED_FORECAST'
    inext='.nc'
else:
    indir='/home/om/cron/ECCOFS_OBS/RIVERS/data/raw/daily/'
    inname='GLOFAS_GRIDDED_FORECAST'
    inext='.grb'
    
    
plotFLAG=False    
#
nday=3
start_date = datetime.datetime.now()
date_times = [start_date  - timedelta(days=i) for i in range(nday)]     
    
############################################################################
#Main Program
############################################################################

def main():

    for iday in date_times:
        idate=iday.strftime('%Y%m%d')
        
        rfile=f'{rpre}{fname}_{idate}.nc'
        ########################################################################
        print(f'Creating Rivers file: {rfile}')
        print('---------------------------')
        ########################################################################
        rutil.create_roms_rivers_file(rfile,routdata)
    
        infile=f'{indir}{inname}_{idate}{inext}'

        download_GLOFAS_forcast_daily(iday,infile)
        print('**********************************')
        if ncflag:
            ds=xr.open_dataset(infile)
            ds=ds.rename({'valid_time':'time'})
        else:
            ds=xr.open_dataset(infile,engine='cfgrib')
        start_time=time.time()
        ds.load()
        end_time=time.time()
        elapsed_time = end_time - start_time
        print(f"input file load time: {elapsed_time} seconds")
    

        dssubset = ds.isel(latitude=lat_indices, longitude=lon_indices)

    
        if plotFLAG:
        
            dsdis=np.log10(ds['dis24'].mean(dim=['forecast_period']))
            dsdis=dsdis.squeeze('forecast_reference_time')
            #print(dsdis)
            plot_2D_LL(dsdis['longitude'],dsdis['latitude'],dsdis,dssubset)
    
        start_time=time.time()
        dis=dssubset.dis24.values
        end_time=time.time()
        elapsed_time = end_time - start_time
        print(f"Disharge variable extraction time: {elapsed_time} seconds")
    
        print(dssubset)
        nriver=dssubset.sizes['points']
        ntime=dssubset.sizes['forecast_period']
        nz=routdata['nz']
        inz=1/nz

        tmp=np.zeros((nz,nriver))+inz
        salt=np.zeros((ntime,nz,nriver))
        rdir=np.zeros((nriver,1))+2
    
        print('Writing River data to file '+rfile)
        ncid = nc.Dataset(rfile, "r+", format="NETCDF4")
    
        t=dssubset.time.values-np.timedelta64(12,'h')
        timeout=(t-np.datetime64(rtime)) / np.timedelta64(1, 'D')
        ncid.variables['river_time'][:]=timeout
    
        ncid.variables['river_Eposition'][:]=rdata['river_Eposition']
        ncid.variables['river_Xposition'][:]=rdata['river_Xposition']
        ncid.variables['river_lat'][:]=dssubset.latitude.values[:]
        ncid.variables['river_lon'][:]=dssubset.longitude.values[:]
        ncid.variables['river_Vshape'][:]=tmp
        ncid.variables['river_direction'][:]=rdir
        ncid.variables['river_salt'][:]=salt
        ncid.variables['river_transport'][:,:]=dis
        
        ncid.createVariable('I_index_glofas_subset','double', ('river'))
        ncid.variables['I_index_glofas_subset'].description = "I index to glofas*.nc files to extract discharge at near coastal sources" 
            
        ncid.createVariable('J_index_glofas_subset', 'double', ('river'))
        ncid.variables['J_index_glofas_subset'].description = "J index to glofas*.nc files to extract discharge at near coastal sources" 
        
        ncid.variables['I_index_glofas_subset'][:]=rdata['I'][:][0]
        ncid.variables['J_index_glofas_subset'][:]=rdata['J'][:][0]
        
        
    
    
        
        ncid.sync()
        ncid.close()
def download_GLOFAS_forcast_daily(iday,infile):
    #SOURCE:https://ewds.climate.copernicus.eu/datasets
    print(f'DOWNLOADING GLOFAS {iday} forecast to {infile}')

    syear=iday.strftime('%Y')
    smonth=iday.strftime('%m')
    sday=iday.strftime('%d')

    dataset = "cems-glofas-forecast"
    request = {
        "system_version": ["operational"],
        "hydrological_model": ["lisflood"],
        "product_type": ["control_forecast"],
        "variable": "river_discharge_in_the_last_24_hours",
        "year": [syear],
        "month": [smonth],
        "day": [sday],
        "leadtime_hour": [
            "24",
            "48",
            "72",
            "96",
            "120",
            "144",
            "168",
            "192"
            ],
        "data_format": "netcdf",
        "download_format": "unarchived",
        "area": [55, -105, -5, -35]
        }

    client = cdsapi.Client()
    client.retrieve(dataset, request,infile)



    
    
    
    
    
def plot_2D_LL(lon,lat,var,subset):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1,projection=proj)
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False

    ax.set_extent([lonmin, lonmax, latmin,latmax])
    
    
    pcid=ax.pcolormesh(lon,lat,var,transform=pc)
    colorbar = plt.colorbar(pcid)
    pcid=ax.scatter(subset['longitude'],subset['latitude'],c='r',s=40,transform=pc,marker='x')
    pcid=ax.scatter(rdata['longitude'],rdata['latitude'],c='r',s=20,transform=pc,marker='o')
    ax.coastlines(resolution='50m')
    ax.add_feature(cfeature.LAKES, edgecolor="black")
    plt.show()    

if __name__ == "__main__":
    print('Running GLOFAS to ROMS point souce forcing file')
    main()
    print('Finished GLOFAS to ROMS point souce forcing file')
    
    
    