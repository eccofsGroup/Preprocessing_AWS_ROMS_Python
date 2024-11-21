"""
GFS_process.py 

Subset, download, and process GFS data from AWS bucket s3://hrrrzarr. 
Variables downloaded are ROMS meteorological forcing variables. Uses kerchunk
to simulate cloud optimized files. 

Created by Elias Hunter, hunter@marine.rutgers.edu, 1/31/2023 
"""

import zarr 
import xarray
import s3fs
import numpy as np
import pandas as pd
import datetime 
import cartopy.crs as ccrs
import cartopy
import xarray as xr
import fsspec
import dask
from kerchunk.grib2 import scan_grib
from kerchunk.combine import MultiZarrToZarr
import netCDF4 as nc
import time,os
from dask.distributed import Client
from glob import glob
import ujson
import smtplib


#Parameters
print('Setting up necessary parmeters')
baseurl = "s3://hrrrzarr"

###############################################################
#Things to edit. 
###############################################################
odir='/home/ubuntu/data/gfs'
json_dir = '/home/ubuntu/python/HRRR/jsons/'
###############################################################


#ECCOFS BOUNDS
mnlon=-116.0
mxlon=-24.0 
mnlat=-6.0 
mxlat=55.0
dlon=1
dlat=1

today = datetime.datetime.utcnow().strftime('%Y%m%d')
tdelta = datetime.timedelta(1)
ttmp= datetime.datetime.utcnow()-tdelta
yesterday=ttmp.strftime('%Y%m%d')

oVariables = {
    'Tair': ("t2m"), 
    'Qair': ("r2"), 
    'Pair': ("prmsl"),
    'rain': ("prate"),
    'wind': ("UV"),
    'swrad': ("SWRAD"), 
    'lwrad': ("LWRAD"),
    'lwrad_down': ("dlwrf"), 
    'swdown': ("dswrf"), 
    'swup': ("uswrf"),
    'lwdown': ("dlwrf"), 
    'lwup': ("ulwrf"),
    'Uwind': ("u10"),
    'Vwind': ("v10")
    
}
    


projection = ccrs.LambertConformal(central_longitude=262.5, 
                                   central_latitude=38.5, 
                                   standard_parallels=(38.5, 38.5),
                                    globe=ccrs.Globe(semimajor_axis=6371229,
                                                     semiminor_axis=6371229))
                                                     
                                                     


rtime=datetime.datetime(2006,1,1,0,0,0,0) #reference time
rformat=rtime.strftime('%Y-%m-%d 00:00:00')
Fillvalue=np.float64(9.999e20)

units = {
    'Tair': ('C'), 
    'Qair': ('percentage'), 
    'Pair': ('millibar'),
    'rain': ('kilogram meter-2 second-1'),
    'Uwind': ('meter second-1'),
    'Vwind': ('meter second-1'),
    'swrad': ('watt meter-2'), 
    'lwrad': ('watt meter-2'),
    'lwrad_down': ('watt meter-2'), 

    }

longname = {
    'Tair': ('temperature 2 m'), 
    'Qair': ('Relative_humidity_height_above_ground 2 m'), 
    'Pair': ('Pressure_reduced_to_MSL_msl'),
    'rain': ("PRATE",'surface'),
    'Uwind': ('Eastward Wind'),
    'Vwind': ('Northward wind'),
    'swrad': ('difference Downward_ minus Upward_Short-Wave_Radiation_Flux_surface'), 
    'lwrad': ('difference Downward_ minus Upward_Long-Wave_Radiation_Flux_surface'),
    'lwrad_down': ('Doward Longwave Raditation Flux'), 

    }
    
ROTCON_P = np.float64(0.622515) 
LON_XX_P = np.float64(-97.5) 
LAT_TAN_P = np.float64(38.5) 

nfcst=121
afilter={'cfVarName':'none'}  
so = {"anon": True}
sc1=np.tile([1.0,2.0,3.0,4.0,5.0,6.0],20)
sc1=np.append(sc1,1.0)

vfilter={'cfVarName':['prmsl','t2m','r2','u10','v10','prate'],"stepType":'instant'}  
rfilter={'shortName':['dswrf','dlwrf','uswrf','ulwrf'],'typeOfLevel':'surface'}   




############################################################################
#Main Program
############################################################################

def main():
    
    st = time.time()
    files = glob(json_dir+'*.json' )
    for f in files: 
        os.remove(f)

    make_grib_jsons()
    make_rad_jsons()
    for key in oVariables.keys():
        print('===========================================')
        #delete any existing reference jsons
  
        match key:
            case 'Tair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key)
                #Convert temperature to Celcius 
                ds=ds-273.15
                write_output(ds,key)                
            case 'Qair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key)
                
                write_output(ds,key)   
            case 'Pair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key)
                ds=ds*0.01
                write_output(ds,key)   
            case 'wind':
                print(f'Processing Variable:{key}')
                uwind=get_3d_variable('Uwind')
            
                vwind=get_3d_variable('Vwind')
             
                
                write_output(uwind,'Uwind')
                write_output(vwind,'Vwind')  
                
            case 'rain':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key)
                write_output(ds,key)  
            case 'swrad':
                print(f'Processing Variable:{key}') 
                dsup=get_RAD_variable('swup')
                dsdown=get_RAD_variable('swdown')
                ds=dsdown.load()-dsup.load()
    
                write_output(ds,key) 
            case 'lwrad':
                print(f'Processing Variable:{key}') 
                dsup=get_RAD_variable('lwup')
                dsdown=get_RAD_variable('lwdown')
                ds=dsdown.load()-dsup.load()


                write_output(ds,key) 
            case 'lwrad_down':
                print(f'Processing Variable:{key}') 
                
                dsdown=get_RAD_variable('lwdown')


                write_output(dsdown,key)   
            

    et = time.time()    
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    
    
def make_grib_jsons():
    
    #Uses kerchubnk to create non-radiation variable metadat jsons. 
    print('Making non-radiation JSONS')
    fs_read = fsspec.filesystem('s3', anon=True, skip_instance_cache=True)
    urls =  ['s3://' + f for f in fs_read.glob(f's3://noaa-gfs-bdp-pds/gfs.{today}/00/atmos/gfs.t00z.pgrb2.0p25.f*')]
    urls = [f for f in urls if not f.endswith('.idx')]
    fs_write = fsspec.filesystem('')
    urls=urls[0:nfcst]
    
  
    
    for ifi in urls:
        print(ifi)
        
        
        out = scan_grib(ifi, common=None,storage_options=so, inline_threshold=200,filter=vfilter) 
        
        for ind,vname in enumerate(vfilter['cfVarName']):
            fcst_t=ifi.split('.')[-1:]
            json_filename=f'{json_dir}GFS_GRIB_{vname}_{fcst_t[0]}.json'
            print(json_filename)
            with open(json_filename, "wb") as f:
                f.write(ujson.dumps(out[ind]).encode())
                
def make_rad_jsons():
    #Uses kerchubnk to create radiation variable metadat jsons. 
    print('Making radiation JSONS')
    
    fs_read = fsspec.filesystem('s3', anon=True, skip_instance_cache=True)
    urls =  ['s3://' + f for f in fs_read.glob(f's3://noaa-gfs-bdp-pds/gfs.{today}/00/atmos/gfs.t00z.pgrb2.0p25.f*')]
    urls = [f for f in urls if not f.endswith('.idx')]
    tmpurl=(f's3://noaa-gfs-bdp-pds/gfs.{yesterday}/18/atmos/gfs.t18z.pgrb2.0p25.f006')
    urls=urls[0:nfcst+1]
    urls[0]=tmpurl
    
  
    
    for ifi in urls:
        print(ifi)
        
        
        out = scan_grib(ifi, common=None,storage_options=so, inline_threshold=200,filter=rfilter) 
        
        for ind,vname in enumerate(rfilter['shortName']):
            fcst_t=ifi.split('.')[-1:]
             
            if '18z' in ifi:
                json_filename=f'{json_dir}GFS_GRIB_{vname}_f000.json'
            else:
                json_filename=f'{json_dir}GFS_GRIB_{vname}_{fcst_t[0]}.json'
            
            print(json_filename)
            with open(json_filename, "wb") as f:
                f.write(ujson.dumps(out[ind]).encode())
                

def write_output(ds,key):
    #Write the netcdf file
    ofile=f'{odir}/GFS_GRIB_AWS_{key}_{today}.nc'
    

#    ds.to_netcdf(ofile)
    ncid = nc.Dataset(ofile, 'w', format='NETCDF4')
    time = ncid.createDimension('time', None)
    lat = ncid.createDimension('lat', ds.sizes['latitude'])
    lon = ncid.createDimension('lon', ds.sizes['longitude'])
    
    
    times = ncid.createVariable('time', 'f8', ('time',))
    times.units = f'days since {rformat}'
    times.long_name='time'
    timeout=(ds.time.values-np.datetime64(rtime)) / np.timedelta64(1, 'D')
    times[:]=timeout
    
    lats = ncid.createVariable('lat', 'f8', ('lat',))
    lats.units='degrees_north'
    lats.long_name='latitude'
    lats.minimum=np.min(ds.latitude.values)
    lats.maximum=np.max(ds.latitude.values)

    lats[:]=ds.latitude.values
    
    lons = ncid.createVariable('lon', 'f8', ('lon',))
    lons.units='degrees_east'
    lons.long_name='longitude'
    lons.minimum=np.min(ds.longitude.values)
    lons.maximum=np.max(ds.longitude.values)
    lons[:]=ds.longitude.values
    
    
    
    
    value = ncid.createVariable(key, 'f8', ('time', 'lat', 'lon',),fill_value=Fillvalue)
    value.units = units[key]
    value.time = 'time'
    value.coordinates = 'lon lat'
    value.long_name =  longname[key]
    value.missing_value =  Fillvalue
    tmp=ds.values
    tmp[tmp==-999]=np.nan
    
    value[:,:,:]=ds.values
    
    extime = datetime.datetime.utcnow().strftime('%Y:%m:%d %H:%M')
    ncid.title='NCEP Global Forecast System (GFS)  model 0.25 degree. https://www.emc.ncep.noaa.gov/emc/pages/numerical_forecast_systems/gfs.php'
    ncid.source='NCEP Global Forecast System (GFS) Data Archive: AWS Open Data Program. https://registry.opendata.aws/noaa-gfs-bdp-pds/ '
    ncid.creator='Elias Hunter (hunter@marine.rutgers.edu)'
    ncid.history=f'Created by Eli Hunter on {extime} using AWS serivces. {__file__} '
    ncid.close()
    print(f'Data written to {ofile}')
    
def gen_json_grib(u,afilter):

   # 's3://noaa-gfs-bdp-pds/gfs.20230127/00/atmos/gfs.t00z.pgrb2.0p25.f000'

    name = u.split('/')[-1:][0]
    outfname = f'{json_dir}{name}.json'
    out = scan_grib(u, common=None,storage_options=so, inline_threshold=200,filter=afilter) 
    
    with open(outfname, "wb") as f:
        f.write(ujson.dumps(out[0]).encode())
        
    
def get_3d_variable(var):        
# get_3d_variable: gets the non-radiation GFS variables 
  # st = time.time()
    print(f'Accessing data for {var}')
    

    reference_jsons = glob(json_dir+f'*{oVariables[var]}*.json')
    reference_jsons= sorted(reference_jsons)
        
    mzz = MultiZarrToZarr(reference_jsons , concat_dims=['valid_time'],
             identical_dims=['latitude', 'longitude'])
    
    d = mzz.translate()
    fs = fsspec.filesystem("reference", fo=d, remote_protocol='s3', remote_options={'anon':True})
    m = fs.get_mapper("")
    ds_fcst = xr.open_dataset(m, engine="zarr", backend_kwargs=dict(consolidated=False), 
                      chunks={'valid_time':1})
    ds_fcst.coords['longitude'] = (ds_fcst.coords['longitude'] + 180) % 360 - 180

    ds_fcst=ds_fcst.sortby('latitude')
    ds_fcst=ds_fcst.sortby('longitude')
    ds_fcst=ds_fcst.sel(latitude=slice(mnlat,mxlat),longitude=slice(mnlon,mxlon))

    ds_fcst=ds_fcst.drop(['time','step'])
    ds_fcst=ds_fcst=ds_fcst.rename({'valid_time':'time'})

    ds_fcst=ds_fcst[oVariables[var]]

  

  #  et = time.time()    
  #  elapsed_time = et - st
  #  print('Execution time:', elapsed_time, 'seconds')

    

    return ds_fcst
    
def get_RAD_variable(var):        
# get_RAD_variables: gets the radiation GFS variables 
    print(f'Accessing data for {var}')
 
    
    print(oVariables[var])
    reference_jsons = glob(json_dir+f'*{oVariables[var]}*.json')
    reference_jsons= sorted(reference_jsons)
   
    mzz = MultiZarrToZarr(reference_jsons , concat_dims=['valid_time'],
             identical_dims=['latitude', 'longitude'])
    
    d = mzz.translate()
    fs = fsspec.filesystem("reference", fo=d, remote_protocol='s3', remote_options={'anon':True})
    m = fs.get_mapper("")
    ds_fcst = xr.open_dataset(m, engine="zarr", backend_kwargs=dict(consolidated=False), 
                      chunks={'valid_time':1})

    ds_fcst.coords['longitude'] = (ds_fcst.coords['longitude'] + 180) % 360 - 180

    ds_fcst=ds_fcst.sortby('latitude')
    ds_fcst=ds_fcst.sortby('longitude')
    ds_fcst=ds_fcst.sel(latitude=slice(mnlat,mxlat),longitude=slice(mnlon,mxlon))
    ds_fcst=ds_fcst.drop(['time','step'])
    ds_fcst=ds_fcst=ds_fcst.rename({'valid_time':'time'})

    ds_fcst=ds_fcst[oVariables[var]]


    vname=oVariables[var]
  #  print(vname)
    ds_fcst2=ds_fcst.copy(deep=True)
    
    
    for ind in range(1,nfcst+1):

        c1=sc1[ind-1]
        if c1==1:
            c2=0
        else:
            c2=sc1[ind-2]   
        

    
        
        tmp1=ds_fcst.isel(time=ind)
        tmp2=ds_fcst.isel(time=ind-1)
    
        tmp=tmp1*c1-tmp2*c2      

        ds_fcst2[ind,:,:]=tmp.values

    
    #This were we get the 5 hour forcast to calculate 6 hour forecast. From 18Z forecast

    afilter2={'shortName':oVariables[var]}
    tmpurl5=f's3://noaa-gfs-bdp-pds/gfs.{yesterday}/18/atmos/gfs.t18z.pgrb2.0p25.f005'
    out = scan_grib(tmpurl5, common=None,storage_options=so, inline_threshold=200,filter=afilter2) 
    fs = fsspec.filesystem("reference", fo=out[0], remote_protocol='s3', remote_options={'anon':True})
    m = fs.get_mapper("")
    ds_fcst5 = xr.open_dataset(m, engine="zarr", backend_kwargs=dict(consolidated=False), 
                      chunks={'valid_time':1})
    ds_fcst5.coords['longitude'] = (ds_fcst5.coords['longitude'] + 180) % 360 - 180

    ds_fcst5=ds_fcst5.sortby('latitude')
    ds_fcst5=ds_fcst5.sortby('longitude')
    ds_fcst5=ds_fcst5.sel(latitude=slice(mnlat,mxlat),longitude=slice(mnlon,mxlon))

    ds_fcst5=ds_fcst5.drop(['time','step'])
    ds_fcst5=ds_fcst5.rename({'valid_time':'time'})
    ds_fcst5=np.squeeze(ds_fcst5)
    tmp1=ds_fcst.isel(time=0)
    tmp2=ds_fcst5.to_array()
    tmp=tmp1*6.0-tmp2*5.0 
    tmp=np.squeeze(tmp)

    ds_fcst2[0,:,:]=tmp.values
    
    stime=ds_fcst2.time.values[0]
    etime=ds_fcst2.time.values[-2]
    ds_fcst2['time']=ds_fcst2.time-np.timedelta64(30,'m')
    
    Ids2=ds_fcst2.interp(time=pd.date_range(stime,etime, freq='1H'))
    
    return Ids2
    
    

if __name__ == "__main__":
    print('Running')
    client = Client()
    main()
 