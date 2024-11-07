"""
HRRR_process.py 

Subset, donwload, and process HRRR data from AWS bucket s3://hrrrzarr. 
Variables dowoaded are ROMS meteorological forcing variables. Radiatio varaibles a

Created by Elias Hunter, hunter@marine.rutgers.edu, 1/31/2023 
"""

import zarr 
import xarray
import s3fs
import numpy as np
import datetime 
import cartopy.crs as ccrs
import cartopy
import xarray as xr
import metpy
import xesmf as xe
import netCDF4 as nc
import time,os
from herbie import Herbie
import smtplib



#Parameters
print('Setting up necessary parmeters')
baseurl = "s3://hrrrzarr"

###############################################################
#Things to edit. 
###############################################################

regrid_coef_file='/home/ubuntu/python/HRRR/HRRR_REGRID.nc'
Herbieregrid_coef_file='/home/ubuntu/python/HRRR/HRRR_REGRID_HERBIE.nc'
odir='/home/ubuntu/data/hrrr'


###############################################################

#ECCOFS BOUNDS
mnlon=-116.0
mxlon=-24.0 
mnlat=-6.0 
mxlat=55.0
d_lon=0.025
d_lat=0.025

icut=10#remove last ten x values

today = datetime.datetime.utcnow().strftime('%Y%m%d')
#today = '20230321'
tdelta = datetime.timedelta(1)
ttmp= datetime.datetime.utcnow()-tdelta
yesterday=ttmp.strftime('%Y%m%d')

Htoday = datetime.datetime.utcnow().strftime('%Y-%m-%d')


oVariables = {
    'Tair': ("TMP",'2m_above_ground'), 
    'Qair': ("RH",'2m_above_ground'), 
    'Pair': ("MSLMA",'mean_sea_level'),
    'rain': ("PRATE",'surface'),
    'wind': ("UV",'10m_above_ground'),
    'swrad': ("SWRAD",'surface'), 
    'lwrad': ("LWRAD",'surface'),
    'lwrad_down': ("DLWRF",'surface'), 
    'swdown': ("DSWRF",'surface'), 
    'swup': ("USWRF",'surface'),
    'lwdown': ("DLWRF",'surface'), 
    'lwup': ("ULWRF",'surface'),
    'Uwind': ("UGRD",'10m_above_ground'),
    'Vwind': ("VGRD",'10m_above_ground')

    }
    


projection = ccrs.LambertConformal(central_longitude=262.5, 
                                   central_latitude=38.5, 
                                   standard_parallels=(38.5, 38.5),
                                    globe=ccrs.Globe(semimajor_axis=6371229,
                                                     semiminor_axis=6371229))
                                                     
                                                     
                                                     


rtime=datetime.datetime(2006,1,1,0,0,0,0)
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




############################################################################
#Main Program
############################################################################

def main():
    
    
    
    for key in oVariables.keys():
        print('===========================================')
        match key:
            case 'Tair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key,True)
             #   #Convert temperature to Celcius 
                ds=ds-273.15
                write_output(ds,key)
                del ds
            case 'Qair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key,True)
                
                write_output(ds,key)
                del ds
            case 'Pair':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key,True)
                ds=ds*0.01
                write_output(ds,key)
                del ds
            case 'wind':
                print(f'Processing Variable:{key}')
                uwind=get_3d_variable('Uwind',False)
             
                vwind=get_3d_variable('Vwind',False)
        
                #Rotation from https://rapidrefresh.noaa.gov/faq/HRRR.faq.html
               


                olon=uwind.lon
                angle2=ROTCON_P*(olon-LON_XX_P)*np.float64(0.017453)
                sinx2 = np.sin(angle2)
                cosx2 = np.cos(angle2)
                u= cosx2*uwind+sinx2*vwind
                v=-sinx2*uwind+cosx2*vwind
                
                
                varint=xe.util.grid_2d(mnlon, mxlon, d_lon, mnlat, mxlat, d_lat)

                if os.path.exists(regrid_coef_file):
                    print('Re-using weights')
                    regridder = xe.Regridder(u, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
                    u = regridder(u)
                else:
                    print('Generating weights')
                    regridder = xe.Regridder(u, varint, "bilinear",periodic=False,unmapped_to_nan=True)
                    regridder.to_netcdf(regrid_coef_file)
                    u = regridder(u)
            

                if os.path.exists(regrid_coef_file):
                    print('Re-using weights')
                    regridder = xe.Regridder(v, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
                    v = regridder(v)
                else:
                    print('Generating weights')
                    regridder = xe.Regridder(v, varint, "bilinear",periodic=False,unmapped_to_nan=True)
                    regridder.to_netcdf(regrid_coef_file)
                    v = regridder(v)

                u=u.transpose('time','y','x')
                v=v.transpose('time','y','x')
                write_output(u,'Uwind')
                write_output(v,'Vwind')  
                del u,v,uwind,vwind
            case 'rain':
                 print(f'Processing Variable:{key}') 
                 ds=get_3d_variable(key,True)
                 ds=get_PRATE_INIT(key,ds)
                 write_output(ds,key) 
                 del ds
            case 'swrad':
                print(f'Processing Variable:{key}') 
                dsup=get_SW_variable('swup')
                dsdown=get_SW_variable('swdown')
                ds=dsdown.load()-dsup.load()
                del dsup,dsdown
             #   print(dsdown.isel(x=10,y=10,time=18).values)
             #   print(dsup.isel(x=10,y=10,time=18).values)
             #   print(ds.isel(x=10,y=10,time=18).values)
   

                #print(ds.values)
                write_output(ds,key)
                del ds
            case 'lwrad':
                print(f'Processing Variable:{key}') 
                dsup=get_3d_variable('lwup',True)
                dsdown=get_3d_variable('lwdown',True)
                ds=dsdown.load()-dsup.load()
                
                del dsup,dsdown
     #          print(dsdown.isel(x=10,y=10,time=18).values)
     #           print(dsup.isel(x=10,y=10,time=18).values)
     #           print(ds.isel(x=10,y=10,time=18).values)
                write_output(ds,key) 
                del ds
            case 'lwrad_down':
                print(f'Processing Variable:{key}') 
                ds=get_3d_variable(key,True)
                write_output(ds,key)   
            
                del ds

    

def write_output(ds,key):
    
    #Write netcdf file. 
    ofile=f'{odir}/HRRR_ZARR_REGRIDED_{key}_{today}.nc'
    

#    ds.to_netcdf(ofile)
    ncid = nc.Dataset(ofile, 'w', format='NETCDF4')
    time = ncid.createDimension('time', None)
    lat = ncid.createDimension('lat', ds.sizes['y'])
    lon = ncid.createDimension('lon', ds.sizes['x'])
    
    
    times = ncid.createVariable('time', 'f8', ('time',))
    times.units = f'days since {rformat}'
    times.long_name='time'
    timeout=(ds.time.values-np.datetime64(rtime)) / np.timedelta64(1, 'D')
    times[:]=timeout
    
    lats = ncid.createVariable('lat', 'f8', ('lat',))
    lats.units='degrees_north'
    lats.long_name='latitude'
    lats.minimum=np.min(ds.lat.values)
    lats.maximum=np.max(ds.lat.values)
    lats.resolution=d_lat
    lats[:]=ds.lat.values[:,0]
    
    lons = ncid.createVariable('lon', 'f8', ('lon',))
    lons.units='degrees_east'
    lons.long_name='longitude'
    lons.minimum=np.min(ds.lon.values)
    lons.maximum=np.max(ds.lon.values)
    lons.resolution=d_lon
    lons[:]=ds.lon.values[0,:]
    
    
    
    
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
    ncid.title='NCEP High Resolution Rapid Refresh (HRRR) model 4km. https://rapidrefresh.noaa.gov/hrrr/'
    ncid.source='NOAA High-Resolution Rapid Refresh (HRRR) Data Archive: AWS Open Data Program.  MesoWest Group. https://mesowest.utah.edu/html/hrrr/'
    ncid.creator='Elias Hunter (hunter@marine.rutgers.edu)'
    ncid.history=f'Created by Eli Hunter on {extime} using AWS serivces. {__file__} '
    ncid.close()
    print(f'Data written to {ofile}')
    
def get_SW_variable(var):  
#Uses Herbie to get radiation variables. 
    st = time.time() 
    H = Herbie(
    today,
    model="hrrr",
    product="sfc",
    fxx=0,
    )
    ds = H.xarray(f'{oVariables[var][0]}:surface')

    for fxx in range(48):
        H = Herbie(
        today,
        model="hrrr",
        product="sfc",
        fxx=fxx+1,
        )
        print(fxx+1)
        tds = H.xarray(f'{oVariables[var][0]}:surface')

        ds = xr.concat([ds,tds], dim='time')
        
    print(ds.dims)
    ds=ds.isel(x=slice(0,-icut))
    print(ds.dims)
    
    varint=xe.util.grid_2d(mnlon, mxlon, d_lon, mnlat, mxlat, d_lat)
    if os.path.exists(Herbieregrid_coef_file):
        print('Re-using weights')
        regridder = xe.Regridder(ds, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
        ds_regrid = regridder(ds)
    else:
        print('Generating weights')
        regridder = xe.Regridder(ds, varint, "bilinear",periodic=False,unmapped_to_nan=True)
        regridder.to_netcdf(Herbieregrid_coef_file)
        ds_regrid = regridder(ds)
        
    ds_regrid=ds_regrid.drop('time')
    ds_regrid=ds_regrid.rename({'valid_time':'time'})

    
    ds_regrid=ds_regrid.to_array(name=oVariables[var][0]).squeeze()

    et = time.time()    
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    

    

    return ds_regrid 
    
def get_3d_variable(var,gridflag):        
    #Get non-radiation variables .
    st = time.time()
    print(f'Accessing data for {var}')

    fcsturl=baseurl+'/sfc/'+today+'/'+today+'_00z_fcst.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]+'/'+oVariables[var][1]
    fcstgurl=baseurl+'/sfc/'+today+'/'+today+'_00z_fcst.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]
    anlurl=baseurl+'/sfc/'+today+'/'+today+'_00z_anl.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]+'/'+oVariables[var][1]
    anlgurl=baseurl+'/sfc/'+today+'/'+today+'_00z_anl.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]

    #print(fcsturl)
    #print(anlurl)
    #Initialize S3 filesystem
    fs = s3fs.S3FileSystem(anon=True)
    ds_fcst = xr.open_mfdataset([s3fs.S3Map(url, s3=fs) for url in [fcsturl, fcstgurl]], engine='zarr')
    ds_anl = xr.open_mfdataset([s3fs.S3Map(url, s3=fs) for url in [anlurl, anlgurl]], engine='zarr')
    
    #Calculate lon/lat and merge analysis and forecast
    ds_anl = ds_anl.rename(projection_x_coordinate="x", projection_y_coordinate="y")
    ds_anl = ds_anl.metpy.assign_crs(projection.to_cf())
    ds_anl = ds_anl.metpy.assign_latitude_longitude()   
    ds_anl = ds_anl.expand_dims(dim='time',axis=0)
    ds_anl = ds_anl.set_coords(('time'))
    ds_anl = ds_anl.drop(('height','pressure'))

    ds_fcst = ds_fcst.rename(projection_x_coordinate="x", projection_y_coordinate="y")
    ds_fcst = ds_fcst.metpy.assign_crs(projection.to_cf())
    ds_fcst = ds_fcst.metpy.assign_latitude_longitude()   

    ds=xr.concat([ds_anl,ds_fcst],'time')    
    print(ds.dims)
    ds=ds.isel(x=slice(0,-icut))
    print(ds.dims)
    ds=ds.chunk(chunks={'x':ds.dims['x'],'y':ds.dims['y'],'time':1})
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})
    ds_TMP=ds[oVariables[var][0]].drop(['x','y'])
    
    varint=xe.util.grid_2d(mnlon, mxlon, d_lon, mnlat, mxlat, d_lat)
    ds_TMP = ds_TMP.where( ~np.isnan(ds_TMP), other=-999)
 
    if gridflag:
        if os.path.exists(regrid_coef_file):
            print('Re-using weights')
            regridder = xe.Regridder(ds_TMP, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
            ds_regrid = regridder(ds_TMP)
        else:
            print('Generating weights')
            regridder = xe.Regridder(ds_TMP, varint, "bilinear",periodic=False,unmapped_to_nan=True)
            regridder.to_netcdf(regrid_coef_file)
            ds_regrid = regridder(ds_TMP)
        

        ds_regrid=ds_regrid.drop(['metpy_crs'])
    else:
        ds_regrid=ds_TMP.drop(['metpy_crs'])       

    et = time.time()    
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')

    

    return ds_regrid 
def get_PRATE_INIT(var,ids):       
    #Get rain rate.
    st = time.time()
    print(f'Accessing data for {var}, first step')
    fcsturl=baseurl+'/sfc/'+yesterday+'/'+yesterday+'_23z_fcst.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]+'/'+oVariables[var][1]
    fcstgurl=baseurl+'/sfc/'+yesterday+'/'+yesterday+'_23z_fcst.zarr/'+oVariables[var][1]+'/'+oVariables[var][0]
        
    #Initialize S3 filesystem
    fs = s3fs.S3FileSystem(anon=True)
    ds_fcst = xr.open_mfdataset([s3fs.S3Map(url, s3=fs) for url in [fcsturl, fcstgurl]], engine='zarr')
    
    #Calculate lon/lat and merge analysis and forecast


    ds_fcst = ds_fcst.rename(projection_x_coordinate="x", projection_y_coordinate="y")
    ds_fcst = ds_fcst.metpy.assign_crs(projection.to_cf())
    ds_fcst = ds_fcst.metpy.assign_latitude_longitude()   
    ds=ds_fcst
    print(ds.dims)
    ds=ds.isel(x=slice(0,-icut))
    print(ds.dims)
    
    ds=ds.chunk(chunks={'x':ds.dims['x'],'y':ds.dims['y'],'time':1})
    ds = ds.rename({"latitude": "lat", "longitude": "lon"})
    ds_TMP=ds.drop(['x','y','forecast_period','forecast_reference_time'])
    
    varint=xe.util.grid_2d(mnlon, mxlon, d_lon, mnlat, mxlat, d_lat)

    ds_TMP = ds_TMP.where( ~np.isnan(ds_TMP), other=Fillvalue)
    
    if os.path.exists(regrid_coef_file):
        print('Re-using weights')
        regridder = xe.Regridder(ds, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
        ds_regrid = regridder(ds_TMP)
    else:
        print('Generating weights')
        regridder = xe.Regridder(ds, varint, "bilinear",periodic=False,unmapped_to_nan=True)
        regridder.to_netcdf(regrid_coef_file)
        ds_regrid = regridder(ds)
        

    ds_regrid=ds_regrid.drop(['metpy_crs'])

    et = time.time()    
    elapsed_time = et - st

    ids.load()
    ids.values[0,:,:]=ds_regrid.PRATE.values[0,:,:]
    
    print('Execution time:', elapsed_time, 'seconds')

    return ids 


if __name__ == "__main__":
    print('Running')
    
    main()
  