"""
acquire_insitu_CFRF_obs_v1.py

Download Insitu obs data from the Commercial Fisheries Research Foundation 
ERDDAP server. 

Author Elias Hunter hunter@marine.rutgers.edu created 2025/09/30
"""

import os,sys,fnmatch,glob
from datetime import datetime,timedelta,date
import os,glob
import copernicusmarine as copernicus_marine
import pandas as pd
import numpy as np
import netCDF4 as nc
import time
#indir='/home/om/cron/ECCOFS_OBS/INSITU_NRT/data/raw/'
#outdir='/home/om/cron/ECCOFS_OBS/INSITU_NRT/data/processed/'
url='erddap.ondeckdata.com'

ERD_ID1='emolt_data'
ERD_ID2='fixed_gear_oceanography'
ERD_ID3='mobile_gear_oceanography'
ERD_ID4='shelf_fleet_profiles_1m_binned'
ERD_ID5='wind_farm_profiles_1m_binned'

dt=timedelta(days=1)
latmin=-3.0
latmax=53.0
lonmin=-100.0
lonmax=-38.0



def main(fconfig):
    
   # indir=fconfig['obs']['insitu']['indir']
   # outdir=fconfig['obs']['insitu']['outdir']
   # nday=fconfig['obs']['insitu']['nday']
   
   
    indir='/home/om/cron/ECCOFS_OBS/CFRF/data/raw'
    outdir='/home/om/cron/ECCOFS_OBS/CFRF/data/processed'
    nday=7
    
    
    
    print('Getting CFRF ERDDAP data')   
    start_date = datetime.today()
    start_date = datetime(start_date.year, start_date.month, start_date.day, tzinfo=start_date.tzinfo)
    date_times = [start_date  - timedelta(days=i) for i in range(nday)] 
    
    print('Getting CFRF Data')
    
    for et in (date_times):
        print('------')
        st=et-dt
        print(st)
        print(et)
        ofile=f'ECCOFS_INSITU_CFRF_{st.strftime('%Y%m%d')}'
        ncfile=os.path.join(outdir,ofile+'.nc')
        get_cfrf(et,st,ncfile,fconfig) 
      #  infile=os.path.join(indir,ofile+'.csv')
     #  ncfile=os.path.join(outdir,ofile+'.nc')
    #    process(infile,ncfile,fconfig)
        
def get_cfrf(etime,stime,ofile,fconfig):
    from erddapy import ERDDAP
    n=0
    ds = ERDDAP(server=f'https://{url}/erddap',protocol="tabledap")
    ds.response = "csv"
    ds.dataset_id =ERD_ID1
    ds.constraints = {
        "time>=": stime,
        "time<=": etime,
        "latitude>=": latmin,
        "latitude<=": latmax,
        "longitude>=": lonmin,
        "longitude<=": lonmax,
        }
    ds.variables = [
        "latitude",
        "longitude",
         "time",
         "depth",
         "temperature",
         ]

    
    df=create_obs_dataframe()
    try:
        ds1 = ds.to_pandas()
      #  print(ds1.columns.to_list())
        # ['latitude (degrees_north)', 'longitude (degrees_east)', 'time (UTC)', 
        #'depth (m)', 'temperature (degree_C)']
        
        df['obs_lat']=ds1['latitude (degrees_north)']
        df['obs_lon']=ds1['longitude (degrees_east)']
        df['obs_time']=ds1['time (UTC)']
        df['obs_depth']=ds1['depth (m)']
        df['obs_value']=ds1['temperature (degree_C)']
        df['obs_provenance']=732.0
        df['obs_type']=6.0
    
       
    except Exception as e:
         print(f'{ds.dataset_id} Failed')
         print(e)

    ds.variables = [
        "latitude",
        "longitude",
        "time",
        "depth",
        "temperature",
             ]
    ds.dataset_id =ERD_ID2
    try:
        ds2 = ds.to_pandas()
      #  ['latitude (degrees_north)', 'longitude (degrees_east)', 'time (UTC)',
      #'depth (m)', 'temperature (degree_C)']
      #  print(ds2.columns.to_list())
        tdf=create_obs_dataframe()
        tdf['obs_lat']=ds2['latitude (degrees_north)']
        tdf['obs_lon']=ds2['longitude (degrees_east)']
        tdf['obs_time']=ds2['time (UTC)']
        tdf['obs_depth']=ds2['depth (m)']
        tdf['obs_value']=ds2['temperature (degree_C)']
        tdf['obs_provenance']=732.0
        tdf['obs_type']=6.0
        df = pd.concat([df, tdf], ignore_index=True)

        
        #n=n+len(ds2.time)
    except Exception as e:
        print(f'{ds.dataset_id} Failed')
       # print(e)
    
    ds.variables = [
        "latitude",
        "longitude",
        "time",
        "pressure",
        "temperature",
        ]        
    ds.dataset_id =ERD_ID3
    try:
         ds3 = ds.to_pandas()
 #       ['latitude (degrees_north)', 'longitude (degrees_east)', 'time (UTC)', 
 #'pressure (dbar)', 'temperature (degree_C)']
         tdf=create_obs_dataframe()
         tdf['obs_lat']=ds3['latitude (degrees_north)']
         tdf['obs_lon']=ds3['longitude (degrees_east)']
         tdf['obs_time']=ds3['time (UTC)']
         tdf['obs_depth']=ds3['pressure (dbar)']
         tdf['obs_value']=ds3['temperature (degree_C)']
         tdf['obs_provenance']=732.0
         tdf['obs_type']=6.0
         df = pd.concat([df, tdf], ignore_index=True)
    except Exception as e:
        print(f'{ds.dataset_id} Failed')
        #print(e)
        
    ds.dataset_id =ERD_ID4
    
    ds.variables = [
        "latitude",
        "longitude",
        "practical_salinity",
        "time",
        "sea_pressure",
        "temperature"]
    try:
         ds4 = ds.to_pandas()
  #      ['latitude (degrees_north)', 'longitude (degrees_east)', 
  #`'practical_salinity (PSU)', 'time (UTC)', 'sea_pressure (dbar)',
  #'temperature (degree_C)']
         tdf=create_obs_dataframe()
         tdf['obs_lat']=ds4['latitude (degrees_north)']
         tdf['obs_lon']=ds4['longitude (degrees_east)']
         tdf['obs_time']=ds4['time (UTC)']
         tdf['obs_depth']=ds4['sea_pressure (dbar)']
         tdf['obs_value']=ds4['temperature (degree_C)']
         tdf['obs_provenance']=732.0
         tdf['obs_type']=6.0
         df = pd.concat([df, tdf], ignore_index=True)
         
         
         tdf=create_obs_dataframe()
         tdf['obs_lat']=ds4['latitude (degrees_north)']
         tdf['obs_lon']=ds4['longitude (degrees_east)']
         tdf['obs_time']=ds4['time (UTC)']
         tdf['obs_depth']=ds4['sea_pressure (dbar)']
         tdf['obs_value']=ds4['practical_salinity (PSU)']
         tdf['obs_provenance']=732.0
         tdf['obs_type']=7.0
         df = pd.concat([df, tdf], ignore_index=True)
         
    except Exception as e:
        print(f'{ds.dataset_id} Failed')
    if not(df.empty):
    
       ds.dataset_id =ERD_ID5
       
       ds.variables = [
           "latitude",
           "longitude",
           "practical_salinity",
           "time",
           "sea_pressure",
           "temperature"]
       try:
            ds5 = ds.to_pandas()
     #      ['latitude (degrees_north)', 'longitude (degrees_east)', 
     #`'practical_salinity (PSU)', 'time (UTC)', 'sea_pressure (dbar)',
     #'temperature (degree_C)']
            tdf=create_obs_dataframe()
            tdf['obs_lat']=ds5['latitude (degrees_north)']
            tdf['obs_lon']=ds5['longitude (degrees_east)']
            tdf['obs_time']=ds5['time (UTC)']
            tdf['obs_depth']=ds5['sea_pressure (dbar)']
            tdf['obs_value']=ds5['temperature (degree_C)']
            tdf['obs_provenance']=732.0
            tdf['obs_type']=6.0
            df = pd.concat([df, tdf], ignore_index=True)
            
            
            tdf=create_obs_dataframe()
            tdf['obs_lat']=ds5['latitude (degrees_north)']
            tdf['obs_lon']=ds5['longitude (degrees_east)']
            tdf['obs_time']=ds5['time (UTC)']
            tdf['obs_depth']=ds5['sea_pressure (dbar)']
            tdf['obs_value']=ds5['practical_salinity (PSU)']
            tdf['obs_provenance']=732.0
            tdf['obs_type']=7.0
            df = pd.concat([df, tdf], ignore_index=True)
            
       except Exception as e:
           print(f'{ds.dataset_id} Failed')
           
       if not(df.empty):
        write_obs_file(df,ofile,fconfig)
    else:
        print(f'NO DATA: {ofile} not written')
# def process(infile,ncfile,fconfig):
#     print('Processing')
#     df=pd.read_csv(infile,dtype={'platform_id':'str'})
    
#     #QAQC
#     df= df[df['value_qc'] == 1]
#     df = df[~((df['variable'] == 'PSAL') & (df['value'] > fconfig['obs']['insitu']['saltcut']))]

    
#     df['PROV'] = df['platform_type'].map(fconfig['obs']['insitu']['prov_list'])
#     df['obs_type'] = df['variable'].map(fconfig['obs']['insitu']['var_list'])

#     df = df.sort_values(by='time')
#     df=df.drop(['variable','platform_type','platform_id','pressure','value_qc','institution','doi','is_depth_from_producer'],axis=1)
#     write_obs_file(df,ncfile,fconfig)
        
                                    
def create_obs_dataframe():
    df = pd.DataFrame(columns=['obs_lat',
            'obs_lon',
            'obs_time',
            'obs_depth',
            'obs_value',
            'obs_provenance'
            'obs_type'
            ])
    return df   
def  write_obs_file(dfobs,ofile,fconfig):
    print('WRITING file: '+ofile )
    
    rtime=np.datetime64(fconfig['obs']['insitu']['rtime'])
    dtime=pd.to_datetime(rtime)
    tunits=dtime.strftime('days since %Y-%m-%d 00:00:00')


    dfobs=dfobs.set_index('obs_time')
  #  print(dfobs)
    t=pd.to_datetime(dfobs.index.values)
 #   print(t)
    t=t.tz_convert(None)
    timeout=(t-rtime) / np.timedelta64(1, 'D') 
  #  print(dfobs)
    dobsx=dfobs.to_xarray()  

    dobsx['obs_lat']=dobsx['obs_lat'].astype(np.float32)
    dobsx['obs_lon']=dobsx['obs_lon'].astype(np.float32)
    dobsx.attrs['Type']='ECOFFS CFRF Observations file'
    dobsx.attrs['Title']='Commercial Fisheries Research Foundation observations in the ECCOFS domain.'
    dobsx.attrs['Source']='https://erddap.ondeckdata.com/erddap/tabledap'
    dobsx.attrs['history']='Created by Eli Hunter (hunter@marine.rutgers.edu),' +datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    dobsx.attrs['Conventions']='CF-1.8'
    

    dobsx['obs_time']=timeout
    dobsx['obs_time'].attrs['units']=tunits
    dobsx['obs_time'].attrs['long_name']='time'
    dobsx['obs_time'].attrs['calender']='gregorian'
    
    dobsx['obs_lat'].attrs['units']='degrees_north'
    dobsx['obs_lat'].attrs['long_name']='latitude'
    
    dobsx['obs_lon'].attrs['units']='degrees_east'
    dobsx['obs_lon'].attrs['long_name']='longitude'
    
    dobsx['obs_depth'].attrs['units']='m'
    dobsx['obs_depth'].attrs['long_name']='depth'
    dobsx['obs_depth'].attrs['orientation']='positive up'
    

    dobsx['obs_provenance'].attrs['long_name']='Provenance'
    dobsx['obs_provenance'].attrs['units']='None'
    dobsx['obs_provenance'].attrs['source']='https://tds.marine.rutgers.edu/erddap/tabledap/ECCOFS_PROVENANCE.html'
    
    dobsx['obs_type'].attrs['long_name']='Provenance'
    dobsx['obs_type'].attrs['units']='unity'
    dobsx['obs_type'].attrs['option_01']="free-surface" 
    dobsx['obs_type'].attrs['option_02']="vertically integrated u-momentum component" 
    dobsx['obs_type'].attrs['option_03']="vertically integrated v-momentum component" 
    dobsx['obs_type'].attrs['option_04']="u-momentum component" 
    dobsx['obs_type'].attrs['option_05']="v-momentum component" 
    dobsx['obs_type'].attrs['option_06']="potential temperature" 
    dobsx['obs_type'].attrs['option_07']="salinity" 
    
    dobsx['obs_value'].attrs['units']='state variable units'
    dobsx['obs_value'].attrs['long_name']='observation value' 
        

    # dobsx['I'].attrs['units']=' '
    # dobsx['I'].attrs['long_name']='ECCOFS ROMS I coordinate' 
    # dobsx['J'].attrs['units']=' '
    # dobsx['J'].attrs['long_name']='ECCOFS ROMS J coordinate' 
    
    dobsx.to_netcdf(path=ofile,unlimited_dims='obs_time')





###########################################################
#After function definitions
###########################################################
###########################################################
# gfile='/home/om/roms/eccofs/grid_eccofs_6km_02.nc'
# gds=xr.open_dataset(gfile)    
# interps=create_IJinterpolator(gds)
# Iinterp=interps[0]
# Jinterp=interps[1]

    
    
#main(sys.argv)
