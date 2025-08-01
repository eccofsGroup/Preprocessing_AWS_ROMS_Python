"""
acquire_insitu_obs_v1.py

Download Insitu obs data from the copernicus marine store for the last 3 days 

Author Elias Hunter hunter@marine.rutgers.edu created 2025/05/21
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


# latmin=-3.0
# latmax=53.0
# lonmin=-100.0
# lonmax=-38.0
# zmin=-10.0
# zmax=6000.0



dt=timedelta(days=1)


# dataset_id="cmems_obs-ins_glo_phybgcwav_mynrt_na_irr"
# dataset_part="latest"
# dataset_version="202311"
# variables=["PSAL", "TEMP"]



#prov_list={'PF':801, 'XX':802, 'ML':803, 'TX':804, 'XB':805, 'GL':806, 'CT':807, 'SM':808, 'BO':809,'FB':821, 'SD':822, 'MO':823, 'TS':824, 'DB':825, 'TG':826}
#var_list={'PSAL':7.0,'TEMP':6.0}

def main(fconfig):
    
    indir=fconfig['obs']['insitu']['indir']
    outdir=fconfig['obs']['insitu']['outdir']
    nday=fconfig['obs']['insitu']['nday']
    
    start_date = date.today()
    date_times = [start_date  - timedelta(days=i) for i in range(nday)] 
    

    
    for t in (date_times):
        print(t)
        st=t-dt
        ofile=f'ECCOFS_NRT_INSITU_CMEMS_{st.strftime('%Y%m%d')}'
        get_nrt(t,st,ofile,fconfig) 
        infile=os.path.join(indir,ofile+'.csv')
        ncfile=os.path.join(outdir,ofile+'.nc')
        process(infile,ncfile,fconfig)
        
def get_nrt(etime,stime,ofile,fconfig):

    
    
    copernicus_marine.subset(
        dataset_id=fconfig['obs']['insitu']['dataset_id'],
        dataset_part=fconfig['obs']['insitu']['dataset_part'],
        dataset_version=fconfig['obs']['insitu']['dataset_version'],
        variables=fconfig['obs']['insitu']['variables'],
        minimum_longitude=fconfig['obs']['insitu']['lonmin'],
        maximum_longitude=fconfig['obs']['insitu']['lonmax'],
        minimum_latitude=fconfig['obs']['insitu']['latmin'],
        maximum_latitude=fconfig['obs']['insitu']['latmax'],
        start_datetime=stime.strftime('%Y-%m-%d'),
        end_datetime=etime.strftime('%Y-%m-%d'),
        minimum_depth=fconfig['obs']['insitu']['zmin'],
        maximum_depth=fconfig['obs']['insitu']['zmax'],
        coordinates_selection_method="strict-inside",
        output_filename = ofile,
        output_directory =  fconfig['obs']['insitu']['indir'],
        overwrite=True,
        disable_progress_bar=True,
        )
    
def process(infile,ncfile,fconfig):
    print('Processing')
    df=pd.read_csv(infile,dtype={'platform_id':'str'})
    
    #QAQC
    df= df[df['value_qc'] == 1]
    df = df[~((df['variable'] == 'PSAL') & (df['value'] > fconfig['obs']['insitu']['saltcut']))]

    
    df['PROV'] = df['platform_type'].map(fconfig['obs']['insitu']['prov_list'])
    df['obs_type'] = df['variable'].map(fconfig['obs']['insitu']['var_list'])

    df = df.sort_values(by='time')
    df=df.drop(['variable','platform_type','platform_id','pressure','value_qc','institution','doi','is_depth_from_producer'],axis=1)
    write_obs_file(df,ncfile,fconfig)
        
                                    
            
                       
   
def  write_obs_file(df,ofile,fconfig):
    print('WRITING file: '+ofile )
    
    rtime=np.datetime64(fconfig['obs']['insitu']['rtime'])
    dtime=pd.to_datetime(rtime)
    tunits=dtime.strftime('days since %Y-%m-%d 00:00:00')

    dfobs=df.rename(columns={'latitude':'obs_lat',
                               'longitude':'obs_lon',
                               'time':'obs_time',
                               'depth':'obs_depth',
                               'value':'obs_value',
                               'PROV':'obs_provenance'})
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
    dobsx.attrs['Type']='ECOFFS in-situ Observations file'
    dobsx.attrs['Title']='Copernicus in-situ observations in the ECCOFS domain.'
    dobsx.attrs['Source']='https://data.marine.copernicus.eu/'
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
    dobsx['obs_provenance'].attrs['units']=' '
    
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
