"""
get_PODAAC_L3S_LEO_PM_N.py

Subset, donwload, and process HRRR data from AWS PODAAC bucket s3://podaac-ops-cumulus-protected/. 

Daily combined SST 

Created by Elias Hunter, hunter@marine.rutgers.edu, 1/23/2024 

"""


import zarr 
import xarray
import s3fs
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import requests
import boto3
import xarray as xr
import dask
import netCDF4 as nc
import time,os
from dask.distributed import Client

outdirec='/home/ubuntu/PODAAC/L3S_LEO_AM_N/'

nday=5
today = datetime.now().date()

date_times = [today - timedelta(days=i) for i in range(nday)]  # This creates a list for the next 5 days

s3_url='s3://podaac-ops-cumulus-protected/L3S_LEO_AM-STAR-v2.80/'
s3 = s3fs.S3FileSystem(anon=False)
encoding = {'sst_dtime': {'units': 'milliseconds'}}
filepart='-ECCOFS-L3S_GHRSST-SST-LEO_AM_N-ACSPO_V2.80-v02.0-fv01.0.nc'
fathon_bucket_name= "fathom-eccofs"
folder = "data/sst/L3S_LEO_AM_N/"
aws_access_key_id="ACCESSKEYID"
aws_secret_access_key="SECRETACCESSKEYID"
def main():
    
    print('Getting Credentials')
    get_credentials()
    
  
    s3 = s3fs.S3FileSystem(anon=False,
                        key=os.environ["AWS_ACCESS_KEY_ID"],
                        secret=os.environ["AWS_SECRET_ACCESS_KEY"],
                        token=os.environ["AWS_SESSION_TOKEN"])
    counter=1
    for dt in date_times:
        print('-----------------------')
        print(counter)
        if counter > 50:
            print('RESET CREDENTIALS')
            try:
                get_credentials()
                s3 = s3fs.S3FileSystem(anon=False,
                key=os.environ["AWS_ACCESS_KEY_ID"],
                secret=os.environ["AWS_SECRET_ACCESS_KEY"],
                token=os.environ["AWS_SESSION_TOKEN"])
                
               
                counter=1
            except:
                print('Failed to get credentials')
                
        sdate=dt.strftime('%Y%m%d')
        search_url=s3_url+sdate+'*LEO_AM_D*.nc'
        remote_files = s3.glob(search_url)
        if len(remote_files)==0:
            print('Not found: '+search_url)
            continue
        else:
            print(remote_files)        
        fileset = [s3.open(file) for file in remote_files]
        
        data = xr.open_mfdataset(fileset, concat_dim="time",combine='nested',engine='h5netcdf' )
      #  print(data)
    #    outdata=data.drop_vars('crs').sel(lat=slice(46,36),lon=slice(-78,-59))
       # outdata=data.sel(lat=slice(46,36),lon=slice(-78,-59))
        outdata=data.sel(lat=slice(52,5),lon=slice(-100,-36))
       
       # print(outdata)
        outdata['quality_level'].attrs['comment'] = 'SST quality levels: 5 corresponds to clearsky pixels'
        ofile=outdirec+sdate+filepart
        print(ofile)
        outdata.to_netcdf(ofile,encoding=encoding)
        
        ofile2=sdate+filepart
        s3boto = boto3.client("s3",
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name="us-east-2")
        s3_key = folder + ofile2
        print(f"Uploading {ofile} → s3://{fathon_bucket_name}/{s3_key}")
        s3boto.upload_file(ofile, fathon_bucket_name, s3_key)
        
def get_credentials():
    temp_creds_url='https://archive.podaac.earthdata.nasa.gov/s3credentials'
    temp_creds_req=requests.get(temp_creds_url).json()
    os.environ["AWS_ACCESS_KEY_ID"] = temp_creds_req["accessKeyId"]
    os.environ["AWS_SECRET_ACCESS_KEY"] = temp_creds_req["secretAccessKey"]
    os.environ["AWS_SESSION_TOKEN"] = temp_creds_req["sessionToken"]
    
if __name__ == "__main__":
    print('Running')
  #  client = Client()
    main()
   