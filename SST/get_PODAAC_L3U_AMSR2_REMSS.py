"""
get_PODAAC_L3U_AMSR2_REMSS.py

Subset, download process AMSR2 data from AWS PODAAC bucket s3://podaac-ops-cumulus-protected/. 

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
import subprocess
from dask.distributed import Client


#THINGS TO EDIT
outdirec='/home/ubuntu/PODAAC/AMSR2_L3/'
nday=7
today = datetime.now().date()
date_times = [today - timedelta(days=i) for i in range(nday)]  # This creates a list for the next 5 days

s3_url='s3://podaac-ops-cumulus-protected/AMSR2-REMSS-L3U-v8.2/'
#s3_url='s3://ppodaac-ops-cumulus-protected/AMSR2-REMSS-L3U_RT-v8.2/'
encoding = {'sst_dtime': {'units': 'milliseconds'}}
filepart='-ECCOFS-REMSS-L3U_GHRSST-SSTsubskin-AMSR2-RSS_daily.nc'
private_key_path = 'PATH_TO_YOUR_PRIVATE KEY'
remote_file_path = 'USER@REMOTESERVER:PATHTOSERVER'

fathon_bucket_name= "fathom-eccofs"
folder = "data/sst/amsr2/"
aws_access_key_id="ACCESSKEYID"
aws_secret_access_key="SECRETACCESSKEYID"

def main():
    
    for dt in date_times:
        print('-----------------------')
        print('Getting Credentials')
        get_credentials()
        sdate=dt.strftime('%Y%m%d')
        #search_url=s3_url+sdate+'*REMSS-L3U_GHRSST-SSTsubskin-AMSR2-RSS_AMSR2_ocean_L3_daily*.nc'
        search_url=s3_url+sdate+'*REMSS*_daily*.nc'
        s3 = s3fs.S3FileSystem(anon=False)
        remote_files = s3.glob(search_url)
        if len(remote_files)==0:
            print('Not found: '+search_url)
            continue
        else:
            print(remote_files)
  
        #fileset = [s3.open(file) for file in remote_files]
        fileset = s3.open(remote_files[0]) 
        
        data = xr.open_dataset(fileset, engine='h5netcdf', decode_times=False )
        
      #  print(data)
    #    outdata=data.drop_vars('crs').sel(lat=slice(46,36),lon=slice(-78,-59))
        outdata=data.sel(lat=slice(0,53),lon=slice(-100,-36))

       # print(outdata)
        outdata['quality_level'].attrs['comment'] = 'SST quality levels: 5 corresponds to clearsky pixels'
        ofile=outdirec+sdate+filepart
        print(ofile)
       # outdata.to_netcdf(ofile,encoding=encoding)
        outdata.to_netcdf(ofile)
        transfer_file(ofile)
        
        ofile2=sdate+filepart
        s3boto = boto3.client("s3",
            aws_access_key_id=aws_access_key_id,
            aws_secret_access_key=aws_secret_access_key,
            region_name="us-east-2")
        s3_key = folder + ofile2
        print(f"Uploading {ofile} → s3://{fathon_bucket_name}/{s3_key}")
        s3boto.upload_file(ofile, fathon_bucket_name, s3_key)
        
        
        os.remove(ofile)
        
        
        
def transfer_file(infile):
    print('Transfering: '+infile  )
  #  ssh = paramiko.SSHClient()
 #   ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
         subprocess.check_call(['scp',infile,remote_file_path],env=os.environ)
        
    #     private_key = paramiko.RSAKey.from_private_key_file(private_key_path)

    # # Connect to the SFTP server using the private key
    #     ssh.connect(hostname, port, username, pkey=private_key)

    # # Open an SFTP session
    #     sftp = ssh.open_sftp()

    # # Upload the file to the remote server
    #     sftp.put(infile, remote_file_path)
         print(f"File {infile} successfully uploaded to {remote_file_path}")

    # # Close the SFTP session
    #     sftp.close()

    # # Close the SSH connection
        
    
    
    except Exception as e:
        print(f"Error: {str(e)}")
        
  #  ssh.close()
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
    try: 
        os.remove(outdirec+'*.*')
    except:
        print('no files to remove')
        
   