import os,glob
import subprocess,shutil
import datetime,time
from datetime import timedelta,date
import numpy as np
import xarray as xr
import boto3

def concat_met(fconfig):
    
    print('Concatenating Merged Met files')
    local_dir=fconfig['transfer']['ldesdir']
    srcdir=fconfig['transfer']['trlist']['met']['srcdir_m']
    prefix=fconfig['transfer']['trlist']['met']['prefix_m']
    varnames=fconfig['transfer']['trlist']['met']['vars']
    nday=fconfig['transfer']['nday']
    nday2=fconfig['transfer']['trlist']['met']['nday']
    start_date = date.today()-timedelta(days=nday)
    
    alltoday=start_date.strftime('%Y%m%d')
    
    for var in varnames:
         ofile=f'{local_dir}{prefix}_{var}_{alltoday}.nc'
        
         dslist=[]
         for dt in range(0,nday2-1):
             day=start_date+timedelta(days=dt)
             today=day.strftime('%Y%m%d')
             file=glob.glob(srcdir+f'*{var}_{today}*.nc')
             
             ds=xr.open_dataset(file[0])
             subset = ds.isel(time=slice(0, 24))
             dslist.append(subset)
         day=start_date+timedelta(days=nday2-1)
         today=day.strftime('%Y%m%d')
         file=glob.glob(srcdir+f'*{var}_{today}*.nc')
         if len(file)>0:
             
             ds=xr.open_dataset(file[0])
             dslist.append(ds)
         else:
             tmp=srcdir+f'*{var}_{today}*.nc'
             print(f'No Met file: {tmp}')
    
         last_time = ds.time.isel(time=-1).values
         
         combined = xr.concat(dslist, dim='time')    
         combined.to_netcdf(ofile)
        
                
         srcdir2=fconfig['transfer']['trlist']['met']['srcdir']
       #  ofile=f'{local_dir}{prefix}_{var}_CONTINUED_{alltoday}.nc'
         ofile=ofile.replace('MERGED','CONTINUE')
         file=glob.glob(srcdir2+f'*{var}_{today}*.nc')
        
         if len(file)>0:
             ds=xr.open_dataset(file[0])
             subset = ds.where(ds.time > last_time, drop=True)
             #print(var)
             #print(subset)
             subset[var].attrs['coordinates'] = "lon lat" 
             subset.to_netcdf(ofile)
         else:
             tmp=srcdir2+f'*{var}_{today}*.nc'
             print(f'No Met file: {tmp}')
    
    
    ##############################################
    #Commented code for GFS only 
    ##############################################
    # print('Concatenating GFS Met files')
    # local_dir=fconfig['transfer']['ldesdir']
    # srcdir=fconfig['transfer']['trlist']['met']['srcdir']
    # prefix=fconfig['transfer']['trlist']['met']['prefix']
    # varnames=fconfig['transfer']['trlist']['met']['vars']
    # nday=fconfig['transfer']['nday']
    # nday2=fconfig['transfer']['trlist']['met']['nday']
    # start_date = date.today()-timedelta(days=nday)
    
    # alltoday=start_date.strftime('%Y%m%d')
    
    # for var in varnames:
    #     ofile=f'{local_dir}{prefix}_{var}_{alltoday}.nc'
    #     dslist=[]
    #     for dt in range(0,nday2-1):
    #         day=start_date+timedelta(days=dt)
    #         today=day.strftime('%Y%m%d')
    #         file=glob.glob(srcdir+f'*{var}_{today}*.nc')
    #         ds=xr.open_dataset(file[0])
    #         subset = ds.isel(time=slice(0, 24))
    #         dslist.append(subset)
    #     day=start_date+timedelta(days=nday2-1)
    #     today=day.strftime('%Y%m%d')
    #     file=glob.glob(fconfig['transfer']['trlist']['met']['srcdir']+f'*{var}*{today}*.nc')
    #     ds=xr.open_dataset(file[0])
    #     dslist.append(ds)
    #     combined = xr.concat(dslist, dim='time')    
    #     combined.to_netcdf(ofile)


def concat_clm(fconfig):
        print('Concatenating Climatology files')
        local_dir=fconfig['transfer']['ldesdir']
        prefix=fconfig['transfer']['trlist']['clm']['prefix']
        nday=fconfig['transfer']['nday']
        nday2=fconfig['transfer']['trlist']['clm']['nday']
        start_date = date.today()-timedelta(days=nday)
        
        today=start_date.strftime('%Y%m%d')
        ofile=f'{local_dir}{prefix}_{today}_clm.nc'
        files=[]
        for dt in range(0,nday2):
            day=start_date+timedelta(days=dt)
            today=day.strftime('%Y%m%d')
            file=glob.glob(fconfig['transfer']['trlist']['clm']['srcdir']+f'*{today}*.nc')
            files.append(file[0])

        combined = xr.open_mfdataset(files,concat_dim='ocean_time',combine='nested')
        
        last = combined.isel(ocean_time=-1)


        
        new_time = last['ocean_time'].values+np.timedelta64(5,'D') # ,your desired time
        last['ocean_time'].values = new_time
        combined = xr.concat([combined, last],dim="ocean_time")

        combined.to_netcdf(ofile)
        
def append_clm(fconfig):
        from netCDF4 import Dataset
        print('Appending Climatology files')
        local_dir=fconfig['transfer']['ldesdir']
        prefix=fconfig['transfer']['trlist']['clm']['prefix']
        nday=fconfig['transfer']['nday']
        nday2=fconfig['transfer']['trlist']['clm']['nday']
        start_date = date.today()-timedelta(days=nday)
        
        today=start_date.strftime('%Y%m%d')
        ofile=f'{local_dir}{prefix}_{today}_clm.nc'
        files=[]
        for dt in range(0,nday2):
            day=start_date+timedelta(days=dt)
            today=day.strftime('%Y%m%d')
            file=glob.glob(fconfig['transfer']['trlist']['clm']['srcdir']+f'*{today}*.nc')
            files.append(file[0])

        
        with xr.open_dataset(files[0]) as ds:
            ds.to_netcdf(ofile, mode="w", unlimited_dims=["ocean_time"])
        
        with Dataset(ofile, "a") as dst:
            for f in files[1:]:
                print(f)
                with Dataset(f, "r") as src:
            # Figure out current size
                    n0 = dst.dimensions["ocean_time"].size
                    n1 = src.dimensions["ocean_time"].size

            # Append time and any time-dependent variables
                    
                    dst["ocean_time"][n0:n0+n1] = src["ocean_time"][:]
                    
                    for vname, var in src.variables.items():
                        if "ocean_time" in var.dimensions and vname != "ocean_time":
                            
                            dst[vname][n0:n0+n1, ...] = var[:]
            n0 = dst.dimensions["ocean_time"].size
            with Dataset(files[-1], "r") as src:
                dst["ocean_time"][n0] = src["ocean_time"][-1]+5.0
                for vname, var in src.variables.items():
                    if "ocean_time" in var.dimensions and vname != "ocean_time":
                        dst[vname][n0, ...] = var[:]
        # for f in files[1:]:
            # with xr.open_dataset(f) as ds:
                
                # print(f)
                # ds.to_netcdf(ofile, mode="a", unlimited_dims=["ocean_time"])
        # last = xr.open_dataset(files[-1])
        # print(last)
        # new_time = last['ocean_time'].values+np.timedelta64(5,'D') # ,your desired time
        # last['ocean_time'].values = new_time
        # last.to_netcdf(ofile, mode="a", unlimited_dims=["ocean_time"])
        #ombined = xr.concat([combined, last],dim="ocean_time")
        #combined.to_netcdf(ofile)
        
        
def move_clm(fconfig):
        print('Moving  Climatology files')
        local_dir=fconfig['transfer']['ldesdir']
        prefix=fconfig['transfer']['trlist']['clm']['prefix']
        nday=fconfig['transfer']['nday']
        nday2=fconfig['transfer']['trlist']['clm']['nday']
        start_date = date.today()-timedelta(days=nday)
        
        today1=start_date.strftime('%Y%m%d')
        for dt in range(0,nday2):
            day=start_date+timedelta(days=dt)
            today=day.strftime('%Y%m%d')
            file=glob.glob(fconfig['transfer']['trlist']['clm']['srcdir']+f'*{today}*.nc')
            ind=dt+1
            ofile=f'{local_dir}{prefix}_{today1}_{ind:02d}_clm.nc'

            print(ofile)
            shutil.copy2(file[0], ofile)
            

        ds = xr.open_dataset(ofile)
    
        new_time = ds['ocean_time'].values+np.timedelta64(5,'D') # ,your desired time
        ds['ocean_time'].values = new_time
        
        ind=ind+1
        ofile=f'{local_dir}{prefix}_{today1}_{ind:02d}_clm.nc'
        ds.to_netcdf(ofile)

        
def delete_remote(fconfig):
    remote_dir=fconfig['transfer']['rdesdir']
    remote_server=fconfig['transfer']['server']
    try:
        subprocess.check_call(['ssh',remote_server, 'rm -rf '+remote_dir+'*.nc'],env=os.environ)
        print(f"Successfully deleted files in {remote_dir}")

    except Exception as e:
        print(f"Error: {str(e)}")
        print(f"File delete: {remote_dir} failed")
    


def delete_local(fconfig):
    local_dir=fconfig['transfer']['ldesdir']
    try:
        for filename in glob.glob(f'{local_dir}*.nc'):
            file_path = os.path.join(local_dir, filename)

            if os.path.isfile(file_path):
                os.remove(file_path)
        print(f"Successfully deleted files in {local_dir}")

    except Exception as e:
        print(f"Error: {str(e)}")
        print(f"File delete: {local_dir} failed")
    
        
def compile_all_files(fconfig):
    nday=fconfig['transfer']['nday']
    local_dir=fconfig['transfer']['ldesdir']    
    start_date = date.today()-timedelta(days=nday)
    today=start_date.strftime('%Y%m%d')
    keys=fconfig['transfer']['trlist'].keys()

   # move_clm(fconfig)
    for key in keys:
        
        sp=fconfig['transfer']['trlist'][key]['searchpattern']
        
        match key:
            
            case 'met':
                concat_met(fconfig)
            case 'clm':
                
                #move_clm(fconfig)
                #concat_clm(fconfig)
                append_clm(fconfig)
                
            case _:
                
                pass
                print(fconfig['transfer']['trlist'][key]['srcdir']+f'{sp}{today}*.nc')
                files=glob.glob(fconfig['transfer']['trlist'][key]['srcdir']+f'{sp}{today}*.nc')
                files=sorted(files)
                    
                for file in files:
                    try:
                        shutil.copy2(file, local_dir)
                        print(f"File {file} successfully moved to {local_dir}")
        
                    except Exception as e:
                        print(f"Error: {str(e)}")
                        print(f"File {file} move to {local_dir} failed")
                    
          
def transfer(fconfig):

    remote_file_path=fconfig['transfer']['server']+':'+fconfig['transfer']['rdesdir']
    print(remote_file_path)
    files=glob.glob(fconfig['transfer']['ldesdir']+'*')
    files=sorted(files)

    for file in files:
        try:
            subprocess.check_call(['scp',file,remote_file_path],env=os.environ)
            print(f"File {file} successfully uploaded to {remote_file_path}")

        except Exception as e:
            print(f"Error: {str(e)}")
            print(f"File {file} upload to {remote_file_path} failed")
            
  #Move to fathom_science bucket.           
    s3 = boto3.client("s3")
    bucket_name = "fathom-eccofs"
    folder = "assimilation_files/"
    
    for filename in files:
   #     filepath = os.path.join(local_folder, filename)
        if os.path.isfile(filename):  # skip subfolders
            parts=filename.split('/')
            s3_key = folder + parts[-1]
            print(f"Uploading {filename} ? s3://{bucket_name}/{s3_key}")
            s3.upload_file(filename, bucket_name, s3_key)
    