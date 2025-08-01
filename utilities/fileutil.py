import os,glob
import subprocess,shutil
import datetime,time
from datetime import timedelta,date
import xarray as xr

def concat_met(fconfig):
    print('Concatenating Met files')
    local_dir=fconfig['transfer']['ldesdir']
    srcdir=fconfig['transfer']['trlist']['met']['srcdir']
    prefix=fconfig['transfer']['trlist']['met']['prefix']
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
            file=glob.glob(fconfig['transfer']['trlist']['met']['srcdir']+f'*{var}*{today}*.nc')
            ds=xr.open_dataset(file[0])
            subset = ds.isel(time=slice(0, 24))
            dslist.append(subset)
        day=start_date+timedelta(days=nday2-1)
        today=day.strftime('%Y%m%d')
        file=glob.glob(fconfig['transfer']['trlist']['met']['srcdir']+f'*{var}*{today}*.nc')
        ds=xr.open_dataset(file[0])
        dslist.append(ds)
        combined = xr.concat(dslist, dim='time')    
        combined.to_netcdf(ofile)
    #combined = xr.open_mfdataset(files,concat_dim='time',combine='nested')
    #combined.to_netcdf(ofile)

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

        combined = xr.open_mfdataset(files,concat_dim='time',combine='nested')
        combined.to_netcdf(ofile)
        
def delete_remote(fconfig):
    remote_dir=fconfig['transfer']['rdesdir']
    remote_server=fconfig['transfer']['server']
    try:
        subprocess.check_call(['ssh',remote_server, 'rm -rf '+remote_dir+'*'],env=os.environ)
        print(f"Successfully deleted files in {remote_dir}")

    except Exception as e:
        print(f"Error: {str(e)}")
        print(f"File delete: {remote_dir} failed")
    


def delete_local(fconfig):
    local_dir=fconfig['transfer']['ldesdir']
    try:
        for filename in os.listdir(local_dir):
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
#     nfiles=fconfig['transfer']['nfiles']
    for key in keys:
        
        sp=fconfig['transfer']['trlist'][key]['searchpattern']
        match key:
            
            case 'met':
                concat_met(fconfig)
            case 'clm':
                concat_clm(fconfig)
                
            case _:
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
            
            
        