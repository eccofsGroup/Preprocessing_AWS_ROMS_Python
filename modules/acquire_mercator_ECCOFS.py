#!/usr/bin/env python
#
import os,sys,fnmatch,glob
from shutil import copyfile
from datetime import datetime,timedelta,date
import subprocess,shutil,os,glob
import copernicusmarine as copernicus_marine
import pandas as pd

# outdir='/home/om/cron/ECCOFS_OBS/MERCATOR/data/raw/'
# tmpdir='/home/om/cron/ECCOFS_OBS/MERCATOR/tmp/'
# days=14
# fdays=7


# latmin=0.0
# latmax=52.0
# lonmin=-100.0
# lonmax=-36.0
# zmin=-10.0;
# zmax=6000.0
# variables=[['so'] ,['thetao'],['uo','vo'],['zos']]
# tmpfiles=['TMP_so.nc','TMP_thetao.nc','TMP_uv.nc','TMP_zos.nc']



# ids=['cmems_mod_glo_phy-so_anfc_0.083deg_P1D-m', 
#      'cmems_mod_glo_phy-thetao_anfc_0.083deg_P1D-m', 
#      'cmems_mod_glo_phy-cur_anfc_0.083deg_P1D-m', 
#      'cmems_mod_glo_phy_anfc_0.083deg_P1D-m']





def main(fconfig):
    print('Getting Mercator Data')    
    
    
    
    outdir=fconfig['force']['mercator']['outdir']
    tmpdir=fconfig['force']['mercator']['tmpdir']
    days=fconfig['force']['mercator']['days']
    fdays=fconfig['force']['mercator']['fdays']


    latmin=fconfig['force']['mercator']['latmin']
    latmax=fconfig['force']['mercator']['latmax']
    lonmin=fconfig['force']['mercator']['lonmin']
    lonmax=fconfig['force']['mercator']['lonmax']
    zmin=fconfig['force']['mercator']['zmin']
    zmax=fconfig['force']['mercator']['zmax']
    variables=fconfig['force']['mercator']['variables']
    tmpfiles=fconfig['force']['mercator']['tmpfiles']


    ids=fconfig['force']['mercator']['ids']

    
    tstart=date.today()-timedelta(days=days)
    tend=date.today()+timedelta(days=fdays)

    days=pd.date_range(tstart,tend)
    for day in days:
        print(day)


        for f in glob.glob(tmpdir+'*.nc'):
            os.remove(f)
        for f in glob.glob(tmpdir+'*.tmp'):
            os.remove(f)
            
        for ind,variable in enumerate(variables):
            print(tmpfiles[ind])
            print(ids[ind])
            day2=day+timedelta(days=0.5)
            
            print(variable)
            copernicus_marine.subset(
                dataset_id = ids[ind],
                variables = variable,
                minimum_longitude = lonmin,
                maximum_longitude = lonmax,
                minimum_latitude = latmin,
                maximum_latitude = latmax,
                minimum_depth = zmin,
                maximum_depth = zmax,
                start_datetime = day.strftime('%Y-%m-%d %H:%M:%S'),
                end_datetime = day2.strftime('%Y-%m-%d %H:%M:%S'),
                output_filename = tmpfiles[ind],
           #     force_download=True,
            #    overwrite_output_data =True,
            #    dataset_version="202406",
                output_directory = tmpdir,
                )
        
        ofilename = day.strftime(outdir+'mercator_doppio_%Y_%m_%d.nc')
        tmpofilename = day.strftime(tmpdir+'mercator_doppio_%Y_%m_%d.nc')
        print(ofilename)
        
        try:
            print('Copying: '+tmpdir+tmpfiles[0]+' to '+tmpofilename)
            shutil.copyfile(tmpdir+tmpfiles[0], tmpofilename)
    
            
            for tmpfile in tmpfiles[1:]:
                print('Adding: '+tmpdir+tmpfile+' to '+tmpofilename)
                subprocess.check_call(['ncks','-A',tmpdir+tmpfile,tmpofilename])
        except Exception as e:
                # Code to handle other exceptions (if any)
            print(f"An error occurred when combining files via ncks: {e}")
            
            if os.path.exists(tmpofilename):
                print(f'Deleting: {tmpofilename}') 
        else:
            shutil.copyfile(tmpofilename,ofilename)
            

    for f in glob.glob(tmpdir+'*.nc'):
        os.remove(f)
    for f in glob.glob(tmpdir+'*.tmp'):
        os.remove(f)
# if __name__ == "__main__":

    
#     print('RUNNING')
#     main(sys.argv)
#     for f in glob.glob(tmpdir+'*.nc'):
#         os.remove(f)
#     for f in glob.glob(tmpdir+'*.tmp'):
#         os.remove(f)
