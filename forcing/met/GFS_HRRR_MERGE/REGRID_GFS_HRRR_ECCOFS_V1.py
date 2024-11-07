"""
REGRID_GFS_HRRR_ECCOFS_v1.py 

Merge GFS and HRR varaibles onto a single grid. 

Created by Elias Hunter, hunter@marine.rutgers.edu, 1/31/2023 
"""

import xarray as xr
import numpy as np
import time,os,shutil,pickle
import netCDF4 as nc
import xesmf as xe
from dask.distributed import Client
import matplotlib.pyplot as plt
import datetime,cv2
import optparse,sys

Variables = {
    'Tair': ("TMP",'2m_above_ground'), 
    'Qair': ("RH",'2m_above_ground'), 
    'Pair': ("MSLMA",'mean_sea_level'),
    'rain': ("PRATE",'surface'),
    'swrad': ("SWRAD",'surface'), 
    'lwrad': ("LWRAD",'surface'),
    'lwrad_down': ("DLWRF",'surface'), 
    'Uwind': ("UGRD",'10m_above_ground'),
    'Vwind': ("VGRD",'10m_above_ground')

    }


###############################################################
#Things to edit. 
###############################################################
direc='/home/om/cron/HRRR_AWS/data/'
odirec='/home/om/cron/HRRR_AWS/data/merged'
regrid_coef_file='/home/om/cron/HRRR_AWS/work/GFS_HRR_GRID_COEFFS.nc'
smooth_coef_file='/home/om/cron/HRRR_AWS/work/SMOOTH_ECCOFS.nc'

###############################################################


today = datetime.datetime.utcnow().strftime('%Y%m%d')


latN=2440
lonN=3680
N1=80
N2=20



def main(argv):
    
    today=opts.date
    
    for key in Variables.keys():
        print(key)
        #GETTING DATA
        Hfile=f'{direc}hrrr/HRRR_ZARR_REGRIDED_{key}_{today}.nc'

        if not(os.path.exists(Hfile)):
            print(f'NO HRRR file : {Hfile}')
            return
        else:
            print(f'{Hfile} exists')
        H_ds=xr.open_dataset(Hfile)
    #H_ds2=H_ds.copy()
        Gfile=f'{direc}gfs/GFS_GRIB_AWS_{key}_{today}.nc'    

        if not(os.path.exists(Hfile)):
            print(f'NO GFS file : {Gfile}')
            return
        else:
            print(f'{Gfile} exists')
            
            
            
        G_ds=xr.open_dataset(Gfile)
        G_ds=G_ds.sel(time=H_ds.time.values)
        Ofile=f'{odirec}/HRRR_GFS_MERGED_{key}_{today}.nc'
        if os.path.exists(Ofile):
            print("Exists, DELETING")
            os.remove(Ofile)
        shutil.copy(Hfile,Ofile)
        
   
      #REGRIDDING GFS to HRRR
        H_Mask = xr.where(np.isnan(H_ds[key].isel(time=0)), 1, 0)
        H_Mask_I= xr.where(~np.isnan(H_ds[key].isel(time=0)), 1, 0)
        print(f'Regridding:{key}_{today}')
        st = time.time()
        if os.path.exists(regrid_coef_file):
            print('Re-using gridding weights')
            regrid_mask = xe.Regridder(G_ds, H_ds, method="bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)
  #      regridder = xe.Regridder(u, varint, "bilinear",periodic=False,unmapped_to_nan=True,filename=regrid_coef_file,reuse_weights=True)

        else:
            print('Generating gridding weights')
            regrid_mask = xe.Regridder(G_ds, H_ds, method="bilinear",periodic=False,unmapped_to_nan=True)
    #    regridder = xe.Regridder(u, varint, "bilinear",periodic=False,unmapped_to_nan=True)
            regrid_mask.to_netcdf(regrid_coef_file)

        dr_out = regrid_mask(G_ds, keep_attrs=True)
        et = time.time()
  # get the execution time
        elapsed_time = et - st
        print('Gridding Execution time:', elapsed_time, 'seconds')



        st = time.time()
        if os.path.exists(smooth_coef_file):
            print('Re-using Smoothing weights')
            H_sc2=xr.open_dataset(smooth_coef_file)
            H_sc2=H_sc2.to_array()
            H_sc2=H_sc2.drop_indexes('variable')
        else:
            print('Generating Smoothing weights')
            distance_transform = cv2.distanceTransform( H_Mask_I.values.astype(np.uint8), cv2.DIST_L2, 5)
            y=(1+np.tanh((distance_transform-N1)/N2))/2
            H_sc2 = xr.DataArray(y, dims=('lat','lon'))
            H_sc2.to_netcdf(smooth_coef_file)
        

        G_sc2=1.0-H_sc2
        G_sc2=G_sc2.where(H_Mask==0.0,1.0)

    
        H_ds=H_ds.fillna(0)
        H_ds_re=H_ds*H_sc2+dr_out*G_sc2
        et = time.time()
  # get the execution time
        elapsed_time = et - st
        print('Smoothing  Execution time:', elapsed_time, 'seconds')    
    
        
        print(f'Saving :{Ofile}')
        ncid = nc.Dataset(Ofile, 'a', format='NETCDF4')
        
        ncvar=ncid[key]
        ncvar[:,:,:]=H_ds_re[key].values
        nchis=ncid.history
        ncid.history=f'{nchis}: Appended by Eli Hunter using gridded GFS data '
       #Created by Eli Hunter on {extime} using AWS serivces. {__file__} '
        ncid.close()
        
        
if __name__ == "__main__":
    
    
    print('Running: Combing GFS and HRRR')
    parser = optparse.OptionParser()
    parser.add_option('-d', '--date',dest='date',help='%Y%m%d',default=today,type='str')
    (opts, args) = parser.parse_args()
    main(sys.argv)
