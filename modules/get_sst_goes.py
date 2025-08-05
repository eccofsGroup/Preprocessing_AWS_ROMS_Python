#!/usr/bin/env python
# coding: utf-8

import sys,os,glob
import pandas as pd
import netCDF4 as nc
import numpy as np

from roms_get_grid import roms_get_grid
from CreateObsFile_SST import CreateObsFile_SST
from add_history import add_history
from accumarrays import accum3d, accum2d

def main(fconfig):

    #base_file_name = './OBS_FILES/sst_goes_nrt_';
    #url = '/p/om/cron/ECCOFS_OBS/PO.DAAC/data/raw/GOES16_L3/';
    
    url = fconfig['obs']['romsobs']['GOES']['URL']
    working_directory =fconfig['obs']['romsobs']['GOES']['WORKING_DIRECTORY']
    prefix = fconfig['obs']['romsobs']['GOES']['prefix']
    output_file  = os.path.join(working_directory, prefix)
    ndays=fconfig['obs']['romsobs']['ndays']
    ref_datum = pd.Timestamp(2011,1,1);           # reference datum in ROMS
    base_file_name = os.path.join(working_directory, prefix)
    end_day = pd.Timestamp.today().normalize() - pd.Timedelta(days=1)
    start_day     = end_day - pd.Timedelta(days=ndays)
    all_files = glob.glob(f"{base_file_name}*.nc")
    all_files = [os.path.basename(f) for f in all_files]

    if all_files:
           # extract the numeric timestamp from filenames
           d = len(prefix)
           times = [int(fn[d:-3]) for fn in all_files]  # strip prefix & ".nc"          for fn in all_files:
           times.sort()

           if (start_day> ref_datum + pd.Timedelta(days=times[-1])):
                 start_day = ref_datum + pd.Timedelta(days=times[-1])

    days = pd.date_range(start_day, end_day, freq='D')
    ref_datum_goes = pd.Timestamp(1981,1,1);       # reference datum in LEO
    mintime = (ref_datum_goes - ref_datum).days    
    dTime      = 6/24 #6 hours


    grd_file =fconfig['obs']['romsobs']['gridfile']
    scoord = fconfig['obs']['romsobs']['scoord']
    g = roms_get_grid(grd_file, scoord)

    L, M = g['lon_rho'].shape
    IC = np.unique(np.concatenate((np.arange(0, L, 4), [L - 1])))
    JC = np.unique(np.concatenate((np.arange(0, M, 4), [M - 1])))
    
    # Create the new grid by selecting the decimated indices.
    # np.ix_ creates an open mesh from the index vectors.
    g1 = {}
    g1['lon_rho'] = g['lon_rho'][np.ix_(IC, JC)]
    g1['lat_rho'] = g['lat_rho'][np.ix_(IC, JC)]
    g1['mask_rho_nan'] = g['mask_rho_nan'][np.ix_(IC, JC)]
    g1['IC'] = IC
    g1['JC'] = JC

    file = fconfig['obs']['romsobs']['GOES']['LLFILE']
    with nc.Dataset(file) as ds:
        lon = ds.variables['lon'][:]  # 1D array of longitudes
        lat = ds.variables['lat'][:]  # 1D array of latitudes
    lond, latd = np.meshgrid(lon, lat)

    # Assuming g1.lon_rho is a 2D array with shape (L, M)
    L, M = g1['lon_rho'].shape
    
    
    # Create grouping arrays, filled with NaN
    latgroup = np.zeros((len(lat), len(lon)), dtype=int)-1
    # latgroup = np.empty((len(lat), len(lon)), dtype=int)
    longroup = latgroup.copy()
    tol = 4 * 0.06  
  #  print([latgroup.min(),latgroup.max()])
  #  print([longroup.min(),longroup.max()])
    for i in range(L):
        for j in range(M):
            # find indices in lon and lat within tolerance of this grid point
            ind_lo = np.where(np.abs(lon - g1['lon_rho'][i, j]) < tol)[0]
            ind_la = np.where(np.abs(lat - g1['lat_rho'][i, j]) < tol)[0]
            if ind_lo.size and ind_la.size:
                
                
                # assign (i+1),(j+1) to match MATLAB's 1-based indexing
                latgroup[np.ix_(ind_la, ind_lo)]  = i 
                longroup[np.ix_(ind_la, ind_lo)] = j 
                #print([i,j])
                #print(latgroup[np.ix_(ind_la, ind_lo)] )
                #print(longroup[np.ix_(ind_la, ind_lo)] )
    M1,N1 = longroup.shape
    #print([M1,N1])
    #mask = (longroup > 0) & (longroup <= M)
    #count = np.sum(mask)
   # print(f"Number of values between {0} and {M}: {count}")
   # print([latgroup.min(),latgroup.max()])
   # print([longroup.min(),longroup.max()])
    for day in days:
        datestr = day.strftime('%Y%m%d')
        print(f"Processing day {datestr}")
        fname_in = f"{url}{datestr}{fconfig['obs']['romsobs']['GOES']['suffix']}"

        if not os.path.exists(fname_in):
            print("  file not found, skipping.")
            continue
        # --- READ & QC ------------------------------------------------
        with nc.Dataset(fname_in) as ds:
#            tmp = ds.variables['sea_surface_temperature'][:] - 273.15
            sst = ds.variables['sea_surface_temperature'][:] - 273.15 - ds.variables['sses_bias'][:]
            # quality flags
            qflag = ds.variables['quality_level'][:]

            # convert 'time' var (s → days since 2011‑01‑01, ROMS reference datum)
            dtime = ds.variables['time'][:] / 86400.0 + mintime 

        # quality control mask
        bad = ( (np.abs(sst) > 100) | np.isnan(sst) | (qflag < 5) )
        sst[bad]   = np.nan; 
        # Compute bin indices
        mindays = (day - ref_datum).days

        K = dtime.shape[0]; M1,N1 = longroup.shape
        timeBin    = np.floor((dtime - mindays) / dTime).astype(int)  # 1d array of indices
        t_bins = int(1/dTime)

        # replicate the 2D group arrays along the time axis to make 3D xBin, yBin
        xBin = np.tile(latgroup[np.newaxis, :, :], (K, 1, 1))  # shape (K, L, M)
        yBin = np.tile(longroup[np.newaxis, :, :], (K, 1, 1))   # shape (K, L, M)
        tBin = np.tile(timeBin[:, None, None], (1, M1, N1))


        good = (~np.isnan(dtime)) 
        # binning times (1d)                                                                                                                       
        time=accum2d(timeBin[good], 0, dtime[good],shape=(t_bins, 1),func=np.mean).ravel()

        # binning in 3d                                                                
        good = (~np.isnan(sst)) & (xBin>0) & (yBin>0) 
        # print(type(good))
        # print(t_bins)
        # print(L)
        # print(M)
        # print(tBin.shape)
        # print([tBin.min(),tBin.max()])
        # print(xBin.shape)
        # print([xBin.min(),xBin.max()])
        # print(yBin.shape)
        # print([yBin.min(),yBin.max()])
        # print(sst.shape)
        # print([sst.min(),sst.max()])
        
        
        sst_mean = accum3d(tBin[good], xBin[good], yBin[good], sst[good], shape=(t_bins, L, M), func=np.mean)
        #sst_std = accum3d(tBin[good], xBin[good], yBin[good], sst[good], shape=(t_bins, L, M), func="std")

        # apply ROMS mask
        mask3d    = np.tile(g1['mask_rho_nan'][np.newaxis,:,:], (t_bins,1,1))
        sst_mean  = sst_mean * mask3d
    
        # create observation file
        fname   = f"{base_file_name}{mindays:04d}.nc"
        
        flag = CreateObsFile_SST(fname, sst_mean, g1, time, 317)
        if flag != 1:
            print(f"  {fname} not created (no data).")
        else:
           add_history(fname, 'GOES19 L3 SST', prepend=True)
           add_history(fname, 'Prepared by Julia Levin (julia@marine.rutgers.edu)', prepend=False)
# if __name__ == "__main__":
#     get_sst_goes()
