#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 10:26:35 2025

@author: julia
"""

import os,sys
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
import warnings
import shutil

sys.path.append('/Users/julia/python')
sys.path.append('/home/julia/python')
from roms_get_grid import roms_get_grid
from roms_lonlat2ij import roms_lonlat2ij
from define_obs_file import define_4dvar_obs_file, obs_provenance_definition_eccofs
from add_history import add_history
from accumarrays import accum2d

def combine_obs():
    OBS_out = "/p/julia/ROMS/eccofs/OBS_python/CombineOBS/OBS_FILES/ssh_sst_cmems_" 
    grd_file = "/p/julia/ROMS/eccofs/Data/grid_eccofs_6km_09_b7.nc"
    grd_file_3km = "/p/julia/ROMS/eccofs/Data/grid_eccofs_3km_08_b7.nc"

    amsr2_file = "/p/julia/ROMS/eccofs/OBS_python/AMSR2/OBS_FILES/sst_amsr2_nrt_"
    leo_file = "/p/julia/ROMS/eccofs/OBS_python/LEO/OBS_FILES/sst_leo_nrt_"
    goes_file = "/p/julia/ROMS/eccofs/OBS_python/GOES1Hour/OBS_FILES/sst_goes_nrt_"
    ssh_file = "/p/julia/ROMS/eccofs/OBS_python/SSH/OBS_FILES_1h_16_3kmtide/along_track_TotalSSH_"
    cmems_file = "/p/julia/ROMS/eccofs/OBS_python/CMEMS/OBS_FILES/cmems_nrt_"

    m1 = np.load('/p/julia/ROMS/eccofs/OBS_python/CombineOBS/sst_statistics.npz')

    scoord = [7, 2, 250, 50, 2, 4]
    g = roms_get_grid(grd_file, scoord)
    g3 = roms_get_grid(grd_file_3km,scoord)
    [Fi,Fj] = roms_lonlat2ij(g3);

    end_day = pd.Timestamp.today().normalize() - pd.Timedelta(days=1)
    start_day     = end_day - pd.Timedelta(days=2)

    days = pd.date_range(start_day, end_day, freq='3D')
    ref_datum = pd.Timestamp(2011,1,1);           # reference datum in ROMS                                                                                             

    #------------------------------------------------------
    #      Process all satellite SST
    #------------------------------------------------------
    # parameters for binning SST
    dTime = 6/24; # 6 hours
    t_epsilon = 5/24/60;  # slightly more than 4 minutes which is 1*dt in eccofs 6km
    
    std_harm = m1['std_harm']
    w = 2*np.pi/365.25;
    
    all_pts = np.column_stack((g['lon_rho'].ravel(), g['lat_rho'].ravel()))
    Fm = NearestNDInterpolator(all_pts, g['mask_rho'].ravel())

# %%

    for day_pd in days:
        day = (day_pd - ref_datum).days
        print(f"=================================Merging observations for {day_pd}==============")
 
        expected_deviation = np.squeeze(
            std_harm[0, :, :] +
            std_harm[1, :, :] * np.cos(w * day) +
            std_harm[2, :, :] * np.sin(w * day) +
            std_harm[3, :, :] * np.cos(2 * w * day) +
            std_harm[4, :, :] * np.sin(2 * w * day)
            )
        
        # specify files for merging
        files = []
        files.extend([f"{amsr2_file}{int(day):04d}.nc", f"{leo_file}{int(day):04d}.nc", f"{goes_file}{int(day):04d}.nc"])
        files.extend([f"{amsr2_file}{int(day+1):04d}.nc", f"{leo_file}{int(day+1):04d}.nc", f"{goes_file}{int(day+1):04d}.nc"])
        files.extend([f"{amsr2_file}{int(day+2):04d}.nc", f"{leo_file}{int(day+2):04d}.nc", f"{goes_file}{int(day+2):04d}.nc"])

        #initiate arrays
        SST = {}
        SST['obs_lon'] = []; SST['obs_lon'] = []; SST['obs_lat'] = []; 
        SST['depth'] = []; SST['obs_value'] = []; SST['obs_Xgrid'] = []; 
        SST['obs_Ygrid'] = []; SST['obs_Zgrid'] = []; SST['obs_provenance']=[];
        SST['obs_time'] = []; SST['obs_depth']=[]; SST['obs_label']=[]; SST['obs_type']=[]

        for file in files:
            if os.path.isfile(file):
                print(f"reading data from {file}")
                with Dataset(file) as nc:
                    for var in SST.keys():
                        SST[var] = np.concatenate([SST[var],nc.variables[var][:]])
        # %%
        # bin size (approx. 25 km at the equator)
        if (len(SST['obs_lon'])>0):
            dLat = 0.25; dLon = 0.25
            minLat = np.min(SST['obs_lat']); minLon = np.min(SST['obs_lon'])

            latBin = np.floor((SST['obs_lat'] - minLat) / dLat).astype(int) 
            lonBin = np.floor((SST['obs_lon']  - minLon) / dLon).astype(int) 

            # Dimensions
            latDimLength = np.max(latBin)+1
            lonDimLength = np.max(lonBin)+1
 
            # Combine 2D indices into a 1D index (like MATLAB's sub2ind)
            varInd = np.ravel_multi_index((latBin, lonBin), dims=(latDimLength, lonDimLength))

            # Ones column and max index
            onesCol = np.ones_like(varInd)
            maxVarInd = np.max(varInd)+1
        
        # Bin at 25 km
            counter = accum2d(varInd,0,onesCol,shape=[maxVarInd, 1],func=np.sum).ravel()
            lat_b = accum2d(varInd,0,SST['obs_lat'],shape=[maxVarInd, 1],func=np.mean).ravel() 
            lon_b = accum2d(varInd,0,SST['obs_lon'],shape=[maxVarInd, 1],func=np.mean).ravel()
            sst_b = accum2d(varInd,0,SST['obs_value'],shape=[maxVarInd, 1],func=np.mean).ravel()
            std_b = accum2d(varInd,0,SST['obs_value'],shape=[maxVarInd, 1],func=np.std).ravel() 
        
            # compute standard deviation at lon_b and lat_b from model std
            mask = np.where(~np.isnan(expected_deviation))
            pts = np.column_stack((g['lon_rho'][mask], g['lat_rho'][mask]))
            vals = expected_deviation[mask]
            Fstd      = LinearNDInterpolator(pts, vals)
        
            indx = np.where(counter>0)[0]
            estd = np.full_like(counter, np.nan, dtype=float)
            estd[indx] = Fstd(lon_b[indx],lat_b[indx] )
        
            # quality control on SST
            ss = np.where((counter<2) | ((counter==2) & (std_b > estd)))[0]

            for s in ss:
                iii = np.where(varInd == s)[0]
                SST['obs_value'][iii] = np.nan

            ss = np.where((counter>2) & (std_b > 1.5 * estd))[0]

            for s in ss:
                iii = np.where(varInd == s)[0]
                iv = np.where(np.abs(SST['obs_value'][iii] - sst_b[s]) > 2.5 * estd[s])[0]
                SST['obs_value'][iii[iv]] = np.nan
            
            #  now  bin onto 100 km grid and remove values more than 5 std away from the
            # mean
            dLat = 1; dLon = 1
            latBin = np.floor((SST['obs_lat'] - minLat) / dLat).astype(int) 
            lonBin = np.floor((SST['obs_lon']  - minLon) / dLon).astype(int) 

            # Dimensions
            latDimLength = np.max(latBin)+1
            lonDimLength = np.max(lonBin)+1
 
            # Combine 2D indices into a 1D index (like MATLAB's sub2ind)
            varInd = np.ravel_multi_index((latBin, lonBin), dims=(latDimLength, lonDimLength))

            # Ones column and max index
            onesCol = np.ones_like(varInd)
            maxVarInd = np.max(varInd)+1
       
            # Bin at 25 km
            counter = accum2d(varInd,0,onesCol,shape=[maxVarInd, 1],func=np.sum).ravel()
            lat_b = accum2d(varInd,0,SST['obs_lat'],shape=[maxVarInd, 1],func=np.mean).ravel() 
            lon_b = accum2d(varInd,0,SST['obs_lon'],shape=[maxVarInd, 1],func=np.mean).ravel()
            sst_b = accum2d(varInd,0,SST['obs_value'],shape=[maxVarInd, 1],func=np.mean).ravel()
            std_b = accum2d(varInd,0,SST['obs_value'],shape=[maxVarInd, 1],func=np.std).ravel() 
 
            indx = np.where(counter>0)[0]
            estd = np.full_like(counter, np.nan, dtype=float)
            estd[indx] = Fstd(lon_b[indx],lat_b[indx] )
        
            ss = np.where((counter>2) & (std_b > 1.5 * estd))[0]

            for s in ss:
                iii = np.where(varInd == s)[0]
                iv = np.where(np.abs(SST['obs_value'][iii] - sst_b[s]) > 5 * estd[s])[0]
                SST['obs_value'][iii[iv]] = np.nan

            # remove tagged observations
            good = np.where(~np.isnan(SST['obs_value']))[0]

            for var in SST.keys():
                SST[var] = SST[var][good];
        
            # design weights for binning
            weights = np.ones_like(good, dtype=float)
            weights[SST['obs_provenance'] == 317] = 0.05  # GOES
            weights[SST['obs_provenance'] == 311] = 10.0  # LEO

            # bin in space and time
            minTime = np.min(SST['obs_time'])

            xind = np.floor(SST['obs_Xgrid']).astype(int); yind = np.floor(SST['obs_Xgrid']).astype(int)
            timeBin = np.floor((SST['obs_time']- minTime) / dTime).astype(int)
            timeDimLength  = np.max(timeBin)+1;
            # linear index of each (xbin, ybin)
            xdimLength = np.max(xind)+1; ydimLength = np.max(yind)+1
            varInd = np.ravel_multi_index((xind, yind), dims=(xdimLength, ydimLength))
            maxVarInd = np.max(varInd)+1;
        
            ind = np.where(~np.isnan(SST['obs_value']))[0];
            counter = accum2d(timeBin[ind], varInd[ind],np.ones_like(ind),shape=[timeDimLength, maxVarInd],func=np.sum)
            good = ~np.isnan(counter)
        
            var_names = ['obs_lat', 'obs_lon', 'obs_Xgrid', 'obs_Ygrid']
            data = {}
            for var in var_names:
                data[var] = accum2d(timeBin[ind], varInd[ind],SST[var][ind],shape=[timeDimLength, maxVarInd],func=np.mean)[good]

            var = 'obs_value'        
            data[var] = accum2d(timeBin[ind], varInd[ind],SST[var][ind]*weights[ind],shape=[timeDimLength, maxVarInd],func=np.sum)[good] 
            data[var] = data[var] / accum2d(timeBin[ind], varInd[ind],weights[ind],shape=[timeDimLength, maxVarInd],func=np.sum)[good]

            time = accum2d(timeBin[ind], 0, SST['obs_time'][ind],[timeDimLength, 1],np.mean).ravel()
 
            # replicate across space dimension
            data['obs_time'] = np.tile(time[:, None], (1, maxVarInd))[good]

            # compute obs_label that contains all the provenances that go into superobs
            prov = np.unique(SST['obs_provenance'])
            Lp = prov.size

            # multipliers: 10**[0,1,2,...]
            mult = 10 ** np.arange(Lp, dtype=int)

            # flag array same shape as provenance_t
            flag = np.zeros_like(SST['obs_provenance'], dtype=int)
            for ix, p in enumerate(prov):
                flag[SST['obs_provenance'] == p] = mult[ix]

            # accumulate flags into super‐obs using the same binning, get sum of multipliers
            flag_b   = accum2d(timeBin[ind], varInd[ind], flag[ind], shape=(timeDimLength, maxVarInd),func=np.mean)[good]
    
            # decode flag_b into labels
            # For each provenance level, extract its digit in flag_b
            flag_ind = np.zeros((flag_b.size, Lp), dtype=int)
            for ix, p in enumerate(prov):
                # digit at place idx is floor(flag_b/mult[idx]) % 10, capped at 1
                digit = (flag_b // mult[ix]) % 10
                flag_ind[:, ix] = p * np.minimum(1, digit)

            # identify rows where sum across prov‐levels is zero (no obs)
            no_obs = np.sum(flag_ind, axis=1) == 0

            # build string labels with concatenated digits, replacing zero‐rows with '0'
            labels = []
            for i, row in enumerate(flag_ind):
                if no_obs[i]:
                    labels.append('0')
                else:
                    # join nonzero digits into string
                    s = ''.join(str(d) for d in row if d != 0)
                    labels.append(s or '0')

            # convert labels back to integers
            data['obs_label'] = np.array([int(s) for s in labels], dtype=int)

            # provenance_b: same as label_b but any label > 999 means that more than one provenance contributed, identifies as  400
            data['obs_provenance'] = data['obs_label'].copy()
            data['obs_provenance'][data['obs_label'] > 999] = 300

# %%
        # add all other arrays

            LL = len(data['obs_lon'])
            data['depth'] = np.zeros(LL)
            data['obs_depth'] = g['N'] * np.ones(LL)
            data['obs_Zgrid'] = data['obs_depth']
            data['obs_type'] = 6 * np.ones(LL)
        
            # fix SST error
            data['obs_error'] = np.ones(LL)
            ind = np.where(data['obs_provenance']==317)[0]  # geostationary only
            data['obs_error'][ind] = 0.9**2 * np.ones_like(ind)
        
            # microwave and all microwave superobs
            ind = np.where((data['obs_provenance'] == 324) | ( (data['obs_provenance'] == 300) & (data['obs_label'] % 10**3==324)))[0] 
            data['obs_error'][ind] = 0.7**2 * np.ones_like(ind)
        
            # infrared and all infrared superobs
            ind = np.where((data['obs_provenance'] == 311) | ( (data['obs_provenance'] == 300) & (data['obs_label'] % 10**3==311)) |
                       ( (data['obs_provenance'] == 300) & (data['obs_label'] % 10**6==311)))[0] 
            data['obs_error'][ind] = 0.6**2 * np.ones_like(ind)
        else:
            data = {}
            for var in SST.keys():
                data[var] = []
 
            
# %%
        
        #---------------------------------------------------------
        #     Add all other observations
        #---------------------------------------------------------
        files = []
        files.extend([f"{ssh_file}{int(day-1):04d}.nc", f"{cmems_file}{int(day-1):04d}.nc"])
        files.extend([f"{ssh_file}{int(day):04d}.nc", f"{cmems_file}{int(day):04d}.nc"])
        files.extend([f"{ssh_file}{int(day+1):04d}.nc", f"{cmems_file}{int(day+1):04d}.nc"])
        files.extend([f"{ssh_file}{int(day+2):04d}.nc", f"{cmems_file}{int(day+2):04d}.nc"])

        for file in files:
            if os.path.isfile(file):
                print(f"reading data from {file}")
                with Dataset(file) as nc:
                    for var in data.keys():
                        data[var] = np.concatenate([data[var],nc.variables[var][:]])
                        
        if (len(data['obs_lon']>0)) :
            # remove values that fall outside of time window (this occasionally happens
            # when doing repeat of of observations with lags, and fall in the water
            good = np.where((data['obs_time']>=day) & (data['obs_time']<day+3) &
                        (Fm(data['obs_lon'], data['obs_lat']) == 1))[0];
            for var in data.keys():
                data[var] = data[var][good]

            ind = np.where(data['obs_type']==1)[0]  # altimetry
            data['obs_error'][ind] = 0.04**2 * np.ones_like(ind)

            # scale error for CMEMS Argo profiles (801), Gliders (806), 
            # CTD profiles (807), Sea Mammals (808), Moorings (823), 
            # thermosalinograph (824)
        
            ind = np.where((data['obs_type']==6) & (data['obs_provenance'] == 823))[0]  # Temp from Moorings
            data['obs_error'][ind] = 0.9**2 * data['obs_error'][ind]

            ind = np.where((data['obs_type']==7) & (data['obs_provenance'] == 823))[0]  # Salt from Moorings
            data['obs_error'][ind] = 0.6**2 * data['obs_error'][ind]

            ind = np.where((data['obs_type']==6) & (data['obs_provenance'] == 807))[0]  # Temp from CTD
            data['obs_error'][ind] = 0.9**2 * data['obs_error'][ind]                 
        
            ind = np.where((data['obs_type']==7) & (data['obs_provenance'] == 807))[0]  # Salt from CTD
            data['obs_error'][ind] = 0.6**2 * data['obs_error'][ind]

            ind = np.where((data['obs_type']==6) & (data['obs_provenance'] == 806))[0]  # Temp from Gliders
            data['obs_error'][ind] = 0.9**2 * data['obs_error'][ind]                 

            ind = np.where((data['obs_type']==7) & (data['obs_provenance'] == 806))[0]  # Salt from Gliders
            data['obs_error'][ind] = 0.6**2 * data['obs_error'][ind]                 

            data['obs_depth'][data['obs_depth']==0] = g['N']
            data['obs_Zgrid'] = data['obs_depth']
        
# %%
        
        # Since we have only one source of in-situ observations, we do not
        # need to identify identical observations like we did in doppio
        # may need to add this when we'd have more than one source
        
        # --------------------------------------------------------
        #      merge into surveys
        #---------------------------------------------------------
            data['obs_time'] = np.round(data['obs_time'], 8)
            survey_time, It, Is = np.unique(data['obs_time'], return_index=True, return_inverse=True)
            # Prepare working arrays                                                                                                                               
            time_tmp = survey_time.copy()
            time_new = np.zeros_like(survey_time)

            #  Merge “close” survey times within t_epsilon                                                                                                         
            for j in range(len(survey_time)):
                if not np.isnan(time_tmp[j]):
                    # find all survey_time within t_epsilon of survey_time[j]                                                                                      
                    close = np.where(np.abs(time_tmp - time_tmp[j]).astype(np.float32) < t_epsilon)[0]
                    # assign the master time to all close entries                                                                                                  
                    time_new[close] = time_tmp[j]
                    # update obs_time for every original observation in each close group                                                                           
                    # report if we merged more than one                                                                                                            
                    if close.size > 1:
                        for idx in close:
                            data['obs_time'][Is == idx] = time_new[idx]

                        merged = survey_time[close]
                        print(f"   merging surveys at times {merged} ")
                        # mark these as done                                                                                                                       
                        time_tmp[close] = np.nan

            # sort the data first by time then by type
            tmp = 1.e9*data['obs_time'] + data['obs_type'];
            ind = np.argsort(tmp)
            for var in data.keys():
                data[var] = data[var][ind];


            # Recompute unique survey times after merging                                                                                                          
            survey_time, It, Is = np.unique( data['obs_time'], return_index=True, return_inverse=True)
            
            # Count how many observations per survey time                                                                                                          
            #    (Is contains the survey index for each obs_time entry)                                                                                            
            Nobs = np.bincount(Is, minlength=survey_time.size)
# %%
            output_fname = f"{OBS_out}{day:04d}.nc"
            survey = len(survey_time)
            define_4dvar_obs_file(output_fname,survey,obs_provenance_definition_eccofs())
            ds = Dataset(output_fname, 'a')
            ds.history = "CMEMS near real time data transformed into 4DVAR format"
            ds.variables['spherical'][0] = 'T'
            ds.variables['Nobs'][:] = Nobs
            ds.variables['survey_time'][:] = survey_time
            ds.variables['obs_variance'][:] = np.ones(7)

            for var in data.keys():
               ds.variables[var][:] = data[var]
            ds.close()
            
            addhistory = (
            "Observations of altimeter, satellite SST and in-situ T and S"
            " merged into one file for assimilation"   
            "Prepared by Julia Levin. "
            )
            add_history(output_fname, addhistory)

          
            # create observation file for 3km grid
            output_fname_3km = f"{OBS_out}{day:04d}_3km.nc"
            # Copy file
            shutil.copy2(output_fname, output_fname_3km)  
            Xgrid = Fi(data['obs_lon'], data['obs_lat']);
            Ygrid = Fj(data['obs_lon'], data['obs_lat']);

            with Dataset(output_fname_3km, 'r+') as nc:
                nc.variables['obs_Xgrid'][:] = Xgrid
                nc.variables['obs_Ygrid'][:] = Ygrid
        else:
            print(f"  no observations for {day_pd}==============")
     
                
#-----------------------------------------------------------------------

if __name__ == "__main__":
   warnings.filterwarnings("ignore")
   combine_obs()
