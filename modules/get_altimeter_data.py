#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import glob
import pandas as pd
from statsmodels.nonparametric.smoothers_lowess import lowess
#------------------------------------------------------------------------
def read_rads2(input_file, *args):
    """
    Read a RADS .asc file and return
      data, dac, lat, lon, time, (and optional extra fields by name).

    Positional args (strings) request extra outputs in the order you list them:
      e.g. read_rads2('file.asc', 'satellite', 'phase', 'cycle')
    """
    
    extras = {};
    extras['col_names'] = []

    hl = 0
    with open(input_file, 'r') as f:
        while hl < 1000:
            pos = f.tell()
            line = f.readline()
            if not line or not line.startswith('#'):
                # rewind so data read sees this line
                f.seek(pos)
                break

            hl += 1
            if len(line) >= 12:
                key = line[2:12]
                val = line[14:].strip()
                if   key == 'Satellite ': extras['satellite'] = val
                elif key == 'Phase     ': extras['phase']     = val
                elif key == 'Cycle     ': extras['cycle']     = int(val)
                elif key == 'Pass      ': extras['pass']  = int(val)
                elif key == 'Equ_time  ': extras['equ_time']  = float(line[14:31])
                elif key == 'Equ_lon   ': extras['equ_lon']   = float(val)
                elif key.startswith('Col'): extras['col_names'].append(val)
                # else: quietly skip unknown header lines

                        # --- 2) read numeric data ---
        # from the first non-# line through EOF
        raw = np.loadtxt(f)    # shape (n_rows, n_cols)
    
    if raw.ndim == 1:
        print(f"  {input_file} does not have any data")
        time = [np.datetime64('1990-01-01'), np.datetime64('1990-01-02')]; 
        lat = []; lon = []; data = []; dac = []
    else:
        sec_offsets = raw[:, 0]; lat = raw[:, 1]; lon = raw[:, 2]
        data = raw[:, 3]; dac = raw[:, 7]

        # --- 4) build actual times (as numpy.datetime64) ---
    # MATLAB: datenum(1985,1,1) + sec/86400
    # Here: 1985-01-01 + sec_offsets seconds
        ref = np.datetime64('1985-01-01')
        time = ref + sec_offsets.astype('timedelta64[s]')
        

        base = np.datetime64('2011-01-01')
        time = (time - base) / np.timedelta64(1, 'D') 

    varargout = [extras.get(name.lower(), None) for name in args]
    return (data, dac, lat, lon, time, *varargout)


#-----------------------------------------------------------------------------
def get_rads_ssh(rads_dir):
    """
    Read all RADS .asc sea surface height files in a directory,
    grouping them by pass number and selecting only those records
    whose time overlaps 2011-01-01 onward.

    Returns
    -------
    rads : list of dict
        One dict per unique pass, with keys:
          'raw_ssh', 'raw_dac', 'raw_lat', 'raw_lon', 'raw_time',
          'pass',   'cycle',   'equ_time', 'col_names'
    """
    
    # 2) find all .asc files
    pattern    = os.path.join(rads_dir, '*.asc')
    all_files  = sorted(glob.glob(pattern))
    all_names  = [os.path.basename(f) for f in all_files]
    # MATLAB used f(4:7) to extract the pass number
    all_passes = [int(name[3:7]) for name in all_names]

    # 3) group by unique pass
    uni_passes = sorted(set(all_passes))
    rads       = []

    for p in uni_passes:
        # files for this pass
        files_p = [all_files[i] for i,pp in enumerate(all_passes) if pp == p]

        # prepare per-pass storage
        raw_ssh, raw_dac = [], []
        raw_lat, raw_lon = [], []
        raw_time         = []
        pass_list        = []
        cycle_list       = []
        equ_time_list    = []
        col_names        = None

        for fpath in files_p:
            # read data + metadata fields
            (cur_ssh, cur_dac, cur_lat, cur_lon, cur_time,
             cur_pass, cur_cycle, cur_equ_time, cur_col_names) = \
                read_rads2(fpath, 'pass','cycle','equ_time','col_names')

            # keep only if time overlaps our window
            #if (cur_time[-1] >= sel_start) :
            raw_ssh.append(cur_ssh)
            raw_dac.append(cur_dac)
            raw_lat.append(cur_lat)
            raw_lon.append(cur_lon)
            raw_time.append(cur_time)
            pass_list.append(cur_pass)
            cycle_list.append(cur_cycle)
            equ_time_list.append(cur_equ_time)
            col_names = cur_col_names

        # assemble this pass’s results
        rads.append({
            'raw_ssh'  : raw_ssh,
            'raw_dac'  : raw_dac,
            'raw_lat'  : raw_lat,
            'raw_lon'  : raw_lon,
            'raw_time' : raw_time,
            'pass'     : pass_list,
            'cycle'    : cycle_list,
            'equ_time' : equ_time_list,
            'col_names': col_names
        })

    return rads
#----------------------------------------------------------------------

def get_altimeter_data(grd, alt_dir, start_day):
    """
    Python port of MATLAB get_altimeter_data.
    
    Parameters
    ----------
    grd : any
      (not used here because the MATLAB in‐grid filter is disabled)
    alt_dir : str
      directory containing .asc RADS files
    start_day : int or float
      numeric day count since 2011-01-01 (same units as AllTime)
    
    Returns
    -------
    data : dict of lists
      keys: AllSSHA, AllDAC, AllLon, AllLat, AllTime, AllPass, AllCycle
    """
    # 1) load and group RADS tracks
    rads = get_rads_ssh(alt_dir)
    AllIPass  = []
    AllICycle = []
    AllPass   = []
    AllCycle  = []
    AllTime   = []

    for ipass in range(len(rads)):
        for icycle in range(len(rads[ipass]['cycle'])):
            AllIPass.append(ipass)
            AllICycle.append(icycle)
            AllPass.append(rads[ipass]['pass'][icycle])
            AllCycle.append(rads[ipass]['cycle'][icycle])
            time = np.array(rads[ipass]['raw_time'][icycle], dtype=np.float64)
            mean_day = np.mean(time)
            AllTime.append(mean_day)        

    # to numpy arrays and sort by time
    AllIPass   = np.array(AllIPass)
    AllICycle  = np.array(AllICycle)
    AllPass    = np.array(AllPass)
    AllCycle   = np.array(AllCycle)
    AllTime    = np.array(AllTime)

    order = np.argsort(AllTime)
    AllTime    = AllTime[order]
    AllIPass   = AllIPass[order]
    AllICycle  = AllICycle[order]
    AllPass    = AllPass[order]
    AllCycle   = AllCycle[order]

    # round-to-day bins
    time_round = np.round(AllTime - 0.5).astype(int)
    unique_days = np.unique(time_round)
    unique_days = unique_days[unique_days >= start_day]
    data = {
        'AllSSHA': [], 'AllDAC': [], 'AllLon': [], 'AllLat': [],
        'AllTime': [], 'AllPass': [], 'AllCycle': []
    }

    # 3) loop over each day and each track
    for day in range(unique_days.min(), unique_days.max()+1):
        inds = np.where(time_round == day)[0]
        for idx in inds:
            ip = AllIPass[idx]
            ic = AllICycle[idx]
            rec = rads[ip]

            ssh = rec['raw_ssh'][ic].copy()
            lon = rec['raw_lon'][ic].copy()
            lat = rec['raw_lat'][ic].copy()
            dac = rec['raw_dac'][ic].copy()
            t0  = rec['raw_time'][ic].copy()
            
           # --- 1) smoothing ---
            if len(ssh) >= 8:
                ll = lon + 1j*lat
                # --- 2) reorder along-track by lon→lat ---
                order2 = np.argsort((ll))
                lon = lon[order2]; lat = lat[order2]
                ssh = ssh[order2]; dac = dac[order2]

 
                # moving average, window = 8
                ssh = pd.Series(ssh).rolling(window=8, center=True).mean().to_numpy()  # Moving average
                # loess
                ssh = lowess(ssh, np.arange(len(ssh)), frac=8/len(ssh), return_sorted=False)

 
                # --- 3) remove jumps in track (>0.48 deg) ---
                ll    = lon + 1j*lat
                diffs = np.abs(np.diff(ll))
                breaks = np.where(diffs > 0.48)[0]
                mask = np.ones(len(lon), dtype=bool)
                for b in breaks:
                    i1 = max(0, b-4)
                    i2 = min(b+5, len(lon)-1)
                    mask[i1:i2+1] = False
                lon = lon[mask]; lat = lat[mask]
                ssh = ssh[mask]; dac = dac[mask]

                # --- 4) drop NaN or zero SSH ---
                good = ~(np.isnan(ssh) | (ssh == 0))
                lon = lon[good]; lat = lat[good]
                ssh = ssh[good]; dac = dac[good]

            # --- 5) filter known bad regions ---
                regions = [
                    (-66.4, -64.5, 44.1, 45.9),   # Bay of Fundy
                    (-70.45,-70.1,41.8, 42.0),    # Grand Banks
                    (-73.5, -71.8,40.9, 41.3),    # Long Island Sound
                    (-75.5, -74.9,38.7, 39.6),    # Delaware Bay
                    (-75.9, -75.65,36.25,36.65)   # NC Coast
                ]
                for lonmin, lonmax, latmin, latmax in regions:
                    m2 = ~((lon>lonmin)&(lon<lonmax)&(lat>latmin)&(lat<latmax))
                    lon = lon[m2]; lat = lat[m2]
                    ssh = ssh[m2]; dac = dac[m2]

                # --- 6) append to output lists ---
                data['AllSSHA'].extend(ssh.tolist())
                data['AllDAC'].extend(dac.tolist())
                data['AllLon'].extend(lon.tolist())
                data['AllLat'].extend(lat.tolist())
                mean_day = AllTime[idx]
                data['AllTime'].extend([mean_day]*len(ssh))
                data['AllPass'].extend([int(AllPass[idx])]*len(ssh))
                data['AllCycle'].extend([int(AllCycle[idx])]*len(ssh))
            else:
                continue


    return data
#-------------------------------------------------------------------------
def get_all_altimetry(grd, start_day,fconfig):

    # --- Get observations
    # ————————————————
    # 1) DEFINE YOUR PATHS
    # ————————————————
    jas3_dir        = fconfig['obs']['romsobs']['SSH']['jas3_dir']
    sentinel3a_dir  = fconfig['obs']['romsobs']['SSH']['sentinel3a_dir']
    sentinel3b_dir  = fconfig['obs']['romsobs']['SSH']['sentinel3b_dir']
    sentinel6a_dir  = fconfig['obs']['romsobs']['SSH']['sentinel6a_dir']
    swot_dir        = fconfig['obs']['romsobs']['SSH']['swot_dir']

    altimeter_data = {
        'time':       [],
        'lon':        [],
        'lat':        [],
        'ssha':       [],
        'pass':       [],
        'cycle':      [],
        'provenance': [],
        'dac':        []
    }
    print("processing Jason 3 data")
    data = get_altimeter_data(grd, jas3_dir, start_day - 6)
    altimeter_data['time'].extend(data['AllTime'])
    altimeter_data['lon'].extend(data['AllLon'])
    altimeter_data['lat'].extend(data['AllLat'])
    altimeter_data['ssha'].extend(data['AllSSHA'])
    altimeter_data['dac'].extend(data['AllDAC'])
    altimeter_data['pass'].extend(data['AllPass'])
    altimeter_data['cycle'].extend(data['AllCycle'])
    altimeter_data['provenance'].extend([403] * len(data['AllSSHA']))

    print("processing Sentinel3a data")
    data = get_altimeter_data(grd, sentinel3a_dir, start_day - 6)
    altimeter_data['time'].extend(data['AllTime'])
    altimeter_data['lon'].extend(data['AllLon'])
    altimeter_data['lat'].extend(data['AllLat'])
    altimeter_data['ssha'].extend(data['AllSSHA'])
    altimeter_data['dac'].extend(data['AllDAC'])
    altimeter_data['pass'].extend(data['AllPass'])
    altimeter_data['cycle'].extend(data['AllCycle'])
    altimeter_data['provenance'].extend([441] * len(data['AllSSHA']))

    print("processing Sentinel3b data")
    data = get_altimeter_data(grd, sentinel3b_dir, start_day - 6)
    altimeter_data['time'].extend(data['AllTime'])
    altimeter_data['lon'].extend(data['AllLon'])
    altimeter_data['lat'].extend(data['AllLat'])
    altimeter_data['ssha'].extend(data['AllSSHA'])
    altimeter_data['dac'].extend(data['AllDAC'])
    altimeter_data['pass'].extend(data['AllPass'])
    altimeter_data['cycle'].extend(data['AllCycle'])
    altimeter_data['provenance'].extend([442] * len(data['AllSSHA']))

    print("processing Sentinel6a data")
    data = get_altimeter_data(grd, sentinel6a_dir, start_day - 6)
    altimeter_data['time'].extend(data['AllTime'])
    altimeter_data['lon'].extend(data['AllLon'])
    altimeter_data['lat'].extend(data['AllLat'])
    altimeter_data['ssha'].extend(data['AllSSHA'])
    altimeter_data['dac'].extend(data['AllDAC'])
    altimeter_data['pass'].extend(data['AllPass'])
    altimeter_data['cycle'].extend(data['AllCycle'])
    altimeter_data['provenance'].extend([404] * len(data['AllSSHA']))
    print(min(altimeter_data['time']))
    print(max(altimeter_data['time']))
    np.savez_compressed(fconfig['obs']['romsobs']['SSH']['altfile'], **altimeter_data)

    return altimeter_data

