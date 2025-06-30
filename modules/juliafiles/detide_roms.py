from netCDF4 import Dataset, num2date
from utide import solve
import numpy as np

#  Open the NetCDF file
nc = Dataset('/home/julia/ROMS/eccofs/OBS/SSH/eccofs6km_xtr.nc')

#  Read coordinate/time variables
#  
time = nc.variables['ocean_time']
times = num2date(time[:], units=time.units)

#  Read the data variable of interest (e.g. “temp”)
lat_tide  = nc.variables['lat_rho'][:]
mask_tide  = nc.variables['mask_rho'][:]

L,M = lat_tide.shape

#data = np.load('tidal_coef_6km.npz', allow_pickle=True)
#coef = data['coef']
coef = np.empty((L, M), dtype=object)
coef[:, :] = np.nan
for i in range(L):
    for j in range(M):
        # compute harmonic coefficients
        if mask_tide[i,j]:
            zeta_tide = nc.variables['zeta'][:,i,j]     
            zeta_tide = zeta_tide - np.nanmean(zeta_tide) 
            coef[i,j] = solve(
                    times,               # array of datetimes
                    zeta_tide,               # demeaned data
                    lat=lat_tide[i,j],    # e.g. 41.5
                    constit=['M2','S2','N2','K2','K1','O1','P1','Q1','S1','M4','MN4','MS4'],
                    nodal=True,      # apply nodal corrections (default)
                    trend=False,     # don’t fit a linear trend
                    method='ols',    # ordinary least squares
                    conf_int='MC',   # Monte Carlo confidence intervals
                    verbose=False,
                    Rayleigh_min=0.95  # drop constituents with low separability
                )
    print(f"-----------------Saving to file after processing row  {i}----------")
    np.savez_compressed('tidal_coef_6km.npz', coef = coef)
