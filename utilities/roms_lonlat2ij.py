import numpy as np
from scipy.interpolate import LinearNDInterpolator
from netCDF4 import Dataset

def roms_lonlat2ij(grd):
    """
    Create interpolants to convert lon/lat to fractional ROMS i,j coordinates
    on the FORTRAN rho‐points grid.

    Parameters
    ----------
    grd : str or dict-like
        - If str: path to a ROMS netCDF file containing 'lon_rho' & 'lat_rho'.
        - If dict-like: must have keys 'lon_rho' and 'lat_rho' as 2D arrays.

    Returns
    -------
    Fi, Fj : LinearNDInterpolator
        Fi(lon, lat) → fractional i-index (0-based)
        Fj(lon, lat) → fractional j-index (0-based)
    """
    # 1) Read lon/lat arrays (and transpose to match MATLAB’s lon_rho')
    if isinstance(grd, str):
        ds  = Dataset(grd)
        lon = ds.variables['lon_rho'][:].T
        lat = ds.variables['lat_rho'][:].T
        ds.close()
    else:
        lon = np.asarray(grd['lon_rho']).T
        lat = np.asarray(grd['lat_rho']).T

    # 2) Build a mesh of i,j indices (0-based)
    M, N = lon.shape
    I2, J2 = np.meshgrid(np.arange(M), np.arange(N), indexing='ij')
    #    I2[i,j] == i,   J2[i,j] == j

    # 3) Flatten lon/lat & corresponding i/j
    pts     = np.column_stack((lon.ravel(), lat.ravel()))
    i_vals  = I2.ravel()
    j_vals  = J2.ravel()

    # 4) Create scattered interpolants
    Fi = LinearNDInterpolator(pts, i_vals)
    Fj = LinearNDInterpolator(pts, j_vals)

    return Fi, Fj