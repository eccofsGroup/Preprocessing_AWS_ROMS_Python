import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from define_eccofs_4dvar_obs_file import define_4dvar_obs_file

def nc_addhist(ncfile, attval):
    """
    Add text to the global history attribute of a NetCDF file.
    
    The new history entry is prepended to any existing history,
    and a date/time stamp is automatically inserted.
    
    Parameters
    ----------
    ncfile : str
        Path to the NetCDF file.
    attval : str
        The text to be added to the global history attribute.
    
    Raises
    ------
    ValueError
        If attval is not a string.
    """
    if not isinstance(attval, str):
        raise ValueError("The history attribute value must be a string.")
    
    # Open the NetCDF file in append mode.
    with Dataset(ncfile, 'a') as ds:
        # Try to get the existing "history" attribute; if it doesn't exist, use an empty string.
        try:
            old_hist = ds.getncattr('history')
        except AttributeError:
            old_hist = ''
        
        # Create a timestamp string.
        now_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        if not old_hist:
            new_history = f"{now_str}:  {attval}"
        else:
            new_history = f"{now_str}:  {attval}\n{old_hist}"
        
        # Set the global "history" attribute with the new combined history.
        ds.setncattr('history', new_history)
        
def write_to_netcdf(nc_filename, varname, data):
    # Open the NetCDF file in append ('a') mode.
    with Dataset(nc_filename, 'a') as ds:
        # Write the provided data to the variable.
        ds.variables[varname][:] = data

def CreateObsFile_SST(fname, sst, grd, time, Provenance):
    """
    Create a NetCDF observation file of AMSR2 SST using the provided sst, grid, and time data.

    Parameters
    ----------
    fname : str
        The output filename for the NetCDF file.
    sst : numpy.ndarray
        3D array containing sea surface temperature values. Expected shape:
        (npasses, Ny, Nx) where npasses is the number of survey records.
    grd : dict
        A dictionary containing grid information, for example:
          - 'lon_rho': 2D array of longitudes,
          - 'lat_rho': 2D array of latitudes,
          - 'IC': array-like of indices (1-indexed in MATLAB; here, adjust accordingly),
          - 'JC': array-like of indices.
    time : array-like
        Array or list of survey times (length = number of passes).

    Returns
    -------
    Flag : int
        Returns 1 if observations were found and the file was created, 0 otherwise.
    
    Notes
    -----
    This function is a translation of a MATLAB routine that
    extracts observations from the SST input and writes a NetCDF file.
    
    Assumes that the following functions are available:
        - nc_addhist(fname, history_string)
      - nc_write(fname, variable_name, data)
    
    Replace or implement these functions as needed.
    """
    
    Flag = 0
    xi = grd['lon_rho']
    yi = grd['lat_rho']
    
    npasses = len(time)
    
    # Initialize lists to accumulate observations.
    XT_list = []
    YT_list = []
    ZT_list = []
    ODepth_list = []
    OType_list = []
    OTime_list = []
    OValue_list = []
    OError_list = []
    LatT_list = []
    LonT_list = []
    survey_time_list = []
    Nobs_list = []
    
    ckk = 0
    got_observations = False
    
    # Loop over the survey passes
    for rec in range(npasses):
        # Get the observation field for this survey record.
        # sst[rec, :, :] should be a 2D array.
        Ti = np.squeeze(sst[rec, :, :])
        Ei = np.ones(Ti.shape)
        Mp, Lp = Ti.shape
        
        # Create grid index arrays from the decimated grid stored in grd.
        # MATLAB: [Y, X] = ndgrid(grd.IC-1, grd.JC-1)
        # In Python, first convert IC and JC to NumPy arrays and subtract 1.
        IC = np.array(grd['IC']) - 1
        JC = np.array(grd['JC']) - 1
        # Use 'ij' indexing to mimic MATLAB’s ndgrid.
        Y_grid, X_grid = np.meshgrid(IC, JC, indexing='ij')
        X = X_grid.ravel()
        Y = Y_grid.ravel()
        
        # For this example, assume the full lon/lat fields (xi,yi)
        # are provided as 2D arrays that match the original grid.
        Lon = xi.ravel()
        Lat = yi.ravel()
        
        # Flatten the observation arrays.
        Ti_vec = Ti.ravel()
        Ei_vec = Ei.ravel()
        
        # Remove any points where the observation equals 0 or is NaN.
        bad_inds = np.where((Ti_vec == 0) | np.isnan(Ti_vec))[0]
        if bad_inds.size > 0:
            Ti_vec = np.delete(Ti_vec, bad_inds)
            Ei_vec = np.delete(Ei_vec, bad_inds)
            X = np.delete(X, bad_inds)
            Y = np.delete(Y, bad_inds)
            Lon = np.delete(Lon, bad_inds)
            Lat = np.delete(Lat, bad_inds)
        
        # Only process if at least five valid observations exist.
        if Ti_vec.size >= 5:
            got_observations = True
            ckk += 1
            survey_time_list.append(time[rec])
            Nobs_list.append(Ti_vec.size)
            # Z is a vector of zeros (same length as valid X values)
            Z = np.zeros_like(X)
            XT_list.append(X)
            YT_list.append(Y)
            LonT_list.append(Lon)
            LatT_list.append(Lat)
            ZT_list.append(Z)
            # For observation depth, use zeros.
            ODepth_list.append(np.zeros_like(X))
            # Set observation type equal to 6 (as in MATLAB code).
            OType_list.append(6 * np.ones_like(X, dtype=int))
            # Associate the survey time with every observation.
            OTime_list.append(np.full_like(X, time[rec], dtype=float))
            # Append the observation values and error (variance)
            OValue_list.append(Ti_vec)
            OError_list.append(Ei_vec**2)
    
    # If any observations were found, proceed to write the NetCDF file.
    if got_observations:
        XT = np.concatenate(XT_list) if XT_list else np.array([])
        YT = np.concatenate(YT_list) if YT_list else np.array([])
        LonT = np.concatenate(LonT_list) if LonT_list else np.array([])
        LatT = np.concatenate(LatT_list) if LatT_list else np.array([])
        ZT = np.concatenate(ZT_list) if ZT_list else np.array([])
        OTime = np.concatenate(OTime_list) if OTime_list else np.array([])
        OType = np.concatenate(OType_list) if OType_list else np.array([])
        ODepth = np.concatenate(ODepth_list) if ODepth_list else np.array([])
        OValue = np.concatenate(OValue_list) if OValue_list else np.array([])
        OError = np.concatenate(OError_list) if OError_list else np.array([])
        survey_time = np.array(survey_time_list)
        Nobs = np.array(Nobs_list)
        
        # Create observation variance: a 7x1 array with the 6th element equal to variance of OValue.
        obs_variance = np.ones((7, 1))
        obs_variance[5, 0] = np.var(OValue)
        
        # Create the NetCDF file.
        OBS_out = fname
        is3d = 1  # Not used further in this Python version.
        survey = survey_time.size
        Nstate = 7
        titulo = 'AMSR2 SST'
        addhistory = 'Prepared by Julia Levin (julia@marine.rutgers.edu)'
        
        # The variable 'obs_provenance_definition' must be defined in your context;
        # here we use a placeholder.
        obs_provenance_definition = None  
        
        # Call the function to define the NetCDF file structure.
        define_4dvar_obs_file(OBS_out, survey)
 
        nc_addhist(OBS_out, titulo)
        nc_addhist(OBS_out, addhistory)
        
        # Write the variables to the file.
        write_to_netcdf(OBS_out, 'spherical', 'T')
        write_to_netcdf(OBS_out, 'Nobs', Nobs)
        write_to_netcdf(OBS_out, 'survey_time', survey_time)
        write_to_netcdf(OBS_out, 'obs_variance', obs_variance)
        write_to_netcdf(OBS_out, 'obs_type', OType)
        write_to_netcdf(OBS_out, 'obs_time', OTime)
        write_to_netcdf(OBS_out, 'obs_Xgrid', XT)
        write_to_netcdf(OBS_out, 'obs_Ygrid', YT)
        write_to_netcdf(OBS_out, 'obs_Zgrid', ZT)
        write_to_netcdf(OBS_out, 'obs_depth', ODepth)
        write_to_netcdf(OBS_out, 'depth', ODepth)
        write_to_netcdf(OBS_out, 'obs_error', OError)
        
        # New variables:
        write_to_netcdf(OBS_out, 'obs_lat', LatT)
        write_to_netcdf(OBS_out, 'obs_lon', LonT)
        Provenance = 324
        write_to_netcdf(OBS_out, 'obs_provenance', Provenance * np.ones_like(OValue))
        Flag = 1
    else:
        Flag = 0
    return Flag
