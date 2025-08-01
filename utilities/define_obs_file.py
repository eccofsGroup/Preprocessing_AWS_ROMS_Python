from netCDF4 import Dataset
import numpy as np

def obs_provenance_definition_eccofs():
    """
    Provides extra-variable definitions for ECCOFS provenance information,
    analogous to the MATLAB obs_provenance_definition_eccofs function.
    
    Returns
    -------
    extra_var : list of dict
        Each dict has keys:
          - 'Name'      : variable name (str)
          - 'Nctype'    : numpy dtype for the variable
          - 'Dimension' : list of dimension names
          - 'Attribute' : list of {'Name':..., 'Value':...} dicts
    """
    extra_var = []

    # --- 1) obs_provenance (int32) ---
    attrs1 = [
        ('SST_super_observations',      np.int32(300)),
        ('LEO_IR_SST',                  np.int32(311)),
        ('GOES_geostationary_IR_SST',   np.int32(317)),
        ('AMSR2_microwave_SST',         np.int32(324)),
        ('Altimeter_superobs',          np.int32(400)),
        ('Jason_3_SSH',                 np.int32(403)),
        ('Sentinel_6a_SSH',             np.int32(404)),
        ('Sentinel_3a_SSH',             np.int32(441)),
        ('Sentinel_3b_SSH',             np.int32(442)),
        ('CMEMS_CORA superobs',         np.int32(800)),
        ('CORA_ARGO_Vertical_Profile',  np.int32(801)),
        ('CORA_Thermistor_Chain',       np.int32(804)),
        ('CORA_XCTD',                   np.int32(805)),
        ('CORA_Gliders',                np.int32(806)),
        ('CORA_CTD',                    np.int32(807)),
        ('CORA_Sea_Mammals',            np.int32(808)),
        ('CORA_Moorings',               np.int32(823)),
        ('CORA_thermosalinograph',      np.int32(824)),
        ('CORA_Drifter_Buoys',          np.int32(825)),
    ]
    var1 = {
        'Name':      'obs_provenance',
        'Nctype':    np.int32,
        'Dimension': ['datum'],
        'Attribute': [{'Name': n, 'Value': v} for n, v in attrs1]
    }
    extra_var.append(var1)

    # --- 2) obs_label (float64) ---
    var2 = {
        'Name':      'obs_label',
        'Nctype':    np.float64,
        'Dimension': ['datum'],
        'Attribute': []  # no attributes in the MATLAB if(0) block
    }
    extra_var.append(var2)

    return extra_var

#--------------------------------------------------------------------------

def define_4dvar_obs_file(output_ncfile, num_surveys, additional_vars=None):
    """
    Define the 4D-Var observation file (using xbt data) by creating
    a new NetCDF file with the required dimensions, variables, and attributes.
    
    Parameters
    ----------
    output_ncfile : str
        The name (or full path) of the NetCDF file to be created.
    num_surveys : int
        The length of the 'survey' dimension.
    additional_vars : list of dict, optional
        A list where each element is a dictionary containing the definition
        of an additional variable. Each dictionary should include the keys:
          - 'Name': variable name (str)
          - 'Nctype': NetCDF type (e.g., 'S1' for characters or 'f8', 'i4', etc.)
          - 'Dimension': a tuple of dimension names (e.g., ('datum',))
          - 'Attribute': a list of dictionaries with keys 'Name' and 'Value'
    
    Returns
    -------
    None
    """
    # Create a new netCDF file (NETCDF4 format).
    ds = Dataset(output_ncfile, 'w', format='NETCDF4')
    
    # The equivalent of nc_padheader is not needed; the library handles headers.
    
    # Create necessary dimensions.
    ds.createDimension('datum', None)        # unlimited dimension
    ds.createDimension('survey', num_surveys)
    ds.createDimension('weight', 8)            # as in MATLAB code
    ds.createDimension('state_variable', 7)    # as in MATLAB code
    
    #-------------------------------------------------------------------------
    # Define variables.
    #-------------------------------------------------------------------------
    # Define 'spherical' (a character variable). Here we define it as a string type.
    var = ds.createVariable('spherical', str)
    var.long_name = 'grid type logical switch'
    var.option_T = 'spherical'
    var.option_F = 'Cartesian'
    
    # obs_variance: type double, dimensions: state_variable
    var = ds.createVariable('obs_variance', 'f8', ('state_variable',))
    var.long_name = 'global time and space observation variance'
    var.units = 'squared state variable units'
    
    # Nobs: type int, dimensions: survey
    var = ds.createVariable('Nobs', 'i4', ('survey',))
    var.long_name = 'number of observations with the same survey time'
    
    # survey_time: type double, dimensions: survey
    var = ds.createVariable('survey_time', 'f8', ('survey',))
    var.long_name = 'survey time'
    var.units = 'days since 2011-01-01 00:00:00'
    
    # obs_type: type int, dimensions: datum
    var = ds.createVariable('obs_type', 'i4', ('datum',))
    var.long_name = 'model state variable associated with observation'
    var.units = 'unity'
    var.option_01 = 'free-surface'
    var.option_02 = 'vertically integrated u-momentum component'
    var.option_03 = 'vertically integrated v-momentum component'
    var.option_04 = 'u-momentum component'
    var.option_05 = 'v-momentum component'
    var.option_06 = 'potential temperature'
    var.option_07 = 'salinity'
    
    # obs_time: type double, dimensions: datum
    var = ds.createVariable('obs_time', 'f8', ('datum',))
    var.long_name = 'time of observation'
    var.units = 'days since 2011-01-01 00:00:00'
    
    # obs_lon: type double, dimensions: datum
    var = ds.createVariable('obs_lon', 'f8', ('datum',))
    var.long_name = 'longitude of observation'
    var.units = 'degrees_east'
    var.valid_min = -180
    var.valid_max = 180
    
    # obs_lat: type double, dimensions: datum
    var = ds.createVariable('obs_lat', 'f8', ('datum',))
    var.long_name = 'latitude of observation'
    var.units = 'degrees_north'
    var.valid_min = -90
    var.valid_max = 90
    
    # obs_depth: type double, dimensions: datum
    var = ds.createVariable('obs_depth', 'f8', ('datum',))
    var.long_name = 'depth_of_observation'
    var.units = 'meter'
    var.negative = 'downwards'
    
    # depth: type double, dimensions: datum
    var = ds.createVariable('depth', 'f8', ('datum',))
    var.long_name = 'depth_of_observation'
    var.units = 'meter'
    var.negative = 'downwards'
    
    # obs_error: type double, dimensions: datum
    var = ds.createVariable('obs_error', 'f8', ('datum',))
    var.long_name = 'observation error, assigned weight'
    var.units = 'magnitude in terms of state variable type'
    
    # obs_value: type double, dimensions: datum
    var = ds.createVariable('obs_value', 'f8', ('datum',))
    var.long_name = 'observation value'
    var.units = 'magnitude in terms of state variable type'
    
    # obs_Xgrid: type double, dimensions: datum
    var = ds.createVariable('obs_Xgrid', 'f8', ('datum',))
    var.long_name = 'x-grid observation location'
    var.units = '1'
    var.left = 'INT(obs_Xgrid(datum))'
    var.right = 'INT(obs_Xgrid(datum))+1'
    
    # obs_Ygrid: type double, dimensions: datum
    var = ds.createVariable('obs_Ygrid', 'f8', ('datum',))
    var.long_name = 'y-grid observation location'
    var.units = '1'
    var.top = 'INT(obs_Ygrid(datum))+1'
    var.bottom = 'INT(obs_Ygrid(datum))'
    
    # obs_Zgrid: type double, dimensions: datum
    var = ds.createVariable('obs_Zgrid', 'f8', ('datum',))
    var.long_name = 'z-grid observation location'
    var.units = '1'
    var.up = 'INT(obs_Ygrid(datum))+1'
    var.down = 'INT(obs_Ygrid(datum))'
    
    #-------------------------------------------------------------------------
    # Define any additional variables if provided.
    if additional_vars is not None:
        for var_dict in additional_vars:
            name = var_dict.get('Name')
            # Expect 'Nctype' to be a netCDF type string, e.g., 'f8' or 'i4' or 'S1'.
            nctype = var_dict.get('Nctype')
            dims = tuple(var_dict.get('Dimension', ()))
            var = ds.createVariable(name, nctype, dims)
            # If there are attributes, set them.
            if 'Attribute' in var_dict:
                for att in var_dict['Attribute']:
                    attname = att['Name']
                    attvalue = att['Value']
                    var.setncattr(attname, attvalue)
    
    #-------------------------------------------------------------------------
    # Set global attributes.
    ds.setncattr('type', 'Roms observation file')
    ds.setncattr('state_units', 'Free-surface (m), 2D momentum (m/s), 3D momentum (m/s), potential temperature (Celsius), salinity (PSU)')
    
    ds.close()

    pass
