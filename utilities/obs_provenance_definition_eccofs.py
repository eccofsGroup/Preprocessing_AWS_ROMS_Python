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
