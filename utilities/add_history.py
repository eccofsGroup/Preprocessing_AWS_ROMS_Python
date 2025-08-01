from netCDF4 import Dataset
from datetime import datetime


def add_history(nc_file, entry, prepend=True):
    """
    Add a timestamped entry to the global 'history' attribute of a netCDF file.

    Parameters
    ----------
    nc_file : str
        Path to the netCDF file.
    entry : str
        Text to add.
    prepend : bool, optional
        If True, new entries go on top; if False, at the bottom.
    """
    # Open for read/write
    with Dataset(nc_file, 'r+') as ds:
        # Get existing history (empty string if none)
        old_hist = getattr(ds, 'history', '')
        
        # Create timestamped entry (UTC ISO format)
        stamp = datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')
        new_line = f"{stamp} {entry}"
        
        # Combine
        if prepend:
            combined = new_line + '\n' + old_hist
        else:
            combined = old_hist + ('\n' if old_hist else '') + new_line
        
        # Write back
        ds.history = combined
