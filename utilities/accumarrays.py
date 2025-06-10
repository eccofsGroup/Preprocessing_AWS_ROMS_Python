import pandas as pd
import numpy as np

# Helper function for 3D binning via pandas
def accum3d(index1, index2, index3, values, shape, func):
    df = pd.DataFrame({
        'd': index1,
        'v': index2,
        't': index3,
        'val': values
    })

    grouped = df.groupby(['d', 'v', 't'])['val'].agg(func).unstack(fill_value=np.nan)
    full = np.full(shape, np.nan)
    for (d, v, t), value in grouped.stack().items():
        full[int(d), int(v), int(t)] = value
    return full

# Helper function for 2D binning via pandas
def accum2d(index1, index2, values, shape, func):
    df = pd.DataFrame({
        'd': index1,
        'v': index2,
        'val': values
    })
    grouped = df.groupby(['d', 'v'])['val'].agg(func).unstack(fill_value=np.nan)
    full = np.full(shape, np.nan)
    for (d, v), value in grouped.stack().items():
        full[int(d), int(v)] = value
    return full

