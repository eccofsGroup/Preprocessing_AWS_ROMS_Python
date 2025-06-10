import numpy as np
from scipy.ndimage import convolve

def shapiro2d(arr, order=1, mode='mirror'):
    """
    2D Shapiro filter (separable 3-point smoother), repeated `order` times.
    
    Parameters
    ----------
    arr   : 2D array_like
      Input field.
    order : int
      Number of times to apply the 2-step (x then y) filter.
    mode  : str
      Boundary condition passed to convolve (e.g. 'mirror', 'nearest', 'wrap').
    
    Returns
    -------
    smoothed : 2D ndarray
      Filtered field.
    """
    # 1-D 3-point kernel [1/4, 2/4, 1/4]
    k1 = np.array([1, 2, 1], dtype=float) / 4.0

    out = arr.astype(float)
    for _ in range(order):
        # horizontal pass
        out = convolve(out, k1[None, :], mode=mode)
        # vertical pass
        out = convolve(out, k1[:, None], mode=mode)
    return out
