import numpy as np
import netCDF4 as nc
import warnings


def nc_isvar(ncfile, varname):
    """Return True if variable exists in the netCDF dataset."""
    return varname in ncfile.variables

#--------------------------------------------------------------------------------------

def set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N, igrid, h, zeta=None, report=1):
    """
    Compute ROMS grid depth from vertical stretched variables.
    
    Parameters
    ----------
    Vtransform : int
        Vertical transformation equation:
          1 = original ROMS transformation,
          2 = ROMS-UCLA transformation.
    Vstretching : int
        Vertical stretching function:
          1 = original (Song and Haidvogel, 1994)
          2 = Shchepetkin (2005)
          3 = Geyer BBL refinement
          4 = Shchepetkin (2010)
    theta_s : float
        S-coordinate surface control parameter.
    theta_b : float
        S-coordinate bottom control parameter.
    hc : float
        Width (m) of the surface or bottom boundary layer where higher vertical
        resolution is required.
    N : int
        Number of vertical levels.
    igrid : int
        C-grid type:
          1  => density (RHO) points,
          2  => streamfunction (PSI) points,
          3  => U-velocity points,
          4  => V-velocity points,
          5  => W-velocity points.
    h : 2D numpy array
        Bottom depth at RHO-points (in m, positive); shape (Lp, Mp).
    zeta : 2D numpy array, optional
        Free-surface elevation at RHO-points; if not provided, assumed zero.
    report : int or bool, optional
        If True (or 1), prints detailed information (default is 1).
    
    Returns
    -------
    z : 3D numpy array
        Depths (in m, negative) at the requested C-grid location.
    
    Notes
    -----
    For Vtransform == 1, the transformation is:
    
        z0 = ( s(k) - C(k) ) * hc + C(k)*H
        z = z0 + zeta * ( 1 + z0 / H )
    
    For Vtransform == 2, the transformation is:
    
        z0 = [hc*s(k) + C(k)*H] / (hc+H)
        z = zeta + (zeta+H)*z0
    
    The vertical stretching arrays (s, C) are computed using the provided
    stretching function. The free-surface zeta is averaged to match the grid
    type when igrid = 2,3,4.
    """
    # Check Vtransform and Vstretching values.
    if Vtransform < 1 or Vtransform > 2:
        raise ValueError(f"Illegal parameter Vtransform = {Vtransform}")
    if Vstretching < 1 or Vstretching > 4:
        raise ValueError(f"Illegal parameter Vstretching = {Vstretching}")
    if Vtransform == 1 and hc > np.min(h):
        raise ValueError(f"Critical depth exceeds minimum bathymetry value: "
                         f"Vtransform = {Vtransform}, hc = {hc}, hmin = {np.min(h)}")
    
    # If zeta is not provided, assume zero free-surface.
    if zeta is None:
        zeta = np.zeros_like(h)
    # Default reporting.
    if report is None:
        report = 1

    Np = N + 1
    Lp, Mp = h.shape
    L = Lp - 1
    M = Mp - 1

    # Optionally report Vtransform and grid type.
    if report:
        print("")
        if Vtransform == 1:
            print(f"Vtransform  = {Vtransform}   original ROMS")
        elif Vtransform == 2:
            print(f"Vtransform  = {Vtransform}   ROMS-UCLA")
        
        if igrid == 1:
            print(f"   igrid    = {igrid}   at horizontal RHO-points")
        elif igrid == 2:
            print(f"   igrid    = {igrid}   at horizontal PSI-points")
        elif igrid == 3:
            print(f"   igrid    = {igrid}   at horizontal U-points")
        elif igrid == 4:
            print(f"   igrid    = {igrid}   at horizontal V-points")
        elif igrid == 5:
            print(f"   igrid    = {igrid}   at horizontal RHO-points")
    
    # Determine whether to compute s-curve at W-points or RHO-points.
    # If igrid == 5, use kgrid = 1, else kgrid = 0.
    if igrid == 5:
        kgrid = 1
    else:
        kgrid = 0

    # Call the stretching function to compute the s and C arrays.
    # (Assumes the stretching function is defined in your environment.)
    s, C = stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, bool(report))
    
    #--------------------------------------------------------------------------
    # Average bathymetry and free-surface to match the requested grid.
    #--------------------------------------------------------------------------
    if igrid in [1, 5]:
        hr = h.copy()
        zetar = zeta.copy()
    elif igrid == 2:
        hp = 0.25 * (h[0:L, 0:M] + h[1:Lp, 0:M] + h[0:L, 1:Mp] + h[1:Lp, 1:Mp])
        zetap = 0.25 * (zeta[0:L, 0:M] + zeta[1:Lp, 0:M] + zeta[0:L, 1:Mp] + zeta[1:Lp, 1:Mp])
    elif igrid == 3:
        hu = 0.5 * (h[0:L, :] + h[1:Lp, :])
        zetau = 0.5 * (zeta[0:L, :] + zeta[1:Lp, :])
    elif igrid == 4:
        hv = 0.5 * (h[:, 0:M] + h[:, 1:Mp])
        zetav = 0.5 * (zeta[:, 0:M] + zeta[:, 1:Mp])
    
    #--------------------------------------------------------------------------
    # Compute depths z at the requested grid location.
    #--------------------------------------------------------------------------
    if Vtransform == 1:
        if igrid == 1:
            z = np.empty((Lp, Mp, N))
            for k in range(N):
                z0 = (s[k] - C[k]) * hc + C[k] * hr
                z[..., k] = z0 + zetar * (1.0 + z0 / hr)
        elif igrid == 2:
            z = np.empty((L, M, N))
            for k in range(N):
                z0 = (s[k] - C[k]) * hc + C[k] * hp
                z[..., k] = z0 + zetap * (1.0 + z0 / hp)
        elif igrid == 3:
            z = np.empty((L, Mp, N))
            for k in range(N):
                z0 = (s[k] - C[k]) * hc + C[k] * hu
                z[..., k] = z0 + zetau * (1.0 + z0 / hu)
        elif igrid == 4:
            z = np.empty((Lp, M, N))
            for k in range(N):
                z0 = (s[k] - C[k]) * hc + C[k] * hv
                z[..., k] = z0 + zetav * (1.0 + z0 / hv)
        elif igrid == 5:
            # For igrid==5, use all N+1 levels.
            z = np.empty((Lp, Mp, Np))
            z[..., 0] = -hr
            for k in range(1, Np):
                z0 = (s[k] - C[k]) * hc + C[k] * hr
                z[..., k] = z0 + zetar * (1.0 + z0 / hr)
    elif Vtransform == 2:
        if igrid == 1:
            z = np.empty((Lp, Mp, N))
            for k in range(N):
                z0 = (hc * s[k] + C[k] * hr) / (hc + hr)
                z[..., k] = zetar + (zetar + hr) * z0
        elif igrid == 2:
            z = np.empty((L, M, N))
            for k in range(N):
                z0 = (hc * s[k] + C[k] * hp) / (hc + hp)
                z[..., k] = zetap + (zetap + hp) * z0
        elif igrid == 3:
            z = np.empty((L, Mp, N))
            for k in range(N):
                z0 = (hc * s[k] + C[k] * hu) / (hc + hu)
                z[..., k] = zetau + (zetau + hu) * z0
        elif igrid == 4:
            z = np.empty((Lp, M, N))
            for k in range(N):
                z0 = (hc * s[k] + C[k] * hv) / (hc + hv)
                z[..., k] = zetav + (zetav + hv) * z0
        elif igrid == 5:
            z = np.empty((Lp, Mp, Np))
            for k in range(Np):
                z0 = (hc * s[k] + C[k] * hr) / (hc + hr)
                z[..., k] = zetar + (zetar + hr) * z0

    return z

    #-------------------------------------------------------------------------------

def stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report=False):
    """
    Compute the ROMS vertical coordinate stretching function.
    
    [s, C] = stretching(Vstretching, theta_s, theta_b, hc, N, kgrid, report)
    
    Parameters
    ----------
    Vstretching : int
        Vertical stretching function:
          - 1: Original (Song and Haidvogel, 1994)
          - 2: Shchepetkin (UCLA-ROMS, 2005)
          - 3: R. Geyer BBL refinement
          - 4: Shchepetkin (UCLA-ROMS, 2010)
    theta_s : float
        S-coordinate surface control parameter.
    theta_b : float
        S-coordinate bottom control parameter.
    hc : float
        Width (m) of the surface or bottom boundary layer where higher vertical
        resolution is required.
    N : int
        Number of vertical levels.
    kgrid : int
        Depth grid type logical switch:
          - 0: function at vertical RHO-points
          - 1: function at vertical W-points
    report : bool, optional
        If True, prints detailed information (default is False).
    
    Returns
    -------
    s : ndarray
        S-coordinate independent variable, a 1D vector in the range [-1, 0]
        at vertical RHO- or W-points.
    C : ndarray
        Non-dimensional, monotonic vertical stretching function C(s), a 1D
        array in the range [-1, 0].
    """
    
    # Check Vstretching
    if Vstretching < 1 or Vstretching > 4:
        raise ValueError(f"Illegal parameter Vstretching = {Vstretching}")
    
    ds = 1.0 / N
    
    # Set vertical grid depending on kgrid:
    if kgrid == 1:
        # W-points: Nlev = N + 1, levels 0, 1, …, N
        Nlev = N + 1
        lev = np.arange(0, N + 1)
        s = (lev - N) * ds
    else:
        # RHO-points: Nlev = N, levels shifted by 0.5: 0.5, 1.5, ..., N-0.5
        Nlev = N
        lev = np.arange(1, N + 1) - 0.5
        s = (lev - N) * ds
    
    #--------------------------------------------------------------------------
    # Compute ROMS S-coordinates vertical stretching function
    #--------------------------------------------------------------------------
    if Vstretching == 1:
        # Original vertical stretching function (Song and Haidvogel, 1994)
        if theta_s > 0:
            Ptheta = np.sinh(theta_s * s) / np.sinh(theta_s)
            Rtheta = np.tanh(theta_s * (s + 0.5)) / (2.0 * np.tanh(0.5 * theta_s)) - 0.5
            C = (1.0 - theta_b) * Ptheta + theta_b * Rtheta
        else:
            C = s.copy()
    
    elif Vstretching == 2:
        # A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function.
        alfa = 1.0
        beta = 1.0
        ds = 1.0 / N  # (ensuring consistency)
        if kgrid == 1:
            Nlev = N + 1
            lev = np.arange(0, N + 1)
            s = (lev - N) * ds
        else:
            Nlev = N
            lev = np.arange(1, N + 1) - 0.5
            s = (lev - N) * ds
        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
            if theta_b > 0:
                Cbot = -1.0 + np.sinh(theta_b * (s + 1.0)) / np.sinh(theta_b)
                weight = (s + 1.0) ** alfa * (1.0 + (alfa / beta) * (1.0 - (s + 1.0) ** beta))
                C = weight * Csur + (1.0 - weight) * Cbot
            else:
                C = Csur
        else:
            C = s.copy()
    
    elif Vstretching == 3:
        # R. Geyer BBL vertical stretching function.
        ds = 1.0 / N
        if kgrid == 1:
            Nlev = N + 1
            lev = np.arange(0, N + 1)
            s = (lev - N) * ds
        else:
            Nlev = N
            lev = np.arange(1, N + 1) - 0.5
            s = (lev - N) * ds
        if theta_s > 0:
            exp_s = theta_s      # surface stretching exponent
            exp_b = theta_b      # bottom stretching exponent
            alpha_val = 3.0      # scale factor for hyperbolic functions
            Cbot = np.log(np.cosh(alpha_val * (s + 1) ** exp_b)) / np.log(np.cosh(alpha_val)) - 1.0
            Csur = -np.log(np.cosh(alpha_val * np.abs(s) ** exp_s)) / np.log(np.cosh(alpha_val))
            weight = (1.0 - np.tanh(alpha_val * (s + 0.5))) / 2.0
            C = weight * Cbot + (1.0 - weight) * Csur
        else:
            C = s.copy()
    
    elif Vstretching == 4:
        # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
        ds = 1.0 / N
        if kgrid == 1:
            Nlev = N + 1
            lev = np.arange(0, N + 1)
            s = (lev - N) * ds
        else:
            Nlev = N
            lev = np.arange(1, N + 1) - 0.5
            s = (lev - N) * ds
        if theta_s > 0:
            Csur = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)
        else:
            Csur = - s ** 2
        if theta_b > 0:
            Cbot = (np.exp(theta_b * Csur) - 1.0) / (1.0 - np.exp(-theta_b))
            C = Cbot
        else:
            C = Csur
    
    #--------------------------------------------------------------------------
    # Report S-coordinate parameters.
    #--------------------------------------------------------------------------
    if report:
        print("\n")
        if Vstretching == 1:
            print(f"Vstretching = {Vstretching}   Song and Haidvogel (1994)")
        elif Vstretching == 2:
            print(f"Vstretching = {Vstretching}   Shchepetkin (2005)")
        elif Vstretching == 3:
            print(f"Vstretching = {Vstretching}   Geyer (2009), BBL")
        elif Vstretching == 4:
            print(f"Vstretching = {Vstretching}   Shchepetkin (2010)")
        if kgrid == 1:
            print(f"   kgrid    = {kgrid}   at vertical W-points")
        else:
            print(f"   kgrid    = {kgrid}   at vertical RHO-points")
        print(f"   theta_s  = {theta_s}")
        print(f"   theta_b  = {theta_b}")
        print(f"   hc       = {hc}")
        print("\n S-coordinate curves: k, s(k), C(k)\n")
        if kgrid == 1:
            for k in range(Nlev - 1, -1, -1):
                print(f"    {k:3d}   {s[k]:20.12e}   {C[k]:20.12e}")
        else:
            for k in range(Nlev - 1, -1, -1):
                print(f"    {k+1:3d}   {s[k]:20.12e}   {C[k]:20.12e}")
        print("\n")
    
    return s, C

    #-------------------------------------------------------------------------------------

def roms_get_grid(grd_file, scoord=None, tindex=0, calc_zuv=False):
    """
    Gets the longitude, latitude, mask, depth [and z coordinates] from a ROMS NetCDF grid file.
    
    Parameters:
      grd_file : str or dict
          If a string, the file name of the ROMS grid file. If already a grid dictionary,
          vertical coordinates will be added/updated.
      scoord : list/tuple or str, optional
          S-coordinate parameters, either as:
            • a 4-element vector [theta_s, theta_b, Tcline, N] (assumes Vtransform=1, Vstretching=1),
            • a 6-element vector [theta_s, theta_b, Tcline, N, Vtransform, Vstretching], or
            • a file name (or URL) from which to read these parameters.
      tindex : int or 2D array, optional
          How to obtain zeta:
            • 0: assume zeta = 0,
            • integer: use as a time index to read zeta from the ROMS file, or
            • a 2D array: explicitly provided zeta values.
      calc_zuv : bool, optional
          If True, compute z_u and z_v on the ROMS C-grid.
          
    Returns:
      grd : dict
          Dictionary containing grid and (if computed) vertical coordinate information.
    """
    
    # If grd_file is already a grid structure (dictionary), use it.
    if isinstance(grd_file, dict):
        grd = grd_file
    else:
        # Open the grid file for reading.
        ds = nc.Dataset(grd_file, 'r')
        grd = {'grd_file': grd_file}
        
        # Read variables needed for the grid.
        varlist = ['mask_rho', 'mask_psi', 'mask_u', 'mask_v', 'h',
                   'pm', 'pn', 'f', 'angle', 'dmde', 'dndx', 'visc_factor', 'diff_factor']
        for vname in varlist:
            try:
                grd[vname] = ds.variables[vname][:]
            except KeyError:
                warnings.warn(f"Variable not found: {vname}", UserWarning)
                if vname == 'angle' and 'h' in grd:
                    grd[vname] = np.zeros_like(grd['h'])
                if vname == 'h':
                    warnings.warn("Using initial bath for h", UserWarning)
                    try:
                        grd[vname] = np.squeeze(ds.variables['bath'][:])
                    except KeyError:
                        grd[vname] = None
        
        # Get coordinate variables if available.
        varlist2 = ['x_rho', 'y_rho', 'x_u', 'y_u', 'x_v', 'y_v', 'x_psi', 'y_psi']
        if nc_isvar(ds, 'x_rho'):
            for vname in varlist2:
                try:
                    grd[vname] = ds.variables[vname][:]
                except KeyError:
                    warnings.warn(f"Variable not found: {vname}", UserWarning)
        
        varlist3 = ['lon_rho', 'lat_rho', 'lon_psi', 'lat_psi',
                    'lon_v', 'lat_v', 'lon_u', 'lat_u']
        for vname in varlist3:
            try:
                grd[vname] = ds.variables[vname][:]
            except KeyError:
                # Substitute x/y coordinates if lon/lat not available.
                if vname.startswith('lon'):
                    usevname = vname.replace('lon', 'x')
                else:
                    usevname = vname.replace('lat', 'y')
                try:
                    grd[vname] = grd[usevname]
                except KeyError:
                    pass
                grd['nolatlon'] = 1
        
        # Compute a bounding box if possible.
        try:
            grd['bounding_box'] = [np.nanmin(grd['lon_psi']), np.nanmax(grd['lon_psi']),
                                   np.nanmin(grd['lat_psi']), np.nanmax(grd['lat_psi'])]
        except Exception:
            pass
        
        # Construct perimeter if possible.
        try:
            lon_psi = grd['lon_psi']
            lat_psi = grd['lat_psi']
            perim_lon = np.concatenate([lon_psi[0, :],
                                        lon_psi[:, -1],
                                        lon_psi[-1, ::-1],
                                        lon_psi[::-1, 0]])
            perim_lat = np.concatenate([lat_psi[0, :],
                                        lat_psi[:, -1],
                                        lat_psi[-1, ::-1],
                                        lat_psi[::-1, 0]])
            grd['perimeter'] = np.column_stack((perim_lon, perim_lat))
        except Exception:
            pass
        
        # Read additional optional attributes.
        for vname in ['rdrag', 'rdrag2', 'ZoBot']:
            if nc_isvar(ds, vname):
                grd[vname] = ds.variables[vname][:]
        
        # Set up the mask variables.
        if 'mask_rho' in grd:
            grd['mask_rho_nan'] = grd['mask_rho'].copy()
            grd['mask_rho_nan'][grd['mask_rho_nan'] == 0] = np.nan
        else:
            if 'h' in grd:
                grd['mask_rho'] = np.ones_like(grd['h'])
                grd['mask_rho_nan'] = grd['mask_rho']
                grd['mask_u'] = np.ones_like(grd['h'][:, 1:])
                grd['mask_v'] = np.ones_like(grd['h'][1:, :])
                grd['mask_psi'] = np.ones_like(grd['h'][1:, 1:])
                grd['nomask'] = 1
        
        # Read coastline data if available.
        for var in ['lon_coast', 'lat_coast']:
            try:
                grd[var] = ds.variables[var][:]
            except KeyError:
                pass
        
        # If tindex is provided (and nonzero as a single integer), attempt to read wet/dry masks.
        if tindex != 0 and isinstance(tindex, int):
            if nc_isvar(ds, 'wetdry_mask_rho'):
                try:
                    usewetdry = True  # Here you might use a user preference.
                    for vname in ['wetdry_mask_rho', 'wetdry_mask_psi',
                                  'wetdry_mask_u', 'wetdry_mask_v']:
                        grd[vname] = ds.variables[vname][tindex - 1, ...]
                except KeyError:
                    pass
        
        # If the grid file is a refinement grid, get the associated attributes.
        try:
            grd['refine_factor'] = ds.getncattr('refine_factor')
            for attr in ['parent_Imin', 'parent_Imax', 'parent_Jmin', 'parent_Jmax']:
                grd[attr] = ds.getncattr(attr)
        except Exception:
            pass
        
        ds.close()
    
    # If s-coordinate parameters are provided, calculate vertical coordinates.
    if scoord is not None:
        # Get the bathymetry.
        h = grd.get('h')
        if h is None:
            raise ValueError("Grid 'h' (bathymetry) is required for vertical coordinate computation.")
        Mp, Lp = h.shape
        L = Lp - 1
        M = Mp - 1
        
        # Determine s-coordinate parameters.
        if not isinstance(scoord, str):
            # Assume scoord is vector-like.
            theta_s = scoord[0]
            theta_b = scoord[1]
            Tcline  = scoord[2]
            N       = scoord[3]
            if len(scoord) < 5:
                Vtransform = 1
                Vstretching = 1
            else:
                Vtransform = scoord[4]
                Vstretching = scoord[5]
            hc = Tcline
        else:
            # Read s-coordinate parameters from the provided file.
            ds_s = nc.Dataset(scoord, 'r')
            theta_s = ds_s.variables['theta_s'][:].item()
            theta_b = ds_s.variables['theta_b'][:].item()
            hc      = ds_s.variables['hc'][:].item()
            try:
                Tcline = ds_s.variables['Tcline'][:].item()
            except KeyError:
                Tcline = hc
            Cs_r = ds_s.variables['Cs_r'][:]
            N = Cs_r.shape[0]
            Vtransform = ds_s.variables['Vtransform'][:] if nc_isvar(ds_s, 'Vtransform') else 1
            Vstretching = ds_s.variables['Vstretching'][:] if nc_isvar(ds_s, 'Vstretching') else 1
            ds_s.close()
        
        if Vtransform == 1:
            hc = min(Tcline, np.nanmin(h))
        
        # Compute s and C(s) curves for rho and w points.
        s_rho, Cs_r = stretching(Vstretching, theta_s, theta_b, hc, N, 0, 0)
        s_w, Cs_w   = stretching(Vstretching, theta_s, theta_b, hc, N, 1, 0)
        sc_r = s_rho
        sc_w = s_w
        
        # Determine zeta (defaulting to zeros).
        zeta = np.zeros_like(h)
        if tindex != 0:
            if np.isscalar(tindex):
                if not isinstance(scoord, str):
                    warnings.warn("Cannot process zeta from file when scoord parameters are given as a vector.")
                    try:
                        ds_temp = nc.Dataset(grd['grd_file'], 'r')
                        zeta = ds_temp.variables['zeta'][tindex - 1, ...]
                        ds_temp.close()
                    except KeyError:
                        warnings.warn("zeta not found. Assuming zero.")
                        zeta = np.zeros_like(h)
                else:
                    ds_temp = nc.Dataset(scoord, 'r')
                    try:
                        zeta = ds_temp.variables['zeta'][tindex - 1, ...]
                    except KeyError:
                        warnings.warn("zeta not found in scoord file. Assuming zero.")
                        zeta = np.zeros_like(h)
                    ds_temp.close()
            else:
                # tindex is treated as a 2D array of zeta values.
                if tindex.shape != h.shape:
                    raise ValueError("Input zeta dimensions do not match grid 'h' dimensions.")
                else:
                    zeta = tindex
        grd['zeta'] = zeta
        
        # Compute vertical depths at rho-points.
        rgrid = 1
        # Note: h and zeta are transposed to mimic MATLAB's use of h' and zeta'
        z_r = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N,
                        rgrid, h.T, zeta.T, 0)
        # Permute dimensions: MATLAB's permute(z_r, [3,2,1]) is equivalent to:
        grd['z_r'] = np.transpose(z_r, (2, 1, 0))
        
        # Compute depths at w-points.
        wgrid = 5
        z_w = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N,
                        wgrid, h.T, zeta.T, 0)
        grd['z_w'] = np.transpose(z_w, (2, 1, 0))
        
        # Compute cell thicknesses at rho-points.
        grd['dz'] = np.diff(grd['z_w'], axis=0)
        
        # Compute cell areas and volumes if possible.
        try:
            dA = 1.0 / (grd['pm'] * grd['pn'])
            # Broadcast dA to the shape of dz.
            dA3D = np.broadcast_to(dA, grd['dz'].shape)
            dV = dA3D * grd['dz']
            grd['dA'] = dA
            grd['dV'] = dV
        except KeyError:
            pass
        
        if calc_zuv:
            # Compute z coordinates on u-points.
            ugrid = 3
            z_u = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N,
                            ugrid, h.T, zeta.T, 0)
            grd['z_u'] = np.transpose(z_u, (2, 1, 0))
            
            # Compute u-point cell edges as the average of w-points.
            z_w_perm = np.transpose(z_w, (2, 1, 0))
            z_uw = 0.5 * (z_w_perm[:-1, :, :] + z_w_perm[1:, :, :])
            grd['z_uw'] = z_uw
            
            # Compute z coordinates on v-points.
            vgrid = 4
            z_v = set_depth(Vtransform, Vstretching, theta_s, theta_b, hc, N,
                            vgrid, h.T, zeta.T, 0)
            grd['z_v'] = np.transpose(z_v, (2, 1, 0))
            
            # Compute v-point cell edges.
            z_v_perm = np.transpose(z_v, (2, 1, 0))
            z_vw = 0.5 * (z_v_perm[:, :-1, :] + z_v_perm[:, 1:, :])
            grd['z_vw'] = z_vw
        
        # Save s-coordinate parameters in the grid structure.
        grd['Vtransform'] = Vtransform
        grd['Vstretching'] = Vstretching
        grd['theta_s'] = theta_s
        grd['theta_b'] = theta_b
        grd['Tcline'] = Tcline
        grd['N'] = N
        grd['hc'] = hc
        grd['sc_w'] = sc_w
        grd['Cs_w'] = Cs_w
        grd['sc_r'] = sc_r
        grd['Cs_r'] = Cs_r
        grd['s_w'] = s_w
        grd['s_rho'] = s_rho

    return grd
