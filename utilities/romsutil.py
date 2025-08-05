import os,glob
import numpy as np
import xroms
import xarray as xr
import pandas as pd
import netCDF4 as nc
import datetime as dt
from math import radians, cos, sin, asin, sqrt,atan2


    
def stretching(Vstretching, theta_s, theta_b, hc, N, kgrid):
    
#--------------------------------------------------------------------------
# Compute ROMS S-coordinates vertical stretching function
#--------------------------------------------------------------------------
    Np=N+1;

# Original vertical stretching function (Song and Haidvogel, 1994).

    if Vstretching==1:
        ds=1.0/N;
        if kgrid==1:


            lev=np.arange(0,Np)
            s=(lev-N)*ds
        else:
            lev=np.arange(1,Np)-0.5
            s=(lev-N)*ds
        
        if theta_s > 0:
            Ptheta=np.sinh(theta_s*s)/np.sinh(theta_s)
            Rtheta=np.tanh(theta_s*(s+0.5))/(2.0*np.tanh(0.5*theta_s))-0.5
            C=(1.0-theta_b)*Ptheta+theta_b*Rtheta
        else:
            C=s

# A. Shchepetkin (UCLA-ROMS, 2005) vertical stretching function.

    elif (Vstretching == 2):

        alfa=1.0
        beta=1.0
        ds=1.0/N
        if kgrid == 1:

            lev=np.arange(0,Np)
            s=(lev-N)*ds
        else:

            lev=np.arange(1,Np)-0.5
            s=(lev-N)*ds
  
        if theta_s > 0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if theta_b > 0:
                Cbot=-1.0+np.sinh(theta_b*(s+1.0))/np.sinh(theta_b)
                weigth=(s+1.0)**alfa*(1.0+(alfa/beta)*(1.0-(s+1.0)**beta))
                C=weigth*Csur+(1.0-weigth)*Cbot
            else:
                C=Csur
    
        else:
            C=s
            
            
            #R. Geyer BBL vertical stretching function.

    elif (Vstretching == 3):

        ds=1.0/N
        if kgrid == 1:
            lev=np.arange(0,Np)
            s=(lev-N)*ds
        else:

            lev=np.arange(1,Np)-0.5
            s=(lev-N)*ds
  
        if theta_s > 0:
            exp_s=theta_s      #  surface stretching exponent
            exp_b=theta_b      #  bottom  stretching exponent
            alpha=3            #  scale factor for all hyperbolic functions
            Cbot=np.log(np.cosh(alpha*(s+1)**exp_b))/np.log(np.cosh(alpha))-1
            Csur=-np.log(np.cosh(alpha*np.abs(s)**exp_s))/np.log(np.cosh(alpha))
            weight=(1-np.tanh( alpha*(s+.5)))/2;
            C=weight*Cbot+(1-weight)*Csur;
        else:
            C=s




 # A. Shchepetkin (UCLA-ROMS, 2010) double vertical stretching function
 # with bottom refinement

    elif Vstretching == 4:

        ds=1.0/N;
        if (kgrid == 1):
            lev=np.arange(0,Np)
            s=(lev-N)*ds;
        else:

            lev=np.arange(1,Np)-0.5
            s=(lev-N)*ds;
  
        if (theta_s > 0):
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0);
        else:
            Csur=-s^2;
  
        if (theta_b > 0):
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
  

 # Quadratic formulation to enhance surface exchange.

 # (J. Souza, B.S. Powell, A.C. Castillo-Trujillo, and P. Flament, 2014:
 #  The Vorticity Balance of the Ocean Surface in Hawaii from a
 #  Regional Reanalysis.'' J. Phys. Oceanogr., 45, 424-440)

    elif Vstretching == 5:

        if (kgrid == 1):
            lev=np.arange(0,Np)
            s=-(lev*lev - 2.0*lev*N + lev + N*N - N) / (N*N - N)- 0.01*(lev*lev - lev*N) / (1.0 - N)
            s[0]=-1.0;
        else:
            lev=np.arange(1,Np)-0.5
            s=-(lev*lev - 2.0*lev*N + lev + N*N - N) / (N*N - N)-0.01*(lev*lev - lev*N) / (1.0 - N)
  
        if theta_s > 0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=-s^2;
  
        if theta_b > 0:
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
  




    return (s,C)

def create_roms_rivers_file(riverfile,data):

    fh        = nc.Dataset(riverfile, "w", format="NETCDF4",clobber=True)
    fh.type   = data['type'] 
    fh.title  =  data['title'] 
    fh.grd_file  =  data['gridfile'] 
    fh.history = data['history']
    fh.source=data['source']
    fh.citation=data['citation']
    fh.rivers=data['rivers']
#     dims=grd.sizes
#     atts=grd.attrs
    fh.createDimension('river', data['nriver'])
    fh.createDimension('river_time', None)
    fh.createDimension('s_rho', data['nz'])

        
# # Defining hydrodynamic variables of initialization file.  

    fh.createVariable('river', 'double', ('river'))
    fh.variables['river'].long_name = 'river runoff identification number'

    fh.createVariable('river_Xposition', 'double', ('river'))
    fh.variables['river_Xposition'].long_name = "river XI-position at RHO-points" 
    fh.variables['river_Xposition'].valid_min = 1. 
    fh.variables['river_Xposition'].valid_max = 1665. 
    fh.variables['river_Xposition'].LuvSrc_meaning = "i-index grid cell of U- or V-face source/sink" 
    fh.variables['river_Xposition'].LwSrc_meaning = "i-index grid cell of RHO-centered source/sink" 
                
    fh.createVariable('river_Eposition', 'double', ('river'))
    fh.variables['river_Eposition'].long_name = "river ETA-position at RHO-points" 
    fh.variables['river_Eposition'].valid_min = 1. 
    fh.variables['river_Eposition'].valid_max = 1441. 
    fh.variables['river_Eposition'].LuvSrc_meaning = "j-index grid cell of U- or V-face source/sink" 
    fh.variables['river_Eposition'].LwSrc_meaning = "j-index grid cell of RHO-centered source/sink" 
                
    fh.createVariable('river_direction', 'double', ('river'))
    fh.variables['river_direction'].long_name = "river runoff grid-cell face flag" 
    fh.variables['river_direction'].flag_values = "0, 1, 2" ;
    fh.variables['river_direction'].flag_meanings = "LuvSrc flow across U-face, LuvSrc flow across V-face, LwSrc flow into cell center" 
    
    fh.createVariable('river_lat', 'double', ('river'))
    fh.variables['river_lat'].long_name = "latitude of grid cell in GloFAS" 
    
    fh.createVariable('river_lon', 'double', ('river'))
    fh.variables['river_lon'].long_name = "longitude of grid cell in GloFAS" 
    
    fh.createVariable('river_Vshape', 'double', ('s_rho','river'))
    fh.variables['river_Vshape'].long_name = "river runoff mass transport vertical profile" ;
    fh.variables['river_Vshape'].requires = "must sum to 1 over s_rho" 

    fh.createVariable('river_salt', 'double', ('river_time', 's_rho','river'))
    fh.variables['river_salt'].long_name = "river runoff salinity" 
    fh.variables['river_salt'].time = "river_time" 
    
    fh.createVariable('river_sign', 'double', ('river'))
    fh.variables['river_sign'].long_name = "river runoff sense of flow" 
    fh.variables['river_sign'].flag_values = "-1, 1" 
    fh.variables['river_sign'].LuvSrc_meaning = "flow in negative u,v direction, flow in positive u,v, direction" ;
    fh.variables['river_sign'].LwSrc_meaning = "flag not used" 

    fh.createVariable('river_temp', 'double', ('river_time', 's_rho','river'))
    fh.variables['river_temp'].long_name = "river runoff potential temperature" 
    fh.variables['river_temp'].units = "Celsius" 
    fh.variables['river_temp'].time = "river_time" 
    
    fh.createVariable('river_time', 'double', ('river_time'))
    fh.variables['river_time'].long_name = "river runoff time" 
    fh.variables['river_time'].units = data['tunits']

    fh.createVariable('river_transport', 'double', ('river_time','river'))
    fh.variables['river_transport'].long_name = "river runoff vertically integrated mass transport" 
    fh.variables['river_transport'].units = "meter3 second-1" 
    fh.variables['river_transport'].positive = "LuvSrc=T flow in positive u,v direction, LwSrc=T flow into RHO-cell" 
    fh.variables['river_transport'].negative = "LuvSrc=T flow in negative u,v direction, LwSrc=T flow out of RHO-cell" 
    fh.variables['river_transport'].time = "river_time" 
    
  
  #  t=dsmerc.time.values
  #  timeout=(t-np.datetime64(rtime)) / np.timedelta64(1, 'D')
  #  ncid.variables['ocean_time'][:]=timeout
    
        
        
    fh.close()   
def create_init_file(initfile,grd,tunits):
    """
    create_init_file(initfile,grd)

    Create a child initilization file. 
    """

    print('Creating Initialization file')    
    # Create the initialization netCDF file with global attributes.
    fh        = nc.Dataset(initfile, "w", format="NETCDF4",clobber=True)
    fh.type   =  'Initial conditions file for L1 grid'
    fh.title  =  'Initial conditions' 
    fh.history = ['Created by create_init_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]



    dims=grd.sizes
    atts=grd.attrs
    #fh.createDimension('xpsi', dims['xi_psi'])
   # fh.createDimension('xi_psi', dims['xi_u'])
    fh.createDimension('xi_rho', dims['xi_rho'])
    fh.createDimension('xi_u', dims['xi_u'])
  #  fh.createDimension('xv', dims['xi_v'])
    fh.createDimension('xi_v', dims['xi_rho'])
    
 #   fh.createDimension('epsi', dims['eta_psi'])
   # fh.createDimension('eta_psi', dims['eta_v'])
    fh.createDimension('eta_rho', dims['eta_rho'])
  #  fh.createDimension('eu', dims['eta_u'])
    fh.createDimension('eta_u', dims['eta_rho'])
    fh.createDimension('eta_v', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    fh.createDimension('s_w', atts['sc_w'])
    
  #  fh.createDimension('sc_r', atts['sc_r'])
  # fh.createDimension('sc_w', atts['sc_w'])
    fh.createDimension('time',None)
    
        
# Defining hydrodynamic variables of initialization file.  

    fh.createVariable('lon_rho', 'f', ('eta_rho', 'xi_rho'))
    fh.variables['lon_rho'].long_name = 'longitude of RHO-points"'
    fh.variables['lon_rho'].units = 'degree_east'
    fh.variables['lon_rho'][:]=grd.lon_rho[:]
    
    fh.createVariable('lat_rho', 'f', ('eta_rho', 'xi_rho'))
    fh.variables['lat_rho'].long_name = 'latitude of RHO-points"'
    fh.variables['lat_rho'].units = 'degree_north'
    fh.variables['lat_rho'][:]=grd.lat_rho[:]
    
    fh.createVariable('lon_u', 'f', ('eta_u', 'xi_u'))
    fh.variables['lon_u'].long_name = 'longitude of U-points"'
    fh.variables['lon_u'].units = 'degree_east'
    fh.variables['lon_u'][:]=grd.lon_u[:]
    
    fh.createVariable('lat_u', 'f', ('eta_u', 'xi_u'))
    fh.variables['lat_u'].long_name = 'latitude of U-points"'
    fh.variables['lat_u'].units = 'degree_north'
    fh.variables['lat_u'][:]=grd.lat_u[:]

    fh.createVariable('lon_v', 'f', ('eta_v', 'xi_v'))
    fh.variables['lon_v'].long_name = 'longitude of V-points"'
    fh.variables['lon_v'].units = 'degree_east'
    fh.variables['lon_v'][:]=grd.lon_v[:]
    
    fh.createVariable('lat_v', 'f', ('eta_v', 'xi_v'))
    fh.variables['lat_v'].long_name = 'latitude of V-points"'
    fh.variables['lat_v'].units = 'degree_north'
    fh.variables['lat_v'][:]=grd.lat_v[:]


    fh.createVariable('spherical', 'short', ('time'))
    fh.variables['spherical'].long_name = 'grid type logical switch'
    fh.variables['spherical'].flag_meanings = 'Cartesian spherical"'
    fh.variables['spherical'].flag_values = '0, 1'

    fh.createVariable('Vtransform','long',('time'));
    fh.variables['Vtransform'].long_name ='vertical terrain-following transformation equation'

    fh.createVariable('Vstretching','long',('time'));
    fh.variables['Vstretching'].long_name ='vertical terrain-following stretching function'

    fh.createVariable('theta_b','double',('time'));
    fh.variables['theta_b'].long_name ='S-coordinate bottom control parameter'
    fh.variables['theta_b'].units ='1'
    
    fh.createVariable('theta_s','double',('time'));
    fh.variables['theta_s'].long_name ='S-coordinate surface control parameter'
    fh.variables['theta_s'].units ='1'
    
    fh.createVariable('Tcline','double',('time'));
    fh.variables['Tcline'].long_name ='S-coordinate surface/bottom layer width'
    fh.variables['Tcline'].units ='meter'
    
    fh.createVariable('hc','double',('time'));
    fh.variables['hc'].long_name ='S-coordinate parameter, critical depth'
    fh.variables['hc'].units ='meter'
    
    fh.createVariable('Cs_r','double',('s_rho'));
    fh.variables['Cs_r'].long_name='S-coordinate stretching curves at RHO-points'
    fh.variables['Cs_r'].valid_min =-1
    fh.variables['Cs_r'].valid_max =0
    
    fh.createVariable('Cs_w','double',('s_w'));
    fh.variables['Cs_w'].long_name='S-coordinate stretching curves at W-points'
    fh.variables['Cs_w'].valid_min =-1
    fh.variables['Cs_w'].valid_max =0
    
    fh.createVariable('s_rho','double',('s_rho'));
    fh.variables['s_rho'].long_name='S-coordinate at RHO-points'
    fh.variables['s_rho'].valid_min =-1
    fh.variables['s_rho'].valid_max =0
    fh.variables['s_rho'].positive ="up" 
    fh.variables['s_rho'].standard_name ="ocean_s_coordinate_g2"
    fh.variables['s_rho'].formula_terms="s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" 
            
    fh.createVariable('s_w','double',('s_w'));
    fh.variables['s_w'].long_name='S-coordinate at W-points'
    fh.variables['s_w'].valid_min =-1
    fh.variables['s_w'].valid_max =0
    fh.variables['s_w'].positive ="up"
    fh.variables['s_w'].standard_name = "ocean_s_coordinate_g2" 
    fh.variables['s_w'].formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"  

    fh.createVariable('ocean_time','double',('time'));
    fh.variables['ocean_time'].long_name='time since initialization'
    fh.variables['ocean_time'].units =tunits
    fh.variables['ocean_time'].field ='ocean_time, scalar, series'
    
    fh.createVariable('salt','f',('time','s_rho','eta_rho','xi_rho'));
    fh.variables['salt'].long_name='salinity'
    fh.variables['salt'].units ='PSU'
    fh.variables['salt'].time ="ocean_time" 
    fh.variables['salt'].coordinates ="lon_rho lat_rho ocean_time" 
    
    fh.createVariable('temp','f',('time','s_rho','eta_rho','xi_rho'));
    fh.variables['temp'].long_name='temperature'
    fh.variables['temp'].units ='C'
    fh.variables['temp'].time ="ocean_time" 
    fh.variables['temp'].coordinates ="lon_rho lat_rho ocean_time" 

    fh.createVariable('u','f',('time','s_rho','eta_u','xi_u'));
    fh.variables['u'].long_name='u-momentum component'
    fh.variables['u'].units ='meter second-1'
    fh.variables['u'].time ="ocean_time" 
    fh.variables['u'].coordinates ="lon_rho lat_rho ocean_time" 

    fh.createVariable('ubar','f',('time','eta_u','xi_u'));
    fh.variables['ubar'].long_name='vertically integrated u-momentum component'
    fh.variables['ubar'].units ='meter second-1'
    fh.variables['ubar'].time ="ocean_time" 
    fh.variables['ubar'].coordinates ="lon_rho lat_rho ocean_time" 
    
    fh.createVariable('v','f',('time','s_rho','eta_v','xi_v'));
    fh.variables['v'].long_name='v-momentum component'
    fh.variables['v'].units ='meter second-1'
    fh.variables['v'].time ="ocean_time" 
    fh.variables['v'].coordinates ="lon_rho lat_rho ocean_time" 

    fh.createVariable('vbar','f',('time','eta_v','xi_v'));
    fh.variables['vbar'].long_name='vertically integrated v-momentum component'
    fh.variables['vbar'].units ='meter second-1'
    fh.variables['vbar'].time = "ocean_time" 
    fh.variables['vbar'].coordinates = "lon_rho lat_rho ocean_time" 
    
    fh.createVariable('zeta','f',('time','eta_rho','xi_rho'));
    fh.variables['zeta'].long_name='free-surface'
    fh.variables['zeta'].units ='meter'
    fh.variables['zeta'].time = "ocean_time" 
    fh.variables['zeta'].coordinates = "lon_rho lat_rho ocean_time"  
    
    
    # Defining Sediment variables of initialization file.  
    
        
        
    fh.close()
        


def create_init_file_COAWST(initfile,grd,tunits):
    """
    create_init_file(initfile,grd)

    Create a child initilization file. 
    """

    print('Creating Initialization file')    
    # Create the initialization netCDF file with global attributes.
    fh        = nc.Dataset(initfile, "w", format="NETCDF4")
    fh.type   =  'Initial conditions file for L1 grid'
    fh.title  =  'Initial conditions' 
    fh.history = ['Created by create_init_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]



    dims=grd.sizes
    atts=grd.attrs
    #fh.createDimension('xpsi', dims['xi_psi'])
    fh.createDimension('xpsi', dims['xi_u'])
    fh.createDimension('xrho', dims['xi_rho'])
    fh.createDimension('xu', dims['xi_u'])
  #  fh.createDimension('xv', dims['xi_v'])
    fh.createDimension('xv', dims['xi_rho'])
    
 #   fh.createDimension('epsi', dims['eta_psi'])
    fh.createDimension('epsi', dims['eta_v'])
    fh.createDimension('erho', dims['eta_rho'])
  #  fh.createDimension('eu', dims['eta_u'])
    fh.createDimension('eu', dims['eta_rho'])
    fh.createDimension('ev', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    
    fh.createDimension('sc_r', atts['sc_r'])
    fh.createDimension('sc_w', atts['sc_w'])
    fh.createDimension('time',None)
    
    if  atts['NNS']+atts['NCS']>0:
        fh.createDimension('Nbed', atts['Nbed'])
        
    if  atts['Nveg']>0:
        fh.createDimension('Nveg', atts['Nveg'])
        
# Defining hydrodynamic variables of initialization file.  

    fh.createVariable('lon_rho', 'double', ('erho', 'xrho'))
    fh.variables['lon_rho'].long_name = 'longitude of RHO-points"'
    fh.variables['lon_rho'].units = 'degree_east'
    fh.variables['lon_rho'][:]=grd.lon_rho[:]
    
    fh.createVariable('lat_rho', 'double', ('erho', 'xrho'))
    fh.variables['lat_rho'].long_name = 'latitude of RHO-points"'
    fh.variables['lat_rho'].units = 'degree_north'
    fh.variables['lat_rho'][:]=grd.lat_rho[:]
    
    fh.createVariable('lon_u', 'double', ('eu', 'xu'))
    fh.variables['lon_u'].long_name = 'longitude of U-points"'
    fh.variables['lon_u'].units = 'degree_east'
    fh.variables['lon_u'][:]=grd.lon_u[:]
    
    fh.createVariable('lat_u', 'double', ('eu', 'xu'))
    fh.variables['lat_u'].long_name = 'latitude of U-points"'
    fh.variables['lat_u'].units = 'degree_north'
    fh.variables['lat_u'][:]=grd.lat_u[:]

    fh.createVariable('lon_v', 'double', ('ev', 'xv'))
    fh.variables['lon_v'].long_name = 'longitude of V-points"'
    fh.variables['lon_v'].units = 'degree_east'
    fh.variables['lon_v'][:]=grd.lon_v[:]
    
    fh.createVariable('lat_v', 'double', ('ev', 'xv'))
    fh.variables['lat_v'].long_name = 'latitude of V-points"'
    fh.variables['lat_v'].units = 'degree_north'
    fh.variables['lat_v'][:]=grd.lat_v[:]


    
    fh.createVariable('spherical', 'short', ('time'))
    fh.variables['spherical'].long_name = 'grid type logical switch'
    fh.variables['spherical'].flag_meanings = 'spherical, Cartesian'
    fh.variables['spherical'].flag_values = '1, 0'

    fh.createVariable('Vtransform','long',('time'));
    fh.variables['Vtransform'].long_name ='vertical terrain-following transformation equation'

    fh.createVariable('Vstretching','long',('time'));
    fh.variables['Vstretching'].long_name ='vertical terrain-following stretching function'

    fh.createVariable('theta_b','double',('time'));
    fh.variables['theta_b'].long_name ='S-coordinate bottom control parameter'
    fh.variables['theta_b'].units ='1'
    
    fh.createVariable('theta_s','double',('time'));
    fh.variables['theta_s'].long_name ='S-coordinate surface control parameter'
    fh.variables['theta_s'].units ='1'
    
    fh.createVariable('Tcline','double',('time'));
    fh.variables['Tcline'].long_name ='S-coordinate surface/bottom layer width'
    fh.variables['Tcline'].units ='meter'
    
    fh.createVariable('hc','double',('time'));
    fh.variables['hc'].long_name ='S-coordinate parameter, critical depth'
    fh.variables['hc'].units ='meter'
    
    fh.createVariable('Cs_r','double',('sc_r'));
    fh.variables['Cs_r'].long_name='S-coordinate stretching curves at RHO-points'
    fh.variables['Cs_r'].units ='1'
    fh.variables['Cs_r'].valid_min =-1
    fh.variables['Cs_r'].valid_max =0
    fh.variables['Cs_r'].field ='Cs_r, scalar'
    
    fh.createVariable('Cs_w','double',('sc_w'));
    fh.variables['Cs_w'].long_name='S-coordinate stretching curves at W-points'
    fh.variables['Cs_w'].units ='1'
    fh.variables['Cs_w'].valid_min =-1
    fh.variables['Cs_w'].valid_max =0
    fh.variables['Cs_w'].field ='Cs_w, scalar'
    
    fh.createVariable('sc_r','double',('sc_r'));
    fh.variables['sc_r'].long_name='S-coordinate at RHO-points'
    fh.variables['sc_r'].units ='1'
    fh.variables['sc_r'].valid_min =-1
    fh.variables['sc_r'].valid_max =0
    fh.variables['sc_r'].field ='sc_r, scalar'
    
    fh.createVariable('sc_w','double',('sc_w'));
    fh.variables['sc_w'].long_name='S-coordinate at W-points'
    fh.variables['sc_w'].units ='1'
    fh.variables['sc_w'].valid_min =-1
    fh.variables['sc_w'].valid_max =0
    fh.variables['sc_w'].field ='sc_w, scalar'

    fh.createVariable('ocean_time','double',('time'));
    fh.variables['ocean_time'].long_name='time since initialization'
    fh.variables['ocean_time'].units =tunits
    fh.variables['ocean_time'].field ='ocean_time, scalar, series'
    
    fh.createVariable('salt','float',('time','sc_r','erho','xrho'));
    fh.variables['salt'].long_name='salinity'
    fh.variables['salt'].units ='PSU'
    fh.variables['salt'].field ='salinity, scalar, series'
    
    fh.createVariable('temp','float',('time','sc_r','erho','xrho'));
    fh.variables['temp'].long_name='temperature'
    fh.variables['temp'].units ='C'
    fh.variables['temp'].field ='temperature, scalar, series'

    fh.createVariable('u','float',('time','sc_r','eu','xu'));
    fh.variables['u'].long_name='u-momentum component'
    fh.variables['u'].units ='meter second-1'
    fh.variables['u'].field ='u-velocity, scalar, series'

    fh.createVariable('ubar','float',('time','eu','xu'));
    fh.variables['ubar'].long_name='vertically integrated u-momentum component'
    fh.variables['ubar'].units ='meter second-1'
    fh.variables['ubar'].field ='ubar-velocity, scalar, series'
    
    fh.createVariable('v','float',('time','sc_r','ev','xv'));
    fh.variables['v'].long_name='v-momentum component'
    fh.variables['v'].units ='meter second-1'
    fh.variables['v'].field ='v-velocity, scalar, series'

    fh.createVariable('vbar','float',('time','ev','xv'));
    fh.variables['vbar'].long_name='vertically integrated v-momentum component'
    fh.variables['vbar'].units ='meter second-1'
    fh.variables['vbar'].field ='vbar-velocity, scalar, series'
    
    fh.createVariable('zeta','float',('time','erho','xrho'));
    fh.variables['zeta'].long_name='free-surface'
    fh.variables['zeta'].units ='meter'
    fh.variables['zeta'].field ='free-surface, scalar, series'
    
    
    # Defining Sediment variables of initialization file.  
    
    
    for mm in range(atts['NCS']):
        ind=mm+1
        sind='{:02d}'.format(ind)
        fh.createVariable('mud_'+sind,'double',('time','sc_r','erho','xrho'))
        fh.variables['mud_'+sind].long_name='suspended cohesive sediment, size class '+sind
        fh.variables['mud_'+sind].units ='kilogram meter-3'
        fh.variables['mud_'+sind].time ='ocean_time'
        fh.variables['mud_'+sind].field ='mud_'+sind+', scalar, series'     
        
        fh.createVariable('mudfrac_'+sind,'double',('time','Nbed','erho','xrho'))
        fh.variables['mudfrac_'+sind].long_name='cohesive sediment fraction, size class '+sind
        fh.variables['mudfrac_'+sind].units ='nondimensional'
        fh.variables['mudfrac_'+sind].time ='ocean_time'
        fh.variables['mudfrac_'+sind].field ='mudfrac_'+sind+', scalar, series'     
        
        fh.createVariable('mudmass_'+sind,'double',('time','Nbed','erho','xrho'))
        fh.variables['mudmass_'+sind].long_name='cohesive sediment mass, size class '+sind
        fh.variables['mudmass_'+sind].units ='kilogram meter-3'
        fh.variables['mudmass_'+sind].time ='ocean_time'
        fh.variables['mudmass_'+sind].field ='mudmass_'+sind+', scalar, series'
        
        
        
    
    for mm in range(atts['NNS']):
        ind=mm+1
        sind='{:02d}'.format(ind)
    
        fh.createVariable('sand_'+sind,'double',('time','sc_r','erho','xrho'))
        fh.variables['sand_'+sind].long_name='suspended noncohesive sediment, size class '+sind
        fh.variables['sand_'+sind].units ='kilogram meter-3'
        fh.variables['sand_'+sind].time ='ocean_time'
        fh.variables['sand_'+sind].field ='sand_'+sind+', scalar, series'     
        
        fh.createVariable('sandfrac_'+sind,'double',('time','Nbed','erho','xrho'))
        fh.variables['sandfrac_'+sind].long_name='noncohesive sediment fraction, size class '+sind
        fh.variables['sandfrac_'+sind].units ='nondimensional'
        fh.variables['sandfrac_'+sind].time ='ocean_time'
        fh.variables['sandfrac_'+sind].field ='sandfrac_'+sind+', scalar, series'     
        
        fh.createVariable('sandmass_'+sind,'double',('time','Nbed','erho','xrho'))
        fh.variables['sandmass_'+sind].long_name='noncohesive sediment mass, size class '+sind
        fh.variables['sandmass_'+sind].units ='kilogram meter-2'
        fh.variables['sandmass_'+sind].time ='ocean_time'
        fh.variables['sandmass_'+sind].field ='sandmass_'+sind+', scalar, series'
        
        fh.createVariable('bedload_Usand_'+sind,'double',('time','eu','xu'))
        fh.variables['bedload_Usand_'+sind].long_name='bed load flux of sand in U-direction, size class '+sind
        fh.variables['bedload_Usand_'+sind].units ='kilogram meter-1 s-1'
        fh.variables['bedload_Usand_'+sind].time ='ocean_time'
        fh.variables['bedload_Usand_'+sind].field ='bedload_Usand_'+sind+', scalar, series'
        
        fh.createVariable('bedload_Vsand_'+sind,'double',('time','ev','xv'))
        fh.variables['bedload_Vsand_'+sind].long_name='bed load flux of sand in V-direction, size class '+sind
        fh.variables['bedload_Vsand_'+sind].units ='kilogram meter-1 s-1'
        fh.variables['bedload_Vsand_'+sind].time ='ocean_time'
        fh.variables['bedload_Vsand_'+sind].field ='bedload_Vsand_'+sind+', scalar, series'
        
    
    
    
    
    if  atts['NNS']+atts['NCS']>0:
        
        fh.createVariable('bed_thickness','double',('time','Nbed','erho','xrho'))
        fh.variables['bed_thickness'].long_name='sediment layer thickness'
        fh.variables['bed_thickness'].units ='meter'
        fh.variables['bed_thickness'].time ='ocean_time'
        fh.variables['bed_thickness'].field ='bed_thickness, scalar, series'  
        
        fh.createVariable('bed_age','double',('time','Nbed','erho','xrho'))
        fh.variables['bed_age'].long_name='sediment layer age'
        fh.variables['bed_age'].units ='day'
        fh.variables['bed_age'].time ='ocean_time'
        fh.variables['bed_age'].field ='bed_age, scalar, series'     
        
        fh.createVariable('bed_porosity','double',('time','Nbed','erho','xrho'))
        fh.variables['bed_porosity'].long_name='sediment layer porosity'
        fh.variables['bed_porosity'].units ='nondimensional'
        fh.variables['bed_porosity'].time ='ocean_time'
        fh.variables['bed_porosity'].field ='bed_porosity, scalar, series'   
        
        fh.createVariable('bed_biodiff','double',('time','Nbed','erho','xrho'))
        fh.variables['bed_biodiff'].long_name='biodiffusivity at bottom of each layer'
        fh.variables['bed_biodiff'].units ='meter2 second-1'
        fh.variables['bed_biodiff'].time ='ocean_time'
        fh.variables['bed_biodiff'].field ='bed_biodiff, scalar, series'   
           
        fh.createVariable('grain_diameter','double',('time','erho','xrho'))
        fh.variables['grain_diameter'].long_name='sediment median grain diameter size'
        fh.variables['grain_diameter'].units ='meter'
        fh.variables['grain_diameter'].time ='ocean_time'
        fh.variables['grain_diameter'].field ='grain_diameter, scalar, series'   
           
        fh.createVariable('grain_density','double',('time','erho','xrho'))
        fh.variables['grain_density'].long_name='sediment median grain density'
        fh.variables['grain_density'].units ='kilogram meter-3'
        fh.variables['grain_density'].time ='ocean_time'
        fh.variables['grain_density'].field ='grain_density, scalar, series'   
           
        fh.createVariable('settling_vel','double',('time','erho','xrho'))
        fh.variables['settling_vel'].long_name='sediment median grain settling velocity'
        fh.variables['settling_vel'].units ='meter second-1'
        fh.variables['settling_vel'].time ='ocean_time'
        fh.variables['settling_vel'].field ='settling_vel, scalar, series'   
           
        fh.createVariable('erosion_stress','double',('time','erho','xrho'))
        fh.variables['erosion_stress'].long_name='sediment median critical erosion stress'
        fh.variables['erosion_stress'].units ='meter2 second-2'
        fh.variables['erosion_stress'].time ='ocean_time'
        fh.variables['erosion_stress'].field ='erosion_stress, scalar, series'   
           
        fh.createVariable('ripple_length','double',('time','erho','xrho'))
        fh.variables['ripple_length'].long_name='bottom ripple length'
        fh.variables['ripple_length'].units ='meter'
        fh.variables['ripple_length'].time ='ocean_time'
        fh.variables['ripple_length'].field ='ripple_length, scalar, series'   
           
        fh.createVariable('ripple_height','double',('time','erho','xrho'))
        fh.variables['ripple_height'].long_name='bottom ripple height'
        fh.variables['ripple_height'].units ='meter'
        fh.variables['ripple_height'].time ='ocean_time'
        fh.variables['ripple_height'].field ='ripple_height, scalar, series'   
           
        fh.createVariable('dmix_offset','double',('time','erho','xrho'))
        fh.variables['dmix_offset'].long_name='dmix erodibility profile offset'
        fh.variables['dmix_offset'].units ='meter'
        fh.variables['dmix_offset'].time ='ocean_time'
        fh.variables['dmix_offset'].field ='dmix_offset, scalar, series'   
           
        fh.createVariable('dmix_slope','double',('time','erho','xrho'))
        fh.variables['dmix_slope'].long_name='dmix erodibility profile slope'
        fh.variables['dmix_slope'].units ='_'
        fh.variables['dmix_slope'].time ='ocean_time'
        fh.variables['dmix_slope'].field ='dmix_slope, scalar, series'     
           
        fh.createVariable('dmix_time','double',('time','erho','xrho'))
        fh.variables['dmix_time'].long_name='dmix erodibility profile time scale'
        fh.variables['dmix_time'].units ='seconds'
        fh.variables['dmix_time'].time ='ocean_time'
        fh.variables['dmix_time'].field ='dmix_time, scalar, series'              
                 
        
    if  atts['Nveg']>0: 
        fh.createVariable('plant_height','double',('time','Nveg','erho','xrho'))
        fh.variables['plant_height'].long_name='plant height'
        fh.variables['plant_height'].units ='meter'
        fh.variables['plant_height'].time ='ocean_time'
        fh.variables['plant_height'].field ='plant_height, scalar, series'  
        
        fh.createVariable('plant_density','double',('time','Nveg','erho','xrho'))
        fh.variables['plant_density'].long_name='plant density'
        fh.variables['plant_density'].units ='plant-meter2'
        fh.variables['plant_density'].time ='ocean_time'
        fh.variables['plant_density'].field ='plant_density, scalar, series' 
        
        fh.createVariable('plant_diameter','double',('time','Nveg','erho','xrho'))
        fh.variables['plant_diameter'].long_name='plant diameter'
        fh.variables['plant_diameter'].units ='meter'
        fh.variables['plant_diameter'].time ='ocean_time'
        fh.variables['plant_diameter'].field ='plant_diameter, scalar, series'  
        
        fh.createVariable('plant_thickness','double',('time','Nveg','erho','xrho'))
        fh.variables['plant_thickness'].long_name='plant thickness'
        fh.variables['plant_thickness'].units ='meter'
        fh.variables['plant_thickness'].time ='ocean_time'
        fh.variables['plant_thickness'].field ='plant_thickness, scalar, series'  
        
        fh.createVariable('marsh_mask','double',('time','Nveg','erho','xrho'))
        fh.variables['marsh_mask'].long_name='marsh mask'
        fh.variables['marsh_mask'].units ='nondimensional'
        fh.variables['marsh_mask'].time ='ocean_time'
        fh.variables['marsh_mask'].field ='marsh_mask, scalar, series'  
        
        
        
        
        
    fh.close()
        
    
def create_bdry_file_COAWST(initfile,grd,tunits):
    """
    create_bdry_file(initfile,grd)

    Create a child boundary conditions file file. 
    """
     
 
    print('Creating Boundary conditions files file')    
 # Create the initialization netCDF file with global attributes.
    fh        = nc.Dataset(initfile, "w", format="NETCDF4")
    fh.type   =  'Boundary conditions file for L1 grid'
    fh.title  =  'Initial conditions' 
    fh.history = ['Created by create_dbry_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]



 #Define the dimensions of the file
    dims=grd.sizes
    atts=grd.attrs
    
    #fh.createDimension('xpsi', dims['xi_psi'])
    fh.createDimension('xpsi', dims['xi_u'])
    fh.createDimension('xrho', dims['xi_rho'])
    fh.createDimension('xu', dims['xi_u'])
  #  fh.createDimension('xv', dims['xi_v'])
    fh.createDimension('xv', dims['xi_rho'])
    
 #   fh.createDimension('epsi', dims['eta_psi'])
    fh.createDimension('epsi', dims['eta_v'])
    fh.createDimension('erho', dims['eta_rho'])
  #  fh.createDimension('eu', dims['eta_u'])
    fh.createDimension('eu', dims['eta_rho'])
    fh.createDimension('ev', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    
    # fh.createDimension('zeta_time', atts['t2d'])
    # fh.createDimension('v2d_time', atts['t2d'])
    # fh.createDimension('v3d_time', atts['t3d'])
    # fh.createDimension('salt_time', atts['t3d'])
    # fh.createDimension('temp_time', atts['t3d'])
 
    fh.createDimension('zeta_time',None)
    fh.createDimension('v2d_time',None)
    fh.createDimension('v3d_time',None)
    fh.createDimension('salt_time', None)
    fh.createDimension('temp_time', None)
 #Define the Variabkes in the file

    fh.createVariable('zeta_time','double',('zeta_time'));
    fh.variables['zeta_time'].long_name='zeta_time'
    fh.variables['zeta_time'].units =tunits
    fh.variables['zeta_time'].field ='zeta_time, scalar, series'

    fh.createVariable('v2d_time','double',('v2d_time'));
    fh.variables['v2d_time'].long_name='v2d_time'
    fh.variables['v2d_time'].units =tunits
    fh.variables['v2d_time'].field ='v2d_time, scalar, series'

    fh.createVariable('v3d_time','double',('v3d_time'));
    fh.variables['v3d_time'].long_name='v3d_time'
    fh.variables['v3d_time'].units =tunits
    fh.variables['v3d_time'].field ='v3d_time, scalar, series'

    fh.createVariable('salt_time','double',('salt_time'));
    fh.variables['salt_time'].long_name='salt_time'
    fh.variables['salt_time'].units =tunits
    fh.variables['salt_time'].field ='salt_time, scalar, series'

    fh.createVariable('temp_time','double',('temp_time'));
    fh.variables['temp_time'].long_name='temp_time'
    fh.variables['temp_time'].units =tunits
    fh.variables['temp_time'].field ='temp_time, scalar, series'
    
    fh.createVariable('zeta_south','f8',('zeta_time','xrho'));
    fh.variables['zeta_south'].long_name='free-surface southern boundary condition'
    fh.variables['zeta_south'].units ='meter'
    fh.variables['zeta_south'].field ='zeta_south, scalar, series'
    
    fh.createVariable('zeta_east','f8',('zeta_time','erho'));
    fh.variables['zeta_east'].long_name='free-surface eastern boundary condition'
    fh.variables['zeta_east'].units ='meter'
    fh.variables['zeta_east'].field ='zeta_east, scalar, series'
    
    fh.createVariable('zeta_west','f8',('zeta_time','erho'));
    fh.variables['zeta_west'].long_name='free-surface western boundary condition'
    fh.variables['zeta_west'].units ='meter'
    fh.variables['zeta_west'].field ='zeta_west, scalar, series'
    
    fh.createVariable('zeta_north','f8',('zeta_time','xrho'));
    fh.variables['zeta_north'].long_name='free-surface northern boundary condition'
    fh.variables['zeta_north'].units ='meter'
    fh.variables['zeta_north'].field ='zeta_north, scalar, series'
    
    fh.createVariable('ubar_south','f4',('v2d_time','xu'));
    fh.variables['ubar_south'].long_name='2D u-momentum southern boundary condition'
    fh.variables['ubar_south'].units ='meter second-1'
    fh.variables['ubar_south'].field ='ubar_south, scalar, series'
    
    fh.createVariable('ubar_east','f4',('v2d_time','eu'));
    fh.variables['ubar_east'].long_name='2D u-momentum eastern boundary condition'
    fh.variables['ubar_east'].units ='meter second-1'
    fh.variables['ubar_east'].field ='ubar_east, scalar, series'
    
    fh.createVariable('ubar_west','f4',('v2d_time','eu'));
    fh.variables['ubar_west'].long_name='2D u-momentum western boundary condition'
    fh.variables['ubar_west'].units ='meter second-1'
    fh.variables['ubar_west'].field ='ubar_west, scalar, series'

    fh.createVariable('ubar_north','f4',('v2d_time','xu'));
    fh.variables['ubar_north'].long_name='2D u-momentum northern boundary condition'
    fh.variables['ubar_north'].units ='meter second-1'
    fh.variables['ubar_north'].field ='ubar_north, scalar, series'
    
    fh.createVariable('vbar_south','f4',('v2d_time','xv'));
    fh.variables['vbar_south'].long_name='2D v-momentum southern boundary condition'
    fh.variables['vbar_south'].units ='meter second-1'
    fh.variables['vbar_south'].field ='vbar_south, scalar, series'
    
    fh.createVariable('vbar_east','f4',('v2d_time','ev'));
    fh.variables['vbar_east'].long_name='2D v-momentum eastern boundary condition'
    fh.variables['vbar_east'].units ='meter second-1'
    fh.variables['vbar_east'].field ='vbar_east, scalar, series'
    
    fh.createVariable('vbar_west','f4',('v2d_time','ev'));
    fh.variables['vbar_west'].long_name='2D v-momentum western boundary condition'
    fh.variables['vbar_west'].units ='meter second-1'
    fh.variables['vbar_west'].field ='vbar_west, scalar, series'

    fh.createVariable('vbar_north','f4',('v2d_time','xv'));
    fh.variables['vbar_north'].long_name='2D v-momentum northern boundary condition'
    fh.variables['vbar_north'].units ='meter second-1'
    fh.variables['vbar_north'].field ='vbar_north, scalar, series'

    fh.createVariable('u_south','f4',('v3d_time','s_rho','xu'));
    fh.variables['u_south'].long_name='3D u-momentum southern boundary condition'
    fh.variables['u_south'].units ='meter second-1'
    fh.variables['u_south'].field ='u_south, scalar, series'
    
    fh.createVariable('u_east','f4',('v3d_time','s_rho','eu'));
    fh.variables['u_east'].long_name='3D u-momentum eastern boundary condition'
    fh.variables['u_east'].units ='meter second-1'
    fh.variables['u_east'].field ='u_east, scalar, series'
    
    fh.createVariable('u_west','f4',('v3d_time','s_rho','eu'));
    fh.variables['u_west'].long_name='3D u-momentum western boundary condition'
    fh.variables['u_west'].units ='meter second-1'
    fh.variables['u_west'].field ='u_west, scalar, series'

    fh.createVariable('u_north','f4',('v3d_time','s_rho','xu'));
    fh.variables['u_north'].long_name='3D u-momentum northern boundary condition'
    fh.variables['u_north'].units ='meter second-1'
    fh.variables['u_north'].field ='u_north, scalar, series'

    fh.createVariable('v_south','f4',('v3d_time','s_rho','xv'));
    fh.variables['v_south'].long_name='3D v-momentum southern boundary condition'
    fh.variables['v_south'].units ='meter second-1'
    fh.variables['v_south'].field ='v_south, scalar, series'
    
    fh.createVariable('v_east','f4',('v3d_time','s_rho','ev'));
    fh.variables['v_east'].long_name='3D v-momentum eastern boundary condition'
    fh.variables['v_east'].units ='meter second-1'
    fh.variables['v_east'].field ='v_east, scalar, series'
    
    fh.createVariable('v_west','f4',('v3d_time','s_rho','ev'));
    fh.variables['v_west'].long_name='3D v-momentum western boundary condition'
    fh.variables['v_west'].units ='meter second-1'
    fh.variables['v_west'].field ='v_west, scalar, series'

    fh.createVariable('v_north','f4',('v3d_time','s_rho','xv'));
    fh.variables['v_north'].long_name='3D v-momentum northern boundary condition'
    fh.variables['v_north'].units ='meter second-1'
    fh.variables['v_north'].field ='v_north, scalar, series'

    fh.createVariable('temp_south','f4',('temp_time','s_rho','xrho'));
    fh.variables['temp_south'].long_name='3D temperature southern boundary condition'
    fh.variables['temp_south'].units ='C'
    fh.variables['temp_south'].field ='temp_south, scalar, series'
    
    fh.createVariable('temp_east','f4',('temp_time','s_rho','erho'));
    fh.variables['temp_east'].long_name='3D temperature eastern boundary condition'
    fh.variables['temp_east'].units ='C'
    fh.variables['temp_east'].field ='temp_east, scalar, series'
    
    fh.createVariable('temp_west','f4',('temp_time','s_rho','erho'));
    fh.variables['temp_west'].long_name='3D temperature western boundary condition'
    fh.variables['temp_west'].units ='C'
    fh.variables['temp_west'].field ='temp_west, scalar, series'

    fh.createVariable('temp_north','f4',('temp_time','s_rho','xrho'));
    fh.variables['temp_north'].long_name='3D temperature northern boundary condition'
    fh.variables['temp_north'].units ='C'
    fh.variables['temp_north'].field ='temp_north, scalar, series'
    

    fh.createVariable('salt_south','f4',('salt_time','s_rho','xrho'));
    fh.variables['salt_south'].long_name='3D salinity southern boundary condition'
    fh.variables['salt_south'].units ='psu'
    fh.variables['salt_south'].field ='salt_south, scalar, series'
    
    fh.createVariable('salt_east','f4',('salt_time','s_rho','erho'));
    fh.variables['salt_east'].long_name='3D salinity eastern boundary condition'
    fh.variables['salt_east'].units ='psu'
    fh.variables['salt_east'].field ='salt_east, scalar, series'
    
    fh.createVariable('salt_west','f4',('salt_time','s_rho','erho'));
    fh.variables['salt_west'].long_name='3D salinity western boundary condition'
    fh.variables['salt_west'].units ='psu'
    fh.variables['salt_west'].field ='salt_west, scalar, series'

    fh.createVariable('salt_north','f4',('salt_time','s_rho','xrho'));
    fh.variables['salt_north'].long_name='3D salinity northern boundary condition'
    fh.variables['salt_north'].units ='psu'
    fh.variables['salt_north'].field ='salt_north, scalar, series'

    fh.close()    
    
def create_bdry_file(initfile,grd,tunits):
    """
    create_bdry_file(initfile,grd)

    Create a child boundary conditions file file. 
    """
     
 
    print('Creating Boundary conditions files file')    
 # Create the initialization netCDF file with global attributes.
    fh        = nc.Dataset(initfile, "w", format="NETCDF4")
    fh.type   =  'Boundary conditions'
    fh.title  =  'Boundary conditions file for ECCOFS grid' 
    fh.history = ['Created by create_dbry_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]



 #Define the dimensions of the file
    dims=grd.sizes
    atts=grd.attrs
    
    fh.createDimension('xu_u', dims['xi_u'])
    fh.createDimension('xi_rho', dims['xi_rho'])
    fh.createDimension('xi_u', dims['xi_u'])
    fh.createDimension('xi_v', dims['xi_rho'])
    
    fh.createDimension('eta_rho', dims['eta_rho'])
    fh.createDimension('eta_u', dims['eta_rho'])
    fh.createDimension('eta_v', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    fh.createDimension('s_w', atts['sc_w'])
    fh.createDimension('tracer', 2)
    

 
    fh.createDimension('bry_time',None)

 #Define the Variabkes in the file

    fh.createVariable('bry_time','double',('bry_time'));
    fh.variables['bry_time'].long_name='bry_time'
    fh.variables['bry_time'].units =tunits
    fh.variables['bry_time'].field ='zeta_time, scalar, series'

    fh.createVariable("spherical", "i4") 
    fh.variables['spherical'].long_name="grid type logical switch"
    fh.variables['spherical'].flag_values="0, 1"
    fh.variables['spherical'].flag_meanings="Cartesian spherical"
    fh.createVariable("Vtransform", "i4") 
    fh.variables['Vtransform'].long_name="vertical terrain-following transformation equation"
    fh.createVariable("Vstretching", "i4") 
    fh.variables['Vstretching'].long_name="vertical terrain-following stretching function"
    fh.createVariable("theta_s", "f8") 
    fh.variables['theta_s'].long_name="S-coordinate surface control parameter"
    fh.createVariable("theta_b", "f8") 
    fh.variables['theta_b'].long_name="S-coordinate bottom control parameter"
    fh.createVariable("Tcline", "f8") 
    fh.variables['Tcline'].long_name="S-coordinate surface/bottom layer width"
    fh.variables['Tcline'].units="meter"
    fh.createVariable("hc", "f8") 
    fh.variables['hc'].long_name="S-coordinate parameter, critical depth"
    fh.variables['hc'].units="meter"
    
    
    
    fh.createVariable("s_rho", "f8",('s_rho'))
    fh.variables["s_rho"].long_name = "S-coordinate at RHO-points" 
    fh.variables["s_rho"].valid_min = -1. 
    fh.variables["s_rho"].valid_max = 0. 
    fh.variables["s_rho"].positive = "up" 
    fh.variables["s_rho"].standard_name = "ocean_s_coordinate_g2" 
    fh.variables["s_rho"].formula_terms = "s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" 
    fh.createVariable("s_w", "f8",('s_w'))
    fh.variables["s_w"].long_name = "S-coordinate at W-points" 
    fh.variables["s_w"].valid_min = -1. 
    fh.variables["s_w"].valid_max = 0. 
    fh.variables["s_w"].positive = "up" 
    fh.variables["s_w"].standard_name = "ocean_s_coordinate_g2" 
    fh.variables["s_w"].formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"     
    fh.createVariable("Cs_r", "f8",('s_rho'))
    fh.variables["Cs_r"].long_name = "S-coordinate stretching function at RHO-point" 
    fh.variables["Cs_r"].valid_min = -1. 
    fh.variables["Cs_r"].valid_max = 0.   
    fh.createVariable("Cs_w", "f8",('s_w'))
    fh.variables["Cs_w"].long_name = "S-coordinate stretching function at W-points" 
    fh.variables["Cs_w"].valid_min = -1. 
    fh.variables["Cs_w"].valid_max = 0. 

  
    
    
    fh.createVariable('zeta_south','f8',('bry_time','xi_rho'));
    fh.variables['zeta_south'].long_name='free-surface southern boundary condition'
    fh.variables['zeta_south'].units ='meter'
    fh.variables['zeta_south'].field ='zeta_south, scalar, series'
    
    fh.createVariable('zeta_east','f8',('bry_time','eta_rho'));
    fh.variables['zeta_east'].long_name='free-surface eastern boundary condition'
    fh.variables['zeta_east'].units ='meter'
    fh.variables['zeta_east'].field ='zeta_east, scalar, series'
    
    fh.createVariable('zeta_west','f8',('bry_time','eta_rho'));
    fh.variables['zeta_west'].long_name='free-surface western boundary condition'
    fh.variables['zeta_west'].units ='meter'
    fh.variables['zeta_west'].field ='zeta_west, scalar, series'
    
    fh.createVariable('zeta_north','f8',('bry_time','xi_rho'));
    fh.variables['zeta_north'].long_name='free-surface northern boundary condition'
    fh.variables['zeta_north'].units ='meter'
    fh.variables['zeta_north'].field ='zeta_north, scalar, series'
    
    fh.createVariable('ubar_south','f4',('bry_time','xi_u'));
    fh.variables['ubar_south'].long_name='2D u-momentum southern boundary condition'
    fh.variables['ubar_south'].units ='meter second-1'
    fh.variables['ubar_south'].field ='ubar_south, scalar, series'
    
    fh.createVariable('ubar_east','f4',('bry_time','eta_u'));
    fh.variables['ubar_east'].long_name='2D u-momentum eastern boundary condition'
    fh.variables['ubar_east'].units ='meter second-1'
    fh.variables['ubar_east'].field ='ubar_east, scalar, series'
    
    fh.createVariable('ubar_west','f4',('bry_time','eta_u'));
    fh.variables['ubar_west'].long_name='2D u-momentum western boundary condition'
    fh.variables['ubar_west'].units ='meter second-1'
    fh.variables['ubar_west'].field ='ubar_west, scalar, series'

    fh.createVariable('ubar_north','f4',('bry_time','xi_u'));
    fh.variables['ubar_north'].long_name='2D u-momentum northern boundary condition'
    fh.variables['ubar_north'].units ='meter second-1'
    fh.variables['ubar_north'].field ='ubar_north, scalar, series'
    
    fh.createVariable('vbar_south','f4',('bry_time','xi_v'));
    fh.variables['vbar_south'].long_name='2D v-momentum southern boundary condition'
    fh.variables['vbar_south'].units ='meter second-1'
    fh.variables['vbar_south'].field ='vbar_south, scalar, series'
    
    fh.createVariable('vbar_east','f4',('bry_time','eta_v'));
    fh.variables['vbar_east'].long_name='2D v-momentum eastern boundary condition'
    fh.variables['vbar_east'].units ='meter second-1'
    fh.variables['vbar_east'].field ='vbar_east, scalar, series'
    
    fh.createVariable('vbar_west','f4',('bry_time','eta_v'));
    fh.variables['vbar_west'].long_name='2D v-momentum western boundary condition'
    fh.variables['vbar_west'].units ='meter second-1'
    fh.variables['vbar_west'].field ='vbar_west, scalar, series'

    fh.createVariable('vbar_north','f4',('bry_time','xi_v'));
    fh.variables['vbar_north'].long_name='2D v-momentum northern boundary condition'
    fh.variables['vbar_north'].units ='meter second-1'
    fh.variables['vbar_north'].field ='vbar_north, scalar, series'

    fh.createVariable('u_south','f4',('bry_time','s_rho','xi_u'));
    fh.variables['u_south'].long_name='3D u-momentum southern boundary condition'
    fh.variables['u_south'].units ='meter second-1'
    fh.variables['u_south'].field ='u_south, scalar, series'
    
    fh.createVariable('u_east','f4',('bry_time','s_rho','eta_u'));
    fh.variables['u_east'].long_name='3D u-momentum eastern boundary condition'
    fh.variables['u_east'].units ='meter second-1'
    fh.variables['u_east'].field ='u_east, scalar, series'
    
    fh.createVariable('u_west','f4',('bry_time','s_rho','eta_u'));
    fh.variables['u_west'].long_name='3D u-momentum western boundary condition'
    fh.variables['u_west'].units ='meter second-1'
    fh.variables['u_west'].field ='u_west, scalar, series'

    fh.createVariable('u_north','f4',('bry_time','s_rho','xi_u'));
    fh.variables['u_north'].long_name='3D u-momentum northern boundary condition'
    fh.variables['u_north'].units ='meter second-1'
    fh.variables['u_north'].field ='u_north, scalar, series'

    fh.createVariable('v_south','f4',('bry_time','s_rho','xi_v'));
    fh.variables['v_south'].long_name='3D v-momentum southern boundary condition'
    fh.variables['v_south'].units ='meter second-1'
    fh.variables['v_south'].field ='v_south, scalar, series'
    
    fh.createVariable('v_east','f4',('bry_time','s_rho','eta_v'));
    fh.variables['v_east'].long_name='3D v-momentum eastern boundary condition'
    fh.variables['v_east'].units ='meter second-1'
    fh.variables['v_east'].field ='v_east, scalar, series'
    
    fh.createVariable('v_west','f4',('bry_time','s_rho','eta_v'));
    fh.variables['v_west'].long_name='3D v-momentum western boundary condition'
    fh.variables['v_west'].units ='meter second-1'
    fh.variables['v_west'].field ='v_west, scalar, series'

    fh.createVariable('v_north','f4',('bry_time','s_rho','xi_v'));
    fh.variables['v_north'].long_name='3D v-momentum northern boundary condition'
    fh.variables['v_north'].units ='meter second-1'
    fh.variables['v_north'].field ='v_north, scalar, series'

    fh.createVariable('temp_south','f4',('bry_time','s_rho','xi_rho'));
    fh.variables['temp_south'].long_name='3D temperature southern boundary condition'
    fh.variables['temp_south'].units ='C'
    fh.variables['temp_south'].field ='temp_south, scalar, series'
    
    fh.createVariable('temp_east','f4',('bry_time','s_rho','eta_rho'));
    fh.variables['temp_east'].long_name='3D temperature eastern boundary condition'
    fh.variables['temp_east'].units ='C'
    fh.variables['temp_east'].field ='temp_east, scalar, series'
    
    fh.createVariable('temp_west','f4',('bry_time','s_rho','eta_rho'));
    fh.variables['temp_west'].long_name='3D temperature western boundary condition'
    fh.variables['temp_west'].units ='C'
    fh.variables['temp_west'].field ='temp_west, scalar, series'

    fh.createVariable('temp_north','f4',('bry_time','s_rho','xi_rho'));
    fh.variables['temp_north'].long_name='3D temperature northern boundary condition'
    fh.variables['temp_north'].units ='C'
    fh.variables['temp_north'].field ='temp_north, scalar, series'
    

    fh.createVariable('salt_south','f4',('bry_time','s_rho','xi_rho'));
    fh.variables['salt_south'].long_name='3D salinity southern boundary condition'
    fh.variables['salt_south'].units ='psu'
    fh.variables['salt_south'].field ='salt_south, scalar, series'
    
    fh.createVariable('salt_east','f4',('bry_time','s_rho','eta_rho'));
    fh.variables['salt_east'].long_name='3D salinity eastern boundary condition'
    fh.variables['salt_east'].units ='psu'
    fh.variables['salt_east'].field ='salt_east, scalar, series'
    
    fh.createVariable('salt_west','f4',('bry_time','s_rho','eta_rho'));
    fh.variables['salt_west'].long_name='3D salinity western boundary condition'
    fh.variables['salt_west'].units ='psu'
    fh.variables['salt_west'].field ='salt_west, scalar, series'

    fh.createVariable('salt_north','f4',('bry_time','s_rho','xi_rho'));
    fh.variables['salt_north'].long_name='3D salinity northern boundary condition'
    fh.variables['salt_north'].units ='psu'
    fh.variables['salt_north'].field ='salt_north, scalar, series'

    fh.close()
def create_clm_file_COAWST(initfile,grd,tunits):
    """
    create_clm_file(initfile)

    Create a child climatology file. 
    """

    print('Creating Climatology file')    
    
    # Create the initialization netCDF file with global attributes.
    fh        = nc.Dataset(initfile, "w", format="NETCDF4")
    fh.type   =  'Climatology file for L1 grid'
    fh.title  =  'Climatology' 
    fh.history = ['Created by create_clm_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]


    dims=grd.sizes
    atts=grd.attrs
    
    #Define the dimensions of the file
    #fh.createDimension('xpsi', dims['xi_psi'])
    fh.createDimension('xpsi', dims['xi_u'])
    fh.createDimension('xrho', dims['xi_rho'])
    fh.createDimension('xu', dims['xi_u'])
  #  fh.createDimension('xv', dims['xi_v'])
    fh.createDimension('xv', dims['xi_rho'])
    
 #   fh.createDimension('epsi', dims['eta_psi'])
    fh.createDimension('epsi', dims['eta_v'])
    fh.createDimension('erho', dims['eta_rho'])
  #  fh.createDimension('eu', dims['eta_u'])
    fh.createDimension('eu', dims['eta_rho'])
    fh.createDimension('ev', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    
    fh.createDimension('ocean_time', None)
    fh.createDimension('zeta_time', None)
    fh.createDimension('v2d_time', None)
    fh.createDimension('v3d_time', None)
    fh.createDimension('salt_time', None)
    fh.createDimension('temp_time', None)



    fh.createDimension('one', 1)
    
    #Define the Variabkes in the file
    fh.createVariable('ocean_time','double',('ocean_time'));
    fh.variables['ocean_time'].long_name='wind field time'
    fh.variables['ocean_time'].units =tunits
    fh.variables['ocean_time'].field ='ocean_time, scalar, series'
    
    fh.createVariable('zeta_time','double',('zeta_time'));
    fh.variables['zeta_time'].long_name='zeta_time'
    fh.variables['zeta_time'].units =tunits
    fh.variables['zeta_time'].field ='zeta_time, scalar, series'

    fh.createVariable('v2d_time','double',('v2d_time'));
    fh.variables['v2d_time'].long_name='v2d_time'
    fh.variables['v2d_time'].units =tunits
    fh.variables['v2d_time'].field ='v2d_time, scalar, series'

    fh.createVariable('v3d_time','double',('v3d_time'));
    fh.variables['v3d_time'].long_name='v3d_time'
    fh.variables['v3d_time'].units =tunits
    fh.variables['v3d_time'].field ='v3d_time, scalar, series'

    fh.createVariable('salt_time','double',('salt_time'));
    fh.variables['salt_time'].long_name='salt_time'
    fh.variables['salt_time'].units =tunits
    fh.variables['salt_time'].field ='salt_time, scalar, series'

    fh.createVariable('temp_time','double',('temp_time'));
    fh.variables['temp_time'].long_name='temp_time'
    fh.variables['temp_time'].units =tunits
    fh.variables['temp_time'].field ='temp_time, scalar, series'
    
    fh.createVariable('lon_rho','float',('erho','xrho'));
    fh.variables['lon_rho'].long_name='lon_rho'
    fh.variables['lon_rho'].units ='degrees'
    fh.variables['lon_rho'].FillValue_ =100000.
    fh.variables['lon_rho'].missing_value =100000.
    fh.variables['lon_rho'].field ='xp, scalar, series'
    
    fh.createVariable('lat_rho','float',('erho','xrho'));
    fh.variables['lat_rho'].long_name='lat_rho'
    fh.variables['lat_rho'].units ='degrees'
    fh.variables['lat_rho'].FillValue_ =100000.
    fh.variables['lat_rho'].missing_value =100000.
    fh.variables['lat_rho'].field ='yp, scalar, series'
    
    fh.createVariable('zeta','float',('zeta_time','erho','xrho'));
    fh.variables['zeta'].long_name='free-surface'
    fh.variables['zeta'].units ='meter'
    fh.variables['zeta'].field ='free-surface, scalar, series'
    
    fh.createVariable('salt','float',('salt_time','s_rho','erho','xrho'));
    fh.variables['salt'].long_name='salinity'
    fh.variables['salt'].units ='PSU'
    fh.variables['salt'].field ='salinity, scalar, series'
    
    fh.createVariable('temp','float',('temp_time','s_rho','erho','xrho'));
    fh.variables['temp'].long_name='temperature'
    fh.variables['temp'].units ='C'
    fh.variables['temp'].field ='temperature, scalar, series'

    fh.createVariable('u','float',('v3d_time','s_rho','eu','xu'));
    fh.variables['u'].long_name='u-momentum component'
    fh.variables['u'].units ='meter second-1'
    fh.variables['u'].field ='u-velocity, scalar, series'
    
    fh.createVariable('v','float',('v3d_time','s_rho','ev','xv'));
    fh.variables['v'].long_name='v-momentum component'
    fh.variables['v'].units ='meter second-1'
    fh.variables['v'].field ='v-velocity, scalar, series'

    fh.createVariable('ubar','float',('v2d_time','eu','xu'));
    fh.variables['ubar'].long_name='vertically integrated u-momentum component'
    fh.variables['ubar'].units ='meter second-1'
    fh.variables['ubar'].field ='ubar-velocity, scalar, series'
    
    fh.createVariable('vbar','float',('v3d_time','ev','xv'));
    fh.variables['vbar'].long_name='vertically integrated v-momentum component'
    fh.variables['vbar'].units ='meter second-1'
    fh.variables['vbar'].field ='vbar-velocity, scalar, series'
    
    fh.close()

def create_clm_file(initfile,grd,tunits):
    """
    create_clm_file(initfile)

    Create a climatology file. 
    """

    print('Creating Climatology file')    
    # Create the initialization netCDF file with global attributes.

    if os.path.isfile(initfile):
        os.remove(initfile)
        print(f"{initfile} exists, Deleting.")
    else:
        print(f"{initfile} does not exist. Creating.")
            
    fh        = nc.Dataset(initfile, "w", format="NETCDF4")
    fh.type   =  'Climatology file for L1 grid'
    fh.title  =  'Climatology' 
    fh.history = ['Created by create_clm_file on '+ dt.datetime.now().strftime('%Y/%m/%d %H:%M:%S')]


    dims=grd.sizes
    atts=grd.attrs
    
    #Define the dimensions of the file
    fh.createDimension('xi_rho', dims['xi_rho'])
    fh.createDimension('xi_u', dims['xi_u'])
    fh.createDimension('xi_v', dims['xi_rho'])
    fh.createDimension('eta_rho', dims['eta_rho'])
    fh.createDimension('eta_u', dims['eta_rho'])
    fh.createDimension('eta_v', dims['eta_v'])
    fh.createDimension('s_rho', atts['sc_r'])
    fh.createDimension('s_w', atts['sc_w'])
    fh.createDimension('time',None)
    
    #Define the Variabes in the file
    fh.createVariable('lon_rho', 'f', ('eta_rho', 'xi_rho'))
    fh.variables['lon_rho'].long_name = 'longitude of RHO-points"'
    fh.variables['lon_rho'].units = 'degree_east'
    fh.variables['lon_rho'][:]=grd.lon_rho[:]
    
    fh.createVariable('lat_rho', 'f', ('eta_rho', 'xi_rho'))
    fh.variables['lat_rho'].long_name = 'latitude of RHO-points"'
    fh.variables['lat_rho'].units = 'degree_north'
    fh.variables['lat_rho'][:]=grd.lat_rho[:]
    
    fh.createVariable('lon_u', 'f', ('eta_u', 'xi_u'))
    fh.variables['lon_u'].long_name = 'longitude of U-points"'
    fh.variables['lon_u'].units = 'degree_east'
    fh.variables['lon_u'][:]=grd.lon_u[:]
    
    fh.createVariable('lat_u', 'f', ('eta_u', 'xi_u'))
    fh.variables['lat_u'].long_name = 'latitude of U-points"'
    fh.variables['lat_u'].units = 'degree_north'
    fh.variables['lat_u'][:]=grd.lat_u[:]
  
    fh.createVariable('lon_v', 'f', ('eta_v', 'xi_v'))
    fh.variables['lon_v'].long_name = 'longitude of V-points"'
    fh.variables['lon_v'].units = 'degree_east'
    fh.variables['lon_v'][:]=grd.lon_v[:]
    
    fh.createVariable('lat_v', 'f', ('eta_v', 'xi_v'))
    fh.variables['lat_v'].long_name = 'latitude of V-points"'
    fh.variables['lat_v'].units = 'degree_north'
    fh.variables['lat_v'][:]=grd.lat_v[:]
  
  
    fh.createVariable('spherical', 'short', ('time'))
    fh.variables['spherical'].long_name = 'grid type logical switch'
    fh.variables['spherical'].flag_meanings = 'Cartesian spherical"'
    fh.variables['spherical'].flag_values = '0, 1'
  
    fh.createVariable('Vtransform','long',('time'));
    fh.variables['Vtransform'].long_name ='vertical terrain-following transformation equation'
  
    fh.createVariable('Vstretching','long',('time'));
    fh.variables['Vstretching'].long_name ='vertical terrain-following stretching function'
  
    fh.createVariable('theta_b','double',('time'));
    fh.variables['theta_b'].long_name ='S-coordinate bottom control parameter'
    fh.variables['theta_b'].units ='1'
    
    fh.createVariable('theta_s','double',('time'));
    fh.variables['theta_s'].long_name ='S-coordinate surface control parameter'
    fh.variables['theta_s'].units ='1'
    
    fh.createVariable('Tcline','double',('time'));
    fh.variables['Tcline'].long_name ='S-coordinate surface/bottom layer width'
    fh.variables['Tcline'].units ='meter'
    
    fh.createVariable('hc','double',('time'));
    fh.variables['hc'].long_name ='S-coordinate parameter, critical depth'
    fh.variables['hc'].units ='meter'
    
    fh.createVariable('Cs_r','double',('s_rho'));
    fh.variables['Cs_r'].long_name='S-coordinate stretching curves at RHO-points'
    fh.variables['Cs_r'].valid_min =-1
    fh.variables['Cs_r'].valid_max =0
    
    fh.createVariable('Cs_w','double',('s_w'));
    fh.variables['Cs_w'].long_name='S-coordinate stretching curves at W-points'
    fh.variables['Cs_w'].valid_min =-1
    fh.variables['Cs_w'].valid_max =0
    
    fh.createVariable('s_rho','double',('s_rho'));
    fh.variables['s_rho'].long_name='S-coordinate at RHO-points'
    fh.variables['s_rho'].valid_min =-1
    fh.variables['s_rho'].valid_max =0
    fh.variables['s_rho'].positive ="up" 
    fh.variables['s_rho'].standard_name ="ocean_s_coordinate_g2"
    fh.variables['s_rho'].formula_terms="s: s_rho C: Cs_r eta: zeta depth: h depth_c: hc" 
            
    fh.createVariable('s_w','double',('s_w'));
    fh.variables['s_w'].long_name='S-coordinate at W-points'
    fh.variables['s_w'].valid_min =-1
    fh.variables['s_w'].valid_max =0
    fh.variables['s_w'].positive ="up"
    fh.variables['s_w'].standard_name = "ocean_s_coordinate_g2" 
    fh.variables['s_w'].formula_terms = "s: s_w C: Cs_w eta: zeta depth: h depth_c: hc"  
  
    fh.createVariable('ocean_time','double',('time'));
    fh.variables['ocean_time'].long_name='time since initialization'
    fh.variables['ocean_time'].units =tunits
    fh.variables['ocean_time'].field ='ocean_time, scalar, series'
    
    fh.createVariable('salt','f',('time','s_rho','eta_rho','xi_rho'));
    fh.variables['salt'].long_name='salinity'
    fh.variables['salt'].units ='PSU'
    fh.variables['salt'].time ="ocean_time" 
    fh.variables['salt'].coordinates ="lon_rho lat_rho ocean_time" 
    
    fh.createVariable('temp','f',('time','s_rho','eta_rho','xi_rho'));
    fh.variables['temp'].long_name='temperature'
    fh.variables['temp'].units ='C'
    fh.variables['temp'].time ="ocean_time" 
    fh.variables['temp'].coordinates ="lon_rho lat_rho ocean_time" 
  
    fh.createVariable('u','f',('time','s_rho','eta_u','xi_u'));
    fh.variables['u'].long_name='u-momentum component'
    fh.variables['u'].units ='meter second-1'
    fh.variables['u'].time ="ocean_time" 
    fh.variables['u'].coordinates ="lon_rho lat_rho ocean_time" 
  
    fh.createVariable('ubar','f',('time','eta_u','xi_u'));
    fh.variables['ubar'].long_name='vertically integrated u-momentum component'
    fh.variables['ubar'].units ='meter second-1'
    fh.variables['ubar'].time ="ocean_time" 
    fh.variables['ubar'].coordinates ="lon_rho lat_rho ocean_time" 
    
    fh.createVariable('v','f',('time','s_rho','eta_v','xi_v'));
    fh.variables['v'].long_name='v-momentum component'
    fh.variables['v'].units ='meter second-1'
    fh.variables['v'].time ="ocean_time" 
    fh.variables['v'].coordinates ="lon_rho lat_rho ocean_time" 
  
    fh.createVariable('vbar','f',('time','eta_v','xi_v'));
    fh.variables['vbar'].long_name='vertically integrated v-momentum component'
    fh.variables['vbar'].units ='meter second-1'
    fh.variables['vbar'].time = "ocean_time" 
    fh.variables['vbar'].coordinates = "lon_rho lat_rho ocean_time" 
    
    fh.createVariable('zeta','f',('time','eta_rho','xi_rho'));
    fh.variables['zeta'].long_name='free-surface'
    fh.variables['zeta'].units ='meter'
    fh.variables['zeta'].time = "ocean_time" 
    fh.variables['zeta'].coordinates = "lon_rho lat_rho ocean_time"  
    
    
        
    fh.close()

def rangebearing(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6378.13 # Radius of earth in kilometers. Use 3956 for miles. Determines return value units.
    dist=c*r
    dLon = lon2 - lon1
    y = sin(dLon) * cos(lat2)
    x = cos(lat1) * sin(lat2) \
        - sin(lat1) * cos(lat2) * cos(dLon)
    bear=atan2(y, x)
     
    return dist,bear



# def uv_rot(u, v, grid,angle,hboundary="extend", hfill_value=None):
#     """Calculate East/Norths Components from u and v components

#     Inputs
#     ------
#     u: DataArray
#         xi component of velocity [m/s]
#     v: DataArray
#         eta component of velocity [m/s]
#     grid: xgcm.grid
#         Grid object associated with u, v
#     hboundary: string, optional
#         Passed to `grid` method calls; horizontal boundary selection
#         for moving to rho grid.
#         From xgcm documentation:
#         A flag indicating how to handle boundaries:
#         * None:  Do not apply any boundary conditions. Raise an error if
#           boundary conditions are required for the operation.
#         * 'fill':  Set values outside the array boundary to fill_value
#           (i.e. a Neumann boundary condition.)
#         * 'extend': Set values outside the array to the nearest array
#           value. (i.e. a limited form of Dirichlet boundary condition.
#     hfill_value: float, optional
#         Passed to `grid` method calls; horizontal boundary fill value
#         selection for moving to rho grid.
#         From xgcm documentation:
#         The value to use in the boundary condition with `boundary='fill'`.

#     Returns
#     -------
#     DataArray of speed calculated on rho/rho grids.
#     Output is `[T,Z,Y,X]`.

#     Notes
#     -----
    

#     Example usage
#     -------------
#     >>> xroms.speed(ds.u, ds.v, grid)
#     """

#     assert isinstance(u, xr.DataArray), "var must be DataArray"
#     assert isinstance(v, xr.DataArray), "var must be DataArray"

#     u = xroms.to_rho(u, grid, hboundary=hboundary, hfill_value=hfill_value)
#     v = xroms.to_rho(v, grid, hboundary=hboundary, hfill_value=hfill_value)

#     u=u.rename({"z_rho_u":"z_rho"})
#     v=v.rename({"z_rho_v":"z_rho"})

#     W=(u+v*1j)*np.exp(1j*angle)
#     var1 = np.real(W)
#     var2 = np.imag(W)
    
#     var1.attrs["name"] = "u_eastward"
#     var1.attrs["long_name"] = "Easterly velocity"
#     var1.attrs["units"] = "m/s"
#     var1.attrs["grid"] = grid
#     var1.name = var1.attrs["name"]

#     var2.attrs["name"] = "v_northward"
#     var2.attrs["long_name"] = "Northerly velocity"
#     var2.attrs["units"] = "m/s"
#     var2.attrs["grid"] = grid
#     var2.name = var2.attrs["name"]
    
#     return var1,var2

def uv_rot_2d(u, v, grid,angle,hboundary="extend", hfill_value=None,reverse=False):
    """Calculate East/Norths Components from u and v components

    Inputs
    ------
    u: DataArray
        xi component of velocity [m/s]
    v: DataArray
        eta component of velocity [m/s]
    grid: xgcm.grid
        Grid object associated with u, v
    hboundary: string, optional
        Passed to `grid` method calls; horizontal boundary selection
        for moving to rho grid.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    hfill_value: float, optional
        Passed to `grid` method calls; horizontal boundary fill value
        selection for moving to rho grid.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.

    Returns
    -------
    DataArray of speed calculated on rho/rho grids.
    Output is `[T,Z,Y,X]`.


    """

    assert isinstance(u, xr.DataArray), "var must be DataArray"
    assert isinstance(v, xr.DataArray), "var must be DataArray"
  #  print('ROTATING')
    if reverse:
     #   print('Geographic TO roms')
        W=(u+v*1j)*np.exp(1j*angle*-1.0)
        u = np.real(W)
        v = np.imag(W)
        var1 = xroms.to_u(u, grid, hboundary=hboundary, hfill_value=hfill_value)
        var2 = xroms.to_v(v, grid, hboundary=hboundary, hfill_value=hfill_value)
    
        var1.attrs["name"] = "ubar"
        var1.attrs["long_name"] = "vertically integrated u-momentum component"
        var1.attrs["units"] = "meter second-1"
        var1.attrs["grid"] = grid
        var1.attrs["location"] = 'edge1'
        var1.attrs["field"] = 'ubar-velocity, scalar, series'
        var1.name = var1.attrs["name"]

        var2.attrs["name"] = "vbar"
        var2.attrs["long_name"] = "vertically integrated v-momentum component"
        var2.attrs["units"] = "meter second-1"
        var2.attrs["grid"] = grid
        var1.attrs["location"] = 'edge2'
        var1.attrs["field"] = 'vbar-velocity, scalar, series'
        var2.name = var2.attrs["name"]
        
    else:
     #   print('ROMS TO Geographic')
        u = xroms.to_rho(u, grid, hboundary=hboundary, hfill_value=hfill_value)
        v = xroms.to_rho(v, grid, hboundary=hboundary, hfill_value=hfill_value)

        W=(u+v*1j)*np.exp(1j*angle)
        var1 = np.real(W)
        var2 = np.imag(W)
    
        var1.attrs["name"] = "ubar_eastward"
        var1.attrs["long_name"] = "Easterly velocity"
        var1.attrs["units"] = "m/s"
        var1.attrs["grid"] = grid
        var1.name = var1.attrs["name"]

        var2.attrs["name"] = "vbar_northward"
        var2.attrs["long_name"] = "Northerly velocity"
        var2.attrs["units"] = "m/s"
        var2.attrs["grid"] = grid
        var2.name = var2.attrs["name"]
    
    return var1,var2

def uv_rot(u, v, grid,angle,hboundary="extend", hfill_value=None,reverse=False):
    """Calculate East/Norths Components from u and v components

    Inputs
    ------
    u: DataArray
        xi component of velocity [m/s]
    v: DataArray
        eta component of velocity [m/s]
    grid: xgcm.grid
        Grid object associated with u, v
    hboundary: string, optional
        Passed to `grid` method calls; horizontal boundary selection
        for moving to rho grid.
        From xgcm documentation:
        A flag indicating how to handle boundaries:
        * None:  Do not apply any boundary conditions. Raise an error if
          boundary conditions are required for the operation.
        * 'fill':  Set values outside the array boundary to fill_value
          (i.e. a Neumann boundary condition.)
        * 'extend': Set values outside the array to the nearest array
          value. (i.e. a limited form of Dirichlet boundary condition.
    hfill_value: float, optional
        Passed to `grid` method calls; horizontal boundary fill value
        selection for moving to rho grid.
        From xgcm documentation:
        The value to use in the boundary condition with `boundary='fill'`.
    Reverse: Logical, optional

    Returns
    -------
    DataArray of speed calculated on rho/rho grids.
    Output is `[T,Z,Y,X]`.

    """

    assert isinstance(u, xr.DataArray), "var must be DataArray"
    assert isinstance(v, xr.DataArray), "var must be DataArray"

    if reverse:
      #  print('Geographic TO roms')
  
        W=(u+v*1j)*np.exp(1j*angle*-1.0)
        u = np.real(W)
        v = np.imag(W)
        var1 = xroms.to_u(u, grid, hboundary=hboundary, hfill_value=hfill_value)
        var2 = xroms.to_v(v, grid, hboundary=hboundary, hfill_value=hfill_value)
        

              
              
        var1.attrs["name"] = "u"
        var1.attrs["long_name"] = "u-momentum component"
        var1.attrs["units"] = "meter second-1"
        var1.attrs["grid"] = grid
        var1.attrs["location"] = 'edge1'
        var1.attrs["field"] = 'u-velocity, scalar, series'
        var1.name = var1.attrs["name"]

        var2.attrs["name"] = "v"
        var2.attrs["long_name"] = "v-momentum component"
        var2.attrs["units"] = "meter second-1"
        var2.attrs["grid"] = grid
        var2.attrs["location"] = 'edge2'
        var2.attrs["field"] = 'v-velocity, scalar, series'
        var2.name = var2.attrs["name"]
    else:
     #   print('ROMS TO Geographic')
        u = xroms.to_rho(u, grid, hboundary=hboundary, hfill_value=hfill_value)
        v = xroms.to_rho(v, grid, hboundary=hboundary, hfill_value=hfill_value)
        W=(u+v*1j)*np.exp(1j*angle)
        var1 = np.real(W)
        var2 = np.imag(W)
    
        var1.attrs["name"] = "u_eastward"
        var1.attrs["long_name"] = "Easterly velocity"
        var1.attrs["units"] = "m/s"
        var1.attrs["grid"] = grid
        var1.name = var1.attrs["name"]

        var2.attrs["name"] = "v_northward"
        var2.attrs["long_name"] = "Northerly velocity"
        var2.attrs["units"] = "m/s"
        var2.attrs["grid"] = grid
        var2.name = var2.attrs["name"]
    


    return var1,var2


