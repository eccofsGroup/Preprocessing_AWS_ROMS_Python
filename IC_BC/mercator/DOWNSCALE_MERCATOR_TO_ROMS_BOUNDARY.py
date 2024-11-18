"""
DOWNSCALE_MERCATOR_TO ROMS.py

Adapted from DOWNSCALE_ROMS_OUTPUT.py
Created by Elias Hunter, hunter@marine.rutgers.edu, 12/19/2023 
"""
import os,glob
import numpy as np
import xroms
import xarray as xr
import pandas as pd
import xesmf as xe
import cartopy 
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import datetime,time
from datetime import timedelta,date
import romsutil as rutil
from  scipy.interpolate import interp1d
import netCDF4 as nc


import warnings
warnings.filterwarnings("ignore")
#Set relevant inout parameters

proj = cartopy.crs.Mercator(central_longitude=-74)
pc = cartopy.crs.PlateCarree()


latmin=-3.0
latmax=53.0
lonmin=-100.0
lonmax=-38.0

#Set time



reftime=datetime.datetime(2011,1,1)
tunits=reftime.strftime('days since %Y-%m-%d')
rtime=np.datetime64(reftime.date())

#Donor grid Info


#Donor inputfiles
datadir='/home/om/cron/ECCOFS_OBS/MERCATOR/data/raw/'
regrid_coef_file='/home/om/cron/ECCOFS_OBS/MERCATOR/data/mercator_bdry_9999.nc'


#Receiver Grid Info
L1grdfile='/home/om/roms/eccofs/grid_eccofs_6km_02.nc' # can be a thredds url
L1theta_s=7.0
L1theta_b=2.0
L1Tcline=250.0
L1N=50
L1Vtransform  =2        #vertical transformation equation
L1Vstretching =4        #vertical stretching function
L1hc=250
# Nbed=1
# Nveg=1
# NCS=0
# NNS=1

(s_w_L1,C_w_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,1)
(s_r_L1,C_r_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,0)




#Output fflags

bdylist={'EAST':True, 'WEST':False, 'NORTH':False, 'SOUTH':True}
#bdylist={'SOUTH':True}


nday=7
start_date = datetime.datetime(2019, 1, 1)
sdate=start_date.strftime('%Y%m%d')
date_times = [start_date  + timedelta(days=i) for i in range(nday)] 

#Receiver Output file
bryfile=f'/home/om/cron/ECCOFS_OBS/MERCATOR/data/bdry/ECCOFS_BDRY_{sdate}_bry.nc'

############################################################################
#Main Program
############################################################################

def main():
    ########################################################################
    print('Initilaizing grids and regridder')
    ########################################################################
    
    #Getting the mercator grid information
    day=date_times[0]
 #   mfilename = day.strftime(datadir+'mercator_doppio_%Y_%m_%d.nc')
 #   dsmgrd=xr.open_dataset(mfilename)
 #   mlat=dsmgrd['latitude']
 #   mlon=dsmgrd['longitude']
 #   mdepth=dsmgrd['depth']
    
    cfgrd=xroms.open_netcdf(L1grdfile)
    cfgrd.attrs['sc_r']=L1N
    cfgrd.attrs['sc_w']=L1N+1
    #plot_points_LL
    filelist=[]
    for day in date_times:
        
        mfilename = day.strftime(datadir+'mercator_doppio_%Y_%m_%d.nc')
        
        filelist.append(mfilename)

    print(filelist)   
    dsmerc=xr.open_mfdataset(filelist)
   # blon,blat = np.meshgrid(dsmerc['longitude'].values,dsmerc['latitude'].values)
    
  #  print(blat)
   # plot_points_LL(cfgrd['lat_rho'],cfgrd['lon_rho'],dsmerc)
    
    ########################################################################
    #Process donwscaling files
    ########################################################################



    start_time=time.time()
    downscale_bdry_file(cfgrd,dsmerc)
    end_time=time.time()
    elapsed_time = end_time - start_time
    print(f"boundary conditions  file processing time: {elapsed_time} seconds")
        
def downscale_bdry_file(cfgrd,dsmerc):
    
    
    print(f'Procesing Boundary Condition file: {bryfile}')
    rutil.create_bdry_file(bryfile,cfgrd,tunits)

#     tmp=dsmerc.rename({'mask_rho':'mask'})
 #   dsmerc=tmp.rename({'lat_rho':'lat'})
 #   dsmerc=tmp.rename({'lon_rho':'lon'})
# ##############################################################################################
# #Initialize regridder. 
# ##############################################################################################
    varnew_south=xroms.subset(cfgrd,Y=slice(0,3))
    varnew_south=varnew_south.rename({'lat_rho':'lat'})
    varnew_south=varnew_south.rename({'lon_rho':'lon'})
    #varnew_south=varnew_south.rename({'mask_rho':'mask'})

    varnew_north=xroms.subset(cfgrd,Y=slice(cfgrd.sizes['eta_rho']-3,cfgrd.sizes['eta_rho']))
    varnew_north=varnew_north.rename({'lat_rho':'lat'})
    varnew_north=varnew_north.rename({'lon_rho':'lon'})
    #varnew_north=varnew_north.rename({'mask_rho':'mask'})


    varnew_west=xroms.subset(cfgrd,X=slice(0,3))
    varnew_west=varnew_west.rename({'lat_rho':'lat'})
    varnew_west=varnew_west.rename({'lon_rho':'lon'})
    #varnew_west=varnew_west.rename({'mask_rho':'mask'})


    varnew_east=xroms.subset(cfgrd,X=slice(cfgrd.sizes['xi_rho']-3,cfgrd.sizes['xi_rho']))
    varnew_east=varnew_east.rename({'lat_rho':'lat'})
    varnew_east=varnew_east.rename({'lon_rho':'lon'})
    #varnew_east=varnew_east.rename({'mask_rho':'mask'})
    start_time=time.time()
    
    regrid_coef_file_tmp=regrid_coef_file.replace('9999','EAST')
    if os.path.exists(regrid_coef_file_tmp):
        print('Re-using weights')
        regridder_E = xe.Regridder( dsmerc,varnew_east,"bilinear",extrap_method="nearest_s2d",locstream_out=False,filename=regrid_coef_file_tmp,reuse_weights=True)  
    else:
        print('Generating weights')
        regridder_E = xe.Regridder(dsmerc,varnew_east, "bilinear",extrap_method="nearest_s2d",locstream_out=False)
        regridder_E.to_netcdf(regrid_coef_file_tmp)

    
    regrid_coef_file_tmp=regrid_coef_file.replace('9999','WEST')
    if os.path.exists(regrid_coef_file_tmp):
        print('Re-using weights')
        regridder_W = xe.Regridder( dsmerc,varnew_west,"bilinear",extrap_method="nearest_s2d",locstream_out=False,filename=regrid_coef_file_tmp,reuse_weights=True)  
    else:
        print('Generating weights')
        regridder_W = xe.Regridder(dsmerc,varnew_west, "bilinear",extrap_method="nearest_s2d",locstream_out=False)
        regridder_W.to_netcdf(regrid_coef_file_tmp)
            
    regrid_coef_file_tmp=regrid_coef_file.replace('9999','NORTH')
    if os.path.exists(regrid_coef_file_tmp):
        print('Re-using weights')
        regridder_N = xe.Regridder( dsmerc,varnew_north,"bilinear",extrap_method="nearest_s2d",locstream_out=False,filename=regrid_coef_file_tmp,reuse_weights=True)  
    else:
        print('Generating weights')
        regridder_N = xe.Regridder(dsmerc,varnew_north, "bilinear",extrap_method="nearest_s2d",locstream_out=False)
        regridder_N.to_netcdf(regrid_coef_file_tmp)
            
    regrid_coef_file_tmp=regrid_coef_file.replace('9999','SOUTH')
    if os.path.exists(regrid_coef_file_tmp):
        print('Re-using weights')
        regridder_S = xe.Regridder( dsmerc,varnew_south,"bilinear",extrap_method="nearest_s2d",locstream_out=False,filename=regrid_coef_file_tmp,reuse_weights=True)  
    else:
        print('Generating weights')
        regridder_S = xe.Regridder(dsmerc,varnew_south, "bilinear",extrap_method="nearest_s2d",locstream_out=False)
        regridder_S.to_netcdf(regrid_coef_file_tmp)
            
            
            
            
    end_time=time.time()
    elapsed_time = end_time - start_time
    print(f'Create regridder: {elapsed_time}')
    
    
# ##############################################################################################
# # 3-D variable interpolation
# ##############################################################################################
    for ind,t in enumerate(dsmerc.time.values):
        print(f'Procesing Boundary Condition time: {t}')
        dsmerc_I=dsmerc.sel(time=t)

        
# ##############################################################################################
# #Regridding
# ##############################################################################################


        for key in bdylist:
            if bdylist[key]:
                print(f'Processing {key} boundary')
            else:
                print(f'Skipping {key} boundary')
                continue
            
            start_time=time.time()    
            if key=="EAST":
                dataout = regridder_E(dsmerc_I,keep_attrs=True)
                dataout['xi_u']=varnew_east['xi_u']
                dataout['eta_v']=varnew_east['eta_v']
                varnew=varnew_east
            if key=="WEST":
                dataout = regridder_W(dsmerc_I,keep_attrs=True)
                dataout['xi_u']=varnew_west['xi_u']
                dataout['eta_v']=varnew_west['eta_v'] 
                varnew=varnew_west
            if key=="SOUTH":
                dataout = regridder_S(dsmerc_I,keep_attrs=True)
                dataout['xi_u']=varnew_south['xi_u']
                dataout['eta_v']=varnew_south['eta_v'] 
                varnew=varnew_south
            if key=="NORTH":
                dataout = regridder_N(dsmerc_I,keep_attrs=True)
                dataout['xi_u']=varnew_north['xi_u']
                dataout['eta_v']=varnew_north['eta_v'] 
                varnew=varnew_north

            
            end_time=time.time()
            elapsed_time = end_time - start_time
            print(f'Regridding {key} boundary: {elapsed_time}')
    



        
# ##############################################################################################
# #Get depths of L0 and L1 grid. 
# ##############################################################################################
        
            dataout_z=dataset_get_Z(dataout,varnew)
    
    
            dim_dict=varnew.dims
##############################################################################################
#Initialize Arrays
##############################################################################################
            print('INITIALZING OUTPUT ARRAYS')
            if key=='EAST' :
                temp=np.full((L1N,dim_dict['eta_rho'],3),np.nan)    
                salt=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
                u_east=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
                v_north=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
            if key=='WEST' :
                temp=np.full((L1N,dim_dict['eta_rho'],3),np.nan)    
                salt=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
                u_east=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
                v_north=np.full((L1N,dim_dict['eta_rho'],3),np.nan)   
            if key=='NORTH' :
                temp=np.full((L1N,3,dim_dict['xi_rho']),np.nan)    
                salt=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
                u_east=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
                v_north=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
            if key=='SOUTH' :
                temp=np.full((L1N,3,dim_dict['xi_rho']),np.nan)    
                salt=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
                u_east=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
                v_north=np.full((L1N,3,dim_dict['xi_rho']),np.nan)   
           

# # tlat=[]
# # tlon=[]
# ##############################################################################################
# #Loading Data for interpolation
# ##############################################################################################
            print('LOADING DATA FOR VERTICAL INTERP.')
            mask=varnew.mask_rho.values
            tmpz=dataout_z.z_rho.values
            tmpz2=dataout['depth']
            temp_I=dataout.thetao.values
            salt_I=dataout.so.values
            u_east_I=dataout.uo.values
            v_north_I=dataout.vo.values

        
#  #       start_time = time.time()
# ##############################################################################################
# # Vertical Interpolation
# ##############################################################################################
          
            start_time=time.time()

            for eta in range(0,dim_dict['eta_rho']):
                for xi in range(0,dim_dict['xi_rho']):
                    maskflag=mask[eta,xi]
                    if maskflag==0.0:
                        continue
                    
        
                #    print('----------------------')
                           
                    all_nan = np.all(np.isnan(temp_I[:,eta,xi]))
                    if ~all_nan:
                     #   print('TEMP NAN')
                        
                        z=-tmpz2.values     
                        data=temp_I[:,eta,xi]
                        Igood=np.where(~np.isnan(data))
                        ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                        ntmp=ifun(tmpz[:,eta,xi])
                        temp[:,eta,xi]=ntmp
                        
                        
                        
                    all_nan = np.all(np.isnan(salt_I[:,eta,xi]))
                    if ~all_nan:
                    #    print('Salt NAN')
                        z=-tmpz2.values
                        data=salt_I[:,eta,xi]
                        Igood=np.where(~np.isnan(data))
                        ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                        ntmp=ifun(tmpz[:,eta,xi])
                        salt[:,eta,xi]=ntmp
                        
        
                    
                    all_nan = np.all(np.isnan(u_east_I[:,eta,xi]))
                    if ~all_nan:
                    #    print('U NAN')
                        z=-tmpz2.values
                        data=u_east_I[:,eta,xi]
                        Igood=np.where(~np.isnan(data))
                        ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                        ntmp=ifun(tmpz[:,eta,xi])
                        u_east[:,eta,xi]=ntmp
                    
                    all_nan = np.all(np.isnan(v_north_I[:,eta,xi]))
                    if ~all_nan:
                       # print('V NAN')
                        z=-tmpz2.values
                        data=v_north_I[:,eta,xi]
                        Igood=np.where(~np.isnan(data))
                        ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                        ntmp=ifun(tmpz[:,eta,xi])
                        v_north[:,eta,xi]=ntmp        
        
        
                    
            end_time=time.time()
            elapsed_time = end_time - start_time
            print(f"processing time: {elapsed_time} seconds") 
         #   eta_rho=dataout['eta_rho'];
         #   xi_rho=dataout['xi_rho'];
    
     
            for eta in range(0,dim_dict['eta_rho']):
                for xi in range(0,dim_dict['xi_rho']):
                    maskflag=mask[eta,xi]
                    if maskflag==0.0:
                        continue

                    all_nan = np.all(np.isnan(temp[:,eta,xi]))
                    
                    if all_nan:
                      
                        tmp=np.squeeze(temp[0,:,:])
                        real_value_indices = np.argwhere(~np.isnan(tmp))
                        
                        distances = np.sqrt((real_value_indices[:, 0] - eta)**2 + 
                            (real_value_indices[:, 1] - xi)**2)
                        cI = real_value_indices[np.argmin(distances)]
                        temp[:,eta,xi]=temp[:,cI[0],cI[1]]
               
                        
                        
                    all_nan = np.all(np.isnan(salt[:,eta,xi]))
                    if all_nan:
                      
                        tmp=np.squeeze(salt[0,:,:])
                        real_value_indices = np.argwhere(~np.isnan(tmp))
                        distances = np.sqrt((real_value_indices[:, 0] - eta)**2 + 
                            (real_value_indices[:, 1] - xi)**2)
                        cI = real_value_indices[np.argmin(distances)]
                        salt[:,eta,xi]=salt[:,cI[0],cI[1]]
                      
                    all_nan = np.all(np.isnan(u_east[:,eta,xi]))
                    if all_nan:
                      
                        tmp=np.squeeze(u_east[0,:,:])
                        real_value_indices = np.argwhere(~np.isnan(tmp))
                        distances = np.sqrt((real_value_indices[:, 0] - eta)**2 + 
                            (real_value_indices[:, 1] - xi)**2)
                        cI = real_value_indices[np.argmin(distances)]
                        u_east[:,eta,xi]=u_east[:,cI[0],cI[1]]
                    
                        
                    all_nan = np.all(np.isnan(v_north[:,eta,xi]))
                    if all_nan:
                      
                        tmp=np.squeeze(v_north[0,:,:])
                        real_value_indices = np.argwhere(~np.isnan(tmp))
                        distances = np.sqrt((real_value_indices[:, 0] - eta)**2 + 
                            (real_value_indices[:, 1] - xi)**2)
                        cI = real_value_indices[np.argmin(distances)]
                        v_north[:,eta,xi]=v_north[:,cI[0],cI[1]]
                 
                
                    
            end_time=time.time()
            elapsed_time = end_time - start_time
            print(f"processing time: {elapsed_time} seconds")


#     ##############################################################################################
#     #Write Output to netcdf file
#     ##############################################################################################
            dataout_z['temp']=(('s_rho', 'eta_rho', 'xi_rho'),temp)
            dataout_z['salt']=(('s_rho', 'eta_rho', 'xi_rho'),salt) 
            dataout_z['u_eastward']=(('s_rho', 'eta_rho', 'xi_rho'),u_east) 
            dataout_z['v_northward']=(('s_rho', 'eta_rho', 'xi_rho'),v_north) 
            dataout_z['zeta']=dataout_z['zeta'].fillna(0)
            (dataout_z,gridout)=xroms.roms_dataset(dataout_z,Vtransform=L1Vtransform)
            uv=rutil.uv_rot(dataout_z.u_eastward, dataout_z.v_northward, gridout,cfgrd.angle,reverse=True)
            ruhis=uv[0]
            rvhis=uv[1]
            dataout_z=xr.merge([dataout_z,ruhis,rvhis])
            dataout_z.load()
            ubar=dataout_z.u.xroms.gridmean(gridout,"Z")
            vbar=dataout_z.v.xroms.gridmean(gridout,"Z")
            
            
            timeout=(t-np.datetime64(rtime)) / np.timedelta64(1, 'D')
            
            
            
            ncid = nc.Dataset(bryfile, "r+", format="NETCDF4")
            ncid.variables['bry_time'][ind]=timeout
            print(temp.shape)
            if key=='EAST':
                ncid.variables['u_east'][ind,:,:]=dataout_z.u.values[:,:,-1]
                ncid.variables['v_east'][ind,:,:]=dataout_z.v.values[:,:,-1]
                ncid.variables['ubar_east'][ind,:]=ubar.values[:,-1]
                ncid.variables['vbar_east'][ind,:]=vbar.values[:,-1]
                ncid.variables['salt_east'][ind,:,:]=salt[:,:,-1]
                ncid.variables['temp_east'][ind,:,:]=temp[:,:,-1]
                ncid.variables['zeta_east'][ind,:]=dataout_z['zeta'].values[:,-1]
            if key=='WEST':
                ncid.variables['u_west'][ind,:,:]=dataout_z.u.values[:,:,0]
                ncid.variables['v_west'][ind,:,:]=dataout_z.v.values[:,:,0]
                ncid.variables['ubar_west'][ind,:]=ubar.values[:,0]
                ncid.variables['vbar_west'][ind,:]=vbar.values[:,0]
                ncid.variables['salt_west'][ind,:,:]=salt[:,:,0]
                ncid.variables['temp_west'][ind,:,:]=temp[:,:,0]
                ncid.variables['zeta_west'][ind,:]=dataout_z['zeta'].values[:,0]
            if key=='NORTH':
                ncid.variables['u_north'][ind,:,:]=dataout_z.u.values[:,-1,:]
                ncid.variables['v_north'][ind,:,:]=dataout_z.v.values[:,-1,:]
                ncid.variables['ubar_north'][ind,:]=ubar.values[-1,:]
                ncid.variables['vbar_north'][ind,:]=vbar.values[-1,:]
                ncid.variables['salt_north'][ind,:,:]=salt[:,-1,:]
                ncid.variables['temp_north'][ind,:,:]=temp[:,-1,:]
                ncid.variables['zeta_north'][ind,:]=dataout_z['zeta'].values[-1,:]
            
            if key=='SOUTH':
                ncid.variables['u_south'][ind,:,:]=dataout_z.u.values[:,0,:]
                ncid.variables['v_south'][ind,:,:]=dataout_z.v.values[:,0,:]
                ncid.variables['ubar_south'][ind,:]=ubar.values[0,:]
                ncid.variables['vbar_south'][ind,:]=vbar.values[0,:]
                ncid.variables['salt_south'][ind,:,:]=salt[:,0,:]
                ncid.variables['temp_south'][ind,:,:]=temp[:,0,:]
                ncid.variables['zeta_south'][ind,:]=dataout_z['zeta'].values[0,:]
                
                
                
#         # ncid.variables['ubar_north'][ind,:]=xromhisL1_N_z.ubar.values[:,-1]
#         # ncid.variables['vbar_north'][ind,:]=xromhisL1_N_z.vbar.values[:,-1]
#         # ncid.variables['zeta_north'][ind,:]=xromhisL1_N_z.zeta.values[-1,:]
#         ncid.variables['u_north'][ind,:,:]=xromhisL1_N_z.u.values[:,:,-1]
#         ncid.variables['v_north'][ind,:,:]=xromhisL1_N_z.v.values[:,:,-1]
#         ncid.variables['salt_north'][ind,:,:]=xromhisL1_N_z.salt.values[:,:]
#         ncid.variables['temp_north'][ind,:,:]=xromhisL1_N_z.temp.values[:,:]

#         # ncid.variables['ubar_south'][ind,:]=xromhisL1_S_z.ubar.values[:,0]
#         # ncid.variables['vbar_south'][ind,:]=xromhisL1_S_z.vbar.values[:,0]
#         # ncid.variables['zeta_south'][ind,:]=xromhisL1_S_z.zeta.values[0,:]
#         ncid.variables['u_south'][ind,:,:]=xromhisL1_S_z.u.values[:,:,0]
#         ncid.variables['v_south'][ind,:,:]=xromhisL1_S_z.v.values[:,:,0]
#         ncid.variables['salt_south'][ind,:,:]=xromhisL1_S_z.salt.values[:,:]
#         ncid.variables['temp_south'][ind,:,:]=xromhisL1_S_z.temp.values[:,:]


#         # ncid.variables['ubar_west'][ind,:]=xromhisL1_W_z.ubar.values[:,0]
#         # ncid.variables['vbar_west'][ind,:]=xromhisL1_W_z.vbar.values[:,0]
#         # ncid.variables['zeta_west'][ind,:]=xromhisL1_W_z.zeta.values[:,0]
#         ncid.variables['u_west'][ind,:,:]=xromhisL1_W_z.u.values[:,:,0]
#         ncid.variables['v_west'][ind,:,:]=xromhisL1_W_z.v.values[:,:,0]
#         ncid.variables['salt_west'][ind,:,:]=xromhisL1_W_z.salt.values[:,:]
#         ncid.variables['temp_west'][ind,:,:]=xromhisL1_W_z.temp.values[:,:]


#         # ncid.variables['ubar_east'][ind,:]=xromhisL1_E_z.ubar.values[:,-1]
#         # ncid.variables['vbar_east'][ind,:]=xromhisL1_E_z.vbar.values[:,-1]
#         # ncid.variables['zeta_east'][ind,:]=xromhisL1_E_z.zeta.values[:,-1]
#         ncid.variables['u_east'][ind,:,:]=xromhisL1_E_z.u.values[:,:,-1]
#         ncid.variables['v_east'][ind,:,:]=xromhisL1_E_z.v.values[:,:,-1]
#         ncid.variables['salt_east'][ind,:,:]=xromhisL1_E_z.salt.values[:,:]
#         ncid.variables['temp_east'][ind,:,:]=xromhisL1_E_z.temp.values[:,:]


            ncid.sync()
            ncid.close()
        
#    ##############################################################################################
#     #Process barotropic data
#     ##############################################################################################
#     dsqckl0=rechunk_eta_xi(dsqckl0)


#     for ind,t in enumerate(dsqckl0.ocean_time.values[0:t2end]):
#         print(t)
#         dsqckl0_I=dsqckl0.sel(ocean_time=t)

#         (xromqckL0,gridqckL0)=xroms.roms_dataset(dsqckl0_I,Vtransform=L0Vtransform)
#         uv=rutil.uv_rot_2d(xromqckL0.ubar, xromqckL0.vbar, gridqckL0,xromqckL0.angle)
#         ru=uv[0]
#         rv=uv[1]
#         xromqckL0=xr.merge([xromqckL0,ru,rv])
#         xromqckL0.load()
#     ##############################################################################################
#     #Regridding
#     # ##############################################################################################
#         xromqckL1_N = regridder_N(xromqckL0,keep_attrs=True)
#         xromqckL1_N=xromqckL1_N.drop_vars(['lon_psi','lat_psi'])
#         xromqckL1_N['xi_u']=varnew_north['xi_u']
#         xromqckL1_N['eta_v']=varnew_north['eta_v']
#         xromqckL1_N=rechunk_eta_xi(xromqckL1_N)
#         (xromqckL1_N,gridqckL1_N)=xroms.roms_dataset(xromqckL1_N,Vtransform=L1Vtransform)
  
        
#         xromqckL1_S = regridder_S(xromqckL0,keep_attrs=True)
#         xromqckL1_S=xromqckL1_S.drop_vars(['lon_psi','lat_psi'])
#         xromqckL1_S['xi_u']=varnew_south['xi_u']
#         xromqckL1_S['eta_v']=varnew_south['eta_v']
#         xromqckL1_S=rechunk_eta_xi(xromqckL1_S)
#         (xromqckL1_S,gridqckL1_S)=xroms.roms_dataset(xromqckL1_S,Vtransform=L1Vtransform)
    
#         xromqckL1_E = regridder_E(xromqckL0,keep_attrs=True)
#         xromqckL1_E=xromqckL1_E.drop_vars(['lon_psi','lat_psi'])
#         xromqckL1_E['xi_u']=varnew_east['xi_u']
#         xromqckL1_E['eta_v']=varnew_east['eta_v']
#         xromqckL1_E=rechunk_eta_xi(xromqckL1_E)
#         (xromqckL1_E,gridqckL1_E)=xroms.roms_dataset(xromqckL1_E,Vtransform=L1Vtransform)
    
#         xromqckL1_W = regridder_W(xromqckL0,keep_attrs=True)
#         xromqckL1_W=xromqckL1_W.drop_vars(['lon_psi','lat_psi'])
#         xromqckL1_W['xi_u']=varnew_west['xi_u']
#         xromqckL1_W['eta_v']=varnew_west['eta_v']
#         xromqckL1_W=rechunk_eta_xi(xromqckL1_W)
#         (xromqckL1_W,gridqckL1_W)=xroms.roms_dataset(xromqckL1_W,Vtransform=L1Vtransform)


#     ##############################################################################################
#     #Rotate 
#     ##############################################################################################
      
#         uv=rutil.uv_rot_2d(xromqckL1_N.ubar_eastward, xromqckL1_N.vbar_northward, gridqckL1_N,xromqckL1_N.angle,reverse=True)
#         ru=uv[0]
#         rv=uv[1]
#         xromqckL1_N=xr.merge([xromqckL1_N,ru,rv])

#         uv=rutil.uv_rot_2d(xromqckL1_S.ubar_eastward, xromqckL1_S.vbar_northward, gridqckL1_S,xromqckL1_S.angle,reverse=True)
#         ru=uv[0]
#         rv=uv[1]
#         xromqckL1_S=xr.merge([xromqckL1_S,ru,rv])
        
#         uv=rutil.uv_rot_2d(xromqckL1_E.ubar_eastward, xromqckL1_E.vbar_northward, gridqckL1_E,xromqckL1_E.angle,reverse=True)
#         ru=uv[0]
#         rv=uv[1]
#         xromqckL1_E=xr.merge([xromqckL1_E,ru,rv])
        
#         uv=rutil.uv_rot_2d(xromqckL1_W.ubar_eastward, xromqckL1_W.vbar_northward, gridqckL1_W,xromqckL1_W.angle,reverse=True)
#         ru=uv[0]
#         rv=uv[1]
#         xromqckL1_W=xr.merge([xromqckL1_W,ru,rv])

        
 
#     # ##############################################################################################
#     # #Write Output to netcdf file
#     # ##############################################################################################
        
#         timeout=(t-np.datetime64(rtime)) / np.timedelta64(1, 'D')
#         ncid = nc.Dataset(bryfile, "r+", format="NETCDF4")
#         ncid.variables['zeta_time'][ind]=timeout
#         ncid.variables['v2d_time'][ind]=timeout

        
#         ncid.variables['ubar_north'][ind,:]=xromqckL1_N.ubar.values[-1,:]
#         ncid.variables['vbar_north'][ind,:]=xromqckL1_N.vbar.values[-1,:]
#         ncid.variables['zeta_north'][ind,:]=xromqckL1_N.zeta.values[-1,:]

#         ncid.variables['ubar_south'][ind,:]=xromqckL1_S.ubar.values[0,:]
#         ncid.variables['vbar_south'][ind,:]=xromqckL1_S.vbar.values[0,:]
#         ncid.variables['zeta_south'][ind,:]=xromqckL1_S.zeta.values[0,:]


#         ncid.variables['ubar_west'][ind,:]=xromqckL1_W.ubar.values[:,0]
#         ncid.variables['vbar_west'][ind,:]=xromqckL1_W.vbar.values[:,0]
#         ncid.variables['zeta_west'][ind,:]=xromqckL1_W.zeta.values[:,0]


#         ncid.variables['ubar_east'][ind,:]=xromqckL1_E.ubar.values[:,-1]
#         ncid.variables['vbar_east'][ind,:]=xromqckL1_E.vbar.values[:,-1]
#         ncid.variables['zeta_east'][ind,:]=xromqckL1_E.zeta.values[:,-1]


#         ncid.sync()
#         ncid.close()        
        
        
def rechunk_eta_xi(ds):
    ds=ds.unify_chunks()   
    ds = ds.chunk({'xi_rho': ds.sizes['xi_rho']})
    ds = ds.chunk({'eta_rho': ds.sizes['eta_rho']})
    ds = ds.chunk({'eta_v': ds.sizes['eta_v']})
    ds = ds.chunk({'xi_u': ds.sizes['xi_u']})
    return ds

def dataset_get_Z(dataout,varnew):

    #Fixing grid
    # xromsds['s_w']=('s_w',s_w_L0)
    # xromsds['Cs_w']=('s_w',C_w_L0)
    # xromsds['s_rho']=('s_rho',s_r_L0)
    # xromsds['Cs_r']=('s_rho',C_r_L0)
    # xromsds.attrs['hc']=L0hc
    # xromsds=xromsds.unify_chunks()   
    # xromsds = xromsds.chunk({'xi_rho': xromsds.sizes['xi_rho']})
    # xromsds = xromsds.chunk({'eta_rho': xromsds.sizes['eta_rho']})
    # xromsds = xromsds.chunk({'eta_v': xromsds.sizes['eta_v']})
    # xromsds = xromsds.chunk({'xi_u': xromsds.sizes['xi_u']})

    #(xromsds,gridhisL1)=xroms.roms_dataset( xromsds,Vtransform=L0Vtransform)
    varnew['zeta']=dataout['zos']
    xromsds_z=varnew.zeta.to_dataset()
    xromsds_z['h']=varnew['h']
    xromsds_z['xi_u']=dataout['xi_u']
    xromsds_z['eta_v']=dataout['eta_v']
        
    xromsds_z['s_w']=('s_w',s_w_L1)
    xromsds_z['Cs_w']=('s_w',C_w_L1)
    xromsds_z['s_rho']=('s_rho',s_r_L1)
    xromsds_z['Cs_r']=('s_rho',C_r_L1)
    xromsds_z.attrs['hc']=L1hc
    xromsds_z['pm']=varnew['pm']
    xromsds_z['pn']=varnew['pn']
    
    (xromsds_z,gridhisL1_z)=xroms.roms_dataset(xromsds_z,Vtransform=L1Vtransform)
    
    

    return xromsds_z       


        

def plot_2D_LL(lon,lat,var):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1,projection=proj)
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines(resolution='50m')

    ax.set_extent([lonmin, lonmax, latmin,latmax])
    
    
    pcid=ax.pcolormesh(lon,lat,var,transform=pc)
    colorbar = plt.colorbar(pcid)
    plt.show()
    
def plot_points_LL(lon,lat,var):
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1,1,1,projection=proj)
    gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
    gl.top_labels = False
    gl.right_labels = False
    ax.coastlines(resolution='50m')

    ax.set_extent([lonmin, lonmax, latmin,latmax])
        
  
    blon,blat = np.meshgrid(var['longitude'].values,var['latitude'].values)
    #pcid=ax.scatter(lon,lat,c=1,s=20,transform=pc,marker='o',vmin=0,vmax=1)
    print(blon)
   # 
    pcid2=ax.plot(blon,blat,'k+',transform=pc)
    pcid=ax.plot(lat,lon,'ro',transform=pc)
  #  colorbar = plt.colorbar(pcid)
    plt.show()
    
    
    
    
if __name__ == "__main__":
    print('Running Downscale')
    main()
    
    print('Finished Downscale')
