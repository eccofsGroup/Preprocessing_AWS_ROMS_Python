"""
DOWNSCALE_MERCATOR_TO ROMS.py

Adapted from DOWNSCALE_ROMS_OUTPUT.py
Created by Elias Hunter, hunter@marine.rutgers.edu, 12/19/2023 
"""
import os
import numpy as np
import xroms
import xarray as xr
import xesmf as xe
import datetime,time
from datetime import timedelta,date
import romsutil as rutil
from  scipy.interpolate import interp1d
import netCDF4 as nc


import warnings
warnings.filterwarnings("ignore")
#Set relevant inout parameters

# proj = cartopy.crs.Mercator(central_longitude=-74)
# pc = cartopy.crs.PlateCarree()


# latmin=-3.0
# latmax=53.0
# lonmin=-100.0
# lonmax=-38.0

#Set time



reftime=datetime.datetime(2011,1,1)
tunits=reftime.strftime('days since %Y-%m-%d')
rtime=np.datetime64(reftime.date())

#Donor grid Info


#Donor inputfiles
# datadir='/home/om/cron/ECCOFS_OBS/MERCATOR/data/raw/'
# regrid_coef_file='/home/om/cron/ECCOFS_OBS/MERCATOR/data/mercator_bdry_9999.nc'


# #Receiver Grid Info
# L1grdfile='/home/om/cron/ECCOFS_OBS/MERCATOR/work/grid_eccofs_3km_08_b7.nc' # can be a thredds url
# L1theta_s=7.0
# L1theta_b=2.0
# L1Tcline=250.0
# L1N=50
# L1Vtransform  =2        #vertical transformation equation
# L1Vstretching =4        #vertical stretching function
# L1hc=250
# Nbed=1
# Nveg=1
# NCS=0
# NNS=1

# (s_w_L1,C_w_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,1)
# (s_r_L1,C_r_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,0)




#Output fflags

# bdylist={'EAST':True, 'WEST':False, 'NORTH':False, 'SOUTH':True}
#bdylist={'SOUTH':True}




#Receiver Output file
#bryfile=f'/home/om/cron/ECCOFS_OBS/MERCATOR/data/bdry/ECCOFS_BDRY_{sdate}_bry.nc'

############################################################################
#Main Program
############################################################################

def main(fconfig):
    ########################################################################
    print('Initilaizing grids and regridder')
    ########################################################################
  #  regrid_coef_file=fconfig['force']['bry']['regrid_coef_file']
 #   today=day.strftime('%Y%m%d')
  #  bryfile=f'{fconfig['force']['bry']['ibrydir']}{fconfig['force']['bry']['inipre']}{today}_BRY.nc'
    L1N=fconfig['force']['L1N']
    L1grdfile=fconfig['force']['L1grdfile']
    datadir=fconfig['force']['datadir']
    days=fconfig['force']['bry']['days']
    hdays=fconfig['force']['bry']['hdays']
    fdays=fconfig['force']['bry']['fdays']
    nday=hdays+fdays
    #Getting the mercator grid information
    middate = date.today()-timedelta(days=days)
    alldates = [middate  + timedelta(days=i) for i in range(days)]  

    for start_date in alldates:
        print(start_date)
        date_times = [start_date  + timedelta(days=i) for i in range(nday)]  
        cfgrd=xroms.open_netcdf(L1grdfile)
        cfgrd.attrs['sc_r']=L1N
        cfgrd.attrs['sc_w']=L1N+1
        #plot_points_LL
        filelist=[]
        for day in date_times:
            
            mfilename = day.strftime(datadir+'mercator_doppio_%Y_%m_%d.nc')
            
            filelist.append(mfilename)
    
      
        start_time=time.time()
        dsmerc=xr.open_mfdataset(filelist)
        dsmerc['zos']= dsmerc['zos'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
        dsmerc['uo']= dsmerc['uo'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
        dsmerc['vo']= dsmerc['vo'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
        dsmerc['so']= dsmerc['so'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
        dsmerc['thetao']= dsmerc['thetao'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
      
        end_time=time.time()
        elapsed_time = end_time - start_time
        print(f"filling processing time: {elapsed_time} seconds")

        ########################################################################
        #Process boudnary files
        ########################################################################
    
    
    
        start_time=time.time()
        downscale_bdry_file(cfgrd,dsmerc,start_date,fconfig)
        end_time=time.time()
        elapsed_time = end_time - start_time
        print(f"boundary conditions  file processing time: {elapsed_time} seconds")
        
def downscale_bdry_file(cfgrd,dsmerc,day,fconfig):

    regrid_coef_file=fconfig['force']['bry']['regrid_coef_file']
    today=day.strftime('%Y%m%d')
    bryfile=f'{fconfig['force']['bry']['brydir']}{fconfig['force']['bry']['brypre']}{today}_bry.nc'
    bdylist=fconfig['force']['bry']['bdylist']

    
   
    L1N=fconfig['force']['L1N']
  #  L1Vstretching=fconfig['force']['L1Vstretching']
    L1Vtransform=fconfig['force']['L1Vtransform']
 #   L1theta_s=fconfig['force']['L1theta_s']
  #  L1theta_b=fconfig['force']['L1theta_b']
  #  L1hc=fconfig['force']['L1hc']
  #  L1Tcline=fconfig['force']['L1Tcline']
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
            print(f'Regrid {key} boundary: {elapsed_time}')
    



        
# ##############################################################################################
# #Get depths of L0 and L1 grid. 
# ##############################################################################################
        
            dataout_z=dataset_get_Z(dataout,varnew,fconfig)
    
    
            dim_dict=varnew.dims
##############################################################################################
#Initialize Arrays
##############################################################################################
         #   print('INITIALZING OUTPUT ARRAYS')
            if key=='EAST' :
                temp=np.full((L1N,dim_dict['eta_rho'],3),0.0)    
                salt=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
                u_east=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
                v_north=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
            if key=='WEST' :
                temp=np.full((L1N,dim_dict['eta_rho'],3),0.0)    
                salt=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
                u_east=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
                v_north=np.full((L1N,dim_dict['eta_rho'],3),0.0)   
            if key=='NORTH' :
                temp=np.full((L1N,3,dim_dict['xi_rho']),0.0)    
                salt=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
                u_east=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
                v_north=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
            if key=='SOUTH' :
                temp=np.full((L1N,3,dim_dict['xi_rho']),0.0)    
                salt=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
                u_east=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
                v_north=np.full((L1N,3,dim_dict['xi_rho']),0.0)   
           

# # tlat=[]
# # tlon=[]
# ##############################################################################################
# #Loading Data for interpolation
# ##############################################################################################
#            print('LOADING DATA FOR VERTICAL INTERP.')
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
         #   print(temp.shape)
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
                
                


            ncid.sync()
            ncid.close()
        
        
def rechunk_eta_xi(ds):
    ds=ds.unify_chunks()   
    ds = ds.chunk({'xi_rho': ds.sizes['xi_rho']})
    ds = ds.chunk({'eta_rho': ds.sizes['eta_rho']})
    ds = ds.chunk({'eta_v': ds.sizes['eta_v']})
    ds = ds.chunk({'xi_u': ds.sizes['xi_u']})
    return ds

def dataset_get_Z(dataout,varnew,fconfig):

    L1N=fconfig['force']['L1N']
    L1Vstretching=fconfig['force']['L1Vstretching']
    L1Vtransform=fconfig['force']['L1Vtransform']
    L1theta_s=fconfig['force']['L1theta_s']
    L1theta_b=fconfig['force']['L1theta_b']
    L1hc=fconfig['force']['L1hc']
    (s_w_L1,C_w_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,1)
    (s_r_L1,C_r_L1)=rutil.stretching(L1Vstretching,L1theta_s,L1theta_b,L1hc,L1N,0)
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


        
    
    
# if __name__ == "__main__":
#     print('Running Downscale')
#     main()
    
#     print('Finished Downscale')
