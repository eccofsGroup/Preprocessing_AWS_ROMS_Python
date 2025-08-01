"""
DOWNSCALE_MERCATOR_TO ROM_CLM.py

Adapted from DOWNSCALE_ROMS_OUTPUT.py. Generate a Climatology file 
for a ROMS gird from Mercator. 
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
import psutil

import warnings
warnings.filterwarnings("ignore")
#Set relevant inout parameters


reftime=datetime.datetime(2011,1,1)
tunits=reftime.strftime('days since %Y-%m-%d')
rtime=np.datetime64(reftime.date())




############################################################################
#Main Program
############################################################################

def main(fconfig):
    ########################################################################
    print('Initilaizing grids and regridder')
    ########################################################################
    
    
    
    datadir=fconfig['force']['datadir']
    L1grdfile=fconfig['force']['L1grdfile']
    L1N=fconfig['force']['L1N']
            
     
    cfgrd=xroms.open_netcdf(L1grdfile)
    cfgrd.attrs['sc_r']=L1N
    cfgrd.attrs['sc_w']=L1N+1
    #Getting the mercator grid information
    nday=fconfig['force']['clm']['nday']
    fdays=fconfig['force']['clm']['fdays']
    start_date = date.today()-timedelta(days=nday)
    date_times = [start_date  + timedelta(days=i) for i in range(nday+fdays)]  
    
    

    for day in date_times:
        mfilename = day.strftime(datadir+'mercator_doppio_%Y_%m_%d.nc')
        print(mfilename)

        
        dsmerc=xr.open_dataset(mfilename)
        
    
                
    ########################################################################
    #Process donwscaling files
    ########################################################################
       
        start_time=time.time()
        downscale_clm_file(cfgrd,dsmerc,day,fconfig)
        end_time=time.time()
        elapsed_time = end_time - start_time
        print(f"Initilization file processing time: {elapsed_time} seconds")


def memory_usage():
    process = psutil.Process(os.getpid())
    mem_info = process.memory_info()
    return mem_info.rss / 1024 ** 2  # Memory usage in MB


        
        
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


        
def downscale_clm_file(cfgrd,dsmerc,day,fconfig):
    today=day.strftime('%Y%m%d')
    clmfile=f'{fconfig['force']['clm']['clmdir']}{fconfig['force']['clm']['clmpre']}{today}_clm.nc'
    rutil.create_clm_file(clmfile,cfgrd,tunits)
    
    regrid_coef_file=fconfig['force']['clm']['regrid_coef_file']

    L1N=fconfig['force']['L1N']
    
   
    L1N=fconfig['force']['L1N']
    L1Vstretching=fconfig['force']['L1Vstretching']
    L1Vtransform=fconfig['force']['L1Vtransform']
    L1theta_s=fconfig['force']['L1theta_s']
    L1theta_b=fconfig['force']['L1theta_b']
    L1hc=fconfig['force']['L1hc']
    L1Tcline=fconfig['force']['L1Tcline']


    tmp=dsmerc
 #   mlat=dsmerc['latitude']
 #   mlon=dsmerc['longitude']
    varnew =cfgrd
    varnew=varnew.rename({'lat_rho':'lat'})
    varnew=varnew.rename({'lon_rho':'lon'})
    
    start_time=time.time()

    dsmerc['zos']= dsmerc['zos'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
    dsmerc['uo']= dsmerc['uo'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
    dsmerc['vo']= dsmerc['vo'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
    dsmerc['so']= dsmerc['so'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
    dsmerc['thetao']= dsmerc['thetao'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")
  #  dsmerc['vbar']= dsmerc['vbar'].interpolate_na(dim="longitude", method="nearest",limit=None,fill_value="extrapolate").interpolate_na(dim="latitude", method="nearest",limit=None,fill_value="extrapolate")

    end_time=time.time()
    elapsed_time = end_time - start_time
    print(f'FILL time: {elapsed_time}')
#    varnew=varnew.rename({'mask_rho':'mask'})
    start_time=time.time()
    
    if os.path.exists(regrid_coef_file):
        print('Re-using weights')
        regridder = xe.Regridder( tmp,varnew, "bilinear",extrap_method="nearest_s2d",filename=regrid_coef_file,reuse_weights=True)  
    else:
        print('Generating weights')
        regridder = xe.Regridder( tmp,varnew, "bilinear",extrap_method="nearest_s2d")
        regridder.to_netcdf(regrid_coef_file)
            

    
    end_time=time.time()
    elapsed_time = end_time - start_time
    print(f'REGRIDDER time: {elapsed_time}')
    
    

    
    print('RUNNING horizontal regridding')
 
    dataout = regridder(dsmerc,keep_attrs=True)
 
    dataout['xi_u']=varnew['xi_u']
    dataout['eta_v']=varnew['eta_v']
    dataout = dataout.rename({"time": "ocean_time"})
    print('Preparing for vertical interpolation')
   

    dataout=rechunk_eta_xi(dataout)    
    varnew=rechunk_eta_xi(varnew)
    dataout_z=dataset_get_Z(dataout,varnew,fconfig)
    # print(dataout_z['s_w'].values)
    # print(dataout_z['s_rho'].values)
    # print(dataout_z['Cs_w'].values)
    # print(dataout_z['Cs_r'].values)
    dim_dict=dataout_z.dims
    
    temp = np.empty((L1N,dim_dict['eta_rho'],dim_dict['xi_rho']))
    temp[:]=0.0
    salt = np.empty((L1N,dim_dict['eta_rho'],dim_dict['xi_rho']))
    salt[:]=0.0
    u_east = np.empty((L1N,dim_dict['eta_rho'],dim_dict['xi_rho']))
    u_east[:]=0.0
    v_north = np.empty((L1N,dim_dict['eta_rho'],dim_dict['xi_rho']))
    v_north[:]=0.0
    
    
    
    print('Running vertical interpolation, this may take a few minutes')
    
    mask=cfgrd.mask_rho.values
    tmpz=dataout_z.z_rho.values
    tmpz2=dataout['depth']
    temp_I=dataout.thetao.values
    salt_I=dataout.so.values
    u_east_I=dataout.uo.values
    v_north_I=dataout.vo.values
    
    
    start_time=time.time()


    for eta in range(0,dim_dict['eta_rho']):
        for xi in range(0,dim_dict['xi_rho']):
            maskflag=mask[eta,xi]
            if maskflag==0.0:
                continue
            

        #    print('----------------------')
                   
            all_nan = np.all(np.isnan(temp_I[0,:,eta,xi]))
            if ~all_nan:
              #  print('TEMP NAN')
                
                z=-tmpz2.values     
                data=temp_I[0,:,eta,xi]
                Igood=np.where(~np.isnan(data))
                ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                ntmp=ifun(tmpz[0,:,eta,xi])
                temp[:,eta,xi]=ntmp
                
                
                
            all_nan = np.all(np.isnan(salt_I[0,:,eta,xi]))
            if ~all_nan:
            #    print('Salt NAN')
                z=-tmpz2.values
                data=salt_I[0,:,eta,xi]
                Igood=np.where(~np.isnan(data))
                ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                ntmp=ifun(tmpz[0,:,eta,xi])
                salt[:,eta,xi]=ntmp
                

            
            all_nan = np.all(np.isnan(u_east_I[0,:,eta,xi]))
            if ~all_nan:
            #    print('U NAN')
                z=-tmpz2.values
                data=u_east_I[0,:,eta,xi]
                Igood=np.where(~np.isnan(data))
                ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                ntmp=ifun(tmpz[0,:,eta,xi])
                u_east[:,eta,xi]=ntmp
            
            all_nan = np.all(np.isnan(v_north_I[0,:,eta,xi]))
            if ~all_nan:
               # print('V NAN')
                z=-tmpz2.values
                data=v_north_I[0,:,eta,xi]
                Igood=np.where(~np.isnan(data))
                ifun=interp1d(z[Igood],data[Igood],kind='linear',bounds_error=False,fill_value=(data[Igood][-1],data[Igood][0]))
                ntmp=ifun(tmpz[0,:,eta,xi])
                v_north[:,eta,xi]=ntmp        


            
    end_time=time.time()
    elapsed_time = end_time - start_time
    print(f"processing time: {elapsed_time} seconds") 


    
    
  #   print('Rotating velocities to receiver grid coordinate system. ')   
    dataout_z['temp']=(('s_rho', 'eta_rho', 'xi_rho'),temp)
    dataout_z['salt']=(('s_rho', 'eta_rho', 'xi_rho'),salt) 
    dataout_z['u_eastward']=(('s_rho', 'eta_rho', 'xi_rho'),u_east) 
    dataout_z['v_northward']=(('s_rho', 'eta_rho', 'xi_rho'),v_north) 
  #   dataout_z['vbar_northward']=xromqckL1['vbar_northward']
  #   dataout_z['ubar_eastward']=xromqckL1['ubar_eastward']
    dataout_z['zeta']=dataout_z['zeta'].fillna(0)
    (dataout_z,gridout)=xroms.roms_dataset(dataout_z,Vtransform=L1Vtransform)
        

    uv=rutil.uv_rot(dataout_z.u_eastward, dataout_z.v_northward, gridout,cfgrd.angle,reverse=True)
    ruhis=uv[0]
    rvhis=uv[1]
    dataout_z=xr.merge([dataout_z,ruhis,rvhis])
    dataout_z.load()
    ubar=dataout_z.u.xroms.gridmean(gridout,"Z")
    vbar=dataout_z.v.xroms.gridmean(gridout,"Z")


    print('Writing climatology data to file'+clmfile)
    ncid = nc.Dataset(clmfile, "r+", format="NETCDF4")
    
    t=dsmerc.time.values
    timeout=(t-np.datetime64(rtime)) / np.timedelta64(1, 'D')
    ncid.variables['ocean_time'][:]=timeout

    ncid.variables['ubar'][0,:,:]=ubar.values[:,:]
    ncid.variables['vbar'][0,:,:]=vbar.values[:,:]
    ncid.variables['zeta'][0,:,:]=dataout['zos'].values[:,:]
    ncid.variables['u'][0,:,:,:]=dataout_z.u.values[:,:,:]
    ncid.variables['v'][0,:,:,:]=dataout_z.v.values[:,:,:]
    ncid.variables['salt'][0,:,:,:]=salt
    ncid.variables['temp'][0,:,:,:]=temp
    ncid.variables['Vtransform'][:]=L1Vtransform
    ncid.variables['Vstretching'][:]=L1Vstretching
    ncid.variables['theta_b'][:]=L1theta_b
    ncid.variables['theta_s'][:]=L1theta_s
    ncid.variables['Tcline'][:]=L1Tcline
    ncid.variables['hc'][:]=L1hc
    
    

    ncid.variables['s_rho'][:]=dataout_z['s_rho'].values
    ncid.variables['s_w'][:]=dataout_z['s_w'].values
    ncid.variables['Cs_r'][:]=dataout_z['Cs_r'].values
    ncid.variables['Cs_w'][:]=dataout_z['Cs_w'].values
    
    ncid.sync()
    ncid.close()

# def plot_2D_LL(lon,lat,var):
#     fig = plt.figure(figsize=(10,6))
#     ax = fig.add_subplot(1,1,1,projection=proj)
#     gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False
#     ax.coastlines(resolution='50m')

#     ax.set_extent([lonmin, lonmax, latmin,latmax])
    
    
#     pcid=ax.pcolormesh(lon,lat,var,transform=pc)
#     colorbar = plt.colorbar(pcid)
#     plt.show()
    
# def plot_points_LL(lon,lat,var):
#     fig = plt.figure(figsize=(10,6))
#     ax = fig.add_subplot(1,1,1,projection=proj)
#     gl = ax.gridlines(draw_labels=True, x_inline=False, y_inline=False)
#     gl.top_labels = False
#     gl.right_labels = False
#     ax.coastlines(resolution='50m')

#     ax.set_extent([lonmin, lonmax, latmin,latmax])
        
        
#     pcid=ax.scatter(lon,lat,c=var,s=20,transform=pc,marker='o',vmin=0,vmax=30)
#     colorbar = plt.colorbar(pcid)
#     plt.show()
    
    
    
    
# if __name__ == "__main__":
#     print('Running Downscale')
#     main()
    
#     print('Finished Downscale')
