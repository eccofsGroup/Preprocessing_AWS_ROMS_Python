#!/usr/bin/env python
import yaml,os,sys,glob
import optparse
import datetime
from datetime import datetime,timedelta,date
import time
import pprint,traceback
import sched
import warnings
warnings.filterwarnings("ignore")

# import DOWNSCALE_ROMS_FORECAST

def main(opts):
    print('OBS driver initializing')
    st=time.time()
    cfgfile=opts.configfile
    fconfig=read_config(cfgfile)
    # #pprint.pprint(fconfig)
    for path in fconfig['modulepath']:
        if path not in sys.path:
            sys.path.append(path)

    # print('Get GLOFAS Rivers Data')
    # status_river=ECCOFS_RIVERS_driver(fconfig)
    

    # print('Get Mercator Data')
    # status_merc=mercator_aquire_driver(fconfig)
    
    
    # # print('Extract Initial conditions from mercator')
    # # status_ini=ECCOFS_INI_driver(fconfig)
    
    print('Extract climatology data from mercator')
    status_clm=ECCOFS_CLM_driver(fconfig)
    
    # print('Extract Boundary conditions from mercator')
    # status_bry=ECCOFS_BRY_driver(fconfig)
    
    
    # print('In situ driver running')
    # status_insitu=insitu_driver(fconfig)
    
    
    # print('Rads Driver Running')
    # status_ssh=rads_driver(fconfig)
    
    # print('Checking SST Availability')
    # status_sst=sst_driver(fconfig)
    # #print(f'There are a total of  {filecount} SST files')

    # print('Checking GLIDER Availability')
    # glider_driver(fconfig)


    # print('RUNNING OBS PreProcessing')
    # status=obs_pre_driver(fconfig)

    # print('RUNNING CombineOBS')
    # status=obs_combine_driver(fconfig)


    # #create scheduler
    # delay=fconfig['obs']['delay']
    # scheduler = sched.scheduler(time.time, time.sleep)


    # if not(status_insitu):
    #     print(f'In situ data aquisition failed, waiting {delay} seconds')
    #     scheduler.enter(delay, 2, insitu_driver,(fconfig,))
    # else:
    #     print('In situ data aquisition succeeded')
            
            
    # if not(status_ssh):
    #     print(f'SSH data aquisition failed, waiting {delay} seconds')
    #     scheduler.enter(delay,1, rads_driver,(fconfig,))
    # else:
    #     print('SSH data aquisition succeeded')

    
    # if not(status_sst):
    #     print(f'SST data aquisition failed, waiting {delay} seconds')
    #     scheduler.enter(delay,1, sst_driver,(fconfig,))
    # else:
    #     print('SST data aquisition succeeded')


    # scheduler.run()
    # et=time.time()
    # elt=et-st
    # print(f'TOTAL processing time: {elt} seconds')
        
    # file_delete_driver(fconfig)
    # file_compile_driver(fconfig)
    # file_transfer_driver(fconfig)
    # et=time.time()
    # elt=et-st
    # print(f'TOTAL processing time with file transfer: {elt} seconds')

def file_compile_driver(fconfig):
    import fileutil as fileu
    
    fileu.compile_all_files(fconfig)
    

def file_delete_driver(fconfig):
    import fileutil as fileu
    
    fileu.delete_remote(fconfig)
    fileu.delete_local(fconfig)
    
def file_transfer_driver(fconfig):
    import fileutil as fileu
    
    fileu.transfer(fconfig)
    
    
    
def obs_combine_driver(fconfig):
    import combine_obs as combine
    
        
    print('Processing Combine')
    try:    
        st=time.time()
        combine.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'Combine processing time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
            
    return status        
def ECCOFS_RIVERS_driver(fconfig):
    import GET_GLOFAS_FOR_ECCOFS_DAILY as ECCOFS_river
    try:    
        st=time.time()
        ECCOFS_river.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'Rivers file time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    return status   
    
def ECCOFS_BRY_driver(fconfig):
    import DOWNSCALE_MERCATOR_TO_ROMS_BOUNDARY as ECCOFS_bry
    try:    
        st=time.time()
        ECCOFS_bry.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'BRY file time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    return status   
   
def ECCOFS_INI_driver(fconfig):
    import DOWNSCALE_MERCATOR_TO_ROMS_INI as ECCOFS_ini
    try:    
        st=time.time()
        ECCOFS_ini.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'IC file time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    return status   

def ECCOFS_CLM_driver(fconfig):
    import DOWNSCALE_MERCATOR_TO_ROMS_CLM as ECCOFS_clm
    try:    
        st=time.time()
        ECCOFS_clm.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'CLM file time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    return status   
   
def mercator_aquire_driver(fconfig):
    import acquire_mercator_ECCOFS as get_merc
    try:    
        st=time.time()
        get_merc.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'Mercator Acquisition time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    return status        
def obs_pre_driver(fconfig):
    import get_sst_amsr2 as amsr2
    import get_sst_goes as GOES
    import get_sst_leo as LEO
    import get_ssh as SSH
    import get_cmems as CMEMS
    
        
    print('Processing CMEMS')
    try:    
        st=time.time()
        CMEMS.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'CMEMS processing time: {elt} seconds')
        
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
            
            
    print('Processing SSH')
    try:    
        st=time.time()
        SSH.main(fconfig)
        status=True  
        et=time.time()
        elt=et-st
        print(f'SSH processing time: {elt} seconds')
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
            
    print('Processing AMSR2')
    try:    
        amsr2.main(fconfig)
        status=True  
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
    
    print('Processing GOES')
    try:    
        GOES.main(fconfig)
        status=True  
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 
     
    print('Processing LEO')
    try:    
        LEO.main(fconfig)
        status=True  
    except Exception as e:
            print(f"An error occurred: {e}")
            traceback.print_exc()  
            status=False 

                           
                   
    return status     
def rads_driver(fconfig):
    import acquire_eccofs_ssh_rads_v1 as sshaquire
    print('--------------------')
    
    try:
        sshaquire.main(fconfig)
        status=True    
    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()  
        status=False   
    
    return status       
def sst_driver(fconfig):
    import acquire_GOES19_sst_v1 as goesacquire
    print('--------------------')
    
    try:
        goesacquire.main(fconfig)
        status=True    
    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()  
        status=False   
    
    direc=fconfig['obs']['sst']['direc']
    subs=fconfig['obs']['sst']['subdir']
    ndays=fconfig['obs']['sst']['ndays']
    dt=timedelta(days=ndays)
    start_date = date.today()-dt
    fcount=0
    for ind,sub in enumerate(subs):
        sdirec=os.path.join(direc,sub,'*.nc')
        filelist=glob.glob(sdirec)
        filelist=sorted(filelist)
        lastfiles=filelist[-ndays:]
        nfiles=0
        for il in lastfiles:
            sdate=il.split('/')[-1][0:8]
            fdate=datetime.strptime(sdate, "%Y%m%d").date()
            if fdate>start_date:
                    nfiles=nfiles+1
                    fcount=fcount+1
                    #print(il)
        print(f' {subs[ind]} has {nfiles} files')
    print(f'There are a total of  {fcount} SST files')
    return status
def insitu_driver(fconfig):
    import acquire_eccofs_insitu_obs_v1 as isaquire
    
    
    print('--------------------')
    try:
        isaquire.main(fconfig)
        status=True    
    except Exception as e:
        print(f"An error occurred: {e}")
        traceback.print_exc()  
        status=False   

    
    return status
def glider_driver(fconfig):
    print('--------------------')
    print('glider obs placeholder')
def read_config(infile):
    # Reading YAML from a file
    with open(infile, 'r') as file:
        #with open('./test.yml', 'r') as file:
        fconfig = yaml.safe_load(file)

    
    
    return fconfig
    

    

if __name__ == "__main__":
    parser = optparse.OptionParser()
    parser.add_option('-c', '--configfile',dest='configfile',help='Configuration YML file',default='eccofs_driver_config.yml',type='str')
    (opts, args) = parser.parse_args()
    
    print('RUNNING DRIVER')
    main(opts)
    print('DONE')
