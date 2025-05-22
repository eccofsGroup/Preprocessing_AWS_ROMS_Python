#!/usr/bin/env python
import yaml,os,sys,glob
import optparse
import datetime
from datetime import datetime,timedelta,date
import time
import pprint,traceback
import sched
# import DOWNSCALE_ROMS_FORECAST

def main(opts):
    print('OBS driver initializing')
    cfgfile=opts.configfile
    fconfig=read_config(cfgfile)
    #pprint.pprint(fconfig)
    for path in fconfig['modulepath']:
        if path not in sys.path:
            sys.path.append(path)

    
    
    print('In situ driver running')
    status_insitu=insitu_driver(fconfig)
    
    
    print('Rads Driver Running')
    status_ssh=rads_driver(fconfig)
    
    print('Checking SST Availability')
    filecount=sst_driver(fconfig)
    print(f'There are a total of  {filecount} SST files')

    print('Checking GLIDER Availability')
    glider_driver(fconfig)




#
#
#
#
#
#
#
#
#
#
#




    # create scheduler
    delay=fconfig['obs']['delay']
    scheduler = sched.scheduler(time.time, time.sleep)


    if not(status_insitu):
        print(f'In situ data aquisition failed, waiting {delay} seconds')
        scheduler.enter(delay, 2, insitu_driver,(fconfig,))
    else:
        print('In situ data aquisition succeeded')
            
            
    if not(status_ssh):
        print(f'SSH data aquisition failed, waiting {delay} seconds')
        scheduler.enter(delay,1, rads_driver,(fconfig,))
    else:
        print('SSH data aquisition succeeded')

    scheduler.run()


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
    print('--------------------')
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

    return fcount
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
