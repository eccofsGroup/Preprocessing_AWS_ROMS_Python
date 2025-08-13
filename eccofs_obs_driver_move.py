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

    file_delete_driver(fconfig)
    file_compile_driver(fconfig)
    file_transfer_driver(fconfig)
    et=time.time()
    elt=et-st
    print(f'TOTAL time file transfer: {elt} seconds')

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
