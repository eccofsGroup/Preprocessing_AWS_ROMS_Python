#!/usr/bin/env python
#
import os,sys,fnmatch,glob
from shutil import copyfile
from datetime import datetime,timedelta 
import subprocess
import pandas as pd


rads2asc='/home/om/cron/RADS4/local/bin/rads2asc4'
outdir='/home/om/cron/ECCOFS_OBS/SSH/data/raw/'
days=90
latlon='-100,-38,5,52'
variables='time,lat,lon,sla,dry_tropo,wet_tropo,iono,inv_bar,tide_ocean_got410,tide_load_got410,ssb,mss_dtu15'


satlist=['j3/b','j3/c','j3/d','3a/a','3b/b','6a/a','sw/b']
dirlist=['jason3/','jason3/','jason3/','sentinel3a/','sentinel3b/','sentinel6a/','swot/']




# dry_tropo: dry tropo ecmwf
# wet_tropo: wet tropo rad wet tropo ecmwf other
# wet_tropo: wet tropo ecmwf c2
# iono: iono gim iono nic09 c2 e2 g1 sa
# iono: iono gim iono nic09 j1 j2 j3 n1 tx 9
# inv bar: inv bar mog2d  #NOT USED
#tide ocean got410
#tide load got410
#ssb: ssb cls j1 j2 j3 n1 tx
#ssb ssb cls j1 j2 j3 n1 tx
#mss : mss dtu15


def main(fconfig):
     print('MAIN')
     rads2asc=fconfig['obs']['ssh']['rads2asc']
     odirdac=fconfig['obs']['ssh']['odirdac']
     odirnodac=fconfig['obs']['ssh']['odirnodac']
     days=fconfig['obs']['ssh']['days']
     latlon=fconfig['obs']['ssh']['latlon']
     variables=fconfig['obs']['ssh']['variables']
     satlist=fconfig['obs']['ssh']['satlist']
     dirlist=fconfig['obs']['ssh']['dirlist']
     radsxmlnodac=fconfig['obs']['ssh']['radsxmlnodac']
     radsxmldac=fconfig['obs']['ssh']['radsxmldac']
     tend=datetime.now()
     tstart=tend-timedelta(days=days)
     
 
     
     timeend=tend.strftime( "%Y%m%d%H%M%S" )
     timestart=tstart.strftime( "%Y%m%d%H%M%S" )
#      print(timestart)
#      print(timeend)
     print('GETTING SL WITH DAC')
     for il,sat in enumerate(satlist):
             print(sat)
             subprocess.check_call([rads2asc,
                   '-S', sat,
                   '--ymd',timestart+','+timeend,
                   '-R',latlon,
                   '-V',variables,
                   '-X',radsxmldac,
                   '-o',odirdac+dirlist[il]])
            
     
     for il,sat in enumerate(satlist):
             print(sat)
             subprocess.check_call([rads2asc,
                   '-S', sat,
                   '--ymd',timestart+','+timeend,
                   '-R',latlon,
                   '-V',variables,
                   '-X',radsxmlnodac,
                   '-o',odirnodac+dirlist[il]])
            

