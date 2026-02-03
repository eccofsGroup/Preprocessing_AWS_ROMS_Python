#!/bin/bash
#
source /home/ubuntu/.bashrc_conda
conda activate HRRR 
echo STARTING

scp /home/ubuntu/data/gfs/GFS_GRIB_AWS*.nc hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/data/gfs/
scp /home/ubuntu/cron/latest_gfs.log hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/

rm /home/ubuntu/data/gfs/GFS_GRIB_AWS*.nc

ipython /home/ubuntu/python/HRRR/GFS_process_v2.py

scp /home/ubuntu/data/gfs/GFS_GRIB_AWS*.nc hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/data/gfs/
scp /home/ubuntu/cron/latest_gfs.log hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/


echo FINISHED
