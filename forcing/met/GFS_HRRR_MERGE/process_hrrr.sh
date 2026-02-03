#!/bin/bash
#
source /home/ubuntu/.bashrc_conda
conda activate HRRR 
echo STARTING

scp /home/ubuntu/data/hrrr/HRRR_ZARR_REGRIDED*.nc hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/data/hrrr/
scp /home/ubuntu/cron/latest_hrrr.log hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/

rm /home/ubuntu/data/hrrr/HRRR_ZARR_REGRIDED*.nc

ipython /home/ubuntu/python/HRRR/HRRR_process_v2.py

scp /home/ubuntu/data/hrrr/HRRR_ZARR_REGRIDED*.nc hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/data/hrrr/
scp /home/ubuntu/cron/latest_hrrr.log hunter@boardwalk.marine.rutgers.edu:/home/om/cron/HRRR_AWS/


echo FINISHED
