#!/bin/bash
#

source ~/.bashrc
conda activate ECCOFS

echo STARTING ECCOFS driver
cd /home/om/cron/ECCOFS_OBS/REALTIME/work/

python -u /home/om/cron/ECCOFS_OBS/REALTIME/work/eccofs_obs_driver.py 
echo MOVING FILES
python -u /home/om/cron/ECCOFS_OBS/REALTIME/work/eccofs_obs_driver_move.py 



