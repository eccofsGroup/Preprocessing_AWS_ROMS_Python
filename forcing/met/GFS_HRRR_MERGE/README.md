Files to download, regrid , and merge GFS and HRRR forecasts for use as forcing inpout to ECCOFS. 

Files:
GFS_process_v2.py: Downloads GFS forecast. All the data are accessed via Kerchunk/scangrib and  fsspec. Radiation variables are treated differently because of averaging.  

HRRR_process_v2.py: Downloads HRRR forecast. All the non-radiation variales are accessed directly via s3fs. Radtiation variables are accessed using Herbie, 

REGRID_GFS_HRRR_ECCOFS_V1.py: Merges GFS and HRRR data. GFS is regridding to high resolution and blended with HRRR via weighted averaging. 


At the top of all three of these files there are directories that need to be set. Then run GFS_process_v2 and HRRR_process_v2 first. Then REGRID_GFS_HRRR_ECCOFS_V1. 

The conda environments I used are in GFS_HRRR_DOWNLOAD.yml (for GFS_process_v2 and HRRR_process_v2) and HRRR_GFS_MERGE.yml (for REGRID_GFS_HRRR_ECCOFS_V1). The reason there are two enviroments is that they are run on different machines at Rutgers. You can creatre one environnment to run it all. 



