#####################################################################################################################
#####################################################################################################################
# sample TROPOMI NO2 observations at African seaports

import os
import glob
import numpy as np
import pandas as pd
import xarray as xr

#####################################################################################################################
# import AF TROPOMI NO2 data extracted from raw TROPOMI L2 swath observations (see "TROPOMI_NO2_1_extract_raw_data.ipynb")

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_extracted')
NO2_AF_files = sorted(glob.glob('NO2_AF_*.nc')) 

print("Number of files:",len(NO2_AF_files))
print("First file:",NO2_AF_files[0])
print("Last  file:",NO2_AF_files[-1])

#####################################################################################################################
# read all AF observations from imported files
NO2_AF_data = [xr.open_dataset(data) for data in NO2_AF_files]

# convert to pandas dataframe
NO2_AF_data = [data.to_dataframe() for data in NO2_AF_data]

# explore the imorted data
for i in range(len(NO2_AF_data)):
    print("number of observations:",NO2_AF_data[i]['date'][0].split('-')[0:2],len(NO2_AF_data[i]))

# check a sample data
NO2_AF_data[0].head()

#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# now start regriding over each city
# 1> define an area covering the city (use a box for now)
# 2> subset the EU dataset (zoom in, quality flag, wind conditions)
# 3> regrid

#####################################################################################################################
# build a "regridding" funtion before looking at any city
# but this regreding approach does not consider numbers of observations in each grid! need to update this!

def regrid_TROPOMI(data):
    """ Aim: given a pandas dataframe (lon,lat,NO2,Pre,...), output a new dataframe (lon_regrid,lat_regrid,NO2_regrid) 
        Note: 1> every time need to set "out_lon" and "out_lat" for different grids
              2> the values of "NO2_regrid" depend on the chosen averaging method """
    data.insert(2, "lon_new", np.nan) # insert an empty column for the new lon 
    data.insert(3, "lat_new", np.nan) # insert an empty column for the new lat 
    for i in range(len(data)):        # for each point in your raw data, find its nearest point on the target grid
        data['lon_new'][i] = out_lon[np.argmin(abs(out_lon-data['lon'][i]))]
        data['lat_new'][i] = out_lat[np.argmin(abs(out_lat-data['lat'][i]))]
    data_regrid = data.groupby(['lon_new','lat_new'],as_index=False).mean() # get mean NO2 at the same grid centre
    data_regrid = data_regrid[['lon_new','lat_new','NO2']]
    return data_regrid
    
    
#####################################################################################################################
# CapeTown

# get all data over CapeTown (here lat/lon limits are boundaries, not grid centres!)
NO2_CapeTown = [data[(data['lon'] >= 17.8)  &
                     (data['lon'] <= 19)    &
                     (data['lat'] >= -34.6) &
                     (data['lat'] <= -33)] for data in NO2_AF_data]

# only get good data under all wind conditions (Flag >= 0.75)
NO2_CapeTown_good = [data[data['Flag'] >= 0.75] for data in NO2_CapeTown]

NO2_CapeTown_good = [data.reset_index(drop=True) for data in NO2_CapeTown_good] # need to reset indices everytime after slicing 

# only get good data only under calm wind conditions
NO2_CapeTown_calm = [data[(abs(data['wind_east']  <= 2)) & 
                          (abs(data['wind_north'] <= 2))] for data in NO2_CapeTown_good]

NO2_CapeTown_calm = [data.reset_index(drop=True) for data in NO2_CapeTown_calm]


# regrid data over CapeTown (list grid centres instead of the boundaries)
out_lon = np.arange(17.8+0.05/2,19,0.05)
out_lat = np.arange(-34.6+0.05/2,-33,0.05)
print("out_lon:",out_lon,
      "out_lat:",out_lat,sep="\n")

NO2_CapeTown_good_regrid = [regrid_TROPOMI(data) for data in NO2_CapeTown_good]
NO2_CapeTown_calm_regrid = [regrid_TROPOMI(data) for data in NO2_CapeTown_calm]

# save the regridding results

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_regrid')

# Gongda: these entries can all be compressed into a single line for-loop, as the only unique entries are in the number in square brackets and the file date. It would be something like this:
# yymmdd=['201909','201910','201911' ...]
# for i in range(10):
#    NO2_CapeTown_good_regrid[i].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_'+yymmdd[i]+'.csv",index=False,sep=',')

NO2_CapeTown_good_regrid[0].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[1].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[2].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[3].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[4].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[5].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[6].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[7].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[8].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_CapeTown_good_regrid[9].to_csv("NO2_CapeTown_regrid_(all winds Flag 0.75)_202006.csv",index=False,sep=',')


NO2_CapeTown_calm_regrid[0].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[1].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[2].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[3].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[4].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[5].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[6].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[7].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[8].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_CapeTown_calm_regrid[9].to_csv("NO2_CapeTown_regrid_(calm winds Flag 0.75)_202006.csv",index=False,sep=',')

# delete CapeTown data to release RAM
del (NO2_CapeTown,NO2_CapeTown_good,NO2_CapeTown_calm,NO2_CapeTown_good_regrid,NO2_CapeTown_calm_regrid)

#####################################################################################################################
# Mombasa

# get all data over Mombasa (here lat/lon limits are boundaries, not grid centres!)
NO2_Mombasa = [data[(data['lon'] >= 39.50) &
                    (data['lon'] <= 39.90)  &
                    (data['lat'] >= -4.15) &
                    (data['lat'] <= -3.9)] for data in NO2_AF_data]

# only get good data under all wind conditions (Flag >= 0.75)
NO2_Mombasa_good = [data[data['Flag'] >= 0.75] for data in NO2_Mombasa]

NO2_Mombasa_good = [data.reset_index(drop=True) for data in NO2_Mombasa_good] # need to reset indices everytime after slicing 
    
# only get good data only under calm wind conditions
NO2_Mombasa_calm = [data[(abs(data['wind_east']  <= 2)) & 
                         (abs(data['wind_north'] <= 2))] for data in NO2_Mombasa_good]

NO2_Mombasa_calm = [data.reset_index(drop=True) for data in NO2_Mombasa_calm]

# regrid data over Mombasa (list grid centres instead of the boundaries)
out_lon = np.arange(39.50+0.05/2,39.90,0.05)
out_lat = np.arange(-4.15+0.05/2,-3.9,0.05)

NO2_Mombasa_good_regrid = [regrid_TROPOMI(data) for data in NO2_Mombasa_good]
NO2_Mombasa_calm_regrid = [regrid_TROPOMI(data) for data in NO2_Mombasa_calm]

# save the regridding results

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_regrid')

NO2_Mombasa_good_regrid[0].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[1].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[2].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[3].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[4].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[5].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[6].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[7].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[8].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Mombasa_good_regrid[9].to_csv("NO2_Mombasa_regrid_(all winds Flag 0.75)_202006.csv",index=False,sep=',')


NO2_Mombasa_calm_regrid[0].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[1].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[2].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[3].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[4].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[5].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[6].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[7].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[8].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Mombasa_calm_regrid[9].to_csv("NO2_Mombasa_regrid_(calm winds Flag 0.75)_202006.csv",index=False,sep=',')

# delete Mombasa data to release RAM
del (NO2_Mombasa,NO2_Mombasa_good,NO2_Mombasa_calm,NO2_Mombasa_good_regrid,NO2_Mombasa_calm_regrid)


#####################################################################################################################
# Suez

# get all data over Suez (here lat/lon limits are boundaries, not grid centres!)
NO2_Suez = [data[(data['lon'] >= 32.25) &
                 (data['lon'] <= 32.85) &
                 (data['lat'] >= 30.85) &
                 (data['lat'] <= 31.75)] for data in NO2_AF_data]

# only get good data under all wind conditions (Flag >= 0.75)
NO2_Suez_good = [data[data['Flag'] >= 0.75] for data in NO2_Suez]

NO2_Suez_good = [data.reset_index(drop=True) for data in NO2_Suez_good] # need to reset indices everytime after slicing 

# only get good data only under calm wind conditions
NO2_Suez_calm = [data[(abs(data['wind_east']  <= 2)) & 
                      (abs(data['wind_north'] <= 2))] for data in NO2_Suez_good]

NO2_Suez_calm = [data.reset_index(drop=True) for data in NO2_Suez_calm]

# regrid data over Suez (list grid centres instead of the boundaries)
out_lon = np.arange(-0.60+0.05/2,0.30,0.05)
out_lat = np.arange(51.25+0.05/2,51.75,0.05)
print("out_lon:",out_lon,
      "out_lat:",out_lat,sep="\n")

NO2_Suez_good_regrid = [regrid_TROPOMI(data) for data in NO2_Suez_good]
NO2_Suez_calm_regrid = [regrid_TROPOMI(data) for data in NO2_Suez_calm]

# save the regridding results

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_regrid')

NO2_Suez_good_regrid[0].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Suez_good_regrid[1].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Suez_good_regrid[2].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Suez_good_regrid[3].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Suez_good_regrid[4].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Suez_good_regrid[5].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Suez_good_regrid[6].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Suez_good_regrid[7].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Suez_good_regrid[8].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Suez_good_regrid[9].to_csv("NO2_Suez_regrid_(all winds Flag 0.75)_202006.csv",index=False,sep=',')


NO2_Suez_calm_regrid[0].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Suez_calm_regrid[1].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Suez_calm_regrid[2].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Suez_calm_regrid[3].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Suez_calm_regrid[4].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Suez_calm_regrid[5].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Suez_calm_regrid[6].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Suez_calm_regrid[7].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Suez_calm_regrid[8].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Suez_calm_regrid[9].to_csv("NO2_Suez_regrid_(calm winds Flag 0.75)_202006.csv",index=False,sep=',')

# delete Suez data to release RAM
del (NO2_Suez,NO2_Suez_good,NO2_Suez_calm,NO2_Suez_good_regrid,NO2_Suez_calm_regrid)

#####################################################################################################################
# Lagos

# get all data over Lagos (here lat/lon limits are boundaries, not grid centres!)
NO2_Lagos = [data[(data['lon'] >= 1.10)  &
                  (data['lon'] <= 1.35)  &
                  (data['lat'] >= 6.05)  &
                  (data['lat'] <= 6.20)] for data in NO2_AF_data]

# only get good data under all wind conditions (Flag >= 0.75)
NO2_Lagos_good = [data[data['Flag'] >= 0.75] for data in NO2_Lagos]

NO2_Lagos_good = [data.reset_index(drop=True) for data in NO2_Lagos_good] # need to reset indices everytime after slicing 
    

# only get good data only under calm wind conditions
NO2_Lagos_calm = [data[(abs(data['wind_east']  <= 2)) & 
                       (abs(data['wind_north'] <= 2))] for data in NO2_Lagos_good]

NO2_Lagos_calm = [data.reset_index(drop=True) for data in NO2_Lagos_calm]

# regrid data over Lagos (list grid centres instead of the boundaries)
out_lon = np.arange(-0.60+0.05/2,0.30,0.05)
out_lat = np.arange(51.25+0.05/2,51.75,0.05)
print("out_lon:",out_lon,
      "out_lat:",out_lat,sep="\n")

NO2_Lagos_good_regrid = [regrid_TROPOMI(data) for data in NO2_Lagos_good]
NO2_Lagos_calm_regrid = [regrid_TROPOMI(data) for data in NO2_Lagos_calm]

# save the regridding results

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_regrid')

NO2_Lagos_good_regrid[0].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Lagos_good_regrid[1].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Lagos_good_regrid[2].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Lagos_good_regrid[3].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Lagos_good_regrid[4].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Lagos_good_regrid[5].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Lagos_good_regrid[6].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Lagos_good_regrid[7].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Lagos_good_regrid[8].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Lagos_good_regrid[9].to_csv("NO2_Lagos_regrid_(all winds Flag 0.75)_202006.csv",index=False,sep=',')


NO2_Lagos_calm_regrid[0].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_201909.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[1].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_201910.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[2].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_201911.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[3].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_201912.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[4].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202001.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[5].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202002.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[6].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202003.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[7].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202004.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[8].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202005.csv",index=False,sep=',')
NO2_Lagos_calm_regrid[9].to_csv("NO2_Lagos_regrid_(calm winds Flag 0.75)_202006.csv",index=False,sep=',')

# delete Lagos data to release RAM
del (NO2_Lagos,NO2_Lagos_good,NO2_Lagos_calm,NO2_Lagos_good_regrid,NO2_Lagos_calm_regrid)
