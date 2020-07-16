#####################################################################################################################
#####################################################################################################################
# sample TROPOMI NO2 observations over EU domain 
# regrid to GEOS-Chem nested EU grids (0.25x0.3125)

# do this regriding for each month in seperate Python jobs as it takes a very long time
# test codes in Jupyter Notebook, submit Python jobs using ".py"

import os
import numpy as np
import pandas as pd
import xarray as xr

#####################################################################################################################
# import EU TROPOMI NO2 data extracted from raw TROPOMI L2 swath observations (see "TROPOMI_NO2_1_extract_raw_data.ipynb")
# use 2020-01 as the example

os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_extracted')
NO2_EU_file_202001 = 'NO2_EU_202001.nc'

#####################################################################################################################
# read all EU observations within this month
NO2_EU_extracted_202001 = xr.open_dataset(NO2_EU_file_202001)

# convert to pandas dataframe
NO2_EU_all_202001 = NO2_EU_extracted_202001.to_dataframe()
print("total observations EU (flag >= 0.5):",len(NO2_EU_all_202001))
NO2_EU_all_202001.head()

#####################################################################################################################
# get good data under all wind conditions (Flag >= 0.75)
NO2_EU_good_202001 = NO2_EU_all_202001[NO2_EU_all_202001['Flag'] >= 0.75]

# need to reset indices everytime after slicing the dataframe
NO2_EU_good_202001 = NO2_EU_good_202001.reset_index(drop=True) 

# number of good data under all wind conditions
print("good observations EU (flag >= 0.75) all winds:",NO2_EU_good_202001['date'][0].split('-')[0:2],len(NO2_EU_good_202001))

#####################################################################################################################
# get good data only under calm wind conditions
NO2_EU_calm_202001 = NO2_EU_good_202001[(abs(NO2_EU_good_202001['wind_east']  <= 2)) &
                                        (abs(NO2_EU_good_202001['wind_north'] <= 2))]

# reset indices
NO2_EU_calm_202001 = NO2_EU_calm_202001.reset_index(drop=True)

# number of good data under all wind conditions
print("good observations EU (flag >= 0.75) calm winds:",NO2_EU_calm_202001['date'][0].split('-')[0:2],
      len(NO2_EU_calm_202001),round(len(NO2_EU_calm_202001)/len(NO2_EU_good_202001),2)*100,'%')

#####################################################################################################################
# build a "regridding" funtion

# fisrt decide the target grids centres after regridding
out_lon = np.arange(-15,40+0.3125,0.3125)  # (lon_min,lon_max+resolution,lon_resolution)
out_lat = np.arange(32.75,61.25+0.25,0.25) # (lat_min,lat_max+resolution,lat_resolution)

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
    
# use a small dataset to see if your "regrid_TROPOMI" does what you expect
sample_df = NO2_EU_good_202001.loc[210:212,]
sample_df = sample_df.reset_index(drop=True)
sample_df

sample_regrid = regrid_TROPOMI(sample_df)
sample_regrid

#####################################################################################################################
# regrid NO2 over EU domain (GEOS-Chem nested EU grids)

# first decide the grids centres after regridding
out_lon = np.arange(-15,40+0.3125,0.3125) 
out_lat = np.arange(32.75,61.25+0.25,0.25)

# compare with GEOS-Chem EU nested grids
# http://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_horizontal_grids#0.25_x_0.3125_EU_nested_grid
print(out_lon,out_lat,sep="\n")

NO2_EU_good_regrid_202001 = regrid_TROPOMI(NO2_EU_good_202001)
NO2_EU_calm_regrid_202001 = regrid_TROPOMI(NO2_EU_calm_202001)


#####################################################################################################################
# convert format and save output
os.chdir('/rds/projects/s/shiz-shi-aphh/TROPOMI_DATA/TROPOMI_NO2_regrid')

# good data all winds
NO2_EU_good_regrid_202001 = xr.Dataset.from_dataframe(NO2_EU_good_regrid_202001)
NO2_EU_good_regrid_202001.attrs = {'Data summary':'regrid TROPOMI NO2 in EU under all winds (quality flag >= 0.75)',
                                   'Year month'  : '202001'}
NO2_EU_good_regrid_202001['NO2'].attrs = {'unit':'1e15 molecules_percm2'}
NO2_EU_good_regrid_202001['lon_new'].attrs = {'unit':'Degrees_east'}
NO2_EU_good_regrid_202001['lat_new'].attrs = {'unit':'Degrees_north'}

NO2_EU_good_regrid_202001.to_netcdf('NO2_EU_good_regrid_202001.nc')

# good data calm winds
NO2_EU_calm_regrid_202001 = xr.Dataset.from_dataframe(NO2_EU_calm_regrid_202001)
NO2_EU_calm_regrid_202001.attrs = {'Data summary':'regrid TROPOMI NO2 in EU under calm winds (quality flag >= 0.75)',
                                   'Year month'  : '202001'}
NO2_EU_calm_regrid_202001['NO2'].attrs = {'unit':'1e15 molecules_percm2'}
NO2_EU_calm_regrid_202001['lon_new'].attrs = {'unit':'Degrees_east'}
NO2_EU_calm_regrid_202001['lat_new'].attrs = {'unit':'Degrees_north'}

NO2_EU_calm_regrid_202001.to_netcdf('NO2_EU_calm_regrid_202001.nc')
