#####################################################################################################################
#####################################################################################################################
# See how TROPOMI NO2 responds to the Suez Canal blockage
# When downloading the data, look at a larger domain (Suez and its surrounding + Mediterranean Sea)

import os
import glob
import numpy  as np
import pandas as pd
from netCDF4 import Dataset
import xarray as xr

''' 
Note on this Suez Canal blockage

Blockage period: 23-29 March 2021
Data download period: 5 January - 26 April 2021
Domain (lon_min,lat_min,lon_max,lat_max): -20,5,60,50
Corresponding hour windows for data donwload: [6,7,8,9,10,11,12,13,14]

First test: sample weekly data before, during and after the blockage, get maps and time serires plot
Second test: get daily maps and combine with GeoViews
'''

#####################################################################################################################
# build a function to read oversampled TROPOMI NO2 as pandas dataframes
def read_oversampled_NO2(TROPOMI_oversampled_NO2_output_file):
    '''read the output file for oversampled TROPOMI NO2'''
    df = pd.read_csv(TROPOMI_oversampled_NO2_output_file,sep="\s+",header=None)
    df = df.iloc[:,2:7]
    df.columns = ['lat','lon','NO2','Count','NO2_uncertainty']
    return df

#####################################################################################################################
# the spatial coverage may not be consistent on different days or during different weeks

# read all the data from the weekly results
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Oversample_output')
Oversampled_NO2_files = sorted(glob.glob("Oversample_output_Suez_NO2_week*"), key=lambda x: int(x.split("_")[-2]))
print(*Oversampled_NO2_files,sep="\n")
oversampled_data = [read_oversampled_NO2(file) for file in Oversampled_NO2_files]

# use all the data ever sampled to decide the max dimension
lat_min = []
lat_max = []
lon_min = []
lon_max = []

for i in range(len(oversampled_data)):
    lat_min.append(oversampled_data[i].lat.min())
    lat_max.append(oversampled_data[i].lat.max())
    lon_min.append(oversampled_data[i].lon.min())
    lon_max.append(oversampled_data[i].lon.max())
    
lat_min = min(lat_min)
lat_max = max(lat_max)
lon_min = min(lon_min)
lon_max = max(lon_max)

# check the full dimension
print("lat_min:",lat_min)
print("lat_max:",lat_max)
print("lon_min:",lon_min)
print("lon_max:",lon_max)

# With the dimension above and the resolution, we can create a consistent domain ("the full grid")
# so that we can combine the data from different days/weeks together

# first list all the lats and lons: use (min,max+1/2 resolutions, resolution) to keep the max value in Python
# just round the floats created by Python to be safe
# as the "pd.merge" step later will require the values of "keys" to be excatly the same
Res = 0.05
domain_lat = np.arange(lat_min,lat_max+Res/2,Res,dtype=None)
domain_lon = np.arange(lon_min,lon_max+Res/2,Res,dtype=None) 

domain_lat = np.round(domain_lat,3)
domain_lon = np.round(domain_lon,3)

# build a function to create a "full grid" by listing the full combinations of lats and lons in the domain
def expand_grid(lat,lon):
    '''list all combinations of lats and lons using expand_grid(lat,lon)'''
    test = [(A,B) for A in lat for B in lon]
    test = np.array(test)
    test_lat = test[:,0]
    test_lon = test[:,1]
    full_grid = pd.DataFrame({'lat': test_lat, 'lon': test_lon})
    return full_grid

# create the "full grid"
domain_grid = expand_grid(domain_lat,domain_lon)
print(domain_grid)

################################################################################################
# Now we can read each single dataset and match it with the full grid

# Step 1> select the oversampled data
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Oversample_output')

# change input time to read daily data or weekly data
time = 'week_1'

Oversampled_NO2_file  = "Oversample_output_Suez_NO2_"+str(time)+"_0.05"

# check the selected data
print(Oversampled_NO2_file)

# Step 2> feed the oversampled data into this data cleaning routine

# read oversampled NO2 data
NO2_data  = read_oversampled_NO2(Oversampled_NO2_file)

# combine the data with the full domain grids
NO2_data = pd.merge(domain_grid,NO2_data,how='left', on=['lat','lon'])
NO2_data = NO2_data.sort_values(by=['lat','lon'], ascending=[True, True])

# reshape the variables from 1D in the dataframe to the map dimension
NO2 = NO2_data['NO2'].values.reshape(len(domain_lat),len(domain_lon))
NO2_uncertainty = NO2_data['NO2_uncertainty'].values.reshape(len(domain_lat),len(domain_lon))
Count = NO2_data['Count'].values.reshape(len(domain_lat),len(domain_lon))

# convert to xarray for plotting
NO2_xarray = xr.DataArray(NO2, coords=[('lat', domain_lat),('lon', domain_lon)])
NO2_uncertainty_xarray = xr.DataArray(NO2_uncertainty, coords=[('lat', domain_lat),('lon', domain_lon)])
Count_xarray = xr.DataArray(Count, coords=[('lat', domain_lat),('lon', domain_lon)]) 

# but it is complicated to save out the results one by one for multiple days or weeks
################################################################################################
################################################################################################
# So here we use the list comprehensions to process multiple files

#################
# weekly data
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Oversample_output')

# select the files and sort them numerically
Oversampled_NO2_files_weekly = sorted(glob.glob("Oversample_output_Suez_NO2_week*"), key=lambda x: int(x.split("_")[-2]))
print(*Oversampled_NO2_files_weekly,sep="\n")

# read oversampled data and match with the "full grid"
Oversampled_NO2_week = [read_oversampled_NO2(file) for file in Oversampled_NO2_files_weekly]
Oversampled_NO2_week = [pd.merge(domain_grid,data,how='left', on=['lat','lon']) for data in Oversampled_NO2_week]
Oversampled_NO2_week = [data.sort_values(by=['lat','lon'], ascending=[True, True]) for data in Oversampled_NO2_week]

# convert the data to the xarray format for plotting
NO2_week = [data['NO2'].values.reshape(len(domain_lat),len(domain_lon)) for data in Oversampled_NO2_week]
NO2_week_xr = [xr.DataArray(data, coords=[('lat', domain_lat),('lon', domain_lon)]) for data in NO2_week]

#################
# daily data
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Oversample_output')

# select the files and sort them numerically
Oversampled_NO2_files_daily = sorted(glob.glob("Oversample_output_Suez_NO2_day*"), key=lambda x: int(x.split("_")[-2]))
print(*Oversampled_NO2_files_daily,sep="\n")

# read oversampled data and match with the "full grid"
Oversampled_NO2_day = [read_oversampled_NO2(file) for file in Oversampled_NO2_files_daily]
Oversampled_NO2_day = [pd.merge(domain_grid,data,how='left', on=['lat','lon']) for data in Oversampled_NO2_day]
Oversampled_NO2_day = [data.sort_values(by=['lat','lon'], ascending=[True, True]) for data in Oversampled_NO2_day]

# convert the data to the xarray format for plotting
NO2_day = [data['NO2'].values.reshape(len(domain_lat),len(domain_lon)) for data in Oversampled_NO2_day]
NO2_day_xr = [xr.DataArray(data, coords=[('lat', domain_lat),('lon', domain_lon)]) for data in NO2_day]

################################################################################################
# Start making maps to have a quick look at the results 

# avoid setting "%matplotlib inline" as it is time consuming when we need to produce many figures
import matplotlib.pyplot as plt
import cartopy.crs as crs
import geopandas as gpd

# read shape file (Global high resolution shoreline database from NOAA: https://www.ngdc.noaa.gov/mgg/shorelines/)
# use "full reolution" here to avoid misrepresentation of land and water
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/shapefiles/gshhg-shp-2.3.7/GSHHS_shp/f")
world_shore = gpd.read_file("GSHHS_f_L1.shp")

################################################################################################
# build a function to quickly generate maps without a legend to save space on a slide

def quick_plot(input_xr,plot_domain,var_min,var_max,output_figure_name):
    '''
    Input a xarray data array, define the map domain, provide the min and max of the values on map. Provide a outputfile name.
    '''
    # set the figure size, the aspect ratio is set to be 2:1 due to the sampling region
    fig = plt.figure(figsize=[20,10])
    
    # set the map projection and domain: https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html#cartopy-projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(plot_domain)
    
    # plot the value on map
    im = input_xr.plot(ax=ax,cmap='jet',vmin=var_min,vmax=var_max)
    
    # add shapefile
    ax.add_geometries(world_shore.geometry, crs=ccrs.PlateCarree(),edgecolor='black',facecolor='none')
    
    # remove the colorbar and tile
    plt.delaxes(fig.axes[1])
    ax.set_title('')

    # save out
    fig.savefig(output_figure_name, dpi=100,bbox_inches='tight')

    # close the figure to avoid taking CPU memory
    plt.close()
    
################################################################################################
# build a function to generatet the bar for the figures above

def plot_color_bar(input_xr,plot_domain,label,var_min,var_max,output_figure_name):
    '''
    Draw the figure in the same way as above, but remove the plot rather than the colorbar.
    '''
    fig = plt.figure(figsize=[20,10])
    cbar_keys = {'shrink': 1, 'pad' : 0.05,'orientation':'horizontal','label':label} 
    
    # set the map projection: https://scitools.org.uk/cartopy/docs/v0.15/crs/projections.html#cartopy-projection
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(plot_domain)
    
    # plotthe value on map
    im = input_xr.plot(ax=ax,cmap='jet',cbar_kwargs=cbar_keys,vmin=var_min,vmax=var_max)
    
    # set color bar label size
    plt.rcParams.update({'font.size':25})
    ax.xaxis.label.set_size(25)
    
    # remove the plot
    plt.delaxes(fig.axes[0])
    
    # save out
    fig.savefig(output_figure_name, dpi=100,bbox_inches='tight')
    
    # close the figure to avoid taking CPU memory
    plt.close()

################################################################################################
# check again the data for plotting

print("weekly data:",len(NO2_week_xr))
print("daily data:",len(NO2_day_xr))

# generate corresponding output file names

# weekly maps
Suez_weeks = list(range(1,17))
Suez_weeks = [str('Suez_NO2_map_week_') + str(week_number) for week_number in Suez_weeks]
print(*Suez_weeks,sep="\n")

# daily maps
Suez_days = list(range(1,29))
Suez_days = [str('Suez_NO2_map_day_') + str(date_number) for date_number in Suez_days]
print(*Suez_days,sep="\n")

################################################################################################
# output multiple plots together

os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Figures')

# maps during the blockage
# week 12
# day 8-14

# plot weekly data

# plot over the big domain [lon_min,lon_max,lat_min,lat_max]
Suez_domain_big = [-20,60,5,50]

for i in range(len(NO2_week_xr)):
    quick_plot(NO2_week_xr[i],Suez_domain_big,0,2,Suez_weeks[i]+str('_big'))

# plot over the small domain [lon_min,lon_max,lat_min,lat_max]
Suez_domain_small = [26,60,10,35]

for i in range(len(NO2_week_xr)):
    quick_plot(NO2_week_xr[i],Suez_domain_small,0,2,Suez_weeks[i]+str('_small'))

# generate the color bar at the end
plot_color_bar(NO2_week_xr[0],Suez_domain_small,'NO$_2$ tropospheric column [$10^{15}$ molec. cm$^{-2}$]',0,2,"Suez_NO2_color_bar")

# plot daily data

# plot over the small domain [lon_min,lon_max,lat_min,lat_max]
Suez_domain_small = [26,60,10,35]

for i in range(len(NO2_day_xr)):
    quick_plot(NO2_day_xr[i],Suez_domain_small,0,2,Suez_days[i]+str('_small'))

################################################################################################
################################################################################################
# Use GeoViews to combine the maps together in time series

# load GeoViews package
import geoviews as gv
import geoviews.feature as gf
import cartopy.crs as crs

# it is important to check your geoviews version, some commands may not work in a wrong version
# this script is written under version 1.9.1
print(gv.__version__)

# there are two backends ('bokeh', 'matplotlib') for the GeoViews
# later we will use "bokeh" for interactive plots

################################################################################################
# weekly maps

# list all the weeks
Suez_weeks = ['01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16']
print(*Suez_weeks,sep="\n")

# combine the xarray data arrays from weekly results

# make a copy first
weekly_data = NO2_week_xr.copy()

# add the variable name
weekly_data = [data.rename('NO2') for data in weekly_data]

# add a time dimension to the data
for i in range(len(NO2_week_xr)):
    NO2_week_xr[i] = NO2_week_xr[i].assign_coords(week=Suez_weeks[i])
    NO2_week_xr[i] = NO2_week_xr[i].expand_dims('week')
    
# combine the data together
NO2_week_xr_combined = xr.concat(NO2_week_xr,'week')

# you can zoom in and change maps, so normally there is no need to make a small map
# but if you have to reduce the file size, you can subset over the small domain
# weekly_data = [data.sel(lat=slice(10,35),lon = slice(26,60)) for data in weekly_data]

# check the results
NO2_week_xr_combined

# output the plots

# first move to the output directory
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Figures')

# turn on "bokeh" backend to enable interactive map
gv.extension('bokeh')

# extract data from the combined xarray
gv_data  = gv.Dataset(NO2_week_xr_combined,['lon','lat','week'],'NO2',crs=crs.PlateCarree())

# use the data to generate the geoviews image
gv_image = gv_data.to(gv.Image)

# decide features of the output figure
gv_image_out = gv_image.opts(cmap='jet', clim=(0,2), colorbar=True, width=800, height=500) * gf.coastline 

# save out the interactive map
renderer = gv.renderer('bokeh')
renderer.save(gv_image_out, 'weekly_maps')

################################################################################################
# daily maps

# list all the dates

def list_dates_between(start_date,end_date):
    '''Select TROPOMI files within the start date ('yyyymmdd') and end date ('yyyymmdd')'''
    # list all the dates between the start and the end
    from datetime import date, timedelta
    start_date = date(int(start_date[0:4]),int(start_date[4:6]),int(start_date[6:8]))
    end_date   = date(int(end_date[0:4]),int(end_date[4:6]),int(end_date[6:8]))
    delta = end_date - start_date 
    sampling_dates = []
    for i in range(delta.days + 1):
        sampling_dates.append((start_date + timedelta(days=i)).strftime('%Y%m%d'))
    
    # print out all the sampling dates
    return sampling_dates

# list all the dates
Suez_days = list_dates_between("20210316","20210412")
print("number of days:",len(Suez_days))
print(*Suez_days,sep="\n")

# combine the xarray data arrays from daily results

# make a copy first
daily_data = NO2_day_xr.copy()

# add the variable name
daily_data = [data.rename('NO2') for data in daily_data]

# add a time dimension to the data
for i in range(len(NO2_day_xr)):
    NO2_day_xr[i] = NO2_day_xr[i].assign_coords(date=Suez_days[i])
    NO2_day_xr[i] = NO2_day_xr[i].expand_dims('date')
    
# combine the data together
NO2_day_xr_combined = xr.concat(NO2_day_xr,'date')

# check the results
NO2_day_xr_combined

# output the plots

# first move to the output directory
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_2_Suez_Canal/Figures')

# turn on "bokeh" backend to enable interactive map
gv.extension('bokeh')

# extract data from the combined xarray
gv_data  = gv.Dataset(NO2_day_xr_combined,['lon','lat','date'],'NO2',crs=crs.PlateCarree())

# use the data to generate the geoviews image
gv_image = gv_data.to(gv.Image)

# decide features of the output figure
gv_image_out = gv_image.opts(cmap='jet', clim=(0,2), colorbar=True, width=800, height=500) * gf.coastline 

# save out the interactive map
renderer = gv.renderer('bokeh')
renderer.save(gv_image_out, 'daily_maps')

# For now, the default coastline from GeoViews is used
# If you can crop and create your own shapefile, you should be able to use high resolution shorelines from NOAA
# Think about how to do this with geopandas
#####################################################################################################################
#####################################################################################################################