###############################################################################################
###############################################################################################
# Here I am trying to plot any oversampled TROPOMI data over any domain using one single script (still under developing)
# Overampled TROPOMI data can be: NO2/HCHO/SO2/CO/O3 VCD, the corresponding VCD uncertainties and number of oversampled pixels at any resolution

# You can jump to the functions:
# Line 30-36 function 1> read_oversampled_TROPOMI(): read oversampled TROPOMI data from the txt file and convert to pandas dataframe
# Line 59-67 function 2> expand_grid(): given the range of lat and lon, return the full cominations of all lats and lons
# Line 195-261 function 3> plot_TROPOMI(): input a few parameters (e.g. domain,species, and time), return the corresponding plots
# Apart from these functions, some other steps are still needed along the process

# Here I use oversampled TROPOMI NO2 over Africa at the 0.25x0.25 resolution to establish the functions. 
# Line 30-278: establish the functions
# Line 280-315: example applications of the functions

import os
import glob
import pandas as pd
import numpy as np
import xarray as xr

import matplotlib.pyplot as plt
from gamap_colormap import WhGrYlRd
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
###############################################################################################
# move to the directory for oversampled output txt files
os.chdir("/rds/projects/s/shiz-shi-aphh/TROPOMI_NO2_2_Oversampling/Output")

# build a function to convert oversampled TROPOMI txt file into pandas dataframe
def read_oversampled_TROPOMI(TROPOMI_oversampled_NO2_file):
    '''read the txt file from the TROPOMI oversampling algorithm'''
    df = pd.read_csv(TROPOMI_oversampled_NO2_file,sep="\s+",header=None)
    df = df.iloc[:,2:7]
    df.columns = ['lat','lon','NO2','Count','NO2_uncertainty']
    return df

# explore the data and try plotting using the coarser resolution
NO2_AF_025_files = sorted(glob.glob('TROPOMI_oversampling_output_AF_NO2_0.75_*_0.25'))
NO2_AF_025 = [read_oversampled_TROPOMI(file) for file in NO2_AF_025_files]
###############################################################################################
# check a sample data frame
print(NO2_AF_025[0])

# the number of data points are not always fixed, at some grids there is no data
# this will surely be the case for even finer resolutions
for i in range(len(NO2_AF_025)):
    print("Number of rows:",NO2_AF_025[i].shape[0])   
    
# create the full list of lats and lons of grids centres within the AF domain at 0.25 resolution
# before oversampling, the boundaries for Africa are chosen as (lat[-40,48], lon[-20,50])
AF_lat_025 = np.arange(-40-0.25/2,48+0.25/2+0.25, 0.25,dtype=None) # use "(min,max+resolution,resolution)" to keep the max
AF_lon_025 = np.arange(-20-0.25/2,50+0.25/2+0.25, 0.25,dtype=None) 

print("Total number of lats on map should be:",len(AF_lat_025))
print("Total number of lons on map should be:",len(AF_lon_025))
print("Total number of data points in 'AF_0.25x0.25' should be:", len(AF_lat_025)*len(AF_lon_025))
###############################################################################################
# build a function to list all combinations of lats and lons in the domain (equivalent to "expand.grid" in R)
def expand_grid(lat,lon):
    '''list all combinations of lon and lat using expand_grid(lon,lat)'''
    test = [(A,B) for A in lat for B in lon]
    test = np.array(test)
    test_lat = test[:,0]
    test_lon = test[:,1]
    full_grid = pd.DataFrame({'lat': test_lat, 'lon': test_lon})
    return full_grid
 
# create the full 0.25x0.25 grids for Africa domain
AF_025_grid = expand_grid(AF_lat_025,AF_lon_025)
print(AF_025_grid)

# outer join the full grid and the oversampled data
# so the domain coverged is always consistent
# later re-think about this joining method, for now focus on plotting
NO2_AF_full_025 = [pd.merge(AF_025_grid,data,how='outer', on=['lat','lon']) for data in NO2_AF_025]

# check a sample result
print(NO2_AF_full_025[0])

# after the previous steps, both "lat" and "lon" should be sorted ascending already, but just to make sure
NO2_AF_full_025 = [df.sort_values(by=['lat','lon'],ascending=[True, True]) for df in NO2_AF_full_025]

# check a sample result
# so lat is repeated and lon is cycled
print(NO2_AF_full_025[0])

# since the lats and lons are created already
print(len(AF_lat_025),len(AF_lon_025))

# now just reshape the variable to match the map dimension
NO2_025 = [df['NO2'].values.reshape((354,282)) for df in NO2_AF_full_025]
NO2_count_025 = [df['Count'].values.reshape((354,282)) for df in NO2_AF_full_025]
NO2_uncertainty_025 = [df['NO2_uncertainty'].values.reshape((354,282)) for df in NO2_AF_full_025]

# convert the "lat","lon","NO2","count" and "Uncertainty" into xarray data array format
# this is because xarray data can surely be plotted by cartopy
# the TROPOMI data is now ready for plotting
NO2_xarray_025 = [xr.DataArray(data, coords=[('lat', AF_lat_025),('lon', AF_lon_025)]) for data in NO2_025]
NO2_count_xarray_025 = [xr.DataArray(data, coords=[('lat', AF_lat_025),('lon', AF_lon_025)]) for data in NO2_count_025]
NO2_uncertainty_xarray_025 = [xr.DataArray(data, coords=[('lat', AF_lat_025),('lon', AF_lon_025)]) for data in NO2_uncertainty_025]
###############################################################################################
# import shapefiles
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_1_AF_SEAPORT/AF_shapefiles")
AF_shapefile = r'Africa.shp'
AF_map = ShapelyFeature(Reader(AF_shapefile).geometries(),ccrs.PlateCarree(), edgecolor='black',facecolor='none')
###############################################################################################
# read locations of seaports and sort by lat
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_1_AF_SEAPORT/AF_surface')
AF_seaports = pd.read_csv("AF_seaports_coordinates.csv")
AF_seaports = AF_seaports.sort_values(by='lat').reset_index(drop=True)
AF_seaports
###############################################################################################
# move to directory for output plots
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/TROPOMI/project_1_AF_SEAPORT/')

# make a single plot
fig = plt.figure(figsize=[10,10])
ax= plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([-20, 50, -40, 48], crs=ccrs.PlateCarree()) # zoom in [lon,lon,lat,lat]
ax.add_feature(AF_map) # add country boundaries in Africa
ax.coastlines() # add coastlines for other domains

# plot the NO2 concentration map
NO2_xarray_025[0].plot(ax=ax,cmap='jet',vmax=2,vmin=0,
                       cbar_kwargs={'shrink': 0.5, 
                                    'pad' : 0.05,
                                    'orientation':'horizontal',
                                    'label':'NO$_2$ tropospheric column [$10^{15}$ molec. cm$^{-2}$]',
                                    'ticks':(0,0.5,1,1.5,2)})
                                    
# mark locations of seaports
ax.scatter(x=AF_seaports['lon'], y=AF_seaports['lat'],facecolors='none',edgecolors='black',linewidths=2,s =100)

# set color bar label size
plt.rcParams.update({'font.size':12})

# provide a "transparent grid line" so we can mark the lat and lon using labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
gl.xlines = False
gl.ylines = False
gl.xlabels_top = False
gl.ylabels_left = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER
gl.xlabel_style = {'size': 15, 'color': 'black'}
gl.ylabel_style = {'size': 15, 'color': 'black'}

# add text based on coordinates
# ax.annotate('2019-08', xy=(0,0), xytext=(-20,-39),fontsize = 20,color='white')

# add text regardless of coordinates
plt.text(0.01, 0.01,'2019-08',size = 20,color='white',transform = ax.transAxes)

# add title
ax.set_title('TROPOMI NO$_2$ over Africa 0.25x0.25', size = 15)

# display the plot
plt.show()

# save out with reduced figure margins
fig.savefig('Africa_NO2_025_201908_all_boundaries_test.png', dpi=300,bbox_inches='tight')
###############################################################################################
# before making multiple plots, first have a rough idea of the ranges of each varialbe
# add more statistics to understand the coverage, like how many percent of grids are covered.

# NO2 VCD
print('#'*50,"NO2 VCD",sep='\n')

for i in range(len(NO2_xarray_025)):
    print("min:",round(np.nanmean(NO2_xarray_025[i]),1),
          "max:",round(np.nanmax(NO2_xarray_025[i]),1),
          "mean:",round(np.nanmean(NO2_xarray_025[i]),1),
          "median:",round(np.nanmedian(NO2_xarray_025[i]),1))

# NO2 VCD uncertainty
print('#'*50,"NO2 VCD uncertainty",sep='\n')

for i in range(len(NO2_uncertainty_xarray_025)):
    print("min:",round(np.nanmin(NO2_uncertainty_xarray_025[i]),1),
          "max:",round(np.nanmax(NO2_uncertainty_xarray_025[i]),1),
          "mean:",round(np.nanmean(NO2_uncertainty_xarray_025[i]),1),
          "median:",round(np.nanmedian(NO2_uncertainty_xarray_025[i]),1))

# Number of oversampled pixels
print('#'*50,"Number of oversampled pixels",sep='\n')

for i in range(len(NO2_count_xarray_025)):
    print("min:",round(np.nanmin(NO2_count_xarray_025[i]),1),
          "max:",round(np.nanmax(NO2_count_xarray_025[i]),1),
          "mean:",round(np.nanmean(NO2_count_xarray_025[i]),1),
          "median:",round(np.nanmedian(NO2_count_xarray_025[i]),1))
###############################################################################################
# build a function to plot TROPOMI data over any domain in the world at any resolution and during any sampling period
# the speceis can be any TROPOMI species (NO2, HCHO, CO, O3 and SO2), but the codes need to be slight edited
# the variable can be VCD, VCD uncertainty or number of oversampled pixels

def plot_TROPOMI(domain_name,domain,resolution,time,var_name,var_xarray,var_min,var_max):    
    '''Iput the following info about one variable, return the corresponding colored concentration map(s). 
       1> domain_name: this string will appear at the start of the output file name
       2> domain: define the output figure domain [lon_min,lon_max,lat_min,lat_max]
       3> resolution: this string will appear in the output file name
       4> time: this string will appear at the bottom-left corner of the figure and in the outputfile name
       5> var_name: the name of the input species (Can be "NO2","NO2_uncertainty" or "NO2_count". Likewise for other species.)
       6> var_xarray: the xarray data array of the input species
       7> var_min: the min on the color bar
       8> var_max: the max on the color bar
    '''
    # firt go through some figure settings based on input species
    if var_name == "NO2":
        var_label = 'NO$_2$ tropospheric column [$10^{15}$ molec. cm$^{-2}$]'
        var_ticks = (0,0.5,1,1.5,2)
    elif var_name == "NO2_uncertainty":
        var_label = 'NO$_2$ tropospheric column uncertainty [$10^{15}$ molec. cm$^{-2}$]'
        var_ticks = (0,0.01,0.02,0.03,0.04,0.05)
    elif var_name == "NO2_count":
        var_label = 'Number of oversampled pixels'
        var_ticks = (0,200,400,600,800)
    else:
        print("Invalid variable name! Available names are: NO2,NO2_uncertainty,NO2_count")
    
    cbar_keys = {'shrink': 0.5, 'pad' : 0.05,'orientation':'horizontal','label':var_label,'ticks':var_ticks}
    
    # start plotting
    fig = plt.figure(figsize=[10,10])
    ax= plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent(domain, crs=ccrs.PlateCarree()) # zoom in to the chosen domain [lon,lon,lat,lat]
    ax.add_feature(AF_map) # add shapefile (here we are using country boundaries in Africa)
    ax.coastlines() # add coastlines for other domains
    
    # color the map based on values of the chosen variable
    var_xarray.plot(ax=ax,cmap='jet',vmin=var_min,vmax=var_max,cbar_kwargs=cbar_keys)

    # set color bar label size
    plt.rcParams.update({'font.size':12})
    
    # mark locations of surface sites (Here we are using African seaports)
    # you can also add colored circles to represent surface data values
    ax.scatter(x=AF_seaports['lon'], y=AF_seaports['lat'],facecolors='none',edgecolors='black',linewidths=2,s =100)

    # provide a "transparent grid line" so we can mark the lat and lon using labels
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True)
    gl.xlines = False
    gl.ylines = False
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 15, 'color': 'black'}
    gl.ylabel_style = {'size': 15, 'color': 'black'}

    # add "time" input at the bottom-left corner
    plt.text(0.01, 0.05,domain_name,size = 20,color='white',weight='bold',transform = ax.transAxes)
    plt.text(0.01, 0.01,time,size = 20,color='white',weight='bold',transform = ax.transAxes)
    
    # save out with reduced figure margins
    fig.savefig(domain_name+'_'+var_name+'_'+resolution+'_'+time+'.png', dpi=300,bbox_inches='tight')
    
# next version: think about how to store and access a list of domain names and domains (use "re"?) 
###############################################################################################
# now producce a map for all periods + daomains/locations

# function usage: plot_TROPOMI(domain_name,domain,resolution,time,var_name,var_xarray,var_min,var_max)

# create a full list of sampling months
months = ['2019-08','2019-09','2019-10','2019-11','2019-12',
          '2020-01','2020-02','2020-03','2020-04','2020-05','2020-06','2020-07']

# plot over whole Africa
for i in range(len(months)):
    plot_TROPOMI('Africa',[-20, 50, -40, 48],'0.25x0.25',months[i],'NO2',NO2_xarray_025[i],0,2)
    plot_TROPOMI('Africa',[-20, 50, -40, 48],'0.25x0.25',months[i],'NO2_uncertainty',NO2_uncertainty_xarray_025[i],0,0.05)
    plot_TROPOMI('Africa',[-20, 50, -40, 48],'0.25x0.25',months[i],'NO2_count',NO2_count_xarray_025[i],0,800)
###############################################################################################
# plot over each seaport

# check the seaport list again
AF_seaports

# for now, use "CapeTown","Mombassa","Lagos","Suez"
# later write codes to use subplots for each seaport
# https://github.com/geoschem/GEOSChem-python-tutorial/blob/main/Chapter01_NetCDF_xarray.ipynb
# https://stackoverflow.com/questions/40996175/loading-a-rds-file-in-pandas
    
# CapeTown
for i in range(len(months)):
    plot_TROPOMI('CapeTown',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2',NO2_xarray_025[i],0,2)
    plot_TROPOMI('CapeTown',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2_uncertainty',NO2_uncertainty_xarray_025[i],0,0.05)
    plot_TROPOMI('CapeTown',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2_count',NO2_count_xarray_025[i],0,800)
    
# Mombassa
for i in range(len(months)):
    plot_TROPOMI('Mombassa',[39.65-0.25, 39.65+0.25, -4.05-0.25, -4.05+0.25],'0.25x0.25',months[i],'NO2',NO2_xarray_025[i],0,1)
    plot_TROPOMI('Mombassa',[39.65-0.25, 39.65+0.25, -4.05-0.25, -4.05+0.25],'0.25x0.25',months[i],'NO2_uncertainty',NO2_uncertainty_xarray_025[i],0,0.05)
    plot_TROPOMI('Mombassa',[39.65-0.25, 39.65+0.25, -4.05-0.25, -4.05+0.25],'0.25x0.25',months[i],'NO2_count',NO2_count_xarray_025[i],0,800)
    
# Lagos
for i in range(len(months)):
    plot_TROPOMI('Lagos',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2',NO2_xarray_025[i],0,1.5)
    plot_TROPOMI('Lagos',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2_uncertainty',NO2_uncertainty_xarray_025[i],0,0.05)
    plot_TROPOMI('Lagos',[18.45-0.25, 18.45+0.25, -33.91-0.25, -33.91+0.25],'0.25x0.25',months[i],'NO2_count',NO2_count_xarray_025[i],0,800)
   
# Suez
for i in range(len(months)):
    plot_TROPOMI('Suez',[32.37-0.25, 32.37+0.25, 31.21-0.25, 31.21+0.25],'0.25x0.25',months[i],'NO2',NO2_xarray_025[i],0,2)
    plot_TROPOMI('Suez',[32.37-0.25, 32.37+0.25, 31.21-0.25, 31.21+0.25],'0.25x0.25',months[i],'NO2_uncertainty',NO2_uncertainty_xarray_025[i],0,0.05)
    plot_TROPOMI('Suez',[32.37-0.25, 32.37+0.25, 31.21-0.25, 31.21+0.25],'0.25x0.25',months[i],'NO2_count',NO2_count_xarray_025[i],0,800)
 
# End
###############################################################################################
# This plotting routine needs to be smarter:
# 1> When plotting over 20-50 cities or seaports over 12 months, the pollution levels vary, so the vmax and vmin for the color bar will always change! 
#    Is there a way to quickly decide the vmax and vmin for the color bar for each set of plots (one set = one speices over one domain during one period)?
# 2> Instead of plotting each location seperately, can I use sth like a list to loop over all locations?
# 3> Is there a package that allows you to have more controls of the figure? Is matplotlib enough? I am thinking about adding/editing all features grammatically (like 'ggplot2' in R).
#    I feel that advanced users prefer things like "imshow"?
# 4> Try subsetting the 0.01x0.01 data over the target domain first, then plot. I guess it will be faster than zooming in everytime.
# 5> Need to use subplots to reduce the number of output plots
