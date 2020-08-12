##################################################################################
##################################################################################
# Read daily surface VOCs from GEOS-Chem 12.0.0
# Save variables of interest into one single netcdf file

import os
import glob
import xarray as xr

# load daily model output files during 2016-11 (with default MEIC emissions)
os.chdir("/rds/projects/2018/maraisea-glu-01/RDS/GEOSChem/GEOS-Chem_12.0.0/Run_Directory/Beijing_VOCs/merra2_05x0625_tropchem_as_2016_November_default_MEIC/APHH_VOCs_daily_outputs")

# GEOSChem.SpeciesConc are in the unit of "mol/mol"
# need to use the number density ("Met_AIRNUMDEN") from "GEOSChem.StateMet" to derive mass concentrations
Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))
StateMet = sorted(glob.glob("GEOSChem.StateMet*.nc4"))
print("number of files:",len(Species),Species[0],Species[-1],sep="  ")
print("number of files:",len(StateMet),StateMet[0],StateMet[-1],sep="  ")

# open daily netcdf files
GC_species  = [xr.open_dataset(file) for file in Species]
GC_statemet = [xr.open_dataset(file) for file in StateMet]
##################################################################################
# Extract species of interest: Benzene/Toluene/Ethane/Acetylene/Acetone/Formaldehyde/Propane/Xylenes/nbutane/ibutane/Ethylene/Ethanol
# Full GEOS-Chem species info (shortname,formula,molec weight, phase...) at http://wiki.seas.harvard.edu/geos-chem/index.php/Species_in_GEOS-Chem
# C2H2 is listed in the model wiki, but not found in the output file, maybe it is added in the newer model versions?

# "time=0" removes the uneccesary dimension, "lev=0" selects surface data only

surface_BENZ = [data['SpeciesConc_BENZ'].isel(time=0,lev=0) for data in GC_species] # Benz (C6H6)
surface_TOLU = [data['SpeciesConc_TOLU'].isel(time=0,lev=0) for data in GC_species] # Toluene (C7H8)
surface_C2H6 = [data['SpeciesConc_C2H6'].isel(time=0,lev=0) for data in GC_species] # Ethane (C2H6)
surface_ACET = [data['SpeciesConc_ACET'].isel(time=0,lev=0) for data in GC_species] # Acetone (CH3C(O)CH3) 
surface_CH2O = [data['SpeciesConc_CH2O'].isel(time=0,lev=0) for data in GC_species] # Formaldehyde (CH2O,HCHO)
surface_C3H8 = [data['SpeciesConc_C3H8'].isel(time=0,lev=0) for data in GC_species] # Propane (C3H8)
surface_EOH  = [data['SpeciesConc_EOH'].isel(time=0,lev=0)  for data in GC_species] # Ethanol (C2H5OH)
##################################################################################
# read air number density from "GEOSChem.StateMet" and fix the unit
# Important note: the unit of air number density in GEOS-Chem output is wrong!
# It should be cm-3, instead of m-3
# Ways to prove this: 1> if it is m-3, then values are not of the same magnitude with those derived from ideal gas law
#                     2> if it is m-3, the mass concentrations in the end will be around 10^6 times lower than what it should be 
# The "Met_AIRNUMDEN" is also not included in the "StateMet" diagnostics on GEOS-Chem wiki
# In order to output "Met_AIRNUMDEN", you need to manually add "Met_AIRNUMDEN" to "History.rc" when setting up the model

# extract air number density
surface_airnumberdensity = [data['Met_AIRNUMDEN'].isel(time=0,lev=0) for data in GC_statemet] 

# fix the unit
for i in range(len(surface_airnumberdensity)):
    surface_airnumberdensity[i].attrs['long_name'] = 'Dry air number density'
    surface_airnumberdensity[i].attrs['units'] = 'cm-3'
##################################################################################    
# combine relevant variables
China_daily_VOCs = []

for i in range(len(GC_species)):
    China_daily_VOCs.append(xr.merge([surface_BENZ[i],
                                      surface_TOLU[i],
                                      surface_C2H6[i],
                                      surface_ACET[i],
                                      surface_CH2O[i],
                                      surface_C3H8[i],
                                      surface_EOH[i],
                                      surface_airnumberdensity[i]]))
                                      
# add global attributes from raw GEOS-Chem output file
for i in range(len(China_daily_VOCs)):
    China_daily_VOCs[i].attrs = GC_species[i].attrs
    China_daily_VOCs[i].attrs['title'] = 'Daily GEOS-Chem surface VOCs concentrations in China'
    
# now the info needed is combined on each day
print(China_daily_VOCs[0])

# combine data from multiple days
China_VOCs = xr.concat(China_daily_VOCs,'time')
print(China_VOCs)

# save the outputs to netcdf
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/geoschem")
China_VOCs.to_netcdf('China_daily_surface_VOCs_geoschem.nc')

# check the netcdf file generated
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/geoschem")
China_daily_surface_VOCs = xr.open_dataset('China_daily_surface_VOCs_geoschem.nc')
print(test)

# plot on map for a sanity check
%matplotlib inline
import matplotlib.pyplot as plt
from gamap_colormap import WhGrYlRd 
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/domain/CHN_ADM')
China_shape = r'CHN_ADM1.shp'
china_map = ShapelyFeature(Reader(China_shape).geometries(),ccrs.PlateCarree(), edgecolor='black',facecolor='none')

t= plt.axes(projection=ccrs.PlateCarree())
t.add_feature(china_map)
t.set_extent([100, 125, 30, 45], crs=ccrs.PlateCarree()) # [lon_min,lon_max,lat_min,lat_max]

VOC_test1 = China_daily_surface_VOCs['SpeciesConc_BENZ'].isel(time=0)
VOC_test1.plot(ax=t,cmap=WhGrYlRd)
##################################################################################
# For unit conversions, here I just use NO2 as an example
# For VOCs species, you may want to convert to total mass or total carbon

# NO2 molecule weight: 46 g mol-1
# Avogadro constant: 6.02214086 × 10^23 mol-1
# GEOS-Chem species unit: mol mol-1
# air number density unit: cm-3

surface_NO2_mass_concentration  = surface_NO2*surface_airnumberdensity*46/(6.022*1e11)
##################################################################################
##################################################################################
# In addtion: conversions between pandas dataframe (csv/xlsx) and xarray data arrays(netcdf) sometimes can be useful
# Although it is not suggested for this task

import pandas as pd

os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/geoschem")
China_daily_surface_VOCs = xr.open_dataset('China_daily_surface_VOCs_geoschem.nc')
print(China_daily_surface_VOCs)

# read BENZ at a random location on a random day using the xarray arrays
surface_BENZ = China_daily_surface_VOCs['SpeciesConc_BENZ']
surface_BENZ_day1 = surface_BENZ.isel(time=0)
surface_BENZ_day1_location1 = surface_BENZ_day1.sel(lat = 39.5, lon=116.25)
print(surface_BENZ_day1_location1)
print("values:",surface_BENZ_day1_location1.values)

# convert the xarray data array to pandas data frame
def xarray_to_pandas(data):
    data = data.to_dataframe()
    data.reset_index(inplace=True)
    return data

# use BENZ as an example
BENZ_df = xarray_to_pandas(surface_BENZ) 
print("number of rows:",len(BENZ_df))
BENZ_df.head()

# read BENZ at the same location and on the same day using the pandas dataframe
test = BENZ_df[(BENZ_df['time'] == '2016-11-01 12:00:00') & 
               (BENZ_df['lat'] == 39.5) & 
               (BENZ_df['lon'] == 116.25)]
test

# now combine multiple variables into a single pandas dataframe
surface_BENZ = China_daily_surface_VOCs['SpeciesConc_BENZ']
surface_TOLU = China_daily_surface_VOCs['SpeciesConc_TOLU']
surface_C2H6 = China_daily_surface_VOCs['SpeciesConc_C2H6']
surface_ACET = China_daily_surface_VOCs['SpeciesConc_ACET']
surface_CH2O = China_daily_surface_VOCs['SpeciesConc_CH2O']
surface_C3H8 = China_daily_surface_VOCs['SpeciesConc_C3H8']
surface_EOH = China_daily_surface_VOCs['SpeciesConc_EOH']
surface_AIRNUMDEN = China_daily_surface_VOCs['Met_AIRNUMDEN']

# no idea why the for loop worked, but the results are not saved!
for data in [surface_BENZ,surface_TOLU,surface_C2H6,surface_ACET,surface_CH2O,surface_C3H8,surface_EOH,surface_AIRNUMDEN]:
        data= xarray_to_pandas(data)
        print(data)
        
# for now, use this to perform the same function for all items and save the results into a single pandas dataframe
surface_VOCs = [surface_BENZ,surface_TOLU,surface_C2H6,surface_ACET,surface_CH2O,surface_C3H8,surface_EOH,surface_AIRNUMDEN]
surface_VOCs = list(map(xarray_to_pandas, surface_VOCs))
surface_VOCs = pd.concat(surface_VOCs,axis=1)
surface_VOCs = surface_VOCs.loc[:,~surface_VOCs.columns.duplicated()]  # remove duplicated columns
surface_VOCs.head()

# save output to csv
surface_VOCs.to_csv("China_surface_VOCs.csv",index=False,sep=',')
##################################################################################
# convert pandas data frame to xarray 
surface_VOCs = xr.Dataset.from_dataframe(surface_VOCs)

# edit attributes
surface_VOCs.attrs = {'title':'Daily surface VOCs in China'}

# edit attributes
surface_VOCs.attrs = {'title':'Daily surface VOCs in China'}
surface_VOCs['lat'].attrs = {'unit':'Degrees_east'}
surface_VOCs['lon'].attrs={'unit':'Degrees_north'}

# check the results
print(surface_VOCs)

# save outputs to netcdf
surface_VOCs.to_netcdf('China_surface_VOCs.nc')

# End
##################################################################################    
