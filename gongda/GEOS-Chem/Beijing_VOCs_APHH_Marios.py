##################################################################################
##################################################################################
# Read daily surface VOCs at IAP (Beijing) during APHH Campaign (2016-11) from GEOS-Chem 12.0.0
# Save variables of interest into one single netcdf file for use in Mario's Beijing VOCs paper

# For full GEOS-Chem species info (shortname,formula,molec weight, phase...)
# http://wiki.seas.harvard.edu/geos-chem/index.php/Species_in_GEOS-Chem
# C2H2 is listed in the model wiki, but not found in the output file 
# maybe it is added in the newer model versions?

import os
import glob
import xarray as xr
##################################################################################
# load daily model output files during 2016-11 (default MEIC emissions)
os.chdir("/rds/projects/2018/maraisea-glu-01/RDS/GEOSChem/GEOS-Chem_12.0.0/Run_Directory/Beijing_VOCs/merra2_05x0625_tropchem_as_2016_November_default_MEIC/APHH_VOCs_daily_outputs")

Species  = sorted(glob.glob("GEOSChem.SpeciesConc*.nc4"))
StateMet = sorted(glob.glob("GEOSChem.StateMet*.nc4"))

print("number of files:",len(Species),"first file:",Species[0],"last file:",Species[-1],sep=" ")
print("number of files:",len(StateMet),"first file:",StateMet[0],"last file:",StateMet[-1],sep=" ")
##################################################################################
# read model outputs from each day
GC_species  = [xr.open_dataset(file) for file in Species]
GC_statemet = [xr.open_dataset(file) for file in StateMet]

# Extract species of interest: 
# Benzene/Toluene/Ethane/Acetylene/Acetone/Formaldehyde/Propane/Xylenes/nbutane/ibutane/Ethylene/Ethanol
surface_BENZ = [data['SpeciesConc_BENZ'].isel(time=0,lev=0) for data in GC_species] # Benz (C6H6)
surface_TOLU = [data['SpeciesConc_TOLU'].isel(time=0,lev=0) for data in GC_species] # Toluene (C7H8)
surface_C2H6 = [data['SpeciesConc_C2H6'].isel(time=0,lev=0) for data in GC_species] # Ethane (C2H6)
surface_ACET = [data['SpeciesConc_ACET'].isel(time=0,lev=0) for data in GC_species] # Acetone (CH3C(O)CH3) 
surface_CH2O = [data['SpeciesConc_CH2O'].isel(time=0,lev=0) for data in GC_species] # Formaldehyde (CH2O,HCHO)
surface_C3H8 = [data['SpeciesConc_C3H8'].isel(time=0,lev=0) for data in GC_species] # Propane (C3H8)
surface_EOH  = [data['SpeciesConc_EOH'].isel(time=0,lev=0)  for data in GC_species] # Ethanol (C2H5OH)

# check the units of each species using a sample file
for data in [surface_BENZ[0],surface_TOLU[0],surface_C2H6[0],surface_ACET[0],
             surface_CH2O[0],surface_C3H8[0],surface_EOH[0]]:
    print(data.units)
    
# so these species are in the unit of "mol/mol"
# need to use the number density ("Met_AIRNUMDEN") to derive mass concentrations
# or multiply by "1e9" to convert to "ppbv"
##################################################################################
# read air number density from "GEOSChem.StateMet" and fix the unit

# Important note: the unit of air number density in GEOS-Chem output is wrong!
# It should be cm-3, instead of m-3
# Ways to prove this: 1> if it is m-3, then values are not of the same magnitude with those derived from ideal gas law
#                     2> if it is m-3, the mass concentrations in the end will be around 10^6 times lower than what it should be 

# The "Met_AIRNUMDEN" is also not included in the "StateMet" diagnostics on GEOS-Chem wiki
# In order to output "Met_AIRNUMDEN", you need to manually add "Met_AIRNUMDEN" to "History.rc" when setting up the model

surface_airnumberdensity = [data['Met_AIRNUMDEN'].isel(time=0,lev=0) for data in GC_statemet] # air number density

for i in range(len(surface_airnumberdensity)):
    surface_airnumberdensity[i].attrs['long_name'] = 'Dry air number density'
    surface_airnumberdensity[i].attrs['units'] = 'cm-3'
    
# combine all relevant variables
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
    
# add global attributes from the raw GEOS-Chem output file
for i in range(len(China_daily_VOCs)):
    China_daily_VOCs[i].attrs = GC_species[i].attrs
    China_daily_VOCs[i].attrs['title'] = 'Daily GEOS-Chem surface VOCs concentrations in China'
    
# now all the data fields needed are combined on each day
print(China_daily_VOCs[0])

# combine data from multiple days into one single file
China_VOCs_201611 = xr.concat(China_daily_VOCs,'time')
print(China_VOCs_201611)

# save outputs to netcdf
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/1_co_author_papers/VOCs_Marios")
China_VOCs_201611.to_netcdf('China_daily_surface_VOCs_geoschem_base_MEIC.nc')
##################################################################################
# check the netcdf file generated
os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/1_co_author_papers/VOCs_Marios")
China_daily_surface_VOCs = xr.open_dataset('China_daily_surface_VOCs_geoschem_base_MEIC.nc')
print(China_daily_surface_VOCs)

# plot on map for a sanity check
%matplotlib inline
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

# a test plot
os.chdir('/rds/projects/2018/maraisea-glu-01/Study/Research_Data/BTH/domain/CHN_ADM')
China_shape = r'CHN_ADM1.shp'
china_map = ShapelyFeature(Reader(China_shape).geometries(),ccrs.PlateCarree(), edgecolor='black',facecolor='none')

t= plt.axes(projection=ccrs.PlateCarree())
t.add_feature(china_map)
t.set_extent([100, 125, 30, 45], crs=ccrs.PlateCarree()) # [lon_min,lon_max,lat_min,lat_max]

VOC_test1 = China_daily_surface_VOCs['SpeciesConc_BENZ'].isel(time=0)
VOC_test1.plot(ax=t,cmap='jet')
##################################################################################
##################################################################################
# Codes below explains how to analyse the outputs from GEOS-Chem

# Convert from mol/mol to mass concentrations, using NO2 as an example
# NO2 molecule weight: 46 g mol-1
# Avogadro constant: 6.02214086 × 10^23 mol-1
# GEOS-Chem species unit: mol mol-1
# air number density unit: cm-3
surface_NO2_mass_concentration  = surface_NO2*surface_airnumberdensity*46/(6.022*1e11)

# Cnovert from mol/mol to ppbv (since both are mixing ratios, just apply the factor):
surface_NO2_ppbv = surface_NO2*1e9
##################################################################################
# In addtion: conversions between pandas dataframe (csv/xlsx) and xarray data arrays(netcdf) sometimes can be useful
# Although it is not suggested here for this task.
import pandas as pd

os.chdir("/rds/projects/2018/maraisea-glu-01/Study/Research_Data/1_co_author_papers/VOCs_Marios")
China_daily_surface_VOCs = xr.open_dataset('China_daily_surface_VOCs_geoschem_base_MEIC.nc')

# currently it is a xarray object
print(China_daily_surface_VOCs)

# now read BENZ at IAP using xarray
BENZ = China_daily_surface_VOCs['SpeciesConc_BENZ']
BENZ_IAP = BENZ.sel(lat = 39.5, lon=116.25)
BENZ_IAP = BENZ_IAP*1e9
print("#"*50,"surface BENZ at IAP:",BENZ_IAP,sep="\n")
print("#"*50,"number of days:",len(BENZ_IAP.values),sep="\n")
print("#"*50,"data type:",type(BENZ_IAP),sep="\n")

# convert the xarray to a pandas dataframe
BENZ_IAP = BENZ_IAP.to_dataframe().reset_index()
BENZ_IAP

# build a function to read other species at IAP
def read_VOCs_IAP(gc_species):
    data = gc_species
    data = data.sel(lat=39.5,lon=116.25)
    data = data*1e9
    data = data.to_dataframe().reset_index()
    return data

# read other species at IAP
TOLU_IAP = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_TOLU'])
C2H6_IAP = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_C2H6'])
ACET_IAP = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_ACET'])
CH2O_IAP = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_CH2O'])
C3H8_IAP = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_C3H8'])
EOH_IAP  = read_VOCs_IAP(China_daily_surface_VOCs['SpeciesConc_EOH'])

# combine all species into a single pandas dataframe
VOCs_IAP = [BENZ_IAP,TOLU_IAP,C2H6_IAP,ACET_IAP,CH2O_IAP,C3H8_IAP,EOH_IAP]
VOCs_IAP = pd.concat(VOCs_IAP,axis=1)

# remove duplicated columns
VOCs_IAP = VOCs_IAP.loc[:,~VOCs_IAP.columns.duplicated()] 
VOCs_IAP.head()

# save output to csv
VOCs_IAP.to_csv("IAP_surface_VOCs.csv",index=False,sep=',')
##################################################################################
# Alternatively, use "map" or "list comprehensions" to perform the same function to a list of objects
# not needed here, but it is good to know the "map" and "list comprehensions" in Python

# read all species seperately
surface_BENZ = China_daily_surface_VOCs['SpeciesConc_BENZ']
surface_TOLU = China_daily_surface_VOCs['SpeciesConc_TOLU']
surface_C2H6 = China_daily_surface_VOCs['SpeciesConc_C2H6']
surface_ACET = China_daily_surface_VOCs['SpeciesConc_ACET']
surface_CH2O = China_daily_surface_VOCs['SpeciesConc_CH2O']
surface_C3H8 = China_daily_surface_VOCs['SpeciesConc_C3H8']
surface_EOH  = China_daily_surface_VOCs['SpeciesConc_EOH']
surface_AIRNUMDEN = China_daily_surface_VOCs['Met_AIRNUMDEN']

# perform the same functions to all
surface_VOCs = [surface_BENZ,surface_TOLU,surface_C2H6,surface_ACET,surface_CH2O,surface_C3H8,surface_EOH,surface_AIRNUMDEN]
# using "map"
surface_VOCs = list(map(read_VOCs_IAP,surface_VOCs))
# or using "list comprehension"
surface_VOCs = [read_VOCs_IAP(data) for data in surface_VOCs]

# merge the data frames and remove duplicated columns
surface_VOCs = pd.concat(surface_VOCs,axis=1)
surface_VOCs = surface_VOCs.loc[:,~surface_VOCs.columns.duplicated()] 
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
