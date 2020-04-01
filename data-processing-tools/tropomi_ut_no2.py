#!/usr/bin/python

''' Code to test feasibility of obtaining NO2 in the upper tropspheree
    from TROPOMI.

    First step is to read and find TROPOMI data points with optically
    thick clouds above 450 hPa.'''


# import relevant packages:
import os
import sys
import glob
import netCDF4 as nc4
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib.colors import LogNorm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
import matplotlib.colors as colors
import leastsq
from cloud_slice_ut_no2 import cldslice

from gcpy.gcpy.constants import K_BOLTZMANN as kB
from gcpy.gcpy.constants import R_GAS_UNIV as Runi
from gcpy.gcpy.constants import G as g
from gcpy.gcpy.constants import AVOGADRO as na
from gcpy.gcpy.constants import MW_AIR as mmair

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# Select season to process:
seas='son'

# Define target grid:
StrRes='4x5'

# Define cloud product to be used (either fresco or dlr-ocra):
cldprd='dlr-ocra'   #'fresco'    # dlr-ocra

log=open("log_"+seas+"_"+StrRes+"_"+cldprd, "w")
sys.stdout=log

print("Processing data for "+seas+" at "+StrRes+" with "+\
      cldprd+" cloud product",flush=True)

# Define grid:
if StrRes=='1x1': dellat,dellon=1,1
if StrRes=='2x25': dellat,dellon=2,2.5
if StrRes=='4x5': dellat,dellon=4,5
#dellat=1#2#2.5#2#0.5#4#0.5
#dellon=1#2.5#3.125#2.5#0.5#5#0.5
out_lon=np.arange(-180,181,dellon)
out_lat=np.arange(-90.,91.,dellat)
# Convert output lats and long to 2D:
X, Y = np.meshgrid(out_lon, out_lat,indexing='ij') 
# Dimensions of output data:
xdim=len(out_lon)
ydim=len(out_lat)
nval=int(700*dellat*dellon)  #10000 Scale with resolution

#Define output arrays for time period of interest (month for now):
gno2vmr=np.zeros((xdim,ydim)) # NO2 VMR in pptv
gcnt=np.zeros((xdim,ydim)) # No. of data points (orbits)
gerr=np.zeros((xdim,ydim))
num=np.zeros(5)
cnt=0
maxpnts=0

#Define conversion factor for no2 VMR:
#den2mr=np.divide((np.multiply(g,mmair)),na)
den2mr=np.divide((np.multiply(g,mmair)),na)

#out_lon_edge=latlontools.latlon_est_bnds(out_lon)
#out_lat_edge=latlontools.latlon_est_bnds(out_lat)

#Initialize:
first=0

# Define pressure range:
pmin=180     #380   #290   #180   #190
pmax=450     #450   #380   #290   #440

# Define string to represent the layer range:
prange=str(pmin)+'-'+str(pmax)

# Define information for each season:
if seas=='jja': 
    strmon=['06','07','08']
    stryr=['2019','2019','2019']
    yrrange='2019'
if seas=='son': 
    strmon=['09','10','11']
    stryr=['2019','2019','2019']
    yrrange='2019'
if seas=='djf': 
    strmon=['12','01','02']
    stryr=['2019','2020','2020']
    yrrange='2019-2020'

# get string of S5P TROPOMI NO2 filenames in directory of interest:
tomidir='/data/uptrop/nobackup/tropomi/Data/'
no2type='OFFL'
no2dir=tomidir+'NO2_'+no2type+'/'+stryr[0]+'/'
files=glob.glob(no2dir+strmon[0]+'/'+'S5P_'+no2type+'_L2__NO2____'+\
                stryr[0]+strmon[0]+'*')
mon=[strmon[1],strmon[2]]
yy=[stryr[1],stryr[2]]
for i,mon in enumerate(mon):
    no2dir='/data/uptrop/nobackup/tropomi/Data/NO2_'+no2type+'/'+yy[i]+'/'
    for filename in glob.glob(no2dir+mon+'/'+'S5P_'+no2type+'_L2__NO2____'+\
                              yy[i]+mon+'*'):
        files.append(filename)
files=sorted(files)

# Get string of S5P TROPOMI cloud product file names:
if cldprd=='dlr-ocra':   
    clddir=tomidir+'CLOUD_OFFL/'+stryr[0]+'/'
    cldfiles=glob.glob(clddir+strmon[0]+'/'+'S5P_OFFL_L2__CLOUD__'+\
                    stryr[0]+strmon[0]+'*')
    mon=[strmon[1],strmon[2]]
    yy=[stryr[1],stryr[2]]
    for i,mon in enumerate(mon):
        clddir=tomidir+'CLOUD_OFFL/'+yy[i]+'/'
        for filename in glob.glob(clddir+mon+'/'+'S5P_OFFL_L2__CLOUD__'+\
                                  yy[i]+mon+'*'):
            cldfiles.append(filename)
            cldfiles=sorted(cldfiles)
    print('No. of cloud files: ',len(cldfiles),flush=True)
    # Check for uneven number of files:
    if len(cldfiles) != len(files):
        print('NO2 files = ',len(files),flush=True)
        print('CLOUD files = ',len(cldfiles),flush=True)
        print('unequal number of files',flush=True)    
        sys.exit()

# Check point:
print('No. of files: ',len(files),flush=True)

# Save % retained valid data in each file:
postfilt=np.zeros(len(files))

# loop over files:
for f,files in enumerate(files):

    print('Processing:',files[-86:],flush=True)

    # Count valid data for each orbit:
    #NOrb=0.

    fh=Dataset(files,mode='r')
    #print(fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].variables['cloud_fraction_crb'])

    # These actually only need to be defined once (consider changing in 
    # the future).
    # Get NO2 conversion factor:
    no2sfac=fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                '_column'].multiplication_factor_to_convert_to_molecules_percm2

    # Get QA flag scale factor:
    qasfac=fh.groups['PRODUCT'].variables['qa_value'].scale_factor

    # NO2 fill/missing value:
    fillval=fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                                           '_column']._FillValue

    # Extract data of interest (lon, lat, clouds, NO2 total column & error):

    # Geolocation data:
    glons=fh.groups['PRODUCT'].variables['longitude'][:]#[:][0,:,:]
    tlons=glons[0,:,:]
    glats=fh.groups['PRODUCT'].variables['latitude'][:]#[:][0,:,:]
    tlats=glats[0,:,:]

    # Columns data:
    # (1) Tropospheric vertical column:
    gtropno2 =fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                                             '_column'][:]
    ttropno2=gtropno2.data[0,:,:]
    # (2) Total vertical column:
    gtotno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
              variables['nitrogendioxide_total_column'][:]
    ttotno2=gtotno2.data[0,:,:] 
    # (3) Total slant column:
    gscdno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
            variables['nitrogendioxide_slant_column_density'][:]
    tscdno2=gscdno2.data[0,:,:]    
    # (4) Stratospheric vertical column:
    gstratno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
              variables['nitrogendioxide_stratospheric_column'][:]
    tstratno2=gstratno2.data[0,:,:]
    # (5) Summed vertical column:
    gsumno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
             variables['nitrogendioxide_summed_total_column'][:]
    tsumno2=gtotno2.data[0,:,:] 

    # AMF for cloudy portion of pixel from cloud top pressure to
    # the tropopause:
    gtropcldamf=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
              variables['air_mass_factor_cloudy'][:]
    ttropcldamf=gtropcldamf.data[0,:,:] 
    # Stratospheric AMF:
    gstratamf=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
               variables['air_mass_factor_stratosphere'][:]
    tstratamf=gstratamf.data[0,:,:] 
    # Total AMF:
    gtotamf=fh.groups['PRODUCT'].variables['air_mass_factor_total'][:]
    ttotamf=gtotamf.data[0,:,:] 
    # Get sum of the two for cloudy AMF from cloud top pressure to the 
    # top of the atmosphere:
    #tcldamf=np.add(ttropcldamf,tstratamf)    #np.multiply(0.5,np.add(ttropcldamf,tstratamf)) #  np.add(ttropcldamf,tstratamf)
    tcldamf=np.add(ttropcldamf,tstratamf)

    # QA value:
    qaval=fh.groups['PRODUCT'].variables['qa_value'][0,:,:]

    # Precisions:
    tropno2err =fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                                           '_column_precision'][0,:,:]
    totno2err=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
               variables['nitrogendioxide_total_column_precision'][0,:,:]
    stratno2err=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
                variables['nitrogendioxide_stratospheric_column_precision'][0,:,:]

    # Cloud input data if fresco cloud product is being used:
    if ( cldprd=='fresco' ):
        tcldfrac=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                  variables['cloud_fraction_crb'][:]
        cldfrac=tcldfrac.data[0,:,:]
        #tcldfrac._FillValue
        #tcldalb=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
        #         variables['cloud_albedo_crb'][:]
        #cldalb=tcldalb.data[0,:,:]
        gcldpres=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                  variables['cloud_pressure_crb'][:]
        tcldpres=np.ma.getdata(gcldpres[0,:,:]) # extract data from masked array

    # Aerosol absorbing index:
    taai=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['aerosol_index_354_388'][:]
    aai=taai.data[0,:,:]
    aai=np.where(aai>1e30, np.nan, aai)

    # Solar zenith angle (degrees):
    tsza=fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].\
          variables['solar_zenith_angle'][:]
    sza=tsza[0,:,:]

    # close file:
    fh.close()

    # Cloud input data if dlr-ocra cloud product is used:
    if ( cldprd=='dlr-ocra' ):
        
        # Read data:
        fh=Dataset(cldfiles[f],mode='r')

        # Check that date is the same as the no2 file:
        strdate=cldfiles[f]
        strdate=strdate[-66:-51]
        if strdate!=files[-66:-51]:
            print('NO2 file, Cloud file: '+files[-66:-51]+", "+strdate,\
                  flush=True)
            print('EXITING: Files are not for the same date!',flush=True)
            sys.exit()
        
        # Get cloud fraction and cloud top pressure:
        tcldfrac=fh.groups['PRODUCT'].variables['cloud_fraction'][:]
        cldfrac=tcldfrac.data[0,:,:]
        gcldpres=fh.groups['PRODUCT'].variables['cloud_top_pressure'][:]
        tcldpres=np.ma.getdata(gcldpres[0,:,:]) # extract data from masked array

        # QA value:
        cldqa=fh.groups['PRODUCT'].variables['qa_value'][0,:,:]

        # Set poor quality cloud data to nan:
        cldfrac=np.where(cldqa<0.5, np.nan, cldfrac )
        tcldpres=np.where(cldqa<0.5, np.nan, tcldpres )
        
        # Check that data shapes are equal:
        if cldfrac.shape!=sza.shape:
            print('Cloud product and NO2 indices ne!',flush=True)
            print(cldfrac.shape,sza.shape,flush=True)
            print('Skipping this swath',flush=True)
            fh.close()
            continue

        # Close file:
        fh.close()

    # Apply scale factor (QA value only):
    #qaval=np.multiply(qaval,qasfac)

    # Convert missing total column NO2 values to NAN before regridding:
    #ttotno2[gtotno2.mask[0,:,:]]=float("nan")
    # Account for files where mask is missing (only appears to be one):
    if len(gscdno2.mask.shape)==0:
        tscdno2=np.where(tscdno2==fillval, np.nan, tscdno2)
    else:
        tscdno2[gscdno2.mask[0,:,:]]=float("nan") 

    # Find relevant data only:
    # Filter to only include very cloudy scenes at high altitude
    # No. of valid points before filter:
    inicnt=np.count_nonzero(~np.isnan(tscdno2))
    # Filter out scenes with cloud fraction < 0.7:
    tscdno2=np.where(cldfrac<0.7, np.nan, tscdno2)
    # Account for cloud fractions that are NANs, but the NO2 column
    # is not:
    if cldprd=='dlr-ocra':
        nanind=np.where(np.isnan(cldfrac))[0]
        tscdno2=np.where(np.isnan(cldfrac), np.nan, tscdno2)
    # Filter out scenes outside the UT (200-450 hPa):
    tscdno2=np.where(tcldpres>pmax*1e2, np.nan, tscdno2)
    tscdno2=np.where(tcldpres<pmin*1e2, np.nan, tscdno2)
    # Filter out low quality retrieval scenes (0.45 suggested
    # by Henk Eskes, KNMI, if targeting cloudy scenes):
    tscdno2=np.where(qaval<0.45, np.nan, tscdno2)
    # Filter out scenes with AAI > 1 (could be misclassified as clouds)
    tscdno2=np.where(aai>1., np.nan, tscdno2)
    # Cloud albedo threshold (not used)
    #tscdno2=np.where(cldalb<=0.8, np.nan, tscdno2)

    # No. of point retained after filtering:
    tcnt=np.count_nonzero(~np.isnan(tscdno2))
    # Save % valid points retained:
    postfilt[f]=100.*(tcnt/inicnt)

    # Trim the data to include only those relevant:
    # This also reshapes the data from a 2D to a 1D array:
    lons=tlons[~np.isnan(tscdno2)]
    lats=tlats[~np.isnan(tscdno2)]
    cldpres=tcldpres[~np.isnan(tscdno2)]
    cldamf=tcldamf[~np.isnan(tscdno2)]
    sumno2=np.multiply(tsumno2[~np.isnan(tscdno2)],no2sfac)
    stratno2=tstratno2[~np.isnan(tscdno2)]
    scdno2=tscdno2[~np.isnan(tscdno2)]

    # Define output data:
    gno2=np.zeros((xdim,ydim,nval))
    gstrat=np.zeros((xdim,ydim,nval))
    gcldp=np.zeros((xdim,ydim,nval))
    #gcnt=np.zeros(X.shape,nval)
    cntloop=np.zeros((xdim,ydim))

    # Regrid by dump and average approach:
    for w, scdno2 in enumerate(scdno2):
        if ( scdno2!=np.nan ):
            
            # Find nearest gridsquare in output grid:
            p,q = np.argmin(abs(out_lon-lons[w])),\
                  np.argmin(abs(out_lat-lats[w]))

            # Convert total column NO2 from SCD to VCD:
            tvcdno2=np.divide(scdno2,cldamf[w])

            # Convert columns from mol/m2 to molec/cm2:
            tvcdno2=np.multiply(tvcdno2,no2sfac) 
            tstrat=np.multiply(stratno2[w],no2sfac) 

            # Error check:
            #if tstrat>tvcdno2:
            #    # Check by how much:
            #    print('Stratosphere, total: ',\
            #          int(tstrat*1e-13),int(tvcdno2*1e-13))
            #    print('EXITING: Stratospheric column > total column')
            #    sys.exit()

            # Check that summed column from file is greater than partial
            # column estimated using the slant column over very cloudy
            # scenes and the cloudy troposphere AMF and stratosphere AMF
            #diff=sumno2[w]-tvcdno2
            #if ( w>100 and w<110 ): 
            #    print(tvcdno2*1e-15,sumno2[w]*1e-15,diff*1e-15)
           
            # Get relevant data in each output grid square:
            # Also convert no2 to molec/cm2:
            gno2[p,q,int(cntloop[p,q])]=tvcdno2
            gstrat[p,q,int(cntloop[p,q])]=tstrat
            gcldp[p,q,int(cntloop[p,q])]=cldpres[w]

            # Increment indices:
            cntloop[p,q]=cntloop[p,q]+1

    # Estimate daily mean VMRs:
    for i in range(xdim):
        for j in range(ydim):
            if gno2[i,j,0]!=0:
                # Get NO2 partial columns above cloudy scenes:
                tcolno2=gno2[i,j,:]
                # Trim to remove trailing zeros:
                tcolno2=tcolno2[np.nonzero(tcolno2)]
                # Convert NO2 from molec/cm2 to molec/m2:
                tcolno2=np.multiply(tcolno2,1e4)

                # Skip if fewer than 10 points:
                if len(tcolno2)<10: 
                    num[0]=num[0]+1
                    continue

                # Get stratospheric columns:
                strat=gstrat[i,j,:]
                # Trim to remove trailing zeros:
                strat=strat[np.nonzero(strat)]
                # Convert NO2 from molec/cm2 to molec/m2:
                strat=np.multiply(strat,1e4)

                # Get cloud top pressure and convert from Pa to hPa:
                tcld=gcldp[i,j,:]*1e-2
                # Trim to remove trailing zeros:
                tcld=tcld[np.nonzero(tcld)]

                # Skip scenes with non-uniform stratosphere using the
                # same threshold as is used for GEOS-Chem:
                if (np.std(strat)/np.mean(strat))>=0.02:
                    #print(i,j,(np.std(strat)/np.mean(strat)))
                    continue

                # Find and remove outliers (2-sigma):
                #hiind=np.where(tcolno2>(np.mean(tcolno2)+2*np.std(tcolno2)))
                #loind=np.where(tcolno2<(np.mean(tcolno2)-2*np.std(tcolno2)))
                #tcolno2[hiind]=0.0
                #tcolno2[loind]=0.0
                #tcld[hiind]=0.0
                #tcld[loind]=0.0
                # Remove additional zeros as a result of outlier filtering:
                #tcolno2=tcolno2[np.nonzero(tcolno2)]
                #tcld=tcld[np.nonzero(tcld)]

                # Get number of points:
                npnts=len(tcld)
                if npnts>maxpnts:
                    maxpnts=npnts
                    print(maxpnts,flush=True)

                # Get mean cloud top pressure:
                pmean=np.mean(tcld)
                # Calculate Gaussian weight:
                gwgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(pmean,315.))),(np.multiply(2,np.square(135.))))))
                    
                # Use cloud_slice_ut_no2 function to get NO2 mixing
                # ratio from cloud-slicing:
                if npnts<=30:
                    csval=cldslice(tcolno2,tcld)
                    if np.isnan(csval[0]) or np.isnan(csval[1]):
                        num[csval[2]-1]=num[csval[2]-1]+1
                        continue
                    #NObr=NOrb+npnts
                    # Cloud-slice estimated NO2 VMRs:
                    # Error weighted:
                    #gno2vmr[i,j]=gno2vmr[i,j]+np.divide(csval[0],(csval[1]**2))
                    #gerr[i,j]=gerr[i,j]+np.divide(1,(csval[1]**2))
                    # Gaussian weighted:
                    gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                    gerr[i,j]=gerr[i,j]+gwgt
                    gcnt[i,j]=gcnt[i,j]+1
                    cnt=cnt+1
                elif npnts>30:
                    # Order the data by cloud height:
                    # Get indices that would sort an array:
                    sind=np.argsort(tcld)

                    # Code error check:
                    if len(tcld)!=len(tcolno2): 
                        print('arrays not =',flush=True)
                        print(len(tcld),len(tcolno2),flush=True)
                    
                    # Sorted arrays (although, not sorting these any more):
                    sortcld=tcld#[sind]
                    sortno2=tcolno2#[sind]

                    niter=round(npnts/17)
                    nloop=list(range(niter))
                    # Loop over data:
                    for w in nloop:
                        if w==0:
                            csval=cldslice(sortno2[::niter],sortcld[::niter])
                        else:
                            csval=cldslice(sortno2[w::niter],\
                                           sortcld[w::niter])
                        if np.isnan(csval[0]) or np.isnan(csval[1]):
                            num[csval[2]-1]=num[csval[2]-1]+1
                            continue
                        # Cloud-slice estimated NO2 VMRs:   
                        # Error weighted mean:
                        #gno2vmr[i,j]=gno2vmr[i,j]+np.divide(csval[0],(csval[1]**2))
                        #gerr[i,j]=gerr[i,j]+np.divide(1,(csval[1]**2))
                        # Gaussian weighted mean:
                        gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                        gerr[i,j]=gerr[i,j]+gwgt
                        gcnt[i,j]=gcnt[i,j]+1
                        cnt=cnt+1

print('No. of valid data points: ',cnt,flush=True)
print('Max no. of data points in a gridsquare: ',np.amax(gcnt),flush=True)
print('No. of data lost because: ',flush=True)
print('(1) Too few points: ',num[0],flush=True)
print('(2) Low cloud height variability: ',num[1],flush=True)
print('(3) Large error: ',num[2],flush=True)
print('(4) Significantly less then zero: ',num[3],flush=True)
print('(5) Outlier (NO2 > 200 pptv): ',num[4],flush=True)

print('Mean % points retained: ', np.mean(postfilt),flush=True)

print('Maximum number of points in grid: ', maxpnts,flush=True)

# Next:
# Work on looping over files to get monthly mean UT NO2 (and 
# check again whether there's a conversion error).

#sys.exit()

# Get average:
gno2vmr=np.divide(gno2vmr, gerr, where=gcnt!=0)
#gno2vmr=np.divide(gno2vmr, gcnt, where=gcnt!=0)
gno2vmr[gcnt==0]=np.nan
#gno2vmr[gcnt==0]=np.nan
#gerr=np.divide(1,np.sqrt(gerr), where=gerr!=0)
gerr=np.divide(gerr, gcnt, where=gcnt!=0)
gerr[gcnt==0]=np.nan
#gerr[gcnt==0]=np.nan
gcnt[gcnt==0]=np.nan

# Save the data to NetCDF:
ncout=Dataset('./Data/tropomi-ut-no2-'+cldprd+'-'+StrRes+'-'+seas+'-'+\
              yrrange+'-v1.nc',mode='w',format='NETCDF4') 
#print(ncfile)

#ncout.createDimension('time', TDim)  
ncout.createDimension('lat', ydim) 
ncout.createDimension('lon', xdim) 

# create longitude axis:
lon = ncout.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longiitude'
lon[:] = out_lon

# Create latitude axis:
lat = ncout.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
lat[:] = out_lat

# Save data: UT NO2 VMR (gno2vmr), UT NO2 error (gerr), No. of data points (gcnt)
utno2 = ncout.createVariable('utno2', np.float32, ('lon','lat'))
utno2.units = 'pptv'
utno2.long_name = 'NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
utno2[:] = gno2vmr

utno2err = ncout.createVariable('utno2err', np.float32, ('lon','lat'))
utno2err.units = 'pptv'
utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
utno2err[:] = gerr

nobs = ncout.createVariable('nobs', np.float32, ('lon','lat'))
nobs.units = 'unitless'
nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
nobs[:] = gcnt

ncout.close()

# Plot the data:
m=Basemap(resolution='l',projection='merc',\
          lat_0=0,lon_0=0,llcrnrlon=-180,\
          llcrnrlat=-75,urcrnrlon=180,urcrnrlat=80)
xi,yi=m(X,Y)
plt.subplot(1, 3, 1)
cs=m.pcolor(xi,yi,np.squeeze(gno2vmr), vmin=0, vmax=150, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('NO2 VMRs')

plt.subplot(1, 3, 2)
cs=m.pcolor(xi,yi,np.squeeze(gerr), vmin=0, vmax=20, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('NO2 error')

plt.subplot(1, 3, 3)
cs=m.pcolor(xi,yi,np.squeeze(gcnt), vmin=0., vmax=50, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Number of points')

#plt.show()

#plt.savefig('./Images/test-plot-tropomi-ut-no2-'+Seas+'-v1.ps', \
#            format='ps')

# Close the log file:
log.close()
