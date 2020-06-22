#!/usr/bin/python

''' Code to apply cloud-slicing to TROPOMI total columns of NO2 to obtain
    seasonl mean NO2 mixing ratios in the global upper troposphere.

    The code includes calculation of a geometric tropospheric column,
    bias correction of the TROPOMI stratospheric and tropospheric columns
    based  on assessment of TROPOMI against Pandora, and a choice to 
    use either FRESCO or DLR cloud pressure height and fraction.

    Uses the project-specific external functions cloudslice and alt2pres.
    '''

# import relevant packages:
import os
import sys
import glob
import netCDF4 as nc4
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from cloud_slice_ut_no2 import cldslice
from convert_height_to_press import alt2pres

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# JR: seas, StrRes, cldprd, cldthld can all be input arguments:
# Select season to process:
seas='mam'    # options are jja, son, djf, mam; default can be jja

# Define target grid:
StrRes='1x1'    # options are 1x1, 2x2.5, 4x5; default is 1x1.

# Define cloud product to be used (either fresco or dlr-ocra):
cldprd='dlr-ocra'   #'fresco'    # dlr-ocra

# Define cloud threshold string:
strcldthld='07'   # options are 07, 08, 09, 10; default is 07

# Define cloud threshold string:
if strcldthld=='07': cldthld=0.7
if strcldthld=='08': cldthld=0.8
if strcldthld=='09': cldthld=0.9
if strcldthld=='10': cldthld=1.0

# Log file:
log=open("log_"+seas+"_"+StrRes+"_"+cldprd+"_"+strcldthld+"_v7", "w")
sys.stdout=log

# Track progress:
print("Processing data for "+seas+" at "+StrRes+" with "+\
      cldprd+" cloud product",flush=True)

# Define grid:
if StrRes=='1x1': dellat,dellon=1,1
if StrRes=='2x25': dellat,dellon=2,2.5
if StrRes=='4x5': dellat,dellon=4,5
out_lon=np.arange(-180,180+dellon,dellon)
out_lat=np.arange(-90,90+dellat,dellat)
# Convert output lats and long to 2D:
X, Y = np.meshgrid(out_lon, out_lat,indexing='ij') 
# Dimensions of output data:
xdim=len(out_lon)
ydim=len(out_lat)
#JR: this is the same issue as the ut_no2_gc_test.py file, where instead of this approach, a variable can be initiated and additional data can be appended to it.
nval=int(700*dellat*dellon)  #Size scales with resolution

#Define final seasonal mean output arrays for the season of interest:
gno2vmr=np.zeros((xdim,ydim)) # NO2 VMR in pptv
gcnt=np.zeros((xdim,ydim)) # No. of data points (orbits)
gerr=np.zeros((xdim,ydim)) # Weighted error

# Initialize additional diagnostics:
num=np.zeros(7)  # Track reasons for data loss.
cnt=0
maxpnts=0

# Define pressure range:
pmin=180 
pmax=450 

# Define information for season selected:
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
if seas=='mam': 
    strmon=['03','04','05']
    stryr=['2020','2020','2020']
    yrrange='2020'

# get string of S5P TROPOMI NO2 filenames in directory of interest:
tomidir='/data/uptrop/nobackup/tropomi/Data/'
no2type='OFFL'
no2dir=tomidir+'NO2_'+no2type+'/'+stryr[0]+'/'
files=glob.glob(no2dir+strmon[0]+'/'+'S5P_'+no2type+'_L2__NO2____'+\
                stryr[0]+strmon[0]+'01*')
mon=[strmon[1],strmon[2]]
yy=[stryr[1],stryr[2]]
for i,mon in enumerate(mon):
    no2dir='/data/uptrop/nobackup/tropomi/Data/NO2_'+no2type+'/'+yy[i]+'/'
    for filename in glob.glob(no2dir+mon+'/'+'S5P_'+no2type+'_L2__NO2____'+\
                              yy[i]+mon+'01*'):
        files.append(filename)
files=sorted(files)

# Get string of S5P TROPOMI cloud product file names:
if cldprd=='dlr-ocra':   
    clddir=tomidir+'CLOUD_OFFL/'+stryr[0]+'/'
    cldfiles=glob.glob(clddir+strmon[0]+'/'+'S5P_OFFL_L2__CLOUD__'+\
                    stryr[0]+strmon[0]+'01*')
    mon=[strmon[1],strmon[2]]
    yy=[stryr[1],stryr[2]]
    for i,mon in enumerate(mon):
        clddir=tomidir+'CLOUD_OFFL/'+yy[i]+'/'
        for filename in glob.glob(clddir+mon+'/'+'S5P_OFFL_L2__CLOUD__'+\
                                  yy[i]+mon+'01*'):
            cldfiles.append(filename)
            cldfiles=sorted(cldfiles)
    print('No. of cloud files: ',len(cldfiles),flush=True)
    # Check for inconsistent number of files:
    if len(cldfiles) != len(files):
        print('NO2 files = ',len(files),flush=True)
        print('CLOUD files = ',len(cldfiles),flush=True)
        print('unequal number of files',flush=True)    
        sys.exit()

# Check point:
print('No. of NO2 files: ',len(files),flush=True)

# Initialize array to save % retained valid data in each file:
postfilt=np.zeros(len(files))

# loop over files:
for f,files in enumerate(files):

    print('Processing:',files[-86:],flush=True)

    # ======= READ THE DATA ========

    fh=Dataset(files,mode='r')

    # no2sfac, qasfac, and fillval only need to be read in once, so could
    # just be read in on the first iteration.
    # NO2 conversion factor:
    no2sfac=fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                '_column'].multiplication_factor_to_convert_to_molecules_percm2

    # QA flag scale factor:
    qasfac=fh.groups['PRODUCT'].variables['qa_value'].scale_factor

    # NO2 fill/missing value:
    fillval=fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric'\
                                           '_column']._FillValue

    # Extract data of interest:

    # Geolocation data:
    glons=fh.groups['PRODUCT'].variables['longitude'][:]
    tlons=glons[0,:,:]
    glats=fh.groups['PRODUCT'].variables['latitude'][:]
    tlats=glats[0,:,:]

    # Column data:
    # (1) Total slant column:
    gscdno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
            variables['nitrogendioxide_slant_column_density'][:]
    tscdno2=gscdno2.data[0,:,:] 
    # (2) Stratospheric vertical column:
    gstratno2=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
              variables['nitrogendioxide_stratospheric_column'][:]
    stratno2_og=gstratno2.data[0,:,:]

    # Precisions/Uncertainties of column data:
    # (1) Total slant column uncertainty:
    gscdno2err=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
                variables['nitrogendioxide_slant_column_density_precision'][:]
    tscdno2err=gscdno2err.data[0,:,:]  
    # (2) Stratospheric vertical column uncertainty:
    stratno2err=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
                variables['nitrogendioxide_stratospheric_column_precision'][0,:,:]
    # Stratospheric AMF:
    gstratamf=fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS'].\
               variables['air_mass_factor_stratosphere'][:]
    tstratamf=gstratamf.data[0,:,:] 

    # QA value:
    qaval=fh.groups['PRODUCT'].variables['qa_value'][0,:,:]

    # Aerosol absorbing index:
    taai=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['aerosol_index_354_388'][:]
    aai=taai.data[0,:,:]
    aai=np.where(aai>1e30, np.nan, aai)

    # Solar zenith angle (degrees):
    tsza=fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].\
          variables['solar_zenith_angle'][:]
    sza=np.ma.getdata(tsza[0,:,:])

    # Viewing zenith angle (degrees):
    tvza=fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS'].\
          variables['viewing_zenith_angle'][:]
    vza=np.ma.getdata(tvza[0,:,:])

    # ======= CALCULATE THE GEOMETRIC TROPOSPHERIC COLUMN ========

    # Calculate the geometric AMF:
    tamf_geo=np.add((np.reciprocal(np.cos(np.deg2rad(sza)))),\
                    (np.reciprocal(np.cos(np.deg2rad(vza)))))

    # Get VCD under cloud conditions. This is done as the current
    # tropospheric NO2 VCD product includes influence from the prior
    # below clouds:
    # Calculate the stratospheric slant columns:
    tscdstrat=np.multiply(stratno2_og,tstratamf)
    # Calculate the tropospheric slant columns:
    ttropscd=np.subtract(tscdno2,tscdstrat)
    # Calculate the tropospheric vertical column using the geometric AMF:
    tgeotropvcd=np.divide(ttropscd,tamf_geo)

    # ======= APPLY PANDORA-INFORMED BIAS CORRECTION ========

    # Bias correct stratosphere based on comparison of TROPOMI to Pandora Mauna Loa:
    tstratno2=np.where(stratno2_og!=fillval, ((stratno2_og-(7.3e14/no2sfac))/0.82), np.nan)
    
    # Bias correct troposphere based on comparison of TROPOMI to Pandora Izana:
    tgeotropvcd=np.where(tgeotropvcd!=fillval, tgeotropvcd/2., np.nan)

    # Get the total column as the sum of the bias-corrected components:
    tgeototvcd=np.add(tgeotropvcd,tstratno2)

    # Calculate updated stratospheric NO2 error after bias correcting.
    # Determine by scaling it by the relative change in stratospheric vertical
    # colum NO2 after applying a bias correction:    
    tstratno2err=np.where(stratno2err!=fillval, \
                          np.multiply(stratno2err,np.divide(tstratno2,stratno2_og)),np.nan)

    # Calculate error by adding in quadrature individual
    # contributions:
    ttotvcd_geo_err=np.sqrt(np.add(np.square(tstratno2err),\
                                   np.square(tscdno2err)))
    # Estimate the tropospheric NO2 error as the total error
    # weighted by the relative contribution of the troposphere
    # to the total column, as components that contribute to the 
    # error are the same:
    ttropvcd_geo_err=np.multiply(ttotvcd_geo_err,\
                                 (np.divide(tgeotropvcd,tgeototvcd)))

    # ======= GET CLOUD INFORMATION  ========

    # Cloud input data if fresco cloud product is being used:
    if ( cldprd=='fresco' ):
        # Cloud fraction:
        tcldfrac=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                  variables['cloud_fraction_crb'][:]
        cldfrac=tcldfrac.data[0,:,:]
        # Cloud top pressure:
        gcldpres=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                  variables['cloud_pressure_crb'][:]
        tcldpres=gcldpres[0,:,:]   

        # Get scene and surface pressure to diagnose clouds misclassified as snow/ice:
        # Apparent scene pressure:
        gscenep=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                 variables['apparent_scene_pressure'][:]
        tscenep=gscenep.data[0,:,:]
        # Surface pressure:
        gsurfp=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
                variables['surface_pressure'][:]
        tsurfp=gsurfp.data[0,:,:]

        # Snow/ice flag:
        gsnow=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
               variables['snow_ice_flag'][:]
        tsnow=gsnow.data[0,:,:]
        # Convert all valid snow/ice free flag values (0,255) to 0.
        # Ocean:
        tsnow=np.where(tsnow==255, 0, tsnow)
        # Coastlines (listed as potential "suspect" in the ATBD document p. 67):
        tsnow=np.where(tsnow==252, 0, tsnow)
        # Less then 1% snow/ice cover:
        tsnow=np.where(tsnow<1, 0, tsnow)
        # Snow/ice misclassified as clouds:
        tsnow=np.where(((tsnow>80)&(tsnow<104)&(tscenep>(0.98*tsurfp))),\
                       0, tsnow)
        # Set clouds over snow/ice scenes to nan:
        cldfrac=np.where(tsnow!=0, np.nan, cldfrac )
        tcldpres=np.where(tsnow!=0, np.nan, tcldpres ) 

    # close KNMI NO2 file:
    fh.close()

    # Cloud input data if dlr-ocra cloud product is used:
    if ( cldprd=='dlr-ocra' ):
        
        # Read data:
        fd=Dataset(cldfiles[f],mode='r')

        # Check that date is the same as the no2 file:
        strdate=cldfiles[f]
        strdate=strdate[-66:-51]
        if strdate!=files[-66:-51]:
            print('NO2 file, Cloud file: '+files[-66:-51]+", "+strdate,\
                  flush=True)
            print('EXITING: Files are not for the same date!',flush=True)
            sys.exit()
        
        # Cloud fraction:
        tcldfrac=fd.groups['PRODUCT'].variables['cloud_fraction'][:]
        cldfrac=tcldfrac.data[0,:,:]
        # Cloud top height (m):
        gcldhgt=fd.groups['PRODUCT'].variables['cloud_top_height'][:]
        tcldhgt=np.ma.getdata(gcldhgt[0,:,:])

        # Define pressure array of zeros:
        tcldpres=np.zeros(tcldhgt.shape)
        
        # Calculate pressure assuming dry atmosphere using external
        # conversion code (convert_height_to_press.py). There's a cloud
        # top pressure entry in the data file, but this is obtained using
        # ECMWF pressure and might have errors. Diego (DLR cloud product PI
        # recommended I use cloud altitude rather than pressure data):
        hgtind=np.where((tcldhgt!=fillval))
        tcldpres[hgtind]=alt2pres(tcldhgt[hgtind])

        # QA value:
        cldqa=fd.groups['PRODUCT'].variables['qa_value'][0,:,:]

        # Snow/ice flag (combined NISE and climatology, so misclassification
        # issues in FRESCO cloud product addressed):
        gsnow=fd.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
               variables['snow_ice_flag'][:]
        tsnow=gsnow.data[0,:,:]

        # Set clouds over snow/ice scenes to nan:
        cldfrac=np.where(tsnow!=0, np.nan, cldfrac )
        tcldpres=np.where(tsnow!=0, np.nan, tcldpres ) 

        # Set poor quality cloud data to nan:
        cldfrac=np.where(cldqa<0.5, np.nan, cldfrac )
        tcldpres=np.where(cldqa<0.5, np.nan, tcldpres ) 
        
        # Check that data shapes are equal:
        if cldfrac.shape!=sza.shape:
            print('Cloud product and NO2 indices ne!',flush=True)
            print(cldfrac.shape,sza.shape,flush=True)
            print('Skipping this swath',flush=True)
            fd.close()
            continue

        # Close DLR CLOUD file:
        fd.close()

    # ======= FILTER FOR RELEVANT DATA  ========

    # Find relevant data only:
    # Filter to only include very cloudy scenes at high altitude
    # No. of valid total NO2 column points before apply filtering:
    inicnt=np.count_nonzero(~np.isnan(tgeototvcd))

    # Filter out scenes with cloud fractions (and hence cloud pressures) that are nan:
    tgeototvcd=np.where(np.isnan(cldfrac), np.nan, tgeototvcd)

    # Filter out scenes with cloud fraction < cloud threshold:
    tgeototvcd=np.where(cldfrac<cldthld, np.nan, tgeototvcd)

    # Filter out scenes with cloud heights outside the UT range of interest (180-450 hPa):
    tgeototvcd=np.where(tcldpres>pmax*1e2, np.nan, tgeototvcd)
    tgeototvcd=np.where(tcldpres<pmin*1e2, np.nan, tgeototvcd)

    # Filter out low quality data (0.45 threshold suggested TROPOMI NO2
    # PI Henk Eskes from the KNMI:
    tgeototvcd=np.where(qaval<0.45, np.nan, tgeototvcd)

    # Filter out scenes with AAI > 1 (could be misclassified as clouds)
    tgeototvcd=np.where(aai>1., np.nan, tgeototvcd)

    # No. of points retained after filtering:
    tcnt=np.count_nonzero(~np.isnan(tgeototvcd))

    # Save % valid points retained to print out average at end of routine:
    postfilt[f]=100.*(tcnt/inicnt)

    # ======= DATA TRIMMING ========

    # Trim the data to include only those relevant:
    # This also reshapes the data from a 2D to a 1D array:
    lons=tlons[~np.isnan(tgeototvcd)]
    lats=tlats[~np.isnan(tgeototvcd)]
    cldpres=tcldpres[~np.isnan(tgeototvcd)]
    stratno2=tstratno2[~np.isnan(tgeototvcd)]
    amf_geo=tamf_geo[~np.isnan(tgeototvcd)]
    geototvcd=tgeototvcd[~np.isnan(tgeototvcd)]

    # ======= DEFINE MONTHLY MEAN OUTPUT ARRAYS ========

    # Define output data for each file:
    gno2=np.zeros((xdim,ydim,nval))
    gstrat=np.zeros((xdim,ydim,nval))
    gcldp=np.zeros((xdim,ydim,nval))
    cntloop=np.zeros((xdim,ydim))

    # ======= CLUSTER THE DATA INTO FINAL GRIDSQUARES ========

    # Data clustering before applying cloud-slicing:
    for w, geototvcd in enumerate(geototvcd):
        if ( geototvcd!=np.nan and geototvcd!=0.0 ):
            
            # Skip over pixels where total column is less than stratospheric
            # column. This addresses positive bias in the cloud-sliced results
            # at low concentrations of UT NO2:
            if ( geototvcd < stratno2[w] ): 
                continue
            
            # Find nearest gridsquare in output grid:
            p,q = np.argmin(abs(out_lon-lons[w])),\
                  np.argmin(abs(out_lat-lats[w]))

            # Convert columns from mol/m2 to molec/cm2:
            tvcdno2=np.multiply(geototvcd,no2sfac) 
            tstrat=np.multiply(stratno2[w],no2sfac) 
           
            # Get relevant data in each output grid square:
            # JR: can you apply the same data append approach as you used for the ut_no2_gc_test.py code?
            gno2[p,q,int(cntloop[p,q])]=tvcdno2
            gstrat[p,q,int(cntloop[p,q])]=tstrat
            gcldp[p,q,int(cntloop[p,q])]=cldpres[w]

            # Increment indices:
            cntloop[p,q]=cntloop[p,q]+1

    # ======= APPLY CLOUD-SLICING RETRIEVAL TO THE CLUSTERS ========

    # Estimate daily mean VMRs from the clustered data:
    for i in range(xdim):
        for j in range(ydim):

            # Only loop over grids with relevant data (identified as 
            # vectors where the first entry is zero):
            if gno2[i,j,0] == 0:
                continue

            # Get NO2 partial columns above cloudy scenes:
            tcolno2=gno2[i,j,:]
            # Trim to remove trailing zeros:
            tcolno2=np.trim_zeros(tcolno2)  
            # Convert NO2 from molec/cm2 to molec/m2 for input to the
            # cloud-slicing algorithm:
            tcolno2=np.multiply(tcolno2,1e4)

            # Skip if fewer than 10 points:
            if len(tcolno2)<10: 
                num[0]=num[0]+1
                continue

            # Get stratospheric columns:
            strat=gstrat[i,j,:]
            # Trim to remove trailing zeros:
            strat=np.trim_zeros(strat)   

            # Get cloud top pressure and convert from Pa to hPa:
            tcld=gcldp[i,j,:]
            # Trim to remove trailing zeros:
            tcld=np.trim_zeros(tcld)    
            # Convert from Pa to hPa for intput to the cloud-slicing algorithm:
            tcld=np.multiply(tcld,1e-2)
                
            # Error check that the sizes of the arrays are equal:
            if ( len(tcld)!=len(tcolno2) ):
                print( 'Array sizes ne: cloud height and partial column')
                print(len(tcld),len(tcolno2))
                print(tcld)
                print(tcolno2)
                print(i,j)
                sys.exit()

            # Skip scenes with non-uniform stratosphere using the
            # same threshold as is used for GEOS-Chem:
            if (np.std(strat)/np.mean(strat))>0.02:
                num[6]=num[6]+1
                continue

            # Get number of points:
            npnts=len(tcld)
            if npnts>maxpnts:
                maxpnts=npnts
                print(maxpnts,flush=True)
                    
            # Use cloud_slice_ut_no2 function to get NO2 mixing
            # ratio from cloud-slicing:
            if ((npnts>=20)&(npnts<100)):
                csval=cldslice(tcolno2,tcld)
                if np.isnan(csval[0]) or np.isnan(csval[1]):
                    num[csval[2]-1]=num[csval[2]-1]+1
                    continue

                # Get mean cloud top pressure:
                pmean=csval[3]
                # Calculate Gaussian weight:
                gwgt=np.exp((-(pmean-315)**2)/(2*135**2))

                # Sum Gaussian weighted mean UT NO2:
                gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                gerr[i,j]=gerr[i,j]+gwgt
                gcnt[i,j]=gcnt[i,j]+1
                cnt=cnt+1

            elif (npnts>=100):

                # Define number of iterations:
                niter=round(npnts/40)
                nloop=list(range(niter))

                # Loop over data:
                for w in nloop:

                    # Cloud-slicing:
                    csval=cldslice(tcolno2[w::niter],tcld[w::niter])

                    # Mean cloud pressure:
                    pmean=csval[3]

                    # Gaussian weight:
                    gwgt=np.exp((-(pmean-315)**2)/(2*135**2)) 

                    if np.isnan(csval[0]) or np.isnan(csval[1]):
                        num[csval[2]-1]=num[csval[2]-1]+1
                        continue
                    else:

                        # Sum Gaussian weighted meann
                        gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                        gerr[i,j]=gerr[i,j]+gwgt
                        gcnt[i,j]=gcnt[i,j]+1
                        cnt=cnt+1

# Print diagnostic information at the end of the routine:
print('Max no. of cloud-sliced retrievals in a grid: ',np.amax(gcnt),flush=True)
print('No. of data lost because: ',flush=True)

print('(1) Too few points: ',num[0],flush=True)
print('(2) Low cloud height range: ',num[1],flush=True)
print('(3) Low cloud height std dev: ',num[2],flush=True)
print('(4) Large error: ',num[3],flush=True)
print('(5) Significantly less then zero: ',num[4],flush=True)
print('(6) Outlier (NO2 > 200 pptv): ',num[5],flush=True)
print('(7) Non-uniform stratosphere: ',num[6],flush=True)
print('(8) Successful retrievals: ',cnt,flush=True)
print('(9) Total possible points: ',(np.sum(num)+cnt),flush=True)

print('Mean % points retained: ', np.mean(postfilt),flush=True)

# Get seasonal means:
gno2vmr=np.divide(gno2vmr, gerr, where=gcnt!=0)
gno2vmr[gcnt==0]=np.nan
gerr=np.divide(gerr, gcnt, where=gcnt!=0)
gerr[gcnt==0]=np.nan
gcnt[gcnt==0]=np.nan

# Save the data to NetCDF:
ncout=Dataset('./test.nc',mode='w',format='NETCDF4') 
#ncout=Dataset('./Data/tropomi-ut-no2-'+cldprd+'-'+strcldthld+'-'+StrRes+'-'+\
#              seas+'-'+yrrange+'-v7.nc',mode='w',format='NETCDF4') 
 
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
cs=m.pcolor(xi,yi,np.squeeze(gno2vmr), vmin=0, vmax=100, cmap='jet')
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
cs=m.pcolor(xi,yi,np.squeeze(gcnt), vmin=0., vmax=30, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Number of points')

plt.show()

# Close the log file:
log.close()
