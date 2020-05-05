#!/usr/bin/python

''' Read and regrid cloud information from the official TROPOMI CLOUD product 
    and from the FRESCO cloud product in the TROPOMI NO2 data product file.

    Data are regridded to a 2x2.5 degree global grid and saved as monthly 
    means in NetCDF format for scenes with TROPOMI FRESCO cloud fraction >
    0.7 and clouds top (centroid) pressure in the upper troposphere 
    (180-450 hPa).

    (1) Initial version created on 15 March 2020 by Eloise Marais
        (eloise.marais@le.ac.uk).

    '''

# Import relevant packages:
import glob
import sys
import os
import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
from scipy import stats
import leastsq
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
#import random
from gcpy.gcpy.plot.colormap import WhGrYlRd

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# Define output resolution:
outres='1x1'

# Define month of interest as string 2 characters in length:
StrMM='01'

#DEFINE GRID:
# Define model grid:
minlat=-88.
maxlat=88.
minlon=-178.
maxlon=178.
if (outres=='1x1'): dellat,dellon=1,1
if (outres=='2x25'): dellat,dellon=2,2.5
if (outres=='4x5'): dellat,dellon=4,5
out_lon=np.arange(minlon,maxlon,dellon)#(-180,181,dellon)
out_lat=np.arange(minlat,maxlat,dellat)#(-90.,91.,dellat)
# Convert output lats and longs to 2D:
X, Y = np.meshgrid(out_lon, out_lat,indexing='ij') 
# Dimensions of output data:
xdim=len(out_lon)
ydim=len(out_lat)

# DEFINE DATA ARRAYS for comparing the two datasets where KNMI cloud
# fraction exceeds 0.7 and cloud top height is between 180 and 450 hPa:
# Define output arrays:
gknmi_cf=np.zeros(X.shape)
gknmi_ct=np.zeros(X.shape)
gknmi_cb=np.zeros(X.shape)
gknmi_cnt=np.zeros(X.shape)
gdlr_cf=np.zeros(X.shape)
gdlr_ct=np.zeros(X.shape)
gdlr_cb=np.zeros(X.shape)
gdlr_od=np.zeros(X.shape)
gdlr_cnt=np.zeros(X.shape)
#gdiff_cf=np.zeros(X.shape)

# Define data arrays to output cloud fraction binned into 0.05 
# increments from 0.7 to 1.0 centred at 0.725, 0.775 etc.:
# Consider also splitting this into latitude bands.
cldbin=np.arange(0.7,1.0,0.05)
latbin=np.arange(-90,90,15)
knmi_cf_freq=np.zeros((len(latbin),len(cldbin)))
dlr_cf_freq=np.zeros((len(latbin),len(cldbin)))

# Coverage count:
nobs_fresco=0
nobs_dlr=0

# Input parameter (to selet month of interest):
# Define string of year and first 3 letters of month name based on above entry:
if StrMM=='01': StrYY,MMName='2020','jan'
if StrMM=='02': StrYY,MMName='2020','feb'
if StrMM=='03': StrYY,MMName='2020','mar'
if StrMM=='04': StrYY,MMName='2020','apr'
if StrMM=='05': StrYY,MMName='2020','may'
if StrMM=='06': StrYY,MMName='2019','jun'
if StrMM=='07': StrYY,MMName='2019','jul'
if StrMM=='08': StrYY,MMName='2019','aug'
if StrMM=='09': StrYY,MMName='2019','sep'
if StrMM=='10': StrYY,MMName='2019','oct'
if StrMM=='11': StrYY,MMName='2019','nov'
if StrMM=='12': StrYY,MMName='2019','dec'

# Define directory:
s5pdir='/data/uptrop/nobackup/tropomi/Data/'

# Define days:
ndays=31
dd=list(range(1,ndays+1))
strdd=[''] * ndays
cnt=0
for d in dd:
    strdd[cnt]='0'+str(d) if d<10 else str(d)
    cnt=cnt+1

# Get DLR data file names:
tdfile=glob.glob(s5pdir+'CLOUD_OFFL/'+StrYY+'/'+StrMM+'/S5P_OFFL_L2__CLOUD__'+\
                 StrYY+StrMM+'*')
tdfile=sorted(tdfile)

# Get FRESCO file names:
tffile=glob.glob(s5pdir+'NO2_OFFL/'+StrYY+'/'+StrMM+'/S5P_OFFL_L2__NO2____'+\
                 StrYY+StrMM+'*')
tffile=sorted(tffile)

# Check that number of files are equal. If not, exit the programme:
if len(tdfile) != len(tffile):
    print('DLR files = ',len(tdfile))
    print('FRESCO files = ',len(tffile))
    print('unequal number of files')    
    sys.exit()

# Get number of files:
nfiles=len(tdfile)

# Define fill value:
fillval=9.96921e+36

# NOTE: Cloud pressure is given in Pa!!!

# Loop over files:
for f in range(nfiles):

    # Track progress:
    print('===> Processing: ', tdfile[f])

    # Get orbit/swath number:
    file0=tffile[f]
    forb=file0[104:109]
    file1=tdfile[f]
    dorb=file1[106:111]
    
    # Check orbit/swath is the same. If not, skip this iteration:
    if ( forb!=dorb ): continue

    # Read DLR cloud data:
    fh=Dataset(tdfile[f],mode='r')

    # Extract data of interest:
    tdlons=fh.groups['PRODUCT'].variables['longitude'][:][0,:,:]
    tdlats=fh.groups['PRODUCT'].variables['latitude'][:][0,:,:]
    tfrc=fh.groups['PRODUCT'].variables['cloud_fraction'][:]
    tdfrc=tfrc.data[0,:,:]
    ttop=fh.groups['PRODUCT'].variables['cloud_top_pressure'][:]
    tdtop=ttop.data[0,:,:]
    tbase=fh.groups['PRODUCT'].variables['cloud_base_pressure'][:]
    tdbase=tbase.data[0,:,:]
    tqval=fh.groups['PRODUCT'].variables['qa_value'][:]
    tdqval=tqval.data[0,:,:]  # for more accurate CAL product
    toptd=fh.groups['PRODUCT'].variables['cloud_optical_thickness'][:]
    tdoptd=toptd.data[0,:,:]
    tsnow=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['snow_ice_flag_nise'][:]
    tdsnow=tsnow.data[0,:,:]    
    # Convert all valid snow/ice free flag values (0,255) to 0.
    # Ocean:
    tdsnow=np.where(tdsnow==255, 0, tdsnow)
    # Coastlines (listed as potentially "suspect" in the ATBD document p. 67):
    tdsnow=np.where(tdsnow==252, 0, tdsnow)
    # Get data indices:
    md,nd = tdlons.shape

    # Count data:
    dlr_ind=np.where((tdqval<0.5)&(tdfrc>=0.7)&(tdtop>=18000)&\
                     (tdtop<=45000)&(tdsnow==0))[0]
    nobs_dlr=nobs_dlr+len(dlr_ind)

    # Set missing/poor quality/irrelevant data to NAN:
    # Apply cloud fraction filter, but set it to 0.7 rather
    #     than 0.9 to account for variability around this threshold
    #     in the two cloud products.
    #inicnt=np.count_nonzero(~np.isnan(tdfrc))
    #print(inicnt)
    #tdfrc=np.where(tdfrc<0.7, np.nan, tdfrc)
    # Fill values (do for cloud fraction and cloud top pressure, as missing
    # values for cloud fraction may not be the same as missing values for
    # other data, as these are obtained with different algorithms):
    # necessarily the same as the fill values for the other data):
    tdfrc=np.where(tdfrc==fillval, np.nan, tdfrc)
    tdtop=np.where(tdtop==fillval, np.nan, tdtop)

    # Set cloud fraciton to nan for scenes with missing cloud top
    # pressure data:
    tdfrc=np.where(np.isnan(tdtop), np.nan, tdfrc)    

    # Poor data quality:
    tdfrc=np.where(tdqval<0.5, np.nan, tdfrc)

    # Snow/ice cover:
    tdfrc=np.where(tdsnow!=0, np.nan, tdfrc)

    # Apply filter to remaining data:
    tdtop=np.where(np.isnan(tdfrc), np.nan, tdtop)
    tdbase=np.where(np.isnan(tdfrc), np.nan, tdbase)
    tdoptd=np.where(np.isnan(tdfrc), np.nan, tdoptd)

    # Error check:
    if ( np.nanmax(tdtop)==fillval ):
        print('Not all missing values converted to NANs')
        sys.exit()
    if ( np.nanmax(tdbase)==fillval ):
        print('Not all missing values converted to NANs')
        sys.exit()
    if ( np.nanmax(tdoptd)==fillval ):
        print('Not all missing values converted to NANs')
        sys.exit()

    # close DLR file:
    fh.close()

    # Read in FRESCO cloud data:
    fh=Dataset(tffile[f],mode='r')

    # Extract data of interest:
    tflons=fh.groups['PRODUCT'].variables['longitude'][:][0,:,:]
    tflats=fh.groups['PRODUCT'].variables['latitude'][:][0,:,:]
    tfrc=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['cloud_fraction_crb'][:]
    tffrc=tfrc.data[0,:,:]
    talb=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['cloud_albedo_crb'][:]
    tfalb=talb.data[0,:,:]
    ttop=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['cloud_pressure_crb'][:]
    tftop=ttop.data[0,:,:]
    tqval=fh.groups['PRODUCT'].variables['qa_value'][:]
    tfqval=tqval.data[0,:,:] 
    tsnow=fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA'].\
          variables['snow_ice_flag'][:]
    tfsnow=tsnow.data[0,:,:]
    # Convert all valid snow/ice free flag values (0,255) to 0.
    # Ocean:
    tfsnow=np.where(tfsnow==255, 0, tfsnow)
    # Coastlines (listed as potential "suspect" in the ATBD document p. 67):
    tfsnow=np.where(tfsnow==252, 0, tfsnow)
    # Get data indices:
    mf,nf=tflons.shape

    # close FRESCO file:
    fh.close()

    # Count data:
    fr_ind=np.where((tfqval<0.45)&(tffrc>=0.7)&(tftop>=18000)&\
                    (tftop<=45000)&(tfsnow==0))[0]
    nobs_fresco=nobs_fresco+len(fr_ind)

    m=md
    n=nd
    # Skip files if the number of indices are not equal:
    if mf!=md:
        print('Indices ne for ', tdfile[f])
        continue
        #m=min(md,mf)
        #n=min(nd,nf)

    # Set missing/poor quality/irrelevant data to NAN:
    # Apply cloud fraction filter, but set it to 0.7 rather
    #     than 0.9 to account for variability around this threshold
    #     in the two cloud products.
    #inicnt=np.count_nonzero(~np.isnan(tffrc))
    #print(inicnt)
    tffrc=np.where(tffrc<0.7, np.nan, tffrc)
    tffrc=np.where(tffrc==fillval, np.nan, tffrc)
    # QA Flags. Threshold of 0.45 suggested by Henk Eskes in email 
    # exchange on 18 Jan 2020:
    tffrc=np.where(tfqval<0.45, np.nan, tffrc)
    # Apply cloud top pressure filter to only consider clouds above
    #   500 hPa and below 150 hPa (more generous than the 450-200 hPa
    #   range to account for variability around this threshold in the
    #   two cloud products:
    tffrc=np.where(tftop>450*1e2, np.nan, tffrc)
    tffrc=np.where(tftop<180*1e2, np.nan, tffrc)
    # Snow/ice cover:
    tffrc=np.where(tfsnow!=0, np.nan, tffrc)

    # Apply filter to remaining data:
    tftop=np.where(np.isnan(tffrc), np.nan, tftop)

    # Gather data on frequency of cloud fraction > 0.7:
    # This is done before filtering for scenes with knmi cloud frac
    # > 0.7 to also include all relevant DLR scenes:
    # loop over cloud and latitude band bins:
    for w in range(len(cldbin)):
        for n in range(len(latbin)):
             
            # (1) KNMI:
            # Get indices for relevant data:
            fbin=np.where((tffrc>=(cldbin[w]-0.025))&\
                          (tffrc<(cldbin[w]+0.025))&\
                          (tflats>=(latbin[n]-7.5))&(tflats<(latbin[n]+7.5))&\
                          (tftop>=18000)&(tftop<=45000)&(~np.isnan(tffrc)))[0]
            if len(fbin)>0:
                knmi_cf_freq[n,w]=knmi_cf_freq[n,w]+len(fbin)

            # (2) DLR-OCRA:
            # Get indices for relevant data:
            dbin=np.where((tdfrc>=(cldbin[w]-0.025))&\
                          (tdfrc<(cldbin[w]+0.025))&\
                          (tdlats>=(latbin[n]-7.5))&(tdlats<(latbin[n]+7.5))&\
                          (tdtop>=18000)&(tdtop<=45000)&(~np.isnan(tdfrc)))[0]
            if len(dbin)>0:
                dlr_cf_freq[n,w]=dlr_cf_freq[n,w]+len(dbin)

    # REGRID THE DATA:
    for i in range(m):
        for j in range(n):

            # Skip where FRESCO cloud product is NAN:
            if (np.isnan(tffrc[i,j]) or np.isnan(tdfrc[i,j])): continue

            # Error checks:
            # Check that lat and lon values are the same (except for
            # those at 180/-180 that can have different signs):
            if (tdlons[i,j] != tflons[i,j]) and \
               (abs(tdlons[i,j])!=180.0):
                print('Longitudes not the same')
                print(tdlons[i,j],tflons[i,j])
                p,q = np.argmin(abs(out_lon-tdlons[i,j])),\
                      np.argmin(abs(out_lat-tdlats[i,j]))
                print(p)
                print(out_lon[p])
                sys.exit()
            if tdlats[i,j] != tflats[i,j]:
                print('Latitudes not the same')
                sys.exit()           

            # Find corresponding gridsquare:
            p,q = np.argmin(abs(out_lon-tdlons[i,j])),\
                  np.argmin(abs(out_lat-tdlats[i,j]))            

            # Add data to output arrays:
            #if np.isnan(tffrc[i,j]): continue
            gknmi_cf[p,q]=gknmi_cf[p,q]+tffrc[i,j]
            gknmi_ct[p,q]=gknmi_ct[p,q]+np.divide(tftop[i,j],100)
            gknmi_cnt[p,q]=gknmi_cnt[p,q]+1.0
            #if np.isnan(tdfrc[i,j]): continue
            gdlr_cf[p,q]=gdlr_cf[p,q]+tdfrc[i,j]
            gdlr_ct[p,q]=gdlr_ct[p,q]+np.divide(tdtop[i,j],100)
            gdlr_cb[p,q]=gdlr_cb[p,q]+np.divide(tdbase[i,j],100)
            gdlr_od[p,q]=gdlr_od[p,q]+tdoptd[i,j]
            gdlr_cnt[p,q]=gdlr_cnt[p,q]+1.0

# Get means and differences (only means for cloud base height):
# (1) Cloud fraction:
gknmi_cf=np.where(gknmi_cnt==0, np.nan, np.divide(gknmi_cf,gknmi_cnt))
gdlr_cf=np.where(gdlr_cnt==0, np.nan, np.divide(gdlr_cf,gdlr_cnt))
gdiff_cf=np.subtract(gknmi_cf,gdlr_cf)
                  
# (2) Cloud top pressure:
gknmi_ct=np.where(gknmi_cnt==0, np.nan, np.divide(gknmi_ct,gknmi_cnt))
gdlr_ct=np.where(gdlr_cnt==0, np.nan, np.divide(gdlr_ct,gdlr_cnt))
gdiff_ct=np.subtract(gknmi_ct,gdlr_ct)

# (3) Cloud albedo/optical depth:
gdlr_od=np.where(gdlr_cnt==0, np.nan, np.divide(gdlr_od,gdlr_cnt))

# (4) Cloud base pressure (DLR only):
gdlr_cb=np.divide(gdlr_cb,gdlr_cnt)
gdlr_cb[gdlr_cnt==0.0]=np.nan    

# Print number of observations to screen:
print('No. of FRESCO obs for '+MMName+' = ',nobs_fresco)
print('No. of DLR obs for '+MMName+' = ',nobs_dlr)

# Write data to file:
outdir='/data/uptrop/Projects/UpTrop/python/Data/'
outfile=outdir+'fresco-dlr-cloud-products-'+MMName+'-'+StrYY+'-'+\
        outres+'-v1.nc'
ncfile=Dataset(outfile, 'w', format='NETCDF4_CLASSIC' )
ncfile.createDimension('xdim',len(out_lon))
ncfile.createDimension('ydim',len(out_lat))
ncfile.createDimension('frdim',len(cldbin))
ncfile.createDimension('lbdim',len(latbin))

# Global attributes:
ncfile.title='FRESCO and DLR TROPOMI monthly mean cloud properties for '+\
    StrMM+' '+StrYY
ncfile.subtitle='Data written to file in '+StrMM+' '+StrYY
ncfile.anything='Verions used are v010107 for the DLR product and v010302 for the FRESCO product'

# Longitudes:
lon = ncfile.createVariable('lon',np.dtype('float'),('xdim',))
lon.units='degrees_east'
lon.long_name='longitude'
lon[:] = out_lon

# Latitudes:
lat = ncfile.createVariable('lat',np.dtype('float'),('ydim',))
lat.units='degrees_north'
lat.long_name='latitude'
lat[:] = out_lat

# Cloud fraction bins:
cbin = ncfile.createVariable('cld_fr_bin',np.dtype('float'),('frdim',))
cbin.units='unitless'
cbin.long_name='Cloud fraction bin centre'
cbin[:] = cldbin

# Latitude bands bins:
lbbin = ncfile.createVariable('lat_band_bin',np.dtype('float'),('lbdim',))
lbbin.units='unitless'
lbbin.long_name='Latitude band bin centre'
lbbin[:] = latbin

# DLR cloud fraction:
data1 = ncfile.createVariable('dlr_cld_frac',np.dtype('float'),('xdim','ydim'))
data1.units='unitless'
data1.long_name='DLR CAL cloud fraction'
data1[:] = gdlr_cf

# DLR cloud top pressure:
data2 = ncfile.createVariable('dlr_cld_top_pres',np.dtype('float'),('xdim','ydim'))
data2.units='hPa'
data2.long_name='DLR CAL cloud top pressure'
data2[:] = gdlr_ct

# DLR cloud base pressure:
data3 = ncfile.createVariable('dlr_cld_base_pres',np.dtype('float'),('xdim','ydim'))
data3.units='hPa'
data3.long_name='DLR CAL cloud base pressure'
data3[:] = gdlr_cb

# DLR cloud optical depth:
data4 = ncfile.createVariable('dlr_cld_optd',np.dtype('float'),('xdim','ydim'))
data4.units='unitless'
data4.long_name='DLR CAL cloud optical thickness'
data4[:] = gdlr_od

# FRESCO cloud fraction:
data5 = ncfile.createVariable('knmi_cld_frac',np.dtype('float'),('xdim','ydim'))
data5.units='unitless'
data5.long_name='FRESCO cloud fraction'
data5[:] = gknmi_cf

# Frequency of FRESCO cloud fractions > 0.7:
data6 = ncfile.createVariable('knmi_cf_freq',np.dtype('float'),\
                              ('lbdim','frdim',))
data6.units='unitless'
data6.long_name='Freqeuncy distribution of FRESCO cloud fractions > 0.7'
data6[:] = knmi_cf_freq

# FRESCO cloud top pressure:
data7 = ncfile.createVariable('knmi_cld_top_pres',np.dtype('float'),('xdim','ydim'))
data7.units='hPa'
data7.long_name='FRESCO cloud top pressure'
data7[:] = gknmi_ct

# FRESCO data points in each gridsquare:
data8 = ncfile.createVariable('knmi_points',np.dtype('float'),('xdim','ydim'))
data8.units='unitless'
data8.long_name='Number of FRESCO data points'
data8[:] = gknmi_cnt

# DLR data points in each gridsquare:
data9 = ncfile.createVariable('dlr_points',np.dtype('float'),('xdim','ydim'))
data9.units='unitless'
data9.long_name='Number of DLR data points'
data9[:] = gdlr_cnt

# Frequency of DRL-OCRA cloud fractions > 0.7:
data10 = ncfile.createVariable('dlr_cf_freq',np.dtype('float'),('lbdim','frdim'))
data10.units='unitless'
data10.long_name='Freqeuncy distribution of DLR cloud fractions > 0.7'
data10[:] = dlr_cf_freq

# close the file.
ncfile.close()

# PLOT THE DATA:
m=Basemap(resolution='l',projection='merc',lat_0=0,lon_0=0,\
          llcrnrlon=minlon,llcrnrlat=-70,\
          urcrnrlon=maxlon,urcrnrlat=70)
xi,yi=m(X,Y)

# (1) Cloud fraction:
plt.figure(1)
plt.subplot(3, 1, 1)
cs=m.pcolor(xi,yi,np.squeeze(gknmi_cf), vmin=0.7, vmax=1.3, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('FRESCO Cloud fraction')

plt.subplot(3, 1, 2)
cs=m.pcolor(xi,yi,np.squeeze(gdlr_cf), vmin=0.7, vmax=1.3, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('DLR Cloud fraction')

plt.subplot(3, 1, 3)
cs=m.pcolor(xi,yi,np.squeeze(gdiff_cf), vmin=-0.2, vmax=0.2, cmap='bwr')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Difference (FRESCO-DLR)')

plt.savefig('./Images/fresco-vs-dlr-cloud-frac-'+MMName+'-'+StrYY+'-v2.ps', \
            format='ps')

# (2) Cloud top pressure:
plt.figure(2)
plt.subplot(3, 1, 1)
cs=m.pcolor(xi,yi,np.squeeze(gknmi_ct), vmin=150, vmax=500, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('FRESCO Cloud top pressure [hPa]')

plt.subplot(3, 1, 2)
cs=m.pcolor(xi,yi,np.squeeze(gdlr_ct), vmin=150, vmax=500, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('DLR Cloud top pressure [hPa]')

plt.subplot(3, 1, 3)
cs=m.pcolor(xi,yi,np.squeeze(gdiff_ct), vmin=-30, vmax=30, cmap='bwr')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Difference (FRESCO-DLR)')

plt.savefig('./Images/fresco-vs-dlr-cloud-top-press-'+MMName+'-'+StrYY+\
            '-v1.ps', format='ps')

# (3) Cloud optical depth/albedo (DLR only):
plt.figure(3)
cs=m.pcolor(xi,yi,np.squeeze(gdlr_od), vmin=0, vmax=200, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('DLR cloud optical thickness')

plt.savefig('./Images/dlr-cloud-optical-depth-'+MMName+'-'+StrYY+'-v1.ps', \
            format='ps')

# (4) Cloud base pressure (DLR only):
plt.figure(4)
cs=m.pcolor(xi,yi,np.squeeze(gdlr_cb), vmin=150, vmax=500, cmap='jet')
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('DLR base pressure [hPa]')

plt.savefig('./Images/dlr-cloud-base-press-'+MMName+'-'+StrYY+'-v1.ps', \
            format='ps')

# (5) Number of points (both):
plt.figure(4)
plt.subplot(2, 1, 1)
cs=m.pcolor(xi,yi,np.squeeze(gknmi_cnt), vmin=0, vmax=500, cmap=WhGrYlRd)
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('FRESCO No. of obs')

plt.subplot(2, 1, 2)
cs=m.pcolor(xi,yi,np.squeeze(gdlr_cnt), vmin=0, vmax=500, cmap=WhGrYlRd)
m.drawcoastlines()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('DLR No. of obs')

plt.savefig('./Images/fresco-dlr-number-of-obs-'+MMName+'-'+StrYY+'-v1.ps', \
            format='ps')

plt.show()

