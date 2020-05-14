#!/usr/bin/python

''' Use synthetic partial columns from GEOS-Chem to obtain cloud-sliced 
    NO2 in the upper troposphere and compare this to UT NO2 obtained if
    simply average the NO2 mixing ratios from the model over the same 
    pressure range (the "truth"). Both are obtained by Gaussian 
    weighting toward the pressure center.

    GEOS-Chem partial columns are obtained over Europe, North America, 
    and China at the GEOS-FP meteorology native resolution (0.25x0.3125)
    (latxlon) for June-August 2016-2017.

    Input options to process the data include the region, the horizontal 
    resolution, and model simulation years.
    '''

# Import relevant packages:
import glob
import sys
import os
import numpy as np
import netCDF4 as nc4
from netCDF4 import Dataset
from scipy import stats
from bootstrap import rma
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits import basemap
from cloud_slice_ut_no2 import cldslice
import leastsq

from gcpy.gcpy.constants import AVOGADRO as na
from gcpy.gcpy.constants import G as g
from gcpy.gcpy.constants import MW_AIR as mmair

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# Decide on region:
Reg='CH' #'NA', 'EU', or 'CH'

# Define target grid:
if Reg=='NA': 
    minlat=6. 
    maxlat=58. 
    minlon=-135. 
    maxlon=-60. 
    dirreg='na'
    factor=40   # used to determine size of gno2
if Reg=='EU':   
    minlat=30.  
    maxlat=62.
    minlon=-20.
    maxlon=40.
    dirreg='eu_naei'
    factor=30
if Reg=='CH':
    minlat=10.
    maxlat=54. 
    minlon=65. 
    maxlon=135. 
    dirreg='ch' 
    factor=100 

# Define information for grid:
StrRes='4x5'   

# Define years of interest:
stryr=['2016','2017']
if len(stryr)==1: yrrange=stryr
if len(stryr)==2: yrrange='2016-2017'

# Name of log file to output code prorgess:
log=open("log_"+Reg+"_"+StrRes+"_"+yrrange+"_v3", "w")
sys.stdout=log

# Define grid information:
if StrRes=='8x10': dellat,dellon=8,10
if StrRes=='4x5': dellat,dellon=4,5
if StrRes=='2x25': dellat,dellon=2,2.5
if StrRes=='1x1': dellat,dellon=1,1
out_lon=np.arange(minlon,maxlon+dellon,dellon)#(-180,181,dellon)
out_lat=np.arange(minlat,maxlat+dellat,dellat)#(-90.,91.,dellat)
# Convert output lats and long to 2D:
X, Y = np.meshgrid(out_lon, out_lat,indexing='ij') 
# Dimensions of output data:
xdim=len(out_lon)
ydim=len(out_lat)
nval=int(factor*dellat*dellon)  #Scale array size with resolution

# Define output arrays:
gno2vcd=np.zeros(X.shape)
gno2vmr=np.zeros(X.shape)
gcldfr=np.zeros(X.shape)
gcld=np.zeros(X.shape)
gerr=np.zeros(X.shape)
trueno2=np.zeros(X.shape)  # "true" cloudy UT NO2
gcutno2=np.zeros(X.shape)  # "true" all-sky UT NO2
truewgt=np.zeros(X.shape)  # "true" all-sky UT NO2 Gaussian weights
gascnt=np.zeros(X.shape)   # Count all-sky
gcnt=np.zeros(X.shape)
num=np.zeros(7)

#Initialize:
cnt=0
#numcnt=0
maxcnt=0

# Define pressure range:
pmin=180     #380   #290   #180   #190
pmax=450     #450   #380   #290   #440

# Define string to represent the layer range:
prange=str(pmin)+'-'+str(pmax)

#Define factor to convert slope of NO2 mixing ratio versus pressure
# to VMR:
den2mr=np.divide((np.multiply(g,mmair)),na)

# Get files (test June 2016 for now)
# 2016:
gcdir='/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_'+dirreg+'_iccw/'
files=glob.glob(gcdir+'nc_sat_files_47L/ts_12_15.'+Reg+'.'+stryr[0]+'06*')
#mon=["06","07","08","09"]
mon=["07","08"]
for i in mon:
    for filename in glob.glob(gcdir+'nc_sat_files_47L/ts_12_15.'+Reg+\
                              '.'+stryr[0]+i+'*'):
        files.append(filename)
#2017:
mon=["06","07","08"]
for i in mon:
    for filename in glob.glob(gcdir+'nc_sat_files_47L/ts_12_15.'+Reg+\
                              '.'+stryr[1]+i+'*'):
        files.append(filename)

#mon=["06","07","08","09"]
#mon=["07","08"]
#for i in mon:
#    for filename in glob.glob(gcdir+'nc_sat_files/ts_12_15.'+Reg+\
#                              '.'+stryr[1]+i+'*'):
#        files.append(filename)

# Order the files:
files=sorted(files)

# Check:
print('Number of files:',len(files),flush=True)
#log.write('Number of files:'+str(len(files))+'\n')

# Loop over files:
for f in files:
    print(f,flush=True)

    # Read dataset:
    fh=Dataset(f,mode='r')

    # Extract data of interest:
    # (Add tropopause height to this in the future)
    tlon,tlat,tgcno2,tcldfr,tcldhgt,tadn,tbxhgt,tpedge=\
          fh.variables['LON'],fh.variables['LAT'],\
          fh.variables['IJ-AVG-S__NO2'],fh.variables['TIME-SER__CF'],\
          fh.variables['TIME-SER__CThgt'],fh.variables['TIME-SER__AIRDEN'],\
          fh.variables['BXHGHT-S__BXHEIGHT'],fh.variables['PEDGE-S__PSURF']
    tlon=tlon[:]
    tlat=tlat[:]
    tgcno2=tgcno2[:]
    tcldfr=tcldfr[:]
    tcldhgt=tcldhgt[0,:,:]
    tadn=tadn[:]   # in molec/cm3
    tbxhgt=tbxhgt[:]
    tpedge=tpedge[:]

    # Convert box height from m to cm:
    tbxhgt=tbxhgt*1e2

    # Define output data for this day:
    gno2=np.zeros((xdim,ydim,nval))
    allcldfr=np.zeros((xdim,ydim,nval))
    stratno2=np.zeros((xdim,ydim,nval))
    gcldp=np.zeros((xdim,ydim,nval))
    gtrueno2=np.zeros((xdim,ydim,nval))
    cntloop=np.zeros((xdim,ydim))
    
    # Get column values:
    for y in range(len(tlat)):
        for x in range(len(tlon)):

            # Find nearest gridsquare in output grid:
            p,q = np.argmin(abs(out_lon-tlon[x])),\
                  np.argmin(abs(out_lat-tlat[y]))

            # Calculate corresponding mid-pressure values:
            tpmid=np.zeros(len(tpedge[:,y,x]))
            # Get mid-pressure values, except for highest layer:
            for k in range(len(tpedge[:,y,x])-1):
                tpmid[k]=np.multiply(0.5,(tpedge[k,y,x]+tpedge[k+1,y,x]))
            # Estimate mid-pressure for highest value (doesn't need to
            # be accurate, as surpassing the range of interset):
            tpmid[46]=np.multiply(0.5,(tpedge[46,y,x]+(tpedge[46,y,x]-0.1)))

            # Get indices that fall between 450 and 180 hPa for estimating
            # "true' all-sky UT NO2 and partial columns:
            lind=np.where((tpmid>=pmin)&(tpmid<=pmax))[0]

            # Get Guassian weights that allocated higher weights to points
            # closest to the pressure centre (315 hPa):
            # Equation is:
            #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
            #         135 hPa is the standard deviation.
            twgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(tpmid[lind],315.))),(np.multiply(2,np.square(135.))))))

            # Get UT NO2 under "true" all-sky conditions:
            gcutno2[p,q]=gcutno2[p,q] + \
                          np.sum(tgcno2[lind,y,x]*twgt*1e3)
            gascnt[p,q]=gascnt[p,q]+1.0
            truewgt[p,q]=truewgt[p,q]+np.sum(twgt)

            # Find where cloud fraction in UT exceeds 0.7 after calculating
            # true all-sky NO2:
            # (Keep for testing effect of thick clouds on cloud-sliced UT NO2):
            #if (np.sum(tcldfr[lind,y,x])<=0.7): continue

            # Get model level of cloud top height closest to lowest 
            # pressure-altitude of interest (pmin):
            lcld = np.argmin(abs(tcldhgt[y,x]-tpmid))

            # Skip if cloud top height ouside pressure range of interest:
            lmin,lmax=np.amin(lind),np.amax(lind)
            if ((lcld<lmin)or(lcld>lmax)): continue

            # This error check is probably redundant, but included in case:
            if lcld==0: 
                print("No cloud detected!!!",flush=True)
                sys.exit()

            # Get partial NO2 column in molec/m2 from cloud top height 
            # to highest model level (output up to level 47):
            no22d=np.sum(tgcno2[lmin:,y,x]*1e-5*tadn[lmin:,y,x]*\
                         tbxhgt[lmin:,y,x])

            # Get stratospheric column from 180 hPa aloft 
            # (In future, update to find the model tropopause layer):
            stratind=np.where(tpmid<180.)[0]
            stratcol=np.sum(tgcno2[stratind,y,x]*1e-5*tadn[stratind,y,x]*\
                            tbxhgt[stratind,y,x])

            # Add relevant data to array:
            gno2[p,q,int(cntloop[p,q])]=no22d
            stratno2[p,q,int(cntloop[p,q])]=stratcol
            gcldp[p,q,int(cntloop[p,q])]=tcldhgt[y,x]
            gtrueno2[p,q,int(cntloop[p,q])]=np.mean(tgcno2[lmin:lmax+1,y,x])
            allcldfr[p,q,int(cntloop[p,q])]=np.sum(tcldfr[lmin:lmax+1,y,x])
            
            # Increment indices:
            cntloop[p,q]=cntloop[p,q]+1

    for i in range(xdim):
        for j in range(ydim):
            if gno2[i,j,0]!=0:

                # Define vectors of relevant data:
                tcolno2=gno2[i,j,:]
                tstratno2=stratno2[i,j,:]
                tfrc=allcldfr[i,j,:]
                tcld=gcldp[i,j,:]
                tmrno2=gtrueno2[i,j,:]

                # Skip if fewer than 10 points:
                if len(tcolno2)<10: 
                    num[0]=num[0]+1
                    continue                

                # Trim to remove trailing zeros:
                tcolno2=tcolno2[np.nonzero(tcolno2)]
                tstratno2=tstratno2[np.nonzero(tstratno2)]
                tfrc=tfrc[np.nonzero(tfrc)]
                tcld=tcld[np.nonzero(tcld)]
                tmrno2=tmrno2[np.nonzero(tmrno2)]

                # Remove non-uniform stratosphere:
                if (np.std(tstratno2)/np.mean(tstratno2))>0.02:
                    num[6]=num[6]+1
                    continue

                # Check if there are any repeat values:
                u, ind = np.unique(tcld, return_inverse=True)
                if (len(ind)!=len(tcld)):
                    print('Repeat cloud values',flush=True)
                    sys.exit()

                # Get number of points:
                npnts=len(tcld)
                if npnts>maxcnt:
                    maxcnt=npnts
                    print(maxcnt,flush=True)

                # Get mean cloud top pressure:
                pmean=np.mean(tcld)

                # Calculate Gaussian weight:
                gwgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(pmean,315.))),(np.multiply(2,np.square(135.))))))

                # Use cloud_slice_ut_no2 function to get cloud-sliced 
                # UT NO2 mixing ratios.
                # Treat data differently depending on whether there are no
                # more than or more than 30 points:
                if npnts<=30:
                    csval=cldslice(tcolno2,tcld)
                    if np.isnan(csval[0]) or np.isnan(csval[1]):
                        num[csval[2]-1]=num[csval[2]-1]+1
                        continue

                    # Apply Gaussian weights to cloud-sliced UT NO2:
                    gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                    gcldfr[i,j]=gcldfr[i,j]+np.mean(tfrc)
                    gerr[i,j]=gerr[i,j]+gwgt
               
                    # "True" cloudy UT NO2 (converted from ppbv to pptv):
                    trueno2[i,j]=trueno2[i,j]+\
                                  np.multiply(np.mean(tmrno2*1e3),gwgt)
                    gcnt[i,j]=gcnt[i,j]+1
                    cnt=cnt+1

                elif npnts>30:
                    # Code error check:
                    if len(tcld)!=len(tmrno2): 
                        print('arrays not =',flush=True)
                        print(len(tcld),len(tmrno2),flush=True)

                    # Determine the number of iterations:
                    niter=round(npnts/17)
                    nloop=list(range(niter))

                    # Loop over iterations:
                    for w in nloop:
                        if w==0:
                            csval=cldslice(tcolno2[::niter],tcld[::niter])
                        else:
                            csval=cldslice(tcolno2[w::niter],tcld[w::niter])

                    # Skip if approach didn't work (i.e. cloud-sliced
                    # UT NO2 is NaN):
                    if np.isnan(csval[0]) or np.isnan(csval[1]):
                        num[csval[2]-1]=num[csval[2]-1]+1
                        continue

                    # Gaussian-weighted mean:
                    gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                    gerr[i,j]=gerr[i,j]+gwgt

                    # "True" cloudy UT NO2 (converted from ppbv to pptv):
                    if ( w==0 ):
                        trueno2[i,j]=trueno2[i,j]+\
                            np.multiply(np.mean(tmrno2[::niter]*1e3),gwgt)
                        gcldfr[i,j]=gcldfr[i,j]+np.mean(tfrc[::niter])
                    else:
                        trueno2[i,j]=trueno2[i,j]+\
                            np.multiply(np.mean(tmrno2[w::niter]*1e3),gwgt)
                        gcldfr[i,j]=gcldfr[i,j]+np.mean(tfrc[w::niter])

                    gcnt[i,j]=gcnt[i,j]+1
                    cnt=cnt+1

# No. of data points:                            
print('No. of valid data points: ',cnt,flush=True)
print('Max no. of data points in a gridsquare: ',np.amax(gcnt),flush=True)

# Track reasons for data loss:
print('(1) Too few points: ',num[0],flush=True)
print('(2) Low cloud height range: ',num[1],flush=True)
print('(3) Low cloud height std dev: ',num[2],flush=True)
print('(4) Large error: ',num[3],flush=True)
print('(5) Significantly less then zero: ',num[4],flush=True)
print('(6) Outlier (NO2 > 200 pptv): ',num[5],flush=True)
print('(7) Non-uniform stratosphere: ',num[6],flush=True)

# Get average:
gno2vmr=np.divide(gno2vmr, gerr, where=gcnt!=0)
gcldfr=np.divide(gcldfr, gcnt, where=gcnt!=0)
trueno2=np.divide(trueno2, gerr, where=gcnt!=0)
gcld=np.divide(gcld, gcnt, where=gcnt!=0)
gerr=np.divide(gerr, gcnt, where=gcnt!=0)
gcutno2=np.divide(gcutno2,truewgt,where=gascnt!=0)
trueno2[gcnt==0]=np.nan
gerr[gcnt==0]=np.nan
gcutno2[gascnt==0]=np.nan
gno2vmr[gcnt==0]=np.nan
gcldfr[gcnt==0]=np.nan

# Plot the data:
m=Basemap(resolution='l',projection='merc',\
          lat_0=0,lon_0=0,llcrnrlon=minlon,\
          llcrnrlat=minlat,urcrnrlon=maxlon,urcrnrlat=maxlat)
xi,yi=m(X,Y)
plt.subplot(2, 3, 1)
cs=m.pcolor(xi,yi,np.squeeze(gno2vmr), vmin=0, vmax=80, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Cloud-sliced NO2 VMRs')

plt.subplot(2, 3, 2)
cs=m.pcolor(xi,yi,np.squeeze(gcldfr), vmin=0.5, vmax=2.0, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('Cloud fraction')

plt.subplot(2, 3, 3)
cs=m.pcolor(xi,yi,np.squeeze(gcutno2), vmin=0, vmax=80, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('True NO2 VMRs under all-sky conditions')

plt.subplot(2, 3, 4)
plt.plot(trueno2,gno2vmr, 'o', color='black',markersize=6)
r=stats.pearsonr(trueno2[~np.isnan(gno2vmr)],\
                 gno2vmr[~np.isnan(gno2vmr)])
#print('Correlation = ', r[0])
result=rma(trueno2[~np.isnan(gno2vmr)],gno2vmr[~np.isnan(gno2vmr)],\
           len(trueno2[~np.isnan(gno2vmr)]),1000)
print(result,flush=True)
xvals=np.arange(0,100,5)
yvals=result[1]+xvals*result[0]
plt.plot(xvals, yvals, '-')
plt.xlim(-4,80)
plt.ylim(-4,80)
plt.xlabel('True NO2 (cloudy)')
plt.ylabel('Cloud-sliced NO2')
print('===== True (cloudy) vs cloud-sliced UT NO2 ====')
print('R = ',r[0],flush=True)
print('Slope = ',result[0])
print('Slope Err = ',result[2],flush=True)
print('Intercept = ',result[1],flush=True)
print('Intercept Err = ',result[3],flush=True)
add2plt=("y = {a:6.2f}x + {b:6.3f}".\
         format(a=result[0],b=result[1]))
plt.text(2,75,add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R = {a:6.2f}".format(a=r[0]))
plt.text(2,65, add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)

plt.subplot(2, 3, 5)
plt.plot(gcutno2,gno2vmr, 'o', color='black',markersize=6)
r=stats.pearsonr(gcutno2[(~np.isnan(gno2vmr))&(~np.isnan(gcutno2))],\
                 gno2vmr[(~np.isnan(gno2vmr))&(~np.isnan(gcutno2))])
result=rma(gcutno2[(~np.isnan(gno2vmr))&(~np.isnan(gcutno2))],\
           gno2vmr[(~np.isnan(gno2vmr))&(~np.isnan(gcutno2))],\
                   len(gcutno2[(~np.isnan(gno2vmr))&(~np.isnan(gcutno2))]),\
                   1000)
xvals=np.arange(0,100,5)
yvals=result[1]+xvals*result[0]
plt.plot(xvals, yvals, '-')
plt.xlim(-4,80)
plt.ylim(-4,80)
plt.xlabel('True NO2 (all-sky)')
plt.ylabel('Cloud-sliced NO2')
add2plt=("y = {a:6.2f}x + {b:6.3f}".\
         format(a=result[0],b=result[1]))
plt.text(2,75,add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R = {a:6.2f}".format(a=r[0]))
plt.text(2,65, add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)

plt.subplot(2, 3, 6)
plt.plot(gcutno2,trueno2, 'o', color='black',markersize=6)
r=stats.pearsonr(gcutno2[(~np.isnan(trueno2))&(~np.isnan(gcutno2))],\
                 trueno2[(~np.isnan(trueno2))&(~np.isnan(gcutno2))])
result=rma(gcutno2[(~np.isnan(trueno2))&(~np.isnan(gcutno2))],\
           trueno2[(~np.isnan(trueno2))&(~np.isnan(gcutno2))],\
                   len(gcutno2[(~np.isnan(trueno2))&(~np.isnan(gcutno2))]),\
                   1000)
xvals=np.arange(0,100,5)
yvals=result[1]+xvals*result[0]
plt.plot(xvals, yvals, '-')
plt.xlim(-4,80)
plt.ylim(-4,80)
plt.xlabel('True NO2 (all-sky)')
plt.ylabel('True NO2 (cloudy)')
add2plt=("y = {a:6.2f}x + {b:6.3f}".\
         format(a=result[0],b=result[1]))
plt.text(2,75,add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)
add2plt=("R = {a:6.2f}".format(a=r[0]))
plt.text(2,65, add2plt, fontsize=8,\
         ha='left', va='center')#, transform=ax.transAxes)
plt.show()

# Save the data to NetCDF:
ncout=Dataset('./Data/gc-v12-1-0-ut-no2-'+Reg.lower()+'-jja-'+yrrange+'-'+\
              StrRes+'-'+prange+'-v3.nc',mode='w',format='NETCDF4')

# Create data dimensions:
ncout.createDimension('lat', ydim) 
ncout.createDimension('lon', xdim) 

# create lon axis:
lon = ncout.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
lon[:] = out_lon

# create lat axis:
lat = ncout.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latgitude'
lat[:] = out_lat

# create data axes:
# (1) Cloud-sliced NO2:
csutno2 = ncout.createVariable('csutno2', np.float32, ('lon','lat'))
csutno2.units = 'pptv'
csutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained using cloud-slicing'
csutno2[:] = gno2vmr

# (2) Cloud-sliced NO2 error:
utno2err = ncout.createVariable('utno2err', np.float32, ('lon','lat'))
utno2err.units = 'pptv'
utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
utno2err[:] = gerr

# (3) Number of observations in each gridsquare:
nobs = ncout.createVariable('nobs', np.float32, ('lon','lat'))
nobs.units = 'unitless'
nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
nobs[:] = gcnt

# (4) Mean cloud pressure for season between 450-180 hPa:
utcld = ncout.createVariable('utcld', np.float32, ('lon','lat'))
utcld.units = 'hPa'
utcld.long_name = 'Mean cloud pressure between 450 and 180 hPa'
utcld[:] = gcld

# (5) Mean NO2 mixing ratio at 450-180 hPa for scenes with clouds:
cldutno2 = ncout.createVariable('cldutno2', np.float32, ('lon','lat'))
cldutno2.units = 'pptv'
cldutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained if clouds are present'
cldutno2[:] = trueno2

# (6) Mean NO2 mixing ratio at 450-180 hPa under all conditions (all-sky):
askutno2 = ncout.createVariable('askutno2', np.float32, ('lon','lat'))
askutno2.units = 'pptv'
askutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained under all conditions (all-sky)'
askutno2[:] = gcutno2

# (7) Cloud fraction:
utcldfrc = ncout.createVariable('utcldfrc', np.float32, ('lon','lat'))
utcldfrc.units = 'unitless'
utcldfrc.long_name = 'GEOS-FP cloud fraction obtained as sum of 3D cloud fractions across range of interest (180-450 hPa)'
utcldfrc[:] = gcldfr

# Close the file:
ncout.close()

# Close the log file:
log.close()
