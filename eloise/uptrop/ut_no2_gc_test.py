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


# Validation against mode
# E runs model over EU, CH, NA. EU is quickest.
# Spits out a netCDF with NO2 from retreval, sat clear, sat cloudy
# Runs about 20 minutes per month


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
import leastsq  # E has replaced this with error-bar having version(so better)

from gcpy.gcpy.constants import AVOGADRO
from gcpy.gcpy.constants import G
from gcpy.gcpy.constants import MW_AIR

# Turn off warnings:
np.warnings.filterwarnings('ignore')

# Decide on region:
REGION= 'EU' #'NA', 'EU', or 'CH'

# Define information for grid:
STR_RES= '4x5'

# Define pressure range:
# These are fixed, so can be constants, as opposed to inputs.
P_MIN=180     
P_MAX=450     


# Define years of interest:
stryr=['2016','2017']
if len(stryr)==1: yrrange=stryr[0]
if len(stryr)==2: yrrange='2016-2017'


class ProcessedData:

    def __init__(self, X):
        # Define output arrays:
        self.gno2vcd = np.zeros(X.shape)
        self.gno2vmr = np.zeros(X.shape)
        self.gcldfr = np.zeros(X.shape)
        self.gcld = np.zeros(X.shape)
        self.gerr = np.zeros(X.shape)
        self.trueo3 = np.zeros(X.shape)  # ozone mixing ratio coincident with cloud-sliced ut self
        self.trueno2 = np.zeros(X.shape)  # "true" cloudy UT NO2
        self.gcutno2 = np.zeros(X.shape)  # "true" all-sky UT NO2
        self.truewgt = np.zeros(X.shape)  # "true" all-sky UT NO2 Gaussian weights
        self.gascnt = np.zeros(X.shape)  # Count all-sky
        self.gcnt = np.zeros(X.shape)

        self.loss_count = np.zeros(7)  # TODO: Refactor to dict. These are counts of loss due to specific errors.

        # Initialize:
        self.cnt = 0
        # numcnt=0
        self.maxcnt = 0

        # Define string to represent the layer range:
        self.prange = str(P_MIN) + '-' + str(P_MAX)

        # Define factor to convert slope of NO2 mixing ratio versus pressure
        # to VMR:
        self.den2mr = np.divide((np.multiply(G, MW_AIR)), AVOGADRO)

    def regrid_and_process(self, x, y, no2):
        # Find nearest gridsquare in output grid:
        p, q = np.argmin(abs(out_lon - no2.tlon[x])), \
               np.argmin(abs(out_lat - no2.tlat[y]))

        self.truewgt[p, q] += np.sum(no2.twgt)

        # Sum up "true" all-sky UT NO2:
        self.gcutno2[p, q] += np.sum(no2.tgcno2[no2.askind, y, x] * no2.twgt * 1e3)
        self.gascnt[p, q] += 1.0

        no2.prepare_no2(no2, x, y)

        # Add relevant data to array:
        gno2[p, q, int(cntloop[p, q])] = no2.no22d  # stratcol+tropcol   #no22d
        stratno2[p, q, int(cntloop[p, q])] = no2.stratcol
        gcldp[p, q, int(cntloop[p, q])] = tcldhgt[y, x]
        go3[p, q, int(cntloop[p, q])] = np.mean(tgco3[no2.lmin:no2.lmax + 1, y, x])
        gtrueno2[p, q, int(cntloop[p, q])] = np.mean(tgcno2[no2.lmin:no2.lmax + 1, y, x])
        allcldfr[p, q, int(cntloop[p, q])] = np.sum(tcldfr[no2.lmin:no2.lmax + 1, y, x])
        # Increment indices:
        cntloop[p, q] += 1

    def print_data_loss_report(self):
        # No. of data points:
        # print('No. of valid data points: ',cnt,flush=True)
        print('Max no. of data points in a gridsquare: ', np.amax(self.gcnt), flush=True)
        # Track reasons for data loss:
        print('(1) Too few points: ', self.loss_count[0], flush=True)
        print('(2) Low cloud height range: ', self.loss_count[1], flush=True)
        print('(3) Low cloud height std dev: ', self.loss_count[2], flush=True)
        print('(4) Large error: ', self.loss_count[3], flush=True)
        print('(5) Significantly less then zero: ', self.loss_count[4], flush=True)
        print('(6) Outlier (NO2 > 200 pptv): ', self.loss_count[5], flush=True)
        print('(7) Non-uniform stratosphere: ', self.loss_count[6], flush=True)
        print('(8) Successful retrievals: ', self.cnt, flush=True)
        print('(9) Total possible points: ', (np.sum(self.loss_count) + self.cnt), flush=True)

    def process_file(self,file_path):
        global tlon, tlat, tgcno2, tcldfr, tcldhgt, tadn, tbxhgt, tpedge, tpause, tgco3, gno2, go3, allcldfr, stratno2, gcldp, gtrueno2, cntloop, maxcnt, cnt
        no2_data = N02_Data(file_path)
        # Get column values:
        for y in range(len(no2_data.tlat)):
            for x in range(len(no2_data.tlon)):
                self.regrid_and_process()
        for i in range(xdim):
            for j in range(ydim):
                if no2_data.gno2[i, j, 0] != 0:
                    self.process_grid_square(no2_data, i, j)

        self.get_average()

    def get_average(self):
        # Get average:
        gno2vmr = np.divide(self.gno2vmr, self.gerr, where=self.gcnt != 0)
        gcldfr = np.divide(self.gcldfr, self.gcnt, where=self.gcnt != 0)
        trueno2 = np.divide(self.trueno2, self.gerr, where=self.gcnt != 0)
        trueo3 = np.divide(self.trueo3, self.gerr, where=self.gcnt != 0)
        gcld = np.divide(self.gcld, self.gcnt, where=self.gcnt != 0)
        gerr = np.divide(self.gerr, self.gcnt, where=self.gcnt != 0)
        gcutno2 = np.divide(self.gcutno2, self.truewgt, where=self.gascnt != 0)
        trueno2[self.gcnt == 0] = np.nan
        trueo3[self.gcnt == 0] = np.nan
        gerr[self.gcnt == 0] = np.nan
        gcutno2[self.gascnt == 0] = np.nan
        gno2vmr[self.gcnt == 0] = np.nan
        gcldfr[self.gcnt == 0] = np.nan

    def process_grid_square(self, no2_data, i, j):
        # Define vectors of relevant data:
        tcolno2 = no2_data.gno2[i, j, :]
        tstratno2 = no2_data.stratno2[i, j, :]
        tfrc = no2_data.allcldfr[i, j, :]
        tcld = no2_data.gcldp[i, j, :]
        tmrno2 = no2_data.gtrueno2[i, j, :]
        to3 = no2_data.go3[i, j, :]
        # Skip if fewer than 10 points:
        if len(tcolno2) < 10:
            print("Grid square {}, {} failed with condition 0".format(i, j))
            self.loss_count[0] += 1
            return

            # Trim to remove trailing zeros:
        tcolno2 = np.trim_zeros(tcolno2)  # [np.nonzero(tcolno2)]
        tstratno2 = np.trim_zeros(tstratno2)  # [np.nonzero(tstratno2)]
        tfrc = np.trim_zeros(tfrc)  # [np.nonzero(tfrc)]
        tcld = np.trim_zeros(tcld)  # [np.nonzero(tcld)]
        tmrno2 = np.trim_zeros(tmrno2)  # [np.nonzero(tmrno2)]
        to3 = np.trim_zeros(to3)  # [np.nonzero(tmrno2)]
        # Remove non-uniform stratosphere:
        if (np.std(tstratno2) / np.mean(tstratno2)) > 0.02:
            self.loss_count[6] += 1
            pass  # continue
        # Check if there are any repeat values:
        u, ind = np.unique(tcld, return_inverse=True)
        if (len(ind) != len(tcld)):
            print('Repeat cloud values', flush=True)
            sys.exit()
        # Get number of points:
        npnts = len(tcld)
        if npnts > self.maxcnt:
            maxcnt = npnts
            print(maxcnt, flush=True)

        # TODO: Check if the treshold ratios will need to be tweaked
        # No, these are fixed. The other numbers were just experiments. I've deleted those now for clarity.
        # Use cloud_slice_ut_no2 function to get cloud-sliced
        # UT NO2 mixing ratios.
        # Treat data differently depending on whether there are no
        # 20-100 points:
        if ((npnts >= 20) & (npnts < 100)):
            csval = cldslice(tcolno2, tcld)
            if np.isnan(csval[0]) or np.isnan(csval[1]):
                self.loss_count[csval[2] - 1] += 1
                pass  # continue

            # Get mean cloud top pressure:
            pmean = csval[3]

            # Calculate Gaussian weight:
            gwgt = np.exp(np.multiply((-1.0), np.divide((np.square(np.subtract(pmean, 315.))),
                                                        (np.multiply(2, np.square(135.))))))

            # Apply Gaussian weights to cloud-sliced UT NO2:
            self.gno2vmr[i, j] += np.multiply(csval[0], gwgt)
            self.gcldfr[i, j] += np.mean(tfrc)
            self.gerr[i, j] += gwgt

            # "True" cloudy UT NO2 (converted from ppbv to pptv):
            self.trueno2[i, j] += np.multiply(np.mean(tmrno2 * 1e3), gwgt)

            # Ozone coincident with cloud-sliced NO2:
            self.trueo3[i, j] += np.multiply(np.mean(to3), gwgt)

            self.gcnt[i, j] += 1
            self.cnt += 1

        elif (npnts >= 100):
            #    # Code error check:
            #    if len(tcld)!=len(tmrno2):
            #        print('arrays not =',flush=True)
            #        print(len(tcld),len(tmrno2),flush=True)

            #    # Determine the number of iterations:
            niter = round(npnts / 40)
            nloop = list(range(niter))

            # TODO: Get what this is from E. Looks like some kind of cumulative function?
            # This was my attempt at subsampling the data: if there are more than 100 satellite pixels in a 1 deg x 1 deg 
            # array I split these up be dividing the number of observations by 40 to obtain niter (the number of iterations)
            # and looping over these iterations.  
            #    # Loop over iterations:
            for w in nloop:
                if w == 0:

                    # Cloud-slicing:
                    csval = cldslice(tcolno2[::niter], tcld[::niter])
                    # Get mean cloud top pressure:
                    pmean = csval[3]
                    # Calculate Gaussian weight:
                    # This could appear as a separate function called "Gaussian weight" that requires as input the average
                    # cloud pressure pmean, the center cloud pressure (315) and the cloud standard deviation (135).
                    gwgt = np.exp(np.multiply((-1.0), np.divide((np.square(np.subtract(pmean, 315.))),
                                                                (np.multiply(2, np.square(135.))))))
                else:

                    # Cloud-slicing:
                    csval = cldslice(tcolno2[w::niter], tcld[w::niter])
                    # Get mean cloud top pressure:
                    pmean = csval[3]
                    # Calculate Gaussian weight:
                    gwgt = np.exp(np.multiply((-1.0), np.divide((np.square(np.subtract(pmean, 315.))),
                                                                (np.multiply(2, np.square(135.))))))

            # Skip if approach didn't work (i.e. cloud-sliced
            # UT NO2 is NaN):
            if np.isnan(no2_data.csval[0]) or np.isnan(no2_data.csval[1]):
                self.loss_count[no2_data.csval[2] - 1] += 1
                pass  # continue

            # Gaussian-weighted mean:
            self.gno2vmr[i, j] += np.multiply(csval[0], gwgt)
            self.gerr[i, j] += gwgt

            # "True" cloudy UT NO2 (converted from ppbv to pptv):
            if (w == 0):
                self.trueno2[i, j] += np.multiply(np.mean(tmrno2[::niter] * 1e3), gwgt)
                self.gcldfr[i, j] += np.mean(tfrc[::niter])
            else:
                self.trueno2[i, j] += np.multiply(np.mean(tmrno2[w::niter] * 1e3), gwgt)
                self.gcldfr[i, j] += np.mean(tfrc[w::niter])  #

            self.gcnt[i, j] += 1
            self.cnt += 1

    def plot_data(self):
        # TODO: Make variables members
        # Plot the data:
        m = Basemap(resolution='l', projection='merc', \
                    lat_0=0, lon_0=0, llcrnrlon=minlon, \
                    llcrnrlat=minlat, urcrnrlon=maxlon, urcrnrlat=maxlat)
        xi, yi = m(X, Y)
        plt.subplot(2, 3, 1)
        cs = m.pcolor(xi, yi, np.squeeze(gno2vmr), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Cloud-sliced NO2 VMRs')
        plt.subplot(2, 3, 2)
        cs = m.pcolor(xi, yi, np.squeeze(trueo3), vmin=40.0, vmax=120, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('O3 VMR [ppbv]')
        plt.subplot(2, 3, 3)
        cs = m.pcolor(xi, yi, np.squeeze(gcutno2), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('True NO2 VMRs under all-sky conditions')
        plt.subplot(2, 3, 4)
        plt.plot(trueno2, gno2vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(trueno2[~np.isnan(gno2vmr)], \
                           gno2vmr[~np.isnan(gno2vmr)])
        # print('Correlation = ', r[0])
        result = rma(trueno2[~np.isnan(gno2vmr)], gno2vmr[~np.isnan(gno2vmr)], \
                     len(trueno2[~np.isnan(gno2vmr)]), 1000)
        print(result, flush=True)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (cloudy)')
        plt.ylabel('Cloud-sliced NO2')
        print('===== True (cloudy) vs cloud-sliced UT NO2 ====')
        print('R = ', r[0], flush=True)
        print('Slope = ', result[0])
        print('Slope Err = ', result[2], flush=True)
        print('Intercept = ', result[1], flush=True)
        print('Intercept Err = ', result[3], flush=True)
        add2plt = ("y = {a:6.2f}x + {b:6.3f}". \
                   format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.subplot(2, 3, 5)
        plt.plot(gcutno2, gno2vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(gcutno2[(~np.isnan(gno2vmr)) & (~np.isnan(gcutno2))], \
                           gno2vmr[(~np.isnan(gno2vmr)) & (~np.isnan(gcutno2))])
        result = rma(gcutno2[(~np.isnan(gno2vmr)) & (~np.isnan(gcutno2))], \
                     gno2vmr[(~np.isnan(gno2vmr)) & (~np.isnan(gcutno2))], \
                     len(gcutno2[(~np.isnan(gno2vmr)) & (~np.isnan(gcutno2))]), \
                     1000)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (all-sky)')
        plt.ylabel('Cloud-sliced NO2')
        add2plt = ("y = {a:6.2f}x + {b:6.3f}". \
                   format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.subplot(2, 3, 6)
        plt.plot(gcutno2, trueno2, 'o', color='black', markersize=6)
        r = stats.pearsonr(gcutno2[(~np.isnan(trueno2)) & (~np.isnan(gcutno2))], \
                           trueno2[(~np.isnan(trueno2)) & (~np.isnan(gcutno2))])
        result = rma(gcutno2[(~np.isnan(trueno2)) & (~np.isnan(gcutno2))], \
                     trueno2[(~np.isnan(trueno2)) & (~np.isnan(gcutno2))], \
                     len(gcutno2[(~np.isnan(trueno2)) & (~np.isnan(gcutno2))]), \
                     1000)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (all-sky)')
        plt.ylabel('True NO2 (cloudy)')
        add2plt = ("y = {a:6.2f}x + {b:6.3f}". \
                   format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, \
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.show()

    def save_to_netcdf(self, out_path):
        # TODO: Make variables members
        # Save the data to NetCDF:
        ncout = Dataset(out_path, mode='w', format='NETCDF4')
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
        csutno2 = ncout.createVariable('csutno2', np.float32, ('lon', 'lat'))
        csutno2.units = 'pptv'
        csutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained using cloud-slicing'
        csutno2[:] = gno2vmr
        # (2) Cloud-sliced NO2 error:
        utno2err = ncout.createVariable('utno2err', np.float32, ('lon', 'lat'))
        utno2err.units = 'pptv'
        utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2err[:] = gerr
        # (3) Number of observations in each gridsquare:
        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = gcnt
        # (4) Mean cloud pressure for season between 450-180 hPa:
        utcld = ncout.createVariable('utcld', np.float32, ('lon', 'lat'))
        utcld.units = 'hPa'
        utcld.long_name = 'Mean cloud pressure between 450 and 180 hPa'
        utcld[:] = gcld
        # (5) Mean NO2 mixing ratio at 450-180 hPa for scenes with clouds:
        cldutno2 = ncout.createVariable('cldutno2', np.float32, ('lon', 'lat'))
        cldutno2.units = 'pptv'
        cldutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained if clouds are present'
        cldutno2[:] = trueno2
        # (6) Mean NO2 mixing ratio at 450-180 hPa under all conditions (all-sky):
        askutno2 = ncout.createVariable('askutno2', np.float32, ('lon', 'lat'))
        askutno2.units = 'pptv'
        askutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained under all conditions (all-sky)'
        askutno2[:] = gcutno2
        # (7) Cloud fraction:
        utcldfrc = ncout.createVariable('utcldfrc', np.float32, ('lon', 'lat'))
        utcldfrc.units = 'unitless'
        utcldfrc.long_name = 'GEOS-FP cloud fraction obtained as sum of 3D cloud fractions across range of interest (180-450 hPa)'
        utcldfrc[:] = gcldfr
        # (8) O3 sampled coincident with cloud-slicing retrieval:
        uto3 = ncout.createVariable('uto3', np.float32, ('lon', 'lat'))
        uto3.units = 'ppbv'
        uto3.long_name = 'GEOS-Chem ozone obtained coincident with cloud-sliced NO2'
        uto3[:] = trueo3
        # Close the file:
        ncout.close()


class N02_Data():

    def __init__(self, file_path):
        print(file_path, flush=True)
        # Read dataset:
        fh = Dataset(file_path, mode='r')
        # Extract data of interest:
        # (Add tropopause height to this in the future)
        tlon, tlat, tgcno2, tcldfr, tcldhgt, tadn, tbxhgt, tpedge, tpause, tgco3 = \
            fh.variables['LON'], fh.variables['LAT'], \
            fh.variables['IJ-AVG-S__NO2'], fh.variables['TIME-SER__CF'], \
            fh.variables['TIME-SER__CThgt'], fh.variables['TIME-SER__AIRDEN'], \
            fh.variables['BXHGHT-S__BXHEIGHT'], fh.variables['PEDGE-S__PSURF'], \
            fh.variables['TR-PAUSE__TP-PRESS'], fh.variables['IJ-AVG-S__O3']
        self.tlon = tlon[:]
        self.tlat = tlat[:]
        self.tgcno2 = tgcno2[:]
        self.tcldfr = tcldfr[:]
        self.tcldhgt = tcldhgt[0, :, :]
        self.tadn = tadn[:]  # in molec/cm3
        self.tbxhgt = tbxhgt[:]
        self.tpedge = tpedge[:]
        self.tpause = tpause[0, :, :]
        self.tgco3 = tgco3[:]
        # Convert box height from m to cm:
        self.tbxhgt = tbxhgt * 1e2
        # Define output data for this day:
        self.gno2 = np.zeros((xdim, ydim, nval))
        self.go3 = np.zeros((xdim, ydim, nval))
        self.allcldfr = np.zeros((xdim, ydim, nval))
        self.stratno2 = np.zeros((xdim, ydim, nval))
        self.gcldp = np.zeros((xdim, ydim, nval))
        self.gtrueno2 = np.zeros((xdim, ydim, nval))
        self.cntloop = np.zeros((xdim, ydim))

    def filter_pixel(self, x, y):
        # Calculate corresponding mid-pressure values:
        self.tpmid = np.zeros(len(self.tpedge[:, y, x]))
        # Get mid-pressure values, except for highest layer:
        for k in range(len(self.tpedge[:, y, x]) - 1):
            self.tpmid[k] = np.multiply(0.5, (self.tpedge[k, y, x] + self.tpedge[k + 1, y, x]))
        # Estimate mid-pressure for highest value (doesn't need to
        # be accurate, as surpassing the range of interset):
        self.tpmid[46] = np.multiply(0.5, (self.tpedge[46, y, x] + (self.tpedge[46, y, x] - 0.1)))  # Why 46?
        # Get model layer of tropopause:
        self.tppind = np.argmin(abs(self.tpmid - self.tpause[y, x]))
        # Get indices that fall between 450 and 180 hPa for estimating
        # "true' all-sky UT NO2 and partial columns:
        lind = np.where((self.tpmid >= P_MIN) & (self.tpmid <= P_MAX))[0]
        # Get UT NO2 under "true" all-sky conditions:
        # Make sure this is below the tropopause:
        # If below tropopause, use full extent (180-450 hPa):
        if (lind[len(lind) - 1] <= self.tppind):
            askind = lind
        # If above tropopause, trim to tropopause-450 hPa:
        if (lind[len(lind) - 1] > self.tppind):
            askind = lind[np.where(lind <= self.tppind)[0]]
        # If tropopause below 450 hPa, skip entirely:
        if (self.tpause[y, x] > P_MAX):
            return # continue
        # Get Guassian weights that allocated higher weights to points
        # closest to the pressure centre (315 hPa):
        # Equation is:
        #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
        #         135 hPa is the standard deviation.
        self.twgt = np.exp(np.multiply((-1.0), np.divide((np.square(np.subtract(self.tpmid[askind], 315.))),
                                                    (np.multiply(2, np.square(135.))))))

    def prepare_no2(self, x, y):

        # Find where cloud fraction in UT exceeds 0.7 after calculating
        # true all-sky NO2:
        # (Keep for testing effect of thick clouds on cloud-sliced UT NO2):
        # if (np.sum(tcldfr[lind,y,x])<=0.7): continue
        # Get model level of cloud top height closest to lowest
        # pressure-altitude of interest (P_MIN):
        self.lcld = np.argmin(abs(self.tcldhgt[y, x] - self.tpmid))
        # Skip if cloud top height ouside pressure range of interest:
        lmin, lmax = np.amin(self.lind), np.amax(self.lind)
        if (self.lcld < lmin) or (self.lcld > lmax):
            pass  # continue
        # This error check is probably redundant, but included in case:
        if self.lcld == 0:
            print("No cloud detected!!!", flush=True)
            sys.exit()
        # Get partial NO2 column in molec/m2 from cloud top height
        # to highest model level (output up to level 47):
        # print(tgcno2[lmin:tppind,y,x])
        # print(tgcno2[lmin:tppind,y,x]*1.5)
        self.no22d = np.sum(tgcno2[lmin:, y, x] * 1e-5 * tadn[lmin:, y, x] * \
                       tbxhgt[lmin:, y, x])
        self.tropcol = np.sum(tgcno2[lmin:self.tppind, y, x] * 1e-5 * \
                         tadn[lmin:self.tppind, y, x] * tbxhgt[lmin:self.tppind, y, x])
        # Get stratospheric column from 180 hPa aloft:
        # Previous approach (remove when model simulations done):
        # tppind=np.where(tpmid<180.)[0]
        self.stratcol = np.sum(tgcno2[self.tppind:, y, x] * 1e-5 * tadn[self.tppind:, y, x] * \
                          tbxhgt[tppind:, y, x])


def get_file_list(gcdir):
    global files, i
    # TODO: Replace with os.path
    files = glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + '.' + stryr[0] + '06*')
    mon = ["07", "08"]
    for i in mon:
        for filename in glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + \
                                  '.' + stryr[0] + i + '*'):
            files.append(filename)
    # 2017:
    mon = ["06", "07", "08"]
    for i in mon:
        for filename in glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + \
                                  '.' + stryr[1] + i + '*'):
            files.append(filename)
    return files


if __name__ == "__main__":

    # Name of log file to output code prorgess:
    # TODO: Maybe replace this with Python logging? See rest of code.
    log = open("log_" + REGION + "_" + STR_RES + "_" + yrrange + "_v4", "w")
    sys.stdout = log

    # Define target grid:
    if REGION == 'NA':
        minlat = 6.
        maxlat = 58.
        minlon = -135.
        maxlon = -60.
        dirreg = 'NA'
        factor = 40  # used to determine size of gno2
    if REGION == 'EU':
        minlat = 30.
        maxlat = 62.
        minlon = -20.
        maxlon = 40.
        dirreg = 'eu_naei'
        factor = 30
    if REGION == 'CH':
        minlat = 10.
        maxlat = 54.
        minlon = 65.
        maxlon = 135.
        dirreg = 'ch'
        factor = 100

    # NOTE FOR LATER: This snippet turns up in fresco_cld_err; a candidate for the library.
        # Define grid information:
    if STR_RES == '8x10': dellat, dellon = 8, 10
    if STR_RES == '4x5': dellat, dellon = 4, 5
    if STR_RES == '2x25': dellat, dellon = 2, 2.5
    if STR_RES == '1x1': dellat, dellon = 1, 1
    out_lon = np.arange(minlon, maxlon + dellon, dellon) 
    out_lat = np.arange(minlat, maxlat + dellat, dellat) 
    # Convert output lats and long to 2D:
    X, Y = np.meshgrid(out_lon, out_lat, indexing='ij')
    # Dimensions of output data:
    xdim = len(out_lon)
    ydim = len(out_lat)
    nval = int(factor * dellat * dellon)  # Scale array size with resolution



    # Get files (test June 2016 for now)
    # 2016:
    gcdir = '/data/uptrop/Projects/DEFRA-NH3/GC/geosfp_' + dirreg + '_iccw/'
    files = get_file_list()

    # Order the files:
    files = sorted(files)

    # Check:
    print('Number of files:', len(files), flush=True)

    # Loop over files:
    for file_path in files:
        process_file(file_path)

    print_data_loss_report()
    plot_data()
    out_path = './Data/gc-v12-1-0-ut-self-' + REGION.lower() + '-jja-' + yrrange + '-' + STR_RES + '-' + prange + '-v4.nc'
    save_to_netcdf(out_path)

    # Close the log file:
    log.close()
