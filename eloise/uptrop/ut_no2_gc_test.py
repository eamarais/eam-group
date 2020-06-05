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
import argparse

from gcpy.gcpy.constants import AVOGADRO
from gcpy.gcpy.constants import G
from gcpy.gcpy.constants import MW_AIR

# Turn off warnings:
np.warnings.filterwarnings('ignore')



# Define pressure range:
# These are fixed, so can be constants, as opposed to inputs.
P_MIN=180     
P_MAX=450     


# Define years of interest:
YEARS_TO_PROCESS=['2016', '2017']

# TODO: Go through code and replace pass statements with exceptions. Add try-catch.
class ProcessingException(Exception):
    pass

class ProcessedData:

    # ----Initialisation methods----
    def __init__(self, region, str_res):

        self.define_grid(region, str_res)
        grid_shape = self.X.shape

        # Define output arrays:
        self.g_no2_vcd = np.zeros(grid_shape)
        self.g_no2_vmr = np.zeros(grid_shape)
        self.g_cld_fr = np.zeros(grid_shape)
        self.g_cld = np.zeros(grid_shape)
        self.g_err = np.zeros(grid_shape)
        self.true_o3 = np.zeros(grid_shape)  # ozone mixing ratio coincident with cloud-sliced ut self
        self.true_no2 = np.zeros(grid_shape)  # "true" cloudy UT NO2
        self.g_cut_no2 = np.zeros(grid_shape)  # "true" all-sky UT NO2
        self.true_wgt = np.zeros(grid_shape)  # "true" all-sky UT NO2 Gaussian weights
        self.g_as_cnt = np.zeros(grid_shape)  # Count all-sky
        self.g_cnt = np.zeros(grid_shape)

        # Define output data for this day:
        out_shape = (self.xdim, self.ydim, self.nval)
        self.g_no2 = np.zeros(out_shape)
        self.g_o3 = np.zeros(out_shape)
        self.g_cld_fr = np.zeros(out_shape)
        self.strat_no2 = np.zeros(out_shape)
        self.g_cld_p = np.zeros(out_shape)
        self.g_true_no2 = np.zeros(out_shape)
        self.cnt_loop = np.zeros(out_shape)

        #NOTE: Lots of these didn't seem to appear in the code
        self.loss_count = {
            "too_few_points": 0,
            "low_cloud_height_range": 0,
            "low_cloud_height_std": 0,
            "large_error": 0,
            "much_less_than_zero": 0,
            "no2_outlier": 0,
            "non_uni_strat": 0,
        }

        # Initialize:
        self.cnt = 0
        self.maxcnt = 0

        # Define string to represent the layer range:
        self.prange = str(P_MIN) + '-' + str(P_MAX)

        # Define factor to convert slope of NO2 mixing ratio versus pressure
        # to VMR:
        self.den2mr = np.divide((np.multiply(G, MW_AIR)), AVOGADRO)

    def define_grid(self, region, str_res):
        # Define target grid:
        if region == 'NA':
            self.minlat = 6.
            self.maxlat = 58.
            self.minlon = -135.
            self.maxlon = -60.
            self.dirreg = 'NA'
            self.factor = 40  # used to determine size of g_no2
        if region == 'EU':
            self.minlat = 30.
            self.maxlat = 62.
            self.minlon = -20.
            self.maxlon = 40.
            self.dirreg = 'eu_naei'
            self.factor = 30
        if region == 'CH':
            self.minlat = 10.
            self.maxlat = 54.
            self.minlon = 65.
            self.maxlon = 135.
            self.dirreg = 'ch'
            self.factor = 100
        else:
            print("Invalid region; valid regions are 'NA','EU','CH'.")
            sys.exit(1)
        # NOTE FOR LATER: This snippet turns up in fresco_cld_err; a candidate for the library.
        # Define grid information:
        if str_res == '8x10':
            self.dellat, self.dellon = 8, 10
        if str_res == '4x5':
            self.dellat, self.dellon = 4, 5
        if str_res == '2x25':
            self.dellat, self.dellon = 2, 2.5
        if str_res == '1x1':
            self.dellat, self.dellon = 1, 1
        else:
            print("Invalid resolution: valid resolutions are 8x10, 4x5, 2x25 (two pt 5) and 1x1")
            sys.exit(1)
        self.out_lon = np.arange(self.minlon, self.maxlon + self.dellon, self.dellon)
        self.out_lat = np.arange(self.minlat, self.maxlat + self.dellat, self.dellat)
        # Convert output lats and long to 2D:
        self.X, self.Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')
        # Dimensions of output data:
        self.xdim = len(self.out_lon)
        self.ydim = len(self.out_lat)
        self.nval = int(self.factor * self.dellat * self.dellon)  # Scale array size with resolution

    # ----Processing methods----
    def process_geochem_day(self, file_path):
        """
        I _think_ this performs some preprocessing on a partiuclar set of no2 pixels, then saves the
        result in a grid square, then performs an update on that grid square.
        """
        this_geoschem_day = GeosChemDay(file_path)
        # Get column values:
        for y in range(len(this_geoschem_day.t_lat)):
            for x in range(len(this_geoschem_day.t_lon)):
                this_geoschem_day.prepare_no2_pixel(x, y)
                self.regrid_and_process(x, y, this_geoschem_day)
        for i in range(self.xdim):
            for j in range(self.ydim):
                if self.g_no2[i, j, 0] != 0:
                    self.process_grid_square(i, j)

    def regrid_and_process(self, x, y, no2):
        # Find nearest gridsquare in output grid:
        lon = np.argmin(abs(self.out_lon - no2.t_lon[x]))
        lat = np.argmin(abs(self.out_lat - no2.t_lat[y]))

        self.true_wgt[lon, lat] += np.sum(no2.twgt)

        # Sum up "true" all-sky UT NO2:
        self.g_cut_no2[lon, lat] += np.sum(no2.t_gc_no2[no2.askind, y, x] * no2.twgt * 1e3)
        self.g_as_cnt[lon, lat] += 1.0

        # Add relevant data to array:
        # QUESTION: Does cnt_loop need to be a grid
        # These are all jagged arrays; may not be the same side
        self.g_no2[lon, lat, int(self.cnt_loop[lon, lat])] = no2.no2_2d  # strat_col+trop_col   #no2_2d
        self.strat_no2[lon, lat, int(self.cnt_loop[lon, lat])] = no2.strat_col
        self.g_cld_p[lon, lat, int(self.cnt_loop[lon, lat])] = no2.t_cld_hgt[y, x]
        self.g_o3[lon, lat, int(self.cnt_loop[lon, lat])] = np.mean(no2.t_gc_o3[no2.level_min:no2.level_max + 1, y, x])
        self.g_true_no2[lon, lat, int(self.cnt_loop[lon, lat])] = np.mean(no2.t_gc_no2[no2.level_min:no2.level_max + 1, y, x])
        self.g_cld_fr[lon, lat, int(self.cnt_loop[lon, lat])] = np.sum(no2.t_cld_fr[no2.level_min:no2.level_max + 1, y, x])
        # Increment indices:
        self.cnt_loop[lon, lat] += 1

    def process_grid_square(self, i, j):
        # Define vectors of relevant data:
        # These should all have the same length
        t_col_no2 = self.g_no2[i, j, :]
        t_strat_no2 = self.strat_no2[i, j, :]
        t_fr_c = self.g_cld_fr[i, j, :]
        t_cld = self.g_cld_p[i, j, :]
        t_mr_no2 = self.g_true_no2[i, j, :]
        t_o3 = self.g_o3[i, j, :]
        # Skip if fewer than 10 points:
        if len(t_col_no2) < 10:
            print("Grid square {}, {} failed with condition 0".format(i, j))
            self.loss_count["too_few_points"] += 1
            return

        #TODO: Can change this
            # Trim to remove trailing zeros:
        t_col_no2 = np.trim_zeros(t_col_no2)
        t_strat_no2 = np.trim_zeros(t_strat_no2)
        t_fr_c = np.trim_zeros(t_fr_c)
        t_cld = np.trim_zeros(t_cld)
        t_mr_no2 = np.trim_zeros(t_mr_no2)
        t_o3 = np.trim_zeros(t_o3)
        # Remove non-uniform stratosphere:
        if (np.std(t_strat_no2) / np.mean(t_strat_no2)) > 0.02:
            self.loss_count["non_uni_strat"] += 1
            return  # continue
        # Check if there are any repeat values:
        u, ind = np.unique(t_cld, return_inverse=True)
        if (len(ind) != len(t_cld)):
            #TODO: Find out why this is bad
            print('Repeat cloud values', flush=True)
            sys.exit()
        # Get number of points:
        n_pnts = len(t_cld)
        if n_pnts > self.maxcnt:
            max_cnt = n_pnts
            print(max_cnt, flush=True)

        # Use cloud_slice_ut_no2 function to get cloud-sliced
        # UT NO2 mixing ratios.
        # Treat data differently depending on whether there are no
        # 20-100 points:

        if n_pnts >= 20 and n_pnts <100:
            self.add_slice(i,j,t_cld,t_col_no2,t_fr_c, t_mr_no2)
        elif n_pnts >= 100:
            num_slices = 40
            stride = round(n_pnts / num_slices)
            nloop = list(range(stride))
            for w in nloop:
                # Cloud-slicing:
                subset_t_col_no2 = t_col_no2[w::stride]
                subset_t_cld = t_cld[w::stride]
                subset_t_fr_c = t_fr_c[w::stride]
                subset_t_mr_no2 = t_mr_no2[w::stride]
                self.add_slice(i, j, subset_t_cld, subset_t_col_no2, subset_t_fr_c, subset_t_mr_no2)

    def add_slice(self, i, j, t_cld, t_col_no2, t_fr_c, t_mr_no2):
        """Applies and adds a cloud slice from the given data"""
        utmrno2, utmrno2err, num, mean_cld_pres = cldslice(t_col_no2, t_cld)
        # Calculate Gaussian weight:
        g_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
        # Skip if approach didn't work (i.e. cloud-sliced UT NO2 is NaN):
        # TODO: Check with E if we need to carry on with the loop in this case, or drop it
        if np.isnan(utmrno2) or np.isnan(utmrno2err):
            self.loss_count[num - 1] += 1
            # return
        # Gaussian-weighted mean for each pass of cldslice:
        gaussian_mean = np.mean(t_mr_no2 * 1e3) * g_wgt
        self.true_no2[i, j] += gaussian_mean
        self.g_cld_fr[i, j] += np.mean(t_fr_c)
        self.g_no2_vmr[i, j] += utmrno2 * g_wgt
        self.g_err[i, j] += g_wgt
        self.g_cnt[i, j] += 1
        self.cnt += 1

    def apply_gaussian_weight(self):
        """
        Apply gaussian weighting to data
        """
        self.g_no2_vmr = np.divide(self.g_no2_vmr, self.g_err, where=self.g_cnt != 0)
        self.g_cld_fr = np.divide(self.g_cld_fr, self.g_cnt, where=self.g_cnt != 0)
        self.true_no2 = np.divide(self.true_no2, self.g_err, where=self.g_cnt != 0)
        self.true_o3 = np.divide(self.true_o3, self.g_err, where=self.g_cnt != 0)
        self.g_cld = np.divide(self.g_cld, self.g_cnt, where=self.g_cnt != 0)
        self.g_err = np.divide(self.g_err, self.g_cnt, where=self.g_cnt != 0)
        self.g_cut_no2 = np.divide(self.g_cut_no2, self.true_wgt, where=self.g_as_cnt != 0)
        self.true_no2[self.g_cnt == 0] = np.nan
        self.true_o3[self.g_cnt == 0] = np.nan
        self.g_err[self.g_cnt == 0] = np.nan
        self.g_cut_no2[self.g_as_cnt == 0] = np.nan
        self.g_no2_vmr[self.g_cnt == 0] = np.nan
        self.g_cld_fr[self.g_cnt == 0] = np.nan

    # ----Reporting and saving methods----
    def print_data_report(self):
        # No. of data points:
        # print('No. of valid data points: ',cnt,flush=True)
        print('Max no. of data points in a gridsquare: ', np.amax(self.g_cnt), flush=True)
        # Track reasons for data loss:
        print('(1) Too few points: ', self.loss_count["too_few_points"], flush=True)
        print('(2) Low cloud height range: ', self.loss_count["low_cloud_height_range"], flush=True)
        print('(3) Low cloud height std dev: ', self.loss_count["low_cloud_height_std"], flush=True)
        print('(4) Large error: ', self.loss_count["large_error"], flush=True)
        print('(5) Significantly less then zero: ', self.loss_count["much_less_than_zero"], flush=True)
        print('(6) Outlier (NO2 > 200 pptv): ', self.loss_count["no2_outlier"], flush=True)
        print('(7) Non-uniform stratosphere: ', self.loss_count["non_uni_strat"], flush=True)
        print('(8) Successful retrievals: ', self.cnt, flush=True)
        print('(9) Total possible points: ', (np.sum(self.loss_count.values()) + self.cnt), flush=True)

    def plot_data(self):
        # Plot the data:
        m = Basemap(resolution='l', projection='merc',
                    lat_0=0, lon_0=0, llcrnrlon=self.minlon,
                    llcrnrlat=self.minlat, urcrnrlon=self.maxlon, urcrnrlat=self.maxlat)
        xi, yi = m(self.X, self.Y)
        plt.subplot(2, 3, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.g_no2_vmr), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Cloud-sliced NO2 VMRs')
        plt.subplot(2, 3, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.true_o3), vmin=40.0, vmax=120, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('O3 VMR [ppbv]')
        plt.subplot(2, 3, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.g_cut_no2), vmin=0, vmax=80, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('True NO2 VMRs under all-sky conditions')
        plt.subplot(2, 3, 4)
        plt.plot(self.true_no2, self.g_no2_vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.true_no2[~np.isnan(self.g_no2_vmr)],
                           self.g_no2_vmr[~np.isnan(self.g_no2_vmr)])
        # print('Correlation = ', r[0])
        result = rma(self.true_no2[~np.isnan(self.g_no2_vmr)], self.g_no2_vmr[~np.isnan(self.g_no2_vmr)],
                     len(self.true_no2[~np.isnan(self.g_no2_vmr)]), 1000)
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
        add2plt = ("y = {a:6.2f}x + {b:6.3f}".format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8,
                 ha='left', va='center')  # , transform=ax.transAxes)
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8,
                 ha='left', va='center')  # , transform=ax.transAxes)
        plt.subplot(2, 3, 5)
        plt.plot(self.g_cut_no2, self.g_no2_vmr, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.g_cut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_cut_no2))],
                           self.g_no2_vmr[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_cut_no2))])
        result = rma(self.g_cut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_cut_no2))],
                     self.g_no2_vmr[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_cut_no2))],
                     len(self.g_cut_no2[(~np.isnan(self.g_no2_vmr)) & (~np.isnan(self.g_cut_no2))]),
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
        plt.plot(self.g_cut_no2, self.true_no2, 'o', color='black', markersize=6)
        r = stats.pearsonr(self.g_cut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_cut_no2))],
                           self.true_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_cut_no2))])
        result = rma(self.g_cut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_cut_no2))],
                     self.true_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_cut_no2))],
                     len(self.g_cut_no2[(~np.isnan(self.true_no2)) & (~np.isnan(self.g_cut_no2))]),
                     1000)
        xvals = np.arange(0, 100, 5)
        yvals = result[1] + xvals * result[0]
        plt.plot(xvals, yvals, '-')
        plt.xlim(-4, 80)
        plt.ylim(-4, 80)
        plt.xlabel('True NO2 (all-sky)')
        plt.ylabel('True NO2 (cloudy)')
        add2plt = ("y = {a:6.2f}x + {b:6.3f}".format(a=result[0], b=result[1]))
        plt.text(2, 75, add2plt, fontsize=8, ha='left', va='center')
        add2plt = ("R = {a:6.2f}".format(a=r[0]))
        plt.text(2, 65, add2plt, fontsize=8, ha='left', va='center')
        plt.show()

    def save_to_netcdf(self, out_path):
        # TODO: Make variables members
        # Save the data to NetCDF:
        ncout = Dataset(out_path, mode='w', format='NETCDF4')
        # Create data dimensions:
        ncout.createDimension('lat', self.ydim)
        ncout.createDimension('lon', self.xdim)
        # create lon axis:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longitude'
        lon[:] = self.out_lon
        # create lat axis:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latgitude'
        lat[:] = self.out_lat
        # create data axes:
        # (1) Cloud-sliced NO2:
        csutno2 = ncout.createVariable('csutno2', np.float32, ('lon', 'lat'))
        csutno2.units = 'pptv'
        csutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained using cloud-slicing'
        csutno2[:] = self.g_no2_vmr
        # (2) Cloud-sliced NO2 error:
        utno2err = ncout.createVariable('utno2err', np.float32, ('lon', 'lat'))
        utno2err.units = 'pptv'
        utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2err[:] = self.g_err
        # (3) Number of observations in each gridsquare:
        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = self.g_cnt
        # (4) Mean cloud pressure for season between 450-180 hPa:
        utcld = ncout.createVariable('utcld', np.float32, ('lon', 'lat'))
        utcld.units = 'hPa'
        utcld.long_name = 'Mean cloud pressure between 450 and 180 hPa'
        utcld[:] = self.g_cld
        # (5) Mean NO2 mixing ratio at 450-180 hPa for scenes with clouds:
        cldutno2 = ncout.createVariable('cldutno2', np.float32, ('lon', 'lat'))
        cldutno2.units = 'pptv'
        cldutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained if clouds are present'
        cldutno2[:] = self.true_no2
        # (6) Mean NO2 mixing ratio at 450-180 hPa under all conditions (all-sky):
        askutno2 = ncout.createVariable('askutno2', np.float32, ('lon', 'lat'))
        askutno2.units = 'pptv'
        askutno2.long_name = 'UT NO2 mixing ratio (180-450 hPa) obtained under all conditions (all-sky)'
        askutno2[:] = self.g_cut_no2
        # (7) Cloud fraction:
        utcldfrc = ncout.createVariable('utcldfrc', np.float32, ('lon', 'lat'))
        utcldfrc.units = 'unitless'
        utcldfrc.long_name = 'GEOS-FP cloud fraction obtained as sum of 3D cloud fractions across range of interest (180-450 hPa)'
        utcldfrc[:] = self.g_cld_fr
        # (8) O3 sampled coincident with cloud-slicing retrieval:
        uto3 = ncout.createVariable('uto3', np.float32, ('lon', 'lat'))
        uto3.units = 'ppbv'
        uto3.long_name = 'GEOS-Chem ozone obtained coincident with cloud-sliced NO2'
        uto3[:] = self.true_o3
        # Close the file:
        ncout.close()


class GeosChemDay:

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
        self.t_lon = tlon[:]
        self.t_lat = tlat[:]
        self.t_gc_no2 = tgcno2[:]
        self.t_cld_fr = tcldfr[:]
        self.t_cld_hgt = tcldhgt[0, :, :]
        self.t_adn = tadn[:]  # in molec/cm3
        self.t_bx_hgt = tbxhgt[:]
        self.t_p_edge = tpedge[:]
        self.t_pause = tpause[0, :, :]
        self.t_gc_o3 = tgco3[:]
        # Convert box height from m to cm:
        self.t_bx_hgt = tbxhgt * 1e2

        #Get outputs ready here for tidyness:
        self.no2_2d = None
        self.trop_col = None
        self.strat_col = None
        self.gcutno2 = None
        self.gascnt

        self.level_min = None
        self.level_max = None
        self.askind = None


    def prepare_no2_pixel(self, x, y):

        # Calculate corresponding mid-pressure values:
        tp_mid = np.zeros(len(self.t_p_edge[:, y, x]))
        # Get mid-pressure values, except for highest layer:
        for k in range(len(self.t_p_edge[:, y, x]) - 1):
            tp_mid[k] = np.multiply(0.5, (self.t_p_edge[k, y, x] + self.t_p_edge[k + 1, y, x]))
        # Estimate mid-pressure for highest value (doesn't need to
        # be accurate, as surpassing the range of interset):
        # TODO: Check with E why 46. Middle of tpmid?
        tp_mid[46] = np.multiply(0.5, (self.t_p_edge[46, y, x] + (self.t_p_edge[46, y, x] - 0.1)))
        # Get model layer of tropopause:
        tppind = np.argmin(abs(tp_mid - self.t_pause[y, x]))
        # Get indices that fall between 450 and 180 hPa for estimating
        # "true' all-sky UT NO2 and partial columns:
        lind = np.where((tp_mid >= P_MIN) & (tp_mid <= P_MAX))[0]
        # Get UT NO2 under "true" all-sky conditions:
        # Make sure this is below the tropopause:
        # If below tropopause, use full extent (180-450 hPa):
        if lind[len(lind) - 1] <= tppind:
            self.askind = lind
        # If above tropopause, trim to tropopause-450 hPa:
        if lind[len(lind) - 1] > tppind:
            self.askind = lind[np.where(lind <= tppind)[0]]
        # If tropopause below 450 hPa, skip entirely:
        if self.t_pause[y, x] > P_MAX:
            return  # continue
        # Get Guassian weights that allocated higher weights to points
        # closest to the pressure centre (315 hPa):
        # Equation is:
        #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
        #         135 hPa is the standard deviation.
        self.twgt = np.exp(np.multiply((-1.0), np.divide((np.square(np.subtract(tp_mid[self.askind], 315.))),
                                                         (np.multiply(2, np.square(135.))))))

        # Find where cloud fraction in UT exceeds 0.7 after calculating
        # true all-sky NO2:
        # (Keep for testing effect of thick clouds on cloud-sliced UT NO2):
        # if (np.sum(t_cld_fr[lind,y,x])<=0.7): continue
        # Get model level of cloud top height closest to lowest
        # pressure-altitude of interest (P_MIN):
        self.lcld = np.argmin(abs(self.t_cld_hgt[y, x] - tp_mid))
        # Skip if cloud top height ouside pressure range of interest:
        self.level_min, self.level_max = np.amin(lind), np.amax(lind)
        if (self.lcld <self.level_min) or (self.lcld > self.level_max):
            pass  # continue
        # This error check is probably redundant, but included in case:
        if self.lcld == 0:
            print("No cloud detected!!!", flush=True)
            sys.exit()
        # Get partial NO2 column in molec/m2 from cloud top height
        # to highest model level (output up to level 47):
        # print(t_gc_no2[self.level_min:tppind,y,x])
        # print(t_gc_no2[self.level_min:tppind,y,x]*1.5)
        self.no2_2d = np.sum(self.t_gc_no2[self.level_min:, y, x] * 1e-5 * self.t_adn[self.level_min:, y, x] *
                             self.t_bx_hgt[self.level_min:, y, x])
        # Get stratospheric column from 180 hPa aloft:
        # Previous approach (remove when model simulations done):
        # tppind=np.where(tpmid<180.)[0]
        self.strat_col = np.sum(self.t_gc_no2[tppind:, y, x] * 1e-5 * self.t_adn[tppind:, y, x] *
                                self.t_bx_hgt[tppind:, y, x])


def get_file_list(gcdir, REGION, YEARS_TO_PROCESS):
    # TODO: Rework after you've looked at the files
    files = glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + '.' + YEARS_TO_PROCESS[0] + '06*')
    mon = ["07", "08"]
    for i in mon:
        for filename in glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + \
                                  '.' + YEARS_TO_PROCESS[0] + i + '*'):
            files.append(filename)
    # 2017:
    mon = ["06", "07", "08"]
    for i in mon:
        for filename in glob.glob(gcdir + 'nc_sat_files_47L/ts_12_15.' + REGION + \
                                  '.' + YEARS_TO_PROCESS[1] + i + '*'):
            files.append(filename)
    return sorted(files)


if __name__ == "__main__":

    # Decide on region:
    REGION = 'EU'  # 'NA', 'EU', or 'CH'

    # Define information for grid:
    STR_RES = '4x5'

    parser = argparse.ArgumentParser()
    parser.add_argument("gc_dir")
    parser.add_argument("out_dir")
    parser.add_argument('resolution', help="Can be 8x10, 4x5, 2x25 or 1x1")
    parser.add_argument("-p", "--plot")
    args = parser.parse_args()

    if len(YEARS_TO_PROCESS) == 1:
        yrrange = YEARS_TO_PROCESS[0]
    if len(YEARS_TO_PROCESS) == 2:
        yrrange = '2016-2017'

    # Name of log file to output code prorgess:
    # TODO: Maybe replace this with Python logging? See rest of code.
    log = open("log_" + REGION + "_" + STR_RES + "_" + yrrange + "_v4", "w")
    sys.stdout = log

   # Get files (test June 2016 for now)
    # 2016:
    gc_dir = args.gc_dir
    files = get_file_list(gc_dir)
    print('Number of files:', len(files), flush=True)

    rolling_total = ProcessedData(REGION, STR_RES)

    # Loop over files:
    for file_path in files:
        rolling_total.process_geochem_day(file_path)

    rolling_total.apply_gaussian_weight()
    rolling_total.print_data_report()
    rolling_total.plot_data()
    out_path = './Data/gc-v12-1-0-ut-self-' + REGION.lower() + '-jja-' + yrrange + '-' + STR_RES + '-' + rolling_total.prange + '-v4.nc'
    rolling_total.save_to_netcdf(out_path)

    # Close the log file:
    log.close()
