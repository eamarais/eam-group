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
from netCDF4 import Dataset
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import argparse

from constants import AVOGADRO
from constants import G
from constants import MW_AIR

from bootstrap import rma
from cloud_slice_ut_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM


# Turn off warnings:
np.warnings.filterwarnings('ignore')




# Define pressure range:
# These are fixed, so can be constants, as opposed to inputs.
P_MIN=180     
P_MAX=450     


# Define years of interest:
YEARS_TO_PROCESS=['2016', '2017']


class ProcessingException(Exception):
    pass


class CloudSliceException(Exception):
    pass


class InvalidRegionException(Exception):
    pass


class InvalidResolutionException(Exception):
    pass

class DomainIssueException(Exception):
    pass



class ProcessedData:
    # ----Initialisation methods----
    def __init__(self, region, str_res, do_temperature_correction=False, do_error_weighting=False):

        self.temperature_correction = do_temperature_correction
        self.error_weight = do_error_weighting

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
        self.cloud_slice_count = 0
        self.maxcnt = 0

        # Define string to represent the layer range:
        self.prange = str(P_MIN) + '-' + str(P_MAX)

        # Define factor to convert slope of NO2 mixing ratio versus pressure
        # to VMR:
        self.den2mr = np.divide((np.multiply(G, MW_AIR)), AVOGADRO)

    def define_grid(self, region, str_res):
        # Define target grid:
        if region == 'NA':
            self.minlat = 4.
            self.maxlat = 60.
            self.minlon = -137.
            self.maxlon = -58.
            self.dirreg = '_na_'
        elif region == 'EU':
            self.minlat = 28.
            self.maxlat = 64.
            self.minlon = -22.
            self.maxlon = 42.
            self.dirreg = '_eu_naei_'
        elif region == 'CH':
            self.minlat = 8.
            self.maxlat = 58.
            self.minlon = 63.
            self.maxlon = 138.
            self.dirreg = '_ch_'
        else:
            print("Invalid region; valid regions are 'NA','EU','CH'.")
            raise InvalidRegionException
        # NOTE FOR LATER: This snippet turns up in fresco_cld_err; a candidate for the library.
        # Define grid information:
        if str_res == '8x10':
            self.dellat, self.dellon = 8, 10
        elif str_res == '4x5':
            self.dellat, self.dellon = 4, 5
        elif str_res == '2x25':
            self.dellat, self.dellon = 2, 2.5
        elif str_res == '1x1':
            self.dellat, self.dellon = 1, 1
        else:
            print("Invalid resolution: valid resolutions are 8x10, 4x5, 2x25 (two pt 5) and 1x1")
            raise InvalidResolutionException
        self.out_lon = np.arange(self.minlon, self.maxlon + self.dellon, self.dellon)
        self.out_lat = np.arange(self.minlat, self.maxlat + self.dellat, self.dellat)
        # Convert output lats and long to 2D:
        self.X, self.Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')
        # Dimensions of output data:
        self.xdim = len(self.out_lon)
        self.ydim = len(self.out_lat)
        # Get maximum possible number of native resolution model gridsquares
        # in the coarser grid (to use in process_grid_square to check that 
        # don't exceed this):
        self.max_limit=(self.dellon/0.3125)*(self.dellat/0.25)


    # ----Processing methods----
    def process_geoschem_day(self, file_path):

        # I plonked this here to force these arrays to be reinitated each day with each new
        # file. Please relocate if this could be more eloquently done. - No, this is good; makes a lot more sense here.
        # Define output data for this day:
        out_shape = (self.xdim, self.ydim)  # This feel gross. A 3-d list of appendable lists.
        self.g_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_o3 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.all_cld_fr = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.strat_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_cld_p = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.g_true_no2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]

        this_geoschem_day = GeosChemDay(file_path,
                                        error_weight=self.error_weight,
                                        temperature_correction=self.temperature_correction)
        # Get column values:
        for y in range(len(this_geoschem_day.t_lat)):
            for x in range(len(this_geoschem_day.t_lon)):
                this_geoschem_day.prepare_no2_pixel(x, y)
                self.regrid_and_process(x, y, this_geoschem_day)
        for i in range(self.xdim):
            for j in range(self.ydim):
                if any(value != 0 for value in self.g_no2[i][j]):
                    self.process_grid_square(i, j)

    def regrid_and_process(self, x, y, no2):

        # If no valid data, skip:
        if ( len(no2.askind) == 0 ):
            return

        # Find nearest gridsquare in output grid:
        lon = int(np.argmin(abs(self.out_lon - no2.t_lon[x])))
        lat = int(np.argmin(abs(self.out_lat - no2.t_lat[y])))

        self.true_wgt[lon, lat] += np.sum(no2.twgt)

        self.g_cut_no2[lon, lat] += np.sum(no2.t_gc_no2[no2.askind, y, x] * no2.twgt * 1e3)
        self.g_as_cnt[lon, lat] += 1.0

        # Add relevant data to array:
        # QUESTION: Does cnt_loop need to be a grid
        # These are all jagged arrays; may not be the same side

        # Only proceed beyond this step if there are clouds at 450-180 hPa:
        if (no2.lcld <no2.level_min) or (no2.lcld > no2.level_max):
            # Should this be return rather than pass? This checks for clouds between min and mas pressure. If none, move to next pixel.
            #print("Cloud top outside pressure range in pixel {},{}".format(x,y))
            return

        self.g_no2[lon][lat].append(no2.no2_2d)  # strat_col+trop_col   #no2_2d
        self.strat_no2[lon][lat].append(no2.strat_col)
        self.g_cld_p[lon][lat].append(no2.t_cld_hgt[y, x])
        self.g_o3[lon][lat].append(
            np.mean(no2.t_gc_o3[no2.level_min:no2.level_max + 1, y, x]))
        self.g_true_no2[lon][lat].append(
            np.mean(no2.t_gc_no2[no2.level_min:no2.level_max + 1, y, x]))
        self.all_cld_fr[lon][lat].append(
            np.sum(no2.t_cld_fr[no2.level_min:no2.level_max + 1, y, x]))
        pass

    def process_grid_square(self, i, j):
        # Define vectors of relevant data:
        # These should all have the same length
        t_col_no2 = np.asarray(self.g_no2[i][j],dtype=np.float)
        t_strat_no2 = np.asarray(self.strat_no2[i][j],dtype=np.float)
        t_fr_c = np.asarray(self.all_cld_fr[i][j],dtype=np.float)
        t_cld = np.asarray(self.g_cld_p[i][j],dtype=np.float)
        t_mr_no2 = np.asarray(self.g_true_no2[i][j],dtype=np.float)
        t_o3 = np.asarray(self.g_o3[i][j],dtype=np.float)

        # Skip if fewer than 10 points:
        if len(t_col_no2) < 10:
            #print("Grid square {}, {} failed with condition too_few_points".format(i, j))
            self.loss_count["too_few_points"] += 1
            return

        # Remove non-uniform stratosphere:
        if (np.std(t_strat_no2) / np.mean(t_strat_no2)) > 0.02:
            self.loss_count["non_uni_strat"] += 1
            return  # continue

        # Get number of points:
        n_pnts = len(t_cld)
        if n_pnts > self.maxcnt:
            self.maxcnt = n_pnts
            print(self.maxcnt, flush=True)

        # Check if number of gridsquares exceeds max possible
        # (indicates error in coarser grid domain):
        if ( self.maxcnt > self.max_limit ):
            print('No. of model grids exceeds max possible')
            raise DomainIssueException

        # Use cloud_slice_ut_no2 function to get cloud-sliced
        # UT NO2 mixing ratios.
        # Treat data differently depending on whether there are no
        # 20-100 points:
        #try:
        if n_pnts >= 20 and n_pnts <100:
            self.add_slice(i,j,t_cld,t_col_no2, t_mr_no2, t_fr_c)
        elif n_pnts >= 100:
            num_slices = 40
            stride = round(n_pnts / num_slices)
            nloop = list(range(stride))
            for w in nloop:
                subset_t_col_no2 = t_col_no2[w::stride]
                subset_t_cld = t_cld[w::stride]
                subset_t_mr_no2 = t_mr_no2[w::stride]
                subset_t_fr_c = t_fr_c[w::stride]
                self.add_slice(i, j, subset_t_cld, subset_t_col_no2, subset_t_mr_no2, subset_t_fr_c)

    def add_slice(self, i, j, t_cld, t_col_no2, t_mr_no2, t_fr_c):
        """Applies and adds a cloud slice from the given data"""
        utmrno2, utmrno2err, stage_reached, mean_cld_pres = cldslice(t_col_no2, t_cld)
        # Calculate Gaussian weight:
        if self.error_weight:
            g_wgt = 1.0 / (utmrno2err ** 2)
        else:
            g_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
        # Skip if approach didn't work (i.e. cloud-sliced UT NO2 is NaN):
        # Yes, drop out after the reason for data loss is added to loss_count.
        if np.isnan(utmrno2) or np.isnan(utmrno2err):
            self.loss_count[CLOUD_SLICE_ERROR_ENUM[stage_reached]] += 1
            #print("Cloud-slice exception {} in pixel i:{} j:{}".format(
            #    CLOUD_SLICE_ERROR_ENUM[stage_reached], i, j
            #))
        else:
            # Gaussian-weighted mean for each pass of cldslice:
            gaussian_mean = np.mean(t_mr_no2 * 1e3) * g_wgt
            self.true_no2[i, j] += gaussian_mean
            self.g_no2_vmr[i, j] += utmrno2 * g_wgt
            self.g_err[i, j] += g_wgt
            self.g_cnt[i, j] += 1
            self.g_cld_fr[i, j] += np.mean(t_fr_c)
            self.cloud_slice_count += 1

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
        # print('No. of valid data points: ',cloud_slice_count,flush=True)
        print('Max no. of data points in a gridsquare: ', np.amax(self.g_cnt), flush=True)
        # Track reasons for data loss:
        print('(1) Too few points: ', self.loss_count["too_few_points"], flush=True)
        print('(2) Low cloud height range: ', self.loss_count["low_cloud_height_range"], flush=True)
        print('(3) Low cloud height std dev: ', self.loss_count["low_cloud_height_std"], flush=True)
        print('(4) Large error: ', self.loss_count["large_error"], flush=True)
        print('(5) Significantly less then zero: ', self.loss_count["much_less_than_zero"], flush=True)
        print('(6) Outlier (NO2 > 200 pptv): ', self.loss_count["no2_outlier"], flush=True)
        print('(7) Non-uniform stratosphere: ', self.loss_count["non_uni_strat"], flush=True)
        print('(8) Successful retrievals: ', self.cloud_slice_count, flush=True)
        print('(9) Total possible points: ', (sum(self.loss_count.values()) + self.cloud_slice_count), flush=True)

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
    def __init__(self, file_path, error_weight=False, temperature_correction=False):
        print(file_path, flush=True)

        self.error_weight = error_weight
        self.temperature_correction = temperature_correction

        # Read dataset:
        fh = Dataset(file_path, mode='r')
        # Extract data of interest:
        # (Add tropopause height to this in the future)
        tlon, tlat, tgcno2, tcldfr, tcldhgt, tadn, tbxhgt, tpedge, tpause, tgco3, tdegk = \
            fh.variables['LON'], fh.variables['LAT'], \
            fh.variables['IJ-AVG-S__NO2'], fh.variables['TIME-SER__CF'], \
            fh.variables['TIME-SER__CThgt'], fh.variables['TIME-SER__AIRDEN'], \
            fh.variables['BXHGHT-S__BXHEIGHT'], fh.variables['PEDGE-S__PSURF'], \
            fh.variables['TR-PAUSE__TP-PRESS'], fh.variables['IJ-AVG-S__O3'], \
            fh.variables['DAO-3D-S__TMPU']
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
        self.t_deg_k = tdegk[:]
        # Convert box height from m to cm:
        self.t_bx_hgt = self.t_bx_hgt * 1e2

        #Get outputs ready here for tidyness:
        self.no2_2d = None
        self.trop_col = None
        self.strat_col = None
        self.gcutno2 = None
        self.gascnt = None

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
        # Data output from the model includes 47 vertical layers. This means that only 46 pressure centres can be calculated as the calculation requires pressure edges.
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
            #print("Tropopause less than P_MAX in geoschem pixel x:{}, y:{}".format(x,y))
            return  # continue
        # Get Guassian weights that allocated higher weights to points
        # closest to the pressure centre (315 hPa):
        # Equation is:
        #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
        #         135 hPa is the standard deviation.
        # The "shorthand" formula (np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))) can be used here too
        if self.error_weight:
            self.twgt = np.ones(len(self.askind))
        else:
            self.twgt = np.exp((-(tp_mid[self.askind] - 315) ** 2) / (2 * 135 ** 2))

        # Find where cloud fraction in UT exceeds 0.7 after calculating
        # true all-sky NO2:
        # (Keep for testing effect of thick clouds on cloud-sliced UT NO2):
        # if (np.sum(t_cld_fr[lind,y,x])<=0.7): continue
        # Get model level of cloud top height closest to lowest
        # pressure-altitude of interest (P_MIN):
        self.lcld = np.argmin(abs(self.t_cld_hgt[y, x] - tp_mid))
        # Skip if cloud top height ouside pressure range of interest:
        self.level_min, self.level_max = np.amin(lind), np.amax(lind)

        if (self.temperature_correction):

            # Equation is from the TROPOMI product ATBD (p. 32, Eqn 18)
            # (product document abbrevation: S5P-KNMI-L2-0005-RP)
            self.temp_corr = 1 - (3.16e-3 * (self.t_deg_k[self.level_min:, y, x] - 220.)) + \
                        (3.39e-6 * ((self.t_deg_k[self.level_min:, y, x] - 220) ** 2))
        else:
            # Set to 1 so that no scaling is applied:
            # (might be a more eloquent way to do this)
            self.temp_corr = np.ones(len(self.t_gc_no2[self.level_min:, y, x]))

        # Get partial NO2 column in molec/m2 from cloud top height
        # to highest model level (output up to level 47):
        # print(t_gc_no2[self.level_min:tppind,y,x])
        # print(t_gc_no2[self.level_min:tppind,y,x]*1.5)
        self.no2_2d = np.sum(self.t_gc_no2[self.level_min:, y, x]
                             * 1e-5
                             * self.temp_corr
                             * self.t_adn[self.level_min:, y, x]
                             * self.t_bx_hgt[self.level_min:, y, x])
        # Get stratospheric column from 180 hPa aloft:
        # Previous approach (remove when model simulations done):
        # tppind=np.where(tpmid<180.)[0]
        self.strat_col = np.sum(self.t_gc_no2[tppind:, y, x]
                                * 1e-5 * self.t_adn[tppind:, y, x]
                                * self.t_bx_hgt[tppind:, y, x])


def get_file_list(gcdir, REGION, YEARS_TO_PROCESS):
    # Define target grid:
    if REGION == 'NA':
        dirreg = '_na_'
    elif REGION == 'EU':
        dirreg = '_eu_naei_'
    elif REGION == 'CH':
        dirreg = '_ch_'
    else:
        print("Invalid region; valid regions are 'NA','EU','CH'.")
        raise InvalidRegionException
        # NOTE FOR LATER: This snippet turns up in fresco_cld_err; a candidate for the library.
    gcdir=gcdir+'geosfp'+dirreg+'iccw/'
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
    #REGION = 'EU'  # 'NA', 'EU', or 'CH'

    # Define information for grid:
    #STR_RES = '4x5'

    parser = argparse.ArgumentParser()
    # Shorten directory name to up to "GC/", then define the subdirectory
    # as 'geosfp' + dirreg + 'iccw/' in get_file_list.
    parser.add_argument("--gc_dir", default='/data/uptrop/Projects/DEFRA-NH3/GC/')
    parser.add_argument("--out_path", default='/home/j/jfr10/eos_library/uptrop_comparison/test.nc2')
    parser.add_argument('--resolution', default="4x5", help="Can be 8x10, 4x5, 2x25 or 1x1")
    parser.add_argument('--region', default="EU", help="Can be EU, NA, or CH")
    parser.add_argument("-p", "--plot", type=bool)
    parser.add_argument("--do_temp_correct", type=bool)
    parser.add_argument("--do_error_weight", type=bool)
    args = parser.parse_args()

    if len(YEARS_TO_PROCESS) == 1:
        yrrange = YEARS_TO_PROCESS[0]
    if len(YEARS_TO_PROCESS) == 2:
        yrrange = '2016-2017'

    # Name of log file to output code prorgess:
    # TODO: Maybe replace this with Python logging? See rest of code.
    #log = open("log_" + REGION + "_" + STR_RES + "_" + yrrange + "_v4", "w")
    #sys.stdout = log

   # Get files (test June 2016 for now)
    # 2016:
    gc_dir = args.gc_dir
    STR_RES = args.resolution
    REGION = args.region
    files = get_file_list(gc_dir, REGION, YEARS_TO_PROCESS)
    print('Number of files:', len(files), flush=True)

    rolling_total = ProcessedData(REGION, STR_RES,
                                  do_temperature_correction=args.do_temp_correct,
                                  do_error_weighting=args.do_error_weight)

    # Loop over files:
    for file_path in files:
        rolling_total.process_geoschem_day(file_path)

    rolling_total.apply_gaussian_weight()
    rolling_total.print_data_report()
    rolling_total.plot_data()
    out_path = args.out_path
    rolling_total.save_to_netcdf(out_path)

    # Close the log file:
    #log.close()

# Define years of interest:
stryr=['2016','2017']
if len(stryr)==1: yrrange=stryr[0]
if len(stryr)==2: yrrange='2016-2017'

# Apply temperature correction?
# (Can be added as an input argument; default is False)
temperature_correction=False
error_weight=False

# Name of log file to output code prorgess:
log=open("log_"+Reg+"_"+StrRes+"_"+yrrange+"_v4", "w")
sys.stdout=log

if ( temperature_correction ):
    print('Temperature correction is being applied', flush=True)

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
trueo3=np.zeros(X.shape)   # ozone mixing ratio coincident with cloud-sliced ut no2
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

# Order the files:
files=sorted(files)

# Check:
print('Number of files:',len(files),flush=True)

# Loop over files:
for f in files:
    print(f,flush=True)

    # Read dataset:
    fh=Dataset(f,mode='r')

    # Extract data of interest:
    # (Add tropopause height to this in the future)
    tlon,tlat,tgcno2,tcldfr,tcldhgt,tadn,tbxhgt,tpedge,tpause,tgco3,tdegk=\
          fh.variables['LON'],fh.variables['LAT'],\
          fh.variables['IJ-AVG-S__NO2'],fh.variables['TIME-SER__CF'],\
          fh.variables['TIME-SER__CThgt'],fh.variables['TIME-SER__AIRDEN'],\
          fh.variables['BXHGHT-S__BXHEIGHT'],fh.variables['PEDGE-S__PSURF'],\
          fh.variables['TR-PAUSE__TP-PRESS'],fh.variables['IJ-AVG-S__O3'],\
          fh.variables['DAO-3D-S__TMPU']
    tlon=tlon[:]
    tlat=tlat[:]
    tgcno2=tgcno2[:]
    tcldfr=tcldfr[:]
    tcldhgt=tcldhgt[0,:,:]
    tadn=tadn[:]   # in molec/cm3
    tbxhgt=tbxhgt[:]
    tpedge=tpedge[:]
    tpause=tpause[0,:,:]
    tgco3=tgco3[:]
    tdegk=tdegk[:]    # 3D temperature in deg K

    # Convert box height from m to cm:
    tbxhgt=tbxhgt*1e2

    # Define output data for this day:
    gno2=np.zeros((xdim,ydim,nval))
    go3=np.zeros((xdim,ydim,nval))
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

            # Get model layer of tropopause:
            tppind=np.argmin(abs(tpmid-tpause[y,x])) 

            # Get indices that fall between 450 and 180 hPa for estimating
            # "true' all-sky UT NO2 and partial columns:
            lind=np.where((tpmid>=pmin)&(tpmid<=pmax))[0]

            # Get UT NO2 under "true" all-sky conditions:
            # Make sure this is below the tropopause:
            # If below tropopause, use full extent (180-450 hPa):
            if ( lind[len(lind)-1]<=tppind ): askind=lind
            # If above tropopause, trim to tropopause-450 hPa:
            if ( lind[len(lind)-1]>tppind ): 
                askind=lind[np.where(lind<=tppind)[0]]
            # If tropopause below 450 hPa, skip entirely:
            if ( tpause[y,x]>pmax ): continue

            # Get Guassian weights that allocated higher weights to points
            # closest to the pressure centre (315 hPa):
            # Equation is:
            #   w = exp(-(p-315)^2/2*135^2 ) where 315 hPa is the centre and
            #         135 hPa is the standard deviation.
            if ( error_weight ):
                twgt=np.ones(len(askind))            
            else:
                twgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(tpmid[askind],315.))),(np.multiply(2,np.square(135.))))))

            # Sum up "true" all-sky UT NO2:
            gcutno2[p,q]=gcutno2[p,q] + \
                          np.sum(tgcno2[askind,y,x]*twgt*1e3)
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

            # Calculate temperature correction to test its effect on cloud-sliced
            # UT NO2 (no temperature correction is applied to the AMF used to 
            # calculated TROPOMI partial NO2 columns, as 3D temperature isn't included
            # in the TROPOMI file. This serves as a test of the effect of this
            # temperature correction on the synthetic results):
            if ( temperature_correction ):

                # Equation is from the TROPOMI product ATBD (p. 32, Eqn 18)
                # (product document abbrevation: S5P-KNMI-L2-0005-RP)
                temp_corr=1-(3.16e-3*(tdegk[lmin:,y,x] - 220.)) + \
                           (3.39e-6*((tdegk[lmin:,y,x] - 220)**2))
            else:
                # Set to 1 so that no scaling is applied:
                # (might be a more eloquent way to do this)
                temp_corr=np.ones(len(tgcno2[lmin:,y,x]))

            # Get partial NO2 column in molec/m2 from cloud top height 
            # to highest model level (output up to level 47):
            no22d=np.sum(tgcno2[lmin:,y,x]*1e-5*temp_corr*tadn[lmin:,y,x]*\
                         tbxhgt[lmin:,y,x])                

            # No longer used, so can be commented out or deleted
            #tropcol=np.sum(tgcno2[lmin:tppind,y,x]*1e-5*\
            #               tadn[lmin:tppind,y,x]*tbxhgt[lmin:tppind,y,x])

            # Get stratospheric column from 180 hPa aloft:
            # Previous approach (remove when model simulations done):
            #tppind=np.where(tpmid<180.)[0]

            stratcol=np.sum(tgcno2[tppind:,y,x]*1e-5*tadn[tppind:,y,x]*\
                            tbxhgt[tppind:,y,x])

            # Add relevant data to array:
            gno2[p,q,int(cntloop[p,q])]=no22d     #stratcol+tropcol   #no22d
            stratno2[p,q,int(cntloop[p,q])]=stratcol
            gcldp[p,q,int(cntloop[p,q])]=tcldhgt[y,x]
            go3[p,q,int(cntloop[p,q])]=np.mean(tgco3[lmin:lmax+1,y,x])
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
                to3=go3[i,j,:]

                # Skip if fewer than 10 points:
                if len(tcolno2)<10: 
                    num[0]=num[0]+1
                    continue                

                # Trim to remove trailing zeros:
                tcolno2=np.trim_zeros(tcolno2)#[np.nonzero(tcolno2)]
                tstratno2=np.trim_zeros(tstratno2)#[np.nonzero(tstratno2)]
                tfrc=np.trim_zeros(tfrc)#[np.nonzero(tfrc)]
                tcld=np.trim_zeros(tcld)#[np.nonzero(tcld)]
                tmrno2=np.trim_zeros(tmrno2)#[np.nonzero(tmrno2)]
                to3=np.trim_zeros(to3)#[np.nonzero(tmrno2)]

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

                # Use cloud_slice_ut_no2 function to get cloud-sliced 
                # UT NO2 mixing ratios.
                # Treat data differently depending on whether there are no
                # more than or more than 30 points:
                if ((npnts>=20)&(npnts<100)): #<=30:
                    csval=cldslice(tcolno2,tcld)
                    if np.isnan(csval[0]) or np.isnan(csval[1]):
                        num[csval[2]-1]=num[csval[2]-1]+1
                        continue

                    # Get mean cloud top pressure:
                    pmean=csval[3]

                    # Calculate weights:
                    if ( error_weight ):
                        gwgt=1.0/(csval[1]**2)
                    else:
                        gwgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(pmean,315.))),(np.multiply(2,np.square(135.))))))

                    # Apply Gaussian weights to cloud-sliced UT NO2:
                    gno2vmr[i,j]=gno2vmr[i,j]+np.multiply(csval[0],gwgt)
                    gcldfr[i,j]=gcldfr[i,j]+np.mean(tfrc)
                    gerr[i,j]=gerr[i,j]+gwgt
               
                    # "True" cloudy UT NO2 (converted from ppbv to pptv):
                    trueno2[i,j]=trueno2[i,j]+\
                                  np.multiply(np.mean(tmrno2*1e3),gwgt)

                    # Ozone coincident with cloud-sliced NO2:
                    trueo3[i,j]=trueo3[i,j]+\
                                 np.multiply(np.mean(to3),gwgt)
                    
                    gcnt[i,j]=gcnt[i,j]+1
                    cnt=cnt+1

                elif (npnts>=100):
                #    # Code error check:
                #    if len(tcld)!=len(tmrno2): 
                #        print('arrays not =',flush=True)
                #        print(len(tcld),len(tmrno2),flush=True)

                #    # Determine the number of iterations:
                    niter=round(npnts/40)
                    nloop=list(range(niter))

                #    # Loop over iterations:
                    for w in nloop:
                        if w==0:
                            
                            # Cloud-slicing:
                            csval=cldslice(tcolno2[::niter],tcld[::niter])

                            # Get mean cloud top pressure:
                            pmean=csval[3]

                            # Calculate weights:
                            if ( error_weight ):
                                gwgt=1.0/(csval[1]**2)
                            else:
                                gwgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(pmean,315.))),(np.multiply(2,np.square(135.))))))
                        else:
                            
                            # Cloud-slicing:
                            csval=cldslice(tcolno2[w::niter],tcld[w::niter])

                            # Get mean cloud top pressure:
                            pmean=csval[3]

                            # Calculate weights:
                            if ( error_weight ):
                                gwgt=1.0/(csval[1]**2)
                            else:
                                gwgt=np.exp(np.multiply((-1.0),np.divide((np.square(np.subtract(pmean,315.))),(np.multiply(2,np.square(135.))))))

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
                        gcldfr[i,j]=gcldfr[i,j]+np.mean(tfrc[w::niter])#

                    gcnt[i,j]=gcnt[i,j]+1
                    cnt=cnt+1

# No. of data points:                            
#print('No. of valid data points: ',cnt,flush=True)
print('Max no. of data points in a gridsquare: ',np.amax(gcnt),flush=True)

# Track reasons for data loss:
print('(1) Too few points: ',num[0],flush=True)
print('(2) Low cloud height range: ',num[1],flush=True)
print('(3) Low cloud height std dev: ',num[2],flush=True)
print('(4) Large error: ',num[3],flush=True)
print('(5) Significantly less then zero: ',num[4],flush=True)
print('(6) Outlier (NO2 > 200 pptv): ',num[5],flush=True)
print('(7) Non-uniform stratosphere: ',num[6],flush=True)
print('(8) Successful retrievals: ',cnt,flush=True)
print('(9) Total possible points: ',(np.sum(num)+cnt),flush=True)

# Get average:
gno2vmr=np.divide(gno2vmr, gerr, where=gcnt!=0)
gcldfr=np.divide(gcldfr, gcnt, where=gcnt!=0)
trueno2=np.divide(trueno2, gerr, where=gcnt!=0)
trueo3=np.divide(trueo3, gerr, where=gcnt!=0)
gcld=np.divide(gcld, gcnt, where=gcnt!=0)
gerr=np.divide(gerr, gcnt, where=gcnt!=0)
gcutno2=np.divide(gcutno2,truewgt,where=gascnt!=0)
trueno2[gcnt==0]=np.nan
trueo3[gcnt==0]=np.nan
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
cs=m.pcolor(xi,yi,np.squeeze(trueo3), vmin=40.0, vmax=120, cmap='jet')
m.drawparallels(np.arange(-80.,81.,45.),labels=[1,0,0,0],fontsize=8)
m.drawmeridians(np.arange(-180.,181.,45.),labels=[0,0,0,1],fontsize=8)
m.drawcoastlines()
m.drawcountries()
cbar = m.colorbar(cs, location='bottom', pad="10%")
plt.title('O3 VMR [ppbv]')

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
              StrRes+'-'+prange+'-v4.nc',mode='w',format='NETCDF4')

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

# (8) O3 sampled coincident with cloud-slicing retrieval:
uto3 = ncout.createVariable('uto3', np.float32, ('lon','lat'))
uto3.units = 'ppbv'
uto3.long_name = 'GEOS-Chem ozone obtained coincident with cloud-sliced NO2'
uto3[:] = trueo3

# Close the file:
ncout.close()

# Close the log file:
log.close()
