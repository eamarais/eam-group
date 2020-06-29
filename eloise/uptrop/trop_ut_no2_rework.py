import glob
import argparse
import sys
import os
from os import path
from netCDF4 import Dataset
import datetime as dt
import re

import numpy as np
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from dateutil import rrule as rr

# Import hack
sys.path.append(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        '..'))

from uptrop.convert_height_to_press import alt2pres
from uptrop.cloud_slice_ut_no2 import cldslice, CLOUD_SLICE_ERROR_ENUM


class CloudFileDateMismatch(Exception):
    pass


class CloudFileShapeMismatch(Exception):
    pass


class GridAggregator:
    def __init__(self, dellat, dellon):
        self.postfilt = []
        self.out_lon = np.arange(-180, 180 + dellon, dellon)
        self.out_lat = np.arange(-90, 90 + dellat, dellat)

        self.xdim = len(self.out_lon)
        self.ydim = len(self.out_lat)

        # Seasonal mean output arrays
        self.gno2vmr = np.zeros((self.xdim, self.ydim))  # NO2 VMR in pptv
        self.gcnt = np.zeros((self.xdim, self.ydim))  # No. of data points (orbits)
        self.gerr = np.zeros((self.xdim, self.ydim))  # Weighted error

        self.file_count = 0

        self.loss_count = {
            "too_few_points": 0,
            "low_cloud_height_range": 0,
            "low_cloud_height_std": 0,
            "large_error": 0,
            "much_less_than_zero": 0,
            "no2_outlier": 0,
            "non_uni_strat": 0,
        }

        # Members from add_trop_data_to_gridsquare
        self.gno2 = None
        self.gstrat = None
        self.gcldp = None
        self.cntloop = None

        self.cloud_slice_count = 0

    def add_trop_data_to_gridsquare(self, trop_data):

        # These are reinitialised with every new product
        self.gno2 = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gstrat = [[[] for n in range(self.ydim)] for m in range(self.xdim)]
        self.gcldp = [[[] for n in range(self.ydim)] for m in range(self.xdim)]

        self.cntloop = [[0 for n in range(self.ydim)] for m in range(self.xdim)]

        for geo_tot_value, trop_lat, trop_lon, strat_no2_val, cloud_pressure \
             in zip(trop_data.geototvcd, trop_data.lats, trop_data.lons, trop_data.stratno2, trop_data.cldpres):
            # Skip over pixels where total column is less than stratospheric
            # column. This addresses positive bias in the cloud-sliced results
            # at low concentrations of UT NO2:
            if trop_data.geototvcd < strat_no2_val:
                continue

            # Find nearest gridsquare in output grid:
            p = np.argmin(abs(self.out_lon - trop_lon))
            q = np.argmin(abs(self.out_lat - trop_lat))

            # Convert columns from mol/m2 to molec/cm2:
            tvcdno2 = np.multiply(geo_tot_value, trop_data.no2sfac)
            tstrat = np.multiply(strat_no2_val, trop_data.no2sfac)

            # Get relevant data in each output grid square:
            # Something is going on here according to Pycharm; check in debugger.
            self.gno2[p, q].append(tvcdno2)
            self.gstrat[p, q].append(tstrat)
            self.gcldp[p, q].append(cloud_pressure)

            # Increment indices:
            self.cntloop[p, q] += 1

        # Save % valid points retained to print out average at end of routine:
        self.postfilt.append(100. * (trop_data.tcnt / trop_data.inicnt))
        self.file_count += 1

    def apply_cloud_slice(self, maxpnts, n_slices=40):
        # Estimate daily mean VMRs from the clustered data:
        for i in range(self.xdim):
            for j in range(self.ydim):

                tcolno2 = self.gno2[i, j]
                strat = self.gstrat[i, j]
                tcld = self.gcldp[i, j]

                # Only loop over grids with relevant data (identified as
                # vectors where the first entry is zero):
                if tcolno2[0] == 0:
                    continue

                # Skip if fewer than 10 points:
                if len(tcolno2) < 10:
                    self.loss_count["too_few_points"] += 1
                    continue

                # Convert from Pa to hPa for intput to the cloud-slicing algorithm:
                tcolno2 = np.multiply(tcolno2, 1e4)
                tcld = np.multiply(tcld, 1e-2)

                # Error check that the sizes of the arrays are equal:
                if (len(tcld) != len(tcolno2)):
                    print('Array sizes ne: cloud height and partial column')
                    raise UnequalColumnException

                # Skip scenes with non-uniform stratosphere using the
                # same threshold as is used for GEOS-Chem:
                if (np.std(strat) / np.mean(strat)) > 0.02:
                    self.loss_count["non_uni_strat"] += 1
                    continue

                # Get number of points:
                npnts = len(tcld)
                if npnts > maxpnts:
                    maxpnts = npnts
                    print(maxpnts, flush=True)

                # Use cloud_slice_ut_no2 function to get NO2 mixing
                # ratio from cloud-slicing:
                # TODO: Confirm with E that we skip here if npnts is between 10 and 20
                if ((npnts >= 20) & (npnts < 100)):
                    self.add_slice(i,j,tcld,tcolno2)

                elif (npnts >= 100):
                    # Define number of iterations:
                    stride = round(npnts / n_slices)
                    nloop = list(range(stride))
                    for w in nloop:
                        subset_t_col_no2 = tcolno2[w::stride]
                        subset_t_cld = tcld[w::stride]
                        self.add_slice(i, j, subset_t_cld, subset_t_col_no2)


    def add_slice(self, i, j, t_cld, t_col_no2):
        """Applies and adds a cloud slice from the given data"""
        # TODO: Check with E that the chunk of unreachable code at the end of cldslice is supposed to be there
        utmrno2, utmrno2err, stage_reached, mean_cld_pres = cldslice(t_col_no2, t_cld)

        # Calculate weights:
        gaus_wgt = np.exp((-(mean_cld_pres - 315) ** 2) / (2 * 135 ** 2))
        # Skip if approach didn't work (i.e. cloud-sliced UT NO2 is NaN):
        # Drop out after the reason for data loss is added to loss_count.
        if np.isnan(utmrno2) or np.isnan(utmrno2err):
            # Cloud-slicing led to nan, but due to rma regression rather
            # than data filtering (these are rare):
            if (stage_reached == 0):
                print("Cloud-sliced NO2 NAN for pixel i:{} j:{}".format(i, j))
                return
            self.loss_count[CLOUD_SLICE_ERROR_ENUM[stage_reached]] += 1
            # print("Cloud-slice exception {} in pixel i:{} j:{}".format(
            #    CLOUD_SLICE_ERROR_ENUM[stage_reached], i, j))
        else:
            self.gno2vmr[i, j] += np.multiply(utmrno2, gaus_wgt)
            self.gerr[i, j] += gaus_wgt
            self.gcnt[i, j] += 1
            self.cloud_slice_count += 1

    def calc_seasonal_means(self):
        self.mean_gno2vmr = np.divide(self.gno2vmr, self.gerr, where=self.gcnt != 0)
        self.mean_gerr = np.divide(self.gerr, self.gcnt, where=self.gcnt != 0)
        self.mean_gno2vmr[self.gcnt == 0] = np.nan
        self.mean_gerr[self.gcnt == 0] = np.nan
        self.gcnt[self.gcnt == 0] = np.nan   # Watch out for this rewriting of gcnt in the future

    def print_report(self):
        print('Max no. of data points in a gridsquare: ', np.amax(self.gcnt), flush=True)
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

    def save_to_netcdf(self, out_file):
        ncout = Dataset(out_file, mode='w',format="NETCDF4")

        ncout.createDimension('lat', self.ydim)
        ncout.createDimension('lon', self.xdim)

        # create longitude axis:
        lon = ncout.createVariable('lon', np.float32, ('lon',))
        lon.units = 'degrees_east'
        lon.long_name = 'longiitude'
        lon[:] = self.out_lon

        # Create latitude axis:
        lat = ncout.createVariable('lat', np.float32, ('lat',))
        lat.units = 'degrees_north'
        lat.long_name = 'latitude'
        lat[:] = self.out_lat

        # Save data: UT NO2 VMR (gno2vmr), UT NO2 error (gerr), No. of data points (gcnt)
        utno2 = ncout.createVariable('utno2', np.float32, ('lon', 'lat'))
        utno2.units = 'pptv'
        utno2.long_name = 'NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2[:] = self.mean_gno2vmr

        utno2err = ncout.createVariable('utno2err', np.float32, ('lon', 'lat'))
        utno2err.units = 'pptv'
        utno2err.long_name = 'Standard error of the NO2 mixing ratios in the UT (180-450 hPa) obtained using cloud-slicing'
        utno2err[:] = self.mean_gerr

        nobs = ncout.createVariable('nobs', np.float32, ('lon', 'lat'))
        nobs.units = 'unitless'
        nobs.long_name = 'Number of observations in each gridsquare used to obtain cloud-sliced UT NO2 mixing ratios'
        nobs[:] = self.gcnt

        ncout.close()

    def plot_data(self):
        # Plot the data:
        m = Basemap(resolution='l', projection='merc',
                    lat_0=0, lon_0=0, llcrnrlon=-180,
                    llcrnrlat=-75, urcrnrlon=180, urcrnrlat=80)
        X, Y = np.meshgrid(self.out_lon, self.out_lat, indexing='ij')
        xi, yi = m(X, Y)
        plt.subplot(1, 3, 1)
        cs = m.pcolor(xi, yi, np.squeeze(self.mean_gno2vmr), vmin=0, vmax=100, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('NO2 VMRs')

        plt.subplot(1, 3, 2)
        cs = m.pcolor(xi, yi, np.squeeze(self.mean_gerr), vmin=0, vmax=20, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('NO2 error')

        plt.subplot(1, 3, 3)
        cs = m.pcolor(xi, yi, np.squeeze(self.gcnt), vmin=0., vmax=30, cmap='jet')
        m.drawparallels(np.arange(-80., 81., 45.), labels=[1, 0, 0, 0], fontsize=8)
        m.drawmeridians(np.arange(-180., 181., 45.), labels=[0, 0, 0, 1], fontsize=8)
        m.drawcoastlines()
        m.drawcountries()
        cbar = m.colorbar(cs, location='bottom', pad="10%")
        plt.title('Number of points')

        plt.show()


class TropomiData:

    #This time, I'm not putting the read in the init - this is _just_ for member initialisation
    def __init__(self, file_path, pmax, pmin):

        self.file_name = path.basename(file_path)
        self.date = get_date(self.file_name)
        self.pmax = pmax
        self.pmin = pmin

        # Members straight from trop body
        self.no2sfac = None
        self.qasfac = None
        self.fillval = None
        self.tlons = None
        self.tlats = None
        self.tscdno2 = None
        self.stratno2_og = None
        self.tscdno2err = None
        self.stratno2err = None
        self.tstratamf = None
        self.qaval = None
        self.aai = None
        self.sza = None
        self.vza = None

        # Members from geometric column
        self.tamf_geo = None
        self.tgeotropvcd = None

        # Members from bias correction
        self.tstratno2 = None
        self.tgeototvcd = None
        self.ttropvcd_geo_err = None  # This one doesn't seem to be used

        # Members from filtering
        self.inicnt = None
        self.tcnt = None

        # Data members from trimming + reshaping
        self.lons = None
        self.lats = None
        self.cldpres = None # Should this be here?
        self.stratno2 = None
        self.amf_geo = None
        self.geototvcd = None

        # Do initialisation
        self.read_trop_file(file_path)

    def read_trop_file(self, file_path):
        fh = Dataset(file_path, mode='r')

        # no2sfac, qasfac, and fillval only need to be read in once, so could
        # just be read in on the first iteration.
        # NO2 conversion factor:
        no2sfac = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric' \
                                                 '_column'].multiplication_factor_to_convert_to_molecules_percm2

        # QA flag scale factor:
        qasfac = fh.groups['PRODUCT'].variables['qa_value'].scale_factor

        # NO2 fill/missing value:
        fillval = fh.groups['PRODUCT'].variables['nitrogendioxide_tropospheric' \
                                                 '_column']._FillValue
        # Extract data of interest:

        # Geolocation data:
        glons = fh.groups['PRODUCT'].variables['longitude'][:]
        tlons = glons[0, :, :]
        glats = fh.groups['PRODUCT'].variables['latitude'][:]
        tlats = glats[0, :, :]

        # Column data:
        # (1) Total slant column:
        gscdno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                      variables['nitrogendioxide_slant_column_density'][:]
        tscdno2 = gscdno2.data[0, :, :]
        # (2) Stratospheric vertical column:
        gstratno2 = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['nitrogendioxide_stratospheric_column'][:]
        stratno2_og = gstratno2.data[0, :, :]

        # Precisions/Uncertainties of column data:
        # (1) Total slant column uncertainty:
        gscdno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                         variables['nitrogendioxide_slant_column_density_precision'][:]
        tscdno2err = gscdno2err.data[0, :, :]
        # (2) Stratospheric vertical column uncertainty:
        stratno2err = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                          variables['nitrogendioxide_stratospheric_column_precision'][0, :, :]
        # Stratospheric AMF:
        gstratamf = fh.groups['PRODUCT']['SUPPORT_DATA']['DETAILED_RESULTS']. \
                        variables['air_mass_factor_stratosphere'][:]
        tstratamf = gstratamf.data[0, :, :]

        # QA value:
        qaval = fh.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Aerosol absorbing index:
        taai = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                   variables['aerosol_index_354_388'][:]
        aai = taai.data[0, :, :]
        aai = np.where(aai > 1e30, np.nan, aai)

        # Solar zenith angle (degrees):
        tsza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['solar_zenith_angle'][:]
        sza = np.ma.getdata(tsza[0, :, :])

        # Viewing zenith angle (degrees):
        tvza = fh.groups['PRODUCT']['SUPPORT_DATA']['GEOLOCATIONS']. \
                   variables['viewing_zenith_angle'][:]
        vza = np.ma.getdata(tvza[0, :, :])

        # Setting members
        self.no2sfac = no2sfac
        self.qasfac = qasfac
        self.fillval = fillval
        self.tlons = tlons
        self.tlats = tlats
        self.tscdno2 = tscdno2
        self.stratno2_og = stratno2_og
        self.tscdno2err = tscdno2err
        self.stratno2err = stratno2err
        self.tstratamf = tstratamf
        self.qaval = qaval
        self.aai = aai
        self.sza = sza
        self.vza = vza

    def calc_geo_column(self):
        # Calculate the geometric AMF:
        tamf_geo = np.add((np.reciprocal(np.cos(np.deg2rad(self.sza)))),
                          (np.reciprocal(np.cos(np.deg2rad(self.vza)))))

        # Get VCD under cloud conditions. This is done as the current
        # tropospheric NO2 VCD product includes influence from the prior
        # below clouds:
        # Calculate the stratospheric slant columns:
        tscdstrat = np.multiply(self.stratno2_og, self.tstratamf)
        # Calculate the tropospheric slant columns:
        ttropscd = np.subtract(self.tscdno2, tscdstrat)
        # Calculate the tropospheric vertical column using the geometric AMF:
        tgeotropvcd = np.divide(ttropscd, tamf_geo)

        # Setting members
        self.tamf_geo = tamf_geo
        self.tgeotropvcd = tgeotropvcd

    def apply_bias_correction(self):
        # Bias correct stratosphere based on comparison of TROPOMI to Pandora Mauna Loa:
        tstratno2 = np.where(self.stratno2_og != self.fillval,
                             ((self.stratno2_og - (7.3e14 / self.no2sfac)) / 0.82),
                             np.nan)

        # Bias correct troposphere based on comparison of TROPOMI to Pandora Izana:
        tgeotropvcd = np.where(self.tgeotropvcd != self.fillval,
                               self.tgeotropvcd / 2.,
                               np.nan)

        # Get the total column as the sum of the bias-corrected components:
        tgeototvcd = np.add(tgeotropvcd, tstratno2)

        # Calculate updated stratospheric NO2 error after bias correcting.
        # Determine by scaling it by the relative change in stratospheric vertical
        # colum NO2 after applying a bias correction:
        tstratno2err = np.where(self.stratno2err != self.fillval,
                                np.multiply(self.stratno2err, np.divide(self.tstratno2, self.stratno2_og)),
                                np.nan)

        # Calculate error by adding in quadrature individual
        # contributions:
        ttotvcd_geo_err = np.sqrt(np.add(np.square(tstratno2err),
                                         np.square(self.tscdno2err)))
        # Estimate the tropospheric NO2 error as the total error
        # weighted by the relative contribution of the troposphere
        # to the total column, as components that contribute to the
        # error are the same:
        ttropvcd_geo_err = np.multiply(ttotvcd_geo_err,
                                       (np.divide(tgeotropvcd, tgeototvcd)))

        self.tstratno2 = tstratno2
        self.tgeototvcd = tgeototvcd
        self.tgeotropvcd = tgeotropvcd  # Filter applied to member defined in geo_column
        self.ttropvcd_geo_err = ttropvcd_geo_err

    def cloud_filter_and_preprocess(self, cloud_data, cldthld):

        # Do date check
        if self.date != cloud_data.date:
            print('NO2 file: {}, Cloud file: {}'.format(self.date, cloud_data.date), flush=True)
            print('EXITING: Files are not for the same date!', flush=True)
            raise CloudFileDateMismatch

        # Check that data shapes are equal:
        if cloud_data.cldfrac.shape != self.sza.shape:
            print('Cloud product and NO2 indices ne!', flush=True)
            print(cloud_data.cldfrac.shape, self.sza.shape, flush=True)
            print('Skipping this swath', flush=True)
            raise CloudFileShapeMismatch

        tgeototvcd = self.tgeototvcd
        # Filter to only include very cloudy scenes at high altitude
        # No. of valid total NO2 column points before apply filtering:
        self.inicnt = np.count_nonzero(~np.isnan(tgeototvcd))

        # Filter out scenes with cloud fractions (and hence cloud pressures) that are nan:
        tgeototvcd = np.where(np.isnan(cloud_data.cldfrac), np.nan, tgeototvcd)

        # Filter out scenes with cloud fraction < cloud threshold:
        tgeototvcd = np.where(cloud_data.cldfrac < cldthld, np.nan, tgeototvcd)

        # Filter out scenes with cloud heights outside the UT range of interest (180-450 hPa):
        tgeototvcd = np.where(cloud_data.tcldpres > self.pmax * 1e2, np.nan, tgeototvcd)
        tgeototvcd = np.where(cloud_data.tcldpres < self.pmin * 1e2, np.nan, tgeototvcd)

        # Filter out low quality data (0.45 threshold suggested TROPOMI NO2
        # PI Henk Eskes from the KNMI:
        tgeototvcd = np.where(self.qaval < 0.45, np.nan, tgeototvcd)

        # Filter out scenes with AAI > 1 (could be misclassified as clouds)
        tgeototvcd = np.where(self.aai > 1., np.nan, tgeototvcd)
        self.tgeototvcd = tgeototvcd

        # No. of points retained after filtering:
        self.tcnt = np.count_nonzero(~np.isnan(self.tgeototvcd))

        # Trim the data to include only those relevant:
        # This also reshapes the data from a 2D to a 1D array:
        self.lons = self.tlons[~np.isnan(tgeototvcd)]
        self.lats = self.tlats[~np.isnan(tgeototvcd)]
        self.stratno2 = self.tstratno2[~np.isnan(tgeototvcd)]
        self.amf_geo = self.tamf_geo[~np.isnan(tgeototvcd)]
        self.geototvcd = self.tgeototvcd[~np.isnan(tgeototvcd)]

        self.cldpres = cloud_data.tcldpres[~np.isnan(tgeototvcd)]


class CloudData:

    def __init__(self, file_path, data_type, fillval = None):
        # Set from filename
        self.file_name = path.basename(file_path)
        self.date = get_date(self.file_name)  # Assuming for now that ocra and S5P have the same timestamping

        # Set from file contents
        self.cldfrac = None
        self.tcldpres = None
        self.tsnow = None

        if data_type == 'fresco':
            self.read_fresco_file(file_path)
        elif data_type == 'dl_ocra':
            self.read_ocra_file(file_path, fillval)

    def read_fresco_file(self, file_path):
        fh = Dataset(file_path)
        # Cloud fraction:
        tcldfrac = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_fraction_crb'][:]
        self.cldfrac = tcldfrac.data[0, :, :]
        # Cloud top pressure:
        gcldpres = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                       variables['cloud_pressure_crb'][:]
        self.tcldpres = gcldpres[0, :, :]

        # Get scene and surface pressure to diagnose clouds misclassified as snow/ice:
        # Apparent scene pressure:
        gscenep = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                      variables['apparent_scene_pressure'][:]
        tscenep = gscenep.data[0, :, :]
        # Surface pressure:
        gsurfp = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                     variables['surface_pressure'][:]
        tsurfp = gsurfp.data[0, :, :]

        # Snow/ice flag:
        gsnow = fh.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                    variables['snow_ice_flag'][:]
        tsnow = gsnow.data[0, :, :]
        # Convert all valid snow/ice free flag values (0,255) to 0.
        # Ocean:
        tsnow = np.where(tsnow == 255, 0, tsnow)
        # Coastlines (listed as potential "suspect" in the ATBD document p. 67):
        tsnow = np.where(tsnow == 252, 0, tsnow)
        # Less then 1% snow/ice cover:
        tsnow = np.where(tsnow < 1, 0, tsnow)
        # Snow/ice misclassified as clouds:
        self.tsnow = np.where(((tsnow > 80) & (tsnow < 104) & (tscenep > (0.98 * tsurfp))),
                         0, tsnow)
        # Set clouds over snow/ice scenes to nan:
        self.cldfrac = np.where(self.tsnow != 0, np.nan, self.cldfrac)
        self.tcldpres = np.where(self.tsnow != 0, np.nan, self.tcldpres)
        fh.close()

    def read_ocra_file(self, file_path, fillval):

        fd = Dataset(file_path)
        # Cloud fraction:
        tcldfrac = fd.groups['PRODUCT'].variables['cloud_fraction'][:]
        cldfrac = tcldfrac.data[0, :, :]
        # Cloud top height (m):
        gcldhgt = fd.groups['PRODUCT'].variables['cloud_top_height'][:]
        tcldhgt = np.ma.getdata(gcldhgt[0, :, :])

        # Define pressure array of zeros:
        tcldpres = np.zeros(tcldhgt.shape)

        # Calculate pressure assuming dry atmosphere using external
        # conversion code (convert_height_to_press.py). There's a cloud
        # top pressure entry in the data file, but this is obtained using
        # ECMWF pressure and might have errors. Diego (DLR cloud product PI
        # recommended I use cloud altitude rather than pressure data):
        hgtind = np.where((tcldhgt != fillval))
        tcldpres[hgtind] = alt2pres(tcldhgt[hgtind])

        # QA value:
        cldqa = fd.groups['PRODUCT'].variables['qa_value'][0, :, :]

        # Snow/ice flag (combined NISE and climatology, so misclassification
        # issues in FRESCO cloud product addressed):
        gsnow = fd.groups['PRODUCT']['SUPPORT_DATA']['INPUT_DATA']. \
                    variables['snow_ice_flag'][:]
        self.tsnow = gsnow.data[0, :, :]

        # Set clouds over snow/ice scenes to nan:
        cldfrac = np.where(self.tsnow != 0, np.nan, cldfrac)
        tcldpres = np.where(self.tsnow != 0, np.nan, tcldpres)

        # Set poor quality cloud data to nan:
        self.cldfrac = np.where(cldqa < 0.5, np.nan, cldfrac)
        self.tcldpres = np.where(cldqa < 0.5, np.nan, tcldpres)

        # Close DLR CLOUD file:
        fd.close()


#   TODO: Move these into a seperate file for reuse maybe
def get_tropomi_file_list(trop_dir, date_range):
    out = []
    for date in date_range:
        out += (get_tropomi_files_on_day(trop_dir, date))
    return sorted(out)


def get_ocra_file_list(ocra_dir, date_range):
    out = []
    for date in date_range:
        out += (get_ocra_files_on_day(ocra_dir, date))


def get_tropomi_files_on_day(tomidir, date):
    # Converts the python date object to a set string representation of time
    # In this case, zero-padded year, month and a datestamp of the Sentinel format
    # See https://docs.python.org/3/library/datetime.html#strftime-and-strptime-format-codes
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    tomi_glob_string = path.join(tomidir, 'NO2_OFFL', year, month,'S5P_OFFL_L2__NO2____'+ datestamp + '*')
    tomi_files_on_day = glob.glob(tomi_glob_string)
    print('Found {} tropomi files for {}: '.format(len(tomi_files_on_day), date))
    tomi_files_on_day = sorted(tomi_files_on_day)
    return tomi_files_on_day


def get_ocra_files_on_day(tomidir,date):
    # Get string of day:
    year = date.strftime(r"%Y")
    month = date.strftime(r"%m")
    datestamp = date.strftime(r"%Y%m%dT")
    cld_glob_string = path.join(tomidir, "CLOUD_OFFL", year, month,
                                   'S5P_OFFL_L2__CLOUD__' + datestamp + '*')
    cldfile = glob.glob(cld_glob_string)[0]
    # Order the files:
    cldfile = sorted(cldfile)
    return cldfile


def get_date(file_name):
    # A regular expression that gets Sentinel datestamps out of filenames
    # See https://regex101.com/r/QNG11l/1
    date_regex = r"\d{8}T\d{6}"
    date_string = re.findall(date_regex, file_name)[0]
    # A line for converting Sentinel string reps to datetime
    return dt.datetime.strptime(date_string, r"%Y%m%dT%H%M%S")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("trop_dir")
    parser.add_argument("out_dir")
    parser.add_argument("--season", default='jja', help="Can be jja, son, djf, mam")
    parser.add_argument("--grid_res", default='1x1', help="Can be 1x1, 2x25, 4x5")
    parser.add_argument("--cloud_product", default = "fresco", help="can be fresco or dl-ocra")
    parser.add_argument("--cloud_threshold", default = "07", help="can be 07, 08, 09, 10")
    parser.add_argument("--ocra_fill", help="Fill value for Ocra data")
    args = parser.parse_args()

    if args.season == "jja":
        start_date = dt.datetime(year=2019, month=6, day=1)
        end_date = dt.datetime(year=2019, month=8, day=31)
        yrrange = '2019'
    elif args.season == "son":
        start_date = dt.datetime(year=2019, month=9, day=1)
        end_date = dt.datetime(year=2019, month=11, day=30)
        yrrange = '2019'
    elif args.season == "djf":
        start_date = dt.datetime(year=2019, month=12, day=1)
        end_date = dt.datetime(year=2020, month=2, day=29)  # Beware the leap year here
        yrrange = '2019-2020'
    elif args.season == "mam":
        start_date = dt.datetime(year=2020, month=3, day=1)
        end_date = dt.datetime(year=2020, month=6, day=29)
        yrrange = '2020'
    else:
        print("Invalid season; can be jja, son, djf, mam")
        sys.exit(1)

    if args.grid_res == '1x1':
        dellat, dellon = 1, 1
    elif args.grid_res == '2x25':
        dellat, dellon = 2, 2.5
    elif args.grid_res == '4x5':
        dellat, dellon = 4, 5
    else:
        print("Invalid grid; values can be 1x1, 2x25, 4x5")
        sys.exit(1)

    date_range = rr.rrule(rr.DAILY, dtstart=start_date, until=end_date)

    trop_files = get_tropomi_file_list(args.trop_dir, date_range)
    if args.cloud_product == "fresco":
        cloud_files = trop_files
    elif args.cloud_product == "dl-ocra":
        cloud_files = get_ocra_file_list(args.trop_dir, date_range)
    else:
        print("Invalid cloud product; can be fresco or dl-ocra")
        sys.exit(1)

    grid_aggregator = GridAggregator(dellat, dellon)

    # Define pressure ranges:
    PMIN=180
    PMAX=450

    for trop_file, cloud_file in zip(trop_files, cloud_files):
        # Looks like this is currently set to process a string of files as
        # opposed to a single file:
        trop_data = TropomiData(trop_file,PMAX,PMIN)
        cloud_data = CloudData(cloud_file, data_type=args.cloud_product, fillval=args.ocra_fill)

        trop_data.calc_geo_column()
        trop_data.apply_bias_correction()
        trop_data.cloud_filter_and_preprocess(cloud_data, args.cloud_threshold)

        grid_aggregator.add_trop_data_to_gridsquare(trop_data)
        grid_aggregator.apply_cloud_slice(args.max_points)


    out_file = path.join(args.out_dir, 'tropomi-ut-no2-'+args.cloud_product
                         + '-' + args.cloud_threshold
                         + '-' + args.grid_res
                         + '-' + args.season
                         + '-' + yrrange+'-v7.nc')
    grid_aggregator.print_report()
    grid_aggregator.save_to_netcdf(args.out_dir)
    grid_aggregator.plot_data()



