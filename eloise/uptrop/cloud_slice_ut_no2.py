#!/usr/bin/python

# Import relevant packages:
import glob
import sys
import os
import numpy as np
from bootstrap import rma

from constants import AVOGADRO as na
from constants import G as g
from constants import MW_AIR as mmair

CLOUD_SLICE_ERROR_ENUM = {
    1: "too_few_points",
    2: "low_cloud_height_range",
    3: "low_cloud_height_std",
    4: "large_error",
    5: "much_less_than_zero",
    6: "no2_outlier",
    7: "non_uni_strat"
}

def cldslice(pcolno2,cldtophgt):

    """ 
    Compute upper troposphere NO2 using partial columns above
    cloudy scenes. 

    Determine NO2 mixing ratio by regressing NO2 partial columns
    against cloud-top heights over cloudy scenes.

    INPUT: vectors of partial columns in molec/m2 and corresponding 
           cloud top heights in hPa.

    OUTPUT: NO2 volumetric mixing ratio, corresponding estimated error on the
            cloud-sliced NO2 value, a number to identify which filtering
            criteria led to loss of data in the case that the cloud-sliced
            NO2 value ia nan, and the mean cloud pressure of data retained
            after 10th and 90th percentile filtering.
    """

    # Initialize:
    utmrno2=0.0
    utmrno2err=0.0
    error_state=0

    # Define factor to convert slope of NO2 partial column vs pressure
    # to VMR:
    den2mr=np.divide((np.multiply(g,mmair)),na)

    # Get 10th and 90th percentiles of data population:
    p10=np.percentile(pcolno2,10)
    p90=np.percentile(pcolno2,90)

    # Remove outliers determined as falling outside the 10th and 90th 
    # percentile range. Not include this or instead using 5th and 95th leads
    # to overestimate in cloud-sliced UT NO2 compared to the "truth":
    sind=np.where((pcolno2>p10)&(pcolno2<p90))[0]
    # Trim the data to remove ouliers:
    pcolno2=pcolno2[sind]
    cldtophgt=cldtophgt[sind]

    # Cloud pressure mean:
    mean_cld_pres=np.mean(cldtophgt)

    # Get number of points in vector:
    npoints=len(cldtophgt)

    # Only consider data with more than 5 points for reasonably
    # robust statistics. This step is added to account for data loss
    # removing outliers:
    if npoints<=10:
        error_state=1
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Get cloud top height standard deviation:
    stdcld=np.std(cldtophgt)
    # Get cloud top height range:
    diffcld=np.nanmax(cldtophgt)-np.nanmin(cldtophgt)

    # Only consider scenes with a dynamic range of clouds:
    # (i) Cloud range:
    if diffcld<=140:
        error_state=2
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # (ii) Cloud standard deviation:
    if stdcld<=30:
        error_state=3
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)


    # Get regression statistics:
    # Partial NO2 column (molec/m2) vs cloud top height (hPa):
    # 300 iterations of regression chosen to compromise between
    # statistics and computational efficiency:
    result=rma(cldtophgt*1e2,pcolno2,len(pcolno2),300)

    # Remove data with relative error > 100%:
    if np.absolute(np.divide(result[2], result[0]))>1.0:
        error_state=4
        utmrno2=np.nan
        utmrno2err=np.nan
        return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Account for negative values:
    # Set points with sum of slope and error less than zero to nan.
    # This is to account for noise in the data that hover near zero.
    if result[0]<0 and (not np.isnan(utmrno2)):
        if (np.add(result[0],result[2])<0):
            error_state=5
            utmrno2=np.nan
            utmrno2err=np.nan
            return (utmrno2, utmrno2err, error_state, mean_cld_pres)

        slope=result[0]
        #slope=np.multiply(slope,sf)
        slope_err=result[2]
        #slope_err=np.multiply(slope_err,sf)
        # Convert slope to mol/mol:
        utmrno2=np.multiply(slope,den2mr)
        # Convert error to mol/mol:
        utmrno2err=np.multiply(slope_err,den2mr)
        # Convert UT NO2 from mol/mol to ppt:
        utmrno2=np.multiply(utmrno2,1e+12)
        # Convert UT NO2 error from mol/mol to ppt
        utmrno2err=np.multiply(utmrno2err,1e+12)

        # Finally, remove outliers in the cloud-sliced NO2
        # 200 pptv threshold is chosen, as far from likely range.
        # Scale factor applied to TROPOMI UT NO2 to account for
        # positive bias in free tropospheric NO2:
        if utmrno2>200:
            error_state=6
            utmrno2=np.nan
            utmrno2err=np.nan
            return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    # Output goes here:
    return (utmrno2, utmrno2err, error_state, mean_cld_pres)

    
