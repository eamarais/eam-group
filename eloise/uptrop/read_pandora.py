#!/usr/bin/python

# Import relevant packages:
import glob
import sys
import os

import numpy as np
import pandas as pd

def readpandora(filename):

    """
    """
    
    # Initialize:
    nvals=250000
    yyyy=np.zeros(nvals)
    mon=np.zeros(nvals)
    day=np.zeros(nvals)
    utc_hh=np.zeros(nvals)
    mins=np.zeros(nvals)
    jday=np.zeros(nvals)
    sza=np.zeros(nvals)
    no2=np.zeros(nvals)
    no2err=np.zeros(nvals)
    qaflag=np.zeros(nvals)
    fitflag=np.zeros(nvals)

    cnt=0    # line counts

    #f = open(filename, "rb")

    f=open(filename, 'r', encoding='utf8', errors='ignore')
    #sys.exit()
    Lines=f.readlines()

    #print(Lines[14:16])

    # Get site latitude:
    txt=Lines[14]
    splttxt=txt.split()
    lat=float(splttxt[3])

    # Get site longitude:
    txt=Lines[15]
    splttxt=txt.split()
    lon=float(splttxt[3])

    print('Pandora lon and lat: ',lon,lat)

    # Initialize:
    First=0

    # Get data:
    # loop over lines:
    for w, Lines in enumerate(Lines):

        # Rename line and decode:
        txt=Lines

        # Find relevant data. This approach seems convoluted, but it's
        # to account for different data entry locations and descriptions in
        # the total and tropospheric column data files.
        # Get entries for relevant data: # Julian Day:
        if (txt[10:42]=='UT date and time for center of m' ):
            dateind=int(txt[6:8])-1
        # Julian Day:
        if (txt[10:42]=='Fractional days since 1-Jan-2000' ):
            jdayind=int(txt[6:8])-1
        # SZA:
        if (txt[10:42]=='Solar zenith angle for center of' ):
            szaind=int(txt[6:8])-1
        # NO2 column:
        # (a) Total:
        if (txt[10:42]=='Nitrogen dioxide total vertical ' ):
            no2ind=int(txt[6:8])-1
        # (b) Tropospheric:
        if (txt[10:42]==' Nitrogen dioxide tropospheric v' ):
            no2ind=int(txt[7:9])-1
        # NO2 column error:
        # (a) Total:
        if (txt[10:51]=='Uncertainty of nitrogen dioxide total ver' ):
            errind=int(txt[6:8])-1
        # (b) Tropospheric:
        if (txt[10:51]==' Uncertainty of nitrogen dioxide troposph' ): 
            errind=int(txt[7:9])-1
        # Data quality flag:
        if (txt[10:42]==' L2 data quality flag for nitrog' ):
            qaflagind=int(txt[7:9])-1
        # Level 2 fit flag:
        # There are two entries (min and max) of this in the tropospheric
        # column file. The minimum is being used here.
        if (txt[10:40]==' Level 2 Fit data quality flag' ):
            fitflagind=int(txt[7:9])-1     

        # Find divide between headers and data:
        if (txt[:10]=='----------'): 
            First=First+1
            continue

        # Skip header lines:
        if First<2: continue

        # Split line:
        splttxt=txt.split()
        # Extract relevant data and recast data type:
        tdate,tjday,tsza,tno2,tno2err,tqaflag,tfitflag=str(splttxt[dateind]),\
                     float(splttxt[jdayind]),float(splttxt[szaind]),\
                     float(splttxt[no2ind]),float(splttxt[errind]),\
                     float(splttxt[qaflagind]),float(splttxt[fitflagind])

        # Add to final data array:
        jday[cnt]=float(tjday)
        sza[cnt]=float(tsza)
        no2[cnt]=float(tno2)
        no2err[cnt]=float(tno2err)
        qaflag[cnt]=float(tqaflag)
        fitflag[cnt]=float(tfitflag)

        # Split date into YYYY MM DD and UTC HH MIN and add to final array:
        yyyy[cnt]=float(tdate[0:4])
        mon[cnt]=float(tdate[4:6])
        day[cnt]=float(tdate[6:8])
        utc_hh[cnt]=float(tdate[9:11])
        mins[cnt]=float(tdate[11:13])

        # Increment:
        cnt=cnt+1
    
    # Trim the arrays to remove trailing zeros:
    mon=mon[np.nonzero(yyyy)]
    day=day[np.nonzero(yyyy)]
    utc_hh=utc_hh[np.nonzero(yyyy)]
    mins=mins[np.nonzero(yyyy)]
    jday=jday[np.nonzero(yyyy)]
    sza=sza[np.nonzero(yyyy)]
    no2=no2[np.nonzero(yyyy)]
    no2err=no2err[np.nonzero(yyyy)]
    qaflag=qaflag[np.nonzero(yyyy)]
    fitflag=fitflag[np.nonzero(yyyy)]
    yyyy=yyyy[np.nonzero(yyyy)]

    # Convert arrays to data frames:
    df=pd.DataFrame({'year':yyyy, 'month':mon, 'day':day, 'hour_utc':utc_hh,\
                     'minute':mins, 'jday':jday, 'sza':sza, 'no2':no2,\
                     'no2err':no2err, 'qaflag':qaflag, 'fitflag':fitflag})

    # Convert latitude and longitude to dataframe:
    loc={'lat':lat, 'lon':lon} 

    # Output:
    return (loc, df)
