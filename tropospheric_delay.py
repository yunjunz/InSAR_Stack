#! /usr/bin/env python
# Author: Heresh Fattahi

import pyaps as pa
import numpy as np
import sys
import h5py
import os, imp, sys, glob
import argparse
import configparser
import  datetime
import time



print '------------------------------------------------'
print 'You are using PyAPS from %s'%pa.__file__
print '------------------------------------------------'

def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='extracts the overlap geometry between master bursts')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True,
            help='An h5 timeseries file that includes a list of dates')
    parser.add_argument('-d', '--dem', type=str, dest='dem', default=None,
            help='dem')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True,
            help='output h5 file that contains timeseries of tropospheric delay')
    parser.add_argument('-t', '--time', type=str, dest='time', required=True,
            help='SAR acquisition time (0 6 12 18 for ERA-I)')
    parser.add_argument('-l', '--look_angle', type=float, dest='lookAngle', required=True,
            help='look angle to project the zenith delay to radar line of sight')
    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps

def getDates(inps):

    h5=h5py.File(inps, 'r')
    dsDates = h5.get('dateList')
    dateList1 = list(dsDates[:])
    dateList = []
    for i in range(len(dateList1)):
        dateList.append(dateList1[i].decode("utf-8")[0:10].replace('-',''))
    return dateList

def download_delay(inps):
    
    dirName = os.path.dirname(os.path.abspath(inps.input))
    ECMWF_dir = os.path.join(dirName, 'ECMWF')
    if not os.path.exists(ECMWF_dir):
       os.makedirs(ECMWF_dir)
    
    dateList = getDates(inps.input)    
    print(dateList)
    pa.ECMWFdload(dateList, inps.time, ECMWF_dir ) 
    ecmwf_files={}
    for d in dateList:
        emcwfName = os.path.join(ECMWF_dir, 'ERA-Int_' + d + '_' + inps.time + '.grb')
        ecmwf_files[d]=emcwfName
    
    return ecmwf_files

def compute_delay(inps, ecmwf_file):

    dateList = getDates(inps.input)    
    aps = pa.PyAPS_rdr(ecmwf_file[dateList[0]],inps.dem, grib='ECMWF', demfmt='HGT')
    h5out = h5py.File(inps.output,'w')
    ds = h5out.create_dataset('ecmwf', shape=(len(dateList),aps.ny,aps.nx),dtype=np.float32)
    numDates = len(dateList)

    for i in range(numDates):
        d = dateList[i]
        print('computing delay for: ', d)
        aps = pa.PyAPS_rdr(ecmwf_file[d],inps.dem, grib='ECMWF', demfmt='HGT')    
        delay = np.zeros((aps.ny, aps.nx))
        aps.getdelay(delay, inc=inps.lookAngle)
        ds[i,:,:] = delay

    h5out.close()

def main(iargs=None):
    inps = cmdLineParse(iargs)
    ecmwf_files = download_delay(inps)
    print(ecmwf_files)
    compute_delay(inps, ecmwf_files)

if __name__ == '__main__' :
  ''' 
  invert a network of the pair's mis-registrations to
  estimate the mis-registrations wrt the Master date.
  '''

  main()


'''
print 'Testing Download Methods'
print 'Testing ECMWF Downloads'
date1 = '20071014'
date2 = '20071129'

pa.ECMWFdload([date1, date2],'6','./ECMWF/')

print 'Testing ECMWF in Radar geometry, with a FLT dem'
dem = '../geom_master/z.14alks_4rlks_masked.rdr'
lat = '../geom_master/lat.14alks_4rlks.rdr'
lon = '../geom_master/lon.14alks_4rlks.rdr'

aps1 = pa.PyAPS_rdr('ECMWF/ERA-Int_' + date1 + '_6.grb',dem, grib='ECMWF', demfmt='HGT')
aps2 = pa.PyAPS_rdr('ECMWF/ERA-Int_' + date2 + '_6.grb',dem, grib='ECMWF', demfmt='HGT')

print 'With Lat Lon files'
phs1 = np.zeros((aps1.ny,aps1.nx))
phs2 = np.zeros((aps2.ny,aps2.nx))

aps1.getdelay(phs1, inc=38.85, wvl=0.2360571)
aps2.getdelay(phs2, inc=38.85, wvl=0.2360571)

#aps1.getgeodelay(phs1, inc=38.85, wvl=0.2360571, lat=lat, lon=lon)
#aps2.getgeodelay(phs2, inc=38.85, wvl=0.2360571, lat=lat, lon=lon)

LLphs = phs2-phs1

#plt.imshow(LLphs)
#plt.colorbar()
#plt.show()

h5=h5py.File('tropoPhase.h5','w')
ds=h5.create_dataset('aps',data=LLphs, shape=LLphs.shape,dtype=np.float32)
h5.close()
'''
