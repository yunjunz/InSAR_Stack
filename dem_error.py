#! /usr/bin/env python
# Author: Heresh Fattahi

import numpy as np
import sys
import h5py
import os, imp, sys, glob
import argparse
import configparser
import  datetime
import time


def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='extracts the overlap geometry between master bursts')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True,
            help='An h5 timeseries file that includes a list of dates')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True,
            help='output h5 file that contains timeseries of tropospheric delay')
    parser.add_argument('-b', '--baseline', type=str, dest='baseline', required=True,
            help='baseline list')
    parser.add_argument('-d', '--data', type=str, dest='data', required=True,
            help='dataset name')
    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps

def getBaselineList(inps, numDates):

    bPerp = np.zeros((numDates,1))
    i = 0
    for line in open(inps.baseline):
       c = line.split()
       if len(c) < 2 or line.startswith('%') or line.startswith('#'):
          next #ignore commented lines or those without variables
       else:
          bPerp[i,0]=str.replace(c[1],'\n','').strip()
          i+=1

    return bPerp

def designMatrix(inps):

    h5=h5py.File(inps.input, 'r')
    dsDates = h5.get('dateList')
    dateList = list(dsDates[:])
    numDates = len(dateList)
    bPerp = getBaselineList(inps, numDates)
    t = np.zeros((numDates,1))
    for i in range(numDates):
        d = dateList[i].decode("utf-8")
        d = datetime.datetime(*time.strptime(d,"%Y-%m-%d %H:%M:%S")[0:6]) #.strftime('%Y%m%d')
        d = d.year + (d.month-1)/12.0 + (d.day-1)/365.0
        t[i,0] = d

    A = np.hstack((.5*t**2, t, np.ones((numDates,1)), bPerp))

    h5.close()

    return A, bPerp

def designMatrix_v2(inps):

    h5=h5py.File(inps.input, 'r')
    dsDates = h5.get('dateList')
    dateList = list(dsDates[:])
    numDates = len(dateList)
    bPerp = getBaselineList(inps, numDates)
    t = np.zeros((numDates,1))
    for i in range(numDates):
        d = dateList[i].decode("utf-8")
        d = datetime.datetime(*time.strptime(d,"%Y-%m-%d %H:%M:%S")[0:6]) #.strftime('%Y%m%d')
        d = d.year + (d.month-1)/12.0 + (d.day-1)/365.0
        t[i,0] = d

    stepFunc = np.ones_like(bPerp)
    t_eq = 2007.0 + (11.-1.)/12. + (4-1)/365.0
    stepFunc[t<t_eq] = 0.0
    A = np.hstack((.5*t**2, t, np.ones((numDates,1)), stepFunc, bPerp))

    h5.close()

    return A, bPerp

def estimateDemErr(inps, A, Bp):
    h5=h5py.File(inps.input, 'r')

    if inps.data:
       k=inps.data
       data = h5.get(inps.data)
    else:
       k=list(h5.keys())[0]
       print(k)
       data = h5.get(k)

    Nz, Ny, Nx = data.shape
    Npar = A.shape[1]
    A1 = np.linalg.pinv(A)
    A1 = np.array(A1,np.float32)

    ##########

    h5out = h5py.File(inps.output,'w')
    ds = h5out.create_dataset(k, shape=(Nz,Ny,Nx),dtype=np.float32)
    ds_demErr = h5out.create_dataset('dem_error', shape=(Ny,Nx),dtype=np.float32)
    dateList = h5.get('dateList')
    dsDateList = h5out.create_dataset('dateList', data=dateList, dtype=dateList.dtype)

    dZ = np.zeros((1,Nx))
    for i in range(Ny):
        print(i, 'out of ',Ny)
        L = data[:,i,:]
        X = np.dot(A1, L)
        dZ[0,:]=X[-1,:]
        r = L-np.dot(Bp,dZ)
        ds[:,i,:] = r
        ds_demErr[i,:] = dZ
    h5out.close()
    h5.close()


def main(iargs=None):
    inps = cmdLineParse(iargs)
    A, Bp = designMatrix_v2(inps)    
    estimateDemErr(inps, A, Bp)

if __name__ == '__main__' :
  ''' 
  invert a network of the pair's mis-registrations to
  estimate the mis-registrations wrt the Master date.
  '''

  main()
