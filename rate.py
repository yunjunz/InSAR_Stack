#! /usr/bin/env python
# Author: Heresh Fattahi

import h5py
import argparse
import configparser
import  datetime
import time
import numpy as np

#################################################################
def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='extracts the overlap geometry between master bursts')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True,
            help='Directory with the pair directories that includes dense offsets for each pair')
    parser.add_argument('-d', '--data', type=str, dest='data', default=None,
            help='Directory with the pair directories that includes dense offsets for each pair')
    parser.add_argument('-o', '--output', type=str, dest='output', required=True,
            help='Directory with the pair directories that includes dense offsets for each pair')
    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args=iargs)
    return inps

def designMatrix(h5File):

    h5=h5py.File(h5File, 'r')
    dsDates = h5.get('dateList')
    dateList = list(dsDates[:])
    A = np.ones((len(dateList),2))
    for i in range(len(dateList)):
        d = dateList[i].decode("utf-8")
        d = datetime.datetime(*time.strptime(d,"%Y-%m-%d %H:%M:%S")[0:6]) #.strftime('%Y%m%d')
        d = d.year + (d.month-1)/12.0 + (d.day-1)/365.0
        A[i,0] = d
    h5.close()  
    time_std = np.sum((A[:,0] - np.mean(A[:,0]))**2)
    return A, time_std

def estimateVelocity(inps):
    A, time_std = designMatrix(inps.input)
    h5=h5py.File(inps.input, 'r')

    if inps.data:
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
    dsv = h5out.create_dataset('velocity', shape=(Ny,Nx),dtype=np.float32)
    dsv_std = h5out.create_dataset('std_velocity', shape=(Ny,Nx),dtype=np.float32)

    for i in range(Ny):
        print(i, 'out of ',Ny)
        L = data[:,i,:]        
        ts = np.dot(A1, L)
        r = L-np.dot(A,ts)
        std_v = np.sqrt(np.sum(r**2,0)/(Nz-2.)/time_std) 
        dsv[i,:] = ts[0,:]
        dsv_std[i,:] = std_v
    h5out.close()
    h5.close()
 
def main(iargs=None):
    inps = cmdLineParse(iargs)
    estimateVelocity(inps)


if __name__ == '__main__' :
  ''' 
  invert a network of the pair's mis-registrations to
  estimate the mis-registrations wrt the Master date.
  '''

  main()

