#!/usr/bin/env python3
# Author: Heresh Fattahi
# 2017

import os, imp, sys, glob
import argparse
import configparser
import  datetime
import time
import numpy as np
import h5py
from insarPair import insarPair
from insarStack import insarStack


PI = np.pi

#################################################################
def createParser():
    '''
    Create command line parser.
    '''

    parser = argparse.ArgumentParser( description='extracts the overlap geometry between master bursts')
    parser.add_argument('-i', '--input', type=str, dest='input', required=True,
            help='Input h5 file that includes pairs of InSAR data')
    parser.add_argument('-t', '--timeseries_output', type=str, dest='output', required=True,
            help='Name of output time-series file')
    parser.add_argument('-o', '--observation', type=str, nargs = '+', dest='observation', default=None,
            help='name of the observation dataset to invert if several datasets exist. If not specified then the first observation dataset is inverted')
    parser.add_argument('-p', '--platform_track', type=str, dest='platformTrack', default=None,
            help='name of the platform-track to invert if several platform-tracks exist.')
    parser.add_argument('-m', '--inversion_method', type=str, dest='method', default='lsq',
            help='inversion method : LSQ , WLSQ')
    parser.add_argument('-a', '--azimuth_looks', type=float, dest='azLooks', default=14.0,
            help='number of azLooks. required only for converting coherence to phase variance in the WLSQ approach')
    parser.add_argument('-r', '--range_looks', type=float, dest='rngLooks', default=4.0,
            help='number of azLooks. required only for converting coherence to phase variance in the WLSQ approach')
    parser.add_argument('-w', '--wavelength', type=float, dest='wavelength', required=True,
            help='output directory to save dense-offsets for each date with respect to the stack Master date')
    parser.add_argument('-d', '--reference_date', type=str, dest='referenceDate', default=None,
            help='Reference date for the time-series. If not specified, the first date is considered as reference date.')
    parser.add_argument('-s', '--scale', type=float, dest='scale', default=None,
                help='scale parameter to scale the observation vector. If not provided, the input is assumed as phase and the scale is defined as : wavelength/4./PI which converts phase to range')

    return parser

def cmdLineParse(iargs = None):
    '''
    Command line parser.
    '''

    parser = createParser()
    inps = parser.parse_args(args=iargs)

    return inps

def plotNetwork(h5File):
  tbase,dateList1,dateDict, masters, slaves = date_list(h5File)
  dateList =[datetime.datetime(*time.strptime(d,"%Y-%m-%d %H:%M:%S")[0:6]) for d in dateList1]

  for d in dateList1[1:]:
      if d not in slaves:
         print("A discontinuity in the interferogram network detected")
         print(d)
         
         #idx = dateList1.index(d)
         #for m in masters:
             

  import matplotlib.pyplot as plt
  offset = np.ones(len(dateList))

  fig1 = plt.figure(figsize=(10,4))
  ax = fig1.add_subplot(1,1,1)
#  ax.plot(dateList[:-1],2*offset[:-1],'ro',ms=10)
#  ax.plot(dateList[1:],offset[1:],'bo',ms=10)

  numPairs = len(masters)
  for i in range(numPairs):
      ndxt1 = dateList1.index(masters[i])
      ndxt2 = dateList1.index(slaves[i])
      ax.plot([dateList[ndxt1], dateList[ndxt2]], [2,1],'-ko',ms=10,mfc='red')


  ax.set_ylim([0, 3])
  if dateList[0].month>1:
     minX = datetime.datetime(dateList[0].year,dateList[0].month-1,dateList[0].day)
  else:
     minX = datetime.datetime(dateList[0].year-1,12,dateList[0].day)
  
  if dateList[-1].month<12:
     print(dateList[-1].month+1)
     maxX = datetime.datetime(dateList[-1].year,dateList[-1].month+1,dateList[-1].day)
  else:
     maxX = datetime.datetime(dateList[-1].year+1,1,dateList[-1].day)

  ax.set_xlim([minX,maxX])
  ax.set_xlabel('Time')
  #ax.set_ylabel('LOS range-change [meters]')
  for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
             item.set_fontsize(20)

  plt.savefig("Network.png")
  #plt.show()
  #print(stop)

def date_list(h5file):
  h5=h5py.File(h5file,'r')
  ds = h5['/platform-track'].get('pairs_idx')
  pairs = ds[:,:]
  numPiars = pairs.shape[0]
  dateList = []
  tbase = []
  masters = [None]*numPiars
  slaves = [None]*numPiars
  for i in range(numPiars):
      master = pairs[i,0].decode("utf-8")
      slave = pairs[i,1].decode("utf-8")
      if master not in dateList: dateList.append(master)
      if slave not in dateList: dateList.append(slave)
      masters[i]=master
      slaves[i]=slave

  dateList.sort()
  d1 = datetime.datetime(*time.strptime(dateList[0],"%Y-%m-%d %H:%M:%S")[0:6])
  for ni in range(len(dateList)):
    d2 = datetime.datetime(*time.strptime(dateList[ni],"%Y-%m-%d %H:%M:%S")[0:6])
    diff = d2-d1
    tbase.append(diff.days)

  dateDict = {}
  for i in range(len(dateList)): dateDict[dateList[i]] = tbase[i]

  return tbase,dateList,dateDict, masters, slaves

def design_matrix(h5File, referenceDate):
  tbase,dateList,dateDict, masters, slaves = date_list(h5File)
  numDates = len(dateDict)
  numPairs = len(masters)
  A = np.zeros((numPairs,numDates))
  B = np.zeros_like(A)
  C = np.zeros((numPairs,numDates-1))
  tbase = np.array(tbase)
  for ni in range(numPairs):
     ndxt1 = dateList.index(masters[ni])
     ndxt2 = dateList.index(slaves[ni])
     A[ni,ndxt1] = -1
     A[ni,ndxt2] = 1
     B[ni,ndxt1:ndxt2] = tbase[ndxt1+1:ndxt2+1]-tbase[ndxt1:ndxt2]
     C[ni,ndxt1:ndxt2] = 1
  #print('A',A)
  #print('%%%%%%%%%%%%%%% %%%%%')  
  refIndex = dateList.index(referenceDate)
  A = np.hstack((A[:,0:refIndex], A[:,(refIndex+1):]))
  #A = A[:,1:]
  B = B[:,:-1]

  print(C)
  return A, B, C



def Coherence2SigPhi(coherence,inps):

    # Calculating the theoretical variance of the 
    # phase based on the coherence of
    # the sub-band interferogran

    N = inps.azLooks*inps.rngLooks
    epsilon = 0.0001
    coherence[coherence==0] = epsilon

    # this is for sub-bands only (Gomba et al, 2015)    
    Sig_phi2 = 3 * (1-coherence**2)/ (2.0*N*coherence**2)
    Sig_phi2[Sig_phi2==0] = epsilon
    return Sig_phi2

def invert(inps, h5File, outName):
    tbase,dateList,dateDict, masters, slaves = date_list(h5File)
    numPairs = len(masters)
    ###############################
    # specifying the reference date for the timeseries
    if inps.referenceDate is None:
       referenceDate = dateList[0]
    else:
       referenceDate = inps.referenceDate
    print('The reference date for timeseries: ', referenceDate)
    print (dateList)
    refIndex = dateList.index(referenceDate)    
    timeIndex = [i for i in range(len(dateList))]
    timeIndex.remove(refIndex)
    ##############################

    A, B, C = design_matrix(h5File, referenceDate)

    print('shape of A : ', A.shape)
    print('Rank of A : ', np.linalg.matrix_rank(A))
    print('number of dates: ', len(dateList))
    h5 = h5py.File(h5File,'r')
    if inps.platformTrack is None:
       inps.platformTrack = list(h5.keys())[0]

    if inps.observation is None:
       print(list(h5[inps.platformTrack]['observations'].keys()))
       inps.observation = list(h5[inps.platformTrack]['observations'].keys())
    
    obsPath = '/' + inps.platformTrack + '/' + 'observations'

    data = h5[obsPath].get(inps.observation[0])
    Nz, Ny, Nx = data.shape
    if inps.method == 'sequential':
        print('inversion method sequential')
        A = C
        timeIndex = [i for i in range(len(dateList)-1)]

        Npar = A.shape[1]
    else:
        Npar = A.shape[1]+1

    A1 = np.linalg.pinv(A)
    A1 = np.array(A1,np.float32)
    print(Npar)
    print(timeIndex)
    
    ##########
    h5residual = h5py.File('residual.h5','w')
    h5out = h5py.File(outName,'w')    
    group = h5residual.create_group('platform-track')
    obsGroup = group.create_group('observations')

    if  inps.scale is None:
        scale = inps.wavelength/4.0/PI
    else:
        scale = inps.scale
    
    print('scale', scale)
    for observation in inps.observation:
        print ('inverting : ', observation , ' pairs from ', inps.platformTrack )
        dsr = obsGroup.create_dataset(observation, shape=(Nz, Ny, Nx),
                        dtype=np.float32) 
        dst = h5out.create_dataset('temporal_coherence_'+observation, shape=(Ny,Nx),dtype=np.float32)
        ds = h5out.create_dataset(observation, shape=(Npar,Ny,Nx),dtype=np.float32)
        print(ds.shape)
        data = h5[obsPath].get(observation)
        #phase2range = inps.wavelength/4.0/PI
        for i in range(Ny):
            print(i, 'out of ',Ny)
            L = data[:,i,:]
            ts = np.dot(A1, L)
            L_residual = L - np.dot(A,ts)
            dsr[:,i,:] =  L_residual
            dst[i,:] = np.absolute(np.sum(np.exp(1j*L_residual),0))/Nz
            #ts = np.vstack((np.zeros((1,ts.shape[1]), dtype=np.float32), ts))
            ds[timeIndex,i,:] = ts*scale
            #ds[:,i,:] = ts*phase2range
            

    dateListE = [d.encode("ascii", "ignore") for d in dateList]
    dateListE = np.array(dateListE)
    dsDateList = h5out.create_dataset('dateList', data=dateListE, dtype=dateListE.dtype)

    h5residual.close()
    h5out.close()
    #h5tempCoh.close()
    h5.close()

def invert_wlsq(inps, h5File, outName, observation=None):
    # inversion using weighted least squares
    print('design matrix')
    tbase,dateList,dateDict, masters, slaves = date_list(h5File)
    numPairs = len(masters)
    ###############################
    # specifying the reference date for the timeseries
    if inps.referenceDate is None:
       referenceDate = dateList[0]
    else:
       referenceDate = inps.referenceDate
    print('The reference date for timeseries: ', referenceDate)
    print (dateList)
    refIndex = dateList.index(referenceDate)
    timeIndex = [i for i in range(len(dateList))]
    timeIndex.remove(refIndex)
    ##############################

    A, B, C = design_matrix(h5File, referenceDate)

    if inps.method == 'WLSQ_sequential':
       A = C
       timeIndex = [i for i in range(len(dateList)-1)]
       print ('WLSQ_sequential')
    h5 = h5py.File(h5File,'r')
    if inps.platformTrack is None:
       inps.platformTrack = list(h5.keys())[0]

    if inps.observation is None:
       print(list(h5[inps.platformTrack]['observations'].keys()))
       inps.observation = list(h5[inps.platformTrack]['observations'].keys())

    obsPath = '/' + inps.platformTrack + '/' + 'observations'

    data = h5[obsPath].get(inps.observation[0])
    Nz, Ny, Nx = data.shape
    Npar = A.shape[1]

    qualityPath = '/' + inps.platformTrack + '/' + 'quality'
    Sigma = h5[qualityPath].get('sig_iono') 

    print('prepare output files')
    ##########
    observation = inps.observation[0] 
    h5residual = h5py.File('residual.h5','w')
    group = h5residual.create_group('platform-track')
    obsGroup = group.create_group('observations')
    dsr = obsGroup.create_dataset(observation, shape=(Nz, Ny, Nx), dtype=np.float32)

    ##########
    h5out = h5py.File(outName,'w')
    ds = h5out.create_dataset(observation, shape=(Npar,Ny,Nx),dtype=np.float32)
    dsq = h5out.create_dataset('quality',shape=(Npar,Ny,Nx),dtype=np.float32)
    dst = h5out.create_dataset('temporal_coherence_'+observation, shape=(Ny,Nx),dtype=np.float32)
    print('kronoker product of design matrix for one raw')
    #I = np.eye(Nx)
    Nx_slice=10
    ii_idx = range(Nx_slice*Npar)
    I = np.eye(Nx_slice)
    Ak = np.kron(I,A)
    phase2range = inps.wavelength/4.0/PI
    range_x = range(0, int(Nx/Nx_slice)*Nx_slice, Nx_slice)

    for i in range(Ny):
    #for i in range(100): 
      print(i, 'out of ',Ny)
      #for j in range(0,Nx,Nx_slice):
      for j in range_x:
        print(j)
        L = data[:,i,j:j+Nx_slice].flatten(1)
        Sig_phi = Sigma[:,i,j:j+Nx_slice].flatten(1)
        #Sig_phi2 = Coherence2SigPhi(C,inps)
        W = np.diag(1./Sig_phi**2)

        Cm = np.linalg.inv(np.dot(np.dot(Ak.T, W),Ak))
        B = np.dot(Cm, (np.dot(Ak.T,W)))

        ts = np.dot(B, L)
        L_residual = L - np.dot(Ak,ts)
        #ts = ts.reshape([Nx,Npar]).T
        ts = ts.reshape([Nx_slice,Npar]).T
        #Cm = np.sqrt(Cm[range(Nx*Npar),range(Nx*Npar)]).reshape([Nx,Npar]).T
        Cm = np.sqrt(Cm[ii_idx, ii_idx]).reshape([Nx_slice,Npar]).T
        ds[timeIndex,i,j:j+Nx_slice] = ts*phase2range
        dsq[timeIndex,i,j:j+Nx_slice] = Cm*phase2range
        #ds[1:,i,j:j+Nx_slice] = ts*phase2range
        #dsq[1:,i,j:j+Nx_slice] = Cm*phase2range

        dsr[:,i,j:j+Nx_slice] =  L_residual.reshape((Nx_slice, Nz)).T
        dst[i,j:j+Nx_slice] = np.absolute(np.sum(np.exp(1j*L_residual),0))/Nz
        #ds[:,i,:] = ts #*phase2range

    dateListE = [d.encode("ascii", "ignore") for d in dateList]
    dateListE = np.array(dateListE)
    dsDateList = h5out.create_dataset('dateList', data=dateListE, dtype=dateListE.dtype)

    h5residual.close()
    h5out.close()
    h5.close()


def invert_wlsq_coherence(inps, h5File, outName, observation=None):
    # inversion using weighted least squares
    print('design matrix')
    tbase,dateList,dateDict, masters, slaves = date_list(h5File)
    numPairs = len(masters)
    A,B = design_matrix(h5File)

    h5 = h5py.File(h5File,'r')
    data = h5['/platform-track/observations'].get(observation)
    coherence = h5['/platform-track/quality'].get('coherence')
    Nz, Ny, Nx = data.shape
    Npar = A.shape[1]
    #A1 = np.linalg.pinv(A)
    #A1 = np.array(A1,np.float32)
    print('prepare output files')
    ##########
    h5residual = h5py.File('residual.h5','w')
    group = h5residual.create_group('platform-track')
    obsGroup = group.create_group('observations')
    dsr = obsGroup.create_dataset(observation, shape=(Nz, Ny, Nx),
                        dtype=np.float32)

    ##########
    h5tempCoh = h5py.File('temporal_coherence.h5','w')
    dst = h5tempCoh.create_dataset('temporal_coherence', shape=(Ny,Nx),dtype=np.float32)
    ##########

    h5out = h5py.File(outName,'w')
    ds = h5out.create_dataset(observation, shape=(len(dateList),Ny,Nx),dtype=np.float32)
    dsq = h5out.create_dataset('quality',shape=(len(dateList),Ny,Nx),dtype=np.float32)
    
    print('kronoker product of design matrix for one raw')
    #I = np.eye(Nx)
    Nx_slice=10
    I = np.eye(Nx_slice)
    Ak = np.kron(I,A)
    #phase2range = inps.wavelength/4.0/PI
    for i in range(Ny):
    #for i in range(10): 
      print(i, 'out of ',Ny)
      for j in range(0,Nx,Nx_slice):
        print(j) 
        L = data[:,i,j:j+Nx_slice].flatten(1)
        C = coherence[:,i,j:j+Nx_slice].flatten(1)
        Sig_phi2 = Coherence2SigPhi(C,inps)
        W = np.diag(1./Sig_phi2)
        
        Cm = np.linalg.inv(np.dot(np.dot(Ak.T, W),Ak))
        B = np.dot(Cm, (np.dot(Ak.T,W)))
        
        ts = np.dot(B, L)
        L_residual = L - np.dot(Ak,ts)
        #ts = ts.reshape([Nx,Npar]).T
        ts = ts.reshape([Nx_slice,Npar]).T
        #Cm = np.sqrt(Cm[range(Nx*Npar),range(Nx*Npar)]).reshape([Nx,Npar]).T
        Cm = np.sqrt(Cm[range(Nx_slice*Npar),range(Nx_slice*Npar)]).reshape([Nx_slice,Npar]).T
        ds[1:,i,j:j+Nx_slice] = ts
        dsq[1:,i,j:j+Nx_slice] = Cm

        dsr[:,i,j:j+Nx_slice] =  L_residual.reshape((Nx_slice, Nz)).T
        dst[i,j:j+Nx_slice] = np.absolute(np.sum(np.exp(1j*L_residual),0))/Nz
        #ds[:,i,:] = ts #*phase2range
        
    dateListE = [d.encode("ascii", "ignore") for d in dateList]
    dateListE = np.array(dateListE)
    dsDateList = h5out.create_dataset('dateList', data=dateListE, dtype=dateListE.dtype)

    h5residual.close()
    h5out.close()
    h5tempCoh.close()
    h5.close()


def main(iargs=None):

  inps = cmdLineParse(iargs)
  h5File = os.path.abspath(inps.input)
  outDir = os.path.dirname(inps.output)
  outDir = os.path.abspath(outDir)
  if not os.path.exists(outDir):
      os.makedirs(outDir)

  plotNetwork(h5File)
  
  timeseriesFile = os.path.join(outDir, os.path.basename(inps.output))

  if inps.method == 'WLSQ' or inps.method == 'WLSQ_sequential':
      print('inversion using Weighted Least Squares')
      h5Timeseries = invert_wlsq(inps, h5File, timeseriesFile, observation=inps.observation)

  else:
      print('inversion using Least Squares')
      h5Timeseries = invert(inps, h5File, timeseriesFile)

  #writeDateOffsets(inps, h5Timeseries)

if __name__ == '__main__' :
  ''' 
  invert a network of the pair's mis-registrations to
  estimate the mis-registrations wrt the Master date.
  '''

  main()









