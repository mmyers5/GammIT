import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import os

def meantime():
  '''
  Get outputs from get_meancov and get_sumtime for grb
  ---INPUTS---
  grb: str
    see get_sumtime
  ---OUTPUTS---
  will print exposure time, frame time, and number of frames
  then mean coverage and standard deviation of coverage
  '''
  grb = raw_input("GRB: ")
  expsumtime,framsumtime,n,grbDir,chan = get_sumtime(grb)
  x0,y0 = get_x0y0(grb,chan)
  meanCov,stdCov,driz = get_meancov(grbDir,x0,y0)
  print "RESULTS!!"
  print "exposure sum time: {exp}".format(exp=expsumtime)
  print "frame sum time: {fram}".format(fram=framsumtime)
  print "number of frames: {num}".format(num=n)
  print "mean coverage: {mean} err: {std}".format(mean=meanCov,std=stdCov)
  print "drizzle factor: {driz}".format(driz=driz)
  
def get_sumtime(grb):
  '''
  Get the total exposure times and frame times for a given mosaic
  ---INPUTS---
  grb: str
    the GRB name. Run from Mosaicking directory or else!
  ---OUTPUTS---
  expsumtime: flt
    the total exposure time for all drizzled frames
  framsumtime: flt
    the total frame time for all drizzled frames
  n: int
    the number of frames drizzled
  '''
  listDir,grbDir,chan = get_list(grb)
  if listDir==1: 
    print "No can do!"
    return 1
  f        = open(listDir, 'r')
  fList    = f.readlines()
  f.close()
  expsumtime  = 0.
  framsumtime = 0.
  n           = 0
  for cbcd in fList:
    cbcdFile = cbcd.strip('\n')
    exptime,framtime = get_time(cbcdFile)
    expsumtime  += exptime
    framsumtime += framtime
    n           += 1

  return expsumtime,framsumtime,n,grbDir,chan

def get_list(grb):
  '''
  Check the cbcd list that indicates the files used mosaic
  ---INPUTS---
  grb: str
    see get_sumtime
  ---OUTPUTS---
  listDir: str
    the path to the list containing the files used in mosaic 
  grbDir: str
    the path to the mosaic materials used for mosaicking
  chan: int
    the channel number queried
  '''
  print os.system('ls '+grb)
  rFile    = raw_input("r string: ")
  chan     = raw_input("channel number: ")
  if rFile=='combo':
    grbDir   = "{GRB}/{r}/ch{num}".format(GRB=grb,r=rFile,num=chan)
  else:
    grbDir   = "{GRB}/{r}/ch{num}/bcd".format(GRB=grb,r=rFile,num=chan)
  testDir  = os.system('ls '+grbDir+'/*list.txt')
  if testDir!=0:
    print "No can do!"
    return 1
  listFile = raw_input("list to use: ")   # file name, not full path
  if listFile=='':
    listFile = 'newcbcdlist.txt'
  listDir  = grbDir+'/'+listFile
  testDir  = test(listDir,grbDir)
  if listDir!=testDir:
    listDir=testDir
  
  return listDir,grbDir,chan

def test(listDir,grbDir):
  '''
  Makes sure that the list of cbcd files has full paths
  '''
  f = open(listDir,'r')
  cbcdFiles = f.readlines()
  f.close()
  newFilename = listDir
  for line in cbcdFiles:
    if line.startswith('SPITZER'):
      newFilename = listDir.split('.txt')[0]+'_addFull.txt'
      newFile     = open(newFilename,'a')
      break
  if newFilename!=listDir:
    for line in cbcdFiles:
      newLine = "/scratch2/GRBH-MZSFR/Mosaicking/{grb}/{line}".\
                format(grb=grbDir,line=line)
      newFile.write(newLine)
    f.close()
  return newFilename

def get_time(cbcdFile):
  '''
  Given a cbcd, get the exposure and frame time
  ---INPUTS---
  cbcdFile: str
    the path to the cbcd file used in mosaicking
  ---OUTPUTS---
  exptime: float
    the exposure time
  framtime: float
    the frame time
  '''
  hdulist  = pyfits.open(cbcdFile)
  exptime  = hdulist[0].header['EXPTIME']
  framtime = hdulist[0].header['FRAMTIME']
  
  return exptime, framtime

def get_driz(fname):
  '''
  Get the drizzle factor used for a mosaic
  '''
  hdulist = pyfits.open(fname)
  drizfac = hdulist[0].header['DRIZ_FAC']
  return drizfac

def get_meancov(mosaicDir,x0,y0):
  '''
  Get the mean coverage and std for a given mosaic
  ---INPUTS---
  mosaicDir: str
    the directory where the mosaic materials are located
  x0: int
    the x image position of grb centroid
  y0: int
    the y image position of grb centroid
  ---OUTPUTS---
  meanCov: float
    the mean coverage around a box centered somewhere
  stdCov: float
    the standard deviation of meanCov
  maxCov: float
    the maximum coverage number for whole mosaic
  '''
  print os.system('ls '+mosaicDir+'/Results* -d')
  results = raw_input("Results Directory: ")
  if results == '':
    results = 'ResultsCLa'
  fname   = "{mos}/{res}/Combine-mosaic/mosaic_cov.fits".\
            format(mos=mosaicDir,res=results)
  driz    = get_driz(fname)
  meanCov,stdCov = mean_std(fname,x0,y0)
  return meanCov,stdCov,driz

def get_x0y0(grb,chan):
  '''
  Get the x0 and y0 image positions for a given grb
  ---INPUTS---
  grb: str
    the grb name
  '''
  chanDict = {'1': 'a',\
              '2': 'b',\
              '3': 'c',\
              '4': 'd'}
  chanApert= chanDict[chan]
  grbApert = '/scratch2/GRBH-MZSFR/Analysis/{grb}/{grb}_phys{let}.reg'.\
             format(grb=grb,let=chanApert)
  testDir  = os.system('ls '+grbApert)
  if testDir!=0:
    print "No can do!"
    return 1
  regFile  = open(grbApert,'r')
  reg      = regFile.readlines()
  regFile.close()
  for line in reg:
    if line.startswith('circle'):
      inside  = line.split('(')
      outside = inside[1].split(')')
      nums    = outside[0].split(',')
      x0 = int(round(float(nums[0])))
      y0 = int(round(float(nums[1])))
      break
  return x0,y0
  
def mean_std(fname,x0,y0,l0=6):
  '''
  Take a square of length l0 centered around (x0,y0) and calculate the mean
  within
  ---INPUTS---
  fname: str
    the filename of the fits file 
  x0: int
    the x-coordinate of the centroid
  y0: int
    the y-coordinate of the centroid
  l0: int
    the half length of the square in arcseconds
  ---OUTPUTS---
  meanCov: float
    the mean of the values within the box
  stdCov: float
    the standard deviation of the values within the box
  maxCov: float
  the max of coverage values for whole mosaic
  '''
  f = pyfits.open(fname)
  data = f[0].data
  header = f[0].header
  f.close()
  Y = np.shape(data)[0]
  X = np.shape(data)[1]
  lp = int((l0/0.4)/2.)
  y,x = np.ogrid[0:Y, 0:X]

  xLow = x >= x0 - lp
  xHigh = x <= x0 + lp
  yLow = y >= y0 - lp
  yHigh = y <= y0 + lp
  mask = xLow*xHigh*yLow*yHigh
  #data[mask]=0.

  meanCov = np.mean(data[mask])
  stdCov = np.std(data[mask])
  #fits.writeto(fname.strip('.fits')+'_test.fits',data,header,clobber=True)
  return meanCov,stdCov

def plotato(fname):
  f = pyfits.open(fname)
  data = f[0].data
  header = f[0].header
  f.close()
  plt.imshow(data, origin='lower')
  plt.show()
