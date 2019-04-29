from subprocess import Popen,PIPE
import numpy as np
from astropy.io import fits
import os, re
from astropy.cosmology import WMAP9 as cosmo

def calc_phot(pxscal,correction,bkgSubtractedCounts,srcUnc,srcPx,bkgUnc,bkgPx):
  '''
  Calculates flux density per frequency in microjansky
  ===INPUTS===
  pxscal: float
    pixel scale in arcseconds per pixel
  correction: float
    aperture correction
  bkgSubtractedCounts: float
    the background subtracted counts in the aperture
  srcUnc: float
    counts in aperture of square uncertainty map
  srcPx: int
    number of pixels in aperture
  bkgUnc: float
    counts in annulus of square uncertainty map
  bkgPx: int
    numbers of pixels in annulus (including clipped pixels)
  ===OUTPUTS===
  flx_uJy: float
    flux density per frequency in microjansky
  flx_unc_uJy: float
    flux density per frequency uncertainty in microjansky
  '''
  MJy_steradian    = 23.50443 #uJy/square-arc-second
  uJy_px  = MJy_steradian*pxscal**2
  flx_uJy = correction * bkgSubtractedCounts * uJy_px
  flx_unc_uJy_Sq = correction * 3 * (srcUnc*(uJy_px**2)/1000) + \
                   ((srcPx/bkgPx)**2 * bkgUnc*(uJy_px**2)/1000)
  flx_unc_uJy = np.sqrt(flx_unc_uJy_Sq)
  return flx_uJy,flx_unc_uJy

def mab_Mab(inFile,outFile):
  inDat = np.genfromtxt(inFile,names=True,dtype=None)
  outWrite = open(outFile,'w')
  outWrite.write('GRB\tz\tMab_H\tMab_H_unc\n')
  for entry in inDat:
    grb = entry['GRB']
    mab = entry['m_H']
    mab_unc = entry['m_H_unc']
    z = entry['z']
    dL = cosmo.luminosity_distance(z).value*1e6
    Mab = mab - 5*np.log10(dL/10)
    outWrite.write('{grb}\t{z}\t{Mab}\t{Mab_unc}\n'.\
                   format(grb=grb,z=z,Mab=Mab,Mab_unc=mab_unc))
  outWrite.close()


def AB2uJy(ab):
  flx = 3631.*10**(-ab*2/5.)*1e6
  return flx

def AB2uJy_unc(ab,abunc):
  flx = AB2uJy(ab)
  flx_unc = abs(AB2uJy(ab)-AB2uJy(ab-abunc))
  return flx,flx_unc
  

def uJy2AB_unc(flx,unc):
  '''
  Given flux and uncertainty, return AB magnitude
  ===INPUTS===
  flx: float
    flux density per frequency in microjansky
  unc: float
    flux density per frequency uncertainty in microjansky
  ===OUTPUTS===
  mAB: float
    AB magnitude
  muncAB: float
    AB magnitude uncertainty
  '''
  mAB    = uJy2AB(flx)
  muncAB = abs(uJy2AB(flx)-uJy2AB(flx+unc))
  return mAB,muncAB

def uJy2AB(flx):
  '''Given flux in microjansky, return absolute magnitude.'''
  mAB  = (-5/2.)*np.log10(flx*1e-6/3631)
  return mAB
  
def sub_cnts(srcCnts,srcPix,bkgCnts,bkgPix):
    '''
    calculate background-subtracted counts in aperture
    ===INPUTS===
    srcCnts: float
      counts in source region of aperture
    srcPix: int
      number of pixels in source region of aperture
    bkgCnts: float
      counts in background region of aperture
    bkgPix: int
      number of pixels in background region of aperture
    ===OUTPUTS===
    sCnts: float
      background-subtracted counts. See Tanmoy's paper
    '''
    sCnts = srcCnts-(bkgCnts*srcPix/bkgPix)
    return sCnts

def get_phys(inFile,ra,dec,phys=True,LTV=False):
  '''
  Convert incoming (ra,dec) coordinates to physical coordinates
  ===INPUTS===
  inFile: str
    .fits file registered with the (ra,dec) coordinates
  radec: tuple
    right ascension and declination in degrees
  phys: bool
    if True, will return physical coordinates. Default is True.
  ===OUTPUTS===
  XY: tuple
    physical x-and-y-coordinate
  '''
  header = fits.open(inFile)[0].header
  # differences from physical and image coordinates
  # physical = image - [LTV1 or LTV2]
  try:
    LTV1, LTV2   = header['LTV1'], header['LTV2']
  except KeyError:
    LTV1, LTV2 = 0, 0

  # dump degree coordinates into a garbage file
  coords = open('temp.txt','w')
  coords.write('RA\tDEC\n{ra}\t{dec}'.format(ra=ra,dec=dec))
  coords.close()
  stdout = Popen('funsky {inFile} {temp}'.\
          format(inFile=inFile,temp='temp.txt'),\
                 shell=True,stdout=PIPE).stdout 
  output = stdout.read()
  os.remove('temp.txt')

  coords = re.split(' +',output)[1:3]
  xy = (float(coords[0]),float(coords[1]))  # image coordinates
  XY = (xy[0] - LTV1,xy[1]-LTV2)
  if LTV == True:
    return LTV1,LTV2
  elif phys == True:
    return XY
  elif phys == False:
    return xy

def get_cnts(inFile,XY,ap,oFile):
  '''
  Get counts in a region
  ===INPUTS===
  inFile: str
    .fits file of image in which to count values
  XY: float
    physical coordinates of centroid
  ap: str
    aperture designation (2-2-6 or 3-3-7)
  oFile: str
    garbage file into which funcnts will dump counts info
  ===OUTPUTS===
  Creates oFile. Use read_cnts to read it
  '''
  # define aperture. values are tuples of pixel sizes
  apDict = {'2-2-6':(6,18),'3-3-7':(9,21)}
  # definte funcnts inputs
  aperture = 'circle({X},{Y},{aperture})'.\
             format(X=XY[0],Y=XY[1],aperture=apDict[ap][0])
  annulus  = 'annulus({X},{Y},{c1},{c2})'.\
             format(X=XY[0],Y=XY[1],c1=apDict[ap][0],c2=apDict[ap][1])
  # funcnts counts regions and puts into oFile
  os.system('funcnts {inFile} "{aperture}" "{annulus}" > {oFile}'.\
            format(inFile=inFile,aperture=aperture,annulus=annulus,oFile=oFile))

def read_cnts(cntsFile):
  '''
  Read the output of funcnts
  ===INPUTS===
  cntsFile: str
    the counts file into which funcnts writes
  ===OUTPUTS===
  srcCnts: float
    the number of counts in the source region
  srcPix: float
    the number of pixels in the source region
  bkgCnts: float
    the number of counts in the background region
  bkgPix: float
    the number of pixels in the background region
  numClip: int
    the number of clipped pixels, if available
  '''
  myFile = open(cntsFile,'r')
  cnts   = np.asarray(myFile.readlines())
  myFile.close()
  # indexing specific to output
  source = np.where(cnts=='# source_data\n')[0][0]+3
  bkgrnd = np.where(cnts=='# background_data\n')[0][0]+3
  inFile = cnts[1].split('/')[1].strip('\n')
  GRB    = inFile.split('_')[0]
  chan   = inFile.split('_')[1]
  dirFile= '{GRB}/{inFile}'.format(GRB=GRB,inFile=inFile)
  srcCnt = float(re.split(' +',cnts[source])[2])
  srcPix = float(re.split(' +',cnts[source])[3].strip('\n'))
  bkgCnt = float(re.split(' +',cnts[bkgrnd])[1])
  bkgPix = float(re.split(' +',cnts[bkgrnd])[2].strip('\n'))

  # check if there are clipped pixels
  header = fits.open(dirFile)[0].header
  if 'CLIPNUM' in header.keys():
    numClip = header['CLIPNUM']
  else:
    numClip = 0

  return srcCnt,srcPix,bkgCnt,bkgPix,numClip
  
def get_acorr(ap,chan):
  '''
  Get Spitzer standard aperture correction
  ===INPUTS===
  ap: str
    3-3-7 or 2-2-6
  chan: int
    channel of Spitzer image (1,2,3,or 4)
  '''
  apert = {'2-2-6':[1.215,1.233,1.366,1.568],'3-3-7':[1.125,1.120,1.135,1.221]}
  aCorr = apert[ap][chan-1]

  return aCorr 
