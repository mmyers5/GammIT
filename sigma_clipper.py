import astropy.io.fits as fits
import numpy as np
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
from collections import OrderedDict

def pipe_mask(GRB,ch,x0,y0,ap,sub=False,sigHigh=2,sigLow=2):
  '''
  Generates sigma clipped images for a mosaic and its associated files
  ===INPUTS===
  GRB: str
    gamma-ray burst...
  ch: int
    channel of Spitzer image
  x0: int
    image x-coordinate of centroid of GRB
  y0: int
    image y-coordinate of centroid of GRB
  ap: str
    aperture and annulus size, of the from '2-2-6' or '3-3-7'
  sub: bool
    if the mosaic is a modeled image, set to True to perform clip on subbed
    image. Default is False.
  ===OUTPUTS===
  saves clipped images, new header titled 'CLIPNUM' to indicate
  number of pixels clipped. See make_mask to see output file names
  '''
  fileTags = OrderedDict({'flx':'maic.wcs','unc':'munc_sq.wcs','psf':'psf'})
  if sub == True:
    fileTags['flx']='sub.wcs'
  pref = '{GRB}/{GRB}_ch{ch}'.format(GRB=GRB,ch=ch)
  primo = '{pref}_{flx}.fits'.format(pref=pref,flx=fileTags['flx'])
  aperture = {'2-2-6':(6.5,18.5),'3-3-7':(9.5,21.5)}
  
  mask,maskArr,clippedNum,sLow,sHigh=make_mask(primo,x0,y0,aperture[ap][0],aperture[ap][1],
                                               sigHigh=sigHigh,sigLow=sigLow)
  for tag in fileTags:
    oldFile = '{pref}_{tag}.fits'.format(pref=pref,tag=fileTags[tag])
    newFile = '{pref}_{ap}_{tag}.fits'.format(pref=pref,ap=ap,tag=fileTags[tag])
    if tag == 'psf':
      psfMask = make_mask(oldFile,65,65,aperture[ap][0],aperture[ap][1],psf=True,\
                         sigHigh=sigHigh,sigLow=sigLow)
      apply_mask(oldFile,psfMask,maskArr,clippedNum,newFile)
    else:
      apply_mask(oldFile,mask,maskArr,clippedNum,newFile)
  return sLow
  
def apply_mask(fname,mask,maskArr,clippedNum,outFile):
  '''
  Applies a clipping mask to a .fits image
  ===INPUTS===
  fname: str
    .fits file location
  mask: array
    subset of .fits image to which masked values are applied
    consists of True/False values
  maskArr: array
    fill values applied to mask. must be of same shape as mask
    consists of 0/1 values, gets multiplied to original image
  clippedNum: int
    number of clipped pixels
  outFile: str
    clipped .fits file locations
  ===OUTPUTS===
  saves clipped image to outFile, new header titled 'CLIPNUM' to indicate
  number of pixels clipped
  '''
  f = fits.open(fname)
  fDat = f[0].data
  fHdr = f[0].header
  f.close()
  print mask.shape, maskArr.shape
  fDat[mask]=fDat[mask]*maskArr
  fHdr['CLIPNUM']=clippedNum
  fits.writeto(outFile,fDat,fHdr,clobber=True)

def make_mask(fname,x0,y0,r,R,psf=False,test=False,sigHigh=2,sigLow=2):
  '''
  r1 - > r source
  r2 - > R annulus
  Performs k-sigma clipping in background annulus of an image
  ===INPUTS===
    fname: str
      name of .fits file to which sigma-clipping as applied
    x0: int
      image x-coordinate of center of source region
    y0: int
      image y-coordinate of center of source region
    r: int
      radius of source region in pixels
    R: int
      radius of annulus region in pixels
    psf: bool
      if the input image is a psf, set to True to adjust x/y
      Default centroid is at (65,65), can adjust if needed.
  ===OUTPUTS===
  mask: array
    subset of .fits image to which masked values are applied
    consists of True/False values
  maskArr: array
    fill values applied to mask. must be of same shape as mask
    consists of 0/1 values, gets multiplied to original image
  clippedNum: int
    number of clipped pixels
  '''
  f=fits.open(fname)               # open the .fits file
  data=f[0].data                   # get image array
  header=f[0].header               # get image header
  f.close()                        # close the .fits file

  Y=np.shape(data)[0]              # "height" of image in pixels
  X=np.shape(data)[1]              # "width" of image in pixels
  reg=np.zeros((Y,X))              # array of width and height of img
  y,x=np.ogrid[-y0:Y-y0, -x0:X-x0] # arrays centered around x0,y0

  mask1 = x*x + y*y <= R*R       # True inside annulus circle region
  mask2 = x*x + y*y >= r*r       # True outside source region
  mask=mask1*mask2                 # True in annulus, False in source and else
  reg[mask]=data[mask]             # put image values in annulus region
  if psf==True:
    return mask
  
  # Median sigma-clip, 2-sigma for upper tail, 3-sigma if fail
  try:
    filtMed=sigma_clip(reg[mask],\
           sigma_lower=sigLow,\
           sigma_upper=sigHigh,\
           cenfunc=np.median,\
           iters=None)               # input array with clipped elements True
  except TypeError:
    sigLow +=1
    sigHigh+=1
    filtMed=sigma_clip(reg[mask],\
            sigma_lower=sigLow,\
            sigma_upper=sigHigh,\
            cenfunc=np.median,\
            iters=None)               # input array with clipped elements True

  nanFilt=np.ma.filled(filtMed,\
           fill_value=np.nan)             # input array with clipped elements 0
  clippedNum = np.sum(np.isnan(nanFilt))
  fillFilt = np.nan_to_num(nanFilt)
  maskArr = np.nan_to_num(fillFilt/fillFilt)
  if test==True:
    return mask,maskArr,clippedNum,data.sigLow,sigHigh
  else:
    return mask,maskArr,clippedNum,sigLow,sigHigh
