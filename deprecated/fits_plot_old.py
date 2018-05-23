from __future__ import division # confidence high 
from astropy.cosmology import WMAP9 as cosmo
from subprocess import Popen,PIPE
from collections import OrderedDict
import os,re,sys,math,pyfits,scipy.ndimage
import numpy as np
import matplotlib.pyplot as plt

MAX_REJECT = 0.5 
MIN_NPIXELS = 5 
GOOD_PIXEL = 0 
BAD_PIXEL = 1 
KREJ = 2.5 
MAX_ITERATIONS = 5 

def zscale (image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None): 
    """Implement IRAF zscale algorithm 
    nsamples=1000 and contrast=0.25 are the IRAF display task defaults 
    bpmask and zmask not implemented yet 
    image is a 2-d numpy array 
    returns (z1, z2) 
    """ 
 
    # Sample the image 
    samples = zsc_sample (image, nsamples, bpmask, zmask) 
    npix = len(samples) 
    samples.sort() 
    zmin = samples[0] 
    zmax = samples[-1] 
    # For a zero-indexed array 
    center_pixel = (npix - 1) // 2 
    if npix%2 == 1: 
        median = samples[center_pixel] 
    else: 
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1]) 
 
    # 
    # Fit a line to the sorted array of samples 
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT)) 
    ngrow = max (1, int (npix * 0.01)) 
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, KREJ, ngrow, 
                                             MAX_ITERATIONS) 
 
    if ngoodpix < minpix: 
        z1 = zmin 
        z2 = zmax 
    else: 
        if contrast > 0: zslope = zslope / contrast 
        z1 = max (zmin, median - (center_pixel - 1) * zslope) 
        z2 = min (zmax, median + (npix - center_pixel) * zslope) 
    return z1, z2 
 
def zsc_sample (image, maxpix, bpmask=None, zmask=None): 
     
    # Figure out which pixels to use for the zscale algorithm 
    # Returns the 1-d array samples 
    # Don't worry about the bad pixel mask or zmask for the moment 
    # Sample in a square grid, and return the first maxpix in the sample 
    nc = image.shape[0] 
    nl = image.shape[1] 
    stride = max (1.0, math.sqrt((nc - 1) * (nl - 1) / float(maxpix))) 
    stride = int (stride) 
    samples = image[::stride,::stride].flatten() 
    return samples[:maxpix] 
     
def zsc_fit_line (samples, npix, krej, ngrow, maxiter): 
 
    # 
    # First re-map indices from -1.0 to 1.0 
    xscale = 2.0 / (npix - 1) 
    xnorm = np.arange(npix) 
    xnorm = xnorm * xscale - 1.0 
 
    ngoodpix = npix 
    minpix = max (MIN_NPIXELS, int (npix*MAX_REJECT)) 
    last_ngoodpix = npix + 1 
 
    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad 
    badpix = np.zeros(npix, dtype="int32") 
 
    # 
    #  Iterate 
 
    for niter in range(maxiter): 
 
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix): 
            break 
         
        # Accumulate sums to calculate straight line fit 
        goodpixels = np.where(badpix == GOOD_PIXEL) 
        sumx = xnorm[goodpixels].sum() 
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum() 
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum() 
        sumy = samples[goodpixels].sum() 
        sum = len(goodpixels[0]) 
 
        delta = sum * sumxx - sumx * sumx 
        # Slope and intercept 
        intercept = (sumxx * sumy - sumx * sumxy) / delta 
        slope = (sum * sumxy - sumx * sumy) / delta 
         
        # Subtract fitted line from the data array 
        fitted = xnorm*slope + intercept 
        flat = samples - fitted 
 
        # Compute the k-sigma rejection threshold 
        ngoodpix, mean, sigma = zsc_compute_sigma (flat, badpix, npix) 
 
        threshold = sigma * krej 
 
        # Detect and reject pixels further than k*sigma from the fitted

        lcut = -threshold 
        hcut = threshold 
        below = np.where(flat < lcut) 
        above = np.where(flat > hcut) 
 
        badpix[below] = BAD_PIXEL 
        badpix[above] = BAD_PIXEL 
         
        # Convolve with a kernel of length ngrow 
        kernel = np.ones(ngrow,dtype="int32") 
        badpix = np.convolve(badpix, kernel, mode='same') 
 
        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0]) 
         
        niter += 1 
 
    # Transform the line coefficients back to the X range [0:npix-1] 
    zstart = intercept - slope 
    zslope = slope * xscale 
 
    return ngoodpix, zstart, zslope 
 
def zsc_compute_sigma (flat, badpix, npix): 
 
    # Compute the rms deviation from the mean of a flattened array. 
    # Ignore rejected pixels 
 
    # Accumulate sum and sum of squares 
    goodpixels = np.where(badpix == GOOD_PIXEL) 
    sumz = flat[goodpixels].sum() 
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum() 
    ngoodpix = len(goodpixels[0]) 
    if ngoodpix == 0: 
        mean = None 
        sigma = None 
    elif ngoodpix == 1: 
        mean = sumz 
        sigma = None 
    else: 
        mean = sumz / ngoodpix 
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1))
        if temp < 0: 
            sigma = 0.0 
        else: 
            sigma = math.sqrt (temp) 
 
    return ngoodpix, mean, sigma 

def map_coords(x0,y0,w0,h0,deg):
  '''
  Map image coordinates to rotated image
  ===INPUTS===
  x0: float
    image x-coordinate
  y0: float
    image y-coordinate
  w0: int
    the width of the image
  h0: int
    the height of the image
  deg: float
    the counter-clockwise rotation angle of image, in degrees
  '''

  cx0 = w0/2              # center of image, x-coordinate
  cy0 = h0/2              # center of image, y-coordinate
  rad = math.radians(deg) # degree in radians
  sin = np.sin(rad)
  cos = np.cos(rad)
  M = np.array([[cos , sin , (1-cos)*cx0-sin*cy0],\
                [-sin, cos , sin*cx0+(1-cos)*cy0]])  # rotation matrix

  w = int(h0*np.abs(sin)+w0*np.abs(cos))   # new width of image
  h = int(h0*np.abs(cos)+w0*np.abs(sin))   # new height of image

  M[0,2] += w/2-cx0       # translated cx0 
  M[1,2] += h/2-cy0       # translated cy0
  v = [x0,y0,1]           # vectorized coordinates

  calc = np.dot(M,v)      

  return calc[0],calc[1]
  
def make_stamp(GRB,ch,X0,Y0,z,r,circol,letters='k'):
  '''
  Make a stamp of a GRB field, 6" a side
  ===INPUTS===
  GRB: str
    the grb
  ch: int
    Spitzer channel
  X0: float
    physical x-coordinate
  Y0: float
    physical y-coordinate
  z: float
    redshift
  circol: string
    the color of stamp
    red (r): HST
    yellow (y): Ground-based optical/NIR
    blue (b): Swift/UVOT/X-ray
  ===OUTPUTS===
  saves a file in plots/stamps. Be careful not to overwrite anything!
  '''
  inFits    = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.format(GRB=GRB,ch=ch)
  inhdulist = pyfits.open(inFits)
  imgHdr    = inhdulist[0].header
  imgData   = inhdulist[0].data
  inhdulist.close()
  NAXIS1 = imgHdr['NAXIS1']
  NAXIS2 = imgHdr['NAXIS2']
  deg    = imgHdr['PA']
  LTV1   = imgHdr['LTV1']
  LTV2   = imgHdr['LTV2']
  x0  = X0 + LTV1
  y0  = Y0 + LTV2

  # flip, apply zscale, and rotate
  imgFlip = np.flipud(imgData)
  z1,z2   = zscale(imgFlip,nsamples=600,contrast=0.25)
  imgRot  = scipy.ndimage.rotate(imgFlip,deg,prefilter=True)
  x,y     = map_coords(x0,NAXIS2-y0,NAXIS1,NAXIS2,deg)    #NAXIS2-y0 because of flip

  # declare figure objects
  fig = plt.figure()
  ax  = fig.gca()
  box    = 3.0/0.4  # 6" a side
  border = 0.5/0.4  # inset 0.5"
  ax.set_xlim((x-box,x+box))
  ax.set_ylim((y-box,y+box))
  ax.axis('off')
  circle  = plt.Circle((x,y),r/0.4,color=circol,lw=3.0,fill=False)
  ax.imshow(imgRot,vmin=z1, vmax=z2, cmap='gray_r',\
            interpolation='nearest')
  ax.add_artist(circle)

  # add stamp information
  GRBtag = '{GRB}\nIRAC - ch{ch}\nz = {z}'.format(GRB=GRB,ch=ch,z=z)
  ax.text(x-box+border,y+box-border,GRBtag,\
          ha='left',va='top',fontsize=22,fontweight='bold',color=letters)

  kpc_arcsec = cosmo.kpc_proper_per_arcmin(z).value/60.
  scale = '{num} kpc'.format(num=str(round(kpc_arcsec,1)))
  ax.hlines(y-box+2*border,x-box+border,x-box+border+1/0.4,lw=3.0)
  ax.text(x-box+border+0.5/0.4,y-box+2.2*border,scale,\
          ha='center',va='bottom',fontsize=22,fontweight='bold',color=letters)
  ax.text(x-box+border+0.5/0.4,y-box+1.7*border,'1"',\
          ha='center',va='top',fontsize=22,fontweight='bold',color=letters)

  if letters=='k':
    figFile = 'plots/stamps/{GRB}_ch{ch}_stamp.eps'.format(GRB=GRB,ch=ch)
  else:
    figFile = 'plots/stamps/{GRB}_ch{ch}_stamp_{col}.eps'.\
              format(GRB=GRB,ch=ch,col=letters)
  if os.path.isfile(figFile):
    os.remove(figFile)
  plt.savefig(figFile,format='eps',transparent=True,\
              pad_inches=0,dpi=plt.gcf().get_dpi())
  plt.close(fig)  

def hst_coords(GRB,ch,img):
  '''
  Get HST coordinates given image coordinates from Spitzer
  ===INPUTS===
  GRB: str
    the grb
  ch: int
    Spitzer image channel
  img: dict
    a dictionary for coordinates. keys are integers, indexed from 0
    values are tuples of (x,y)
  '''
  imgFile   = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.format(GRB=GRB,ch=ch)
  hstFile   = '{GRB}/{GRB}_hst.wcs.fits'.format(GRB=GRB)
#  hstFile   = '{GRB}/{GRB}_aglow.fits'.format(GRB=GRB)
  inhdulist = pyfits.open(imgFile)
  imgHdr    = inhdulist[0].header
  inhdulist.close()
  LTV1   = imgHdr['LTV1']
  LTV2   = imgHdr['LTV2']
  inhdulist = pyfits.open(hstFile)
  hstHdr    = inhdulist[0].header
  inhdulist.close()
  hstCoords = OrderedDict()        # the hst image coordinates
  deg = OrderedDict()              # the hst degree coordinates

  # for each region passed in
  for obj in img:
    # grab degree coordinates
    temp = open('temp.txt','w')
    x = img[obj][0]
    y = img[obj][1]
    temp.write('{x}\t{y}'.format(x=x,y=y))
    temp.close()
    stdout = Popen('funsky -r {imgFile} {temp}'.\
            format(imgFile=imgFile,temp='temp.txt'),\
                   shell=True,stdout=PIPE).stdout
    output = stdout.read()
    degree = re.split(' +',output)[1:3]
    deg[obj] = [float(degree[0]),float(degree[1])]

    # grab hst image coordinates
    temp = open('temp.txt','w')
    temp.write('{ra}\t{dec}'.format(ra=deg[obj][0],dec=deg[obj][1]))
    temp.close()
    stdout = Popen('funsky {hstFile} {temp}'.\
                   format(hstFile=hstFile,temp='temp.txt'),\
                   shell=True,stdout=PIPE).stdout
    output = stdout.read()
    image = re.split(' +',output)[1:3]
    hstCoords[obj] = [float(image[0]),float(image[1])]

  # calculate pixel scale
  scals = []

  for i in deg.keys():
    for j in deg.keys():
      if i==j:
        continue
      xTop = abs(deg[i][0]-deg[j][0]) * 3600
      xBot = abs(hstCoords[i][0]-hstCoords[j][0])  
      yTop = abs(deg[i][1]-deg[j][1]) * 3600
      yBot = abs(hstCoords[i][1]-hstCoords[j][1])  
      scals.append([xTop/xBot,yTop/yBot])
  scalAll = np.asarray(scals)
  pxscal = round(np.mean([scalAll[:,0],scalAll[:,1]]),2)
  os.remove('temp.txt')

  return hstCoords,pxscal

def seek_HST(GRB,hstFile):
  '''
  Given a .txt file indicating which HST images went with the models
  and GRB, return the HST image file filter
  ===INPUTS===
  GRB: str
    the grb name
  hstFile:
    the hstFile where images are keyed to GRBs. Necessary headers,
    GRB, HST
  ===OUTPUTS===
  hstFilt: str
    the hst filter associated with the hst image
  '''
  inDat   = np.genfromtxt(hstFile,names=True,dtype=None,delimiter='\t')
  hstFilt = inDat['HST'][inDat['GRB']==GRB]
  return hstFilt[0]

def make_model_stamp(GRB,ch,xy,incirc,ap=None):
  '''
  Generate stamps for modeled images
  ===INPUTS===
  GRB: str
    the grb
  ch: int
    Spitzer channel 
  xy: tuple
    image x and y coordinates
  ap: tuple
    to plot aperture in panel 2, set to (2.4,7.2) or (3.6,8.4)
  ===OUTPUTS===
  Writes a file to plots/stamps
  '''
  # Need original image, modeled image, HST image, and model .inp file
  imgFits = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.format(GRB=GRB,ch=ch)
  inhdulist = pyfits.open(imgFits)
  imgHdr    = inhdulist[0].header
  imgData   = inhdulist[0].data
  inhdulist.close()
  subFits = '{GRB}/{GRB}_ch{ch}_sub.wcs.fits'.format(GRB=GRB,ch=ch)
  inhdulist = pyfits.open(subFits)
  subData   = inhdulist[0].data
  inhdulist.close()
  hstFits = '{GRB}/{GRB}_hst.wcs.fits'.format(GRB=GRB,ch=ch)
#  hstFits = '{GRB}/{GRB}_aglow.fits'.format(GRB=GRB,ch=ch)
  inhdulist = pyfits.open(hstFits)
  hstHdr = inhdulist[0].header
  hstData   = inhdulist[0].data
  hNAXIS1=hstHdr['NAXIS1']
  hNAXIS2=hstHdr['NAXIS2']
  inhdulist.close()
  NAXIS1 = imgHdr['NAXIS1']
  NAXIS2 = imgHdr['NAXIS2']
  deg    = imgHdr['PA']
  modFile= '{GRB}/{GRB}_ch{ch}_model.inp'.format(GRB=GRB,ch=ch)
  modOpen= open(modFile,'r')
  models = modOpen.readlines()
  modOpen.close()

  # get coordinates of modeled regions in Spitzer
  coords = OrderedDict()  # fully rotated, image coordinates in Spitzer
  img    = OrderedDict()  # plain image coordinates, not rotated
  i = 0
  for line in np.arange(len(models)):
    if models[line].startswith('# Component number: ') and \
       'sky' not in models[line+1]:

      coordStr   = models[line+2]
      coordSplit = coordStr.split(' ')[2:4]
      x0 = float(coordSplit[0])
      y0 = float(coordSplit[1])
      img[str(i)] = [x0,y0]
      x,y = map_coords(x0,NAXIS2-y0,NAXIS1,NAXIS2,deg)
      coords[str(i)]=[x,y]
      i+=1
    # add the afterglow circle
    xyr = map_coords(xy[0],NAXIS2-xy[1],NAXIS1,NAXIS2,deg)
    img['aglow'] = [xyr[0],xyr[1]]
  # get coordinates of modeled regions in HST
  hstCoords,hstpx = hst_coords(GRB,ch,img)
  for hstObj in hstCoords:
    x0 = hstCoords[hstObj][0]
    y0 = hstCoords[hstObj][1]
    x,y = map_coords(x0,hNAXIS2-y0,hNAXIS1,hNAXIS2,0)
    hstCoords[hstObj][0]=x
    hstCoords[hstObj][1]=y

  # define figure parameters
  fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(12,4))
  plt.subplots_adjust(wspace=0.01)
  pxscal = 0.4
  box = 11/pxscal # 22" a side
  border = 0.5/pxscal # 0.5" inset
  i   = 1

  # plot original mosaic, modeled mosaic, and hst
  for image in [imgData,subData,hstData]:
    if i == 1:
      #ax = ax1
      ax = ax3
    elif i == 2:
      ax = ax2
    elif i == 3:
      #ax = ax3
      ax = ax1
      coords = hstCoords
      pxscal = hstpx
      box = 11/pxscal # 22" a side
      border = 0.5/pxscal # 0.5" inset
      deg = 0.

    # manipulate image (flip, scale, rotate)
    imageFlip = np.flipud(image)
    z1,z2     = zscale(imageFlip,nsamples=600,contrast=0.25)
    imageRot  = scipy.ndimage.rotate(imageFlip,deg,prefilter=False)
    ax.imshow(imageRot,vmin=z1,vmax=z2,cmap='gray_r',interpolation='nearest')

    # add circles for reach region
    for obj in coords.keys():
      if obj == '0':
        circol = '#ff7f0e'
      elif obj == 'aglow':
        circol = incirc
      else:
        circol = '#1f77b4'
      circle = plt.Circle((coords[obj][0],coords[obj][1]),1/pxscal,\
                          color=circol,lw=2.0,fill=False)
      ax.add_artist(circle)
      if obj != '0' and obj!='aglow':
        ax.text(float(coords[obj][0]),float(coords[obj][1])+1.5/pxscal,obj,\
                ha='center',va='bottom',fontsize=10,fontweight='bold',color=circol,\
                clip_on=True)

    # add stamp label to panels
    if i == 1:
      GRBtag = '{GRB} - ch{ch}\nOriginal Mosaic'.format(GRB=GRB,ch=ch)
      ax.text(coords['0'][0]-box+border,coords['0'][1]+box-border, GRBtag,\
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')
    if i == 2:
      GRBtag = '{GRB} - ch{ch}\nSubtracted Image'.format(GRB=GRB,ch=ch)
      ax.text(coords['0'][0]-box+border,coords['0'][1]+box-border, GRBtag,\
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')
    if i == 3:
      hstFilt= seek_HST(GRB,'hstfilts.txt')
      GRBtag = '{GRB} - ch{ch}\n{hstFilt}'.format(GRB=GRB,ch=ch,hstFilt=hstFilt)
      ax.text(coords['0'][0]-box+border,coords['0'][1]+box-border, GRBtag,\
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')
    
    if i == 2 and ap!=None:
      for region in ap:
        circle = plt.Circle((coords['0'][0],coords['0'][1]),region/pxscal,\
                           color='cyan',lw=1.0,fill=False)
        return coords
        ax.add_artist(circle)

    # more figure options
    ax.set_xlim((float(coords['0'][0])-box,float(coords['0'][0])+box))
    ax.set_ylim((float(coords['0'][1])-box,float(coords['0'][1])+box))
    ax.axis('off')
    i+=1
  
  if ap!=None:
    figFile = 'plots/stamps/{GRB}_ch{ch}_{ap}_modelstamp.eps'.\
              format(GRB=GRB,ch=ch,ap=int(ap[0]))
  else:
    figFile = 'plots/stamps/{GRB}_ch{ch}_modelstamp.eps'.\
              format(GRB=GRB,ch=ch)
  if os.path.isfile(figFile):
    os.remove(figFile)
  #plt.savefig(figFile,format='eps',dpi=plt.gcf().get_dpi(),transparent=True,\
  #            pad_inches=0)
  #plt.close(fig)
  plt.show(fig)
