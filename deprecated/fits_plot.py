from __future__ import division # confidence high 
from astropy.cosmology import WMAP9 as cosmo
from subprocess import Popen,PIPE
from collections import OrderedDict
import os,re,sys,math,pyfits,scipy.ndimage
import numpy as np
import matplotlib.pyplot as plt
import fast_phot as fap

# global variables for zscale functions
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
  """Map image coordinates to rotated image.
  :param x0:
    Image x-coordinate.

  :param y0:
    Image y-coordinate.

  :param w0:
    The full image width in pixels.

  :param h0:
    The full image height in pixels.

  :param deg:
    The counter-clockwise rotation angle of the image in degrees.

  :returns calc:
    An (x,y) tuple of rotated image coordinates.
  """

  # get proper rotation matrix, M
  cx0 = w0/2              # center of image, x-coordinate
  cy0 = h0/2              # center of image, y-coordinate
  rad = math.radians(deg) # degree in radians
  sin = np.sin(rad)
  cos = np.cos(rad)
  M = np.array([[cos , sin , (1-cos)*cx0-sin*cy0],\
                [-sin, cos , sin*cx0+(1-cos)*cy0]])

  w = int(h0*np.abs(sin)+w0*np.abs(cos))   # new width of image
  h = int(h0*np.abs(cos)+w0*np.abs(sin))   # new height of image

  M[0,2] += w/2-cx0       # translated cx0 
  M[1,2] += h/2-cy0       # translated cy0
  v = [x0,y0,1]           # vectorized coordinates

  calc = np.dot(M,v)      

  return calc[0],calc[1]

def get_file(stringkeys, prestring='{GRB}/{GRB}_ch{ch}', poststring='maic.wcs.fits'):  
  """Create filenames given basic information. Files are structured as follows:
  prefix_postfix. E.G.:
  {GRB}/{GRB}_{ch}_maic.wcs.fits

  :param prefix:
    A dict determining how to format a prestring using key,value pairs.
  
  :param postfix: (optional, default: None)
    A dict determining how to format a poststring using key,value pairs.

  :param prestring: (optional, default: '{GRB}/{GRB}_ch{ch}')
    The prefix string for creating the filename.

  :param postfix: (optional, default: maic.wcs.fits)
    The postfix string for finishing the filename.

  :returns filename:
    The filename string.
  """

  filename = '_'.join((prestring.format(**stringkeys),
                       poststring.format(**stringkeys)))
  return filename

def read_headers(fitsFile, headers, return_data=False):
  inHduList = pyfits.open(fitsFile)
  inHeader  = inHduList[0].header
  outHeader = []
  for header in headers:
    outHeader.append(inHeader[header])
  if return_data == True: 
    outHeader.append(inHduList[0].data)
  inHduList.close()
  return outHeader

def inp_coords(inpLines, *rotate):
  img    = OrderedDict() # un-rotated coords
  if len(rotate) == 3:
    rotateImg = True
    coords = OrderedDict()

  i = 0
  for line in range(len(inpLines)):
    if inpLines[line].startswith('# Component number: ') and 'sky' not in inpLines[line+1]:
      coordStr   = inpLines[line+2]
      coordSplit = coordStr.split(' ')[2:4]
      x0, y0     = float(coordSplit[0]), float(coordSplit[1])
      img[i]  = [x0, y0]

      if rotateImg == True:
        inv_y0  = rotate[1]-y0
        x, y = map_coords(x0, y0, *rotate)
        coords[i] = [x,y]

      i+=1

  if rotateImg == True:
    return img, coords
  else:
    return img

def make_stamp(GRB, ch, X0, Y0, z, r, circol, tcol='k', savefile=False, **scale_kwargs):
  """Make a stamp of a GRB field.

  :param GRB:
    The GRB string name. Used for creating filename.

  :param ch:
    The Spitzer channel. Used for creating filename.

  :param X0:
    The physical x-coordinate of the image.

  :param Y0:
    The physical y-coordinate of the image.

  :param z:
    Redshift of the image.

  :param r:
    The radius of the afterglow circle to be plotted.

  :param circol:
    The color assigned to the reference image used for the afterglow.
    -red, r: HST
    -yellow, y: Ground-based optical/NIR
    -blue, b: Swift/UVOT/X-ray

  :param tcol: (optional, default: 'k')
    The color used for the text in the stamps.

  :param savefile: (optional, default: False)
    Whether or not to save the figure.
    
  :param scale_kwargs: (optional)
    Keywords passed to zscale.
    -nsamples
    -contrast
    -bpmask
    -zmask

  :returns:
    A file in plots/stamps of the format...
  """
  colDict = {'r':'r', 'y':'m', 'b': 'c'} # red,red yellow,green blue,cyan
  
  # create filenames
  inFits  = get_file({'GRB':GRB, 'ch':ch})
  outFigs = get_file({'GRB':GRB, 'ch':ch, 'tcol': tcol},
                     prestring='plots/stamps/{GRB}_ch{ch}',
                     poststring='stamp_{tcol}.eps')
  print inFits
  
  # open and get basic image information
  inhdulist = pyfits.open(inFits)
  imgHdr    = inhdulist[0].header
  imgData   = inhdulist[0].data
  inhdulist.close()
  NAXIS1 = imgHdr['NAXIS1']
  NAXIS2 = imgHdr['NAXIS2']
  deg    = imgHdr['PA']
  try:
    LTV1, LTV2 = imgHdr['LTV1'], imgHdr['LTV2']
  except:
    LTV1, LTV2 = 0, 0
  x0  = X0 + LTV1
  y0  = Y0 + LTV2
  # flip, apply zscale, and rotate
  z1, z2  = zscale(imgData, **scale_kwargs)
  imgRot  = scipy.ndimage.rotate(imgData, -deg, prefilter=True)
  x,y     = map_coords(x0, y0, NAXIS1, NAXIS2, -deg)    #NAXIS2-y0 because of flip
  
  # declare figure objects
  fig = plt.figure()
  ax  = fig.gca()
  box = 6/0.4       # 6" a side, given pxscale = 0.4"/px
  border = box/12   # for placing texts in border
  ax.set_xlim((x-box/2,x+box/2))
  ax.set_ylim((y-box/2,y+box/2))
  ax.axis('off')
  
  # make plots
  ax.imshow(imgRot,vmin=z1, vmax=z2, cmap='gray_r',\
            interpolation='nearest',origin='lower')
  circle  = plt.Circle((x,y), r/0.4, color=colDict[circol], lw=4.0, fill=False)
  ax.add_artist(circle)

  # add stamp text information
  GRBtag = '{GRB}\nIRAC - ch{ch}\nz = {z}'.format(GRB=GRB,ch=ch,z=z)
  ax.text(x-box/2+border,y+box/2-border,GRBtag,\
          ha='left',va='top',fontsize=22,fontweight='bold',color=tcol)

  kpc_arcsec = cosmo.kpc_proper_per_arcmin(z).value/60.
  scale = '{num} kpc'.format(num=str(round(kpc_arcsec,1)))
  ax.hlines(y-box/2+2*border,x-box/2+border,x-box/2+border+1/0.4,lw=3.0,color=tcol)
  ax.text(x-box/2+border+0.5/0.4,y-box/2+2.1*border,scale,\
          ha='center',va='bottom',fontsize=22,fontweight='bold',color=tcol)
  ax.text(x-box/2+border+0.5/0.4,y-box/2+1.8*border,'1"',\
          ha='center',va='top',fontsize=22,fontweight='bold',color=tcol)

  # save the file
  if savefile==True:
    if os.path.isfile(outFigs):
      os.remove(outFigs)
    plt.savefig(outFigs,format='eps',transparent=True,\
                 pad_inches=0,dpi=plt.gcf().get_dpi())
    plt.close(fig)  

  return fig

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
      if (i==j) or (i=='aglow' or j=='aglow'):
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

def make_model_stamp(GRB,ch,xy,incirc,z,ap=None, c1=0.25, c2=0.25, savefile=False,
addnum=False, hstScale=0):
  """Generate stamps for modeled images
  
  :param GRB:
    The GRB string.

  :param ch:
    The channel integer of desired stamp.

  :param xy:
    The tuple containing the x and y image coordinates of centroid in units of
    degrees.

  :param incirc:
    The color of the marker used for the afterglow in the images.

  :param z:
     The redshift of the host.

  :param ap: (optional, default: None)
    The aperture used for display in the 'subtracted' panel. Either (2.4, 7.2)
    or (3.6, 8.4) for the 2-2-6 or 3-3-7 apertures.

  :param c1: (optional, default: 0.25)

  :param c2: (optional, default: 0.25)

  :param savefile: (optional, default: False)

  :param addnum: (optional, default: False)

  :returns fig:
    The figure object.
  """
  colDict = {'r':'r', 'y':'m', 'b': 'c'} # red,red yellow,green blue,cyan

  # Need original image, modeled image, HST image, and model .inp file
  imgFits = get_file({'GRB':GRB, 'ch':ch})
  imgFits   = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.format(GRB=GRB,ch=ch)
  NAXIS1, NAXIS2, deg, imgData = read_headers(imgFits, ['NAXIS1', 'NAXIS2', 'PA'],
                                              return_data=True)
  subFits   = get_file({'GRB':GRB, 'ch':ch},
                       poststring='sub.wcs.fits')
  NAXIS1, NAXIS2, deg, subData = read_headers(subFits, ['NAXIS1', 'NAXIS2', 'PA'],
                                              return_data=True)
  hstFits   = get_file({'GRB':GRB},
                       prestring = '{GRB}/{GRB}',
                       poststring='hst.wcs.fits')
  hNAXIS1, hNAXIS2, hstData = read_headers(hstFits, ['NAXIS', 'NAXIS2'],
                                           return_data=True)
  modFile = get_file({'GRB': GRB, 'ch':ch},
                     poststring='model.inp')
  with open(modFile) as mf:
    models = mf.readlines()

  # get coordinates of modeled regions in Spitzer
  img, coords = inp_coords(models, NAXIS1, NAXIS2, -deg)
  # add the afterglow circle
  xy_img = fap.get_phys(imgFits, xy[0], xy[1], phys=False)
  img['aglow']    = [xy_img[0], xy_img[1]]
  aglowx, aglowy  = map_coords(img['aglow'][0], img['aglow'][1],
                               NAXIS1, NAXIS2, -deg)
  coords['aglow'] = [aglowx, aglowy]

  # get coordinates of modeled regions in HST
  hstCoords,hstpx = hst_coords(GRB,ch,img)
  # rotate HST regions (not needed, just flip)
  #for hstObj in hstCoords.keys():
  #  hstCoords[hstObj][1] = hNAXIS2-hstCoords[hstObj][1]
  #for hstObj in hstCoords:
  #  x0 = hstCoords[hstObj][0]
  #  y0 = hstCoords[hstObj][1]
  #  x,y = map_coords(x0,hNAXIS2-y0,hNAXIS1,hNAXIS2,0)
  #  hstCoords[hstObj][0]=x
  #  hstCoords[hstObj][1]=y

  # define figure parameters
  fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(12, 4))
  plt.subplots_adjust(wspace=0.01)
  pxscal = 0.4    # for spitzer!

  for i, image in enumerate([imgData, subData, hstData]):
    # define data for each panel
    if i == 2:
      ax = ax0
      plotCoords = hstCoords
      plotPxscal = hstpx + hstScale
      #if i==3: plotPxscal -= .030
      plotDeg  = 0.
      contrast = c2
    elif i == 0 or i ==1:
      if   i == 0: 
        ax = ax1
      elif i == 1: 
        ax = ax2
      plotCoords = coords
      plotPxscal = pxscal
      plotDeg  = deg
      contrast = c1
  
    # define panel sizes
    values = np.array( plotCoords.values() )
    center = values[0,:]
    if i == 0:
      maxLength = (abs(center - values)).max()
    box = int(maxLength/plotPxscal) # (maxLength*pxscal)" a side
    if box < 21/plotPxscal:         # assure there's not too much zoom
      box = 21/plotPxscal
    border = box/24.

    # manipulate image (flip, scale, rotate)
    #imageFlip = np.flipud(image)
    z1, z2    = zscale(image, contrast=contrast)
    imageRot  = scipy.ndimage.rotate(image, -plotDeg, prefilter=True)
    ax.imshow(imageRot, vmin=z1, vmax=z2, cmap='gray_r',
              interpolation='nearest', origin='lower')

    # add circles for reach region
    for obj in plotCoords.keys():
      if obj == 'aglow':
        circol = colDict[incirc]
        if i == 0:
          marker = '+'
      else:
        circol = '#ff7f0e'
        marker='o'

      if (i != 1) and (obj == 0):
        continue
      elif i == 1: 
        for region in ap:
          apcol = '#1f77b4'
          circle = plt.Circle((plotCoords[0][0], plotCoords[0][1]),
                               region/plotPxscal, color=apcol, lw=1.0, fill=False)
          ax.add_artist(circle)
        continue
      elif obj == 'aglow' and i!=2:
        ax.plot(plotCoords[obj][0], plotCoords[obj][1], color=colDict[incirc],
                ms=12, mew=2, marker=marker)
      else:
        circle = plt.Circle((plotCoords[obj][0], plotCoords[obj][1]),
                             1/plotPxscal, color=circol, lw=2.0, fill=False)
        ax.add_artist(circle)

        if addnum == True:
          ax.text(float(plotCoords[obj][0]), float(plotCoords[obj][1]) + 1.5/plotPxscal,obj,
                  ha='center', va='bottom', fontsize=12, fontweight='bold',
                  color=circol, clip_on=True)

    # more figure options
    cx = plotCoords[0][0]
    cy = plotCoords[0][1]
    ymin, ymax = cy - box/2, cy + box/2
    xmin, xmax = cx - box/2, cx + box/2
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.axis('off')

    # add stamp label to panels
    if i == 2:
      hstFilt= seek_HST(GRB,'hstfilts.txt')
      GRBtag = '{GRB} - ch{ch}\n{hstFilt}'.format(GRB=GRB, ch=ch, hstFilt=hstFilt)
      ax.text(xmin + border, ymax - border, GRBtag,
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')
    if i == 0:
      GRBtag = '{GRB} - ch{ch}\nOriginal Mosaic'.format(GRB=GRB,ch=ch)
      ax.text(xmin + border, ymax - border, GRBtag,
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')
    if i == 1:
      GRBtag = '{GRB} - ch{ch}\nSubtracted Image'.format(GRB=GRB,ch=ch)
      ax.text(xmin + border, ymax - border, GRBtag,
              ha='left',va='top',fontsize=16,fontweight='bold',color='k')

    # add stamp line to panels
    kpc_arcsec = cosmo.kpc_proper_per_arcmin(z).value/60.
    scale = '{num} kpc'.format(num=str(round(5*kpc_arcsec,1)))
    ax.hlines(ymin + 2*border, xmin + border, xmin + border + 5/plotPxscal,
              lw=3.0, color='k')
    ax.text(xmin + border + 2.5/plotPxscal, ymin + 2.1*border, scale,
            ha='center', va='bottom', fontsize=12, fontweight='bold', color='k')
    ax.text(xmin + border + 2.5/plotPxscal, ymin + 1.8*border, '5"',
            ha='center', va='top', fontsize=12, fontweight='bold', color='k')
  
  if savefile == True:
    apDict = {'2.4':2,'3.6':3}
    if ap != None:
      apStr = apDict[str(ap[0])]
      figFile = 'plots/stamps/{GRB}_ch{ch}_{ap}_modelstamp.eps'.\
                 format(GRB=GRB,ch=ch,ap=apStr)
    else:
      figFile = 'plots/stamps/{GRB}_ch{ch}_modelstamp.eps'.\
                 format(GRB=GRB,ch=ch)
    if os.path.isfile(figFile):
      os.remove(figFile)
    fig.savefig(figFile,format='eps',dpi=plt.gcf().get_dpi(),transparent=True,\
                pad_inches=0)
    plt.close(fig)
  else:
    return fig
