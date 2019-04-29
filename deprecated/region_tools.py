import numpy as np
import fast_phot as fap
from collections import OrderedDict
import os
import astropy.io.fits as fits
import numpy.ma as ma

def region_analysis(inFile, x0, y0, r, ext=0):
  """Give information about object field.
  
  :param inFile:
    A fits image string.

  :param x0:
    Image x-coordinate in pixels.

  :param y0:
    Image y-coordinate in pixels.

  :param r:
    Radius about which to study the field in pixels.

  :param ext: (optional, default: 0)
    If the fits image has multiple extensions, specify the index extension.

  :returns:
    The minimum pixel value in a circle of radius r centered on (x0,y0)

  :returns:
    The maximum pixel value in a circle of radius r centered on (x0,y0)

  :returns rms:
    The RMS of the field.
  """
  f = fits.open(inFile)
  data = f[ext].data
  f.close()
  Y = np.shape(data)[0]
  X = np.shape(data)[1]
  reg = np.zeros((Y,X))
  y,x = np.ogrid[-y0:Y-y0,-x0:X-x0]
  mask = x*x + y*y >= r*r  # true outside source region
  mx = ma.masked_array(data,mask=mask)
  rms = np.sqrt((np.square(mx)).mean())
  return mx.min(),mx.max(),rms

def make_apertures(GRB,ra,dec):
  """
  Create ds9 regions for the aperture/annulus of a target.

  :param GRB:
    The grb name as a string. Presumably follows my specific directory layout.

  :param xy:
    The right ascension and declination of the GRB in degrees.

  :returns:
    Two region files in ./GRB, labeled 3-3-7.reg or 2-2-6.reg for the respective
    apertures. Place on image to calculate flux counts using that aperture.
    Primarily used as a visual diagnostic of the field.
  """
  # Spitzer apertures are 2.4" and 3.6", corresponding to 2 or 3 native pixels
  apDict = {'2-2-6': ('2.4"', '7.2"'),
            '3-3-7': ('3.6"', '8.4"')}

  for aperture in ('3-3-7','2-2-6'):
    reg = '{GRB}/{ap}.reg'.format(GRB=GRB, ap=aperture)
    # write some stuff
    with open(reg, 'w') as regFile:
      regFile.write('# Region file format: DS9 version 4.1\n')
      regFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
      regFile.write('fk5\n')
      ap = apDict[aperture][0]
      an = apDict[aperture][1]
      regFile.write('circle({ra},{dec},{ap}\n'.format(ra=ra,dec=dec,
                                                      ap=apDict[aperture][0]))
      regFile.write('-circle({ra},{dec},{ap} # background\n'.format(ra=ra,dec=dec,
                                                                    ap=apDict[aperture][0]))
      regFile.write('circle({ra},{dec},{an} # background\n'.format(ra=ra,dec=dec,
                                                                   an=apDict[aperture][1]))

def make_region_from_model(inpFile,outFile,hstcoords=None,color='green'):
  """Will make a ds9 region file in image coordinates given a galfit input file.

  :param inpFile:
    The galfit input file string.

  :param outFile:
    The output region file string.

  :param hstcoords: (optional, default: None)
     What to label each component. Deprecated usage.

  :param color: (optional, default: 'green')
    The color that you want the region circles to be. Superfluous!

  :returns:
    A region file in image coordinates that should have the modeled regions from
    the input file.
  """
  with open(inpFile, 'r') as f:
    models = f.readlines()
  with open(outFile, 'w') as outWrite:
    outWrite.write('# Region file format: DS9 version 3.0\n')
    outWrite.write('global color=green font="helvetica 10 normal" select=1 edit=1 move=1 delete=1 include=1 fixed=0\n')
    i=0
    for line in np.arange(len(models)):
      if models[line].startswith('# Component number: ') and 'sky' not in models[line+1]:
        coordStr   = models[line+2]
        coordSplit = coordStr.split(' ')[2:4]
        x0  = float(coordSplit[0])
        y0  = float(coordSplit[1])
        if hstcoords is None:
          tag = '# color='+color+' '+'text={'+str(i)+'}\n'
        else:
          tag = '# color='+color+' '+'text={'+str(hstcoords[i])+'}\n'
        outWrite.write('image;circle({x}, {y}, 2.5) '.format(x=x0,y=y0))
        outWrite.write(tag)
        i+=1

def sex_obj(inFile, obj, image=False, fwhm=False, spitzer=False):
  """Return centroid information from a .stars file. Will print additional
     things that I, for some reason, didn't think to include at the time.

  :param inFile:
    The input .stars file from runsex.

  :param obj:
    The integer corresponding to the object in .stars.

  :param image: (optional, default: False)
    If True, will return image coordinates. If False, will return world
    coordinates in degrees. 

  :param fwhm: (optional, default: False)
    If True, will return the FWHM in addition to the other return arguments.

  :param spitzer: (optional, default: False)
    If True, will return the axis angle of an extended object in Spitzer
    coordinates (for use with galfit, add 90). Only works if fwhm is True.

  :returns xy:
    The right ascension and declination of the object in degrees or pixels 
    and the RMS position error along major axis.
  """
  stars = np.loadtxt(inFile,skiprows=13)  # will change if SExtractor is configured differently
  row   = np.where(stars[:,0] == obj)[0]  
  if len(row) > 1:
    row = row[0]
  if image == True:
    xy = stars[row, [1,2,12]].flatten()
    xy[2] *= 3600/0.4       # convert a thing in degrees to pixels
  else:
    xy = stars[row, 10:13].flatten()

  # For FWHM inquiries
  if fwhm == True:
    if image==True:
      fwhm  = float(stars[row, 8])
    else:
      fwhm = float(stars[row, 13])
    theta = -1*float(stars[row, 14])
    axrat = float(stars[row, 7])
    if spitzer == True:
      theta = -1*theta
      if image==True:
        fwhm = float(stars[row,13])*3600/0.4
    print 'x(deg/pix): {x}\ny(deg/pix): {y}\nfwhm(deg/pix): {fwhm}\n theta(deg): {theta}\nb/a(pix): {ax}\n'.\
    format(x=xy[0], y=xy[1], fwhm=fwhm, theta=theta, ax=1/axrat)
  return xy

def init_inp(GRB,ch,num):
  """
  Initialize a galfit .inp file.

  :param GRB:
    The object, assuming my specific file structure.

  :param ch:
    The channel of the Spitzer image.

  :param num:
  """
  # labeling system for .inp and .inp files
  # a,b,c,d == ch1,ch2,ch3,ch4
  chInp = {'1':'a0','2':'b0','3':'c0','4':'d0'}

  objects = obj_maker(num)            # {compononent: (source, objnum)}
  ties    = tie_maker(objects)        # {source: tie uncertainty}
  fullObjects = unc_calc(GRB, ch, objects, ties)  # {component: (source,objnum,unc)}

  # Image coordinates for each object in Spitzer
  xy = model_obj_image(GRB,fullObjects,ch)

  # open and write initial parameters to input file by changing template
  # see psf_template.inp
  # get uncertainties and write out constraints file
  constraintFile  = '{GRB}/{GRB}_ch{ch}.constraints'.format(GRB=GRB,ch=ch)
  constraintFile2 = '{GRB}_ch{ch}.constraints'.format(GRB=GRB,ch=ch)
  write_constraints(fullObjects, constraintFile)  

  # fill in psf templates
  with open('/scratch2/GRBH-MZSFR/Analysis/template.inp', 'r') as f:
    intro = f.readlines()
  intro = [lib.replace('GRB', GRB) for lib in intro]
  intro = [lib.replace('CH', str(ch)) for lib in intro]
  intro = [lib.replace('FNUM', chInp[str(ch)]) for lib in intro]
  intro = [lib.replace('constraintFile', constraintFile2) for lib in intro]

  with open('/scratch2/GRBH-MZSFR/Analysis/psf_template.inp','r') as f:
    bodies = f.readlines()
  outInp = open('{GRB}/{inp}.inp'.format(GRB=GRB, inp=chInp[str(ch)]), 'w')
  outInp.write(''.join(intro))

  primary   = bodies.index('# primary\n')+1
  secondary = bodies.index('# secondary\n')+1
  sky = bodies.index('# sky\n')+1
  i = 1     # labeling component numbers

  # write input coordinates, changes pending re: freeing GRB host position? or
  # constrained?
  for key, val in xy.items():
    if i == 1:
      body = bodies[primary:secondary-1]
    else:
      body =  bodies[secondary:sky-1]
    body = [lib.replace('xx', str(val[0])) for lib in body]
    body = [lib.replace('yy', str(val[1])) for lib in body]
    outInp.write('# Component number: {i}\n'.format(i=i))
    outInp.write(''.join(body))
    print "{i}: {bod}".format(i=i, bod=''.join(body))
    i += 1

  # write sky model
  outInp.write('# Component number: {i}\n'.format(i=i))
  outInp.write(''.join(bodies[sky:]))
  outInp.close()

def tie_maker(objects):
  """Make a dictionary for calculating tie uncertainties. Will need spreadsheet
  of all of the xi and eta tie uncertainties from SExtractor. Can get from
  .ccout files.

  :param objects:
    See obj_maker.

  :returns ties:
    For each key, find the uncertainty in the tie between the key and the
    spitzer mosaic. Number is uncertainty in Spitzer pixels.
  """
  # Ensure not to calculate redundant ties
  sources = np.unique([ val[0] for key,val in objects.items() ])
  ties    = OrderedDict().fromkeys(sources)

  for key in ties.keys():
    if key == 'spitzer':
      xi, eta=(0, 0)
    else:
      print 'COMPONENT {key}'.format(key=key)
      xi  = input('xi unc (arcsec): ')
      eta = input('eta unc (arcsec): ')
    # Convert to Spitzer pixels, assuming 0.4"/pixel
    tieUnc    = (xi/0.4, eta/0.4)
    ties[key] = tieUnc

  return ties

def obj_maker(num):
  """Make the dictionary for use in parse_uncs and init_inp
  
  :param num:
    The number of entries in the dictionary that is desired.

  :returns objects:
    A dict of either tuples that correspond to image coordinates in the 
    Spitzer mosaic or SExtractor objects that correspond to either
    the HST or Spitzer images.
    Each key is an integer indicating the component number in the .inp file.
    The values are tuples, the first element of which is 'hst', 'spitzer', or
    'aglow' to indicate where the SExtractor object is located. The second
    element is the object number.
  """
  objects = OrderedDict.fromkeys(range(num))
  prev = '---'
  for key in objects.keys():
    source = raw_input('Component {i} source (spitzer,hst,aglow) {p}: '.\
                       format(i=key,p=prev))
    if source == '':
      source = objects[key-1][0]
    else:
      prev = source
    obj = input('Component {i} object number: '.format(i=key))
    objects[key] = (source,obj)
  return objects

def unc_calc(GRB,ch,objects,ties):
  """
  Calculate the full uncertainty of barycenter, including centroid uncertainty
  and tie uncertainty.
  
  :param GRB:
    The object.
  
  :param ch:
    The channel in Spitzer of observations.

  :param objects:
    See obj_maker.

  :params ties:
    See tie_maker.

  :returns objects:
    A modified objects dictionary. The keys are component numbers. The values
    are tuples. First element is the image source. Second element is the
    SExtractor object. Third element is the total uncertainty in the position of
    the object in Spitzer pixel coordinates, including the tie uncertainty and
    the barycenter uncertainty from SExtractor.
  """
  errSpitz    = parse_uncs(GRB, ch, objects) # centroid error in Spitzer pixels
  fullObjects = OrderedDict()
  for key, val in objects.items():
    tieUnc = ties[val[0]]
    objUnc = errSpitz[key][1]
    totUnc = np.linalg.norm((objUnc,tieUnc[0],tieUnc[1])) # quadrature sum
    fullObjects[key] = (val[0],val[1],totUnc)
  return fullObjects  

def parse_uncs(GRB,ch,objects):
  """Get the barycenter uncertainty of a SExtractor objects.

  :param GRB:
    The object.

  :param ch:
    The object channel observation in Spitzer.

  :param objects:
    See obj_maker.

  :returns errSpitz:
    The values of uncertainty in spitzer pixels for each component.
    Keys are the same as objects. Values are a tuple, the first element of which
    is the object number, the second element of which is the uncertainty in the 
    position in Spitzer pixels.
  """
  # Assume file structures
  prefix = '{GRB}/{GRB}'.format(GRB=GRB)
  suffix = {'hst':'hst',\
            'spitzer':'ch{ch}_maic'.format(ch=ch),\
            'aglow':'aglow'}
  errSpitz = OrderedDict()

  # Look for each component in their respective images
  for key,val in objects.items():
    if isinstance(val[1], (tuple,list)):  # I forget why this is here.
      errDeg = 0
    else:
      fName = '{pref}_{suff}.fits.stars'.format(pref=prefix, suff=suffix[val[0]])
      xy = sex_obj(fName, val[1])
    errSpitz[key] = val[0], xy[-1]*3600/0.4
    #  x, y, errDeg  = sex_obj(fName, val[1])    # units are in degrees
    #errSpitz[key]   = val[0], errDeg*3600/0.4   # units are no longer in degrees
  return errSpitz
    
def write_constraints(objects, outFile):
  """Write the constraints file used for galfit .inp file.
  
  :param objects:
    See obj_maker. If the uncertainty is 0, will not be constrained.

  :param outFile:
    The filename of the constraint file that you will write into.
  """
  with open('template.constraints', 'r') as f:
    tempUnc = f.readlines()

  with open(outFile, 'w') as outWrite:
    outWrite.write('# component\tparameter\tconstraint\tcomment\n')
    outWrite.write('# operation (see below)\trange\n\n')

    for key, val in objects.items():
      compNum = key + 1   # Components are numbered from 1, not 0
      if val[2] == 0:     # If there is no constraint
        continue
      else:
        lim = str(val[2]) # String-ify the constraint value

      body = [lib.replace('component', str(compNum)) for lib in tempUnc]
      body = [lib.replace('--', '-'+lim) for lib in body]
      body = [lib.replace('++', lim) for lib in body]
      outWrite.write(''.join(body))
    
def model_obj_image(GRB, objects, ch):
  """Give image coordinates given a GRB, several objects, and a target mosaic.

  :param GRB:
    The object.

  :param ch:
    The Spitzer band.
  
  :param objects:
    See obj_maker.


  :returns xy:
    xy coordinates of objects in Spitzer, using image coordinates.
  """
  # File assertions
  prefix = '{GRB}/{GRB}'.format(GRB=GRB)
  suffix = {'hst':'hst.wcs',\
            'spitzer':'ch{ch}_maic'.format(ch=ch),\
            'spitzerWcs':'ch{ch}_maic.wcs'.format(ch=ch),\
            'aglow':'aglow'}
  imgFile = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.\
            format(GRB=GRB,ch=ch)

  # Initialize dicts
  radec = OrderedDict()
  xy = OrderedDict()

  # Handle objects dict
  for key, val in objects.items():
    # Initialize .stars files
    starFile = '{pref}_{suff}.fits.stars'.format(pref=prefix, suff=suffix[val[0]])
    if val[0]=='spitzer' and not os.path.isfile(starFile):
      starFile = '{pref}_{suff}.fits.stars'.\
                 format(pref=prefix,suff=suffix['spitzerWcs'])
    
    if isinstance(val[1], int):
      #dictDex = str(val[1])
      radec[key] = sex_obj(starFile, val[1])  # radec in degrees
      if val[0]=='spitzer':
        xy[key] = sex_obj(starFile, val[1], image=True)[:-1]
      else:
        xy[key] = fap.get_phys(imgFile, radec[key][0], radec[key][1],
                             phys=False)

    elif isinstance( val[1], (tuple, list) ):     # If custom coordinates
      xy[key] = val[1]
  return xy

def get_coords_from_inp(inpFile): 
  """Given a galfit .inp file, return image coordinates of each component.
  Deprecated, look at fits_plot.
  
  :param inpFile:
    The galfit .inp file from which to read image coordinates.

  :returns inpXY:
    The xy image coordinates in a list.
  """
  models = open(inpFile, 'r').readlines()
  xyList = []
  for line in np.arange(len(models)):
    if models[line].startswith('# Component number: ') and 'sky' not in models[line+1]:
      coordStr   = models[line+2]
      coordSplit = coordStr.split(' ')[2:4]
      x0  = float(coordSplit[0])
      y0  = float(coordSplit[1])
      xyList.append(np.array((x0,y0)))
  return xyList
