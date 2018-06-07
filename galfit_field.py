import numpy as np
from galfit_config import GalfitConfig
import star_parser as sp
from subprocess import Popen, PIPE
from astropy.io import fits
import os, re
# to update constraints, just read in the old input file, old constriants file,
# and new input file. to change from psf to gaussian or vice versa, provide
# component number, initial input param text file, and new input file.

class GalfitSource(object):

  params = {'psf': ('x_world', 'y_world', 'erra_world'),
            'gaussian': ('x_world', 'y_world', 'erra_world',
                         'fwhm_world', 'theta_world', 'elongation')}

  def __init__(self, number, source_type, star_object=None):
    if star_object != None:
      for param in self.params[source_type]:
        setattr(self, param, star_object.get_object(number, param)[0])
    else:
      self.x_world = number[0]
      self.y_world = number[1]
      self.erra_world = 0.    # maybe default this to size of psf or something

  def tie_unc(self, tied_img=None, target_pxscal=0.4):
    # could be redundant way of getting ccout file, but just in case images
    # don't have wcs in their filename...
    if tied_img == None:
      setattr(self, 'tie_rms', 0)
    else:
      tied_ccout = ''.join(tied_img.split('.wcs')).replace('fits', 'ccout')
      with open(tied_ccout, 'r') as f:
        rms_line = [i for i in f if 'fit rms' in i][0].split()
      rms_index = rms_line.index('rms:')
      xi, eta = float(rms_line[rms_index+1]), float(rms_line[rms_index+2])
      rms = np.linalg.norm((xi, eta))/target_pxscal
      setattr(self, 'tie_rms', rms)
    
  def transform_coords(self, target_img, target_pxscal=0.4):
    if not all(i in self.__dict__.keys() for i in self.params['psf']):
      raise AttributeError('Please provide bare minimum for coordinate transform.')
    # dump degree coordinates into a garbage file
    with open('temp.txt', 'w') as f:
      names   = '\t'.join(self.__dict__.keys())
      values  = '\t'.join( map(str, self.__dict__.values()) )
      f.write('\n'.join((names, values)))

    # read ra, dec
    stdout = Popen('funsky {infits} {temp} {ra} {dec}'.format(infits=target_img, 
                                                              temp='temp.txt',
                                                              ra='x_world:d',
                                                              dec='y_world:d'),
             shell=True, stdout=PIPE).stdout
    output = stdout.read()
    os.remove('temp.txt')
    coords = re.split(' +', output)[1:3]
    xy = ( float(coords[0]), float(coords[1]) )
    xy_err = (self.erra_world * 3600)/target_pxscal

    setattr(self, 'target_x', xy[0])
    setattr(self, 'target_y', xy[1])
    setattr(self, 'target_xy_err', xy_err)

    # see if additional data needed
    if all(i in self.__dict__.keys() for i in self.params['gaussian']):
      setattr(self, 'target_fwhm', (self.fwhm_world * 3600)/target_pxscal)
      # PA = world north - image north...usually.
      header = fits.open(target_img)[0].header
      try:
        PA = header['PA']
      except KeyError:
        PA = 0
      # I think Spitzer images are backwards as hell
      setattr(self, 'target_theta', -1.*(self.theta_world - PA))
      setattr(self, 'target_axis_ratio', 1./(self.elongation))

class GalfitField(GalfitConfig):
  """For getting the field of objects to be used by galfit and transforming
  coordinates between images.

  :param init_file:
    Filename that specifies SEXtracted sources or coordinates from reference 
    (or target) images and the type of star to which each source is modeled.
    Can also specify coordinate locations in wcs using comma-separated value,
    but will default to using this coordinate as a psf, the implication being
    that SExtractor didn't extract this source for whatever reason.
  """
  def __init__(self, GRB, ch, init_file):
    GalfitConfig.__init__(self, GRB, ch)
    self.init_file = init_file
    self.read_init()
    self.transform_init()

  def read_init(self):
    usecols = ('component', 'image', 'number', 'source_type')
    dtype   = ('<i8', '|S5', '|S20', '|S10')
    init = np.genfromtxt(self.init_file, names=True, usecols=usecols,
                         dtype=dtype)

    for image in np.unique(init['image']):
      try:
        star_file = getattr(self, '{}_star'.format(image))
      except AttributeError:
        raise AttributeError('Image type {} does not exist.',
                             'Please use hst, spitz, or aglow'.format(image))
      
      try:
        setattr(self, '{}_star_object'.format(image), sp.Star(star_file))
      except IOError:
        raise IOError('File {} does not eist.',
                      'Please change star file for {}'.format(star_file, image))
      else:
        self.init = init

  def transform_init(self):
    for component in self.init:
      if ',' in component['number']:
        star = GalfitSource(component['number'].split(','), 'psf')
      else:
        star_object = getattr(self, '{}_star_object'.format(component['image']))
        star = GalfitSource(int(component['number']), component['source_type'], 
                            star_object)
      star.transform_coords(self.spitz_img)
      # calculate tie uncertainty, primary contributor to constraint
      # assumes hst is tied to mosaic and mosaic is tied to afterglow
      # obviously a coordinate in its own image gives tie_unc = 0
      if component['image'] == 'aglow':
        tied_img = self.spitz_img
      elif component['image'] == 'hst':
        tied_img = self.hst_img
      else:
        tied_img = None
      star.tie_unc(tied_img)
      # by this point, available star keys are target_xy_err, target_x,
      # target_y, target_fwhm, target_theta, and tie_rms, target_axis_ratio
      setattr(self, 'source{}'.format(component['component']), star)

  def write_init_inp(self):
    # intro
    inpfile, cofile, intro = intro_writer(self.GRB, self.ch)
    if os.path.isfile(inpfile):
      yn = raw_input("Do you want to overwrite the existing file {}? (y/n): ".\
                     format(inpfile))
      if yn.lower() == 'y':
        print("Right-o! Replacing {}".format(inpfile))
      elif yn.lower() == 'n':
        inpfile = raw_input("Please specify a new input file: ")
      else:
        print("I don't know what you are saying. I'm outta here!")
        return
    if os.path.isfile(cofile):
      yn = raw_input("Do you want to overwrite the existing file {}? (y/n): ".\
                     format(cofile))
      if yn.lower() == 'y':
        print("Right-o! Replacing {}".format(inpfile))
      elif yn.lower() == 'n':
        cofile = raw_input("Please specify a new constraints file: ")
      else:
        print("I don't know what you are saying. I'm outta here!")
        return

    # body and constraints
    body, const = [], []
    for component in self.init:
      compnum = component['component']
      stardata = getattr(self, 'source{}'.format(compnum))
      if compnum == 1:
        skip = 1
      else:
        skip = 0

      bodystrheader = '# Component number: {}'.format(compnum)

      if component['source_type'] == 'psf':
        bodystr = psf_writer(stardata.target_x, stardata.target_y, skip)
      elif component['source_type'] == 'gaussian':
        bodystr = gauss_writer(stardata.target_x, stardata.target_y, 
                               stardata.target_fwhm,
                               stardata.target_axis_ratio,
                               stardata.target_theta, skip)
      bodystr = '\n'.join((bodystrheader, bodystr))
      body.append(bodystr)

      if stardata.target_xy_err != 0 or stardata.tie_rms != 0:
        err = np.linalg.norm((stardata.target_xy_err, stardata.tie_rms))
        conststr = constraints_writer(compnum, err)
        const.append(conststr)

    body = '\n'.join(body)
    const = ''.join(const)

    # sky
    skyheader = '# Component number: {}'.format(compnum+1)
    sky = sky_writer()
    sky = '\n'.join((skyheader, sky))
    
    # write inp 
    with open(inpfile, 'w') as f:
      f.write('\n'.join((intro, body, sky))) 

    # write constraints
    with open(cofile, 'w') as f:
      f.write(const)

def intro_writer(GRB, ch):
  chstyle = {1:'a0', 2:'b0', 3:'c0', 4:'d0'}
  with open('template.inp', 'r') as f:
    intro_raw = ''.join(f.readlines())
  cofile = '{}_ch{}.constraints'.format(GRB, ch)
  inpfile = '{}/{}.inp'.format(GRB, chstyle[ch])
  repls = (('GRB', str(GRB)), ('CH', str(ch)), ('FNUM', chstyle[ch]),
           ('constraintFile', cofile))
  intro = reduce(lambda i, j: i.replace(*j), repls, intro_raw)

  # fix cofile
  cofile = '{GRB}/{GRB}_ch{ch}.constraints'.format(GRB=GRB, ch=ch)
  return inpfile, cofile, intro
  
def psf_writer(x, y, skip):
  with open('template_psf.inp', 'r') as f:
    psf_raw = ''.join(f.readlines())
  repls = ('XX', str(x)), ('YY', str(y)), ('SKIP', str(skip))
  psf = reduce(lambda i, j: i.replace(*j), repls, psf_raw)
  return psf

def gauss_writer(x, y, fwhm, axrat, pa, skip):
  with open('template_gauss.inp', 'r') as f:
    gauss_raw = ''.join(f.readlines())
  repls = (('XX', str(x)), ('YY', str(y)), ('SKIP', str(skip)),
           ('FWHM', str(fwhm)), ('AXISRATIO', str(axrat)), ('PA', str(pa)))
  gauss = reduce(lambda i, j: i.replace(*j), repls, gauss_raw)
  return gauss

def sky_writer():
  # Created separately in case I want to change sky things later
  with open('template_sky.inp', 'r') as f:
    sky_raw = ''.join(f.readlines())
  return sky_raw

def constraints_writer(compnum, err, xy=False):
  with open('template.constraints', 'r') as f:
    const_raw = ''.join(f.readlines())
  
  if xy==True:
    xerr = err[0]
    yerr = err[1]
  else:
    xerr, yerr = err, err
  repls = (('component', str(compnum)), ('x--', str(xerr)), ('x++', str(xerr)),
                                        ('y--', str(yerr)), ('y++', str(yerr)))
    
  const = reduce(lambda i, j: i.replace(*j), repls, const_raw)
  return const

def read_inp(inpfile):
  with open(inpfile, 'r') as f:
    inpraw = [line for line in f if '1)' in line][:-1]

  imgcoords = np.zeros((len(inpraw), 2))
  for i in range(len(inpraw)):
    x = float(inpraw[i].split()[1])
    y = float(inpraw[i].split()[2])
    imgcoords[i,0], imgcoords[i,1] = x, y
  return imgcoords

def read_consts(constsfile):
  # assumes each coordinate has equal constraints in xy in both directions
  consts = np.genfromtxt(constsfile, dtype=None,
                         names=('component', 'xy', 'err', 'err2'),
                         usecols=('component', 'err'))
  consts = consts[::2]
  # super annoying structured array
  return consts

def update_constraints(inp0file, con0file, inp1file, con1file):
  inp0 = read_inp(inp0file)
  inp1 = read_inp(inp1file)+10
  con0 = read_consts(con0file)
  const = []
  for compnum in con0['component']:
    inpindex = compnum-1
    err = con0[inpindex][1]
    x0, y0 = inp0[inpindex]
    x1, y1 = inp1[inpindex]
    xerr, yerr = np.abs(x1-x0)+err, np.abs(y1-y0)+err
    newconst = constraints_writer(compnum, (xerr, yerr), xy=True)
    const.append(newconst)
  const = ''.join(const)
  with open(con1file, 'w') as f:
    f.write(const)
