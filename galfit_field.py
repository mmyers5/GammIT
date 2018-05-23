import numpy as np
from galfit_config import GalfitConfig
import star_parser as sp
from subprocess import Popen, PIPE
from astropy.io import fits
import os, re
# Might want to make GalfitSource and GalfitTarget separate to make
# re-calculating Target coordinate constraints easier. i.e. GalfitSource holds
# the position and geometry in the source img, plus additional shit necessary to
# transform. GalfitTarget holds only the position and geometry in the target
# img, so having two GalfitTargets to compare coordinates will allow
# re-calculation of constraints. This makes GalfitSource pretty...slim, so maybe
# don't make it a separate object. Keep it as a star_object, so only have
# GalfitTarget.

class GalfitSource(object):

  params = {'psf': ('x_world', 'y_world', 'erra_world'),
            'gaussian': ('x_world', 'y_world', 'erra_world',
                         'fwhm_world', 'theta_world', 'elongation')}

  def __init__(self, star_object, number, source_type):
    for param in self.params[source_type]:
      setattr(self, param, star_object.get_object(number, param)[0])

  def tie_unc(self, tied_img, target_pxscal=0.4):
    # could be redundant way of getting ccout file, but just in case images
    # don't have wcs...
    tied_ccout = ''.join(tied_img.split('.wcs')).replace('fits', 'ccout')
    with open(tied_ccout, 'r') as f:
      rms_line = [i for i in f if 'fit rms' in i][0].split()
    rms_index = rms_line.index('rms:')
    xi, eta = float(rms_line[rms_index+1]), float(rms_line[rms_index+2])
    rms = np.linalg.norm((xi, eta))/target_pxscal
    setattr(self, 'rms', rms)
    
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
      setattr(self, 'axis_ratio', 1./(self.elongation))

class GalfitField(GalfitConfig):
  """For getting the field of objects to be used by galfit and transforming
  coordinates between images.

  :param init_file:
    Filename that specifies SEXtracted sources or coordinates from reference 
    (or target) images and the type of star to which each source is modeled.
  """
  def __init__(self, GRB, ch, init_file):
    GalfitConfig.__init__(self, GRB, ch)
    self.init_file = init_file
    self.read_init()

  def read_init(self):
    usecols = ('component', 'image', 'number', 'source_type')
    init = np.genfromtxt(self.init_file, names=True, dtype=None,
                         usecols=usecols)

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
