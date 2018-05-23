import numpy as np
import star_parser as sp
from os.path import isfile

class GalfitConfig(object):
  """For initializing/modifying my galfit objects.
  :param GRB:
    The GRB string. Assumes my custom directory layout.
  
  :param ch:
    Spitzer channel integer. Assumes my custom directory layout.
  """
  def __init__(self, GRB, ch):
    self.GRB = GRB
    self.ch  = ch
    
    self.spitz_star = '{GRB}/{GRB}_ch{ch}_maic.fits.stars'.\
                      format(GRB=GRB, ch=ch)
    self.aglow_star = '{GRB}/{GRB}_aglow.fits.stars'.\
                      format(GRB=GRB)
    self.hst_star   = '{GRB}/{GRB}_ch{ch}_hst.wcs.fits.stars'.\
                      format(GRB=GRB, ch=ch)

    self.spitz_img  = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.\
                      format(GRB=GRB, ch=ch)
    self.aglow_img  = '{GRB}/{GRB}_aglow.fits'.\
                      format(GRB=GRB)
    self.hst_img    = '{GRB}/{GRB}_ch{ch}_hst.wcs.fits'.\
                      format(GRB=GRB, ch=ch)

  def read_init(self, finit):
    """Specify the objects from each star file to use. Make sure all of the
    relevant star files are available, and set the star parser object.

    :param finit:
      The galfit initialization file. Should contain the following columns:
      -component: index starting from 1
      -image: hst, spitz, or aglow
      -number: the source number from SExtractor
      -source_type: how the source should be modeled in galfit (currently
      supports psf and gaussian)

    :returns init:
      A structured array using the columns specified above.
    """

    init = np.genfromtxt(initFile, names=True, dtype=None, usecols=('component',
                                                                    'image',
                                                                    'number',
                                                                    'source_type'))
    for image in np.unique(init['image']):
      try:      # ensure star files exist for all that are being demanded
        star_file = getattr(self, '{}_star'.format(image))
      except AttributeError:
        raise AttributeError('Image type {} does not exist. '\
                             'Please use hst, spitz, or aglow'.format(image))
      else:     # get the star object (see star_parser)
        try:
          setattr(self, '{}_star_object'.format(image), sp.Star(star_file))
        except IOError:
          raise IOError('File {} does not exist. '\
                        'Please change star file for {}'.format(star_file,
                                                                image))
    return init

  def check_files(self):
    """Simple diagnostic to check if all star files and image files exist."""
    checklist = ['spitz_star', 'aglow_star', 'hst_star',
                 'spitz_img', 'aglow_img', 'hst_img']
    badfiles  = [(i, getattr(self, i)) for i in checklist\
                 if not isfile( getattr(self,i) )]
    if len(badfiles) == 0:
      print('All of your files exist! Hurray!')
    else:
      print('The following attribute(s) and file(s) do not exist.')
      print(badfiles)
