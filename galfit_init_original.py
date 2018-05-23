import numpy as np
import star_parser as sp

class GalfitConfig(object):
  """For initializing/modifying my galfit objects.
  :param GRB:
    The GRB string. Assumes my custom directory layout.
  
  :param ch:
    Spitzer channel integer. Assumes my custom directory layout.
  
  :param hst: (optional, default: 'hst')
    The string used to specify an hst image. Assumes my custom file naming
    conventions.

  :param spitz: (optional, default: 'maic')
    The string used to specify a spitzer image. Assumes my custom file naming
    conventions.

  :param aglow: (optional, default: 'aglow')
    The string used to specify an afterglow image. Assumes my custom file naming
    conventions.
  """
  def star_file(image):
    """My complicated way of setting star file names. So totally
       unnecessary"""
    image_name = '{}_star_file'.format(image)
    image_sffx = '_{}'.format(image)

    @property
    def prop(self):
      return getattr(self, image_sffx)

    @prop.setter
    def prop(self, sffx):
      new_file_name = '{GRB}/{GRB}_{sffx}'.format(GRB=self.GRB, sffx=sffx)
      setattr(self, image_name, new_file_name)
      setattr(self, image_sffx, sffx)
  
    return prop

  hst = star_file('hst')
  spitz = star_file('spitz')
  aglow = star_file('aglow')


  def __init__(self, GRB, ch,
               hst='hst.wcs.fits.stars',
               spitz='maic.fits.stars',
               aglow='aglow.fits.stars'):
    self.GRB = GRB
    self.ch  = ch
    self._hst, self._spitz, self._aglow = hst, spitz, aglow
    self.hst_star_file = '{GRB}/{GRB}_ch{ch}_{img}'.format(GRB=GRB, ch=ch,
                                                           img=hst)
    self.spitz_star_file = '{GRB}/{GRB}_ch{ch}_{img}'.format(GRB=GRB, ch=ch,
                                                             img=spitz)
    self.aglow_star_file = '{GRB}/{GRB}_{img}'.format(GRB=GRB, ch=ch,
                                                      img=aglow)
    self.constraints_file = '{GRB}/{GRB}_ch{ch}.constraints'.format(GRB=GRB, ch=ch)
    self._template_inp = 'template.inp'
    self._template_psf_inp = 'template_psf.inp'
    self._template_gauss_inp = 'template_gauss.inp'
    self._template_constraints = 'template.constraints'
    self._chInp = {1:'a', 2:'b', 3:'c', 4:'d'}

  def read_init(self, initFile):
    """Specify the objects from each star file to use"""

    init = np.genfromtxt(initFile, names=True, dtype=None, usecols=('component',
                                                                    'image',
                                                                    'number',
                                                                    'source_type'))
    for image in np.unique(init['image']):
      try:
        star_file = getattr(self, '{}_star_file'.format(image))
      except AttributeError:
        raise AttributeError('Image type {} does not exist. '\
                             'Please use hst, spitz, or aglow'.format(image))
      else:
        try:
          setattr(self, '{}_star_object'.format(image), sp.Star(star_file))
        except IOError:
          raise IOError('File {} does not exist. '\
                        'Please change star file for {}'.format(star_file,
                                                                image))
    return init

  def write_init_inp(self, init):
    """write the initial galfit inp file inp file. Send to galfit_writer"""
    with open(self._template_inp, 'r') as f:
      intro_raw = f.readlines()
    repls = ('GRB', self.GRB), ('CH', self.ch),
            ('FNUM', self._chInp[self.ch]),
            ('constraintFile', self.constraints_file)
    intro = reduce(lambda i, j: i.replace(*j), repls, intro_raw)

    for source in init:
      str_component = '# Component number: {}'.format(source['component'])
      str_template  = getattr(self,
                      '_template_{}_inp'.format(source['source_type']))
      with open(str_template, 'r') as f:
        body_raw = f.readlines()

      func = getattr(self, source['source_type'])
      star_data_raw = func(source['image'], source['number'])

      star_data = image1_image2(star_data_raw)

      repls = ('XX', star_data['X_IMAGE']),
              ('YY', star_data['Y_IMAGE']),
              ('FWHM', star_data['FWHM_IMAGE']),
              ('AXISRATIO', 1./star_data['ELONGATION']),
              ('PA', star_data['THETA_IMAGE']),
              ('SKIP', star_data['SKIP'])

      body = reduce(lambda i, j: i.replace(*j), repls, body_raw)
    
  def psf(self, image, number):
    star_object = getattr(self, '{}_star_object'.format(image))
    star_data = star_object.get_object(number, 'x_world', 'y_world',
                                               'erra_world')
    return star_data

  def gaussian(self, image, number):
    star_object = getattr(self, '{}_star_object'.format(image))
    star_data = star_object.get_object(number, 'x_world', 'y_world',
                                               'erra_world', 'fwhm_world',
                                               'theta_image')
    return star_data
