from galfit_field import GalfitField, GalfitSource
import numpy as np
import sys

if len(sys.argv) == 4:
  GRB = str(sys.argv[1])
  ch  = int(sys.argv[2])
  init_file = str(sys.argv[3])
else:
  sys.exit('Usage: python galfit_master.py [GRB] [ch] [init_file]')

field = GalfitField(GRB, ch, init_file)
intro = write_init_inp(GRB, ch) 
for component in field.init:
  star_object = getattr(field, '{}_star_object'.format(component['image']))
  star = GalfitSource(star_object, component['number'], component['source_type')
  star.transform_coords(field.spitz_img)
  # assuming hst tied to mosaic and mosaic tied to afterglow
  if component['image'] == 'aglow' or component['image'] == 'spitz':
    tied_img = field.spitz_img
  elif component['image'] == 'hst':
    tied_img = field.hst_img
  star.tie_unc(tied_img)
  # available keys, star.target_xy_err, star.target_x, star.target_y,
  # star.target_fwhm, star.target_theta, star.rms

def write_source(compnum, template, source_type, xx, yy,):
  str_comp = '# Component number: {}'.format(compnum)
  with open(template, 'r') as f:
    body_raw = f.readlines()
  repls = ('XX', 

def write_init_inp(GRB, ch, cofile=None):
  chstyle = {1:'a', 2:'b', 3:'c', 4:'d'}
  with open('template.inp', 'r' as f):
    intro_raw = f.readlines() 
  if cofile == None:
    cofile = '{GRB}/{GRB}_ch{ch}.constraints'.format(GRB=GRB, ch=ch)
  repls = ('GRB', GRB), ('CH', ch), ('FNUM', chstyle[ch]),
          ('constraintFile', cofile)
  intro = reduce(lambda i, j: i.replace(*j), repls, intro_raw)
  return intro
