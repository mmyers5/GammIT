import os, sys
import numpy as np
import star_parser as sp

def dec2sex(x, y):
  dh, dd  = 360/24., 1.
  dm, dam = dh/60., dd/60.
  ds, das = dm/60., dam/60.
  h, d = int(x/dh), int(y)
  m, arcm = int((x%dh) / dm), int((y%dd) / dam)
  s, arcs = int((x%dh) % dm / ds), int((y%dd) % dam / das)
  ra  = '{}:{}:{:.3f}'.format(h, m, s)
  if np.sign(d) == -1:
    pm = '-'
  else:
    pm = '+'
  dec = '{}{}:{}:{:.3f}'.format(pm, abs(d), arcm, arcs)
  return ra, dec

if len(sys.argv) == 4:
  reference_file = sys.argv[2]
  target_file = sys.argv[3]
  init_file   = sys.argv[1]  
else:
  print("Usage: python imgmatch_manual.py init reference target.")

reference = sp.Star(reference_file)
target = sp.Star(target_file)
init   = np.genfromtxt(init_file, names=True, usecols=('reference',
                                                       'target'))
#outfile = target_file.split('fits')[0] + 'ccmap'
outfile = 'test.ccmap'
ccmap = []
for i in init:
  objtgt = i['target']
  objref = i['reference']
  x_tgt, y_tgt = target.get_object(objtgt, 'x_image', 'y_image')
  x_ref, y_ref = reference.get_object(objref, 'x_world', 'y_world')
  ra, dec = dec2sex(x_ref, y_ref)  # string!

  objtgt = 'SEX_{:.0f}'.format(objtgt)
  objref = 'SEX_{:.0f}'.format(objref)
  ccmapstr = '{:.6f} {:.6f} {} {} {} {}'.format(x_tgt, y_tgt, ra, dec,
                                                objtgt, objref)
  ccmap.append(ccmapstr)

with open(outfile, 'w') as f:
  f.write('\n'.join(ccmap))
