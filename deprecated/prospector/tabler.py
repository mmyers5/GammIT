import numpy as np
from os import listdir
import model_plotato_2 as mplo
from copy import deepcopy as dc

grbs = np.genfromtxt('grblist.txt', names=None, dtype=None)
outfile = 'test.txt'
allRows = []

def string_thetas(thetas):
  stringList = []
  for i in range(thetas.shape[0]):
    prim  = thetas[i,0]
    upp   = thetas[i,1]
    low   = thetas[i,2]
    stringBasic = '${prim:.2f}^[+{upp:.2f}]_[-{low:.2f}]$'.format(prim=prim,
                                                                upp=upp, low=low)
    stringModified  = stringBasic.replace('[','{').replace(']','}')
    stringList.append(stringModified)
  
  return stringList

def alter_chain(obj):
  obj.old_chain = dc(obj.chain)
  # get logmass
  obj.chain[:,:,0]  = np.log10(obj.chain[:,:,0])
  # get logage
  obj.chain[:,:,2]  = np.log10(obj.chain[:,:,2]*1e9)
  
def get_last_mc(grb):

  try:
    grbFiles  = listdir( 'cluster_run/{grb}'.format(grb=grb) )
  except OSError:
    return None
  grbMcs  = sorted([f for f in grbFiles if f.endswith('_mcmc')])
  runMax  = int( ( grbMcs[-1].split('_')[1] ).strip('run') )
  for mc in grbMcs:
    nrun  = int( ( mc.split('_')[1] ).strip('run') )
    if nrun > runMax:
      runMax = nrun 

  grbLastMc = [f for f in grbMcs if 'run{num}'.format(num=runMax) in f][0]
  fullPath  = 'cluster_run/{grb}/{mc}'.format(grb=grb, mc=grbLastMc)
  
  return fullPath

for grb in grbs:
  mcPath = get_last_mc(grb)
  if mcPath is None:
    continue
  obj = mplo.Result(mcPath)
  alter_chain(obj)  

  thetas  = obj.quantilize()
  stringList  = string_thetas(thetas)
  stringRow   = ' & '.join(stringList)
  fullRow     = r'{grb} & {stringRow} \\'.format(grb=grb, stringRow=stringRow)
  allRows.append(fullRow)

with open(outfile, 'w') as f:
  f.write('\n'.join(allRows))
