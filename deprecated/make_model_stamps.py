import numpy as np
import fast_phot as fap
import fits_plot as fip
import sys

"""Usage: python this.py inFile indStart indEnd"""

sargv     = sys.argv
fileIn    = str(sargv[1])
indStart  = int(sargv[2])
indEnd    = int(sargv[3])


fileRead  = np.genfromtxt(fileIn, dtype=None, names=True)
index = fileRead['index']
subs  = fileRead['type']
stampsMake = fileRead[(index>=indStart) & (index<=indEnd) & (subs=='sub')]
for entry in stampsMake:
  fip.make_model_stamp(entry['GRB'], entry['ch'], (entry['ra'], entry['dec']),
                       entry['circol'], entry['z'], 
                       ap=(2.4,7.2), savefile=True)
  fip.make_model_stamp(entry['GRB'], entry['ch'], (entry['ra'], entry['dec']),
                       entry['circol'], entry['z'], 
                       ap=(3.6,8.4), savefile=True)
