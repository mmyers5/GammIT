import numpy as np
import fast_phot as fap
import fits_plot as fip
import sys

"""Usage: python this.py inFile indStart indEnd tcol"""

sargv     = sys.argv
fileIn    = str(sargv[1])
indStart  = int(sargv[2])
indEnd    = int(sargv[3])
tcol      = str(sargv[4])

fileRead  = np.genfromtxt(fileIn, dtype=None, names=True)
index = fileRead['index']
stampsMake = fileRead[(index>=indStart) & (index<=indEnd)]
for entry in stampsMake:
  GRB = entry['GRB']
  ch = entry['ch']
  fileImg = '/scratch2/GRBH-MZSFR/Analysis/{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.\
            format(GRB=GRB, ch=ch)

  X,Y = fap.get_phys(fileImg, entry['ra'], entry['dec'])
  fip.make_stamp(GRB, ch, X, Y, entry['z'], entry['r'], entry['circol'], 
                 savefile=True, tcol=tcol)
