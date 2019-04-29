import fast_phot as fap
import numpy as np
import sys
import os

if len(sys.argv) == 3:
  countsFile = sys.argv[1]
  photFile   = sys.argv[2]
else:
  sys.exit('Usage: calc_all_phot.py [input] [output]')

inCounts  = np.genfromtxt(countsFile,names=True,dtype=None)
photWrite = open(photFile, 'w')
photWrite.write('GRB\tch\ttype\tflux\tflux_unc\tmab\tmab_unc\tap\tsig\n')

for index in np.arange(len(inCounts['GRB'])):
  GRB = inCounts['GRB'][index]
  #GRB = '0'+str(inCounts['GRB'][index])
  correction   = inCounts['corrCorr'][index]
  bkgSubCounts = inCounts['flxSubCnts'][index]
  srcUnc = inCounts['uncCnts'][index]
  srcPx  = inCounts['apPix'][index]
  bkgUnc = inCounts['uncBkg'][index]
  bkgPx  = inCounts['anPix'][index]
  ch = inCounts['ch'][index]
  ftype = inCounts['type'][index]
  ap    = inCounts['ap'][index]
  sig   = inCounts['sig'][index]

  if inCounts['type'][index] == 'unc':
    correction = inCounts['apCorr'][index]
    flx = fap.calc_phot(0.4,correction,0,srcUnc,srcPx,bkgUnc,bkgPx)
    mab = fap.uJy2AB_unc(flx[1]*3,0)
    flx = np.asarray(flx)
    flx[1]*=3
  elif inCounts['type'][index] == 'flx' or inCounts['type'][index] == 'sub':
    flx = fap.calc_phot(0.4,correction,float(bkgSubCounts),srcUnc,srcPx,bkgUnc,bkgPx)
    mab = fap.uJy2AB_unc(flx[0],flx[1])
  else:
    sys.exit('What happened? Type not recognized')
  
  photWrite.write('{GRB}\t{ch}\t{ftype}\t{flx}\t{flxunc}\t{m}\t{munc}\t{ap}\t{sig}\n'.\
           format(GRB=GRB,ch=ch,ftype=ftype,flx=flx[0],flxunc=flx[1],\
                  m=mab[0],munc=mab[1],ap=ap,sig=sig))
photWrite.close()
