import sys, os, glob
from astropy.io import fits as pyfits
from copy import deepcopy as dc

'''Usage: python list_creator.py [grb] [ch] [outdir]
e.g. python list_creator.py 130606A ch1 cdf/lists
or
e.g. python list_creator.py 130606A ch1 cdf/lists rblah
'''
grb = str(sys.argv[1])
ch  = str(sys.argv[2])
outdir  = str(sys.argv[3])

try:
  rarg  = str(sys.argv[4])
except IndexError:
  pre = '/scratch2/GRBH-MZSFR/Mosaicking/{grb}/r*/{ch}'.format(grb=grb, ch=ch)
else:
  pre = '/scratch2/GRBH-MZSFR/Mosaicking/{grb}/{rarg}/{ch}'.format(grb=grb, rarg=rarg, ch=ch)

fitsList  = glob.glob('{pre}/bcd/*{i}*.fits'.format(pre=pre, i='cbcd'))
approvedList  = []
for fitsName in fitsList:
  hdulist = pyfits.open(fitsName)
  exptime = hdulist[0].header['EXPTIME']
  framtime  = hdulist[0].header['FRAMTIME']
  if exptime/framtime >= 0.86:
    approvedList.append(fitsName)

for fitsType in ['cbcd', 'cbunc', 'bimsk']:
  outFile = '{outdir}/{grb}_{ch}_{i}list.txt'.format(outdir=outdir, i=fitsType,
                                                     grb=grb, ch=ch)
  if fitsType == 'cbcd':
    modifiedList = dc(approvedList)
  else:
    modifiedList = [i.replace('cbcd', fitsType) for i in approvedList]

  with open(outFile, 'w') as f:
    f.write('\n'.join(sorted(modifiedList)))
