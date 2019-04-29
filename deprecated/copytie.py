import os,sys
from astropy.io import fits

original = str(sys.argv[1])
subbed   = str(sys.argv[2])
output   = str(sys.argv[3])

originalHeader = fits.open(original)[0].header
subbedData     = fits.open(subbed)[0].data

try:
  fits.writeto(output,subbedData,originalHeader)
except:
  os.remove(output)
  fits.writeto(output,subbedData,originalHeader)
