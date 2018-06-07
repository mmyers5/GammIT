import fits_plot as fip # plotting routines
import fast_phot as fap # fast photometry
import sigma_clipper as sc # sigma-clipping
import region_tools as rt
from collections import OrderedDict
from astropy.io import fits
import numpy as np
import sys, os

def gen_unc(inFile):
  '''
  Takes incoming .fits file and squares every pixel
  then multiplies by 1000 (for purpose of using funcnts with .fits
  file with very small values).
  ===INPUTS===
  inFile: str
    .fits file location
  ===OUTPUTS===
  Generates a file called fname_sq.fits.
  '''
  f = fits.open(inFile)
  data   = f[0].data                      # image array
  header = f[0].header                    # image header
  f.close()
  newData=(data**2)*1000                  # square the array and multiply by 1k

  outFile=(inFile.split('.fits'))[0]+'_sq.fits'     # generate new filename
  fits.writeto(outFile,newData,header,clobber=True) # write out filename

def region_coords(imgFile,regFile,inCoord='deg'):
  '''
  Take a ds9 region file of circles, return image coordinates
  ===INPUTS===
  imgFile: str
    image .fits file to which region coordinates belong
  regFile: str
    ds9 region file (currently must be in degree coordinates)
  inCoord: str
    currently nothing. will implement converting from different coord systems
  ===OUTPUTS===
  coords: dict
    a dictionary keyed by integers starting from 0.
    corresponding values are image x and y coordinates.
  '''
  regions = np.asarray(open(regFile,'r').readlines())        # read in regions
  coords  = {}
  i=0
  for line in regions:
    # assumes every relevant region is a plain circle
    if line.startswith('circle'):
      temp = open('temp.txt','w')                             # garbage file
      deg  = line[line.find('(')+1:line.find(')')].split(',') # get coords in deg
      temp.write('{ra}\t{dec}'.format(ra=deg[0],dec=deg[1]))
      temp.close()

      stdout = Popen('funsky {imgFile} {temp}'.\
               format(imgFile=imgFile,temp='temp.txt'),\
               shell=True,stdout=PIPE).stdout # write funsky output to stdout
      output = stdout.read()
      os.remove('temp.txt')
      coords[str(i)]=re.split(' +',output)[1:3]
      i+=1
  return coords

def get_file(GRB,ch,fType,aperture=None):
  '''
  Generate a file name given a GRB and channel, filetype and apertures.
  Currently assumes format like ./{GRB}/{GRB}_ch{ch}_{fType}
  ===INPUTS===
  GRB: str
  ch: int
  fType: str
    special file types include unc, sub, psf. By default, will assume 
    file intended is a "regular" wcs-registered mosaic
  aperture: str
    will accept '2-2-6' or '3-3-7' for purpose of files that are sigma-clipped
  ===OUTPUTS===
  fName: str
    file name used given the above. Does not check if file exists
  '''
  # prefix assumes files are located in subdirectory of current director
  fileDict = {'flx':'maic.wcs.fits','sub':'sub.wcs.fits','unc':'munc_sq.wcs.fits',\
              'psf':'psf.fits','modreg':'modelstamp.reg','modinp':'model.inp'}
  pref = '{GRB}/{GRB}_ch{ch}'.format(GRB=GRB,ch=ch)

  # make sure input matches known filetypes
  if fType not in fileDict.keys():
    return 'Invalid filetype, please use flx,sub,unc,psf,or modreg'

  # for use with files that have been sigma-clipped
  if aperture == None:
    fName = '{pref}_{tag}'.format(pref=pref,tag=fileDict[fType])
  else:
    fName = '{pref}_{ap}_{tag}'.\
             format(pref=pref,ap=aperture,tag=fileDict[fType])

  # make sure uncertainty square files exist
  if not os.path.isfile('{pref}_munc_sq.fits'.format(pref=pref)):
    gen_unc('{pref}_munc.fits'.format(pref=pref))
    
  # maybe file doesn't exist if .wcs isn't tied
  if not os.path.isfile(fName):
    registered   = '{pref}_maic.wcs.fits'.format(pref=pref)
    unregistered = ''.join('{pref}_{tag}'.format(pref=pref,tag=fileDict[fType]).split('.wcs'))
    copytie(registered,unregistered,fName)

  if not os.path.isfile('{pref}_munc_sq.wcs.fits'.format(pref=pref)):
    registered   = '{pref}_maic.wcs.fits'.format(pref=pref)
    copytie(registered, '{pref}_munc_sq.fits'.format(pref=pref),
                        '{pref}_munc_sq.wcs.fits'.format(pref=pref))

  return fName

def copytie(original,subbed,output):
  '''
  ===INPUTS===
  copy wcs registration to image
  original: str
    filename with wcs coordinates to be copied over
  subbed: str
    filename with modeled regions (default has no wcs coords)
  output: str
    the output filename, subbed's data with original's header
  '''

  originalHeader = fits.open(original)[0].header
  subbedData     = fits.open(subbed)[0].data
  fits.writeto(output,subbedData,originalHeader)


'''
Given a table with appropriate headers, will calculate photometry for each row
of data and generate appropriate plots
USAGE: master_phot.py [input] [output] [startnum]
===INPUTS===
input: str
  the table; headers are (in order): 
  index: int
    an index, or row number
  GRB: str
    grb name
  z: float
    redshift
  type: str
    indicates "flx" for host galaxies that don't need any additional treatment,
    "unc" for non-detections, and "sub" for host galaxies that underwent
    modeling
  ch:
    channel of Spitzer image
  ra: float
    right ascension of host in degrees, 
  dec: float
    declination of host in degrees
  circol: str
    a string to indicate color for stamps, by default will be "r" for red
output: str
  the filename to which output photometry will be written. Headers are rather
  annoying to outline.
startnum: int
  index the entries in input, use index here to start at that index
'''

# ensure proper usage
if len(sys.argv) == 5:
  masterFile  = sys.argv[1]
  outFile     = sys.argv[2]
  startWith   = sys.argv[3]
  sig         = float(sys.argv[4])  # for sigma clipping
  if  not os.path.isfile(masterFile):
    sys.exit('{masterFile} does not exist'.format(masterFile=masterFile))
else:
  sys.exit('Usage: master_phot.py [input] [output] [start index] [sigma]')

# write out the header of output table
seq = ('index','GRB','ch','type','apCorr','corrCorr','apPix','anPix',\
       'flxCnts','flxBkg','flxSubCnts','uncCnts','uncBkg','numClip',\
       'opsfCnts','opsfBkg','opsfSubCnts','psfCnts','psfBkg','psfSubCnts',\
       'opsfapPix','opsfanPix','ap','sig')
if not os.path.isfile(outFile):
  outWrite = open(outFile,'a')
  outWrite.write('{header}\n'.format(header=' '.join(seq)))
else:
  outWrite = open(outFile,'a')
inData   = np.genfromtxt(masterFile,names=True,dtype=None)

# to help interpret aperture. values are in img coordinates
# using 0.4"/pixel
apDict = {'2-2-6':(2.4,7.2),'3-3-7':(3.6,8.4)}

# make sure indices are always consecutively labeled
# truly silly way to index data, FIX  THIS
for idx,i in enumerate(inData['index']): 
  if int(i)<int(startWith):
    continue 

  # grab preliminary data for each row
  GRB = inData['GRB'][idx]
  z   = inData['z'][idx]
  ch  = inData['ch'][idx]
  r   = inData['r'][idx]
  fType  = inData['type'][idx]  
  raDeg  = inData['ra'][idx]
  decDeg = inData['dec'][idx]
  circol = inData['circol'][idx]

  for ap in apDict.keys():
    # fill dictionary with preliminary data for each aperture
    raw = OrderedDict.fromkeys(seq)
    raw['GRB'] = GRB
    raw['ch']  = ch
    raw['type']= fType
    raw['ap']  = ap
    raw['index']  = i
    raw['apCorr'] = fap.get_acorr(ap,ch)
    # easiest case, if the method is only to grab upper limit
    if fType == 'unc':
      raw['corrCorr'] = raw['apCorr']             # no aperture corr needed
      uncFile = get_file(GRB,ch,fType=fType)
      XY = fap.get_phys(uncFile,raDeg,decDeg)     # spitzer physical centroid
      fap.get_cnts(uncFile,XY,ap,'cntsTemp.txt')  # counts output to cntsTemp.txt
      raw['uncCnts'],raw['apPix'],raw['uncBkg'],raw['anPix'],raw['numClip']=fap.read_cnts('cntsTemp.txt')

    # not-so-straightforward case of regular flux or subbed calculations
    # requires psf and k-sigma corrections
    # grab appropriate files 
    elif fType == 'flx' or fType == 'sub':
      files = OrderedDict()                         # files to be used
      files['opsf'] = get_file(GRB,ch,fType='psf')  # non-k-sigma corrected psf

      # grab unc, source flux, and psf files for k-sigma corrections
      for prefix in ['unc','flx','psf']:
        if prefix=='flx' and fType =='sub':         # modeled files adjust filename
          files[prefix]= get_file(GRB,ch,fType='sub',aperture=ap)
        else:
          files[prefix] = get_file(GRB,ch,fType=prefix,aperture=ap)

      # perform counts on unc, source flux, and psf image files
      for img in files:

        # grab physical and image coordinate centroids
        if img=='psf' or img=='opsf':
          XY = (65,65)   # assumes PSFs that are sized 128 by starfinder
          xy = XY
        elif fType=='sub':
          inpFile = get_file(GRB,ch,fType='modinp',aperture=None)
          xy = tuple(rt.get_coords_from_inp(inpFile)[0])
          LTV1,LTV2=fap.get_phys(files['unc'],raDeg,decDeg,phys=False,LTV=True)
          XY = (xy[0] - LTV1,xy[1]-LTV2)
        else:
          XY = fap.get_phys(files['unc'],raDeg,decDeg)
          xy = fap.get_phys(files['unc'],raDeg,decDeg,phys=False)

        
        # sigma_clipper the images
        if img=='flx' and fType=='sub':
          sig=sc.pipe_mask(GRB,ch,int(xy[0]),int(xy[1]),ap,sub=True,sigHigh=sig,sigLow=sig)
        elif img=='flx' and fType!='sub':
          sig=sc.pipe_mask(GRB,ch,int(xy[0]),int(xy[1]),ap,sigHigh=sig,sigLow=sig)
        raw['sig']=sig
        # get counts for unc, source flux, and psf files
        fap.get_cnts(files[img],XY,ap,'cntsTemp.txt')
        if img=='psf' or img == 'opsf':
          raw[img+'Cnts'],raw[img+'apPix'],raw[img+'Bkg'],raw[img+'anPix'],raw['numClip']=fap.read_cnts('cntsTemp.txt')
        else:
          raw[img+'Cnts'],raw['apPix'],raw[img+'Bkg'],raw['anPix'],raw['numClip']=fap.read_cnts('cntsTemp.txt')
      
      # background-subtract the counts in source flux and psf files
      for img in ['flx','opsf','psf']:
        if img=='psf':
          raw[img+'SubCnts']=fap.sub_cnts(raw[img+'Cnts'],raw['opsfapPix'],raw[img+'Bkg'],raw['opsfanPix'])
        elif img =='opsf':
          raw[img+'SubCnts']=fap.sub_cnts(raw[img+'Cnts'],raw['opsfapPix'],raw[img+'Bkg'],raw['opsfanPix']-raw['numClip'])
        else:
          raw[img+'SubCnts']=fap.sub_cnts(raw[img+'Cnts'],raw['apPix'],raw[img+'Bkg'],raw['anPix']-raw['numClip'])

      # calculated corrected aperture correction
      raw['corrCorr']=raw['apCorr']*raw['opsfSubCnts']/raw['psfSubCnts']

    # delete non-essential keys
    if 'psfanPix' in raw:
      del raw['psfanPix']
      del raw['psfapPix']

    # write counts data out to file
    for header in raw:
      outWrite.write('{header}\t'.format(header=raw[header]))
    outWrite.write('\n')
    
    # PLOTTING PROCEDURES
    # make relevant stamps
    #XY = fap.get_phys('{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.\
    #                  format(GRB=GRB,ch=ch),raDeg,decDeg)
    #fip.make_stamp(GRB,ch,XY[0],XY[1],z,r,circol,tcol='white')
    #if fType == 'sub':
    #  fip.make_model_stamp(GRB,ch,xy,circol,z)
    #  fip.make_model_stamp(GRB,ch,xy,circol,z,ap=(2.4,7.2))
    #  fip.make_model_stamp(GRB,ch,xy,circol,z,ap=(3.6,8.4))
    
outWrite.close()
