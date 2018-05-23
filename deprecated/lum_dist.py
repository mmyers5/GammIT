from astropy.cosmology import WMAP9 as cosmo
import numpy as np
import fast_phot as fap
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import LogLocator,AutoMinorLocator

def get_DL(GRB):
  '''
  Given a gamma-ray burst, return the luminosity distance in megaparsec.
  '''
  master = np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/master_file.txt',\
                         names=True,dtype=None,delimiter='\t')
  z = master['z'][master['GRB']==GRB][0]
  DL= cosmo.luminosity_distance(z)
  
  return DL

def rest_wavelength(z,ch=1):
  '''
  Given a redshift and a channel, return the rest wavelength in nm.
  ===INPUTS===
  z: float
    redshift
  ch: int
    the Spitzer channel of observation. Default is 1.
  ===OUTPUTS===
  rest_nm: quantity
    the rest wavelength in nanometers (see astropy.units)
  '''
  chDict = {1:3.6*u.um,2:4.5*u.um,3:5.8*u.um,4:8.0*u.um}
  filt   = chDict[ch]
  rest   = filt/(1+z)
  rest_nm= rest.to(u.nm)

  return rest_nm

def mab_to_Mab(mab,DL):
  '''
  Given apparent magnitude and luminosity distance, return absolute magnitude
  Mab_10pc = mab - 5*np.log10(DL_pc/10)
  return Mab_10pc
  ===INPUTS===
  mab: float
    apparent AB magnitude
  DL: quantity
    luminosity distance in Mpc or parsec
  ===OUTPUTS===
  Mab: float
    absolute AB magnitude
  '''
#  Mab = mab - 5*np.log10(DL.to(u.pc).value/10)
  Mab = mab - 5*np.log10(DL/10)
  return Mab

def kcorr(z,band='H',ch=1):
  '''
  Get k-correction to object at z given band and rest filter
  ===INPUTS===
  z: float
    redshift
  band: str
    the EM-band to which to k-correct
  ch: int
    the Spitzer channel
  ===OUTPUTS===
  kcorr: float
    the k-correction
  '''
  bandDict = {'H':1630*u.nm,'K':2190*u.nm}
  chDict   = {1:3.6*u.um,2:4.5*u.um,3:5.8*u.um,4:8.0*u.um}
  band_nm  = bandDict[band].value
  ch_nm    = chDict[ch].to(u.nm).value
  a = 2.5*(-2.3+2)*np.log10((1+z)*band_nm/ch_nm)
  b = 2.5*np.log10(1+z)
  kcorr = a - b

  return kcorr

def Mab_band(inFile,outFile,mych=1):
  bands   = {1:'H',2:'K'}
  inDat   = np.genfromtxt(inFile,names=True,dtype=None,delimiter='\t')
  outWrite= open(outFile,'w')
  outDict  = {}
  outMab   = {}
  outKcorr = {}
  outWrite.write('GRB\ttype\tMab_{band}\tk-correction\n'.\
                 format(band=bands[mych]))
  for idx,GRB in enumerate(inDat['GRB']):
    ch   = inDat['ch'][idx]
    fType= inDat['type'][idx]
    if ch==mych:
      z = inDat['z'][idx]
      DL  = cosmo.luminosity_distance(z)

      mab = inDat['mab'][idx]
      Mab = mab_to_Mab(mab,DL)
      kCorr = kcorr(z,band=bands[mych],ch=mych)

      Mab_band = Mab - kCorr
      outMab[GRB]  = Mab_band
      outKcorr[GRB]= kCorr
      
      outWrite.write('{GRB}\t{fType}\t{Mab}\t{kCorr}\n'.\
                     format(GRB=GRB,fType=fType,Mab=Mab_band,kCorr=kCorr))
  outWrite.close()
  outDict['Mab']  = outMab
  outDict['kCorr']= outKcorr
  print 'Mean Mab_{band}: {mean}'.format(band=bands[mych],mean=np.mean(outMab.values()))
  print 'Mean k-correction: {mean}'.format(mean=np.mean(outKcorr.values()))
  
  return outDict

def kcorr2(z,inBand,outBand):
  '''
  Get k-correction to object at z given band and rest filter
  ===INPUTS===
  z: float
    redshift
  band: str
    the EM-band to which to k-correct
  ch: int
    the Spitzer channel
  ===OUTPUTS===
  kcorr: float
    the k-correction
  '''
  a = 2.5*(-2.3+2)*np.log10((1+z)*outBand/inBand)
  b = 2.5*np.log10(1+z)
  kcorr = a - b

  return kcorr

def grb_hosts(inFile,outFile,outBand=1.630,inBand=3.6):
  bands = {1.63:'H',2.19:'K'}   # bands in micron
  inDat = np.genfromtxt(inFile,names=True,dtype=None)
  outWrite = open(outFile,'w')
  outWrite.write('GRB\tz\tMab_{band}\tMab_{band}_unc\tkcorr\tkcorr_unc\tsource_type\tsource\n'.\
                 format(band=bands[outBand]))
  prevGRB   = ['null']
  for idx,GRB in enumerate(inDat['GRB']):
    z  = inDat['z'][idx]
    DL = cosmo.luminosity_distance(z)
    mab = inDat['mab'][idx]
    mab_unc = inDat['mab_unc'][idx]
    s = inDat['source'][idx]
    stype = inDat['source_type'][idx]
    
    if GRB not in prevGRB: # requires "this" sample to be first in txt
      prevGRB.append(GRB)
    else:
      continue
    Mab = mab_to_Mab(mab,DL)
    kcorr = kcorr2(z,inBand,outBand)
      
    Mab_band = Mab - kcorr
    kcorr_unc = abs(0.2*(-2.5*np.log10((1+z)*outBand/inBand)))
    Mab_band_unc = mab_unc+kcorr_unc
    
    outWrite.write('{GRB}\t{z}\t{Mab_band}\t{Mab_band_unc}\t{kcorr}\t{kcorr_unc}\t{stype}\t{s}\n'.\
    format(GRB=GRB,z=z,Mab_band=Mab_band,Mab_band_unc=Mab_band_unc,kcorr=kcorr,kcorr_unc=kcorr_unc,stype=stype,s=s))
  outWrite.close()
  allDat = np.genfromtxt(outFile,names=True,dtype=None)
  return allDat
  
def k_plots(inFile,figFile):
  inDat = np.genfromtxt(inFile,names=True,dtype=None,delimiter='\t')
  fig_size = plt.rcParams["figure.figsize"]
  fig_size[0] = 12  # width, inches
  fig_size[1] = 9   # height, inches
  
  # parse the sun!
  myGRBs    = inDat['ID'][inDat['source_type']=='this']
  for idx,item in enumerate(inDat):
    if item['source_type']=='grb' and item['ID'] in myGRBs:
      inDat['source_type'][idx]='null'   # redundant hosts get type null

  # formatting plot
  plt.xlabel('Redshift')
  plt.ylabel('Absolute magnitude (AB)')
  plt.xlim(0,4)
  plt.ylim(-30,-10)
  plt.gca().set_xticks(np.arange(1,5,1))
  plt.gca().set_xticklabels(np.arange(1,5,1))
  plt.gca().xaxis.set_minor_locator(AutoMinorLocator(10))
  plt.gca().yaxis.set_minor_locator(AutoMinorLocator(10))
  plt.gca().invert_yaxis()
  plt.gca().tick_params('both',length=15,width=1,which='major')
  plt.gca().tick_params('both',length=5,width=1,which='minor')

  colors = {'this':'black','grb':'orange','fieldgal':'gray'}
  labels = {'this':'This Study','grb':'Other $\it{Spitzer}$ Studies of GRB Hosts','fieldgal':'Field Galaxies'}
  alphas = {'this':1.0,'grb':1.0,'fieldgal':0.2}
  edges  = {'this':'1.0','grb':'1.0','fieldgal':'0.2'}
  markers= {'this':16,'grb':10,'fieldgal':10}

  for stype in ['fieldgal','grb','this']:
    gals = inDat[inDat['source_type']==stype]

    # upper limits
    ulZ = [i['z'] for i in gals if i['mab_unc']==0 or i['mab_unc']==99]
    ulM = [i['Mab_H'] for i in gals if i['mab_unc']==0 or i['mab_unc']==99]
    plt.errorbar(ulZ,ulM,fmt='v',c=colors[stype],\
                 alpha=alphas[stype],ecolor=edges[stype],markersize=markers[stype])

    # detections
    detZ   = [i['z'] for i in gals if i['mab_unc']!=0 and i['mab_unc']!=99]
    detM   = [i['Mab_H'] for i in gals if i['mab_unc']!=0 and i['mab_unc']!=99]
    detErr = [i['mab_unc'] for i in gals if i['mab_unc']!=0 and i['mab_unc']!=99]
    plt.errorbar(detZ,detM,yerr=detErr,fmt='o',c=colors[stype],\
                 alpha=alphas[stype],ecolor=edges[stype],markersize=markers[stype],label=labels[stype])
    
  handles,labels = plt.gca().get_legend_handles_labels()
  handles = [h[0] for h in handles]
  plt.gca().legend(list(reversed(handles)),list(reversed(labels)),\
             loc='lower right',numpoints=1,scatterpoints=1,frameon=False,\
             prop={'size':16})
    
  # more formatting
  restwav = np.arange(1,4,1)
  newax = plt.gca().twiny()
  newax.set_xlim((0.72,3.6))
  newticks = (3.6/restwav)-1
  newax.set_xticks(newticks)
  newax.set_xticklabels(restwav.astype(int))
  minors = (3.6/np.arange(0.7,3.6,0.1)) - 1
  newax.set_xticks(minors,minor = True)
#  newax.xaxis.set_minor_locator(LogLocator(base=10,numticks=1,subs=np.arange(0.5,3.5,0.1)))
#  newax.xaxis.set_tick_params(which='minor')
#  newax.set_xscale('log')
  plt.gca().tick_params('both',length=15,width=1,which='major')
  plt.gca().tick_params('both',length=5,width=1,which='minor')
  newax.set_xlabel(r'Rest Wavelength $(\mu m)$')
  matplotlib.rcParams.update({'font.size': 18})
  #plt.show()
  plt.savefig(figFile,format='png',transparent=True,\
              pad_inches=0,dpi=plt.gcf().get_dpi())
