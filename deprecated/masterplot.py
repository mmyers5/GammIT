import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp,mannwhitneyu
from scipy import interpolate
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D

def mass_light(inFile,modFile):
  '''
  Make plot of masses...
  '''
  table       = np.genfromtxt(inFile,names=True,dtype=None,comments='#')
  SED = np.genfromtxt(modFile,names=True,dtype=None,comments='#')
  fig_size    = plt.rcParams['figure.figsize']
  fig_size[0] = 12
  fig_size[1] = 9
  myObj = 0
  myObjUL = 0
  SEDObj = 0 
  SEDObjUL = 0
  # plot upper limits and lines from M_70 to M_Max 
  for item in table:
    GRB = item['GRB']
    z   = item['z']
    if item['Mass_ULZ1_upper']==0:
      # my data
      mm = 10**(item['Mass_UL_Z1'])
      plt.plot(z,mm,'v',color='0.8',markersize=12)
    else:
      # my data
      mm_up  = 10**(item['Mass_UL_Z1']+item['Mass_ULZ1_upper'])
      mm_low = 10**(item['Mass_Typical_Z1']-item['Mass_Typical_Z1_lower'])
      plt.errorbar(z,np.mean((mm_up,mm_low)),yerr=mm_up-np.mean((mm_up,mm_low)),lw=12,ecolor='0.8',\
                   capsize=None,capthick=None,markeredgecolor='k',markeredgewidth=1)
  for item in table:
    GRB = item['GRB']
    z   = item['z']
    if item['Mass_ULZ1_upper']==0:
      # SED data
      ms = SED['Mass_Gsol'][SED['GRB']==GRB][0]*1e9
      ms_up  = SED['Mass_Gsol_upper'][SED['GRB']==GRB][0]*1e9
      plt.plot(z,ms+ms_up,'rv',markersize=12)
    else:
      # SED data
      ms = SED['Mass_Gsol'][SED['GRB']==GRB][0]*1e9
      ms_up  = SED['Mass_Gsol_upper'][SED['GRB']==GRB][0]*1e9
      ms_low = SED['Mass_Gsol_lower'][SED['GRB']==GRB][0]*1e9
      plt.errorbar(z,ms,yerr=[(ms_low,ms_up)],fmt='ro',lw=2,capsize=10,capthick=2,markersize=9)

  # format plot

  SEDObj, = plt.plot((0,),(0,),marker='o',color='red',linestyle='none',markersize=9)
  SEDObjUL, = plt.plot((0,),(0,),marker='v',color='red',markersize=12,linestyle='none')
  myObj = mpatches.Patch(color='0.8')
  myObjUL, = plt.plot((0,),(0,),marker='v',color='0.8',markersize=12,linestyle='none')
  plt.ylim(1e6,1e12)
  plt.xlim(0,4)
  plt.xticks(np.arange(1,5,1))
  plt.yscale('log')
  plt.xlabel('Redshift')
  plt.ylabel(r'Stellar Mass ($M_\odot$)')
  plt.legend([myObj,SEDObj],['Mass Ranges','SED Fits'],\
             loc='lower right',numpoints=1,frameon=False)
#  plt.legend(loc='lower right',numpoints=1,scatterpoints=1)
  matplotlib.rcParams.update({'font.size': 18})
  plt.show()

def k_hist(inFile,wav):
  '''
  Makes a histogram/scatter of k correction stuff with previous
  studies. Relevant headers are 'GRB', 'z', and 'Mab_H' or 'Mab_K'
  for wav
  '''
  kcorrs     = np.genfromtxt(inFile,names=True,dtype=None,delimiter='\t')
  fig_size   = plt.rcParams["figure.figsize"]
  fig_size[0]= 12   # width in inches
  fig_size[1]= 9    # height in inches
  fieldGals = np.genfromtxt('field_gals.txt',names=True,dtype=None,comments='#')
  
  
  # parse the sun!
  for item in kcorrs:
    if item['GRB'] in np.unique(kcorrs['GRB']):
      ind = np.where(kcorrs['GRB']==item['GRB'])
      kcorrs = np.delete(kcorrs,ind[0][0:-1])

  theirs = kcorrs[np.where(kcorrs['Source']!='grbh')]
  mine   = kcorrs[np.where(kcorrs['Source']=='grbh')]

  # formatting plot
  plt.xlabel('Redshift')
  plt.ylabel(r'$M_{AB}$ (H-Band)')
  plt.xlim(0,4)
  plt.ylim(-10,-30)
  plt.gca().set_xticks(np.arange(1,5,1))
  plt.gca().invert_yaxis()

  # fields, upper limits
  ulZ   = fieldGals['z'][fieldGals['mab_unc']==0]
  ulM   = fieldGals['Mab_H'][fieldGals['mab_unc']==0]
  #plt.errorbar(ulZ,ulM,yerr=np.array([ulM*0,np.ones(len(ulM))*1.0]),uplims=True,\
  #             fmt='o',c='cyan',alpha=0.5,ecolor='0.50',markersize=8)
  plt.errorbar(ulZ,ulM,fmt='^',c='cyan',alpha=0.5,ecolor='0.50',markersize=8)

  # fields, detections
  detZ   = fieldGals['z'][fieldGals['mab_unc']!=0]
  detM   = fieldGals['Mab_H'][fieldGals['mab_unc']!=0]
  detErr = fieldGals['mab_unc'][fieldGals['mab_unc']!=0]
  plt.errorbar(detZ,detM,yerr=detErr,fmt='o',markersize=8,c='cyan',alpha=0.2,\
               ecolor='0.75',label='field galaxies')

  # previous, upper limits
  ulZ   = theirs['z'][theirs['mab_unc']==0]
  ulM   = theirs['Mab_H'][theirs['mab_unc']==0]
  plt.errorbar(ulZ,ulM,fmt='^',c='magenta',alpha=0.5,ecolor='0.50',markersize=8)

  # previous, detections
  detZ   = theirs['z'][theirs['mab_unc']!=0]
  detM   = theirs['Mab_H'][theirs['mab_unc']!=0]
  detErr = theirs['mab_unc'][theirs['mab_unc']!=0]
  plt.errorbar(detZ,detM,yerr=detErr,fmt='o',markersize=8,c='blue',alpha=0.5,\
               ecolor='0.75',label='previous studies')

  # mine, upper limits
  ulZ   = mine['z'][mine['mab_unc']==0]
  ulM   = mine[wav][mine['mab_unc']==0]
  plt.errorbar(ulZ,ulM,fmt='^',c='black',ecolor='0.50',markersize=10)

  # mine, detections
  detZ   = mine['z'][mine['mab_unc']!=0]
  detM   = mine['Mab_H'][mine['mab_unc']!=0]
  detErr = mine['mab_unc'][mine['mab_unc']!=0]
  plt.errorbar(detZ,detM,yerr=detErr,fmt='o',markersize=9,c='black',\
               ecolor='0.50',label='This Study')
  
  # more formatting
  restwav= 3.6/(1+np.arange(1,5,1))
  newax = plt.gca().twiny()
  newax.set_xlim(plt.gca().get_xlim())
  newax.set_xticks(np.arange(1,5,1))
  newax.set_xticklabels(restwav.astype(int))
  newax.set_xlabel(r'Rest Wavelength $(\mu m)$')
  plt.legend(loc='lower right',numpoints=1,scatterpoints=1)
  matplotlib.rcParams.update({'font.size': 18})
  plt.show()

def test_samps(inFile):
  kcorrs     = np.genfromtxt(inFile,names=True,dtype=None,comments='#')
  fieldGals = np.genfromtxt('field_gals.txt',names=True,dtype=None,comments='#')
  
  # parse the sun!
  for item in kcorrs:
    if item['GRB'] in np.unique(kcorrs['GRB']):
      ind = np.where(kcorrs['GRB']==item['GRB'])
      kcorrs = np.delete(kcorrs,ind[0][0:-1])

  theirs = kcorrs[np.where(kcorrs['Source']!='grbh')]
  theirs = theirs[np.where(theirs['mab_unc']!=0)]
  mine   = kcorrs[np.where(kcorrs['Source']=='grbh')]
  mine   = mine[np.where(mine['mab_unc']!=0)]
  gals   = fieldGals[np.where(fieldGals['mab_unc']!=0)]

  # previous studies
  D,p   = ks_2samp(theirs['Mab_H'],mine['Mab_H'])
  N     = len(mine['Mab_H'])
  print("PREVIOUS STUDIES\nD: {D}\tp: {p}\tN: {N}".\
         format(D=D,p=p,N=N))
  U,p   = mannwhitneyu(theirs['Mab_H'],mine['Mab_H'])
  return U
  print("U: {U}\tp: {p}".format(U=U,p=2*p)) 
  # ks-test, field galaxies
  D,p   = ks_2samp(gals['Mab_H'],mine['Mab_H'])
  print("FIELD GALAXY STUDIES\nD: {D}\tp: {p}\tN: {N}".\
         format(D=D,p=p,N=N))
  # ks-test, boths
  D, p  = ks_2samp(np.concatenate((theirs['Mab_H'],fieldGals['Mab_H'])),\
                   mine['Mab_H'])
  print("ALL THE STUDIES\nD: {D}\tp: {p}\tN: {N}".\
         format(D=D,p=p,N=N))

def z_wav_histo(inFile):
  '''
  Makes a histogram of redshifts and rest wavelengths
  ---INPUTS---
  inFile: str
    the text file with redshift and channel as columns
    labeled as 'z' and 'ch'
  '''
  zWav = np.genfromtxt(inFile,names=True,dtype=None,comments='#')
  restWav = zWav['ch']/(zWav['z']+1)

  fig_size = plt.rcParams["figure.figsize"]
  fig_size[0]=12   # width in inches
  fig_size[1]=9    # height in inches
  plt.subplot(121)
  matplotlib.rcParams.update({'font.size':16})
  n,bins,patches=plt.hist(zWav['z'],bins=14,range=(0,4),histtype='stepfilled',\
    color='white',edgecolor='black',hatch='\\\\\\\\')
  plt.xticks(np.arange(1.0,5.0,1.0))
  plt.yticks(np.arange(2.0,13.0,2.0))
  plt.ylim(0,10.1)
  plt.xlim(0,4)
  plt.xlabel('Redshift',fontsize=16)
  plt.ylabel('Number of Targets Observed',fontsize=16)

  plt.subplot(122)
  n,bins,patches=plt.hist(restWav,bins=14,range=(0,3.5),\
    histtype='stepfilled',color='white',edgecolor='black',hatch='\\\\\\\\')
  plt.xticks(np.arange(1.0,5.0,0.5))
  plt.yticks(np.arange(2.0,13.0,2.0))
  plt.ylim(0,12.1)
  plt.xlim(0.5,3.6)
  plt.xlabel(r'Rest Wavelength $(\mu m)$',fontsize=16)
  plt.ylabel('Number of Targets Observed',fontsize=16)
  plt.show()
