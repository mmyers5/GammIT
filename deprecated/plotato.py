import numpy as np
import matplotlib.pyplot as plt

def plotSED(fMod,fRes,fDat,title):
  '''
  Plot SEDs. Will need to adjust xlim and ylim manually.
  ===
  INPUTS
  ===
  fMod: str
    filename that has the model data, column0 is wlength in angstroms, column1
    is flux density per wavelength
  fDat: str
    filename that has the filter res stuff. same structure as above
  title: str
    the title for the plot
  '''
  mod=np.loadtxt(fMod,skiprows=1)    # read columns
  mLength=mod[:,0]                   # wavelength in Ang
  mFlx=mod[:,1]                      # flux density in erg/s/cm2/A
  mFlx=lamb2Nu(mLength,mFlx)
  
  dat=np.loadtxt(fDat,skiprows=1)
  dLength=dat[:,0]
  dFlx=dat[:,1]
  derr = dat[:,2]
  #dFlx=lamb2Nu(dLength,dFlx)
  
  res=np.loadtxt(fRes,skiprows=1)
  rLength=res[:,0]
  rFlx=res[:,1]
  rFlx=lamb2Nu(rLength,rFlx)

  plt.plot(mLength/10,mFlx)          # plot nanometers vs flux
  plt.errorbar(dLength/10,dFlx,yerr=derr,fmt='o')
  plt.plot(rLength/10,rFlx,'o')

  plt.yscale('log')
  plt.xlim( (150,5000) )
  plt.ylim( (1,210) )
  plt.xlabel('Wavelength (nm)')
  plt.ylabel(r'$F_\nu (\mu Jy)$')
  plt.title(title)

  plt.show()

def lamb2Nu(lambLength,lambFlx):
  c = np.double(2.998e18)
  nu = lambFlx*1e10*(lambLength**2)/c
  return nu

def lamb2Nu2(fLamb,fOut):
  '''
  Convert uJy to flux density per frequency
  ===
  INPUTS
  ===
  fLamb: str
    file with flux density in uJy. column0 is wavelength in angstroms, column1
    is flux density in uJy
  fOut: str
    output file
  '''
  c = np.double(2.998e18)          # speed of light in Ang/s
  dat = np.loadtxt(fLamb, skiprows=1)
  dLength = dat[:,0]    # wavelengths in Ang
  dFlx = dat[:,1]       # flux density in 10^-19 erg/s/cm^2/Angstrom

  fFlx = dFlx*1e10*(dLength**2)/c  # flux density in 10^-29 ergs/s/cm^2/Hz

  f = open(fOut,'w')
  i=0
  f.write('# wl fl (x 10^-29 ergs s^-1 cm^-2 Hz^-1)\n')
  for item in fFlx:
    f.write('%f %f\n' % (dLength[i],item))
    i+=1
  f.close()

def ABMag(F1):
  m1 = 23.93
  m2 = 25
  F2 = 10**((m1-m2)/(-2.5))*F1
  return F2

def writeMe(fIn,fOut):
  dat = np.loadtxt(fIn,skiprows=1)
  dLength = dat[:,0]
  dFlx = dat[:,1]
  dErr = dat[:,2]

  dFlx = ABMag(dFlx)
  dErr = ABMag(dErr)

  f = open(fOut, 'w')
  i = 0
  f.write('# wl fl (x 10^-19 ergs s^-1 cm^-2 Ang^-1, AB-Mag 0-point = 25)\n')
  for flux in dFlx:
    f.write('%f %f %f\n' % (dLength[i],flux,dErr[i]))
    i+=1
  f.close()
  
  
