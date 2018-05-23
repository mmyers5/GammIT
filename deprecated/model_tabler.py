import numpy as np
import lum_dist as ld
from astropy.cosmology import WMAP9 as cosmo
import os
import calc_SFH as cs
import time

def find_delta(GRB,y,var):
  fZeros    = ['051022','051016B','060218','070721B','080928','100621A','110731A']
  var = var[~np.isnan(y)]
  y   = y[~np.isnan(y)]
  if GRB in fZeros:
      var = var[y!=0]
      y   = y[y!=0]
  chi = y
  minChi  = np.min(chi)  # minimum chi
  indMid  = np.where(chi==minChi)[0][0] # index of min
  indGood = np.where(chi<minChi+1)[0] # good indices in chi2, within delta

  lowGood   = indGood[indGood<indMid]
  highGood  = indGood[indGood>indMid]
  if len(lowGood)==1:
    lowGood=lowGood
  else:
    for i in np.arange(len(lowGood)-1,1,-1):
      thisInd = lowGood[i]
      nextInd = lowGood[i-1]
      if (chi[nextInd]<chi[thisInd]) or (nextInd!=thisInd-1):
       lowGood=lowGood[i:]
       break
  if len(highGood)==1:
    highGood=highGood
  else:
    for i in np.arange(0,len(highGood)-1,1):
      thisInd = highGood[i]
      nextInd = highGood[i+1]
      if (chi[nextInd]<chi[thisInd]) or (nextInd!=thisInd+1):
       highGood=highGood[0:i+1]
       break
  indGood = np.concatenate((lowGood,[indMid],highGood))
  #indLow  = np.min(indGood)
  #indHigh = np.max(indGood)

  lowVar  = np.min(var[indGood])
  highVar = np.max(var[indGood])

  #lowVar  = var[indLow]
  #highVar = var[indHigh]
  midVar = var[indMid]


  return midVar,highVar,lowVar,minChi

def make_bigger_table(outfile):
  '''
  tags: SSP_200grid_FZ,CSFH_200grid_FZ,chi2grid_200grid,chi2grid_2d,
  plot_delta,plot_delta_m2

  '''
  fZeros    = ['051022','051016B','060218','070721B','080928','100621A','110731A']
  tags      = ['SSP_200grid_FZ','CSFH_200grid_FZ','plot_delta','plot_delta_m2']
  grbFile   = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/grbs.txt'
  myGRBs    = np.genfromtxt(grbFile,names=True,dtype=None)
  outWrite  = open(outfile,'w')
  outWrite.write('GRB & z & model & SFR & Mass (Gsol) & Age (Myr) & Av & chi\n')
  parDict   = {}
  row       = 0
  # populate parDict
  for tagNum,tag in enumerate(tags):
    if (tagNum == 0) or (tagNum == 1):
      uncFile = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}/outmodelparunc.dat'.\
      format(tag=tag)
      models  = np.genfromtxt(uncFile,dtype=None,skiprows=2)
      GRBs = models['f0']
      z    = models['f2']
      SFR  = (models['f4'],models['f5'],models['f6'])
      mass = (models['f8']*1e9,models['f9']*1e9,models['f10']*1e9) # in sol
      age  = (models['f12']*1e6,models['f13']*1e6,models['f14']*1e6) # in yr
      Av   = (models['f16'],models['f17'],models['f18'])
      chi2 = models['f20']

	    # make pretty
      massval = np.log10(mass[0])
      massupp = np.log10(mass[0]+mass[1])-np.log10(mass[0])
      masslow = np.log10(mass[0])-np.log10(mass[0]-mass[2])
      mass    = (massval,massupp,masslow)
      ageval = np.log10(age[0])
      ageupp = np.log10(age[0]+age[1])-np.log10(age[0])
      agelow = np.log10(age[0])-np.log10(age[0]-age[2])
      age    = (ageval,ageupp,agelow)

      if tagNum == 1:
        parFile = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}/outmodelpar.dat'.format(tag=tag)
        dnt_mod = np.genfromtxt(parFile,dtype=None,skiprows=2)
        tau = dnt_mod['f5']
        unc = np.zeros(len(tau))
        newSFR = []
        for index,subSFR in enumerate(SFR[0]):
          if subSFR == 0: # use the upper-limit on SFR to calc age
            newSFR.append(SFR[1][index])
          else:
            newSFR.append(SFR[0][index])
        newSFR = np.array(newSFR)
        #age = (cs.calc_age2(mass[0]*1e9,newSFR)/1e6,unc,unc)  # in Myr
        age = (np.log10(cs.calc_age2(mass[0]*1e9,newSFR+0.0001)),unc,unc)  # in Myr
    elif (tagNum == 2) or (tagNum == 3):
      tagGRB  = []
      tagAge  = []
      tagMass = []
      tagAv   = []
      tagChi2 = []
      tagSFR  = []

      for GRB in myGRBs['GRB']:
        if tagNum == 2:
          parDir  = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}'.\
      format(tag=tag)
          gridFile  = '{parDir}/{GRB}.txt'.format(parDir=parDir,GRB=GRB)
          try:
            grid    = np.genfromtxt(gridFile,names=True,dtype=None)
          except:
            continue
          else:
            grid    = np.genfromtxt(gridFile,names=True,dtype=None)
          ages    = grid['AGE1']
          masses  = grid['MASS']
          Avs     = grid['AV']
          chi2s   = grid['CHI2']


          age   = find_delta(GRB,chi2s,ages)[0:3]
          mass  = find_delta(GRB,chi2s,masses)[0:3]
          Av    = find_delta(GRB,chi2s,Avs)[0:3]
          minChi= find_delta(GRB,chi2s,ages)[-1]

        elif  tagNum == 3:
          parDir  = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}'.\
      format(tag=tag)
          gridFile  = '{parDir}/{GRB}grid.txt'.format(parDir=parDir,GRB=GRB)
          try:
            grid    = np.genfromtxt(gridFile,names=True,dtype=None)
          except:
            continue
          else:
            grid    = np.genfromtxt(gridFile,names=True,dtype=None)
          ages    = grid['AGE1']
          masses  = grid['MASS']

          gridFile  = '{parDir}/{GRB}chi2.txt'.format(parDir=parDir,GRB=GRB)
          grid  = np.loadtxt(gridFile,skiprows=1)
          grid[grid==0] = 999
          grid  = np.nan_to_num(grid)
          grid[grid==0]= 777
          grid[grid==999]=0
          if GRB in fZeros:
            grid[grid==0] = 999
          minGrid = np.min(grid)
          rowGrid = np.where(grid==minGrid)[0][0] # grid mass, use for age
          colGrid = np.where(grid==minGrid)[1][0] # grid age, use for mass
          chi2s   = (grid[rowGrid,:],grid[:,colGrid])
          gridFile  = '{parDir}/{GRB}Av.txt'.format(parDir=parDir,GRB=GRB)
          grid  = np.loadtxt(gridFile,skiprows=1)
          Avs   = (grid[rowGrid,:],grid[:,colGrid])


          age   = find_delta(GRB,chi2s[0],ages)[0:3]
          mass  = find_delta(GRB,chi2s[1],masses)[0:3]
          Av0   = find_delta(GRB,chi2s[0],Avs[0])[0:3]
          Av1   = find_delta(GRB,chi2s[1],Avs[1])[0:3]
          minChi= find_delta(GRB,chi2s[1],ages)[-1]
          AvUp  = np.sqrt(Av0[1]**2+Av1[1]**2)
          AvLow = np.sqrt(Av0[2]**2+Av1[2]**2)
          Av    = (Av0[0],AvUp,AvLow)

        massval = np.log10(mass[0])
        massupp = np.log10(mass[1])-np.log10(mass[0])
        masslow = np.log10(mass[0])-np.log10(mass[2])
        masser  = (massval,massupp,masslow)
        ageval  = np.log10(age[0])
        ageupp  = np.log10(age[1])-np.log10(age[0])
        agelow  = np.log10(age[0])-np.log10(age[2])
        ageer     = (ageval,ageupp,agelow)
        SFR     = (0,0,0) 

        tagGRB.append(GRB)
        tagAge.append(ageer)
        tagMass.append(masser)
        tagAv.append(Av)
        tagChi2.append(minChi)
        tagSFR.append(SFR)
    
      GRBs  = np.array(tagGRB)
      SFR   = np.array(tagSFR).transpose()
      mass  = np.array(tagMass).transpose()
      age   = np.array(tagAge).transpose()
      Av    = np.array(tagAv).transpose()
      chi2  = np.array(tagChi2)

    modelDict = {'GRBs':GRBs,\
                 'SFR' :SFR,\
                 'mass':mass,\
                 'age' :age,\
                 'Av'  :Av,\
                 'chi2':chi2}
    parDict[tag] = modelDict


  for index,GRB in enumerate(myGRBs['GRB']):
    z = myGRBs['z'][index]
    for row, tag in enumerate(tags):
      if row==0:
        outWrite.write('{GRB} & {z} & {row} &'.format(GRB=GRB,z=z,row=row+1))
      else:
        outWrite.write('& & {row} & '.format(row=row+1))
      if GRB in parDict[tag]['GRBs']:
        idx = np.where(parDict[tag]['GRBs']==GRB)[0][0]
      else:
        outWrite.write('{...} & {...} & {...} & {...} & {...} \\\\\n')
        continue
      for par in ['SFR','mass','age','Av','chi2']:
        if par!='chi2':
          try:
            val   = parDict[tag][par][0][idx]
          except:
            print GRB,par
            val   = 0
            valup = 0
            vallow= 0
          else:
            val   = parDict[tag][par][0][idx]
            valup = parDict[tag][par][1][idx]
            vallow= parDict[tag][par][2][idx]
          if (val!=0) and (val!=np.inf):
            outStr='{val:1.2f}[$^[+{valup:.2f}]_[-{vallow:.2f}]$] & '.\
                   format(val=val,valup=valup,vallow=vallow)
            outStr=outStr.replace('[','{').replace(']','}')
          else:
            outStr='{...} & '
        else:
          val   = parDict[tag][par][idx]
          outStr='{val} \\\\\n'.format(val=val)
        outWrite.write(outStr)
  outWrite.close() 
  return parDict

def make_big_table_withfig(tags,outfile,outfig):
  grbFile = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/grbs.txt'
  myGRBs = np.genfromtxt(grbFile,names=True,dtype=None)
  row = 0
  outWrite = open(outfile,'w')
  outWriteFig = open(outfig,'w')
  outWrite.write('GRB & z & model & SFR & Mass (Gsol) & Age (Myr) & Av & chi\n')
  parDict = {}
  for tag in tags:
    uncFile = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}/outmodelparunc.dat'.format(tag=tag)
    models  = np.genfromtxt(uncFile,dtype=None,skiprows=2)
    GRBs = models['f0']
    
    z    = models['f2']
    SFR  = (models['f4'],models['f5'],models['f6'])
    mass = (models['f8']*1e9,models['f9']*1e9,models['f10']*1e9) # in sol
    age  = (models['f12']*1e6,models['f13']*1e6,models['f14']*1e6) # in yr
    Av   = (models['f16'],models['f17'],models['f18'])
    chi2 = models['f20']

    # make pretty
    massval = np.log10(mass[0])
    massupp = np.log10(mass[0]+mass[1])-np.log10(mass[0])
    masslow = np.log10(mass[0])-np.log10(mass[0]-mass[2])
    mass    = (massval,massupp,masslow)
    ageval = np.log10(age[0])
    ageupp = np.log10(age[0]+age[1])-np.log10(age[0])
    agelow = np.log10(age[0])-np.log10(age[0]-age[2])
    age    = (ageval,ageupp,agelow)

    #if tag == 'double_negative_tau':
    if tag == 'CSFH_200grid':
      parFile = '/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}/outmodelpar.dat'.format(tag=tag)
      dnt_mod = np.genfromtxt(parFile,dtype=None,skiprows=2)
      tau = dnt_mod['f5']
      unc = np.zeros(len(tau))
      newSFR = []
      for index,subSFR in enumerate(SFR[0]):
        if subSFR == 0:
          newSFR.append(SFR[1][index])
        else:
          newSFR.append(SFR[0][index])
      newSFR = np.array(newSFR)
      if tag == 'CSFH_200grid':
        #age = (cs.calc_age2(mass[0]*1e9,newSFR)/1e6,unc,unc)  # in Myr
        age = (np.log10(cs.calc_age2(mass[0]*1e9,newSFR)),unc,unc)  # in Myr
      else:
        age = (cs.calc_age(tau,mass[0]*1e9,newSFR)/1e6,unc,unc)  # in Myr
    modelDict = {'GRBs':GRBs,\
                 'SFR' :SFR,\
                 'mass':mass,\
                 'age' :age,\
                 'Av'  :Av,\
                 'chi2':chi2}
    parDict[tag] = modelDict
  for index,GRB in enumerate(myGRBs['GRB']):
    z = myGRBs['z'][index]
    for row,tag in enumerate(tags):
      figFile='/scratch2/GRBH-MZSFR/Analysis/galbuilds/{tag}/plots/{GRB}hostfit.eps'.\
      format(tag=tag,GRB=GRB)
      if os.path.isfile(figFile):
        texFile = 'Models/{tag}/{GRB}hostfit.eps'.format(tag=tag,GRB=GRB)
        texStr  = '\\begin<overpic>[width=0.5\\textwidth]<{texFile}>\n\
                   \\put (73,60) <{GRB}>\n\
                   \\end<overpic>\n'.\
                   format(tag=tag,GRB=GRB,texFile=texFile)
        texStr = texStr.replace('<','{').replace('>','}')
        if row % 2 == 0:
          texStr = texStr+'&\n'
        else:
          texStr = texStr+'\\\\ [-1.5ex]\n'
      else:
        texStr = 'Fitting failed.'
        if row % 2 == 0:
          texStr = texStr+'&\n'
        else:
          texStr = texStr+'\\\\\n'
      outWriteFig.write(texStr)

      if row==0:
        outWrite.write('{GRB} & {z} & {row} & '.format(GRB=GRB,z=z,row=row+1))
      else:
        outWrite.write('& & {row} & '.format(row=row+1))
      if GRB in parDict[tag]['GRBs']:
        idx = np.where(parDict[tag]['GRBs']==GRB)[0][0]
      else:
        outWrite.write('{...} & {...} & {...} & {...} & {...} \\\\\n')
        continue

      for par in ['SFR','mass','age','Av','chi2']:
        if par!='chi2':
          try:
            val   = parDict[tag][par][0][idx]
          except:
            print GRB,par
            val   = 0
            valup = 0
            vallow= 0
          else:
            val   = parDict[tag][par][0][idx]
            valup = parDict[tag][par][1][idx]
            vallow= parDict[tag][par][2][idx]
          if (val!=0) and (val!=np.inf):
            outStr='{val:1.2f}[$^[+{valup:.2f}]_[-{vallow:.2f}]$] & '.\
                   format(val=val,valup=valup,vallow=vallow)
            outStr=outStr.replace('[','{').replace(']','}')
          else:
            outStr='{...} & '
        else:
          val   = parDict[tag][par][idx]
          outStr='{val} \\\\\n'.format(val=val)
        outWrite.write(outStr)
  outWrite.close() 
  outWriteFig.close()

def append_to_table(outmodel,original_table,outfile):
  inTable = np.genfromtxt(original_table,names=True,dtype=None,delimiter='&')
  outWrite = open(outfile,'w')
  for index in range(len(inTable['GRB'])):
    GRB = inTable['GRB'][index].strip(' ')
    sed = make_mass_string(GRB,outmodel)
    mag  = inTable['mag'][index].strip(' ')
    flux = inTable['flux'][index].strip(' ')
    typical = inTable['typical_mass'][index].strip(' ')
    maximum = inTable['maximum_mass'][index].strip(' ')
    outWrite.write('{GRB} & {flux} & {mag} & {typical} & {maximum} & {sed}'.\
            format(GRB=GRB,flux=flux,mag=mag,typical=typical,maximum=maximum,sed=sed))
  outWrite.close()

def make_mass_string(GRB,outmodel):
  a,b  = get_logmass(GRB,outmodel)
  mass = str(round(a[0],2)).strip('0')
  up   = str(round(a[1],2)).strip('0')
  low  = str(round(a[2],2)).strip('0')
  age  = get_age(GRB,outmodel)
  mass_string = '{mass}[$_[-{low}]^[+{up}]$] & {age:.2f} & {chi}'.\
           format(mass=mass,low=low,up=up,age=age,chi=b)
  mass_string = mass_string.replace('[','{')
  mass_string = mass_string.replace(']','}')
  mass_string = mass_string+r'\\'+'\n'
  return mass_string

def make_mass_table(myGRBs,outmodel,outfile):
  if myGRBs == 'all':
    myGRBs = np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/galbuilds/grbs.txt',dtype=None)
  f = open(outfile,'w')
  for GRB in myGRBs:
    a,b  = get_logmass(GRB,outmodel)
    mass = str(round(a[0],2)).strip('0')
    up   = str(round(a[1],2)).strip('0')
    low  = str(round(a[2],2)).strip('0')
    age  = get_age(GRB,outmodel)
    mass_string = '{mass}[$_[-{low}]^[+{up}]$] & {age} & {chi}'.\
             format(mass=mass,low=low,up=up,age=age,chi=b)
    mass_string = mass_string.replace('[','{')
    mass_string = mass_string.replace(']','}')
    mass_string = mass_string+r'\\'+'\n'
    f.write(mass_string)
  f.close()
  
def make_table_stuff(myGRBs,outmodel):
  f = open('/scratch2/GRBH-MZSFR/Analysis/out_table_temp.txt','w')
  for GRB in myGRBs:
    print GRB
    a,b  = get_logmass(GRB,outmodel)
    mass = str(round(a[0],2)).strip('0')
    up   = str(round(a[1],2)).strip('0')
    low  = str(round(a[2],2)).strip('0')
    age  = get_age(GRB,outmodel)
    string = '{mass}[$_[-{low}]^[+{up}]$] & {age} & {chi}\n'.\
             format(mass=mass,low=low,up=up,age=age,chi=b)
    string = string.replace('[','{')
    string = string.replace(']','}')
    f.write(string)
  f.close()
    

def get_logmass(GRB,outmodel):
  # get model data from galbuilder
  inDat = np.genfromtxt(outmodel,names=True,dtype=None)
  myDat = inDat[inDat['GRB']==GRB][0]
  logmass      = np.log10(myDat['Mass_Gsol'])+np.log10(1e9)
  logmassupper = np.log10(myDat['Mass_Gsol']+myDat['Mass_Gsol_upper'])+\
                 np.log10(1e9)
  logmasslower = np.log10(myDat['Mass_Gsol']-myDat['Mass_Gsol_lower'])+\
                 np.log10(1e9)
  chi2  = myDat['chi2']
  return (logmass,logmassupper-logmass,logmass-logmasslower),chi2

def get_age(GRB,outmodel):
  # get model age from galbuilder
  inDat = np.genfromtxt(outmodel,names=True,dtype=None)
  myDat = inDat[inDat['GRB']==GRB][0]
  age    = (myDat['Age_Myr']+10)/1e3  # start of burst in Gyr
  return age

def get_z(GRB,outmodel):
  # get model z from galbuilder
  inDat = np.genfromtxt(outmodel,names=True,dtype=None)
  myDat = inDat[inDat['GRB']==GRB][0]
  z     = myDat['z']  # redshift
  return z

def get_kcorr(GRB,kcorrFile):
  # get kcorr data
  inDat = np.genfromtxt(kcorrFile,names=True,dtype=None)
  myDat = [i for i in inDat if i['ID']==GRB and i['source_type']=='this']
  z   = myDat[0][1]
  mab = myDat[0][2]
  Mab = myDat[0][4]
  mab_unc = myDat[0][3]
  Mab_unc = mab_unc+abs(-2/5*np.log10((1+z)*(1.63/3.6))*0.2)
  if mab_unc == 0:
    Mab = Mab+Mab_unc
    Mab_unc = 0
  return (mab,mab_unc),(Mab,Mab_unc)

def get_ml_mass(GRB,mlFile):
  # get mass-light data
  inDat = np.genfromtxt(mlFileFile,names=True,dtype=None)
  myDat = inDat[inDat['GRB']==GRB]
  massUL= myDat['Mass_UL_Z1']
  massUL_lower = myDat['Mass_UL_Z1_lower']
  mass_UL_upper= myDat['Mass_UL_Z1_upper']
