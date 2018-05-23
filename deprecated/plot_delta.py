import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import FormatStrFormatter
import model_tabler as mt
from scipy.interpolate import interp1d
#from time import sleep


def plot_delta_m1(inDir):
  grbs =  np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/galbuilds/grbs.txt',\
          names=True,dtype=None)
  grbs =  grbs['GRB']
  plt.figure() 
  fig_size    = plt.rcParams['figure.figsize']
  fig_size[0] = 12
  fig_size[1] = 9
  outFile = '{inDir}/big_table.txt'.format(inDir=inDir)
  outWrite = open(outFile,'w')
  matplotlib.rcParams.update({'font.size':24})

  falseZero = ['051022','051016B','060218','070721B','080928','100621A','110731A']
  for index,grb in enumerate(grbs):
    inFile  = '{inDir}/{grb}.txt'.format(inDir=inDir,grb=grb)
    myArray = np.genfromtxt(inFile,names=True,dtype=None)
    outPlot = '{inDir}/{grb}.eps'.format(inDir=inDir,grb=grb)

    x = myArray['AGE1']/1e9 # Gyr
    y = myArray['CHI2']
    plt.clf()
   
    plt.semilogx(x,y,'k',linewidth=4)
    
    #plt.semilogx(x,y,'r',linewidth=2)
    plt.xlabel(r'Age (Gyr)')
    plt.ylabel(r'$\chi^2_{red}$')

    midX,highX,lowX,minChi  = mt.find_delta(grb,y,x)

    plt.vlines([lowX,highX],[0,0],[max(y),max(y)],colors='r',linestyles='dashed',\
    linewidth=4)
    lowLim=lowX-0.01
    if lowLim<min(x): lowLim=min(x)
    highLim = highX+1
    if highLim>max(x):  highLim=max(x)
    plt.xlim(lowLim,highLim)
    plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    lowLim = minChi-1
    if lowLim<0: lowLim=0
    highLim=minChi+3
    if highLim>max(y):  highLim=max(y)
    plt.ylim(lowLim,highLim)
    plt.title(r', $\chi^2_{min}$='+str(minChi))
    plt.savefig(outPlot,format='eps',transparent=True,pad_inches=0,\
                dpi=plt.gcf().get_dpi())

    texFile = 'Models/Delta/{grb}.eps'.format(grb=grb)
    texStr  = '\\begin<overpic>[width=0.5\\textwidth]<{texFile}>\n\
               \\put (73,60) <{grb}>\n\
               \\end<overpic>\n'.\
               format(grb=grb,texFile=texFile)
 
    texStr = texStr.replace('<','{').replace('>','}')
    if index % 2 == 0:
      texStr = texStr+'&\n'
    else:
      texStr = texStr+'\\\\ [-1.5ex]\n'
    outWrite.write(texStr)
  outWrite.close()

def plot_delta_m2(inDir):
  grbs =  np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/galbuilds/grbs.txt',\
          names=True,dtype=None)
  grbs =  grbs['GRB']
  plt.figure() 
  fig_size    = plt.rcParams['figure.figsize']
  fig_size[0] = 12
  fig_size[1] = 9
  outFile = '{inDir}/big_table.txt'.format(inDir=inDir)
  outWrite = open(outFile,'w')
  matplotlib.rcParams.update({'font.size':24})

  falseZero = ['051022','051016B','060218','070721B','080928','100621A','110731A']
  for index,grb in enumerate(grbs):
    inFile  = '{inDir}/{grb}grid.txt'.format(inDir=inDir,grb=grb)
    gridArray = np.genfromtxt(inFile,names=True,dtype=None)
    inFile  = '{inDir}/{grb}chi2.txt'.format(inDir=inDir,grb=grb)
    z = np.loadtxt(inFile,skiprows=1)
    outPlot = '{inDir}/{grb}.eps'.format(inDir=inDir,grb=grb)

    x = gridArray['AGE1']/1e9 # Gyr
    y = gridArray['MASS']/1e9 #Gsol
    plt.clf()
    try:
      minChi2 = np.min(z[z!=0])
    except:
      print '1: '+grb
      continue
    else:
      minChi2 = np.min(z[z!=0])


    stepsize  = 0.1
    levels    = np.arange(minChi2,1.1+minChi2,stepsize)
    try:
      CS = plt.contour(x,y,z,levels=levels)
    except:
      print '2: '+grb
      continue
    else:
      CS = plt.contour(x,y,z,levels=levels)
      
    plt.clabel(CS)
    plt.xlabel(r'Age (Gyr)')
    plt.ylabel(r'Mass (Gsol)')

    plt.title(r', $\chi^2_{min}$='+str(minChi2))
    plt.savefig(outPlot,format='eps',transparent=True,pad_inches=0,\
                dpi=plt.gcf().get_dpi())

    texFile = 'Models/Delta/{grb}.eps'.format(grb=grb)
    texStr  = '\\begin<overpic>[width=0.5\\textwidth]<{texFile}>\n\
               \\put (73,60) <{grb}>\n\
               \\end<overpic>\n'.\
               format(grb=grb,texFile=texFile)
 
    texStr = texStr.replace('<','{').replace('>','}')
    if index % 2 == 0:
      texStr = texStr+'&\n'
    else:
      texStr = texStr+'\\\\ [-1.5ex]\n'
    outWrite.write(texStr)
  outWrite.close()
