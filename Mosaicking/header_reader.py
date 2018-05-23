'''Will delete filenames in the list that have insufficient exptime'''

from astropy.io import fits as pyfits
import os

fList=['cbunclist.txt','bimsklist.txt']                         # these lists need to be generated first
os.system('ls *cbcd.fits > cbcdlist.txt')                                      # only works in linux OS
os.system('ls *cbunc.fits > cbunclist.txt')
os.system('ls *bimsk.fits > bimsklist.txt')

f = open('cbcdlist.txt','r')
allFiles=f.readlines()
f.close()
f= open('cbcdlist.txt','r')
allFilesCop=f.readlines()
indices=[]
i=0
for fname in allFilesCop:
    hdulist=pyfits.open(fname.strip('\n'))
    exptime=hdulist[0].header['EXPTIME']
    framtime=hdulist[0].header['FRAMTIME']
    if i==0: 
      init=exptime/framtime
      allFiles.pop(0)
      indices.append(0)
      print '%i: initial, %f %f' %(i,exptime,framtime) 
      i+=1
      continue
    if (exptime/framtime)<=0.86:
    #if (exptime/framtime)<=.90:
      index = allFiles.index(fname)
      allFiles.pop(index)
      indices.append(index)
      print '%i: spit, %f %f' %(i,exptime,framtime) 
      i+=1
      continue
    print '%i: accept, %f %f' %(i,exptime,framtime) 
    i+=1
f=open('cbcdlist.txt','w')
f.writelines(allFiles)
f.close()

for aux in fList:
    f = open(aux,'r')
    allFiles=f.readlines()
    f.close()
    for index in indices:
        allFiles.pop(index)
    f=open(aux,'w')
    f.writelines(allFiles)
    f.close()
