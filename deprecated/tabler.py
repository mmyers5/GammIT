import numpy as np
import sys, os
from collections import OrderedDict

if len(sys.argv) == 4:
  inFile  = sys.argv[1]
  outFile = sys.argv[2]
  targetAp= sys.argv[3]
  if not os.path.isfile(inFile):
    sys.exit('{inFile} does not exist'.format(inFile=inFile))
else:
  sys.exit('Usage: tabler.py [input] [output] [aperture]')

inDat    = np.genfromtxt(inFile,names=True,dtype=None)
outWrite = open(outFile,'w')

GRBList = np.unique(inDat['GRB'])
myTable = OrderedDict.fromkeys(GRBList)
current   = inDat['GRB'][0]
tempTable = ['','','','']
aps = ['2-2-6','3-3-7']
for idx,grb in enumerate(inDat['GRB']):
  if grb!=current:
    myTable[current]=' & '.join(tempTable)+r' \\'
    tempTable = ['','','','']
    current   = grb
  if inDat[idx]['usage']=='*':
    ap  = inDat[idx]['ap']
    otherAp = [i for i in aps if i!=ap][0]
    fType   = inDat[idx]['type']
    ch = inDat[idx]['ch']
    if ap == targetAp:
      mab = inDat[idx]['flux']
      mab_unc = inDat[idx]['flux_unc']
    elif otherAp == targetAp:
      mab = [i['flux'] for i in inDat if i['GRB']==grb and i['ch']==ch\
                              and i['type']==fType and i['ap']==otherAp][0]
      mab_unc = [i['flux_unc'] for i in inDat if i['GRB']==grb and i['ch']==ch\
                              and i['type']==fType and i['ap']==otherAp][0]
    mab = round(mab,2)
    mab_unc= round(mab_unc,2)
    if fType == 'unc':
      mab_write = str('<{upper_limit}'.format(upper_limit=mab_unc))
    elif mab < 0:
      mab_write = str('<{upper_limit}'.format(upper_limit=3*mab_unc))
    else:
      mab_write = str('{mab}\pm{mab_unc}'.format(mab=mab,mab_unc=mab_unc))
    if inDat[idx]['usage']=='*' and inDat[idx]['ap']==targetAp:
      mab_write = r'$\mathbf{'+mab_write+'}$'
    else:
      mab_write = '$'+mab_write+'$'
    tempTable[ch-1] = mab_write 

myTable[current]=' & '.join(tempTable)+r' \\ \hline'
for key,val in myTable.items():
  outWrite.write('{grb} & {mab}\n'.format(grb=key,mab=val))
outWrite.close()
