import sys
import os
import numpy as np

myPath  = '/scratch2/GRBH-MZSFR/Analysis/'
grbfile = myPath+'galbuilds/grbs.txt'
myGRBs  = np.genfromtxt(grbfile,dtype=None)
notes   = []
for GRB in myGRBs:
  modelGRB = '{folder}/{grb}/{grb}_ch1_model.inp'.format(folder=myPath,grb=GRB)
  constGRB = '{folder}/{grb}/{grb}_ch1.constraints'.format(folder=myPath,grb=GRB)
  if os.path.isfile(modelGRB):
    sysString = 'less {modelfile}'.format(modelfile=modelGRB)
    os.system(sysString)
  else: continue
  if os.path.isfile(constGRB):
    sysString = 'less {constfile}'.format(constfile=constGRB)
    os.system(sysString)

  noteMsg = 'Make a note of GRB {grb}? '.format(grb=GRB)
  note = raw_input(noteMsg)
  if (note=='no') or (note=='') or (note=='n'):
    continue
  elif (note=='quit') or (note=='q'):
    break
  else:
    noteString = 'GRB {grb}: {note}'.format(grb=GRB,note=note)
    notes.append(noteString)
f=open('model_notes.txt','w')
f.write('\n'.join(notes))
f.close()
