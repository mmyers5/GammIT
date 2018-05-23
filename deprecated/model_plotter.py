import fits_plot as fip
import numpy as np
import matplotlib.pyplot as plt

def model_plotter(inFile):
  apDict = {'2': (2.4, 7.2), '3': (3.6, 8.4)}
  master = np.genfromtxt(inFile, names=True, dtype=None)
  for i in range(master.shape[0]):
    GRB = master[i]['GRB']
    xy  = (master[i]['ra'], master[i]['dec'])
    ch  = master[i]['ch']
    ap  = apDict[str(master[i]['ap'])]
    z   = master[i]['z']
    incirc = master[i]['circol']
    figDict = fip.make_model_stamp(GRB, ch, xy, incirc, z, ap=ap)
    figFile  = 'plots/stamps/{GRB}_ch{ch}_modelstamp.eps'.\
               format(GRB=GRB, ch=ch) 
    finished = 'n'

    while finished!='y':
      if finished =='quit' or finished =='next': break
      else: finished,figDict = query_func(figFile,figDict,GRB,ch,xy,incirc,z,ap)
    plt.close('all')
    if finished == 'quit': return
    elif finished == 'nextone': continue

def query_func(figFile,figDict,GRB,ch,xy,incirc,z, ap):
  n = 'n'
  y = 'y'
  quit='quit'
  nextone='nextone'

  query = input("Savefile? (y/n/quit/nextone):\n{} ".format(figFile))
  if query=='y':
    "Plot saved."
    plt.savefig(figFile, format='eps', dpi=plt.gcf().get_dpi(),
                transparent=True, pad_inches=0)
    return query,figDict
  
  elif query=='n':
    print "Suit yourself."
    return query,figDict

  elif query=='quit':
    print "Quitting."
    return query, figDict

  elif query=='next':
    print "Going to next."
    return query, figDict

  elif type(query)==str:
    try: plt.savefig(query, format='eps', dpi=plt.gcf().get_dpi(),
                     transparent=True, pad_inches=0)
    except IOError: 
      print '{} does not exist.'.format(query)
      skipped.append(figFile)
      return 'n',figDict
    else: return 'y',figDict

  elif type(query)==tuple:
    try: c1,c2 = query
    except ValueError:
      print 'Invalid inputs!'
      return 'n',figDict
    else: 
      figDict = fip.make_model_stamp(GRB, ch, xy, incirc, z, ap=ap, c1=c1, c2=c2)
      return 'n',figDict
