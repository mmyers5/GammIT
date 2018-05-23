import numpy as np
import matplotlib.pyplot as plt
import model_utils
from astropy.cosmology import WMAP9 as cosmo
def lum_dist(infile):
  hosts = np.genfromtxt(infile, names=True, dtype=None)
  plotted = []
  for line in hosts:
    if line['filter'] != 3.6: continue
    GRB = line['ID']
    if GRB in plotted:
      continue
    else:
      plotted.append(GRB)

    z = line['z']

    if line['source_type'] == 'this':
      arrcol='k'
      color = model_utils.colorize('C3')
      zorder=2
      markersize=10
      alpha=1
      label = 'This Study'
      markeredgewidth=0.5
    elif line['source_type'] == 'grb':
      arrcol='k'
      color = model_utils.colorize('C9')
      zorder=1
      markersize=9
      alpha=1
      label = 'Other Studies'
      markeredgewidth=0.5
    else:
      #color = model_utils.colorize('C7')
      arrcol='lightgray'
      color = 'lightgray'
      zorder=0
      markersize=10
      alpha=1.0
      label = 'Field Galaxies'
      markeredgewidth=0
    
    if line['mab'] > 0:
      d_pc  = cosmo.luminosity_distance(z).value*1e6
      M_AB = line['mab'] - 5*( np.log10(d_pc) - 1)
    else:
      M_AB = line['mab']
    
    if line['mab_unc'] == 0:
      marker = 'o'
      markersize=9
    else:
      marker='o'
    M_AB_unc = line['mab_unc']
    
    plt.errorbar(z, M_AB, yerr=M_AB_unc, color=color, marker=marker,
                 zorder=zorder, markersize=markersize, alpha=alpha,
                 markeredgewidth=markeredgewidth, ecolor=arrcol)
    if line['mab_unc'] == 0:
      plt.arrow(z, M_AB, 0, 0.5, head_width=0.05, color=arrcol, zorder=zorder-1,
                facecolor=None, linewidth=0.1)
  plt.gca().set_xlabel('Redshift', fontsize=12)
  plt.gca().set_ylabel('Absolute Magnitude (AB)', fontsize=12)
  plt.gca().invert_yaxis()
  plt.minorticks_on()

  plt.plot(z, M_AB, 'o', mew=0, ms=9, zorder=-2, label='This Study', color=model_utils.colorize('C3'))
  plt.plot(z, M_AB, 'o', mew=0, ms=9, zorder=-2, label='Other Studies', color=model_utils.colorize('C9'))
  #plt.plot(z, M_AB, 'o', mew=0, ms=9, zorder=-2, label='Field Galaxies', color=model_utils.colorize('C7'))
  plt.plot(z, M_AB, 'o', mew=0, ms=9, zorder=-2, label='Field Galaxies', color='lightgray')

  plt.legend(loc='lower right', numpoints=1)
