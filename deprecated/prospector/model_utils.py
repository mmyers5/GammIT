import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import corner
from scipy.stats import mode

"""Utility functions for use with model plotting."""

__all__ = ["scale_units"]

def scale_units(units, *values):
  """Will convert input values to some other unit, given the inputs are in
  either maggies or angstroms.

  :param units:
    Valid units to convert to are 'uJy', 'cgs', 'maggies', 'um', or 'Hz'.
    Using 'Hz' or 'um' implies that values are in units of angstroms, while
    all else should be in units of maggies.
  
  :param values:
    Some sort of sequence of values that need conversion.
  """
  values = np.array(values)
  unitsDict = {'ujy'      : 3631*1e6,
               'cgs'      : 3631*1e-23,
               'maggies'  : 1,
               'um'       : 1e-4,
               'hz'       : 3e18/values}

  try:
    scale = unitsDict[units.lower()]
  except KeyError:
    raise KeyError("Please use one of these units:\n{}".\
              format( unitsDict.keys() ))

  valueScaled = values*scale
  if valueScaled.shape[0] == 1: return valueScaled[0]
  else: return valueScaled

def colorize(CN):
  """For coloring plots since I can't get pyplot to work with these color
  formats.
  
  :param CN:
    String of format C#, where # is from 0-9.

  :returns color:
    A color code.
  """
  color_dict = {'C0'  : '#1f77b4',  # blue
                'C1'  : '#ff7f0e',  # orange
                'C2'  : '#2ca02c',  # green
                'C3'  : '#d62728',  # red
                'C4'  : '#9467bd',  # purple
                'C5'  : '#8c564b',  # brown
                'C6'  : '#e377c2',  # pink
                'C7'  : '#7f7f7f',  # gray
                'C8'  : '#bcbd22',  # yellow
                'C9'  : '#17becf'}  # cyan

  try:
    color = color_dict[CN]
  except KeyError:
    raise KeyError("Please use one of these keys:\n{}".\
            format( color_dict.keys() ))

  return color

def sed_figure_init(overfig):
  if overfig is None:
    sedfig = plt.figure(num=None, figsize=(12,9), dpi=80)
    sedfig.subplots_adjust(hspace=0)
    gs  = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    ax0 = sedfig.add_subplot(gs[0])
    ax1 = sedfig.add_subplot(gs[1], sharex=ax0)
  else:
    sedfig = overfig[0]
    ax0 = overfig[1]
    ax1 = overfig[2]

  return sedfig, ax0, ax1

def y_getter(x, imin, imax, z, zx):
  valid = z[(zx<imax) & (zx>imin)]
  ythresh = mode(valid)
  if ythresh[1][0] > 1:
    yval  = ythresh[0][0]
    ymin, ymax  = 0.90*yval, 2*yval
  else:
    ymin = valid.min()
    ymax = valid.max()
    if ymin > x.min(): ymin = x.min()
    if ymax < x.max(): ymax = x.max()
    ymin, ymax = 0.90*ymin, 2*ymax
  #ymin, ymax = 1e-15, .5e-12
  return ymin, ymax

def x_getter(x,y):
  xmin = 0.90*(x - y).min()
  xmax = 1.10*(x + y).max()
  #xmax = 8
  return xmin, xmax

def axes_formatter(ax0, ax1, z):

  # scale axes
  ax0.set_xscale('log')
  ax1.set_xscale('log')
  ax0.set_yscale('log')

  # plot more readable x axes
  xlist = ax0.get_xticks(minor=True)
  newxlist = []
  for num in xlist:
    if num<1: newxlist.append('{:.1f}'.format(num))
    else: newxlist.append('{:.0f}'.format(num))
  ax0.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
  ax0.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
  #ax0.set_xticklabels(newxlist)

  # add seconds y axis
  ax2 = ax0.twiny()
  ax2.set_xscale('log')
  ax2.set_xticks(xlist)
  newxlist = []
  for num in xlist/(1+z):
    if num<1: newxlist.append('{:.1f}'.format(num))
    else: newxlist.append('{:.0f}'.format(num))
  
  ax2.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
  ax2.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
  #ax2.set_xticklabels(newxlist)

  xbnd = ax0.get_xbound()
  ax2.set_xbound(xbnd[0]/(1+z), xbnd[1]/(1+z))
  ax2.yaxis.set_visible(False)
  return ax0, ax1, ax2

def axes_labeler(units):
  if units.lower() == 'cgs':
    label = r'Specific Flux $\nu f_{\nu}$ $\mathrm{(erg/s/cm^{2})}$'
  elif units.lower() == 'umobs':
    label = r'Wavelength $(\mu m)$'
  elif units.lower() == 'umrest':
    label = r'Rest Wavelength $(\mu m)$'
  elif units.lower() == 'chi':
    label = r'$\chi$'
  return label

def thetas_labeler(label, chain, masslogify=True, dustify=True,
                   agelogify=True):
  for idx, param in enumerate(label):
    if param == 'mass':
      if masslogify == True:
        label[idx]   = r'$\log_{10}\left(M_*/M_{\odot}\right)$'
        chain[:, idx] = np.log10(chain[:, idx])
      else:
        label[idx] = r'Mass $(M_{\odot})$'
    elif param == 'dust2':
      if dustify == True:
        label[idx]  = r'$A_\mathrm{V}$'
        chain[:, idx] /= 0.921
      else:
        label = 'Optical Depth in V' 
    elif param == 'tage':
      if agelogify == False:
        if np.median(chain[:, idx]) < 0.1:
          label[idx]  = r'Age $\mathrm{(Myr)}$'
          chain[:, idx] *= 1e3
        else:
          label[idx]  = r'Age $\mathrm{(Gyr)}$'
      elif agelogify == True:
        label[idx]   = r'$\log_{10}\left(t_\mathrm{age}/\mathrm{yr}\right)$'
        chain[:, idx] = np.log10(chain[:, idx]*1e9)

  return label, chain

def plot_corner_extra(flatChain, extra, cornerBins, cornerRange):
  idx = extra['index']
  val = extra['value']

  chainMask = flatChain[:, idx] > val
  hist_kwargs = {'alpha': 0.4, 'histtype': 'stepfilled'}
  data_kwargs = {'alpha': 0.2}
  cornerfig = corner.corner(flatChain[chainMask,:-1],
                            color=colorize('C2'), bins=cornerBins,
                            range=cornerRange[:-1], show_titles=False,
                            plot_contours=False, no_fill_contours=True,
                            plot_density=False, hist_kwargs=hist_kwargs, 
                            data_kwargs=data_kwargs)
  hist_kwargs = {'alpha': 0.4, 'histtype': 'stepfilled'}
  data_kwargs = {'alpha': 0.2}
  cornerfig = corner.corner(flatChain[~chainMask,:-1],
                            color=colorize('C3'), bins=cornerBins,
                            range=cornerRange[:-1], show_titles=False,
                            plot_contours=False, no_fill_contours=True,
                            plot_density=False, hist_kwargs=hist_kwargs, 
                            data_kwargs=data_kwargs, fig=cornerfig)
  return cornerfig

def theta_cut(flatChain, extra, typicalCut=50):
  idx = extra['index']
  val = extra['value']

  chainMask = flatChain[:, idx] > val

  upper = np.percentile(flatChain[:, -1][chainMask], typicalCut)
  upperIndex  = np.abs(flatChain[:, -1][chainMask] - upper).argmin()
  lower = np.percentile(flatChain[:, -1][~chainMask], typicalCut)
  lowerIndex  = np.abs(flatChain[:, -1][~chainMask] - lower).argmin()

  thetaPair = np.array((flatChain[chainMask][upperIndex, :-1],
                        flatChain[~chainMask][lowerIndex, :-1]))
  return thetaPair
