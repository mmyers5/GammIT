import numpy as np
from sedpy.observate import load_filters,Filter


def do_it_all(filename,outfile):
  catalog   = np.genfromtxt(filename, comments='#', dtype=None, names=True)
  starsInd  = np.where(catalog['col1']=='*')[0]
  stars     = catalog[starsInd]['col2']
  bigstr  = r'{GRB} & {z:.4f} & {Av:.4f} & {filt} & ${mag:2.2f}\pm{unc:1.2f}$ & \cite{huh}\\'
  bigstr2  = r'{GRB} & {z} & {Av} & {filt} & ${mag:2.2f}\pm{unc:1.2f}$ & \cite{huh}\\'
  bigstr3  = r'{GRB} & {z} & {Av} & {filt} & $>{mag:2.2f}$ & \cite{huh}\\'
  bigstr4 = r'{GRB} & {z:.4f} & {Av:.4f} & {filt} & $>{mag:2.2f}$ & \cite{huh}\\'

  with open(outfile,'w') as outf:
    for objid in stars:
      liststr = []
      outArr, zred, Av  = read_text(filename, objid)
      i = 0
      for line in outArr:
        filt, mag, unc = line[0], float(line[1]), float(line[2])
        filt  = filt.split('_')[-1]
        if i == 0:
          if unc==0:
            liststr.append(bigstr4.format(GRB=objid, z=zred, Av=Av, filt=filt,
                                          mag=mag, unc=unc, huh='{}'))
          else:
            liststr.append(bigstr.format(GRB=objid, z=zred, Av=Av, filt=filt,
                                         mag=mag, unc=unc, huh='{}'))
        else:
          if unc==0:
            liststr.append(bigstr3.format(GRB='', z='', Av='', filt=filt,
                                          mag=mag, unc=unc, huh='{}'))
          else:
            liststr.append(bigstr2.format(GRB='', z='', Av='', filt=filt,
                                          mag=mag, unc=unc, huh='{}'))
        i+=1
      biggeststr = '\n'.join(liststr) + '\n'
      outf.write(biggeststr)
  

def read_text(filename, objid):  
  filt_dir  = '/scratch2/GRBH-MZSFR/prospector/sedpy/sedpy/data/filters/'
  catalog = np.genfromtxt(filename, comments='#', dtype=None, names=True)
  ind     = np.where(catalog['col2'] == objid)[0][0]
  stars   = np.where(catalog['col1'] == '*')[0]
  
  if ind == stars[-1]:
    end = len(catalog)
  else:
    end  = stars[ np.where(ind == stars)[0][0] +1 ]
   
  filternames = catalog['col1'][ind+1: end]
  mags  = catalog['col2'][ind+1: end]
  uncs  = catalog['col3'][ind+1: end]
  zred  = float(catalog['col3'][ind])
  Av    = float(catalog['col4'][ind])
  Av_corrmask = catalog['col4'][ind+1:end]

  Av_corrmask = Av_corrmask=='True'
  det_mask    = uncs != '-'
  uncs[~det_mask] = 0
  mags  = mags.astype(np.float)
  uncs  = uncs.astype(np.float)
  corrs = correct_ext(filternames, zred) * (Av/0.73)
  mags[~Av_corrmask]  = mags[~Av_corrmask] + corrs[~Av_corrmask]
  if np.ndim(filternames)>0:
    filt_list = load_filters(filternames, directory=filt_dir)
    wav_eff   = np.array([filt.wave_effective for filt in filt_list])
  elif np.ndim(filters)==0:
    filt_list = Filter(filternames,directory=filt_dir)
    wav_eff   = filt_list.wave_effective
  sort_indices = np.argsort(wav_eff, axis=0)
  mags  = mags[sort_indices]
  uncs  = uncs[sort_indices]
  filternames = filternames[sort_indices]
  return zip(*(filternames, mags, uncs)), zred, Av

def correct_ext(filters, zred, exttype='MW'):
  """Code from Tanmoy Laskar to calculate galactic foreground reddening
    corrections.

  :param filters:
    Filter names. See filt_dir for available filters.

  :param zred:
    Object redshift.

  :param exttype: (optional, default: 'MW')
    Extinction type, could be 'MW', 'LMC', or 'SMC'. Milky way for lyfe!

  :returns:
    The correction to the filter in AB magnitudes.
  """
  filt_dir  = '/scratch2/GRBH-MZSFR/prospector/sedpy/sedpy/data/filters/'

  # Calculate "effective" wavelengths in angstroms
  if np.ndim(filters)>0:
    filt_list = load_filters(filters, directory=filt_dir)
    wav_eff   = np.array([filt.wave_effective for filt in filt_list])
  elif np.ndim(filters)==0:
    filt_list = Filter(filters,directory=filt_dir)
    wav_eff   = filt_list.wave_effective
  # Go from observed wavelength to rest wavelength in micron
  wav = wav_eff*1e-4/(zred+1)
  
  # Take it away, code!
  if (exttype == ''):
    raise ValueError("Extinction curve type must be set first!")
  else:
    if (exttype == 'MW'):
      a = np.array([165., 14., 0.045, 0.002, 0.002, 0.012])
      l = np.array([0.047, 0.08, 0.22, 9.7, 18., 25])
      b = np.array([90., 4., -1.95, -1.95, -1.8, 0.0])
      n = np.array([2., 6.5, 2., 2., 2., 2.])
      K = np.array([2.89, 0.31, 0.16, 0.31, 0.28, 0.76])
    elif (exttype == 'LMC'):                
      a = np.array([175., 19., 0.023, 0.005, 0.006, 0.020])
      l = np.array([0.046, 0.08, 0.22, 9.7, 18., 25.])
      b = np.array([90., 5.5, -1.95, -1.95, -1.8, 0.0])
      n = np.array([2., 4.5, 2., 2., 2., 2.])
      K = np.array([3.00, 0.56, 0.08, 0.77, 0.86, 1.26])
    elif (exttype == 'SMC'):
      a = np.array([185., 27., 0.005, 0.010, 0.012, 0.030])
      l = np.array([0.047, 0.08, 0.22, 9.7, 18., 25])
      b = np.array([90., 5.5, -1.95, -1.95, -1.8, 0.0])
      n = np.array([2., 4.0, 2., 2., 2., 2.])
      K = np.array([2.89, 0.91, 0.02, 1.55, 1.72, 1.89])

    dims = np.ndim(wav)
    if (dims == 0):
      return sum(a/((wav/l)**n + (l/wav)**n + b))
    elif (dims == 1):
      return np.array([sum(a/((w/l)**n + (l/w)**n + b)) for w in wav])
    elif (dims == 2):
      return np.array(wav.shape[0]*[[sum(a/((w/l)**n + (l/w)**n + b)) for w in wav[0]]])
    else:
      raise ValueError('Hmm... Not quite sure what to do with a frequency array that has %i dimensions in extinction.extcurve' %(dims))
