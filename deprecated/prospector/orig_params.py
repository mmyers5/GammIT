import numpy as np
from prospect.models import priors, sedmodel
from prospect.sources import CSPSpecBasis
from sedpy.observate import load_filters,Filter
from astropy.cosmology import WMAP9 as cosmo

# --------------
# RUN_PARAMS
# --------------

run_params = {'verbose':True,
              'debug':False,
              'outfile':'demo_galphot',
              # Fitter parameters
              'nwalkers':128,
              'nburn':[10, 10, 10], 'niter':1280,
              'do_powell': False,
              'do_levenburg': True, 'nmin':20,
              'ftol':0.5e-5, 'maxfev':5000,
              'initial_disp':0.1,
              # Obs data parameters
              'objid':'0',
              'phottable': 'host_photometry.dat',
              'logify_spectrum':False,
              'normalize_spectrum':False,
              # SPS parameters
              'zcontinuous': 1,
              }

# --------------
# OBS
# --------------
def get_z(phottable='demo_photometry.dat',**kwargs):
  """Get redshift from the photometry table
  
  :param phottable: (optional, default: 'demo_photometry.dat')
    File with host photometry. Columns are filters, host/ABmag(or the upper
    limit), redshift/ABmag error, and Av/TrueOrFalse if corrected for
    extinction. Different objects are delimited by rows that stars with an *.

  :param kwargs: (optional)
    Dictionary used by prospector.

  :returns zred:
    The redshift of the object.
  """
  objid = kwargs['objid']
  catalog = np.genfromtxt(phottable, comments='#', dtype=None, names=True)
  ind   = np.where(catalog['col2'] == objid)[0][0]
  zred  = float( catalog['col3'][ind] )
  return zred

def load_obs(phottable='demo_photometry.dat', **kwargs):
  """Primary function in reading the photometry table.

  :param phottable: (optional, default: 'demo_photometry.dat')
    File with host photometry. Columns are filters, host/ABmag(or the upper
    limit), redshift/ABmag error, and Av/TrueOrFalse if corrected for
    extinction. Different objects are delimited by rows that stars with an *.

  :param kwargs: (optional)
    Dictionary used by prospector.

  :returns obs:
    Dictionary with the necessary observational data needed for prospector fits.
  """
  objid = kwargs['objid']
  catalog = np.genfromtxt(phottable, comments='#', dtype=None, names=True)
  ind   = np.where(catalog['col2'] == objid)[0][0]
  stars = np.where(catalog['col1'] == '*')[0]
  if ind==stars[-1]:
    end = len(catalog)
  else:
    end = stars[ np.where(ind == stars)[0][0] + 1 ]

  # Lots of unpacking
  filternames = catalog['col1'][ind+1:end]
  mags  = catalog['col2'][ind+1:end]
  uncs  = catalog['col3'][ind+1:end]
  zred  = float(catalog['col3'][ind])
  Av    = float(catalog['col4'][ind])
  Av_corrmask = catalog['col4'][ind+1:end]
  # Play with some masks and ensure proper data types
  Av_corrmask = Av_corrmask=='True'
  det_mask    = uncs != '-'   # True if detected
  uncs[~det_mask]  = 0        # uncertainty of 0 if not detected
  mags  = mags.astype(np.float)
  uncs  = uncs.astype(np.float)

  # Perform extinction correction on uncorrected galaxies
  corrs = correct_ext(filternames, zred) * (Av/0.73)
  mags[~Av_corrmask] = mags[~Av_corrmask] + corrs[~Av_corrmask]
  
  # Put it all together!
  obs = {}
  obs['filters'] = load_filters(filternames)
  # Convert from AB mags to maggies, use higher flux limit for uncertainty
  obs['maggies']      = np.squeeze( 10**(-mags/2.5) )
  obs['maggies_unc']  = np.abs(obs['maggies'] - np.squeeze( 10**(-(mags-uncs) / 2.5) ))
  obs['phot_mask']    = np.isfinite( np.squeeze(mags) )     # Not really used
  obs['wavelength']   = None
  obs['objid']        = objid
  obs['det_mask'] = det_mask
  obs['zred']     = zred
  obs['Av']       = Av
  obs['Av_corrmask']  = Av_corrmask
  obs['Av_corrs']     =  corrs
  
  return obs

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
            return np.array(wav.shape[0]*[[sum(a/((w/l)**n + (l/w)**n + b)) for
w in wav[0]]])
        else:
            raise ValueError('Hmm... Not quite sure what to do with a frequency\
array that has %i dimensions in extinction.extcurve' %(dims))

# --------------
# SPS Object
# --------------

def load_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    sps = CSPSpecBasis(zcontinuous=zcontinuous,
                       compute_vega_mags=compute_vega_mags)
    return sps

# -----------------
# Gaussian Process
# ------------------

def load_gp(**extras):
    return None, None

# --------------
# MODEL_PARAMS
# --------------

# You'll note below that we have 5 free parameters:
# mass, logzsol, tage, tau, dust2
# They are all scalars.
# Just kidding, the free parameters are mass, tage, and dust2.
#
# The other parameters are all fixed, but we want to explicitly set their
# values, possibly from something differnt than the FSPS defaults

model_params = []

# --- SFH --------
# FSPS parameter.  sfh=4 is a delayed-tau SFH
model_params.append({'name': 'sfh', 'N': 1,
                        'isfree': False,
                        'init': 0,
                        'units': 'type'
                    })

# Normalization of the SFH.  If the ``mass_units`` parameter is not supplied,
# this will be in surviving stellar mass.  Otherwise it is in the total stellar
# mass formed.
model_params.append({'name': 'mass', 'N': 1,
                        'isfree': True,
                        'init': 1e7,
                        'init_disp': 1e6,
                        'units': r'M_\odot',
                        'prior':priors.TopHat(mini=1e5, maxi=1e11)})

# Since we have zcontinuous=1 above, the metallicity is controlled by the
# ``logzsol`` parameter.
model_params.append({'name': 'logzsol', 'N': 1,
                        'isfree': False,
                        'init': -0.3,
                        'init_disp': 0.01,
                        'units': r'$\log (Z/Z_\odot)$',
                        'prior': priors.TopHat(mini=-0.6, maxi=0.0)})

# FSPS parameter
model_params.append({'name': 'tau', 'N': 1,
                        'isfree': False,
                        'init': 1.0,
                        'units': 'Gyr',
                        'prior':priors.LogUniform(mini=0.1, maxi=100)})


# FSPS parameter
model_params.append({'name': 'fburst', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': '',
                        'prior':priors.TopHat(mini=0.0, maxi=1.0)})

# --- Dust ---------
# FSPS parameter
model_params.append({'name': 'dust1', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': '',
                        'prior':priors.TopHat(mini=0.1, maxi=2.0)})

# FSPS parameter
model_params.append({'name': 'dust2', 'N': 1,
                        'isfree': True,
                        'init': 0.35,
                        'reinit': True,
                        'init_disp': 0.3,
                        'units': '',
                        'prior':priors.TopHat(mini=0.0, maxi=2.0)})

# FSPS parameter
model_params.append({'name': 'dust_index', 'N': 1,
                        'isfree': False,
                        'init': -0.7,
                        'units': '',
                        'prior':priors.TopHat(mini=-1.5, maxi=-0.5)})

# FSPS parameter
model_params.append({'name': 'dust1_index', 'N': 1,
                        'isfree': False,
                        'init': -1.0,
                        'units': '',
                        'prior':priors.TopHat(mini=-1.5, maxi=-0.5)})

# FSPS parameter
model_params.append({'name': 'dust_tesc', 'N': 1,
                        'isfree': False,
                        'init': 7.0,
                        'units': 'log(Gyr)',
                        'prior': None})

# FSPS parameter
model_params.append({'name': 'dust_type', 'N': 1,
                        'isfree': False,
                        'init': 2,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'add_dust_emission', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'duste_umin', 'N': 1,
                        'isfree': False,
                        'init': 1.0,
                        'units': 'MMP83 local MW intensity'})

# --- Stellar Pops ------------
# FSPS parameter
model_params.append({'name': 'tpagb_norm_type', 'N': 1,
                        'isfree': False,
                        'init': 2,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'add_agb_dust_model', 'N': 1,
                        'isfree': False,
                        'init': True,
                        'units': 'index'})

# FSPS parameter
model_params.append({'name': 'agb_dust', 'N': 1,
                        'isfree': False,
                        'init': 1,
                        'units': 'index'})

# --- Nebular Emission ------

# FSPS parameter
model_params.append({'name': 'add_neb_emission', 'N': 1,
                     'isfree': False,
                     'init': True})

# Here is a really simple function that takes a **dict argument, picks out the
# `logzsol` key, and returns the value.  This way, we can have gas_logz find
# the value of logzsol and use it, if we uncomment the 'depends_on' line in the
# `gas_logz` parameter definition.
#
# One can use this kind of thing to transform parameters as well (like making
# them linear instead of log, or divide everything by 10, or whatever.) You can
# have one parameter depend on several others (or vice versa).  Just remember
# that a parameter with `depends_on` must always be fixed.  It's also not a
# good idea to have one parameter depend on another parameter that *also*
# depends on something, since dependency resolution order is arbitrary.

def stellar_logzsol(logzsol=0.0, **extras):
    return logzsol


# FSPS parameter
model_params.append({'name': 'gas_logz', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': r'log Z/Z_\odot',
#                        'depends_on': stellar_logzsol,
                        'prior':priors.TopHat(mini=-2.0, maxi=0.5)})

# FSPS parameter
model_params.append({'name': 'gas_logu', 'N': 1,
                        'isfree': False,
                        'init': -2.0,
                        'units': '',
                        'prior':priors.TopHat(mini=-4, maxi=-1)})

# --- Calibration ---------
model_params.append({'name': 'phot_jitter', 'N': 1,
                        'isfree': False,
                        'init': 0.0,
                        'units': 'mags',
                        'prior':priors.TopHat(mini=0.0, maxi=0.2)})

def load_model(**extras):
    # In principle (and we've done it) you could have the model depend on
    # command line arguments (or anything in run_params) by making changes to
    # `model_params` here before instantiation the SedModel object.  Up to you.
    zred  = get_z(phottable=extras['phottable'],objid=extras['objid'])
    # --- Distance ---
    # This is the redshift.  Because we are not separately supplying a ``lumdist``
    # parameter, the distance will be determined from the redshift using a WMAP9
    # cosmology, unless the redshift is 0, in which case the distance is assumed to
    # be 10pc (i.e. for absolute magnitudes)
    model_params.append({'name': 'zred', 'N': 1,
                            'isfree': False,
                            'init': zred,
                            'units': '',
                            'prior':priors.TopHat(mini=0.0, maxi=4.0)})

    # FSPS parameter
    model_params.append({'name': 'tage', 'N': 1,
                            'isfree': True,
                            'init': 0.0,
                            'init_disp': 1.0,
                            'units': 'Gyr',
                            'prior':priors.TopHat(mini=0.005,maxi=cosmo.age(zred).value)})
    # FSPS parameter
    model_params.append({'name': 'tburst', 'N': 1,
                            'isfree': False,
                            'init': 0.0,
                            'units': '',
                            'prior':priors.TopHat(mini=0.0,maxi=cosmo.age(zred).value)})
    return sedmodel.SedModel(model_params)
