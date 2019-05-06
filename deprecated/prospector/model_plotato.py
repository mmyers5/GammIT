import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl
from copy import deepcopy as dc
import prospect.io.read_results as rr
import corner
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import random

"""Functions for plotting meaningful results from prospector output files"""

__all__ = ["flatten_results", "burn_chain", "get_theta_max", "model_spec",
           "obs_phot", "angstrom_micron", "wav_nu", "plot_max_sed",
           "spec_range", "proper_label", "plot_corner", "plot_corner_extra",
           "theta_cut", "quantilize", "plot_burn_chain", "colorize"]

def flatten_results(filename, thresh=1, order='C'):
  """Flatten chain from results file
  :param filename:
    Name and path to output prospector file ending in 'mcmc'.

  :param thresh: (optional, default: 1)
    Number of sigma away from which to burn chain (see burn function).
    If burn unwanted, set to numpy.inf

  :param order: (optional, default: 'C')
    Passed to numpy.reshape to determine order of rows. 'C' will order by
    walkers and their iterations, 'F' will order by iterations and their
    walkers.

  :returns flatchain:
    An array containing all of the walker and iteration model values
    where each column is a theta, with the last column as the corresponding
    lnprobability. By default, the rows are ordered by walkers and their
    respective iterations.
  """

  results, obs, model = rr.results_from(filename)
  iter_burn, walker_noburn  = burn_chain(results,thresh=thresh)   # index below which to burn

  #rows  = (results['run_params']['nwalkers']) * \
  #        (results['run_params']['niter'] - burn)
  rows  = len(walker_noburn) * \
          (results['run_params']['niter'] - iter_burn)
  
  cols  = len(results['theta_labels'])

  chain   = results['chain']
  lnprob  = results['lnprobability']

  flatchain = np.zeros((rows, cols+1))
  flatchain[:,:cols]  = chain[walker_noburn,iter_burn:,:].reshape(rows, cols, order=order)
  flatchain[:,-1]     = lnprob[walker_noburn,iter_burn:].reshape(rows, order=order)

  return flatchain
  
def burn_chain(results, thresh=1, order='C'):
  """Read the lnprobability in a results file, get the mean of last 
  10% of each walker's iterations, and find the index at which the mean
  of the walkers' iterations first falls above the mean minus some delta,
  where delta is the standard deviation of the mean times thresh.
 
  :param results:
    The results dictionary from an mcmc file.
  
  :param thresh: (optional, default: 1)
    The threshold sigma distance from the mean that determines when to burn.
  
  :param order: (optional, default: 'C')
    Passed to numpy.reshape to determine order of rows. 'C' will order by
    walkers and their iterations, 'F' will order by iterations and their
    walkers.

  :returns burn:
    The index below which to burn iterators.
  """
  lnprob  = results['lnprobability']
  iters   = results['run_params']['niter']
  # use the last 10% of iterations from each walker
  iter_few     = lnprob[:, -int(0.10*iters)+1:]

  walker_mean  = np.median(iter_few, axis=1)    # mean of walkers over iters
  walker_sigma = np.std(walker_mean, axis=0)
  mean  = np.median(walker_mean, axis=0)    # mean of all walkers

  # compare to the mean of all of the walkers per iteration 
  iter_mean = np.median(lnprob, axis=0)     # mean of iters over walkers

  # retrieve index where mean over iterations first exceeds (mean-delta)
  iter_burn  = np.where( iter_mean > mean - (thresh*walker_sigma) )[0][0]
  leftovers  = np.median(lnprob[:,iter_burn:], axis=1)
  walker_noburn = np.where(leftovers > mean-thresh*walker_sigma)[0]

  return iter_burn, walker_noburn
    
def get_theta_max(results):
  """Read the results file from prospector and get the thetas that correspond to
  the maximum likelihood.

  :param results:
    The results dictionary from an mcmc file.

  :returns theta_max:
    An array with the maximum theta results.
  """
  ind_max = results['lnprobability'].argmax()
  walker, iteration  = np.unravel_index(ind_max, results['lnprobability'].shape)
  theta_max = results['chain'][walker,iteration,:]

  return theta_max

def model_spec(theta, model, obs, sps, units='cgs'):
  """Get the spectrum and photometry given the input parameters.

  :param theta:
    The model parameters.
  
  :param model:
    The model information from prospector.

  :param obs:
    The input observation dictionary from prospector.

  :param sps:
    The sps information from prospector.

  :param units: (optional, default: 'uJy'):
    To convert the output spectrum and photometry to uJy instead of maggies.

  :returns spec:
    Spectrum that corresponds to wavelengths in obs['wavelength']

  :returns phot:
    The photometry in each observed filter integrated over the spectrum.
  """
  spec, phot, x = model.mean_model(theta,obs=obs,sps=sps)

  if units.lower()=='ujy':
    spec *= 3631*1e6
    phot *= 3631*1e6

  if units.lower()=='cgs':
    spec *= 3631 * 1e-23
    phot *= 3631 * 1e-23
  
  if units.lower()=='maggies':
    return spec,phot

  return spec, phot

def obs_phot(results, units='cgs'):
  """Given results, get the photometric detections.

  :param results:
    The results dictionary from an mcmc file.

  :param units: (optional, default: 'uJy')
    To convert the output spectrum and photometry to uJy instead of maggies.

  :returns orig:
    Observed photometry array.
  
  :returns err:
    One-sigma photometry limits.

  :returns det_mask:
    The mask that specifies detections in results.
  """
  det_mask = results['obs']['det_mask']       # detections mask
  results['obs']['maggies_unc'][~det_mask] = results['obs']['maggies'][~det_mask]

  orig    = results['obs']['maggies']
  err     = results['obs']['maggies_unc']
  wav_err = np.array([f.rectangular_width for f in results['obs']['filters']])

  if units.lower()=='ujy':
    orig    *= 3631*1e6
    err     *= 3631*1e6
  elif units.lower()=='cgs':
    orig    *= 3631*1e-23
    err     *= 3631*1e-23

  return orig, err, wav_err, det_mask

def angstrom_micron(*wavelengths):
  """Given wavelengths in angstroms, return wavelengths in microns.

  :param wavelengths:
    A list of wavelengths in angstroms.

  :returns wav_micron:
    A tuple of wavelengths in microns.
  """
  wav_angstrom = np.array(wavelengths)
  wav_micron   = wav_angstrom * 1e-4

  return tuple(wav_micron)

def wav_nu(wavelengths, units='um'):
  """Given wavelengths in angstroms or microns, return frequency in Hz.

  :param wavelengths:
    A list of wavelengths in microns or angstroms.

  :param units: (optional, default: 'um')
    Units of wavelengths, either a single string for all wavelengths,
    or a list per wavelength set.

  :returns nu_Hz:
    A tuple of frequencies in Hz.
  """
  wav_angstrom = np.array(wavelengths)

  if units.lower()=='um' or units.lower()=='micron':
    c = 3.0e14
  elif units.lower()=='a' or units.lower()=='angstrom':
    c = 3.0e18
  else:
    return 'Yo I need an actual unit bruv...'

  nu_Hz = c/wav_angstrom

  return tuple(nu_Hz)
  
def plot_max_sed(myfile, overfig=None,thresh=1, save_file=None, extra_thetas=None, **kwargs):
  """Given an mcmc file, plot the spectrum with photometry and residual info. Will
  presumably save figure files sometime. Incorporates burn specified by thresh.
  
  :param myfile:
    The mcmc prospector file.
  
  :param thresh: (optional, default: 1)
    Passed to burn_chain.

  :param save_file: (optional, default: None)
    If not None, will save the figure to the string specified by save_file.

  :param extra_thetas: (optional, default: None)

  :param kwargs: (optional)
    A dictionary that has stuff. Mostly used to skip re-getting sps.
    Will also pass to spec_range.
    -num_spec
    -quantiles
    -units
  """
  results, obs, model = rr.results_from(myfile)
  theta_max = get_theta_max(results)
  
  # check keyword arguments
  try:
    sps = kwargs['sps']
  except KeyError:              # if key does not exist
    sps = rr.get_sps(results)
    kwargs['sps'] = sps
  except NameError:             # if dict does not exist
    kwargs = {}
    sps = rr.get_sps(results)
    kwargs['sps'] = sps
  font_kwargs = {'fontsize': 14}

  # wavelengths in angstroms
  spec_wave = dc(sps.wavelengths)
  phot_wave = np.array([f.wave_effective for f in results['obs']['filters']])

  # get maximum ln spectrum, original photometry, and quantiled spectra
  obs['wavelength'] = spec_wave   # must be set to get proper spec wavelengths
  spec_max, phot_max = model_spec(theta_max, model, obs, sps,units='cgs')
  orig, err, wav_err, det_mask = obs_phot(results, units='cgs')
  quant_specs = spec_range(results, model, obs, **kwargs)

  # get fitted photometry
  #obs['wavelength'] = phot_wave   # must be set to get proper phot wavelengths
  #spec, phot_max = model_spec(theta_max, model, obs, sps,units='cgs')
  # angstroms to microns
  spec_wave, phot_wave, wav_err = angstrom_micron(spec_wave, phot_wave, wav_err)
  spec_nu, phot_nu = wav_nu([spec_wave, phot_wave])
  #err_nu  = np.abs( phot_nu*err / (phot_nu*phot_max*np.log(10)) )   # error prop
  err_nu = phot_nu*err

  # plotting time! Declare some stuff first.
  if overfig==None:
    sedfig = plt.figure(num=None, figsize=(12,9), dpi=80)
    sedfig.subplots_adjust(hspace=0)
    gs  = mpl.gridspec.GridSpec(2, 1, height_ratios=[3,1])
    ax0 = sedfig.add_subplot(gs[0])
    ax1 = sedfig.add_subplot(gs[1], sharex=ax0)
  else:
    sedfig = overfig[0]
    ax0 = overfig[1]
    ax1 = overfig[2]
  

  # plot spectra
  if overfig==None:
    ax0.plot(spec_wave, spec_nu*spec_max, color=colorize('C0'))
  else:
    ax0.plot(spec_wave, spec_nu*spec_max, color=colorize('C3'))
  if quant_specs is not None:
    if overfig is None:
      ax0.fill_between(spec_wave, spec_nu*quant_specs[0],
                      spec_nu*quant_specs[1], interpolate=True,
                       alpha=0.5, color=colorize('C7'))
    if overfig is not None:
      ax0.fill_between(spec_wave, spec_nu*quant_specs[0],
                      spec_nu*[1], interpolate=True,
                       alpha=0.2, color=colorize('C3'))
  # handle extra
  if extra_thetas is not None:
    obs['wavelength'] = dc(sps.wavelengths)
    for i in range(extra_thetas.shape[0]):
      color = 'C{}'.format(i+2)
      extra_spec, extra_phot = model_spec(extra_thetas[i,:], model, obs, sps)
      ax0.plot(spec_wave, spec_nu*extra_spec, color=colorize(color))

  # plot photometry, first by getting data from this study
  study   = np.array([j.startswith('my_spitzer') for j in obs['filternames']])
  ax0.errorbar(phot_wave[det_mask*study], 
               phot_nu[det_mask*study]*orig[det_mask*study],
               fmt='o', yerr=err_nu[det_mask*study],
               xerr=wav_err[det_mask*study], 
               color=colorize('C3'), ms=9, label='phot')
  ax0.errorbar(phot_wave[~det_mask*study], 
               phot_nu[~det_mask*study]*orig[~det_mask*study],
               fmt='v', yerr=0*err_nu[~det_mask*study],
               xerr=wav_err[~det_mask*study], 
               color=colorize('C3'), ms=9, label='phot')
  ax0.errorbar(phot_wave[det_mask*~study], 
               phot_nu[det_mask*~study]*orig[det_mask*~study],
               fmt='o', yerr=err_nu[det_mask*~study],
               xerr=wav_err[det_mask*~study], 
               color=colorize('C1'), ms=9, label='phot')
  ax0.errorbar(phot_wave[~det_mask*~study], 
               phot_nu[~det_mask*~study]*orig[~det_mask*~study],
               fmt='v', yerr=0*err_nu[~det_mask*~study],
               xerr=wav_err[~det_mask*~study], 
               color=colorize('C1'), ms=9, label='phot')
  """
  ax0.errorbar(phot_wave[det_mask*study], 
               phot_nu[det_mask*study]*orig[det_mask*study],
               fmt='o', yerr=err_nu[det_mask*study],
               xerr=wav_err[det_mask*study], 
               color=colorize('C3'), ms=9, label='phot')
  ax0.errorbar(phot_wave[~det_mask*study], 
               phot_nu[~det_mask*study]*orig[~det_mask*study],
               fmt='v', yerr=0*err_nu[~det_mask*study],
               xerr=wav_err[~det_mask*study], 
               color=colorize('C3'), ms=9, label='phot')
  ax0.errorbar(phot_wave[det_mask*~study], 
               phot_nu[det_mask*~study]*orig[det_mask*~study],
               fmt='o', yerr=err_nu[det_mask*~study],
               xerr=wav_err[det_mask*~study], 
               color=colorize('C1'), ms=9, label='phot')
  ax0.errorbar(phot_wave[~det_mask*~study], 
               phot_nu[~det_mask*~study]*orig[~det_mask*~study],
               fmt='v', yerr=0*err_nu[~det_mask*~study],
               xerr=wav_err[~det_mask*~study], 
               color=colorize('C1'), ms=9, label='phot')
  """
  # plot residuals
  resid = phot_nu*orig-phot_nu*phot_max
  ax1.plot(phot_wave[~det_mask], resid[~det_mask], marker='v', ls='',
           color=colorize('C0'), ms=9)
  ax1.errorbar(phot_wave[det_mask], resid[det_mask], marker='o', ls='',
           color=colorize('C0'), ms=9, yerr=err_nu[det_mask])
  ax1.axhline(y=0, xmin=0, xmax=(phot_wave+wav_err).max()*1.10, ls='--',
              color=colorize('C7'))
  ax1.set_ylim([-(np.abs(resid).max()*1.10), np.abs(resid).max()*1.10])

  # final adjustments
  # ax0.legend(loc='lower right', numpoints=1)
  ax0.set_xscale('log')
  ax1.set_xscale('log')
  ax0.set_yscale('log')
  xmin = (phot_wave-wav_err).min()*0.90
  xmax = (phot_wave+wav_err).max()*1.10
  ax0.set_xlim([xmin, xmax])
  specValid = spec_max[(spec_wave < xmax) & (spec_wave > xmin)]
  specNuValid = spec_nu[(spec_wave < xmax) & (spec_wave > xmin)]
  ymin = (specNuValid*specValid).min()
  ymax = (specNuValid*specValid).max()
  if ymin > (phot_nu*orig).min():
    ymin = (phot_nu*orig).min()
  if ymax < (phot_nu*orig).max():
    ymax = (phot_nu*orig).max()

  ymin*=0.90
  ymax*=2
  ax0.set_ylim([ymin,ymax])

  
  
  #ax0.set_ylim([(phot_nu*phot_max).min()*0.95,
  #              (phot_nu*phot_max).max()*1.05])
  ax0.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
  ax0.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
  ax1.set_xlabel(r'$\lambda_{obs} (\mu m)$', **font_kwargs)
  ax0.set_ylabel(r'$\nu f_{\nu} [erg/s/cm^{2}]$', **font_kwargs)
  ax1.set_ylabel(r'$\chi$', **font_kwargs)

  # add second axis
  zred  = results['obs']['zred']
  ax3   = ax0.twiny()
  ax3.set_xscale('log')
  ax0_x = ax0.get_xticks(minor=True)
  ax3_x = ax0_x/(1+zred)
  ax3.set_xticks(ax0_x)
  ax3.set_xticklabels(ax3_x)
  xbnd = ax0.get_xbound()
  ax3.set_xbound(xbnd[0]/(1+zred), xbnd[1]/(1+zred) )
  ax3.xaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
  ax3.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
  ax3.set_xlabel(r'$\lambda_{rest} (\mu m)$', **font_kwargs)

  # add text
  num_param  = len(results['initial_theta'])
  fit_value  = sum(resid[det_mask]**2)/(sum(det_mask)-num_param)
  fit_string = 'GRB {}'.format(results['obs']['objid']) + '\n'
  fit_string = fit_string + r'best-fit $\chi^2/N_{dof}={x:.2f}$'.\
               format(phot=r'{phot}', x=fit_value, dof='{dof}')
  fit_string = fit_string + '\n' + r'$z={}$'.format(zred)
  ax0.text(0.02, 0.90, fit_string, horizontalalignment='left', verticalalignment='center',
           transform=ax0.transAxes, color='k', **font_kwargs)


  # save the figure
  if save_file is not None:
    sedfig.savefig(save_file, format='eps', transparent=True, pad_inches=0,
              dpi=80)

  return sedfig,ax0,ax1

def spec_range(results, model, obs, sps, num_spec=100, quantiles=[16,84], units='cgs'):
  """Given prospector outputs, get spectra at given quantiles per wavelength.

  :param results:
    A dictionary from prospector.
  
  :param model:
    A dictionary from prospector.

  :param obs:
    A dictionary from prospector. Make sure it has wavelength from sps.

  :param sps:
    A dictionary from prospector.

  :param num_spec: (optional, default: 100)
    The number of spectra from which to get quantiles. If 0, won't do
    anything. If 'all', will do all unburnt chains.

  :param quantiles: (optional, default:[16,84])
    A list determining the quantiles to get (i.e. 16th percentile and 84th
    percentile).
  
  :param units: (optional, default: 'uJy')
    A string denoting whether or not to convert to uJy or to keep things at
    maggies. Passed to model_spec.

  :returns quant_specs:
    An array of the shape of the number of quantiles desired and
    obs['wavelength'].
  """
  
  # if this is totally uncalled for...
  if num_spec==0:
    return None
  elif num_spec=='all':
    iter_burn, walker_noburn = burn_chain(results)
    num_spec = results['lnprobability'].flatten().shape[0] - burn

  # get num_spec flattened indices of highest lnprobability
  iter_burn, walker_noburn = burn_chain(results)
  lnprob  = results['lnprobability'][walker_noburn,iter_burn:]
  
  rand_ind  = random.sample(xrange(lnprob.flatten().shape[0]),num_spec)

  # get model thetas of high lnprobability
  walkers, iterations = np.unravel_index(rand_ind, lnprob.shape)
  thetas = results['chain'][walkers,iterations,:]
  
  # spectra of each set of model thetas
  #specs = np.zeros((num_spec,obs['wavelength'].shape[0]))
  specs = np.zeros(num_spec)
  for i in range(thetas.shape[0]):
    spec, phot = model_spec(thetas[i], model, obs, sps, units=units)
    specs[i] = spec
  return specs 
  quant_specs = np.percentile(specs, quantiles, axis=0)

  return quant_specs

def proper_label(theta_labels):
  """Given theta labels from prospector, get theta labels that are formatted
  prettily for plots.

  :param theta_label:
    List of strings from results['theta_labels']

  :returns proper_label:
    List of strings that are properly formatted for plots. Units will be automatically
    assigned somehow!
  """
  label_dict = {'mass' : r'$\log_{10}\left(M_*/M_{\odot}\right)$',
                'dust2': r'$A_\mathrm{V}$',
                'tage' : r'Age $\mathrm{(Gyr)}$'}

  proper_label = []
  for label in theta_labels:
    try:
      proper_label.append(label_dict[label])
    except KeyError:
      return "Key doesn't exist: {}".format(label)
  
  return proper_label


def plot_corner(myfile,thresh=3, overfig=None, bins=20, extra=None, save_file=None, *hist_kwargs, **hist2d_kwargs):
  """Given an mcmc file, make a dope corner plot. Needs cleaning!

  :param myfile:
    The mcmc prospector file.

  :param thresh: (optional, default: 1)
    Passed to burn_chain.

  :param bins: (optional, default: 20)
    Passed to corner to determine the number of bins to use for the 1-D
    histogram plots.

  :param extra: (optional, default: None)
    A dictionary specifying whether cuts should be made and displayed in the
    histograms. Keys should describe the theta label and a float
    above which and below which to cut.

  :param save_file: (optional, default: None)
    If not None, will save the figure to the string specified by save_file.
  
  :param hist_kwargs: (optional)
    Additional kwargs passed to corner for 1-D histogram plots.
    -visible
    -linewidth
    -alpha
    -histtype

  :param hist2d_kwargs: (optional)
    Additional kwargs passed to corner.hist2d for 2-D historgam plots.
    -color
    -plot_datapoints
    -plot_density
    -plot_contours
    -no_fill_contours
    -data_kwargs
      -visible
      -alpha
  """
  results, obs, model = rr.results_from(myfile)

  # get burned chain and reshape
  iter_burn, walker_noburn  = burn_chain(results)
  chain = results['chain'][walker_noburn,iter_burn:,:]
  walk_sum, iter_sum, theta_sum = chain.shape
  shapely_chain = chain.reshape(walk_sum*iter_sum, theta_sum)

  # get 'good' units for thetas
  theta_labels = results['theta_labels']
  for idx, label in enumerate(theta_labels):
    if label == 'mass':             # use log(mass)
      mass_changed = True
      shapely_chain[:,0] = np.log10(shapely_chain[:,0])
    if label == 'dust2':
      #dust_changed = True
      shapely_chain[:,1] /=0.921
  #shapely_chain = shapely_chain[shapely_chain[:,2]<0.5,:]
  
  # get 'good' labels for thetas
  proper_labels = proper_label(theta_labels)

  myrange = []      # to directly specify ranges for each label in cornerfig
  for i in range(shapely_chain.shape[1]):
    mymin = np.min(shapely_chain[:, i])
    mymax = np.max(shapely_chain[:, i])
    myrange.append((mymin,mymax))
  #one_sigma_2d = 1-np.exp(-0.5)
  levels = (0.39, 0.86, 0.99)

  # handle extras
  if extra is not None:
    try:
      idx = theta_labels.index(extra['label'])
    except ValueError:
      print "Extra doesn't look right."
    else:
      extra['index'] = idx
      if overfig is None:
        cornerfig = plot_corner_extra(shapely_chain, extra, bins, myrange)
      else:
        cornerfig = plot_corner_extra(shapely_chain, extra, bins, myrange)

      # update kwargs for primary plot
      try:
        hist2d_kwargs.keys()
      except NameError:
        hist2d_kwargs = {}

      # plot "main" histogram
      data_kwargs = {'visible': False}
      hist2d_kwargs['data_kwargs'] = data_kwargs
      cornerfig = corner.corner(shapely_chain, quantiles=[0.16, 0.5, 0.84],
                                fig = cornerfig,
                                bins=bins, range=myrange, color="k",
                                labels=proper_labels, show_titles=True, levels=levels,
                                no_fill_contours=True, plot_density=False,
                                *hist_kwargs, **hist2d_kwargs)

      # plot SED with extra labels
      lnprobability = results['lnprobability'][walker_noburn,iter_burn:].flatten()
      shapely_chain = np.column_stack((shapely_chain, lnprobability))
      sps        = rr.get_sps(results)
      theta_pair = theta_cut(shapely_chain, extra=extra)
      # make sure to pass proper values of mass
      if mass_changed==True:
        theta_pair[:,0] = np.power(10,theta_pair[:,0])
      sed_plot   = plot_max_sed(myfile, extra_thetas=theta_pair, num_spec=0)

  # just your everyday corner plot
  else:
    if overfig is None:
      cornerfig = corner.corner(shapely_chain, quantiles=[0.16, 0.5, 0.84],
                                bins=bins, range=myrange, color="k",
                                labels=proper_labels, show_titles=True, levels=levels,
                                *hist_kwargs, **hist2d_kwargs)
    else:
      cornerfig = corner.corner(shapely_chain, quantiles=[0.16, 0.5, 0.84],
                                bins=bins, range=myrange, color="k",
                                labels=proper_labels, show_titles=True, levels=levels,
                                figure=overfig,*hist_kwargs, **hist2d_kwargs)

  # Add texts
  GRB, z = obs['objid'], obs['zred']
  cornerfig.suptitle('GRB {GRB}\nz={z}'.format(GRB=GRB, z=z), x=0.95, y=0.95,
                     horizontalalignment='right', fontsize=16)
  # save the figure
  if save_file is not None:
    cornerfig.savefig(save_file, format='png', transparent=True, pad_inches=0,
                      dpi=80)
    plt.close(cornerfig)

  return cornerfig

def plot_corner_extra(shapely_chain, extra, bins, myrange):
  """Plot corner plot if there's a cut in the data.

  :param shapely_chain:
    An array of MxN dimension, where M is the iteration/walker and N is the
    theta.
  
  :param extra:
    A dictionary containing the keys index and value, specifying the theta and
    the value at which to cut.

  :param bins:
    The number of bins to use in the histogram. Used to make sure the histogram
    makes 'sense'.

  :param myrange:
    The span of the histogram, same usage as bins.

  :returns cornerfig:
    A figure onto which further fun things can be plotted. Shows corner plots in
    two different colors that correspond to both sides of the cut that was made.
  """

  idx = extra['index']
  val = extra['value']

  # determine where values lie above cutoff
  chain_mask = shapely_chain[:,idx] > val
  temp_chain = dc(shapely_chain)

  # start plotting above cutoff
  cornerfig = corner.corner(shapely_chain[chain_mask],
                            color=colorize('C2'), bins=bins, range=myrange,
                            show_titles=False, plot_contours=False,
                            no_fill_contours=True, plot_density = False,
                            hist_kwargs={'alpha': 0.4, 'histtype': 'stepfilled'},
                            data_kwargs={'alpha': 0.2})

  # plot below cutoff
  cornerfig = corner.corner(shapely_chain[~chain_mask], fig=cornerfig,
                            color=colorize('C3'), bins=bins, range=myrange,
                            show_titles=False, plot_contours=False,
                            no_fill_contours=True, plot_density = False,
                            hist_kwargs={'alpha': 0.4, 'histtype': 'stepfilled'},
                            data_kwargs={'alpha': 0.2})

  return cornerfig
  
def theta_cut(shapely_chain, extra, typical_cut=50):
  """Grabs theta for an SED showing a cut specified by extra. Given both sides of the cut,
  this will find typical theta values based on the lnprobability.

  :param shapely_chain:
    An array with walker/iteration and theta + lnprobability.

  :param extra:
    A dictionary specifying whether cuts should be made. See above.

  :param cut_type: (optional, default: 50)
    A float determining how to use lnprobability to get typical theta values. 
    -50: use the median lnprobability in either regime
    -100: use the maximum lnprobability in either regime

  :returns theta_pair:
    An array of a pair of model thetas that are typical of both sides of the
    cut. The first set is above the cut, second set is below the cut.
  """

  idx = extra['index']
  val = extra['value']
  
  # determine where values lie above cutoff
  chain_mask = shapely_chain[:,idx] > val
  
  # determine index of typical representatives
  upper = np.percentile(shapely_chain[:,-1][chain_mask] , typical_cut)
  upper_index  = np.abs(shapely_chain[:,-1][chain_mask]  - upper).argmin()
  lower = np.percentile(shapely_chain[:,-1][~chain_mask], typical_cut)
  lower_index  = np.abs(shapely_chain[:,-1][~chain_mask] - lower).argmin()

  # grab the theta of typical representatives
  theta_pair = np.array((shapely_chain[chain_mask][upper_index,:-1],
                        shapely_chain[~chain_mask][lower_index,:-1]))

  return theta_pair

def spec_range_tester(myfile,units='cgs'):
  results, obs, model = rr.results_from(myfile)
  sps = rr.get_sps(results)
  obs['wavelength'] = sps.wavelengths[5013]

  flatchain   = flatten_results(myfile, thresh=3)
  total_theta = flatchain.shape[0]
  all_specs   = np.zeros(total_theta)
  
  for idx,theta in enumerate(flatchain):
    all_specs[idx], phot = model_spec(theta[:-1], model, obs, sps, units=units)
  end_num = total_theta / 100 * 100
  num_specs = np.arange(100,end_num,100)
  fileWrite = '041006_std_tester.txt'
  fileOpen  = open(fileWrite,'w')
  f, (ax1, ax2) = plt.subplots(2,1)
  for num_spec in num_specs:
    stds = np.zeros(4)
    spec_means = np.zeros(4)
    for i in range(4):
      rand_ind  = random.sample(xrange(total_theta), num_spec)
      specs     = all_specs[rand_ind]
      stds[i]   = np.std(specs)
      spec_means[i] = np.mean(specs)
    spec_mean = np.mean(spec_means)/1e-29
    std_stds = np.std(stds)/1e-29
    mean_stds = np.mean(stds)/1e-29
    label = r'{} {} {:.2f}$\pm${}$'.\
            format(num_spec,spec_mean, mean_stds, std_stds)
    fileOpen.write(label+'\n')
    ax1.errorbar(num_spec,spec_mean,yerr=mean_stds)
    ax2.errorbar(num_spec,mean_stds,yerr=std_stds)

  ax1.set_xlabel('Number of Spectra Sampled')
  ax2.set_xlabel('Number of Spectra Sampled')
  ax1.set_xlim((0,end_num+100))
  ax2.set_xlim((0,end_num+100))
  ax1.set_ylabel('Mean of Spectra')
  ax2.set_ylabel('Mean of $\sigma$')
  fileOpen.close()
  plt.savefig('041006_std_tester.eps', format='eps', transparent=True, pad_inches=0,
              dpi=80)

def quantilize(shapely_chain):
  thetas  = np.percentile(shapely_chain,[50,84,16],axis=0)
  thetas  = np.array(thetas)
  uppers  = thetas[1,:]-thetas[0,:]
  lowers  = thetas[0,:]-thetas[2,:]
  final_theta = np.array([thetas[0,:],uppers,lowers])
  return final_theta

def plot_burn_chain(myfile, x=3, order='C'):
  results,obs,model =rr.results_from(myfile)
  lnprobability = results['lnprobability']
  
  plt.figure()
  for i in range(np.shape(lnprobability)[0]):
    plt.plot(lnprobability[i,:],'-', alpha=0.5,color='gray')
  total_iter  = np.shape(lnprobability)[1]
  favored_few = lnprobability[:,-int(0.10*total_iter)+1:]
  favored_mean  = np.mean(favored_few,axis=1)
  sigma = np.std(favored_mean,axis=0)
  mean  = np.mean(favored_mean,axis=0)
  
  all_mean  = np.mean(lnprobability,axis=0)
  iter_burn, walker_noburn = burn_chain(results,thresh=x)
  
#  burn  = np.where(all_mean>mean-x*sigma)[0][0]
  plt.vlines(iter_burn,np.min(lnprobability),np.max(lnprobability))  
  for i in range(lnprobability.shape[0]):
    if i in walker_noburn: continue
    else: plt.plot(lnprobability[i,:],'r-')
  return iter_burn, walker_noburn

def colorize(CN):
  color_dict = {'C0': '#1f77b4',
                'C1': '#ff7f0e',
                'C2': '#2ca02c',
                'C3': '#d62728',
                'C4': '#9467bd',
                'C5': '#8c564b',
                'C6': '#e377c2',
                'C7': '#7f7f7f',
                'C8': '#bcbd22',
                'C9': '#17becf'}

  try:
    return color_dict[CN]
  except KeyError:
    print "No can do, bud!"
    return 'k'
