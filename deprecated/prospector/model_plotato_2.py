import numpy as np
import prospect.io.read_results as rr
import random
import model_utils
import matplotlib.pyplot as plt
import corner
from copy import deepcopy as dc

class Result(object):
  """For analyzing results data.
  :param myfile:
    The name and path to output prospector file ending in 'mcmc'
  """
  def __init__(self, myfile):
    # get results, obs, model, and sps from the mcmc file
    self.results, self.obs, self.model  = rr.results_from(myfile)
    self.sps  = rr.get_sps(self.results)

    # unpack some of the more used params
    self.niter        = self.results['run_params']['niter']
    self.nwalkers     = self.results['run_params']['nwalkers']
    self.lnprob       = self.results['lnprobability']
    self.chain        = self.results['chain']

    # ensure that wavelengths are properly defined
    self.obs['wavelength'] = self.sps.wavelengths
  
    # ensure that lnprob are finite
    badWalkers  = np.unique(np.where(np.isfinite(self.lnprob) == False)[0])
    if badWalkers.shape[0] > 0:
      self.nwalkers = self.nwalkers - badWalkers.shape[0]
      self.lnprob   = np.delete(self.lnprob, badWalkers, 0)
      self.chain    = np.delete(self.chain, badWalkers, 0)

  def burn_chain(self, thresh=3, order='C'):
    """Read the lnprobability of the mcmc run, get the median of the last 10% of
    every walker's iterations (primary median) and find the iteration at which 
    the median of the walkers at that iteration first falls above the primary
    median minus some delta, where delta is the standard deviation of the
    primary median times thresh. 

    Will also burn walkers with medians past the burned iterations that don't
    reach the primary median within the same delta.

    :param thresh: (optional, default: 3)
      The factor by which to multiply the standard deviation of the primary
      median.
  
    :param order: (optional, default: 'C')
      Passed to numpy.reshape to determine order of rows. 'C' will order by
      walker and iteration, whereas 'F' will order by iteration and walker.

    :returns iterBurn:
      The index below which to burn iterators.

    :returns walkerNoBurn:
      Indices specifying walkers that will not be burned.
    """
    # use only the last 10% of iterations from each walker
    iterFew = self.lnprob[:, -int(.1*self.niter):]

    walkerMedian  = np.median(iterFew, axis=1)        # median of walkers at each iteration
    totalMedian   = np.median(walkerMedian, axis=0)   # median of all walkers
    walkerSigma   = np.std(walkerMedian, axis=0)

    # compare to the median of all walkers per iteration
    iterMedian  = np.median(self.lnprob, axis=0)      # median of iterations at each walker

    # retrieve index where median over iterations exceeds (median - delta)
    iterBurn  = np.where(iterMedian > totalMedian - (thresh*walkerSigma) )[0][0]
    
    # retrieve unburned chains and keep only certain walkers
    leftovers     = np.median(self.lnprob[:, iterBurn:], axis=1)
    walkerNoBurn  = np.where(leftovers > totalMedian - (thresh*walkerSigma))[0]
   
    #iterBurn=2000. 
    return iterBurn, walkerNoBurn

  def plot_burn_chain(self, index=('mass','dust2','tage','lnprob'), nrows=2, 
                      burnVisible=True, plotStd=True, **burnParams):
    """Crude diagnostic to check convergence of chains.
    
    :param index:
      Either an index of a flattened chain specifying what particular parameter
      to look at, or a sequence of labels that correspond to all of the indices
      if you want to look at all of them. Ability to only pick out specific
      indices not implemented. Also only really tested for four parameters,
      mass, dust2, age, and lnprob.

    :param burnParams: (optional)
      Additional arguments to send to the burn function.

    :param nrows: (optional, default: 2)
      The number of rows desired in the subplot. Limited testing done on this,
      might misbehave. Will reshape the number of columns automatically
      depending on how much data there are. Will also delete empty axes, so if
      you specify more rows than the data will fit, prepare to have one wonky
      figure.

    :param burnVisible: (optional, default: True)
      Show the burned chains or not.

    :returns burnfig:
      The figure that helpfully shows the parameter space explored by each
      walker during each iteration, and where the chains were cut.
    """
    # create the full flattened chain and get certain run_param info
    flatChain = self.flatten_results(thresh=np.inf)
    iterBurn, walkerNoBurn = self.burn_chain(**burnParams)

    # initialize plot kwargs
    kwargsIntact = {'color': 'gray', 'visible': True, 'zorder': 0, 'alpha': 0.2}
    kwargsBurned = {'color': 'red', 'visible': burnVisible, 'zorder': 1, 'alpha': 1.0}

    allIters = np.arange(0, flatChain.shape[0], self.niter)
    itersstd  = np.zeros(self.niter)

    # if only probing one index
    if isinstance(index, int):
      for i in range(self.niter):
        allwalkers  = flatChain[allIters+i,index]
        walkerstd = np.std(allwalkers)
        itersstd[i] = walkerstd
      burnfig, ax = plt.subplots(1,1)
      for walker in range(0, self.nwalkers):
        walkerFlatIndex = walker * self.niter  # assuming chain ordered 'C'
        if walker in walkerNoBurn:
          kwargsUsage = kwargsIntact
        else:
          kwargsUsage = kwargsBurned
        if plotStd == True:
          ax.plot(itersstd, '-', **kwargsUsage)
        else:
          ax.plot(flatChain[walkerFlatIndex: walkerFlatIndex+self.niter, index],
                  '-', **kwargsUsage)
      ax.vlines(iterBurn, 0.9*flatChain[:, index].min(),
                          1.1*flatChain[:, index].max(), **kwargsUsage)

    else:
      # declare subplot things
      if nrows == 1:
        ncols = len(index)
      else:
        ncols  = (len(index) + nrows - 1)/nrows % (len(index) + nrows - 1)
      burnfig, (axes) = plt.subplots(nrows, ncols)
      
      for i in range(len(index)):
        for j in range(self.niter):
          allwalkers  = flatChain[allIters+j,i]
          walkerstd   = np.std(allwalkers)
          itersstd[j] = walkerstd
        # assign proper axes to indices
        if (ncols == 1) or (nrows == 1):
          ax = axes[i]
        else:  
          row = (i/ncols) % nrows
          col = i % ncols
          ax = axes[row, col]

        for walker in range(0, self.nwalkers):
          walkerFlatIndex = walker * self.niter
          if walker in walkerNoBurn:
            kwargsUsage = kwargsIntact
          else:
            kwargsUsage = kwargsBurned
          if plotStd == True:
            ax.plot(itersstd, '-', **kwargsUsage)
          else:
            ax.plot(flatChain[walkerFlatIndex: walkerFlatIndex+self.niter, i],
                    '-', **kwargsUsage)

        # axis logistics
#        ymin = 0.9*flatChain[:, i].min()
#        ymax = 1.1*flatChain[:, i].max()
#        ax.vlines(iterBurn, ymin, ymax, **kwargsBurned)
        ax.set_xlim(0, self.niter)
#        ax.set_ylim(ymin, ymax)
        ax.set_ylabel(index[i])
      # delete empty axes if there are more axes than data
      if nrows * ncols > 0:
        extra = (nrows * ncols) - len(index)
        flat  = axes.flatten()
        for i in np.arange(extra) + 1:
          burnfig.delaxes(flat[-i])

    print iterBurn
    return burnfig

  def flatten_results(self, noln=False, **burnArgs):
    """Flatten chain from the results.

    :param burnArgs:
      Keyword arguments passed to burn_chain to specify how the burn proceeds.

    :param noln:
      If True, will not include the lnlikelihood in the chain.

    :returns chainFlat:
      An array containing all of the walker and iteration model values,
      where each column is a theta, with the last column as the corresponding
      lnprobability. By default, the rows are ordered by walkers and their
      respective iterations.
    """
    # burn the chain
    iterBurn, walkerNoBurn = self.burn_chain(**burnArgs)

    # get dimensions of flattened chain
    rows  = walkerNoBurn.shape[0] * (self.niter - iterBurn)
    cols  = len(self.results['theta_labels'])

    # flatten the chain with burn
    chainFlat = np.zeros((rows, cols))
    chainFlat = self.chain[walkerNoBurn, iterBurn:, :].reshape(rows, cols)

    if noln==False:
      chainLnprob = self.lnprob[walkerNoBurn, iterBurn:].reshape(rows)
      chainFlat = np.column_stack((chainFlat, chainLnprob))

    return chainFlat

  def quantilize(self, quants=[50,84,16], basic=False, **burnArgs):
    chainFlat = self.flatten_results(noln=True, **burnArgs)
    thetas    = np.array(np.percentile(chainFlat, quants, axis=0))
    if basic == True:
      return thetas
    else:
      uppers  = thetas[1, :]-thetas[0, :]
      lowers  = thetas[0, :]-thetas[2, :]
      final_thetas = np.column_stack((thetas[0,:], uppers, lowers))
      return final_thetas

  def theta_max(self):
    """Get thetas that correspond to the maximum likelihood.

    :returns thetaMax:
      An array with the maximum theta.
    """
    indMax  = self.lnprob.argmax()
    walker, iteration = np.unravel_index(indMax, self.lnprob.shape)
    thetaMax  = self.chain[walker, iteration, :]

    return thetaMax

  def model_spec(self, theta, units='cgs'):
    """Given theta values, get the modeled spectrum and synthetic
    photometry. In hindsight, probably doesn't need its own function,
    hindsight being what it is...
  
    :param theta:
    if synth_phot == True:
      
      Desired model parameters.

    :param units: (optional, default: 'cgs')
      Desired units of the output spectrum and photometry.

    :returns specScaled:
      The SED that corresponds to wavelengths in self.obs['wavelength'].

    :returns photScaled:
      The photometry in each observed filted integrated over the spectrum.
    """
    spec, phot, x = self.model.mean_model(theta, obs=self.obs, sps=self.sps)
    specScaled, photScaled = model_utils.scale_units(units, spec, phot)

    return specScaled, photScaled

  def obs_phot(self, unitsPhot='cgs', unitsWave='um'):
    """Get the photometry data.

    :param unitsPhot: (optional, default: 'cgs')
      String specifying to what unit to convert the output spectrum and
      photometry.

    :param unitsWave: (optional, default: 'um')
      String specifying to what unit to convert the output wavelengths.

    :returns photOrigScaled:
      An array of the observed photometry.

    :returns photErrScaled:
      An array of the uncertainties of the observed photometry. In the case of 
      upper limits, will return the 3-sigma upper limits.

    :returns waveErrScaled:
      An array of the error bars for each photometric filter.
    """
    # initialize data
    detMask   = self.obs['det_mask']
    photOrig  = self.obs['maggies']
    photErr   = self.obs['maggies_unc']
    wave      = np.array([ f.wave_effective for f in self.obs['filters'] ])
    waveErr   = np.array([ f.rectangular_width for f in self.obs['filters'] ])

    # consequence of how obs was constructed
    photErr[~detMask]  = self.obs['maggies'][~detMask]

    # rescale
    photOrigScaled, photErrScaled = model_utils.scale_units(unitsPhot, photOrig, photErr)
    waveScaled, waveErrScaled = model_utils.scale_units(unitsWave, wave, waveErr)

    return photOrigScaled, photErrScaled, waveScaled, waveErrScaled

  def spec_range(self, numSpec=0, quantiles=[16,84], units='cgs', **burnParams):
    """Given prospector outputs, get spectra at given quantiles per wavelength.

    :param numSpec: (optional, default: 5000)
      The number of spectra from which to get quantiles. If 0, won't do
      anything. If 'all', will do all unburnt chains.

    :param quantiles: (optional, default: [16,84])
      A list determining the quantiles to calculate (i.e. 16th percentile and
      84th percentile).

    :param units: (optional, default: 'cgs')
      A string denoting what to convert the spectrum to. Passed to model_spec.

    :returns specsQuant:
      An array of the shape of the number of quantiles desired and
      sps.wavelength
    """
    # build chain
    iterBurn, walkerNoBurn = self.burn_chain(**burnParams)
    lnprob  = self.lnprob[walkerNoBurn, iterBurn:]
    if numSpec == 0:
      return None
    elif numSpec == 'all':
      numSpec = self.lnprob.flatten().shape[0]

    # get indices and unpack thetas
    indRand = random.sample(range(lnprob.flatten().shape[0]), numSpec)
    walkers, iterations = np.unravel_index(indRand, lnprob.shape)
    thetas  = self.chain[walkers, iterations, :]

    # get spectra
    specs = np.zeros((numSpec, self.sps.wavelengths.shape[0]))
    for idx, theta in enumerate(thetas):
      spec, phot  = self.model_spec(theta, units=units)
      specs[idx]  = spec

    specsQuant  = np.percentile(specs, quantiles, axis=0)
    return specsQuant

  def plot_max_sed(self, overfig=None, savefile=None, thetaExtra=None,
                   numSpec=0, synth_phot=False, **kwargs):
    """Plot the spectrum with photometry and residual info. Will presumably save
    figure files sometime. Will support phot_mask where, confusingly, True means
    the data point was used.

    :param overfig: (optional, default: None)
      If not None, will attempt to plot the figure to be plotted onto another
      figure. Might not behave well.

    :param savefile: (optional, default: None)
      If not None, will save the figure to the string specified.

    :param thetaExtra: (optional, default: None)
      For a cut in the SED.

    :param kwargs: (optional)
      A dictionary to pass to spec_range and burn_chain, if you don't like the
      defaults. Could probably be implemented better.
    """
    # specify some kwargs for fonts, can epxand later
    font_kwargs  = {'fontsize': 16, 'color': 'k'}

    # get photometric data and model data
    data, dataSig, dataWave, dataWaveSig  = self.obs_phot()
    spec, phot  = self.model_spec(self.theta_max())

    # get wavelength data for spectrum. Can't be helped! Then turn to Hz
    specWave  = model_utils.scale_units('um', self.obs['wavelength'])
    specHz, photHz = model_utils.scale_units('Hz', specWave*1e4, dataWave*1e4)

    # Get n*f_nu for all the things
    specNu, photNu    = specHz*spec, photHz*phot
    dataNu, dataSigNu = photHz*data, photHz*dataSig

    # Plotting time! Declare some things first
    sedfig, ax0, ax1  = model_utils.sed_figure_init(overfig)

    # Plot spectra
    if overfig == None: 
      specColor = model_utils.colorize('C0')
    else: 
      specColor = model_utils.colorize('C3')
    ax0.plot(specWave, specNu, color=specColor)
    # Handle extra thetas
    if thetaExtra != None:
      for i in range(len(thetaExtra)):
        specColor = model_utils.colorize( 'C{}'.format(i+2) )
        specExtra, photExtra = self.model_spec(thetaExtra[i])
        ax0.plot(specWave, specHz*specExtra, color=specColor)

    if numSpec != 0:
      if overfig == None:
        specColor = model_utils.colorize('C7')
      else:
        try:
          specColor = model_utils.colorize(kwargs['CN'])
        except KeyError:
          specColor = model_utils.colorize('C8')
      specsQuant    = self.spec_range(numSpec=numSpec)
      specsQuantNu  = specHz*specsQuant
      ax0.fill_between(specWave, specsQuantNu[0], specsQuantNu[1],
                       interpolate=True, alpha=0.5, color=specColor)
    
    # Plot photometry, and visually separate data from this study
    dataMine = np.array([ i.startswith('my_spitzer') for i in self.obs['filternames']])
    detMask  = self.obs['det_mask']
    dataSigNu[~detMask] *= 0
    for i in range(photNu.shape[0]):
      if self.obs['phot_mask'][i] != True: continue # skip if masked in fitting

      if dataMine[i] == True: 
        photColor = model_utils.colorize('C3')
      else: 
        photColor = model_utils.colorize('C1')

      if detMask[i] == True: 
        photMrk = 'o'
      else: 
        photMrk = 'v'
      
      ax0.errorbar(dataWave[i], dataNu[i], yerr=dataSigNu[i], xerr=dataWaveSig[i],
                   marker=photMrk, color=photColor, ms=9)
    if synth_phot == True:
      ax0.plot(dataWave, photNu, 'ko')
      

    # Plot chi
    chi = (dataNu - photNu)/dataSigNu
    ax1.plot(dataWave[detMask], chi[detMask],'o', 
             color=model_utils.colorize('C0'), ms=9)
    #ylim = np.abs(chi[detMask]-dataSigNu[detMask]).max()*1.10
    #ax1.set_ylim((-ylim, ylim))

    # Final adjustments
    xmin, xmax = model_utils.x_getter(dataWave, dataWaveSig)
    #ymin, ymax = model_utils.y_getter(photNu, xmin, xmax, specNu, specWave)
    ymin, ymax = (dataNu-dataSigNu).min()/3., (dataNu+dataSigNu).max()*3
    ax1.axhline(y=0, xmin=0, xmax=xmax, ls='--', color=model_utils.colorize('C7'))
    ax0.set_xlim([xmin,xmax])
    ax0.set_ylim([ymin,ymax])
    ax0, ax1, ax2 = model_utils.axes_formatter(ax0, ax1, self.obs['zred'])
    # Add text information
    numParams = len(self.results['initial_theta'])
    chiNdof   = sum(chi[detMask]**2) / (sum(detMask*self.obs['phot_mask']) - numParams)
    if (np.isinf(chiNdof)) or (chiNdof < 0):
      chiNdof = sum(chi[detMask]**2)
      chiNdofStr  = 'GRB {}\n'.format(self.obs['objid']) +\
                   '$\chi^2={x:.2e}$\n'.\
                   format(x=chiNdof) +\
                   '$z={}$'.format(self.obs['zred'])
    else:
      chiNdofStr  = 'GRB {}\n'.format(self.obs['objid']) +\
                   '$\chi^2/\mathrm{N}_\mathrm{dof}={x:.3f}$\n'.\
                   format(N=r'{N}', dof=r'{dof}', x=chiNdof) +\
                   '$z={}$'.format(self.obs['zred'])
    ax0.text(0.02, 0.90, chiNdofStr, 
             horizontalalignment='left', verticalalignment='center',
             transform=ax0.transAxes, **font_kwargs)
    # Label axes
    ylabel0 = model_utils.axes_labeler('cgs')
    ylabel1 = model_utils.axes_labeler('chi')
    xlabel1 = model_utils.axes_labeler('umobs')
    xlabel2 = model_utils.axes_labeler('umrest')
    ax0.set_ylabel(ylabel0, **font_kwargs)
    ax1.set_ylabel(ylabel1, **font_kwargs)
    ax1.set_xlabel(xlabel1, **font_kwargs)
    ax2.set_xlabel(xlabel2, **font_kwargs)

    if savefile != None:
      sedfig.savefig(savefile, format='pdf', transparent=True, pad_inces=0, dpi=80)
      plt.close(sedfig)
    else:
      return sedfig, ax0, ax1

  def plot_corner(self, overfig=None, savefile=None, thetaCut=None, thresh=3,
                  cornerBins=20, color='k', agelogify=True, 
                  *hist_kwargs, **hist2d_kwargs):
    """Given an mcmc file, make a corner plot. It's nice.
  
    :param overfig: (optional, default: None)
      If plotting extras, supply overfig. Not necessary unless performing
      multiple cuts. Deprecated.

    :param savefile: (optional, default: None)
      If not None, will save the figure to the string specified.

    :param thetaCut: (optional, default: None)
      A dictionary specifying whether cuts should be made and displayed in the
      histograms. Kets should describe the theta label and a float at which to
      cut.
  
    :param thresh: (optional, default: 3)
      For burning purpsoes. Please see burn_chain.

    :param cornerBins: (optional, default: 20)
      Bins for the histogram plots, for consistency.

    :param hist_kwargs: (optional)
      Additional kwargs passed to corner for 1-D histogram plots.
      -visible
      -linewidth
      -alpha
      -histtype

    :param hist2d_kwargs: (optional)
      Additional kwargs passed to corner.hist2d for 2-D histogram plots.
      -color
      -plot_datapoints
      -plot_density
      -plot_contours
      -data_kwargs
        -visible
        -alpha
    """
    # burn chain and reshape
    chainFlat = self.flatten_results(thresh=thresh)
    thetaFlat = dc(self.results['theta_labels'])
    
    # 'fix' units for thetas
    thetaLabels, thetaChain = model_utils.thetas_labeler(thetaFlat, chainFlat,
                                                         agelogify=agelogify)
 
    # cornerfig specifications
    cornerRange = zip(thetaChain.min(axis=0), thetaChain.max(axis=0))
    cornerLevel = (0.39, 0.86, 0.99)  #2d sigmas
    quantiles   = (0.16, 0.5, 0.84)   #1d sigmas

    # update kwargs for "primary, full" plot
    try:
      hist2d_kwargs.keys()
    except NameError:
      hist2d_kwargs = {}
    # declare a few formatting details
    linecolor = (model_utils.colorize('C9'), 
                 model_utils.colorize('C3'),
                 model_utils.colorize('C9'))
    contour_kwargs = {'colors': ('k',
                                 model_utils.colorize('C1'),
                                 model_utils.colorize('C9')),
                      'linewidths': 2}
    contourf_kwargs = {'colors': ('0.9','0.5','0.3','0.0')}
    corner_kwargs = {'fig': overfig, 'quantiles': quantiles, 'bins': cornerBins,
                     'range': cornerRange[:-1], 'levels': cornerLevel,
                     'labels': thetaLabels, 'show_titles': True,
                     'no_fill_contours': False, 'plot_density': False,
                      'fill_contours': False, 'contourf_kwargs': contourf_kwargs,
                     'linecolor': linecolor, 'contour_kwargs': contour_kwargs,
                     'label_kwargs': {'fontsize':14}, 
                     'title_kwargs': {'fontsize': 14}}

    # handle extras
    if thetaCut is not None:
      cornerfig = model_utils.plot_corner_extra(chainFlat, thetaCut, cornerBins, cornerRange)

      # plot "main" histogram without 2d data
      corner_kwargs['data_kwargs']    = {'visible': False}
      corner_kwargs['no_fill_contours'] = True
      corner_kwargs['plot_density']     = False
      corner_kwargs['fig']    = cornerfig
      cornerfig = corner.corner(thetaChain[:,:-1], *hist_kwargs, **corner_kwargs)

      # plot SED with extra labels
      thetaPair = model_utils.theta_cut(chainFlat, thetaCut)
      # pass proper mass values
      thetaPair[:,0]  = np.power(10, thetaPair[:,0])
      sedPlot = self.plot_max_sed(thetaExtra=thetaPair, numSpec=0)

    # just your everyday corner plot, no cuts, possibly an overfig
    else:
      #cornerfig = corner.corner(chainFlat[:, :-1], *hist_kwargs, **corner_kwargs)
      truths  = dc(self.theta_max())
      truths[0] = np.log10(truths[0])
      truths[1] /=  0.921
      truths[2] = np.log10(truths[2]*1e9)
      cornerfig = corner.corner(chainFlat[:, :-1], truths=truths, *hist_kwargs,
                                **corner_kwargs)
    # add texts
    GRB, z = self.obs['objid'], self.obs['zred']
    cornerfig.suptitle('GRB {GRB}\nz={z}'.format(GRB=GRB, z=z),
                       x=0.95, y=0.95, horizontalalignment='right',
                       fontsize=20)

    # save the figure
    if savefile is not None:
      cornerfig.savefig(savefile, format='pdf', transparent=True, 
                        pad_inches=0, dpi=80)
      plt.close(cornerfig)
    else:
      return cornerfig
