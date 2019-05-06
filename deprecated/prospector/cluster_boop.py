import numpy as np
import os, subprocess
import model_plotato_2 as mplo
from multiprocessing import Process
from multiprocessing.dummy import Pool as ThreadPool

def init_run_params(**run_params):
  """Creates a string by which to run my_prospector.py using python.
  :param run_params: (optional)
    A dictionary of keyword arguments to pass into my_prosector.
    Things you might use:
      nwalkers
      niter
      param_file
      objid
      outfile
      phottable
  """
  initStr = ['python','my_prospector.py']
  for key, val in run_params.items():
    addStr  = '--{key}={val}'.format(key=key, val=val)
    initStr.append(addStr)

  runStr  = ' '.join(initStr)
  return runStr

def init_model_params(**model_params):
  """Creates a dictionary by which to initialize model parameters as befits a
  SSP with a calzetti dust law.
  :param model_params: (optional)
    A dictionary of keyword arguments to use for the model and the dispersion. 
    Defaults are as follows:
      init_tage: 1 (Gyr)
      init_disp_tage: 1 (Gyr)
      init_mass: 1e9 (solar masses)
      init_disp_mass: 1e9 (solar masses)
      init_dust2: 3.5 (must remember why)
      init_disp_dust2: 3.5 (see above)

  :returns template_params:
    A dictionary with the key and value pairs of desired model parameters.
  """
  # age is now dependent on redshift
  template_params  = {'init_mass': 5e10, 'init_disp_mass': 5e10,
                      'init_dust2': 3.5, 'init_disp_dust2': 3.5}
  """
  template_params  = {'init_tage': 3.5, 'init_disp_tage': 3.5,
                      'init_mass': 5e10, 'init_disp_mass': 5e10,
                      'init_dust2': 3.5, 'init_disp_dust2': 3.5}
  """
  for key, val in model_params.items():
    template_params[key] = val
  return template_params

def seek_params(key, val, templateLines):
  """Given a template param file, replace model parameters with provided values.
  :param key:
    A string indicating the parameter being replaced.

  :param val:
    A value indicating the parameter value to use. Can't get any clearer.

  :param templateLines:
    The template stored as a string.

  :returns templateLines:
    Modified input.
  """
  for i,j in enumerate(templateLines):
    if key in j:
      templateLines[i] = templateLines[i].replace(key, str(val))
  return templateLines

def model_boop(outfile, infile='template_params.py', **model_params):
  """Writes a parameter file to be accepted as an input to my_prospector.py.
  :param outfile:
    The desired filename of the parameter file.

  :param infile: (optional, default: 'template_params.py')
    The template parameter file.

  :param model_params: (optional)
    Parameters to change in the template file. See init_model_params for
    defaults.
  """
  template_params = init_model_params(**model_params)
  templateLines = open(infile, 'r').readlines()
  for key, val in template_params.items():
    templateLines = seek_params(key, val, templateLines)

  with open(outfile, 'w') as f:
    f.write(''.join(templateLines))

def get_new_params(objid, run, outdir='cluster_run'):
  """Using my very special file structure, read a parameter file and get new
  model parameters using the median values as initial guesses and the standard
  deviations as dispersion widths.

  :param objid:
    The object name specified in the photometry data file.

  :param run:
    Run number as an integer unique to my file structure. Used to generate the 
    parameter file for reading.

  :param outdir: (optional, default: 'cluster_run')
    The output directory from which the parameter file is to be read. E.G.
    'cluster_run/objid_runnum_blah_mcmc'

  :returns model_params:
    A dictionary specifying the new init values and their respective dispersion.
  """
  paramsDict  = {0: 'mass', 1: 'dust2', 2: 'tage'}
  filesAll    = os.listdir('{outdir}/{objid}'.format(outdir=outdir, objid=objid))
  start = '{objid}_run{num}'.format(objid=objid, num=run)
  mcstr = [f for f in filesAll if (f.startswith(start)) and \
          (f.endswith('mcmc'))][-1]
  mcfile  = '{outdir}/{objid}/{mcstr}'.format(outdir=outdir, objid=objid, mcstr=mcstr)
  
  obj = mplo.Result(mcfile)
  thetas  = obj.quantilize(basic=True)
  model_params  = {}
  for key, val in paramsDict.items():
    median, upper, lower  = thetas[0, key], thetas[1, key], thetas[2, key]
    modelStr  = 'init_{par}'.format(par=val)
    modelStrDisp  = 'init_disp_{par}'.format(par=val)
    #model_params[modelStr]  =  median
    model_params[modelStr]  =  obj.theta_max()[key]
    model_params[modelStrDisp]  = np.max((upper-median, median-lower))
  return model_params

def main_boop(objid, nstart=0, nruns=8, outdir='cluster_run', run_params=None,
              model_params=None):
  """Calls prospector to run a parameter file.
  :param objid:
    Object of interest as specified in photometry data file.
  
  :param nstart: (optional, default: 0)
    Somewhat deprecated usage. Specifies at what run number to start the
    process. Numbering starts from 0.

  :param nruns: (optional, default: 8)
    Definitely deprecated usage. Specifies number of runs to do in prospector.
    Initially written as a way to customize burn-in process.

  :param outdir: (optional, default: 'cluster_run')
    The output directory to which to store the runs for each objid.

  :param run_params: (optional, default: None)
    Optional run parameters. See init_run_params for potential values to
    specify.

  :param model_params: (optional, default: None)
    Optional model parameters. See init_model_params for potential values to
    specify.
  """
  # ensure that the directory structure exists
  outfileTemp = '{outdir}/{objid}/{objid}_run{num}'
  try:
    os.mkdir('{outdir}/{objid}'.format(outdir=outdir,objid=objid))
  except OSError:
    try:
      os.mkdir('{outdir}'.format(outdir=outdir))
      os.mkdir('{outdir}/{objid}'.format(outdir=outdir,objid=objid))
    except OSError:
      pass

  # create default dictionaries
  if run_params == None:
    run_params = {}
  elif not isinstance(run_params, dict):   # can remove
    print 'Your run_params are invalid. Defaulting.'
    run_params = {}
  # IGNORE
  #run_params['nwalkers'] = 128 
  #run_params['niter'] = 1280
  if model_params == None:
    model_params = {}
  elif not isinstance(model_params, dict): # can remove
    print 'Your model_params are invalid. Defaulting.'
    model_params = {}
 
  run_params['objid'] = objid
  runStr  = init_run_params(**run_params)

  for run in range(nstart, nstart+nruns):
    outfile   = outfileTemp.format(outdir=outdir, objid=objid, num=run)
    paramFile = outfile+'_params.py'

    # recalculate new parameters if run is not zeroth run
    if run!=0:
      model_params  = get_new_params(objid, run-1, outdir=outdir)

    # create a parameter file
    model_boop(paramFile, **model_params)

    # update parameter inputs into my_prospector
    run_params['param_file']  = paramFile
    run_params['outfile']     = outfile

    # create call str for my_prospector
    runStr  = init_run_params(**run_params)
    print runStr
    subprocess.call(runStr, shell=True)

def clustered_main_boop(objs, **boop_args):
  """Perform model fits for many objects. Assumes every object uses the same
  arguments.

  :param objs:
    Iterable object (or not) specifying which objects to run as specified by
    the input photometry data file.

  :param boop_args: (optional)
    Arguments to send to main_boop. Possible arguments are:
      nstart
      nruns
      outdir
      run_params, would not recommend use here
      model_params, would not recommend use here
  """
  processes = []
  for obj in objs:
    processes.append(Process(target=main_boop, args=(obj,), kwargs=boop_args))

  for p in processes:
    p.start()
  for p in processes:
    p.join()
  
def final_run(objs, outdir='cluster_run'):
  """Deprecated. Please ignore. Finds the most recent model fits for a set of
  objects and supersamples the resultant posterior distribution.

  :param objs:
    The objects in question as specified by the photometry data file.
  :param outdir: (optional)
    The output directory where the runs are stored.
  """
  boop_args = {'nruns': 1, 'outdir': outdir}
  boop_args['run_params'] = {'niter': 1000, 'nwalkers': 128}  # overkill

  # for each object, get final most recent run, then make process
  processes = []
  for obj in objs:
    nstart = get_final_run(obj, outdir)
    if nstart is None:
      continue

    boop_args['nstart'] = nstart+1
    processes.append(Process(target=main_boop, args=(obj,), kwargs=boop_args))

  for p in processes:
    p.start()
  for p in processes:
    p.join()

def get_final_run(obj, outdir='final_runs', givemc=False):
  """Get the run number of the mcfile string value of the most recent model fit
  that was run for a
  specified object.

  :param obj:
    The object as specified by the photometry data file.

  :param outdir:
    The output directory where the runs for each object are stored.

  :param givemc: (optional, default: False)
    Flag indicating what user wants returns. If False, will return the most
  recent run number. If True, will return the most recent run's mcmc file.
  """
  # query the directory
  try:
    allfiles  = os.listdir('{outdir}/{obj}'.format(outdir=outdir, obj=obj))
  except OSError:
    return None

  # query all the parameter files
  params    = sorted([f for f in allfiles if f.endswith('_mcmc')])

  # get a run number
  nstart    = int( (params[-1].split('_')[1] ).strip('run') )

  # run comparison
  for param in params:
    x = int( ( param.split('_')[1] ).strip('run') )
    if x > nstart:
      nstart = x

  # decide which value to return
  if givemc == False:
    return nstart
  else:
    mcrun  = 'run{}'.format(nstart)
    mcfile = [f for f in params if mcrun in f][0]
    mcfile = '{outdir}/{obj}/{mc}'.format(outdir=outdir,obj=obj,mc=mcfile)
    return mcfile

def final_run_seds(fileIn='grblist.txt'):
  grbs = np.genfromtxt(fileIn, names=None, dtype=str)
  sed_args  = {'numSpec': 5000}
  processes = []
  for grb in grbs:
    mcfile  = get_final_run(grb, givemc=True)
    if mcfile is None:
      continue
    sed_args['savefile'] = 'cluster_plots/{grb}_SED.pdf'.format(grb=grb)
    obj = mplo.Result(mcfile)
    processes.append(Process(target=obj.plot_max_sed, kwargs=sed_args))

  for p in processes:
    p.start()
  for p in processes:
    p.join()

def final_run_corners(fileIn='grblist.txt'):
  grbs = np.genfromtxt(fileIn, names=None, dtype=str)
  corner_args  = {}
  processes = []
  for grb in grbs:
    mcfile  = get_final_run(grb, givemc=True)
    if mcfile is None:
      continue
    corner_args['savefile'] = 'cluster_plots/{grb}_corner.pdf'.format(grb=grb)
    obj = mplo.Result(mcfile)
    processes.append(Process(target=obj.plot_corner, kwargs=corner_args))

  for p in processes:
    p.start()
  for p in processes:
    p.join()
