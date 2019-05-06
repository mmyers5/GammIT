import os
import re
import sys

DESIRED_GRBS = ['040924', '050908', '041006', '061007', '050824', '070318']
CLUSTER_DIR = '/scratch2/GRBH-MZSFR/Analysis/prospector/final_runs_old'
NITER = 500

def most_recent_h5(grb, cluster_dir):
  grb_dir = '{}/{}'.format(cluster_dir, grb)
  all_h5 = []
  re_pattern = '^\w+(run\d+)'
  for item in os.listdir(grb_dir):
    if item.endswith('.h5'):
      all_h5.append(item)

  max_iter = 0
  max_iter_string = None
  for item in all_h5:
    item_run_iter = re.match(re_pattern, item)
    item_iter = re.findall('\d+', item_run_iter.group(1))
    if item_iter > max_iter:
      max_iter = item_iter
      max_iter_string = item_run_iter.string
  return '{}/{}/{}'.format(cluster_dir, grb, max_iter_string)

for grb in DESIRED_GRBS:
  h5_file = most_recent_h5(grb, CLUSTER_DIR)
  sys.argv = [
    '--niter={}'.format(NITER),
    '--restart_from={}'.format(h5_file)]
  execfile('prospector_restart.py')
