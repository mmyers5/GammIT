import numpy as np
import fast_phot as fap
import fits_plot as fip

def plot_em_all(infile='master_file.txt', indices=None, tcol='k'):
  master = np.genfromtxt(infile, names=True, dtype=None)
  if indices == None:
    index_mask = np.unique(master['GRB'], return_index=True)[1]
    indices = tuple(master['index'][index_mask])
  
  for i in indices:
    entry = master[master['index']==i]
    GRB, ch, z = entry['GRB'][0], entry['ch'][0], entry['z'][0],
    r, circol  = entry['r'][0], entry['circol'][0]
    ra, dec = entry['ra'][0], entry['dec'][0]

    infile = '{GRB}/{GRB}_ch{ch}_maic.wcs.fits'.format(GRB=GRB, ch=ch)
    savefile = 'plots/stamps/{GRB}_ch{ch}_stamp.eps'.format(GRB=GRB, ch=ch)  

    X, Y = fap.get_phys(infile, ra, dec, phys=True)

    fip.make_stamp(GRB, ch, X, Y, z, r, circol, tcol=tcol, savefile=True)

def plot_models(infile='master_file.txt', indices=None):
