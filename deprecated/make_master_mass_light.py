import numpy as np
from astropy.cosmology import WMAP9 as cosmo

master = np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/master_file_out_phot_usage.txt',\
                  names=True,dtype=None,delimiter='\t')
ref    = np.genfromtxt('/scratch2/GRBH-MZSFR/Analysis/master_file.txt',\
                  names=True,dtype=None,delimiter='\t')
master_out = open('/scratch2/GRBH-MZSFR/Analysis/master_file_out_phot_usage_ageDL.txt','w')
master_out.write('GRB\tz\tflx\tflx_unc\tage_Gyr\tDL_Mpc\n')

for idx in np.arange(len(master)):
  if master['usage'][idx]=='*' and master['ch'][idx]==1:
    grb = master['GRB'][idx]
    z   = [entry['z'] for entry in ref if entry['GRB']==grb][0]
    flx = master['flux'][idx]
    flx_unc = master['flux_unc'][idx]

    age_Gyr = cosmo.age(z).value
    DL  = cosmo.luminosity_distance(z).value
  
    master_out.write('{GRB}\t{z}\t{flx}\t{flxunc}\t{age}\t{DL}\n'.\
                   format(GRB=grb,z=z,flx=flx,flxunc=flx_unc,age=age_Gyr,DL=DL))
master_out.close()
