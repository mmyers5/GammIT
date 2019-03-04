"""Delete filenames that have insufficient exposure times and create
`cbcdlist.txt`, `cbunclist.txt`, and `bimsklist.txt`.
"""
import os
import argparse
import logging
from itertools import compress
from astropy.io import fits as pyfits

DEFAULT_EXP_TIME = 0.86    # desire exptime/frametime > GOOD_EXP_TIME

# parse args
parser = argparse.ArgumentParser()
parser.add_argument(
    '-g', '--grb',
    default=None,
    help=('Desired GRB')
)
parser.add_argument(
    '-c', '--ch',
    default=1,
    help=('Desired channel')
)
parser.add_argument(
    '-e', '--exp',
    default=DEFAULT_EXP_TIME,
    help=('Good exposure/framtime ratio')
)
args = vars(parser.parse_args())
GOOD_EXP_TIME = args['exp']

DATA_DIR = os.path.join(os.getcwd(), args['grb'])

TAGS = ['cbcd', 'cbunc', 'bimsk']

def get_tagged_fits_files(tag):
    for path, subdirs, files in os.walk(DATA_DIR):
        for name in files:
            if name.endswith('.fits') and 'ch'+args['ch'] in name:
                yield os.path.join(path, name)

def good_exptime(fits_file):
    hdulist = pyfits.open(fits_file)
    exptime = hdulist[0].header['EXPTIME']
    framtime = hdulist[0].header['FRAMTIME']
    ratio = exptime/framtime
    logging.info('Ratio {} for {}'.format(ratio, fits_file))
    return ratio >= GOOD_EXP_TIME

def write(tag, fits_files):
    write_file = os.path.join(
        DATA_DIR, '{}_ch{}_list.txt'.format(
            tag, args['ch']
        )
    )
    with open(write_file, 'w') as wf:
        wf.writelines('\n'.join(fits_files))
    logging.info('File {} written'.format(write_file))

if __name__ == '__main__':
    for tag in TAGS:
        fits_files = list(get_tagged_files(tag))[1:] # remove calibration file
        good_fits_mask = [good_exptime(i) for i in fits_files]
        good_fits_files = list(compress(fits_files, good_fits_mask))
        write(tag, good_fits_files)
