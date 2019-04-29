import argparse
import logging
import os
from itertools import compress
import re

from astropy.io import fits

import utils

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
    '-t', '--thresh',
    default=0.86,
    help=('Good exposure/framtime threshold')
)
args = vars(parser.parse_args())

GRB = args['grb']
CH = args['ch']
THRESH = args['thresh']

def get_tagged_fits_files(tag, grb, ch, grb_dir):
    tag_prefix = '_{tag}.fits'.format(tag=tag)
    ch_prefix = 'ch{ch}'.format(ch=ch)
    for path, subdirs, files in os.walk(grb_dir):
        if ch_prefix not in path:
            continue
        for name in files:
            if name.endswith(tag_prefix):
                yield os.path.join(path, name)

def above_thresh(fits_file, thresh):
    hdulist = fits.open(fits_file)
    exptime = hdulist[0].header['EXPTIME']
    framtime = hdulist[0].header['FRAMTIME']
    ratio = exptime/framtime
    return ratio > thresh

def write(tag, grb_dir, ch, fits_files):
    write_file = os.path.join(
        grb_dir, '{tag}_ch{ch}_list.txt'.format(
            tag=tag, ch=ch
        )
    )
    with open(write_file, 'w') as wf:
        wf.writelines('\n'.join(fits_files))
    logging.info('File {} written'.format(write_file))

def main(grb, ch, thresh=THRESH):
    grb_dir = os.path.join(utils.MOSAIC_DIR, grb)
    cbcd_files = get_tagged_fits_files('cbcd', grb, ch, grb_dir)
    cbcd_files = sorted(list(cbcd_files))
    cbcd_files.pop(0)   # remove calibration file
    good_cbcd_mask = [above_thresh(i, thresh) for i in cbcd_files]
    good_cbcd_files = list(compress(cbcd_files, good_cbcd_mask))

    good_cbunc_files = [re.sub('cbcd', 'cbunc', f) for f in good_cbcd_files]
    good_bimsk_files = [re.sub('cbcd', 'bimsk', f) for f in good_cbcd_files]

    write('cbcd', grb_dir, ch, good_cbcd_files)
    write('cbunc', grb_dir, ch, good_cbunc_files)
    write('bimsk', grb_dir, ch, good_bimsk_files)

if __name__ == '__main__':
    main(GRB, CH)
