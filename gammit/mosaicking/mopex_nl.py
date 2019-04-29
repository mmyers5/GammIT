import os
import argparse

import utils
from gammit.mosaicking import mopex_listmaker

NL_DIR =  os.path.join(utils.MOSAIC_DIR, 'cdf')
DEFAULT_NL = os.path.join(NL_DIR, 'default.nl')

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
args = vars(parser.parse_args())

GRB = args['grb']
CH = args['ch']

def create_nl(grb, ch):
    with open(DEFAULT_NL, 'r') as f:
        nl = f.read()
    grb_dir = os.path.join(utils.MOSAIC_DIR, grb)
    cbcd_list = 'cbcd_ch{ch}_list.txt'.format(ch=ch)
    if not os.path.isfile(cbcd_list):
        mopex_listmaker.main(grb, ch)

    image_stack_file_name = os.path.join(
        grb_dir, 'cbcd_ch{ch}_list.txt'.format(ch=ch)
    )
    sigmalist_file_name = os.path.join(
        grb_dir, 'cbunc_ch{ch}_list.txt'.format(ch=ch)
    )
    dce_status_mask_list = os.path.join(
        grb_dir, 'bimsk_ch{ch}_list.txt'.format(ch=ch)
    )
    pmask_file_name = os.path.join(
        utils.MOPEX_DIR,
        'cal/super_masks/chan{ch}_ormask_bcd.fits'.format(ch=ch)
    )
    output_dir = os.path.join(
        grb_dir, 'results_ch{ch}'.format(ch=ch)
    )
    nl = nl.format(
        image_stack_file_name=image_stack_file_name,
        sigmalist_file_name=sigmalist_file_name,
        dce_status_mask_list=dce_status_mask_list,
        pmask_file_name=pmask_file_name,
        output_dir=output_dir
    )

    nl_file_name = os.path.join(
        NL_DIR, '{grb}_ch{ch}.nl'.format(grb=grb, ch=ch)
    )
    with open(nl_file_name, 'w') as f:
        f.write(nl)

    return nl_file_name

if __name__ == '__main__':
    create_nl(GRB, CH)
