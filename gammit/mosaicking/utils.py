import os

BASE_DIR = '/'.join(os.getcwd().split('/')[:-2])

MOSAIC_DIR = os.getenv(
    'MOSAIC_DIR',
    '{base}/data/mosaics'.format(base=BASE_DIR)
)

EXP_FRAM_THRESH = 0.86  # want exptime/framtime > thresh
