import os

BASE_DIR = '/'.join(os.getcwd().split('/'))

MOSAIC_DIR = os.getenv(
    'MOSAIC_DIR',
    '{base}/data/mosaics'.format(base=BASE_DIR)
)
