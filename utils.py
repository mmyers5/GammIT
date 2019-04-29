import os

BASE_DIR = os.getcwd()

MOSAIC_DIR = os.getenv(
    'MOSAIC_DIR',
    '{base}/data/mosaics'.format(base=BASE_DIR)
)

MOPEX_DIR = os.getenv(
    'MOPEX_DIR',
    os.path.join(BASE_DIR, 'tools/mopex')
)
