# Code Audit

## calc_SFH

Utility script to calculate star formation histories.

## calc_all_phot

Execution script to calculate photometry from a file with pixel counts.

### Input

A file with the following fields:

- GRB
- corrCorr: correction on the aperture correction
- flxSubCnts: pixel count of background-subtracted image
- uncCnts: pixel count of source uncertainty
- apPix: number of pixels in source aperture
- uncBkg: pixel count of background uncertainty
- anPix: pixel count of background annulus
- ch: Spitzer IRAC channel number
- ftype: flux type ("unc", "flx", or "sub")
- ap: aperture size string
- sig: some sort of sigma value

### Output

A tab-delimited file with the following fields:

- GRB
- ch
- type
- flux
- flux_unc
- mab
- mab_unc
- ap
- sig

## copytie

Execution script to copy headers between FITS files, usually used to copy the header from mosaic that has been astrometrically tied to an optical afterglow image and apply the astrometric tie to a mosaic in different channel.

## distribution_maker

Utility script to plot sample distribution by AB magnitude and and redshift.

## fast_phot

Utility script to calculate photometry. Used in `calc_all_phot`.

## fits_plot

Utility script to create stamps from FITS images.

## fits_plot_old

I would assume it's a copy of fits_plot...

## galfit_config

GalfitConfig class which appears to simply verify galfit objects.

## galfit_field

Classes GalfitSource and GalfitField for area modeling.

## galfit_init_original

Old version for `galfit_config`.

## galfit_master

Execution script to model fields.

### Inputs

GRB and ch as well as a file that doesn't appear to do anything. Needs more digging.

## galfit_modifier

Dummy file.

## galfit_writer

Dummy file.

## imgmatch_manual
