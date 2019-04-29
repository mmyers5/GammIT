###### Define which modules to run ######

compute_uncertainties_internally = 0
have_uncertainties = 1
run_fiducial_image_frame = 1
run_mosaic_geom = 0
run_medfilter = 1
run_detect_radhit = 1
run_mosaic_interp = 1
run_detect_outlier = 1
run_mosaic_proj = 1
run_mosaic_covg = 1
run_mosaic_dual_outlier = 1
run_level = 1
run_mosaic_outlier = 1
run_mosaic_box_outlier = 1
run_mosaic_rmask = 1
run_mosaic_reinterp = 0
run_fix_coverage = 0
run_mosaic_coadder = 1
run_mosaic_combiner = 1
run_mosaic_medfilter = 0
create_rmask_mosaic = 0
make_array_corr_files = 0
make_array_corr_mosaic = 0


##### Define output of run #####

run_median_mosaic = 0
run_absolute_minimum_mosaic = 0
create_std_mosaic = 1
create_unc_mosaic = 1
create_outlier_mosaic = 0
create_dual_outlier_mosaic = 0

###### Test for array location corr ######

delete_intermediate_files_array_corr = 0

###### Input file lists ######

IMAGE_STACK_FILE_NAME = /Users/mmyers/GammIT/data/mosaics/120119A/cbcd_ch2_list.txt
SIGMALIST_FILE_NAME = /Users/mmyers/GammIT/data/mosaics/120119A/cbunc_ch2_list.txt
DCE_STATUS_MASK_LIST = /Users/mmyers/GammIT/data/mosaics/120119A/bimsk_ch2_list.txt
PMASK_FILE_NAME = /Users/mmyers/GammIT/tools/mopex/cal/super_masks/chan2_ormask_bcd.fits
RMASK_LIST = 

FIF_FILE_NAME =

###### Output Dir and Files ######

OUTPUT_DIR = /Users/mmyers/GammIT/data/mosaics/120119A/results_ch2
MEDFILTER_DIR = Medfilter-mosaic
INTERP_DIR = Interp-mosaic
DETECT_DIR = Detect-mosaic
DUAL_OUTLIER_DIR = DualOutlier-mosaic
OUTLIER_DIR = Outlier-mosaic
BOX_OUTLIER_DIR = BoxOutlier-mosaic
RMASK_DIR = Rmask-mosaic
REINTERP_DIR = Reinterp-mosaic
COADDER_DIR = Coadd-mosaic
COMBINER_DIR = Combine-mosaic

###### Mask Bit Parameter ######

DCE_Status_Mask_Fatal_BitPattern = 32520
DCE_Status_Mask_Radhit_Bit = 9
PMask_Fatal_BitPattern = 32767
RMask_Fatal_BitPattern = 15

###### Global Parameters ######

USE_REFINED_POINTING = 0
USE_OUTLIER_FOR_RMASK = 1
USE_BOX_OUTLIER_FOR_RMASK = 1
USE_DUAL_OUTLIER_FOR_RMASK = 1
MOSAIC_PIXEL_SIZE_X = -1.1111E-4
MOSAIC_PIXEL_SIZE_Y = 1.1111E-4

###### Other Parameters ######

overwrite_dmask = 0
keep_coadded_tiles = 1
sigma_weighted_coadd = 0
delete_intermediate_files = 0
do_multiprocess = 'on'
DMASK_DIR = Dmask-mosaic
RMASK_MOSAIC_DIR = RmaskMosaic-mosaic
ARRAY_CORR_IMAGE = 
SIGMA_DIR = Sigma-mosaic
ncpu_multiprocess = 1
ARRAY_PIXAREA_IMAGE = 


###### Modules ######

&SNESTIMATORIN
Gain = 1290.0,
Read_Noise = 90.0,
Confusion_Sigma = 0.0,
&END

&FIDUCIALIMAGEFRAMEIN
Edge_Padding = 10,
Projection_Type = 'TAN',
Coordinate_System = 'J2000',
CROTA2 = A,
&END

&MOSAICGEOM
&END

&MEDFILTER
Window_X = 45,
Window_Y = 45,
N_Outliers_Per_Window = 50,
Use_Sbkg_for_Med = 0,
Sbkg_Filt_Size_for_Med = 1,
&END

&DETECT_RADHIT
Segmentation_Threshold = 3.0,
Detection_Max_Area = 3,
Radhit_Threshold = 6.0,
&END

&MOSAICINTIN
INTERP_METHOD = 2,
FINERES = 0.0,
DRIZ_FAC = 0.8,
GRID_RATIO = 2,
ALPHA = -0.5,
&END

&DETECT
Detection_Max_Area = 100,
Detection_Min_Area = 0,
Detection_Threshold = 4.0,
Threshold_Type = 'simple',
Detect_Min_Peak_Fraction = 0.0,
&END

&MOSAICPROJIN
&END

&MOSAICCOVGIN
TILEMAX_X = 2000,
TILEMAX_Y = 2000,
&END

&MOSAICDUALOUTLIERIN
MAX_OUTL_IMAGE = 2,
MAX_OUTL_FRAC = 0.51,
TILE_XSIZ = 1000,
TILE_YSIZ = 1000,
&END

&LEVEL
Threshold_Ratio = 0.5,
&END

&MOSAICOUTLIERIN
THRESH_OPTION = 1,
BOTTOM_THRESHOLD = 0.0,
TOP_THRESHOLD = 0.0,
MIN_PIX_NUM = 3,
TILE_XSIZ = 1000,
TILE_YSIZ = 1000,
&END

&MOSAICBOXOUTLIERIN
BOX_X = 3,
BOX_Y = 3,
BOX_MEDIAN_BIAS = 1,
TILE_XSIZ = 500,
TILE_YSIZ = 500,
&END

&MOSAICRMASKIN
RM_THRESH = 0.8,
BOTTOM_THRESHOLD = 6.0,
TOP_THRESHOLD = 6.0,
MIN_COVERAGE = 4,
MAX_COVERAGE = 100,
REFINE_OUTLIER = 1,
REFINE_OUTLIER_THRESH = 12,
BOX_BOTTOM_THRESHOLD = 6.0,
BOX_TOP_THRESHOLD = 6.0,
BOX_MIN_COVERAGE = 2.0,
&END

&MOSAICREINTIN
&END

&FIX_COVERAGE
Min_Single_Coverage = 0.95,
Min_Block_Coverage = 0.83,
&END

&MOSAICCOADDIN
TILEMAX_X = 1000,
TILEMAX_Y = 1000,
INTEG_TIME_KWD = EXPTIME,
&END

&MOSAICCOMBINER
&END

&MOSAIC_MEDFILTER
Window_X = 45,
Window_Y = 45,
N_Outliers_Per_Window = 500,
Use_Sbkg_for_Med = 0,
Sbkg_Filt_Size_for_Med = 3,
&END

&CREATERMASKMOSAIC
&END

&MAKEARRAYCORRFILES
&END

&MAKEARRAYCORRMOSAIC
&END

#END


