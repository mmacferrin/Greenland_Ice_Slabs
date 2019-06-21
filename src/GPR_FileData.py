# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:41:45 2016

@author: mmacferrin
email: michael.macferrin@colorado.edu

FileData.py -- Handles file locations for all imported and exported files.
Specifically for MacFerrin, et. al. (2019), "Rapid expansion of Greenland's low-permeability ice slabs"

Used by: MacFerrin_IceLenses_MAIN.py
         FirnCore_Manager.py
         InSituGPR_Manager.py
"""

import os

## DATA locations
# FIGSHARE REPOSITORY DATA AT: https://figshare.com/account/home#/projects/47690
# You can download the data there, put the local folder name in the line below, and go from there.
# NOTE: Not all the intermediate data is cached in the repository. Much of it must be generated.
# See contact info above for details.
FIGSHARE_BASE = "/path/to/figshare/folder"
DATA_FIGURES_OUTPUT_FOLDER = os.path.join(FIGSHARE_BASE, "Figures")
GPR_DIR = r'C:\Users\mmacferrin\Dropbox\Research\DATA\ACT-13 GPR Data\MALAGS'

# Pre-prepared picklefiles for the ACT-IceBridge validation, so we don't need to redo them each time.
GPR_DOWNSAMPLED_FOLDER = os.path.join(GPR_DIR, "DOWNSAMPLE_PICKLEFILES")
ACT13_DOWNSAMPLED_PICKLEFILE = os.path.join(GPR_DOWNSAMPLED_FOLDER, "ACT13_DOWNSAMPLED_TO_20130409_01_010_012.pickle")
ACT13_DOWNSAMPLED_COUNT_PICKLEFILE = os.path.join(GPR_DOWNSAMPLED_FOLDER, "ACT13_DOWNSAMPLED_CELL_COUNT.pickle")
ACT13_ICEBRIDGE_SUBSET__MASK_HMASK_VMASK_PICKLEFILE = os.path.join(GPR_DOWNSAMPLED_FOLDER, "20130409_01_010_012_SUBSET_MASK_HMASK_VMASK.pickle")
ACT13_ICEBRDIGE__LAT_LON_E_N_DEPTHS_PICKLEFILE =  os.path.join(GPR_DOWNSAMPLED_FOLDER, "ICEBRIDGE_SUBSET_LAT_LON_N_E_DEPTHS.pickle")
ACT13_SUBSET__LAT_LON_ELEV_E_N_DEPTHS_PICKLEFILE = os.path.join(GPR_DOWNSAMPLED_FOLDER, "GPR_SUBSET_LAT_LON_ELEV_N_E_DEPTHS.pickle")

# ALL of the GPR files, including ones that are on transects with ones that aren't.
GPR_FILEBASES = ["DAT_0118_A1",
                 "DAT_0119_A1",
                 "DAT_0121_A1",
                 "DAT_0122_A1",
                 "DAT_0126_A1",
                 "DAT_0127_A1",
                 "DAT_0128_A1",
                 "DAT_0129_A1",
                 "DAT_0131_A1",
                 "DAT_0132_A1",
                 "DAT_0134_A1",
                 "DAT_0137_A1",
                 "DAT_0138_A1",
                 "DAT_0143_A1",
                 "DAT_0144_A1",
                 "DAT_0149_A1",
                 "DAT_0150_A1",
                 "DAT_0151_A1",
                 "DAT_0152_A1",
                 "DAT_0154_A1"
                ]

# (Index into 20-element array,
#  Filebase,
#  Switch?:  True: reverse the traces, False: Keep them in sequential order
#  0: Don't crop either, 1: Use 1st, crop 2nd, 2: Use 2nd, crop 1st.)
MAIN_TRANSECT_COLLECTION = [(3,  "DAT_0122_A1", True,  0),
                            (0,  "DAT_0118_A1", False, 1),
                            (8,  "DAT_0131_A1", False, 1),
                            (10, "DAT_0134_A1", False, 1),
                            (12, "DAT_0138_A1", False, 0),
                            (17, "DAT_0151_A1", False, 0),
                            (14, "DAT_0144_A1", True,  1)]


KANU_TRANSECT_COLLECTION = [(5, "DAT_0127_A1", False, 0),
                            (6, "DAT_0128_A1", False, 1),
                            (7, "DAT_0129_A1", False, 2)]

NON_TRANSECT_GPR_LINES = ["DAT_0119_A1",
                          "DAT_0121_A1",
                          "DAT_0126_A1",
                          "DAT_0132_A1",
                          "DAT_0137_A1",
                          "DAT_0143_A1",
                          "DAT_0152_A1",
                          "DAT_0154_A1"]

ASCII_DIR = os.path.join(GPR_DIR, r'ASCII')
IMG_DIR = os.path.join(GPR_DIR, 'TraceImages')
RAW_IMG_DIR = os.path.join(IMG_DIR, "Raw")
RESAMPLED_DIR = os.path.join(IMG_DIR, "resampled")
RESAMPLED_COR_DIR = os.path.join(ASCII_DIR, "resampled")
EXPORT_DIR = os.path.join(RESAMPLED_DIR, "exported")

TRANSECT_GPR_MERGED = os.path.join(IMG_DIR, 'ACT_Transect.tif')
TRANSECT_GPR_MERGED_RESAMPLED = os.path.join(RESAMPLED_DIR, 'ACT_Transect_resampled.tif')
KANU_GPR_MERGED_RESAMPLED = os.path.join(RESAMPLED_DIR, 'KANU_resampled.tif')
TRANSECT_COR_MERGED = os.path.join(GPR_DIR, r'ASCII\ACT_Transect.cor')
TRANSECT_COR_MERGED_RESAMPLED = os.path.join(RESAMPLED_COR_DIR, r'ACT_Transect_resampled.cor')
KANU_COR_MERGED_RESAMPLED = os.path.join(RESAMPLED_COR_DIR, r'KANU_resampled.cor')

def corfilename_to_picklefilename(corfile):
    '''Input a .cor file (resampled), output the appropriate resampled_logvariance_picklefile.'''
    if corfile.find("KANU_resampled") > -1:
        f = os.path.splitext(KANU_GPR_MERGED_RESAMPLED)[0] + "_logvariance.pickle"
    elif corfile.find("ACT_Transect_resampled") > -1:
        f = os.path.splitext(TRANSECT_GPR_MERGED_RESAMPLED)[0] + "_logvariance.pickle"
    else:
        assert corfile.find("DAT_") > -1
        f = os.path.splitext(os.path.join(RESAMPLED_DIR, os.path.split(corfile)[-1]))[0] + "_resampled_logvariance.pickle"

    return f

def logvariance_to_detrended_filename(picklefile_name):
    '''Take the name of a "_resampled_logvariance.pickle" file and replace it with a
    "_resampled_logvariance_detrended.pickle" name, returning the new name.'''
    base, ext = os.path.splitext(picklefile_name)
    return base + "_detrended" + ext

# Coordinate files that include transect plus non-transect lines, but NOT individual cores lines that were combined to make a transect.
# lines 149 & 150 are lines that were re-driven at the time, and hence are duplicate data.  Leaving those out.
RESAMPLED_COORDINATE_FILES = [TRANSECT_COR_MERGED_RESAMPLED, KANU_COR_MERGED_RESAMPLED] + [os.path.join(RESAMPLED_COR_DIR, f + '.cor') for f in NON_TRANSECT_GPR_LINES]
RESAMPLED_GPR_LOGVARIANCE_PICKLEFILES = [corfilename_to_picklefilename(f) for f in RESAMPLED_COORDINATE_FILES]
GPR_DETRENDED_PICKLEFILES = [logvariance_to_detrended_filename(f) for f in RESAMPLED_GPR_LOGVARIANCE_PICKLEFILES]

# Firn Cores
CORE_XLSX_FILE = os.path.join(FIGSHARE_BASE, r'firn_cores_2009_2012_2013_2015_2016_temp.xlsx')

# IceBridge Data Folder
ICEBRIDGE_DATA_FOLDER = r'F:\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\accum'
ICEBRIDGE_DATA_H5FILE = os.path.join(ICEBRIDGE_DATA_FOLDER, "IceBridgeDB.h5")
ICEBRIDGE_EXPORT_FOLDER = os.path.join(ICEBRIDGE_DATA_FOLDER, "exported")
# TODO: Fix Folder Names
ICEBRIDGE_ICELENS_QUICKLOOK_FOLDER = r'C:\Users\mmacferrin\Dropbox\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\Ice Lenses Seen'

# TODO: Fix Folder Names
ICEBRIDGE_EXCLUSIONS_FOLDER = r'C:\Users\mmacferrin\Dropbox\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\Exclusion Guide Files'
ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, r'Exclusions_SURFACE_PICKS.txt')
ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, r'SURFACE_STARTING_PICKS_Suggestions.txt')
ICEBRIDGE_SURFACE_INDICES_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Surface_Indices_Picklefiles")
ICEBRIDGE_SURFACE_SLICE_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Surface_Slice_100m_Picklefiles")
ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "SURFACE_MISMATCH_EXCLUSIONS.txt")
ICEBRIDGE_EXCLUSIONS_LAKES_OTHER_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "LAKES_AND_OTHER_EXCLUSIONS.txt")
ICEBRIDGE_SAMPLE_DEPTHS_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Sample Depths_Picklefiles")
ICEBRIDGE_ROLL_CORRECTED_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Roll_Corrected_Picklefiles")
ICEBRIDGE_ROLL_CORRECTION_OUTPUT_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "ROLL_CORRECTION_PARAMETERS.csv")
ICEBRIDGE_DEPTH_CORRECTED_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Depth_Corrected_Picklefiles")
ICEBRIDGE_DEPTH_CORRECTION_OUTPUT_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "DEPTH_CORRECTION_PARAMETERS.csv")
ICEBRIDGE_VALIDATION_OUTPUT_CSV_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "Validation_Outputs.csv")
ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Boolean Array Picklefiles")
ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE = os.path.join(ICEBRIDGE_EXCLUSIONS_FOLDER, "Ice_Layer_Output_Thicknesses.csv")
ICEBRIDGE_SMOOTHED_ICE_LAYER_SHAPEFILE = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "Exported Shapefiles\IceBridge_out.shp")


