# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 13:28:15 2017

@author: mmacferrin
"""

#####################
## USER NOTE: Most of the .nc user files are *not* located in the FigShare repository.
# The files (and/or newer versions of them) are available from the groups that
# develop and process each of these Regional Climate Models.
# The file names are included here so users can see what files I had and what we did
# with them in RCM_Manager.py
#####################

import os
# Get one variables from the Ice Lens Observational Data
from GPR_FileData import ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, FIGSHARE_BASE

DATA_BASEDIR = os.path.join(FIGSHARE_BASE, "RCMs")

ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE = ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE
STATISTICS_OUTPUT_FOLDER = os.path.join(FIGSHARE_BASE, "IceBridge_Accumulation_Radar", "txt")
ICEBRIDGE_THICKNESS_VALIDATION_TRACKS_LIST_FILE = os.path.join(STATISTICS_OUTPUT_FOLDER, 'Modeled_Ice_Validation_Tracks_List.txt')
EMI_VALIDATION_CUTOFF_VALUES_PICKLEFILE = os.path.join(STATISTICS_OUTPUT_FOLDER, r'_EMI_Validation_Cutoff_Results.pickle')
EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE = os.path.join(STATISTICS_OUTPUT_FOLDER, r'ICE_SLAB_Decadal_Outputs.csv')

#############################################################################################
## HIRHAM5
#############################################################################################

HIRHAM_DATA_FOLDER = os.path.join(DATA_BASEDIR,'HIRHAM')

# Mask and grid data
HIRHAM_MASK_DATAFILE = os.path.join(HIRHAM_DATA_FOLDER, "GL2offline_CountryGlacierMasks.nc")
HIRHAM_GRID_DATAFILE = os.path.join(HIRHAM_DATA_FOLDER, "Grid_Cellareas_GL2domain.nc")
HIRHAM_BASIN_DATAFILE = os.path.join(HIRHAM_DATA_FOLDER, "ZwallyMasks_all.nc")

# ERA Interim - "Current climate" -- best to evaluate against observations
HIRHAM_FILENAMES_2014_ERA = {"DAILY_MEANS"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GL2_input_1980_2014_YM.nc"),
                      "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_GL2LIN_Darcy_60m_liqCL_wh1_SMB_1980_2014_YM.nc"),
                      "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_GL2LIN_Darcy_60m_liqCL_wh1_MRROS_1980_2014_YM.nc"),
                      "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_GL2LIN_Darcy_60m_liqCL_wh1_SNMELT_1980_2014_YM.nc"),
                      "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GL2_input_1980_2014_YM.nc")}

# ECEearth GCM "Current Climate" (q2), best to run against future GCM outputs (biases consistent)
HIRHAM_FILENAMES_2010_GCM = {"YEARLY_MEANS" :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6q2_input_1991_2010_yearsum.nc"),
                      "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6q2_SMB_1991_2010.nc"),
                      "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6q2_MRROS_1991_2010.nc"),
                      "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6q2_SNMELT_1991_2010.nc"),
                      "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6q2_input_1991_2010_tas-yearmean.nc")}

# ECEearth GCM mid-Century 2031-2050 (r2), RCP 4.5 scenario
HIRHAM_FILENAMES_2050_GCM_RCP45 = {"YEARLY_MEANS" :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6r2_input_2031_2050_yearsum.nc"),
                            "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6r2_SMB_2031_2050.nc"),
                            "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6r2_MRROS_2031_2050.nc"),
                            "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6r2_SNMELT_2031_2050.nc"),
                            "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6r2_input_2031_2050_tas-yearmean.nc")}

# ECEearth GCM mid-Century 2031-2050 (y2), RCP 8.5 scenario
HIRHAM_FILENAMES_2050_GCM_RCP85 = {"YEARLY_MEANS" :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6y2_input_2031_2050_yearsum.nc"),
                            "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6y2_SMB_2031_2050.nc"),
                            "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6y2_MRROS_2031_2050.nc"),
                            "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6y2_SNMELT_2031_2050.nc"),
                            "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6y2_input_2031_2050_tas-yearmean.nc")}

# ECEearth GCM end-Century 2081-2100 (s2), RCP 4.5 scenario
HIRHAM_FILENAMES_2100_GCM_RCP45 = {"YEARLY_MEANS" :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6s2_input_2081_2100_yearsum.nc"),
                            "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6s2_SMB_2081_2100.nc"),
                            "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6s2_MRROS_2081_2100.nc"),
                            "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6s2_SNMELT_2081_2100.nc"),
                            "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6s2_input_2081_2100_tas-yearmean.nc")}

# ECEearth GCM end-Century 2081-2100 (z2), RCP 8.5 scenario
HIRHAM_FILENAMES_2100_GCM_RCP85 = {"YEARLY_MEANS" :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6z2_input_2081_2100_yearsum.nc"),
                            "YEARLY_SMB"   :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6z2_SMB_2081_2100.nc"),
                            "YEARLY_RUNOFF":os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6z2_MRROS_2081_2100.nc"),
                            "YEARLY_MELT"  :os.path.join(HIRHAM_DATA_FOLDER,"HH5_ECEarth_GR6z2_SNMELT_2081_2100.nc"),
                            "YEARLY_TEMP"  :os.path.join(HIRHAM_DATA_FOLDER,"HH_GR6z2_input_2081_2100_tas-yearmean.nc")}


HIRHAM_DATA_FILENAMES = [os.path.join(HIRHAM_DATA_FOLDER, HIRHAM_MASK_DATAFILE),
                  os.path.join(HIRHAM_DATA_FOLDER, HIRHAM_GRID_DATAFILE)] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2014_ERA.values()] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2010_GCM.values()] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2050_GCM_RCP45.values()] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2050_GCM_RCP85.values()] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2100_GCM_RCP45.values()] + \
                  [os.path.join(HIRHAM_DATA_FOLDER, fname) for fname in HIRHAM_FILENAMES_2100_GCM_RCP85.values()]

HIRHAM_ELEV_PICKLEFILE = os.path.join(HIRHAM_DATA_FOLDER, "HIRHAM_ELEVATIONS.pickle")
HIRHAM_POLYGON_WKTS_TEXTFILE_ALL = os.path.join(HIRHAM_DATA_FOLDER, r"WKT Polygons\HIRHAM_Polygon_WKTs_ALL.txt")
HIRHAM_POLYGON_WKTS_TEXTFILE_GREENLAND = os.path.join(HIRHAM_DATA_FOLDER, r"WKT Polygons\HIRHAM_Polygon_WKTs_GREENLAND.txt")
HIRHAM_ICEBRIDGE_CELLS_SUMMARY_CSV = os.path.join(HIRHAM_DATA_FOLDER, r'_HIRHAM_IceBridge_Cells_Summary.csv')
HIRHAM_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE = os.path.join(HIRHAM_DATA_FOLDER, r'_HIRHAM_IceBridge_Traces_Rows_Cols.pickle')
HIRHAM_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE = os.path.join(HIRHAM_DATA_FOLDER, r'_HIRHAM_Upper_Lower_Validation_indices.pickle')

#############################################################################################
## RACMO 2.3
#############################################################################################
RACMO_basedir = os.path.join(DATA_BASEDIR, "RACMO")
RACMO_EXPORT_SHAPEFILE_FOLDER = RACMO_basedir
RACMO_ICEBRIDGE_SUMMARY_CSV_FILE = os.path.join(RACMO_basedir, r"_RACMO_IceBridge_Cells_Summary.csv")
RACMO_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE = os.path.join(RACMO_basedir, r'_RACMO_IceBridge_Traces_Rows_Cols.pickle')
RACMO_GRID_INTEGERS_PICKLEFILE = os.path.join(RACMO_basedir, "RACMO Grid Integers.pickle")
RACMO_POLYGON_WKTS_TEXTFILE_GREENLAND = os.path.join(RACMO_basedir, r"RACMO_Polygon_WKTS_GrIS.txt")
RACMO_POLYGON_WKTS_TEXTFILE_ALL = os.path.join(RACMO_basedir, r"RACMO_Polygon_WKTS_ALL.txt")
RACMO_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE = os.path.join(RACMO_basedir, r'_RACMO_Upper_Lower_Validation_indices.pickle')

RACMO_FILENAMES_DICT = {}
RACMO_FILENAMES_DICT["RACMO2.3_ERA-Int_1958-2015"] = {"GRID"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_gridinfo.nc"),
                        "RAIN"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_precip_monthly_1958-2015.nc"),
                        "RUNOFF" : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_runoff_monthly_1958-2015.nc"),
                        "SMB"    : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_smb_monthly_1958-2015.nc"),
                        "SNOW"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_snowfall_monthly_1958-2015.nc"),
                        "MELT"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_snowmelt_monthly_1958-2015.nc"),
                        "SUBL"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_subl_monthly_1958-2015.nc"),
                        "TEMP_2M": os.path.join(RACMO_basedir, "RACMO2.3_GRN11_t2m_monthly_1958-2015.nc")}

RACMO_FILENAMES_DICT["RACMO2.1_HadGEM2_1971-2004"] = {"GRID"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_gridinfo.nc"),
                        "RAIN"   : os.path.join(RACMO_basedir, "mon_cum_precip_hist.nc"),
                        "RUNOFF" : os.path.join(RACMO_basedir, "mon_cum_runoff_hist.nc"),
                        "SMB"    : os.path.join(RACMO_basedir, "mon_cum_smb_hist.nc"),
                        "SNOW"   : os.path.join(RACMO_basedir, "mon_cum_snowfall_hist.nc"),
                        "MELT"   : os.path.join(RACMO_basedir, "mon_cum_snowmelt_hist.nc"),
                        "SUBL"   : os.path.join(RACMO_basedir, "mon_cum_subl_hist.nc"),
                        "TEMP_2M": os.path.join(RACMO_basedir, "mon_ave_t2m_hist.nc")}

RACMO_FILENAMES_DICT["RACMO2.1_HadGEM2_RCP4.5_2005-2098"] = {"GRID"   : os.path.join(RACMO_basedir, "RACMO2.3_GRN11_gridinfo.nc"),
                        "RAIN"   : os.path.join(RACMO_basedir, "mon_cum_precip.nc"),
                        "RUNOFF" : os.path.join(RACMO_basedir, "mon_cum_runoff.nc"),
                        "SMB"    : os.path.join(RACMO_basedir, "mon_cum_smb.nc"),
                        "SNOW"   : os.path.join(RACMO_basedir, "mon_cum_snowfall.nc"),
                        "MELT"   : os.path.join(RACMO_basedir, "mon_cum_snowmelt.nc"),
                        "SUBL"   : os.path.join(RACMO_basedir, "mon_cum_subl.nc"),
                        "TEMP_2M": os.path.join(RACMO_basedir, "mon_ave_t2m.nc")}

#############################################################################################
## MAR 3.5.2
#############################################################################################
# The dictionary for ANNUAL outputs of the various MAR datasets.  Might look into implementing
# the monthly outputs later.
MAR_basedir = os.path.join(DATA_BASEDIR, "MAR")
MAR_BASIN_INTEGERS_PICKLEFILE = os.path.join(MAR_basedir, "MAR_BASIN_INTEGERS.pickle")
MAR_POLYGON_WKTS_TEXTFILE_GREENLAND = os.path.join(MAR_basedir, "MAR_POLYGON_WKTs_GREENLAND.txt")
MAR_POLYGON_WKTS_TEXTFILE_ALL       = os.path.join(MAR_basedir, "MAR_POLYGON_WKTs_ALL.txt")
MAR_EXPORT_SHAPEFILE_FOLDER = os.path.join(MAR_basedir, "_EXPORT_SHAPEFILES")
MAR_ICEBRIDGE_SUMMARY_CSV_FILE = os.path.join(MAR_basedir, "_MAR_IceBridge_Cells_Summary.csv")
MAR_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE = os.path.join(MAR_basedir, r'_MAR_IceBridge_Traces_Rows_Cols.pickle')
MAR_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE = os.path.join(MAR_basedir, r'_MAR_Upper_Lower_Validation_indices.pickle')

# The "ERA_20C_1900-2010_20km" dataset doesn't have a yearly output.  I can create it from the monthly files, but am just skipping it at the moment.
# The "NCEPv1_1948-2015_20km" has the same issue, no yearly data set.
# For the datasets that don't originally have yearly data (only month), a picklefile is listed,
# which contains (or will contain) the annual data summarized from reading the monthly data.
# These are an entirely different file format (variable dictionary) which will need to be handled in the code.
MAR_ANNUAL_FILENAMES_DICT = {'CanESM2-histo_1950-2005_25km': os.path.join(MAR_basedir,r'CanESM2-histo_1950-2005_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-CanESM2-histo-1950-2005.pickle'),
                             'CanESM2-rcp45_2006-2100_25km': os.path.join(MAR_basedir,r'CanESM2-rcp45_2006-2100_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-CanESM2-rcp45-2006-2100.pickle'),
                             'CanESM2-rcp85_2006-2100_25km': os.path.join(MAR_basedir,r'CanESM2-rcp85_2006-2100_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-CanESM2-rcp85-2006-2100.pickle'),
                             'ERA-int_1979-2014_10km'      : os.path.join(MAR_basedir,r'ERA-int_1979-2014_10km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5.2-10km-yearly-ERA-Interim-1979-2014.pickle'),
                             'MIROC5-histo_1900-2005_25km' : os.path.join(MAR_basedir,r'MIROC5-histo_1900-2005_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-MIROC5-histo-1900-2005.pickle'),
                             'MIROC5-rcp45_2006-2100_25km' : os.path.join(MAR_basedir,r'MIROC5-rcp45_2006-2100_25km', 'monthly_outputs_interpolated_at_5km' , 'MARv3.5-yearly-MIROC5-rcp45-2006-2100.pickle'),
                             'MIROC5-rcp85_2006-2100_25km' : os.path.join(MAR_basedir,r'MIROC5-rcp85_2006-2100_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-MIROC5-rcp85-2006-2100.pickle'),
                             'NorESM1-histo_1950-2005_25km': os.path.join(MAR_basedir,r'NorESM1-histo_1950-2005_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-NorESM1-histo-1950-2005.pickle'),
                             'NorESM1-rcp45_2006-2100_25km': os.path.join(MAR_basedir,r'NorESM1-rcp45_2006-2100_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-NorESM1-rcp45-2006-2100.pickle'),
                             'NorESM1-rcp85_2006-2100_25km': os.path.join(MAR_basedir,r'NorESM1-rcp85_2006-2100_25km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5-yearly-NorESM1-rcp85-2006-2100.pickle'),
                             'ERA_20C_1900-2010_20km'      : os.path.join(MAR_basedir,r'ERA-20C_1900-2010_20km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5.2-20km-yearly-ERA-20C-1900-2010.pickle'),
                             'NCEPv1_1948-2015_20km'       : os.path.join(MAR_basedir,r'NCEPv1_1948-2015_20km', 'monthly_outputs_interpolated_at_5km', 'MARv3.5.2-20km-yearly-NCEP1-1948-2015.pickle'),
}


# Utility to open up and look at the files.
if __name__ == "__main__":
    pass