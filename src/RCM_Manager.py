EMI_VALIDATION_CUTOFF_VALUES_PICKLEFILE# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 11:37:43 2017

@author: mmacferrin
RCM_Manager.py - class dictionaries for handling gridded RCM data (HIRHAM5, RACMO, MAR).
In this case it's over the Greenland ice sheet, but should be generalized for everything.

In these classes, the following definitions are used:
DATASET - refers to a particular run/version of the model outputs.  For instance in HIRHAM
          we have a 1991-2010 GCM dataset, and a 1980-2014 ERA-Interim dataset.
"""

import os
import netCDF4
import numpy
#import datetime as dt
import matplotlib.pyplot as plt
import osgeo.ogr as ogr
import osgeo.osr as osr
import pickle
import re
import scipy.stats

#from IceBridgeGPR_Manager_v2 import IceBridgeGPR_Manager_v2

from RCM_FileData import HIRHAM_FILENAMES_2010_GCM, \
                         HIRHAM_FILENAMES_2014_ERA, \
                         HIRHAM_FILENAMES_2050_GCM_RCP45, \
                         HIRHAM_FILENAMES_2050_GCM_RCP85, \
                         HIRHAM_FILENAMES_2100_GCM_RCP45, \
                         HIRHAM_FILENAMES_2100_GCM_RCP85, \
                         HIRHAM_GRID_DATAFILE, \
                         HIRHAM_BASIN_DATAFILE, \
                         HIRHAM_MASK_DATAFILE, \
                         HIRHAM_ELEV_PICKLEFILE, \
                         HIRHAM_POLYGON_WKTS_TEXTFILE_GREENLAND, \
                         HIRHAM_POLYGON_WKTS_TEXTFILE_ALL, \
                         RACMO_FILENAMES_DICT, \
                         RACMO_GRID_INTEGERS_PICKLEFILE, \
                         RACMO_POLYGON_WKTS_TEXTFILE_GREENLAND, \
                         RACMO_POLYGON_WKTS_TEXTFILE_ALL, \
                         MAR_ANNUAL_FILENAMES_DICT, \
                         MAR_BASIN_INTEGERS_PICKLEFILE, \
                         MAR_POLYGON_WKTS_TEXTFILE_GREENLAND, \
                         MAR_POLYGON_WKTS_TEXTFILE_ALL, \
                         ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, \
                         HIRHAM_EXPORT_SHAPEFILE_FOLDER, \
                         RACMO_EXPORT_SHAPEFILE_FOLDER, \
                         MAR_EXPORT_SHAPEFILE_FOLDER, \
                         HIRHAM_ICEBRIDGE_CELLS_SUMMARY_CSV, \
                         RACMO_ICEBRIDGE_SUMMARY_CSV_FILE, \
                         MAR_ICEBRIDGE_SUMMARY_CSV_FILE, \
                         HIRHAM_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE, \
                         RACMO_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE, \
                         MAR_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE, \
                         ICEBRIDGE_THICKNESS_VALIDATION_TRACKS_LIST_FILE, \
                         HIRHAM_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE, \
                         RACMO_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE, \
                         MAR_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE, \
                         EMI_VALIDATION_CUTOFF_VALUES_PICKLEFILE, \
                         EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE, \
                         STATISTICS_OUTPUT_FOLDER

########################################################################################################
#########################################  RCM Superclass  #############################################
########################################################################################################
class RCM_Manager(object):
    def __init__(self):
        pass

    def get_manager_name(self):
        return {HIRHAM5_Manager: "HIRHAM",
                RACMO_Manager  : "RACMO" ,
                MAR_Manager    : "MAR"   }[type(self)]

    def get_variable(self, var_name, dataset_key=None, yearspan=None, annual=True):
        # Use upper-case keys.
        var_name = var_name.upper()

        if type(dataset_key) in (list, tuple):
            # If we're given multiple datsets and yearspans, concatenate the values within.


            if var_name in ("SMB", "RUNOFF", "MELT", "TIME", "EVAP_SUBL", "SUBL", "RAIN", "SNOW"):
                if yearspan is None:
                    yearspan = [None] * len(dataset_key)
                elif (len(yearspan) == 2) and (type(yearspan[0]) in (int,numpy.int,numpy.int32)):
                    yearspan = (yearspan,)

                # Sometimes we're just given one date range (say, 1951-1960) with two datset_keys.
                if len(dataset_key) != len(yearspan):
                    assert (len(yearspan) == 1) and (len(yearspan[0]) == 2) and (type(yearspan[0][0]) in (int,numpy.int,numpy.int32))

                    L = [(dset, yearspan[0]) for dset in dataset_key if self._does_yearspan_overlap_dataset_key(dset, yearspan[0])]

                    if len(L) == 1:
                        # It only overlaps one dataset
                        yearspan = (L[0][1],)
                        dataset_key = (L[0][0],)

                    else:

                        dataset_key = tuple([item[0] for item in L])
                        # If the yearspan spans both datasets, divvy it up.
                        yearspan = [None] * len(L)
                        for i,(dset, yspan) in enumerate(L):
                            TIME = self.get_variable("TIME", dataset_key = dset)
                            yearspan[i] = (max(yspan[0], TIME[0]), min(yspan[1], TIME[-1]))

                assert len(dataset_key) == len(yearspan)

                i = 0
                for (dset_key, yspan) in zip(dataset_key, yearspan):
                    if i == 0:
                        var = self.get_variable(var_name, dataset_key=dset_key, yearspan=yspan, annual=annual)
                    else:
                        var = numpy.append(var, self.get_variable(var_name, dataset_key=dset_key, yearspan=yspan, annual=annual), axis=0)
                    i += 1
            else:
                # If it's just grid variables, just get one of them (don't need to retrieve them all).
                assert var_name in ("MASK", "LAT", "LON", "CELLAREA", "TEMP", "BASINS", "ELEV")
                var = self._get_variable_SUBCLASS(var_name, dataset_key = dataset_key[0])

        else:
            var = self._get_variable_SUBCLASS(var_name, dataset_key=dataset_key, yearspan=yearspan, annual=annual)

        return var

    def _does_yearspan_overlap_dataset_key(self, dataset_key, yearspan):
        assert (len(yearspan) == 2) and (type(yearspan[0]) in (int,numpy.int,numpy.int32)) and (yearspan[1] >= yearspan[0])
        dataset_key_years_match_object = re.search("\d{4}-\d{4}", dataset_key)
        assert dataset_key_years_match_object != None
        dataset_key_years_string = dataset_key_years_match_object.group()
        dataset_key_startyear = int(dataset_key_years_string[:4])
        dataset_key_endyear = int(dataset_key_years_string[5:])
        assert dataset_key_endyear >= dataset_key_startyear

        return not ((dataset_key_startyear > yearspan[1]) or (dataset_key_endyear < yearspan[0]))

    def get_lats_lons(self):
        # Returns the center lats and lons
        return self.get_variable("LAT"), self.get_variable("LON")

    def _get_export_shapefile_folder(self):
        return {HIRHAM5_Manager: HIRHAM_EXPORT_SHAPEFILE_FOLDER,
                RACMO_Manager  : RACMO_EXPORT_SHAPEFILE_FOLDER ,
                MAR_Manager    : MAR_EXPORT_SHAPEFILE_FOLDER   }[type(self)]

    def _get_IceBridge_summary_csv_filename(self):
        return {HIRHAM5_Manager: HIRHAM_ICEBRIDGE_CELLS_SUMMARY_CSV,
                RACMO_Manager  : RACMO_ICEBRIDGE_SUMMARY_CSV_FILE  ,
                MAR_Manager    : MAR_ICEBRIDGE_SUMMARY_CSV_FILE    }[type(self)]

    def _get_icebridge_traces_cellnumbers_picklefile(self):
        return {HIRHAM5_Manager: HIRHAM_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE,
                RACMO_Manager  : RACMO_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE  ,
                MAR_Manager    : MAR_ICEBRIDGE_TRACES_CELLNUMBERS_PICKLEFILE    }[type(self)]

    def _get_icebridge_validation_upper_lower_grid_indices_picklefile_name(self):
        return {HIRHAM5_Manager: HIRHAM_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE,
                RACMO_Manager  : RACMO_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE ,
                MAR_Manager    : MAR_UPPER_LOWER_VALIDATION_INDICES_PICKLEFILE    }[type(self)]

        return EMI_VALIDATION_CUTOFF_VALUES_PICKLEFILE

    def _get_WKT_textfile_name(self, mask_greenland=False):
        if mask_greenland:
            return {HIRHAM5_Manager: HIRHAM_POLYGON_WKTS_TEXTFILE_GREENLAND,
                    RACMO_Manager  : RACMO_POLYGON_WKTS_TEXTFILE_GREENLAND,
                    MAR_Manager    : MAR_POLYGON_WKTS_TEXTFILE_GREENLAND    }[type(self)]
        else:
            return {HIRHAM5_Manager: HIRHAM_POLYGON_WKTS_TEXTFILE_ALL,
                    RACMO_Manager  : RACMO_POLYGON_WKTS_TEXTFILE_ALL ,
                    MAR_Manager    : MAR_POLYGON_WKTS_TEXTFILE_ALL    }[type(self)]

    def get_basin_names_lookup(self, short=True):
        '''The integers in each of the basins from "get_basin_integers" map to a basin/region in Greenland.
        This returns a dictionary lookup to the names. "short" gives a dictionary of abbreviations (i.e. "N", "NE", etc).
        If "short" is False, then it gives longer names ("North", "Northeast", etc).'''

        if short:
            return { 1: "N",
                     2: "NE",
                     3: "E",
                     4: "SE",
                     5: "S",
                     6: "SW",
                     7: "W",
                     8: "NW",
                     9: "FI"}

        else: ## short is False, use long names.
            return { 1: "North",
                     2: "Northeast",
                     3: "East",
                     4: "Southeast",
                     5: "South",
                     6: "Southwest",
                     7: "West",
                     8: "Northwest",
                     9: "Flade Isblink"}

    def interpolate_between_grids(self, source_RCM_Manager, var_name, src_dataset_key=None, interpolation_type="nearest neighbor"):
        '''Take a gridded variable from a source RCM manager, and interpolate it to the grid of the source dataset.
        Return the re-gridded data product.'''
        # just deal with lower-case here.
        interpolation_type = interpolation_type.lower()

        # Sanity checks
        assert interpolation_type in ("nearest neighbor", "bilinear")
        assert isinstance(source_RCM_Manager, RCM_Manager)

        # 1) Get the source dataset and variable
        src_variable = source_RCM_Manager.get_variable(var_name, dataset_key=src_dataset_key)
        # FOR NOW, just deal with 2-D "grid" type variables, for those with a third axis (time), it'd just take a bit of extra code, deal with it later if needed.
        assert len(src_variable.shape) == 2
        # 2) Get the source lat/lon
        src_lat, src_lon = source_RCM_Manager.get_lats_lons()
        # 3) Get the destination (self) lat/lon
        dst_lat, dst_lon = self.get_lats_lons()
        # 4) regrid both to polar coordinates.
        src_eastings, src_northings = source_RCM_Manager.return_coordinates_polarstereo(lats=src_lat, lons=src_lon)
        dst_eastings, dst_northings = self.return_coordinates_polarstereo(lats=dst_lat, lons=dst_lon)
        # 5) Create a blank grid for the new variable (assume 2D grid variables for now)
        dst_variable = numpy.empty(dst_lat.shape, dtype=src_variable.dtype)
        # 6) Compute the bounding boxes of the destination pixels, in eastings/northings.
        dst_bbox_northings, dst_bbox_eastings = self._compute_polygon_boundaries(X=dst_eastings, Y=dst_northings)

        if interpolation_type == "nearest neighbor":
            # For each row in the destionation variable, search for the nearest neighbor in rows in the source variable.
            for row_i in range(dst_variable.shape[0]):
                print "Row", row_i+1, "of", dst_variable.shape[0],
                dst_eastings_row = dst_eastings[row_i,:]
                dst_northings_row = dst_northings[row_i,:]
                # Assign extra dimensions to the DEST eastings and northings for array broadcasting.
                dst_eastings_row.shape = tuple(list(dst_eastings_row.shape) + [1])
                dst_northings_row.shape = tuple(list(dst_northings_row.shape) + [1])

                # Subset src_northings and src_eastings for quicker searches -- use a bounding box of the row.
                dst_row_max_easting = numpy.max(dst_bbox_eastings[row_i:row_i+2,:])
                dst_row_min_easting = numpy.min(dst_bbox_eastings[row_i:row_i+2,:])
                dst_row_max_northing = numpy.max(dst_bbox_northings[row_i:row_i+2,:])
                dst_row_min_northing = numpy.min(dst_bbox_northings[row_i:row_i+2,:])

                # These are the src pixels that are within the bounding box of this dest row (boundaries of the pixels, not just the pixel centers)
                src_subset_mask = (src_eastings >= dst_row_min_easting) & (src_eastings <= dst_row_max_easting) & \
                                  (src_northings >= dst_row_min_northing) & (src_northings <= dst_row_max_northing)
                src_subset_rows, src_subset_cols = numpy.where(src_subset_mask)

                src_eastings_subset = src_eastings[src_subset_mask]
                print src_eastings_subset.shape, "Matches.",
                if len(src_eastings_subset) == 0:
                    dst_variable[row_i,:] = 0
                else:
                    src_eastings_subset.shape = 1,src_eastings_subset.shape[0]
                    src_northings_subset = src_northings[src_subset_mask]
                    src_northings_subset.shape = 1,src_northings_subset.shape[0]

                    distances = numpy.sqrt(numpy.power(dst_eastings_row - src_eastings_subset, 2) + numpy.power(dst_northings_row - src_northings_subset, 2))
                    # Change this into a 2D array to use the argmin method, which can only go along once axis.
                    min_locs_flattened = numpy.argmin(distances, axis=1)
                    src_min_rows = src_subset_rows[min_locs_flattened]
                    src_min_cols = src_subset_cols[min_locs_flattened]

                    # Assign the destination variables to those nearest neighbors.
                    dst_variable[row_i,:] = [src_variable[src_min_rows[i], src_min_cols[i]] for i in range(len(src_min_rows))]

                print numpy.unique(dst_variable[row_i,:])

        else:
            raise ValueError("Interpolation type '{0}' not recognized.".format(interpolation_type))

        # 8) Return the new grid.
        return dst_variable

    def return_coordinates_polarstereo(self, lats, lons):
        '''Return the coordinates in North Polar Stereo projection, in (eastings, northings)'''
        lats_shape = lats.shape
        lons_shape = lons.shape
        assert lats_shape == lons_shape
        lats = lats.flatten()
        lons = lons.flatten()
        points = zip(lons, lats)
        ## 2. Convert all points to Polar Stereo grid projection
        gcsSR = osr.SpatialReference()
        gcsSR.SetWellKnownGeogCS("WGS84") # WGS84 geographical coordinates
        npsSR = osr.SpatialReference()
        npsSR.ImportFromEPSG(3413) # NSIDC North Pole Stereographic
        GCS_TO_NPS = osr.CoordinateTransformation(gcsSR, npsSR)

        # Transform to North Polar Stereo
        nps_points = numpy.array(GCS_TO_NPS.TransformPoints(points))
        # Cut off the zero "height" dimension
        eastings  = numpy.array([p[0] for p in nps_points])
        northings = numpy.array([p[1] for p in nps_points])

        eastings.shape = lons_shape
        northings.shape = lats_shape

        return eastings, northings

    def _compute_polygon_boundaries(self, Y=None, X=None):
        '''Compute the coordinate boundaries of polygons defined by each pixel in the grid.  For now just use lat/lon positions,
        in the future use a better projection (or find a projection that actually simulates the HIRHAM5 Grid).'''
        # Get the latitudes & longitudes of our center points.
        if X is None or Y is None:
            lats, lons = self.get_lats_lons()
        else:
            lats, lons = Y, X
        # Make sure lats and lons are same shape, and 2-dimensional
        assert lats.shape == lons.shape and len(lats.shape) == 2

        # Define arrays for the polygon corner lat/lons
        poly_lats = numpy.empty([d+1 for d in lats.shape], dtype=lats.dtype)
        poly_lons = numpy.empty([d+1 for d in lons.shape], dtype=lons.dtype)

        # Define the center lats & lons by simply interpolating between adjacent pixels
        lats_middle = (lats[:-1,:] + lats[1:,:]) / 2
        lons_middle = (lons[:,:-1] + lons[:,1:]) / 2

        # Fill in the middle pixels of each, and then the leading edge that was left out. We're not handling edge values carefully because they'll be masked out anyway, won't matter.
        poly_lats[1:-1,:-1] = lats_middle
        poly_lons[:-1,1:-1] = lons_middle
        # Fill in edges with semi-nonsense values
        # Right edge of poly_lats, minus top & bottom
        poly_lats[1:-1,-1] = lats_middle[:,-1]
        # Top (north) edge of poly_lons, minus left & right
        poly_lons[-1,1:-1] = lons_middle[-1,:]

        # Top edge of poly_lats
        poly_lats[-1,:] = poly_lats[-2,:] + (poly_lats[-2,:] - poly_lats[-3,:])
        # Bottom edge of poly lats
        poly_lats[ 0,:] = poly_lats[ 1,:] - (poly_lats[ 2,:] - poly_lats[ 1,:])
        # Left edge of poly_lons
        poly_lons[:, 0] = poly_lons[:, 1] - (poly_lons[:, 2] - poly_lons[:, 1])
        # Bottom edge of poly lats
        poly_lons[:,-1] = poly_lons[:,-2] + (poly_lons[:,-2] - poly_lons[:,-3])

        # These should now be the polygon boundaries of our grid.  Use them accordingly.
        return poly_lats, poly_lons

    def get_list_of_icebridge_modeled_ice_validation_tracks(self):
        '''The Modeled_Ice_Validation_Tracks_List.txt file contains a list of track names, and some comments.
        Read this list and return it.'''
        f = open(ICEBRIDGE_THICKNESS_VALIDATION_TRACKS_LIST_FILE, 'r')
        lines = [line.strip() for line in f.readlines() if (re.search("\d{8}_\d{2}_\d{3}_\d{3}",line) is not None)]
        f.close()
        return lines

    def compute_EMI_values(self, dataset_key, yearspan=None, include_rain=True, compute_mean=True, vars_given_dict=None):

        MC_ratios, MC_thresholds = self.compute_MC_ratio_and_threshold(dataset_key,
                                                                       yearspan=yearspan,
                                                                       include_rain = include_rain,
                                                                       compute_mean = compute_mean,
                                                                       vars_given_dict = vars_given_dict)

        if vars_given_dict is None:
            ACCUM = self.get_variable("SNOW", dataset_key = dataset_key, yearspan=yearspan)
            if compute_mean:
                ACCUM = numpy.mean(ACCUM, axis=0)
        else:
            ACCUM = vars_given_dict["SNOW"]

        return (MC_ratios - MC_thresholds) * ACCUM


    def compute_MC_ratio_and_threshold(self, dataset_key, yearspan=None, include_rain=True, compute_mean=True, vars_given_dict=None):
        '''Computes the right hand of the Pfeffer '91 equation, and returns the ratio that it's above or below)
        the threshold met for runoff.  Anything greater than 1.0 is "runoff", anything greater than 0 is saturating.

        If we want to feed the variables given rather than have the function retrieve it, we can use the "vars_given_dict" with the following
        keys: "MELT", "SNOW", "RAIN", "LAT", "LON", "ELEV", "TEMP" and the variables as values.'''

        if vars_given_dict is None:
            if type(dataset_key) == str:
                assert dataset_key in self._LIST_dataset_keys
            elif type(dataset_key) in (list,tuple):
                for dset_key in dataset_key:
                    assert dset_key in self._LIST_dataset_keys

            # Get key variables
            M = self.get_variable("MELT", dataset_key = dataset_key, yearspan = yearspan)
            C = self.get_variable("SNOW", dataset_key = dataset_key, yearspan = yearspan)
            if include_rain:
                R = self.get_variable("RAIN", dataset_key = dataset_key, yearspan = yearspan)
            TEMP = self.get_variable("TEMP", dataset_key = dataset_key, yearspan = yearspan)

            LAT, LON = self.get_lats_lons()
            ELEV = self.get_variable("ELEV", dataset_key = dataset_key)

        else:
            M = vars_given_dict["MELT"]
            C = vars_given_dict["SNOW"]
            if include_rain:
                R = vars_given_dict["RAIN"]
            LAT = vars_given_dict["LAT"]
            LON = vars_given_dict["LON"]
            ELEV = vars_given_dict["ELEV"]
            TEMP = vars_given_dict["TEMP"]

        t_annual_mean = numpy.mean(TEMP, axis=0)

        pi = .873 # Refrozen ice density, from Machguth 2016
        # Snow density parameterization from Langen, et al. (2017)
        a_rho = 328.35
        b_rho = -0.049376
        c_rho = 1.0427
        d_rho = -0.11186
        ps = (a_rho + b_rho*ELEV + c_rho*LAT + d_rho*LON) / 1000.0

        # Maybe try with firn as well, ps = 550?
        # A bit of shorthand to simplify the equation
        Z = (pi - ps)/ps
        c = 2.108 # Heat capacity of ice, kJ/(kg *C)
        L = 333.55 # Latent heat of fusion of ice, kJ/kg
        # The t_annual value is being used to calculate firn/snow temperatures.
        # For this calculation, null-out any average temperatures above zero (should be none for the ice sheet).
        t_annual_mean[t_annual_mean > 0.0] = 0.0
        Q = (c/L)*-1.0*t_annual_mean # Expressed in *positive* degrees below zero.

        MC_Threshold = (Q + Z)/(1 + Z)  # Right side of the Pfeffer '91 equation

        if compute_mean:
            M_mean = numpy.mean(M, axis=0)
            C_mean = numpy.mean(C, axis=0)
            if include_rain:
                R_mean = numpy.mean(R, axis=0)
                MC_Ratio = (M_mean + R_mean) / C_mean
            else:
                MC_Ratio = M_mean / C_mean
        else: # not compute_mean
            if include_rain:
                MC_Ratio = (M + R) / C
            else:
                MC_Ratio = M / C

        # Get rid of the 1-axis in here.
        if compute_mean:
            MC_Ratio.shape = [d for d in MC_Ratio.shape if d>1]
            MC_Threshold.shape = MC_Ratio.shape

        return MC_Ratio, MC_Threshold

    def _get_polygon_WKTs(self, mask=None, mask_greenland=True, overlap_nudge=True):
        '''Create a list of WKT strings to create a polygon for each pixel.  The order will start from bottom-left and go across by row,
        same as a grid.flatten() operation would.  This can be masked or not for only points over Greenland if we choose.

        Two mask options.  Either (a) mask out using a custom mask in the (mask) variable, or (b) select mask_greenland=True, and use the Greenland mask.
        If both are selected, filter it out by both.

        In the GIS software, the boundaries between these pixels show up as tiny little lines.  It's annoying.
        If overlap_nudge is True, then return the points "expanded" just a little bit.'''
        # If we can import the instead of computing them all, just do that.
        # For each _get_WKT_textfile_name() function, modify to choose either this overlap or the non-overlap-nudged textfiles.
        fname = self._get_WKT_textfile_name(mask_greenland = (mask_greenland and mask is None), overlap_nudge=overlap_nudge)
        if os.path.exists(fname):
            f = open(fname, 'r')
            lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
            f.close()

            # If we gave a custom mask, use it to subset the variables.
            if mask is not None:
                lines_array = numpy.array(lines)
                lines_array.shape = mask.shape # This SHOULD work if the masks are correctly sized

                if mask_greenland:
                    gr_mask = self.get_variable("MASK")
                    lines = list(lines_array[mask & gr_mask])
                else:
                    lines = list(lines_array[mask])

            return lines

        corner_lats, corner_lons = self._compute_polygon_boundaries()

        # Define corners
        LL_lats = corner_lats[:-1,:-1]
        LR_lats = corner_lats[:-1,1:]
        UL_lats = corner_lats[1:,:-1]
        UR_lats = corner_lats[1:,1:]

        LL_lons = corner_lons[:-1,:-1]
        LR_lons = corner_lons[:-1,1:]
        UL_lons = corner_lons[1:,:-1]
        UR_lons = corner_lons[1:,1:]

        if overlap_nudge:
            # An epsilon small enough to actually change a 32-bit floating point number with a range of 50-110.
            epsilon = 1e-4
            # Nudge the right side and the top side.  This ever-so-slightly dismorphs the shapes, but not enough to notice.
            # Hopefully just enough to make the squares slightly overlap and not have hairline gaps in the output maps.
            UL_lats = UL_lats + epsilon
            UR_lats = UR_lats + epsilon
            LR_lons = LR_lons + epsilon
            UR_lons = UR_lons + epsilon

        if mask is None and mask_greenland:
            mask = self.get_variable("MASK")

        M,N = LL_lats.shape
        WKTs = [None] * ((M * N) if (mask is None) else numpy.count_nonzero(mask))

        count = 0
        for i in range(M):
            for j in range(N):
                if mask is not None and not mask[i,j]:
                    continue
                # Define the polygon in a 5-point closed loop, counter_clockwise.
                WKTs[count] = "POLYGON (({0} {1}, {2} {3}, {4} {5}, {6} {7}, {0} {1}))".format(
                LL_lons[i,j], LL_lats[i,j], LR_lons[i,j], LR_lats[i,j], UR_lons[i,j], UR_lats[i,j], UL_lons[i,j], UL_lats[i,j])
                count += 1

        return WKTs

    def export_shapefile_WKTs(self):
        '''Create a list of polygon WKT's and export them to a text file, so we don't have to keep doing this again and again.'''
        for mask_greenland in True, False:
            for overlap_nudge in True, False:
                fname = self._get_WKT_textfile_name(mask_greenland = mask_greenland, overlap_nudge = overlap_nudge)
                WKTs = self._get_polygon_WKTs(mask_greenland=mask_greenland, overlap_nudge = overlap_nudge)
                print "Writing", os.path.split(fname)[-1]
                f = open(fname, 'w')
                f.writelines([wkt + '\n' for wkt in WKTs])
                f.close()

    def grid_icebridge_thickness_data(self):
        '''Take the raw icebridge shapefile data, and compile it into our grids.'''
        # Open the CSV file, read in the points (lat, lon, thickness.)
        csv_file = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'r')
        lines = csv_file.readlines()
        csv_file.close()

        # Get a header, and strip off any blank lines
        header = dict([(item,i) for i,item in enumerate(lines[0].strip("\n").split(","))])
        lines = [line.strip() for line in lines[1:] if len(line.strip()) > 0]

        lats = numpy.empty(len(lines), dtype=numpy.float)
        lons = numpy.empty(len(lines), dtype=numpy.float)
        ice_contents_20m = numpy.empty(len(lines), dtype=numpy.float)
        # Create a list of (lon,lat) points for all the IceBridge data
        points = [ogr.Geometry(ogr.wkbPoint) for i in range(len(lats))]

        for i,line in enumerate(lines):
            items = line.split(",")
            lats[i] = float(items[header["lat"]])
            lons[i] = float(items[header["lon"]])
            ice_contents_20m[i] = float(items[header["20m_ice_content_m"]])
            point = points[i]
            point.AddPoint(lons[i], lats[i])

        # Get only gridcells that are on the greenland ice sheet.
        greenland_mask = self.get_variable("MASK")
        greenland_indices = numpy.argwhere(greenland_mask)
        WKTs = self._get_polygon_WKTs(mask_greenland = True)
        drainages = self.get_variable("BASINS")

        ######################################
        # Create a Shapefile output.
        ######################################
        # Get the full file path
        shapefile_fname ={HIRHAM5_Manager: "HIRHAM",
                          RACMO_Manager  : "RACMO" ,
                          MAR_Manager    : "MAR"   }[type(self)] + "_IceBridge_Cells_Only.shp"
        shapefile_fname = os.path.join(self._get_export_shapefile_folder(), shapefile_fname)

        # Create the data source and file.
        driver = ogr.GetDriverByName("ESRI Shapefile")

        if os.path.exists(shapefile_fname):
            print "Deleting previous", os.path.split(shapefile_fname)[-1]
            driver.DeleteDataSource(shapefile_fname)

        data_source = driver.CreateDataSource(shapefile_fname)
        print "Creating", os.path.split(shapefile_fname)[-1]

        # Create a WGS84 geographic spatial reference to go with it.
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        ############################################################
        # Create the IceSlabs layer, with Trace_Count, mean_ice, and std_ice fields
        ############################################################
        layer = data_source.CreateLayer("IceSlabs", srs, geom_type=ogr.wkbPolygon)
        # Add the variable field we're interested in
        count_field = ogr.FieldDefn("TraceCount", ogr.OFTInteger)
        layer.CreateField(count_field)
        mean_field = ogr.FieldDefn("Ice_Mean_m", ogr.OFTReal)
        layer.CreateField(mean_field)
        std_field = ogr.FieldDefn("Ice_Std_m", ogr.OFTReal)
        layer.CreateField(std_field)
        drainage_field = ogr.FieldDefn("Drainage", ogr.OFTInteger)
        layer.CreateField(drainage_field)
        gridX_field = ogr.FieldDefn("Grid_Row", ogr.OFTInteger)
        layer.CreateField(gridX_field)
        gridY_field = ogr.FieldDefn("Grid_Col", ogr.OFTInteger)
        layer.CreateField(gridY_field)

        ################################################
        # Create a CSV File for outputting these grid statistics.
        ################################################
        csv_fname = self._get_IceBridge_summary_csv_filename()
        print "Writing", os.path.split(csv_fname)[-1]
        csv_file = open(csv_fname, 'w')
        csv_header = "Row,Column,Drainage_N,Drainage,Trace Count,Mean Thickness m,Std Thickness m\n"
        csv_file.write(csv_header)

        basin_names_dict = self.get_basin_names_lookup(short=True)

        # Loop through gridcells.
        # Any gridcells that contain icebridge points, compute the count, mean, and stddev thicknesses in that gridcell.
        count = 0

        for i, index in enumerate(greenland_indices):
            if (i % 100) == 0:
                print '.',
            polygon = ogr.CreateGeometryFromWkt(WKTs[i])
            min_lon, max_lon, min_lat, max_lat = polygon.GetEnvelope()

            # In order to keep OGR from breaking on 100s of thousands of points,
            # when we see an indication that there might be an intersection, we'll create a multipoint subset with it.
            subset_mask = ((lons >= min_lon) & \
                           (lons <= max_lon) & \
                           (lats >= min_lat) & \
                           (lats <= max_lat))
            N = numpy.count_nonzero(subset_mask)

            if N > 0:
                # Create a multi-point geometry to perform an interesection with the grid.
                multipoint = ogr.Geometry(ogr.wkbMultiPoint)
                lats_subset = lats[subset_mask]
                lons_subset = lons[subset_mask]
                # Add points within our bounding-box subset to this multipoint, rather than doing it on all the points.
                for j in range(N):
                    point = ogr.Geometry(ogr.wkbPoint)
                    point.AddPoint(lons_subset[j],lats_subset[j])
                    multipoint.AddGeometry(point)

                intersection = polygon.Intersection(multipoint)
                N = intersection.GetGeometryCount()
                if N == 0: # Our subset didn't quite work, move along)
                    continue

                # The Flade Isblink ice cap in NE Greenland is not labeled as a Greenland drainage.  Label it here, #9
                if (80.5 <= numpy.mean(lats_subset) <= 82.0) and \
                   (-20.0 <= numpy.mean(lons_subset) <= -10.0):

                    drainages[index[0],index[1]] = 9

                elif drainages[index[0], index[1]] == 0: # If our icebridge points are not on the HIRHAM ice sheet (which is a bit different
                                                         # than the one I used in the IceBridge data), omit it. This only happens in 1-2 places.
                    continue

                # Now get back to the original points that we just pulled from this, so we can get all the ice thicknesses.
                ice_thickness_subset = ice_contents_20m[subset_mask]
                ice_thickness_intersection = numpy.zeros((N,), dtype=numpy.float)
                intersection_points = [intersection.GetGeometryRef(j).GetPoint()[0:2] for j in range(N)]
                for j,point in enumerate(intersection_points):
                    p_matches = (lats_subset == point[1]) & (lons_subset == point[0])
                    assert numpy.count_nonzero(p_matches) == 1
                    ice_thickness_intersection[j] = ice_thickness_subset[p_matches][0]

                print
                print i, index, intersection.GetGeometryCount()

                feature = ogr.Feature(layer.GetLayerDefn())
                feature.SetField("TraceCount",N)
                mean_ice_thickness = numpy.mean(ice_thickness_intersection)
                feature.SetField("Ice_Mean_m",mean_ice_thickness)
                std_ice_thickness = numpy.std(ice_thickness_intersection)
                feature.SetField("Ice_Std_m",numpy.std(std_ice_thickness))
                drainage = drainages[index[0],index[1]]
                feature.SetField("Drainage",drainage)
                feature.SetField("Grid_Row", int(index[0]))
                feature.SetField("Grid_Col", int(index[1]))

                csv_row = "{0},{1},{2},{3},{4},{5},{6}\n".format( \
                            int(index[0]),
                            int(index[1]),
                            drainage,
                            basin_names_dict[drainage],
                            N,
                            mean_ice_thickness,
                            std_ice_thickness)

                csv_file.write(csv_row)

                feature.SetGeometry(polygon)
                layer.CreateFeature(feature)
                feature = None

                count+= 1

        # Destroy the data_source, saving the file
        data_source = None

        csv_file.close()

        print os.path.split(shapefile_fname)[-1], "written."
        print os.path.split(csv_fname)[-1], "written."


    def read_icebridge_grid_summary(self):
        '''Read the _IceBridge_Cells_Summary.csv" file, which has the outputs for the IceBridge data gridded onto HIRHAM5 cells.
        Return a dictionary with the column names as keys and the data as values.'''
        fname = self._get_IceBridge_summary_csv_filename()
        if not os.path.exists(fname):
            print "File", os.path.split(fname)[-1], 'not found.'
            return

        print "Reading", os.path.split(fname)[-1]
        f = open(fname, 'r')
        lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
        f.close()

        header = dict([(name,i) for (i,name) in enumerate(lines[0].split(','))])
        lines = lines[1:]

        N = len(lines)

        # Define a numpy data storage type for these rows.
        dt = numpy.dtype([('row',numpy.int,1),
                          ('col',numpy.int,1),
                          ('drainage_number',numpy.int,1),
                          ('drainage_name',numpy.string_, 2),
                          ('trace_count',numpy.int,1),
                          ('mean_ice_thickness_m',numpy.float,1),
                          ('std_ice_thickness_m',numpy.float,1)])

        data_array = numpy.empty((N,), dtype=dt)
        for i,line in enumerate(lines):
            items = line.split(",")
            data_array[i]['row']                  = int(items[header['Row']])
            data_array[i]['col']                  = int(items[header['Column']])
            data_array[i]['drainage_number']      = int(items[header['Drainage_N']])
            data_array[i]['drainage_name']        = items[header['Drainage']]
            data_array[i]['trace_count']          = int(items[header['Trace Count']])
            data_array[i]['mean_ice_thickness_m'] = float(items[header['Mean Thickness m']])
            data_array[i]['std_ice_thickness_m']  = float(items[header['Std Thickness m']])

        return data_array

    def export_data_to_shapefile(self, variable, variable_name, layer_name, filename, mask=None, variable_2=None, variable_2_name=None):
        '''Take a grid dataset (whatever it is) and export it to a shapefile.  It will create a shapefile
        polygon for each grid cell (on the greenland ice sheet) within the grid, and give that polygon a value.  This is a fairly
        inefficient way of doing this but without a good projection from which to use the data, it's the best I have so far.

        "mask" should either be None (output all the values in the dataset),
                            "greenland" (output just the values on the GrIS),
                            or a numpy ndarray of boolean values for a custom mask.
        datasets can either be a single value, or a whole list/tuple of values.
        data_names should be the same length as datasets, or a string name if the datasets is just a single dataset.
        '''
        foldername = self._get_export_shapefile_folder()

        # Create the data source and file.
        driver = ogr.GetDriverByName("ESRI Shapefile")

        filename = os.path.join(foldername, filename if filename[-4:].lower() == ".shp" else filename + ".shp")
        if os.path.exists(filename):
            print "Deleting previous", os.path.split(filename)[-1]
            driver.DeleteDataSource(filename)

        data_source = driver.CreateDataSource(filename)

        # Create a WGS84 geographic spatial reference to go with it.
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        # Get polygon WKT's
        WKTs = self._get_polygon_WKTs(mask=mask, mask_greenland=False)
        polygons = [ogr.CreateGeometryFromWkt(wkt) for wkt in WKTs]

        # Create the layer
        layer = data_source.CreateLayer(layer_name, srs, geom_type=ogr.wkbPolygon)

        # Add the variable field we're interested in
        field_datatype = {numpy.dtype("float") : ogr.OFTReal,
                          numpy.dtype("float32") : ogr.OFTReal,
                          numpy.dtype("int32") : ogr.OFTInteger,
                          numpy.int   : ogr.OFTInteger,
                          numpy.int32 : ogr.OFTInteger,
                          numpy.bool  : ogr.OFTInteger ,
                          numpy.dtype("bool") : ogr.OFTInteger }[variable.dtype]
        field = ogr.FieldDefn(variable_name, field_datatype)
        layer.CreateField(field)

        if variable_2 is not None and variable_2_name is not None:
            # Add the optional second variable field we're interested in
            field_2_datatype = {numpy.dtype("float") : ogr.OFTReal,
                          numpy.dtype("float32") : ogr.OFTReal,
                          numpy.dtype("int32") : ogr.OFTInteger,
                          numpy.int   : ogr.OFTInteger,
                          numpy.int32 : ogr.OFTInteger,
                          numpy.bool  : ogr.OFTInteger ,
                          numpy.dtype("bool") : ogr.OFTInteger }[variable.dtype]
            field = ogr.FieldDefn(variable_2_name, field_2_datatype)
            layer.CreateField(field)

        if mask is not None:
            assert ((type(mask) == numpy.ndarray) or (isinstance(mask, numpy.ma.core.MaskedArray))) and (mask.dtype == numpy.bool)
            var = variable[mask]
            if variable_2 is not None:
                var_2 = variable_2[mask]
        else:
            var = variable.flatten()
            if variable_2 is not None:
                var_2 = variable_2.flatten()

        for i,poly in enumerate(polygons):
            varnum = var[i]
            # Get variable value into a format that OGR likes (i.e. non-numpy)
            if variable.dtype in (numpy.dtype("float"), numpy.dtype("float32")):
                varnum = float(varnum)
            elif variable.dtype in (numpy.dtype("int"), numpy.dtype("int32")):
                varnum = int(varnum)
            elif variable.dtype in (numpy.dtype("bool"), numpy.bool):
                varnum = int(varnum)

            if variable_2 is not None:
                varnum_2 = var_2[i]
                if variable_2.dtype in (numpy.dtype("float"), numpy.dtype("float32")):
                    varnum_2 = float(varnum_2)
                elif variable_2.dtype in (numpy.dtype("int"), numpy.dtype("int32")):
                    varnum_2 = int(varnum_2)
                elif variable_2.dtype in (numpy.dtype("bool"), numpy.bool):
                    varnum_2 = int(varnum_2)

            # Create the feature:
            feature = ogr.Feature(layer.GetLayerDefn())
            # Set the attribute
            try:
                feature.SetField(variable_name, varnum)

                if variable_2 is not None and variable_2_name is not None:
                    feature.SetField(variable_2_name, varnum_2)

            except NotImplementedError:
                print variable_name, varnum, type(varnum)
                if variable_2 is not None:
                    print variable_2_name, varnum_2, type(varnum_2)
                raise NotImplementedError
            # Create the polygon and put it in the feature
            feature.SetGeometry(poly)
            # Put this feature in the layer
            layer.CreateFeature(feature)
            # Dereference the feature to write it out.
            feature = None

        # Now dereference the data source to save it.
        data_source = None
        print "Exported", os.path.split(filename)[-1]
        return

    def export_EMIs_by_decade(self):
        '''Run through ALL the dataset_keys in this manager, export shapefiles
        for all the approximate decades within each dataset.'''

        manager_name = self.get_manager_name()
        mask = self.get_variable("MASK")

        for dataset_key in self._LIST_dataset_keys:

            decade_yearspans = self._split_dataset_into_decade_yearspans(dataset_key)

            for yearspan in decade_yearspans:
                MRC_ratios, MRC_thresholds = self.compute_MC_ratio_and_threshold(dataset_key, yearspan=yearspan, compute_mean=True)

                print "RATIOS:", numpy.mean(MRC_ratios[mask])
                print "THRESHOLDS:", numpy.mean(MRC_thresholds[mask])
                ACCUM = self.get_variable("SNOW", dataset_key = dataset_key, yearspan = yearspan)
                ACCUM_mean = numpy.mean(ACCUM, axis=0)
                EMIs = (MRC_ratios - MRC_thresholds) * ACCUM_mean

                # Create layer name from dataset and years
                if manager_name == "MAR":
                    spos = re.search('_\d{4}-\d{4}', dataset_key)
                    dataset_key_short = dataset_key[:spos.start()]
                else:
                    dataset_key_short = dataset_key

                layername = dataset_key_short + "_{0}-{1}".format(yearspan[0],yearspan[1])

                variable_name = "ScaledEMI"

                filename = manager_name + "_EMI_" + layername + ".shp"

                self.export_data_to_shapefile(EMIs,
                                              variable_name = variable_name,
                                              layer_name = layername,
                                              filename = filename,
                                              mask=mask)

        return

    def _get_dataset_keys_and_yearspans_for_GCM_decadal_modeling(self):
        manager_name = self.get_manager_name()

        # Each year timestep
        interval = 1

        if manager_name == "HIRHAM":
            # These are each of the datasets that will be parsed.
            dataset_keys = [("2010GCM",),
                            ("2014ERA",),
                            ("2050GCM45",),
                            ("2050GCM85",),
                            ("2100GCM45",),
                            ("2100GCM85",),
                            ]
            # These are the decades of each yearspan over which to model/map.
            yearspans = [tuple([(year,year+9) for year in numpy.arange(1991,2010-9+interval,interval,dtype=numpy.int)]),
                         tuple([(year,year+9) for year in numpy.arange(1981,2014-9+interval,interval,dtype=numpy.int)]),# + [(2004,2013)]),
                         tuple([(year,year+9) for year in numpy.arange(2031,2050-9+interval,interval,dtype=numpy.int)]),
                         tuple([(year,year+9) for year in numpy.arange(2031,2050-9+interval,interval,dtype=numpy.int)]),
                         tuple([(year,year+9) for year in numpy.arange(2081,2100-9+interval,interval,dtype=numpy.int)]),
                         tuple([(year,year+9) for year in numpy.arange(2081,2100-9+interval,interval,dtype=numpy.int)]),
                        ]
            # These are the datasets from which we pull the "baseline" values
            baseline_dataset_keys = [("2010GCM"),
                                     ("2014ERA"),
                                     ("2010GCM",),
                                     ("2010GCM",),
                                     ("2010GCM",),
                                     ("2010GCM",),
                                    ]
            # This is the yearspan we define as previous "steady-state", masking out any pixels seen in the "threshold_baseline_yearspans" that overlap these.
            baseline_yearspans = [(1991,2000),
                                  (1981,1990),
                                  (1991,2000),
                                  (1991,2000),
                                  (1991,2000),
                                  (1991,2000),
                                 ]


        elif manager_name == "MAR":
            dataset_keys = [('CanESM2-histo_1950-2005_25km','CanESM2-rcp45_2006-2100_25km'),
                            ('CanESM2-histo_1950-2005_25km','CanESM2-rcp85_2006-2100_25km'),
                            ('ERA-int_1979-2014_10km',),
                            ('MIROC5-histo_1900-2005_25km','MIROC5-rcp45_2006-2100_25km'),
                            ('MIROC5-histo_1900-2005_25km','MIROC5-rcp85_2006-2100_25km'),
                            ('NorESM1-histo_1950-2005_25km','NorESM1-rcp45_2006-2100_25km'),
                            ('NorESM1-histo_1950-2005_25km','NorESM1-rcp85_2006-2100_25km'),
                            ('ERA_20C_1900-2010_20km',),
                            ('NCEPv1_1948-2015_20km',)
                            ]
            # Subsets of yearspans over which to model.
            yearspans = [tuple([(year,year+9) for year in numpy.arange(1951,2100-9+interval,interval,dtype=numpy.int)]), # CanESM2 4.5, 1951-2100
                         tuple([(year,year+9) for year in numpy.arange(1951,2100-9+interval,interval,dtype=numpy.int)]), # CanESM2 8.5, 1951-2100
                         tuple([(year,year+9) for year in numpy.arange(1981,2014-9+interval,interval,dtype=numpy.int)]), # ERA-Int, 1981-2013
                         tuple([(year,year+9) for year in numpy.arange(1901,2100-9+interval,interval,dtype=numpy.int)]), # MIROC5 4.5, 1901-2100
                         tuple([(year,year+9) for year in numpy.arange(1901,2100-9+interval,interval,dtype=numpy.int)]), # MIROC5 8.5, 1901-2100
                         tuple([(year,year+9) for year in numpy.arange(1951,2100-9+interval,interval,dtype=numpy.int)]), # NorESM1 4.5, 1951-2100
                         tuple([(year,year+9) for year in numpy.arange(1951,2100-9+interval,interval,dtype=numpy.int)]), # NorESM1 8.5, 1951-2100
                         tuple([(year,year+9) for year in numpy.arange(1901,2010-9+interval,interval,dtype=numpy.int)]), # ERA 20C, 1901-2010
                         tuple([(year,year+9) for year in numpy.arange(1951,2015-9+interval,interval,dtype=numpy.int)]), # NCEPv1, 1951-2013
                        ]
            baseline_dataset_keys = [('CanESM2-histo_1950-2005_25km','CanESM2-rcp45_2006-2100_25km'),
                                     ('CanESM2-histo_1950-2005_25km','CanESM2-rcp85_2006-2100_25km'),
                                     ('ERA-int_1979-2014_10km',),
                                     ('MIROC5-histo_1900-2005_25km','MIROC5-rcp45_2006-2100_25km'),
                                     ('MIROC5-histo_1900-2005_25km','MIROC5-rcp85_2006-2100_25km'),
                                     ('NorESM1-histo_1950-2005_25km','NorESM1-rcp45_2006-2100_25km'),
                                     ('NorESM1-histo_1950-2005_25km','NorESM1-rcp85_2006-2100_25km'),
                                     ('ERA_20C_1900-2010_20km',),
                                     ('NCEPv1_1948-2015_20km',)
                                    ]
            baseline_yearspans = [(1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                  (1981,1990),
                                 ]

        elif manager_name == "RACMO":
            dataset_keys = [("RACMO2.3_ERA-Int_1958-2015",),
                            ("RACMO2.1_HadGEM2_1971-2004","RACMO2.1_HadGEM2_RCP4.5_2005-2098"),
                            ]
            # Subsets of yearspans over which to model.
            yearspans = [tuple([(year,year+9) for year in numpy.arange(1961,2015-9+interval,dtype=numpy.int)]), # ERA-Interim, 1961-2010
                         tuple([(year,year+9) for year in numpy.arange(1971,2098-9+interval,dtype=numpy.int)]),
                        ]
            baseline_dataset_keys = [("RACMO2.3_ERA-Int_1958-2015",),
                                     ("RACMO2.1_HadGEM2_1971-2004","RACMO2.1_HadGEM2_RCP4.5_2005-2098"),
                                     ]
            baseline_yearspans = [(1981,1990),
                                  (1981,1990),
                                  ]

        return dataset_keys, yearspans, baseline_dataset_keys, baseline_yearspans

    def _get_dataset_keys_and_yearspans_for_thickness_parameterization(self):
        # THIS FUNCTION IS DEPRECATED>  I'M NOT USING THE ICE THICKNESS PARAMETERIZATIONS ANYMORE. ALL WE NEED IS THE SNOW DATA MASK.
        manager_name = self.get_manager_name()
        if manager_name == "HIRHAM":
            dataset_keys = [("2010GCM",),
                            ("2014ERA",)]
            # The list of "yearspans" dictionaries define, for each icebridge year, which years should be pulled from each dataset in the dataset_keys list.
            yearspans = [{2010:((2000,2009),), 2011:((2001,2010),), 2012: ((2001,2010),), 2013: ((2001,2010),), 2014: ((2001,2010),)},
                         {2010:((2000,2009),), 2011:((2001,2010),), 2012: ((2002,2011),), 2013: ((2003,2012),), 2014: ((2004,2013),)}]
            baseline_yearspans = [((1991,2000),),
                                  ((1980,1988),)]

        elif manager_name == "RACMO":
            dataset_keys = [("RACMO2.3_ERA-Int_1958-2015",)]
            # The list of "yearspans" dictionaries define, for each icebridge year, which years should be pulled from each dataset in the dataset_keys list.
            yearspans = [{2010:((2000,2009),), 2011:((2001,2010),), 2012: ((2002,2011),), 2013: ((2003,2012),), 2014: ((2004,2013),)}]
            baseline_yearspans = [((1976,1985),)]

        elif manager_name == "MAR":
            dataset_keys = [('CanESM2-histo_1950-2005_25km','CanESM2-rcp45_2006-2100_25km'),
                            ('CanESM2-histo_1950-2005_25km','CanESM2-rcp85_2006-2100_25km'),
                            ('ERA-int_1979-2014_10km',),
                            ('MIROC5-histo_1900-2005_25km','MIROC5-rcp45_2006-2100_25km'),
                            ('MIROC5-histo_1900-2005_25km','MIROC5-rcp85_2006-2100_25km'),
                            ('NorESM1-histo_1950-2005_25km','NorESM1-rcp45_2006-2100_25km'),
                            ('NorESM1-histo_1950-2005_25km','NorESM1-rcp85_2006-2100_25km'),
                            ('NCEPv1_1948-2015_20km',),
                            ]
            # The list of "yearspans" dictionaries define, for each icebridge year, which years should be pulled from each dataset in the dataset_keys list.
            yearspans = [{2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2009),), 2011:((2001,2010),), 2012: ((2002,2011),), 2013: ((2003,2012),), 2014: ((2004,2013),)},
                         {2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2005),(2006,2009)), 2011:((2001,2005),(2006,2010)), 2012: ((2002,2005),(2006,2011)), 2013: ((2003,2005),(2006,2012)), 2014: ((2004,2005),(2006,2013))},
                         {2010:((2000,2009),), 2011:((2001,2010),), 2012: ((2002,2011),), 2013: ((2003,2012),), 2014: ((2004,2013),)}]

            baseline_yearspans = [((1976,1985),),
                                  ((1976,1985),),
                                  ((1979,1988),),
                                  ((1976,1985),),
                                  ((1976,1985),),
                                  ((1976,1985),),
                                  ((1976,1985),),
                                  ((1976,1985),)]


        else:
            raise ValueError("Uknown model type: " + manager_name)

        return dataset_keys, yearspans, baseline_yearspans


    def _identify_individual_track_spans(self, track_tracenums, gap_pixel_cutoff = 5, span_pixel_cutoff = 20):
        '''A sub-function used by ::find_ice_thickness_minmax_indices_along_validation_tracks() to
        identify where there are "gaps" in a given icebridge track.
        Methodology:
            1) Ignore gaps smaller than "overlook_cutoff", just join the two contiguous sections together and call it a single section.
            2) Find the length of each section (in pixels).  Fill in any gaps that are smaller than the two sections on either side of them.
            3) Iterate over sections and gaps until they settle.
        Return a dict:
            A list of tuples, designating filtered (both for gaps and outlier points) spans with "start" and "end" pixel indices of track spans.
        '''
        assert gap_pixel_cutoff > 1
        assert span_pixel_cutoff > 1
        # track_gap_flags_v1 is a boolean array indicating whether the break in that particular track section is greater than the overlook_cutoff.
        track_raw_gap_sizes = numpy.append([1], (track_tracenums[1:] - track_tracenums[:-1]))
        track_gap_i_v1 = numpy.where(track_raw_gap_sizes > gap_pixel_cutoff)[0]

        section_lengths = []
        section_start_end_tracenums = []
        gap_lengths = []

        last_tracenum = track_tracenums[0]
        for tgi in track_gap_i_v1:
            tracenum = track_tracenums[tgi]
            assert tracenum > (last_tracenum + gap_pixel_cutoff)

            section_start_end_tracenums.append((last_tracenum, tracenum - track_raw_gap_sizes[tgi]))
            section_lengths.append(tracenum - track_raw_gap_sizes[tgi] - last_tracenum + 1)
            gap_lengths.append(track_raw_gap_sizes[tgi])


            last_tracenum = tracenum

        section_start_end_tracenums.append((last_tracenum, track_tracenums[-1]))
        section_lengths.append(track_tracenums[-1] - last_tracenum + 1)
        gap_lengths.append(None)

        last_sections_count = len(section_lengths) + 1
        # Keep iterating over track sections, eliminating "small" gaps, until no more small gaps (bigger than the sections on one of the sides) appear.
        while len(section_lengths) < last_sections_count:
            last_sections_count = len(section_lengths)
            gaps_to_eliminate = []
            # Starting with the last and the first sections, look ahead by "lookahead_lag" spots
            # to see if there are any sum of gaps that are smaller than both of the sections on either side of it.
            # In this case, get rid of any of the gaps in between.
            for lookahead_lag in range(len(section_lengths)-1, 0, -1):
                for i in range(len(section_lengths) - lookahead_lag):
                    start_section_len = section_lengths[i]
                    end_section_len = section_lengths[i+lookahead_lag]
                    gap_len_sum = numpy.sum(gap_lengths[i:i+lookahead_lag])
                    if (gap_len_sum <= start_section_len) and (gap_len_sum <= end_section_len):
                        gaps_to_eliminate.extend(range(i,i+lookahead_lag))

            gaps_to_eliminate = numpy.unique(gaps_to_eliminate)

            if len(gaps_to_eliminate) > 0:
                new_section_start = section_start_end_tracenums[0][0]
                new_section_start_end_tracenums = []
                new_section_lengths = []
                new_gap_lengths = []
                for i in range(len(section_lengths)):
                    if i not in gaps_to_eliminate:
                        # We keep this gap, end the section
                        new_section_end = section_start_end_tracenums[i][1]
                        new_section_start_end_tracenums.append((new_section_start, new_section_end))
                        new_section_lengths.append(new_section_end - new_section_start + 1)
                        new_gap_lengths.append(gap_lengths[i])
                        if i < (len(section_lengths)-1):
                            new_section_start = section_start_end_tracenums[i+1][0]

                section_start_end_tracenums = new_section_start_end_tracenums
                section_lengths = new_section_lengths
                gap_lengths = new_gap_lengths


        # Now look at short section lengths and eliminate them (len < overlook_cutoff)
        # Look at short sections (1-5 pixels) and eliminate them.  Get rid of the noise!  Add to the gaps if it's in between.
        sections_to_keep = numpy.where(numpy.array(section_lengths) > span_pixel_cutoff)[0]
        new_section_start_end_tracenums = []
        new_section_lengths = []

        for i in range(len(section_lengths)):
            if i in sections_to_keep:
                new_section_start_end_tracenums.append(section_start_end_tracenums[i])
                new_section_lengths.append(section_lengths[i])

        section_start_end_tracenums = new_section_start_end_tracenums
        section_lengths = new_section_lengths

        gap_lengths = [(section_start_end_tracenums[i+1][0] - section_start_end_tracenums[i][1]) for i in range(len(section_lengths)-1)]
        gap_lengths.append(None)

        return new_section_start_end_tracenums


    def _find_ice_thickness_minmax_indices_along_validation_tracks(self, export=True):
        '''Use the ICEBRIDGE_THICKNESS_VALIDATION_TRACKS_LIST_FILE file to
        determine which ice thickness pixels in the output files are being used to determine where ice slabs exist and where they don't.
        This creates a file with a set of pixel coordinates that the function can use.'''
        validation_picklefilename = self._get_icebridge_validation_upper_lower_grid_indices_picklefile_name()
        if (not export) and os.path.exists(validation_picklefilename):
            print "Reading", os.path.split(validation_picklefilename)[1]
            f = open(validation_picklefilename, 'r')
            data = pickle.load(f)
            f.close()
            return data

        flightlines = self.get_list_of_icebridge_modeled_ice_validation_tracks()
        icebridge_csv_data = self.import_icebridge_csv_data()
        icebridge_rows, icebridge_cols = self.get_icebridge_traces_grid_indices()
        icebridge_tracenumbers = icebridge_csv_data["Tracenumber"]
        icebridge_track_names = icebridge_csv_data["Track_name"]
        ELEV = self.get_variable("ELEV")

        main_span_lower_grid_indices_rows = numpy.empty((len(flightlines),), dtype=numpy.int)
        main_span_lower_grid_indices_cols = numpy.empty((len(flightlines),), dtype=numpy.int)
        main_span_upper_grid_indices_rows = numpy.empty((len(flightlines),), dtype=numpy.int)
        main_span_upper_grid_indices_cols = numpy.empty((len(flightlines),), dtype=numpy.int)
        outlier_grid_indices_rows         = numpy.empty((len(flightlines),), dtype=numpy.int)
        outlier_grid_indices_cols         = numpy.empty((len(flightlines),), dtype=numpy.int)

        for i,track_name in enumerate(flightlines):
            print track_name
            track_mask = (icebridge_track_names == track_name)
            track_grid_rows, track_grid_cols = icebridge_rows[track_mask], icebridge_cols[track_mask]
            track_elevs = ELEV[track_grid_rows, track_grid_cols]

            track_tracenums = icebridge_tracenumbers[track_mask]

            track_span_indices = self._identify_individual_track_spans(track_tracenums)
            # First, find the longest "main" span
            span_lengths = [(s[1] - s[0] + 1) for s in track_span_indices]

            longest_span_index = numpy.argmax(span_lengths)
            main_span_indices = track_span_indices[longest_span_index]

            # Now to determine if this is an "uphill" track or a "downhill" track, based on elevation.
            left_indices = numpy.array((main_span_indices[0],)) \
                           if (longest_span_index==0) \
                           else numpy.append(numpy.array(track_span_indices[0:longest_span_index]).flatten(), main_span_indices[0])
            right_indices = numpy.array((main_span_indices[1],)) \
                            if (longest_span_index==(len(track_span_indices)-1)) \
                            else numpy.append(main_span_indices[1], numpy.array(track_span_indices[longest_span_index+1:]).flatten())
            # We need to map the ORIGINAL icebridge indices to the subset of the indices just along this track.  Use the track_tracenums to do this.
            # Reformat shapes for array broadcasting
            track_tracenums.shape = track_tracenums.shape[0], 1
            left_indices.shape    = 1, left_indices.shape[0]
            right_indices.shape   = 1, right_indices.shape[0]

            left_indices_remapped = numpy.where(left_indices == track_tracenums)[0].flatten()
            right_indices_remapped = numpy.where(right_indices == track_tracenums)[0].flatten()
            # Revert back to original shapes
            left_indices = left_indices.flatten()
            right_indices = right_indices.flatten()

            left_elevations  = track_elevs[left_indices_remapped]
            right_elevations = track_elevs[right_indices_remapped]
            left_mean_elevation = numpy.mean(left_elevations)
            right_mean_elevation = numpy.mean(right_elevations)

            if left_mean_elevation > right_mean_elevation:
                lower_main_span_index = right_indices_remapped[0]
                upper_main_span_index = left_indices_remapped[-1]
                if (len(left_indices) > 1):
                    outlier_index = numpy.min(left_indices_remapped)
                else:
                    outlier_index = None
            elif right_mean_elevation > left_mean_elevation:
                lower_main_span_index = left_indices_remapped[-1]
                upper_main_span_index = right_indices_remapped[0]
                if (len(right_indices) > 1):
                    outlier_index = numpy.max(right_indices_remapped)
                else:
                    outlier_index = None
            else:
                # SHould not have a level transect, handling functions not implemented here.
                assert left_mean_elevation == right_mean_elevation
                print "~~~~~~~~~~~~ LEVEL??? ~~~~~~~~~~~~~~~~"
                assert False

            main_span_lower_grid_indices_rows[i] = track_grid_rows[lower_main_span_index]
            main_span_lower_grid_indices_cols[i] = track_grid_cols[lower_main_span_index]
            main_span_upper_grid_indices_rows[i] = track_grid_rows[upper_main_span_index]
            main_span_upper_grid_indices_cols[i] = track_grid_cols[upper_main_span_index]
            outlier_grid_indices_rows[i]         = track_grid_rows[outlier_index] if (outlier_index is not None) else -9999
            outlier_grid_indices_cols[i]         = track_grid_cols[outlier_index] if (outlier_index is not None) else -9999

        # 3: Record the INDEX (row,col) of (a) bottom of main track, (b) top of main track, (c) top of highest outlier section
        data = (main_span_lower_grid_indices_rows,
                main_span_lower_grid_indices_cols,
                main_span_upper_grid_indices_rows,
                main_span_upper_grid_indices_cols,
                outlier_grid_indices_rows,
                outlier_grid_indices_cols)

        if export:
            print "Writing", os.path.split(validation_picklefilename)[1]
            f = open(validation_picklefilename, 'w')
            pickle.dump(data, f)
            f.close()

        return data


    def import_icebridge_csv_data(self):
        '''Read the ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, return a numpy array with the data.'''
        print "Reading", os.path.split(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE)[-1]
        csv_file = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'r')
        lines = [line.strip() for line in csv_file.readlines() if (len(line.strip()) > 0)]
        csv_file.close()
        # Get the header and the lines
        header_list = lines[0].split(",")
        lines = lines[1:]
        dt = numpy.dtype([(header_list[0], numpy.string_, 19),
                          (header_list[1], numpy.int    , 1 ),
                          (header_list[2], numpy.float64, 1 ),
                          (header_list[3], numpy.float64, 1 ),
                          (header_list[4], numpy.float64, 1 ),
                          (header_list[5], numpy.float64, 1 )
                          ])
        # Ingest the data lines
        ib_data = numpy.empty((len(lines),), dtype=dt)
        for i,line in enumerate(lines):
            items = line.split(',')
            ib_data[i] = (items[0], int(items[1]), float(items[2]), float(items[3]), float(items[4]), float(items[5]))

        return ib_data

    def is_inside_polygon(self, points_x, points_y, polygon_x, polygon_y):
        '''A utility function for seeing whether a series of points are within our outside a polygon.
        Returns a boolean array of size points_x that tells whether each point is or isn't inside the polygon.
        polygon_x and polylgon_y can either be an N-length vector, or an MxN length array of M polygons of N points each.
        The output array will be boolean, shaped (len(points),M)'''
        n = len(polygon_x)
        inside = numpy.zeros(points_x.shape, dtype=numpy.bool)

        p1x,p1y = polygon_x[0], polygon_y[0]
        for i in range(n+1):
            p2x = polygon_x[i % n]
            p2y = polygon_y[i % n]
            mask = (points_y > min(p1y, p2y)) & (points_y <= max(p1y, p2y)) & (points_x <= max(p1x, p2x))
            if p1y != p2y:
                xints = (points_y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x:
                inside[mask] = ~inside[mask]
            else:
                mask = mask & (points_x <= xints)
                inside[mask] = ~inside[mask]
            p1x,p1y = p2x,p2y

        return inside

    def get_icebridge_traces_grid_indices(self, export=False):
        if (not export):
            picklefile_fname = self._get_icebridge_traces_cellnumbers_picklefile()
            if os.path.exists(picklefile_fname):
                print "Reading", os.path.split(picklefile_fname)[-1]
                f = open(picklefile_fname, 'r')
                output_pixel_rows, output_pixel_cols = pickle.load(f)
                f.close()
                return output_pixel_rows, output_pixel_cols

        icebridge_csv_data = self.import_icebridge_csv_data()
        output_pixel_rows, output_pixel_cols = self.compute_grid_indices_of_data_points(icebridge_csv_data['lat'], icebridge_csv_data['lon'])

        if export:
            picklefile_fname = self._get_icebridge_traces_cellnumbers_picklefile()
            f = open(picklefile_fname, 'w')
            pickle.dump((output_pixel_rows, output_pixel_cols), f)
            f.close()
            print "Exporting", os.path.split(picklefile_fname)[-1]

        return output_pixel_rows, output_pixel_cols


    def compute_grid_indices_of_data_points(self, input_lats, input_lons):
        '''Given a list of lat/lon points, return the row and column indices where each data point resides.'''

        assert (input_lats.shape == input_lons.shape) and (len(input_lats.shape)==1)
        input_indices = numpy.arange(input_lats.shape[0],dtype=numpy.int)

        input_lats_max = numpy.max(input_lats)
        input_lats_min = numpy.min(input_lats)
        input_lons_max = numpy.max(input_lons)
        input_lons_min = numpy.min(input_lons)

        # Assign blank values to -1
        output_pixel_rows = numpy.zeros((len(input_lats),), dtype=numpy.int) - 1
        output_pixel_cols = numpy.zeros((len(input_lats),), dtype=numpy.int) - 1

        # 2: Get the pixels of the RCM that overlie that track, including a bit before and after the ice in the track.
        RCM_boundary_lats, RCM_boundary_lons = self._compute_polygon_boundaries()

        mask = self.get_variable("MASK")

        for row in range(mask.shape[0]):

            for col in range(mask.shape[1]):
                poly_lats = RCM_boundary_lats[row:row+2, col:col+2].flatten()
                poly_lons = RCM_boundary_lons[row:row+2, col:col+2].flatten()
                # Switch the last two points to make it the outline of a square (instead of angling across the middle)
                poly_lats[2], poly_lats[3] = poly_lats[3], poly_lats[2]
                poly_lons[2], poly_lons[3] = poly_lons[3], poly_lons[2]

                poly_lat_max = numpy.max(poly_lats)
                poly_lat_min = numpy.min(poly_lats)
                poly_lon_max = numpy.max(poly_lons)
                poly_lon_min = numpy.min(poly_lons)

                # If this polygon is outside of where any of the icebridge points are, just skip it and move on.
                if (poly_lat_min > input_lats_max) or (poly_lat_max < input_lats_min) or \
                   (poly_lon_min > input_lons_max) or (poly_lon_max < input_lons_min):
                    continue

                # Skip this pixel if no points lie within its bounding box.
                bbox_mask = (input_lats  > poly_lat_min) & \
                            (input_lats <= poly_lat_max) & \
                            (input_lons  > poly_lon_min) & \
                            (input_lons <= poly_lon_max)
                if ((numpy.count_nonzero(bbox_mask) == 0)):
                    continue
                else:
                    input_indices_subset = input_indices[bbox_mask]


                # Now find the points that are within this pixel.
                inside_mask = self.is_inside_polygon(points_x = input_lons[bbox_mask],
                                                     points_y = input_lats[bbox_mask],
                                                     polygon_x = poly_lons,
                                                     polygon_y = poly_lats)

                count = numpy.count_nonzero(inside_mask)
                if count > 0:
                    output_pixel_rows[input_indices_subset[inside_mask]] = row
                    output_pixel_cols[input_indices_subset[inside_mask]] = col

        assert numpy.count_nonzero(output_pixel_rows == -1) == 0

        return output_pixel_rows, output_pixel_cols

    def _max_minus_outliers(self, data):
        '''Small utility function used by ::export_EMI_thresholded_maps(). Remove the outliers from a dataset and then take the minimum of the remaining points.
        Here, outliers are defined as 1.5 x Inter-quartile range from the median.  Need only remove upper outliers in this case.'''
        q1,median,q3 = numpy.percentile(data, (25,50,75))
        iqr = q3 - q1
        upper_cutoff = median + 1.5*iqr
        return numpy.max(data[data <= upper_cutoff])

    def compute_EMI_of_IceBridge_pixels(self, dataset_key):
        years = self.get_variable("TIME", dataset_key = dataset_key)
        EMI = self.compute_EMI_values(dataset_key = dataset_key, compute_mean = False)

        icebridge_rows, icebridge_cols = self.get_icebridge_traces_grid_indices(export=False)
        icebridge_mask = numpy.zeros(EMI.shape, dtype=numpy.bool)
        icebridge_mask[:,icebridge_rows, icebridge_cols] = True

        EMI_masked = EMI[icebridge_mask]
        EMI_masked.shape = years.shape[0], EMI_masked.size / years.shape[0]

        return years, EMI_masked

    def compute_ice_layers_and_excess_melt_runoff_V3(self,
                                                     csv_output_file=EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE,
                                                     forget_header=False,
                                                     skip_shapefiles=False):
        '''This builds upon the excess melt calculations above.
        Using the decision tree outlined in my research notebook.
        csv_output_file can be either a string to a file name, or an open file object with "write" or "append" permissions.
        '''
        # Code up my "decision tree v2.0" in here.

        # Create a new shapefile every *this many* years in the dataset.
        shapefile_every_N_years = 5

        # Get all my metadata information set up. Which datasets to use, what yearspans,
        # and the same for the baseline periods in each set
        dataset_keys_list,     \
        yearspans_list,        \
        baseline_dataset_keys, \
        baseline_yearspans_list = self._get_dataset_keys_and_yearspans_for_GCM_decadal_modeling()

        # Mask for the Greenland ice sheet.
        ICE_SHEET_MASK = self.get_variable("MASK", dataset_key = dataset_keys_list[0])

        # Cutoffs for excess melt calculations
        median_max_continuous_EM,        \
        median_min_continuous_EM,        \
        median_min_intermittent_EM, \
        mean_snow_95_cutoff = self.find_median_EMI_lower_upper_intermittent_and_snow_cutoffs()

        # On-the-ground area of all the grid cells.
        CELLAREA = self.get_variable("CELLAREA")

        if type(csv_output_file) == str:
            f = open(csv_output_file, 'w')
        elif type(csv_output_file) == file and not csv_output_file.closed:
            f = csv_output_file

        # Header for the CSV file output
        if forget_header:
            TOTAL_DATA_STRING = ""
        else:
            TOTAL_HEADER_STRING = ",".join(["DATASET",
                                          "YEAR_START",
                                          "YEAR_END",
                                          "FIRN_AQUIFER_POTENTIAL_PIXELCOUNT",
                                          "FIRN_AQUIFER_POTENTIAL_AREA_KM2",
                                          "FIRN_AQUIFER_PIXELCOUNT",
                                          "FIRN_AQUIFER_AREA_KM2",
                                          "ICE_SLAB_INTERMITTENT_PIXELCOUNT",
                                          "ICE_SLAB_INTERMITTENT_AREA_KM2",
                                          "ICE_SLAB_CONTINUOUS_PIXELCOUNT",
                                          "ICE_SLAB_CONTINUOUS_AREA_KM2",
                                          "STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2",
                                          "PORESPACE_ATOP_ICE_SLABS_GT_START_0000",
                                          "PORESPACE_ATOP_ICE_SLABS_GT_START_2000",
                                          "ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000",
                                          "ICE_SLAB_RUNOFF_ANNUAL_GT_START_2000",
                                          "ICE_SLAB_RUNOFF_TOTAL_GT_START_0000",
                                          "ICE_SLAB_RUNOFF_TOTAL_GT_START_2000",
                                          "INTERIOR_RUNOFF_AREA_RCM_COMPUTED_KM2",
                                          "ICE_SLAB_RUNOFF_ANNUAL_RCM_COMPUTED_GT",
                                          "ICE_SLAB_RUNOFF_TOTAL_RCM_COMPUTED_GT",
                                          # These variables aren't what we compute,
                                          # they're just what the RCM models say.
                                          "RCM_RUNOFF_AREA_STEADYSTATE_KM2",
                                          "RCM_RUNOFF_AREA_ANNUAL_KM2",
                                          "RCM_RUNOFF_AREA_CUMULATIVE_KM2",
                                          "RCM_RUNOFF_ANNUAL_FROM_STEADYSTATE_AREA_GT",
                                          "RCM_RUNOFF_CUMULATIVE_FROM_STEADYSTATE_AREA_GT",
                                          "RCM_RUNOFF_ANNUAL_FROM_INTERIOR_GT",
                                          "RCM_RUNOFF_CUMULATIVE_FROM_INTERIOR_GT"]
                                          ) + "\n"

            f.write(TOTAL_HEADER_STRING)
            TOTAL_DATA_STRING = TOTAL_HEADER_STRING

        # Loop through all the datasets, every one
        for i,dataset_key in enumerate(dataset_keys_list):

            # Get the baseline (steady-state past) year spans
            # Get the baseline dataset and forward-looking yearspans
            baseline_dataset_key = baseline_dataset_keys[i]
            baseline_yearspan = baseline_yearspans_list[i]
            yearspans = yearspans_list[i]

            # Copute the Exces Melt of all pixels in the baseline, to see what was already runoff zone prior to baseline period.
            steady_state_EM = self.compute_EMI_values(baseline_dataset_key,
                                                      yearspan=baseline_yearspan,
                                                      include_rain=True,
                                                      compute_mean=True)

            # Calculate accumulation to get aquifer areas
            steady_state_ACCUM = self.get_variable("SNOW",
                                                   dataset_key = baseline_dataset_key,
                                                   yearspan=baseline_yearspan)
            steady_state_ACCUM_mean = numpy.mean(steady_state_ACCUM, axis=0)

            # Calculate melt to get aquifer areas that will actually have aquifers.
            # This 200 mmWE annual melt is an approximate value from (Kuipers-Munneke, 2014, GRL)
            steady_state_MELT = self.get_variable("MELT",
                                                  dataset_key = baseline_dataset_key,
                                                  yearspan = baseline_yearspan)
            steady_state_RAIN = self.get_variable("RAIN",
                                                  dataset_key = baseline_dataset_key,
                                                  yearspan = baseline_yearspan)

            # For comparison, compute the runoff as calculated by the RCM, not by our model.
            rcm_steady_state_RUNOFF = self.get_variable("RUNOFF",
                                                        dataset_key = baseline_dataset_key,
                                                        yearspan = baseline_yearspan)

            rcm_steady_state_RUNOFF = numpy.sum(rcm_steady_state_RUNOFF, axis=0)
            rcm_steady_state_runoff_MASK = (rcm_steady_state_RUNOFF > 0.0) & ICE_SHEET_MASK
            rcm_cumulative_runoff_MASK = numpy.copy(rcm_steady_state_runoff_MASK)

            RCM_RUNOFF_AREA_STEADYSTATE_KM2 = numpy.sum(CELLAREA[rcm_steady_state_runoff_MASK])

            # Compute the mean melt and rain for these periods
            steady_state_MELT_mean = numpy.mean(steady_state_MELT, axis=0)
            steady_state_RAIN_mean = numpy.mean(steady_state_RAIN, axis=0)
            steady_state_MELT_gt_200_MASK = (steady_state_MELT_mean + steady_state_RAIN_mean) >= 200

            #####################
            ## MAIN STEM
            #####################

            # Make AQUIFER_MASK if snow > mean_snow_95_cutoff
            steady_state_ACCUM_mask = steady_state_ACCUM_mean >= mean_snow_95_cutoff
            # AQUIFER_POTENTIAL_MASK refers to the areas that get enough accumulation to form an aquifer, regardless of melt
            AQUIFER_POTENTIAL_MASK = ICE_SHEET_MASK & steady_state_ACCUM_mask
            # AQUIFER_MASK refers to the areas that get enough accumulation and annual melt to actually form an aquifer.
            AQUIFER_MASK = AQUIFER_POTENTIAL_MASK & steady_state_MELT_gt_200_MASK

            # Make CUMULATIVE_PORESPACE_ATOP_ICE_SLABS mask, set to zero to initialize (since no perched ice slabs *yet*).
            # Given that some of the ice slabs form SUBSURFACE, let there be two bounds. They start with 0.0 porespace atop
            # them, and then that they start with 2.0 m porespace atop them (approx 3-4 m snow)
            PORESPACE_ATOP_ICE_SLABS_START_0000 = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.float64)
            PORESPACE_ATOP_ICE_SLABS_START_2000 = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.float64)
            # Same for continuous ice slab mask... in steady state, this is zero.
            CONTINUOUS_ICE_SLAB_MASK = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.bool)

            # Make STEADY_STATE_RUNOFF_ZONE_MASK mask from steady state EM values (>median_max_continuous_EM)
            STEADY_STATE_RUNOFF_ZONE_MASK = ICE_SHEET_MASK & (steady_state_EM > median_min_intermittent_EM)
            STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2 = numpy.sum(CELLAREA[STEADY_STATE_RUNOFF_ZONE_MASK])

            STEADY_STATE_ACCUM_ZONE_MASK = ICE_SHEET_MASK & ~STEADY_STATE_RUNOFF_ZONE_MASK


            RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_0000 = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.float64)
            RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_2000 = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.float64)
            RUNOFF_CUMULATIVE_SUM_START_0000 = 0.0
            RUNOFF_CUMULATIVE_SUM_START_2000 = 0.0

            INTERMITTENT_ICE_SLAB_MASK = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.bool)

            ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_TOTAL = 0.0
            ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID  = numpy.zeros(ICE_SHEET_MASK.shape, dtype=numpy.float64)

            RCM_RUNOFF_CUMULATIVE_FROM_STEADYSTATE_AREA_GT = 0.0
            RCM_RUNOFF_CUMULATIVE_FROM_INTERIOR_GT = 0.0

            # Loop through all the yearspans now, for each year.
            for yearspan in yearspans:

                # Since we've changed this to give only change in area since the baseline, we're not bothering to model before the baseline.
                if yearspan[1] < baseline_yearspan[1]:
                    continue

                # Get 10-year means to compute ice slabs
                # Compute 1-year excess melt.
                EM_single_year = self.compute_EMI_values(dataset_key = dataset_key,
                                                         yearspan = (yearspan[-1],yearspan[-1]),
                                                         compute_mean = False)
                # This comes out as a 3D array, we just want 2D
                EM_single_year.shape = EM_single_year.shape[1:]

                EM_10y_mean = self.compute_EMI_values(dataset_key = dataset_key,
                                                                  yearspan = yearspan,
                                                                  compute_mean = True)

                # Get 10-year mean totals for snow, melt, rain

                SNOW_10y = self.get_variable("SNOW", dataset_key = dataset_key,
                                                             yearspan = yearspan)
                SNOW_10y_mean = numpy.mean(SNOW_10y, axis=0)

                # RAIN
                RAIN_10y = self.get_variable("RAIN", dataset_key = dataset_key,
                                                             yearspan = yearspan)
                RAIN_10y_mean = numpy.mean(RAIN_10y, axis=0)

                # MELT
                MELT_10y = self.get_variable("MELT", dataset_key = dataset_key,
                                                             yearspan = yearspan)
                MELT_10y_mean = numpy.mean(MELT_10y, axis=0)

                ###############################################################################
                ## RUNOFF AS COMPUTED BY THE RCM directly, not our empirical parameterization.
                ###############################################################################

                RUNOFF_RCM_COMPUTED_single_year = self.get_variable("RUNOFF",
                                                                    dataset_key = dataset_key,
                                                                    yearspan = (yearspan[-1], yearspan[-1]))
                # This comes out as a 3D array, we just want 2D
                RUNOFF_RCM_COMPUTED_single_year.shape = RUNOFF_RCM_COMPUTED_single_year.shape[1:]

                RCM_RUNOFF_SINGLE_YEAR_MASK = (RUNOFF_RCM_COMPUTED_single_year > 0) & ICE_SHEET_MASK
                rcm_single_year_interior_runoff_MASK = RCM_RUNOFF_SINGLE_YEAR_MASK & ~rcm_steady_state_runoff_MASK

                RCM_RUNOFF_AREA_ANNUAL_KM2 = numpy.sum(CELLAREA[RCM_RUNOFF_SINGLE_YEAR_MASK])
                rcm_cumulative_runoff_MASK = rcm_cumulative_runoff_MASK | RCM_RUNOFF_SINGLE_YEAR_MASK
                RCM_RUNOFF_AREA_CUMULATIVE_KM2 = numpy.sum(CELLAREA[rcm_cumulative_runoff_MASK])
                RCM_RUNOFF_ANNUAL_FROM_STEADYSTATE_AREA_GT = \
                    numpy.sum(RUNOFF_RCM_COMPUTED_single_year[RCM_RUNOFF_SINGLE_YEAR_MASK] * CELLAREA[RCM_RUNOFF_SINGLE_YEAR_MASK] * 1e-6)
                RCM_RUNOFF_CUMULATIVE_FROM_STEADYSTATE_AREA_GT += RCM_RUNOFF_ANNUAL_FROM_STEADYSTATE_AREA_GT

                RCM_RUNOFF_ANNUAL_FROM_INTERIOR_GT = \
                    numpy.sum(RUNOFF_RCM_COMPUTED_single_year[rcm_single_year_interior_runoff_MASK] * CELLAREA[rcm_single_year_interior_runoff_MASK] * 1e-6)

                RCM_RUNOFF_CUMULATIVE_FROM_INTERIOR_GT += RCM_RUNOFF_ANNUAL_FROM_INTERIOR_GT

                ################################################################################
                # BRANCH #1: Wasn't continuous ice slab, but melt is high enough to form it now.
                ################################################################################

                # NEW_CONTINUOUS_ICE_SLAB_MASKS
                # Places where (1) it's on the ice sheet,
                #              (2) not in the long-term runoff zone, (1 & 2 both covered by STEADY_STATE_ACCUM_ZONE_MASK)
                #              (3) not where there were already ice slabs,
                #              (4) not where aquifers already exist instead of ice slabs, and
                #              (5) where melt has increased enough to form continuous ice slabs.
                NEW_CONTINUOUS_ICE_SLAB_MASK = STEADY_STATE_ACCUM_ZONE_MASK \
                                               & ~CONTINUOUS_ICE_SLAB_MASK \
                                               & ~AQUIFER_MASK \
                                               & (EM_10y_mean > median_min_continuous_EM)

                # For brand new ice slabs, set the porespace to zero, this has
                # seen enough melt where porespace is not accumulating on it.
                PORESPACE_ATOP_ICE_SLABS_START_0000[NEW_CONTINUOUS_ICE_SLAB_MASK] = 0.0
                PORESPACE_ATOP_ICE_SLABS_START_2000[NEW_CONTINUOUS_ICE_SLAB_MASK] = 2000.0

                ##################################################################
                # Branch 2, new aquifer areas, NOT where ice slabs already existed.
                ##################################################################

                # Make aquifer mask for this year... if yes, set this to an aquifer.
                ACCUM_MASK_10yr = SNOW_10y_mean >= mean_snow_95_cutoff
                # Mask only to areas that are in the long-term accumulation zone where firn exists.
                AQUIFER_POTENTIAL_MASK = STEADY_STATE_ACCUM_ZONE_MASK & ACCUM_MASK_10yr & ~CONTINUOUS_ICE_SLAB_MASK
                MELT_gt_200_MASK_10yr = (MELT_10y_mean + RAIN_10y_mean) >= 200.0
                AQUIFER_MASK = AQUIFER_POTENTIAL_MASK & MELT_gt_200_MASK_10yr

                ###################################################
                # BRANCH #3: If on a continuous ice slab area
                ###################################################

                # Subtract annual excess melt from available porespace --> new porespace
                # CHECK: Is this the right way to do this?
                # Melt will subtract both the pore space from near-surface snow, PLUS fill in pore-space beneath.
                # But if it's "excess melt" only, then near-surface snow is already filled, thus we're filling in old porespace.

                # NOTE: EM_single_year can be positive or negative.
                # If positive, it is the amount of melt beyond what filled in the annual snow.
                # If negative, it is the amount of potential porespace left behind after melt has done its job, therefore adding to the porespace already there.
                PORESPACE_ATOP_ICE_SLABS_START_0000 = PORESPACE_ATOP_ICE_SLABS_START_0000 - EM_single_year
                PORESPACE_ATOP_ICE_SLABS_START_0000[~STEADY_STATE_ACCUM_ZONE_MASK] = 0.0

                PORESPACE_ATOP_ICE_SLABS_START_2000 = PORESPACE_ATOP_ICE_SLABS_START_2000 - EM_single_year
                PORESPACE_ATOP_ICE_SLABS_START_2000[~STEADY_STATE_ACCUM_ZONE_MASK] = 0.0

                # For areas porespace <= 0, add to runoff, set porespace to zero
                # Multiply porespace by negative one, since runoff is areas with negative porespace availability in this forumulation.
                RUNOFF_ICE_SLAB_single_year_START_0000 = PORESPACE_ATOP_ICE_SLABS_START_0000 * -1
                RUNOFF_ICE_SLAB_single_year_START_2000 = PORESPACE_ATOP_ICE_SLABS_START_2000 * -1
                # Wherever excess melt has exceeded the available porespace, porespace will be negative.
                # Mask out all negative values of runoff (porespace still > 0) and set to zero.
                # ONLY areas considered in the previous continuous ice slab mask will be considered.
                RUNOFF_ICE_SLAB_single_year_START_0000[RUNOFF_ICE_SLAB_single_year_START_0000 < 0] = 0.0
                RUNOFF_ICE_SLAB_single_year_START_0000[~CONTINUOUS_ICE_SLAB_MASK] = 0.0

                RUNOFF_ICE_SLAB_single_year_START_2000[RUNOFF_ICE_SLAB_single_year_START_2000 < 0] = 0.0
                RUNOFF_ICE_SLAB_single_year_START_2000[~CONTINUOUS_ICE_SLAB_MASK] = 0.0

                # This is the sum of the runoff
                # Computed mmWE * km2 * (1e-6 km/mmWE) --> km3 of water, i.e. GT
                RUNOFF_ICE_SLAB_single_year_SUM_START_0000 = numpy.sum((RUNOFF_ICE_SLAB_single_year_START_0000 * CELLAREA * 1.0e-6).flatten())
                RUNOFF_ICE_SLAB_single_year_SUM_START_2000 = numpy.sum((RUNOFF_ICE_SLAB_single_year_START_2000 * CELLAREA * 1.0e-6).flatten())

                # For areas porespace > 0, no runoff, keep porespace. Set negative values to zero
                PORESPACE_ATOP_ICE_SLABS_START_0000[PORESPACE_ATOP_ICE_SLABS_START_0000 < 0] = 0.0
                PORESPACE_ATOP_ICE_SLABS_START_2000[PORESPACE_ATOP_ICE_SLABS_START_2000 < 0] = 0.0

                #################################################################
                # BRANCH #4: Not a previous aquifer, but high enough melt to form intermittent ice slabs
                #################################################################

                # Mask out intermittent ice slab areas
                # If accumulation too high, set it to an aquifer
                # If accumulation is low, set it to an Intermittent Ice Slab
                INTERMITTENT_ICE_SLAB_MASK = STEADY_STATE_ACCUM_ZONE_MASK \
                                             & ~CONTINUOUS_ICE_SLAB_MASK \
                                             & ~NEW_CONTINUOUS_ICE_SLAB_MASK \
                                             & ~AQUIFER_MASK \
                                             & (INTERMITTENT_ICE_SLAB_MASK | (EM_10y_mean > median_min_intermittent_EM))

                #################################################################
                # Combine single-year masks into steady masks.
                #################################################################
                CONTINUOUS_ICE_SLAB_MASK = CONTINUOUS_ICE_SLAB_MASK | NEW_CONTINUOUS_ICE_SLAB_MASK

                RUNOFF_CUMULATIVE_SUM_START_0000 = RUNOFF_CUMULATIVE_SUM_START_0000 + RUNOFF_ICE_SLAB_single_year_SUM_START_0000
                RUNOFF_CUMULATIVE_SUM_START_2000 = RUNOFF_CUMULATIVE_SUM_START_2000 + RUNOFF_ICE_SLAB_single_year_SUM_START_2000

                RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_0000 = RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_0000 + RUNOFF_ICE_SLAB_single_year_START_0000
                RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_2000 = RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_2000 + RUNOFF_ICE_SLAB_single_year_START_2000

                # Now compute what the RCM gives us for runoff over this area, in the SAME ACCUMULATION ZONE ONLY for apples-to-apples.
                # ALSO we're not computing/comparing runoff over aquifers here, so leave those out too.
                RUNOFF_RCM_COMPUTED_single_year[((~STEADY_STATE_ACCUM_ZONE_MASK) | AQUIFER_POTENTIAL_MASK)] = 0.0
                ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID = ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID + RUNOFF_RCM_COMPUTED_single_year
                RUNOFF_RCM_COMPUTED_single_year_SUM = numpy.sum(RUNOFF_RCM_COMPUTED_single_year * CELLAREA * 1e-6)
                ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_TOTAL = numpy.sum(ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID * CELLAREA * 1e-6)

                # Determine area with positive runoff each year (with memory) from the accumulation zone, as already computed by RCMs
                INTERIOR_RUNOFF_AREA_RCM_COMPUTED_KM2 = numpy.sum(CELLAREA[ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID > 0.0])

                # We're only concerned about porespace on top of ice slabs, not porespace accumulating elsewhere. Zero out porespace elsewhere.
                PORESPACE_ATOP_ICE_SLABS_START_0000[~(CONTINUOUS_ICE_SLAB_MASK | INTERMITTENT_ICE_SLAB_MASK)] = 0.0
                PORESPACE_ATOP_ICE_SLABS_START_2000[~(CONTINUOUS_ICE_SLAB_MASK | INTERMITTENT_ICE_SLAB_MASK)] = 0.0

                PORESPACE_ATOP_ICE_SLABS_ANNUAL_SUM_START_0000 = numpy.sum(PORESPACE_ATOP_ICE_SLABS_START_0000)
                PORESPACE_ATOP_ICE_SLABS_ANNUAL_SUM_START_2000 = numpy.sum(PORESPACE_ATOP_ICE_SLABS_START_2000)

                # Define an ice slab values mask, containing both the Continuous and Intermittent ice slab locations
                ICE_SLAB_VALUES = numpy.zeros(CONTINUOUS_ICE_SLAB_MASK.shape, dtype=numpy.int)
                ICE_SLAB_VALUES[CONTINUOUS_ICE_SLAB_MASK]   = 2 # 2 for continuous
                ICE_SLAB_VALUES[INTERMITTENT_ICE_SLAB_MASK] = 1 # 1 for intermittent


                #################################################################
                # Export annual data summary data lines
                #################################################################
                dataset_identifier_string = self.get_manager_name() + \
                            self._return_short_dataset_key(dataset_keys_list[i][-1])
                yearstart, yearend = yearspan

                datastring_line = ("{0},{1:4d},{2:4d},{3:d},{4:f},{5:d},{6:f},{7:d}," + \
                                   "{8:f},{9:d},{10:f},{11:f},{12:f},{13:f},{14:f}," + \
                                   "{15:f},{16:f},{17:f},{18:f},{19:f},{20:f},{21:f}," + \
                                   "{22:f},{23:f},{24:f},{25:f},{26:f},{27:f}\n").format(
                                  dataset_identifier_string,
                                  yearstart,
                                  yearend,
                                  numpy.count_nonzero(AQUIFER_POTENTIAL_MASK),
                                  numpy.sum(CELLAREA[AQUIFER_POTENTIAL_MASK]),
                                  numpy.count_nonzero(AQUIFER_MASK),
                                  numpy.sum(CELLAREA[AQUIFER_MASK]),
                                  numpy.count_nonzero(INTERMITTENT_ICE_SLAB_MASK),
                                  numpy.sum(CELLAREA[INTERMITTENT_ICE_SLAB_MASK]),
                                  numpy.count_nonzero(CONTINUOUS_ICE_SLAB_MASK),
                                  numpy.sum(CELLAREA[CONTINUOUS_ICE_SLAB_MASK]),
                                  STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2,
                                  PORESPACE_ATOP_ICE_SLABS_ANNUAL_SUM_START_0000,
                                  PORESPACE_ATOP_ICE_SLABS_ANNUAL_SUM_START_2000,
                                  RUNOFF_ICE_SLAB_single_year_SUM_START_0000,
                                  RUNOFF_ICE_SLAB_single_year_SUM_START_2000,
                                  RUNOFF_CUMULATIVE_SUM_START_0000,
                                  RUNOFF_CUMULATIVE_SUM_START_2000,
                                  INTERIOR_RUNOFF_AREA_RCM_COMPUTED_KM2,
                                  RUNOFF_RCM_COMPUTED_single_year_SUM,
                                  ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_TOTAL,
                                  RCM_RUNOFF_AREA_STEADYSTATE_KM2,
                                  RCM_RUNOFF_AREA_ANNUAL_KM2,
                                  RCM_RUNOFF_AREA_CUMULATIVE_KM2,
                                  RCM_RUNOFF_ANNUAL_FROM_STEADYSTATE_AREA_GT,
                                  RCM_RUNOFF_CUMULATIVE_FROM_STEADYSTATE_AREA_GT,
                                  RCM_RUNOFF_ANNUAL_FROM_INTERIOR_GT,
                                  RCM_RUNOFF_CUMULATIVE_FROM_INTERIOR_GT
                                 )

                f.write(datastring_line)

                TOTAL_DATA_STRING = TOTAL_DATA_STRING + datastring_line

                #################################################################
                # Export shapefiles of the relevant variables
                #################################################################
                # Only export a shapefile for every shapefile_every_N_years, and in 2013, or the very last year in the run.
                # and only if "skip_shapefiles" isn't set to true.

                if (not skip_shapefiles) and \
                   ((yearspan[1] % shapefile_every_N_years) == 0) or (yearspan[1] == 2013) or (yearspan == yearspans[-1]):
                    # Print the output to the screen for reference.
                    print
                    print datastring_line,

                    # Export shapefiles of:
                    # - firn aquifer potential mask
                    # - firn aquifer mask
                    # - ice slab mask with both continuous and intermittent in there
                    # - pore space atop ice slabs
                    # - runoff ANNUAL
                    # - runoff CUMULATIVE

                    # A helpful string to use for layer naming.
                    # Fill in {0} with "AQUIFER", "AQUIFER_POTENTIAL", "ICE_SLAB", "PORESPACE", etc...
                    layer_template_string = dataset_identifier_string + "_{0}_" + "{0}-{1}".format(yearstart, yearend)

                    layer_name_firn_aquifer_potential = layer_template_string.format("Aquifer_Potential")
                    layer_name_firn_aquifer           = layer_template_string.format("Aquifer")
                    layer_name_ice_slabs              = layer_template_string.format("Ice_Slabs")
                    layer_name_porespace_0            = layer_template_string.format("Porespace_START_0")
                    layer_name_porespace_2            = layer_template_string.format("Porespace_START_2")
                    layer_name_runoff_annual_0        = layer_template_string.format("Runoff_Annual_START_0")
                    layer_name_runoff_annual_2        = layer_template_string.format("Runoff_Annual_START_2")
                    layer_name_runoff_cumulative_0    = layer_template_string.format("Runoff_Cumulative_START_0")
                    layer_name_runoff_cumulative_2    = layer_template_string.format("Runoff_Cumulative_START_2")
                    layer_name_runoff_rcm_comp_annual = layer_template_string.format("Runoff_Annual_RCM_Computed")
                    layer_name_runoff_rcm_comp_cum    = layer_template_string.format("Runoff_Cumulative_RCM_Computed")


                    self.export_data_to_shapefile(variable=ICE_SLAB_VALUES,
                                                  variable_name="Ice_Slabs",
                                                  layer_name=layer_name_ice_slabs,
                                                  filename = layer_name_ice_slabs + ".shp",
                                                  variable_2 = EM_10y_mean,
                                                  variable_2_name = "ExcessMelt",
                                                  mask = (CONTINUOUS_ICE_SLAB_MASK | INTERMITTENT_ICE_SLAB_MASK))

                    self.export_data_to_shapefile(variable=AQUIFER_MASK,
                                                  variable_name="Aquifer",
                                                  layer_name=layer_name_firn_aquifer,
                                                  filename = layer_name_firn_aquifer + ".shp",
                                                  mask = AQUIFER_MASK)

                    self.export_data_to_shapefile(variable=AQUIFER_POTENTIAL_MASK,
                                                  variable_name="Aq_Poten",
                                                  layer_name=layer_name_firn_aquifer_potential,
                                                  filename = layer_name_firn_aquifer_potential + ".shp",
                                                  mask = AQUIFER_POTENTIAL_MASK)

                    self.export_data_to_shapefile(variable=PORESPACE_ATOP_ICE_SLABS_START_0000,
                                                  variable_name="PoreSpace",
                                                  layer_name=layer_name_porespace_0,
                                                  filename = layer_name_porespace_0 + ".shp",
                                                  mask = (PORESPACE_ATOP_ICE_SLABS_START_0000 > 0.0))

                    self.export_data_to_shapefile(variable=PORESPACE_ATOP_ICE_SLABS_START_2000,
                                                  variable_name="PoreSpace",
                                                  layer_name=layer_name_porespace_2,
                                                  filename = layer_name_porespace_2 + ".shp",
                                                  mask = (PORESPACE_ATOP_ICE_SLABS_START_2000 > 0.0))

                    self.export_data_to_shapefile(variable=RUNOFF_ICE_SLAB_single_year_START_0000,
                                                  variable_name="RunoffAnn",
                                                  layer_name=layer_name_runoff_annual_0,
                                                  filename = layer_name_runoff_annual_0 + ".shp",
                                                  mask = (RUNOFF_ICE_SLAB_single_year_START_0000 > 0.0))

                    self.export_data_to_shapefile(variable=RUNOFF_ICE_SLAB_single_year_START_2000,
                                                  variable_name="RunoffAnn",
                                                  layer_name=layer_name_runoff_annual_2,
                                                  filename = layer_name_runoff_annual_2 + ".shp",
                                                  mask = (RUNOFF_ICE_SLAB_single_year_START_2000 > 0.0))

                    self.export_data_to_shapefile(variable=RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_0000,
                                                  variable_name="RunoffCum",
                                                  layer_name=layer_name_runoff_cumulative_0,
                                                  filename = layer_name_runoff_cumulative_0 + ".shp",
                                                  mask = (RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_0000 > 0.0))

                    self.export_data_to_shapefile(variable=RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_2000,
                                                  variable_name="RunoffCum",
                                                  layer_name=layer_name_runoff_cumulative_2,
                                                  filename = layer_name_runoff_cumulative_2 + ".shp",
                                                  mask = (RUNOFF_ICE_SLAB_CUMULATIVE_GRID_START_2000 > 0.0))

                    self.export_data_to_shapefile(variable=RUNOFF_RCM_COMPUTED_single_year,
                                                  variable_name="RunoffAnn",
                                                  layer_name=layer_name_runoff_rcm_comp_annual,
                                                  filename = layer_name_runoff_rcm_comp_annual + ".shp",
                                                  mask = (RUNOFF_RCM_COMPUTED_single_year > 0.0))

                    self.export_data_to_shapefile(variable=ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID,
                                                  variable_name="RunoffCum",
                                                  layer_name=layer_name_runoff_rcm_comp_cum,
                                                  filename = layer_name_runoff_rcm_comp_cum + ".shp",
                                                  mask = (ICE_SLAB_RUNOFF_RCM_COMPUTED_SUM_GRID > 0.0))

        #################################################################
        # Export data to CSV file, close if necessary.
        #################################################################
        if type(csv_output_file) == str:
            f.close()
            print
            print os.path.split(csv_output_file)[-1], "written."


    def find_median_lower_upper_intermittent_limits_of_dataset(self, data):
        '''Given a gridded dataset "data", find the median lower, upper and "intermittent" limits of the validation
        IceBridge tracks gethered from ::_find_ice_thickness_minmax_indices_along_validation_tracks().
        Return the median of the values at each threshold.'''
        # FIND THE MIN-MAX THICKNESS OF WHERE THE OBSERVATIONS ARE: COMPUTE THRESHOLDS FROM THAT AND MASK OUT TO GET EXTENT.
        track_lower_rows,  \
        track_lower_cols,  \
        track_upper_rows,  \
        track_upper_cols,  \
        track_outlier_rows,\
        track_outlier_cols = self._find_ice_thickness_minmax_indices_along_validation_tracks(export=False)

        track_outlier_rows = track_outlier_rows[track_outlier_rows != -9999]
        track_outlier_cols = track_outlier_cols[track_outlier_cols != -9999]
        assert len(track_outlier_rows) == len(track_outlier_cols)

        values_lower = data[track_lower_rows, track_lower_cols]
        values_upper = data[track_upper_rows, track_upper_cols]
        values_intermittent = data[track_outlier_rows, track_outlier_cols]

        return numpy.median(values_lower), numpy.median(values_upper), numpy.median(values_intermittent)

    def _does_regional_pixel_exceed_regional_median_baseline_area_EMI(self, prow, pcol, EMI, baseline_mask):
        '''Helper function for ::export_EMI_thresholded_maps() to help prevent regional "overreach" of which pixels have ice lenses.
        Ice lenses should not suddenly appear with a HIGHER EMI (usually lower elevation) than pixels around it that were in the original mask.
        Occasional false negatives are okay since we'll include all the baseline_mask pixels in our analysis anyway.'''
        # NOTE: DEFUNCT FUNCTION SINCE USING RETHRESHOLDING STRATEGY.  DON'T NEED IT.
        pixel_EMI = EMI[prow,pcol]

        bbox_row_min = max(0,prow-50)
        bbox_row_max = min(prow+50,EMI.shape[0]-1)
        bbox_col_min = max(0,pcol-50)
        bbox_col_max = min(pcol+50,EMI.shape[0]-1)
        bbox_mask = numpy.zeros(EMI.shape, dtype=numpy.bool)
        bbox_mask[bbox_row_min:bbox_row_max+1, bbox_col_min:bbox_col_max+1] = True
        # Are there any EMI pixels nearby?
        regional_EMI_mask = bbox_mask & baseline_mask
        if numpy.count_nonzero(regional_EMI_mask) == 0:
            return False

        regional_baseline_EMIs = EMI[regional_EMI_mask]
        return pixel_EMI > numpy.median(regional_baseline_EMIs)


    def find_median_EMI_lower_upper_intermittent_and_snow_cutoffs(self, export=False, export_image=False):
        '''Looking at the median values for each Reanalysis dataset, find the mean of the median EMI values for the
        lower, upper, and upper-intermittent limits where ice slabs occur in the data.  Return these three values.'''

        cutoffs_picklefile_name = self._get_EMI_cutoff_value_results_picklefile_name()
        reanalysis_dataset_key_list = (('2014ERA',),('ERA-int_1979-2014_10km',),('NCEPv1_1948-2015_20km',),('RACMO2.3_2015',))

        EMI_cutoffs_dict = None

        if (not export) and os.path.exists(cutoffs_picklefile_name):
            print "Reading", os.path.split(cutoffs_picklefile_name)[-1]
            f = open(cutoffs_picklefile_name, 'r')
            EMI_cutoffs_dict, snow_cutoffs_dicts_global = pickle.load(f)
            f.close()

            # These are the four reanalysis datasets we're using.
            lower_medians = [EMI_cutoffs_dict[dkey]["Median Lower Continuous"] for dkey in reanalysis_dataset_key_list]
            upper_medians = [EMI_cutoffs_dict[dkey]["Median Upper Continuous"] for dkey in reanalysis_dataset_key_list]
            intermittent_medians = [EMI_cutoffs_dict[dkey]["Median Intermittent Continuous"] for dkey in reanalysis_dataset_key_list]
            snow_cutoffs_95 = [snow_cutoffs_dicts_global[dkey]["95%"] for dkey in reanalysis_dataset_key_list]

            lower_mean = numpy.mean(lower_medians)
            upper_mean = numpy.mean(upper_medians)
            intermittent_mean = numpy.mean(intermittent_medians)
            snow_cutoffs_95_mean = numpy.mean(snow_cutoffs_95)

            if not export_image:
                return lower_mean, upper_mean, intermittent_mean, snow_cutoffs_95_mean

        if EMI_cutoffs_dict is None:
            # Let's save these values in case we want to plot them with the ranges needed.
            # Keys are the dataset names, values are dictionaries with
            # ["Lower Continuous", "Upper Continuous", "Upper Intermittent"] as keys and
            # the whole population of cutoff values as values.
            EMI_cutoffs_dict = {}
            # Change the snow_cutoffs_dicts from being a manager-specific item to a global item
            snow_cutoffs_dicts_global = {}

            MAR_m = MAR_Manager()
            HIR_m = HIRHAM5_Manager()
            RAC_m = RACMO_Manager()

            for manager in (HIR_m, MAR_m, RAC_m):
                # Load the dataset for the time span needed. In this case we're ignoring the OLS models and just using the snow dictionaries.
                dataset_keys_list, \
                yearspans_dicts_list,    \
                baseline_yearspans_list, \
                ols_models_list, \
                ols_models_EMI_only_list, \
                snow_cutoffs_dicts = manager.compute_ice_thickness_parameterizations(export=False, verbose=False)

                # FIND THE MIN-MAX THICKNESS OF WHERE THE OBSERVATIONS ARE: COMPUTE THRESHOLDS FROM THAT AND MASK OUT TO GET EXTENT.
                track_lower_rows,  \
                track_lower_cols,  \
                track_upper_rows,  \
                track_upper_cols,  \
                track_outlier_rows,\
                track_outlier_cols = manager._find_ice_thickness_minmax_indices_along_validation_tracks(export=False)

                # Some of the tracks have no upper outlier points.  Get rid of outlier that have empty indices (-9999)
                track_outlier_rows = track_outlier_rows[track_outlier_rows != -9999]
                track_outlier_cols = track_outlier_cols[track_outlier_cols != -9999]

                for i,dataset_key in enumerate(dataset_keys_list):
                    yearspans = yearspans_dicts_list[i][2014]
                    print "~~~", dataset_key, yearspans, "~~~"
                    EMI = manager.compute_EMI_values(dataset_key, yearspan=yearspans,include_rain=True,compute_mean=True)

                    lower_EMIs = EMI[track_lower_rows, track_lower_cols]
                    upper_EMIs = EMI[track_upper_rows, track_upper_cols]
                    outlier_EMIs = EMI[track_outlier_rows, track_outlier_cols]
                    EMI_cutoffs_dict[dataset_key] = {}
                    EMI_cutoffs_dict[dataset_key]["Lower Continuous"] = lower_EMIs
                    EMI_cutoffs_dict[dataset_key]["Upper Continuous"] = upper_EMIs
                    EMI_cutoffs_dict[dataset_key]["Upper Intermittent"] = outlier_EMIs

                    mean_lower_EMIs = numpy.mean(lower_EMIs)
                    mean_upper_EMIs = numpy.mean(upper_EMIs)
                    mean_outlier_EMIs = numpy.mean(outlier_EMIs)
                    median_lower_EMIs = numpy.median(lower_EMIs)
                    median_upper_EMIs = numpy.median(upper_EMIs)
                    median_outlier_EMIs = numpy.median(outlier_EMIs)

                    EMI_cutoffs_dict[dataset_key]["Median Lower Continuous"] = median_lower_EMIs
                    EMI_cutoffs_dict[dataset_key]["Median Upper Continuous"] = median_upper_EMIs
                    EMI_cutoffs_dict[dataset_key]["Median Intermittent Continuous"] = median_outlier_EMIs

                    snow_cutoffs_dicts_global[dataset_key] = snow_cutoffs_dicts[i]

                    print "Lower (min,mean,std,median,max): {0:0.3f}, {1:0.3f}, {2:0.3f}, {3:0.3f}, {4:0.3f}".format(numpy.min(lower_EMIs),
                                                                                            mean_lower_EMIs,
                                                                                            numpy.std(lower_EMIs),
                                                                                            median_lower_EMIs,
                                                                                            numpy.max(lower_EMIs))
                    print "Upper (min,mean,std,median,max): {0:0.3f}, {1:0.3f}, {2:0.3f}, {3:0.3f}, {4:0.3f}".format(numpy.min(upper_EMIs),
                                                                                            mean_upper_EMIs,
                                                                                            numpy.std(upper_EMIs),
                                                                                            median_upper_EMIs,
                                                                                            numpy.max(upper_EMIs))
                    print "Outlier (min,mean,std,median,max): {0:0.3f}, {1:0.3f}, {2:0.3f}, {3:0.3f}, {4:0.3f}".format(numpy.min(outlier_EMIs),
                                                                                              mean_outlier_EMIs,
                                                                                              numpy.std(outlier_EMIs),
                                                                                              median_outlier_EMIs,
                                                                                              numpy.max(outlier_EMIs))

        if export:
            # Let's save the cutoff populations to a picklefile
            f = open(cutoffs_picklefile_name, 'w')
            pickle.dump([EMI_cutoffs_dict, snow_cutoffs_dicts_global],f)
            f.close()
            print os.path.split(cutoffs_picklefile_name)[-1], "written."

        # Pull out the populations, export box-and-whisker plots for all of these.
        if export_image:
            lower_continuous_EMIs_lists = [EMI_cutoffs_dict[dkey]["Lower Continuous"] for dkey in reanalysis_dataset_key_list]
            upper_continuous_EMIs_lists = [EMI_cutoffs_dict[dkey]["Upper Continuous"] for dkey in reanalysis_dataset_key_list]
            intermittent_EMIs_lists = [EMI_cutoffs_dict[dkey]["Upper Intermittent"] for dkey in reanalysis_dataset_key_list]

            lower_medians = [EMI_cutoffs_dict[dkey]["Median Lower Continuous"] for dkey in reanalysis_dataset_key_list]
            upper_medians = [EMI_cutoffs_dict[dkey]["Median Upper Continuous"] for dkey in reanalysis_dataset_key_list]
            intermittent_medians = [EMI_cutoffs_dict[dkey]["Median Intermittent Continuous"] for dkey in reanalysis_dataset_key_list]

            dkeys = [dkey[0] for dkey in reanalysis_dataset_key_list]
            labels_dict = {"2014ERA":"HIRHAM5, ERA-Interim",
                           "NCEPv1_1948-2015_20km":"MAR, NCEPv1",
                           "ERA-int_1979-2014_10km":"MAR, ERA-Interim",
                           "RACMO2.3_2015":"RACMO 2.3, ERA-Interim"}
            dlabels = [labels_dict[dset] for dset in dkeys]

            fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True)
            positions = range(len(reanalysis_dataset_key_list)-1,-1,-1)

            ax1.boxplot(lower_continuous_EMIs_lists, vert=False, showfliers=True, positions=positions, labels=dlabels)
            ax1.axvline(x=numpy.mean(lower_medians), color="blue",linestyle="--")
            ax1.set_title("Continuous\nLower\nLimit")
            ax2.boxplot(upper_continuous_EMIs_lists, vert=False, showfliers=True, positions=positions, labels=dlabels)
            ax2.axvline(x=numpy.mean(upper_medians), color="blue",linestyle="--")
            ax2.set_title("Continuous\nUpper\nLimit")
            ax3.boxplot(intermittent_EMIs_lists, vert=False, showfliers=True, positions=positions, labels=dlabels)
            ax3.axvline(x=numpy.mean(intermittent_medians), color="blue",linestyle="--")
            ax3.set_title("Intermittent\nUpper\nLimit")

            ax2.set_xlabel("Excess Melt (mm w.e. yr$^{-1}$)")

            fig.tight_layout()
            figname = os.path.join(STATISTICS_OUTPUT_FOLDER, "EMI_cutoff_boxplots.png")
            fig.savefig(figname, dpi=600)
            print "Exported", os.path.split(figname)[1]
            plt.close()

        # These are the four reanalysis datasets we're using.
        lower_medians = [EMI_cutoffs_dict[dkey]["Median Lower Continuous"] for dkey in reanalysis_dataset_key_list]
        upper_medians = [EMI_cutoffs_dict[dkey]["Median Upper Continuous"] for dkey in reanalysis_dataset_key_list]
        intermittent_medians = [EMI_cutoffs_dict[dkey]["Median Intermittent Continuous"] for dkey in reanalysis_dataset_key_list]
        snow_cutoffs_95 = [snow_cutoffs_dicts_global[dkey]["95%"] for dkey in reanalysis_dataset_key_list]
        print snow_cutoffs_95

        lower_mean = numpy.mean(lower_medians)
        upper_mean = numpy.mean(upper_medians)
        intermittent_mean = numpy.mean(intermittent_medians)
        snow_cutoffs_95_mean = numpy.mean(snow_cutoffs_95)
        return lower_mean, upper_mean, intermittent_mean, snow_cutoffs_95_mean

    def _read_EMI_ICE_SLAB_decadal_outputs_v1(self):
        print "Reading", os.path.split(EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE)[-1]
        csv_file = open(EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE, 'r')
        lines = [line.strip() for line in csv_file.readlines() if len(line.strip()) > 0]
        csv_file.close()

        header = dict([(item,i) for i,item in enumerate(lines[0].split(','))])
        lines = lines[1:]

        max_key_length = numpy.max([len(line.split(',')[0]) for line in lines])
        dt = numpy.dtype([('DATASET'    ,numpy.string_, max_key_length),
                          ('YEAR_START' ,numpy.int, 1),
                          ('YEAR_END'   ,numpy.int, 1),
                          ('NUM_PIXELS_NEWLAYERS' ,numpy.int, 1),
                          ('AREA_NEWLAYERS' ,numpy.float64, 1),
                          ('NUM_PIXELS_TOTAL' ,numpy.int, 1),
                          ('TOTAL_AREA' ,numpy.float64, 1)])

        # Filter out other header lines scattered through the rest of them.
        lines = [line for line in lines if line[0:7] != "DATASET"]

        data_array = numpy.empty((len(lines),), dtype=dt)
        for i,line in enumerate(lines):
            items = line.split(",")
            data_array[i]['DATASET']                  = items[header['DATASET']]
            data_array[i]['YEAR_START']               = int(items[header['YEAR_START']])
            data_array[i]['YEAR_END']                 = int(items[header['YEAR_END']])
            data_array[i]['NUM_PIXELS_NEWLAYERS']     = int(items[header['NUM_PIXELS_NEWLAYERS']])
            data_array[i]['AREA_NEWLAYERS']           = float(items[header['AREA_NEWLAYERS']])
            data_array[i]['NUM_PIXELS_TOTAL']         = int(items[header['NUM_PIXELS_TOTAL']])
            data_array[i]['TOTAL_AREA']               = float(items[header['TOTAL_AREA']])

        return data_array

    def _read_EMI_ICE_SLAB_decadal_outputs_v3(self):
        print "Reading", os.path.split(EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE)[-1]
        csv_file = open(EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE, 'r')
        lines = [line.strip() for line in csv_file.readlines() if len(line.strip()) > 0]
        csv_file.close()

        header = dict([(item,i) for i,item in enumerate(lines[0].split(','))])
        lines = lines[1:]

        max_key_length = numpy.max([len(line.split(',')[0]) for line in lines])
        dt = numpy.dtype([('DATASET'    ,numpy.string_, max_key_length),
                          ('YEAR_START' ,numpy.int, 1),
                          ('YEAR_END'   ,numpy.int, 1),
                          ('FIRN_AQUIFER_POTENTIAL_PIXELCOUNT',numpy.int, 1),
                          ('FIRN_AQUIFER_POTENTIAL_AREA_KM2',numpy.float, 1),
                          ('FIRN_AQUIFER_PIXELCOUNT',numpy.int, 1),
                          ('FIRN_AQUIFER_AREA_KM2',numpy.float, 1),
                          ('ICE_SLAB_INTERMITTENT_PIXELCOUNT',numpy.int, 1),
                          ('ICE_SLAB_INTERMITTENT_AREA_KM2',numpy.float, 1),
                          ('ICE_SLAB_CONTINUOUS_PIXELCOUNT',numpy.int, 1),
                          ('ICE_SLAB_CONTINUOUS_AREA_KM2',numpy.float, 1),
                          ('STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2',numpy.float, 1),
                          ('PORESPACE_ATOP_ICE_SLABS_GT_START_0000',numpy.int, 1),
                          ('PORESPACE_ATOP_ICE_SLABS_GT_START_2000',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_ANNUAL_GT_START_2000',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_TOTAL_GT_START_0000',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_TOTAL_GT_START_2000',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_ANNUAL_RCM_COMPUTED_GT',numpy.float, 1),
                          ('ICE_SLAB_RUNOFF_TOTAL_RCM_COMPUTED_GT',numpy.float, 1),
                        ])
        # Filter out other header lines scattered through the rest of them.
        lines = [line for line in lines if line[0:7] != "DATASET"]

        data_array = numpy.empty((len(lines),), dtype=dt)
        for i,line in enumerate(lines):
            items = line.split(",")
            data_array[i]['DATASET']                  = items[header['DATASET']]
            data_array[i]['YEAR_START']               = int(items[header['YEAR_START']])
            data_array[i]['YEAR_END']                 = int(items[header['YEAR_END']])
            data_array[i]['FIRN_AQUIFER_POTENTIAL_PIXELCOUNT']     = int  (items[header['FIRN_AQUIFER_POTENTIAL_PIXELCOUNT']])
            data_array[i]['FIRN_AQUIFER_POTENTIAL_AREA_KM2']       = float(items[header['FIRN_AQUIFER_POTENTIAL_AREA_KM2']])
            data_array[i]['FIRN_AQUIFER_PIXELCOUNT']               = int  (items[header['FIRN_AQUIFER_PIXELCOUNT']])
            data_array[i]['FIRN_AQUIFER_AREA_KM2']                 = float(items[header['FIRN_AQUIFER_AREA_KM2']])
            data_array[i]['ICE_SLAB_INTERMITTENT_PIXELCOUNT']      = int  (items[header['ICE_SLAB_INTERMITTENT_PIXELCOUNT']])
            data_array[i]['ICE_SLAB_INTERMITTENT_AREA_KM2']        = float(items[header['ICE_SLAB_INTERMITTENT_AREA_KM2']])
            data_array[i]['ICE_SLAB_CONTINUOUS_PIXELCOUNT']        = int  (items[header['ICE_SLAB_CONTINUOUS_PIXELCOUNT']])
            data_array[i]['ICE_SLAB_CONTINUOUS_AREA_KM2']          = float(items[header['ICE_SLAB_CONTINUOUS_AREA_KM2']])
            data_array[i]['STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2'] = float(items[header['STEADY_STATE_RUNOFF_AREA_BELOW_ICE_SLABS_KM2']])
            data_array[i]['PORESPACE_ATOP_ICE_SLABS_GT_START_0000']= float(items[header['PORESPACE_ATOP_ICE_SLABS_GT_START_0000']])
            data_array[i]['PORESPACE_ATOP_ICE_SLABS_GT_START_2000']= float(items[header['PORESPACE_ATOP_ICE_SLABS_GT_START_2000']])
            data_array[i]['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000']  = float(items[header['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000']])
            data_array[i]['ICE_SLAB_RUNOFF_ANNUAL_GT_START_2000']  = float(items[header['ICE_SLAB_RUNOFF_ANNUAL_GT_START_2000']])
            data_array[i]['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000']   = float(items[header['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000']])
            data_array[i]['ICE_SLAB_RUNOFF_TOTAL_GT_START_2000']   = float(items[header['ICE_SLAB_RUNOFF_TOTAL_GT_START_2000']])
            data_array[i]['ICE_SLAB_RUNOFF_ANNUAL_RCM_COMPUTED_GT']= float(items[header['ICE_SLAB_RUNOFF_ANNUAL_RCM_COMPUTED_GT']])
            data_array[i]['ICE_SLAB_RUNOFF_TOTAL_RCM_COMPUTED_GT'] = float(items[header['ICE_SLAB_RUNOFF_TOTAL_RCM_COMPUTED_GT']])

        return data_array


########################################################################################################
#########################################  HIRHAM-5  ###################################################
########################################################################################################
class HIRHAM5_Manager(RCM_Manager):

    _LIST_dataset_keys = ["2010GCM", "2014ERA", "2050GCM45", "2050GCM85", "2100GCM45", "2100GCM85"]
    _LIST_variable_keys = ["SMB", "RUNOFF", "MELT", "TIME", "EVAP_SUBL", "RAIN", "SNOW", "MASK", "LAT", "LON", "CELLAREA", "TEMP", "BASINS", "ELEV", "TEMPMAGNITUDE"]

    def __init__(self):
        # Initialize the parent class
        super(HIRHAM5_Manager, self).__init__()
        # Get filler dictionaires
        fnames_dict, dsets_dict, vars_dict = self._build_empty_filename_and_dataset_dictionaries()
        # Assign them.
        self._DICT_filenames = fnames_dict # Maintains the filename for each dataset and variable in that dataset_key
        self._DICT_datasets = dsets_dict   # Maintains the NetCDF4 dataset object for each variable in that dataset_key
        self._DICT_variables = vars_dict  # Maintains the actual values for each variable in that dataset_key
        self._LAST_DATASET_KEY_USED = None

    def _return_export_shapefile_folder(self):
        return HIRHAM_EXPORT_SHAPEFILE_FOLDER

    def _return_short_dataset_key(self, dataset_key):
        return dataset_key

    def _return_datatype_base(self):
        '''Useful for creating filenames with the base moniker in them, used by the parent class.'''
        return "HIRHAM5"

    def _split_dataset_into_decade_yearspans(self, dataset_key):
        '''Take each dataset, and return a list of year-spans that correspond to approximate "decades" within that dataset.'''
        assert dataset_key in self._LIST_dataset_keys
        if dataset_key == "2010GCM":
            return [(1991,2000),
                    (2001,2010)]
        elif dataset_key == "2014ERA":
            return [(1980,1994),
                    (1995,2004),
                    (2005,2014)]
        elif dataset_key in ("2050GCM45", "2050GCM85"):
            return [(2031,2040),
                    (2041,2050)]
        elif dataset_key in ("2100GCM45", "2100GCM85"):
            return [(2081,2090),
                    (2091,2100)]
        else:
            raise ValueError("Unknown dataset key '{0}' in HIRHAM5_Manager.".format(dataset_key))

    def _build_empty_filename_and_dataset_dictionaries(self):
        '''Build a dictionary from which to get datasets on demand.
        These won't all be immediately opened, but will be opened as needed as they're requested.
        STRUCTURE: self.hirham_datasets {
            ["2010GCM", "2014ERA", "2050GCM45", "2050GCM85", "2100GCM45", "2100GCM85"]
            each of which is a dict with keys of the variables.'''
        _hirham_filenames_dict = {}
        _hirham_datasets_dict = {}
        _hirham_variables_dict = {}

        for dset in self._LIST_dataset_keys:
            _hirham_filenames_dict[dset] = {}
            _hirham_datasets_dict[dset]  = {}
            _hirham_variables_dict[dset] = {}

            if dset == "2010GCM":
                fname_means  = HIRHAM_FILENAMES_2010_GCM["YEARLY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2010_GCM["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2010_GCM["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2010_GCM["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2010_GCM["YEARLY_TEMP"]

            elif dset == "2014ERA":
                fname_means  = HIRHAM_FILENAMES_2014_ERA["DAILY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2014_ERA["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2014_ERA["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2014_ERA["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2014_ERA["YEARLY_TEMP"]

            elif dset == "2050GCM45":
                fname_means  = HIRHAM_FILENAMES_2050_GCM_RCP45["YEARLY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2050_GCM_RCP45["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2050_GCM_RCP45["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2050_GCM_RCP45["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2050_GCM_RCP45["YEARLY_TEMP"]

            elif dset == "2050GCM85":
                fname_means  = HIRHAM_FILENAMES_2050_GCM_RCP85["YEARLY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2050_GCM_RCP85["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2050_GCM_RCP85["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2050_GCM_RCP85["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2050_GCM_RCP85["YEARLY_TEMP"]

            elif dset == "2100GCM45":
                fname_means  = HIRHAM_FILENAMES_2100_GCM_RCP45["YEARLY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2100_GCM_RCP45["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2100_GCM_RCP45["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2100_GCM_RCP45["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2100_GCM_RCP45["YEARLY_TEMP"]

            elif dset == "2100GCM85":
                fname_means  = HIRHAM_FILENAMES_2100_GCM_RCP85["YEARLY_MEANS"]
                fname_smb    = HIRHAM_FILENAMES_2100_GCM_RCP85["YEARLY_SMB"]
                fname_runoff = HIRHAM_FILENAMES_2100_GCM_RCP85["YEARLY_RUNOFF"]
                fname_melt   = HIRHAM_FILENAMES_2100_GCM_RCP85["YEARLY_MELT"]
                fname_temp   = HIRHAM_FILENAMES_2100_GCM_RCP85["YEARLY_TEMP"]

            _hirham_filenames_dict[dset]["TIME"  ] = fname_means
            _hirham_filenames_dict[dset]["EVAP_SUBL"] = fname_means
            _hirham_filenames_dict[dset]["RAIN"  ] = fname_means
            _hirham_filenames_dict[dset]["SNOW"  ] = fname_means
            _hirham_filenames_dict[dset]["SMB"   ] = fname_smb
            _hirham_filenames_dict[dset]["RUNOFF"] = fname_runoff
            _hirham_filenames_dict[dset]["MELT"  ] = fname_melt
            _hirham_filenames_dict[dset]["TEMP"  ] = fname_temp
            _hirham_filenames_dict[dset]["LAT"   ] = HIRHAM_GRID_DATAFILE
            _hirham_filenames_dict[dset]["LON"   ] = HIRHAM_GRID_DATAFILE
            _hirham_filenames_dict[dset]["CELLAREA"] = HIRHAM_GRID_DATAFILE
            _hirham_filenames_dict[dset]["BASINS"] = HIRHAM_BASIN_DATAFILE
            _hirham_filenames_dict[dset]["MASK"  ] = HIRHAM_MASK_DATAFILE

            # Build an empty dictionary, for performance's sake.  Fill in the None's with datasets as they're called.
            for key in self._LIST_variable_keys:
                _hirham_datasets_dict[dset][key] = None
                _hirham_variables_dict[dset][key] = None

        return _hirham_filenames_dict, _hirham_datasets_dict, _hirham_variables_dict


    def _get_variable_SUBCLASS(self, var_name, dataset_key=None, yearspan=None, annual=None):
        '''The Child function for RCM_Manager.get_variable().  This implements the HIRHAM version of it.
        Returns the dataset needed here.

        Yearspans are inclusive on both ends, to return 1981-2010 (2010 included), put a (1981,2010) tuple.'''

        if dataset_key is None:
            dataset_key = self._LAST_DATASET_KEY_USED if (self._LAST_DATASET_KEY_USED != None) else "2014ERA"

        self._LAST_DATASET_KEY_USED = dataset_key

        # Convert all keys to upper-case.
        var_name = var_name.upper()
        # Make sure this variable is in our options.
        assert var_name in self._LIST_variable_keys

        ################################################################################
        # STEP 1: If we haven't opened the dataset yet, do it and fill in the rest of _DICT_datasets for that dataset.
        #       The ELEV datasets are different (don't have a netCDF4 dataset)
        if (self._DICT_datasets[dataset_key][var_name] == None) and (var_name not in ("ELEV")):

            try:
                dataset = netCDF4.Dataset(self._DICT_filenames[dataset_key][var_name], mode="r")
            except RuntimeError:
                print "'{0}'[{1}]".format(dataset_key, var_name), "in file", self._DICT_filenames[dataset_key][var_name], "not found."

            if var_name in ("LAT","LON", "CELLAREA","MASK","BASINS"):
                if var_name in ("LAT","LON", "CELLAREA"):
                    var_list = ("LAT","LON", "CELLAREA")
                elif var_name == "MASK":
                    var_list = ("MASK",)
                elif var_name == "BASINS":
                    var_list = ("BASINS",)

                # Put this dataset in all the dataset entries it shares (in all dataset_keys too!), so we needn't open the dataset again.
                for dset_key in self._LIST_dataset_keys:
                    for _vname in var_list:
                        self._DICT_datasets[dset_key][_vname] = dataset

            elif var_name in ("TIME","EVAP_SUBL", "RAIN", "SNOW", "SMB", "RUNOFF", "MELT", "TEMP"):
                if var_name in ("TIME","EVAP_SUBL", "RAIN", "SNOW"):
                    var_list = ("TIME","EVAP_SUBL", "RAIN", "SNOW")
                elif var_name == "SMB":
                    var_list = ("SMB",)
                elif var_name == "RUNOFF":
                    var_list = ("RUNOFF",)
                elif var_name == "MELT":
                    var_list = ("MELT",)
                elif var_name == "TEMP":
                    var_list = ("TEMP",)
                # Put this dataset in all the dataset entries it shares, so we needn't open the dataset again.
                for _vname in var_list:
                    self._DICT_datasets[dataset_key][_vname] = dataset

        ################################################################################
        # STEP 2: If we haven't opened that variable yet, do it and fill in the variable.
        if self._DICT_variables[dataset_key][var_name] is None:

            # This is the key for the actual variable names in the NetCDF4 file.
            var_key = {"TIME"     :"time"    ,
                       "EVAP_SUBL":"evspsbl" ,
                       "RAIN"     :"rainfall",
                       "SNOW"     :"snfall"  ,
                       "RUNOFF"   :"mrros"   ,
                       "SMB"      :"smb"     ,
                       "MELT"     :"snmelt"  ,
                       "TEMP"     :"tas"     ,
                       "MASK"     :"glacGRL" ,
                       "BASINS"   :"maskall" ,
                       "LAT"      :"lat"     ,
                       "LON"      :"lon"     ,
                       "CELLAREA" :"cellarea",
                       "ELEV"     : None     ,
                       "TEMPMAGNITUDE": None ,
                       } [var_name]

            if var_name == "ELEV":
                # Fill in elevations (read from picklefile, or interpolated from RACMO)
                variable = self._get_elevation_data()
            elif var_name == "RAIN" and dataset_key == "2050GCM85":
                # DUE TO A BUG in the PRECIP variables, the 2050 RCP8.5 dataset gives really wrong (negative) rainfall values.
                # Just use the ones from the 2050 RCP4.5 dataset for now.  We will get the fixed data later.
                variable = self._DICT_datasets["2050GCM45"][var_name].variables[var_key]
            else:
                variable = self._DICT_datasets[dataset_key][var_name].variables[var_key]

            # If it's time, convert "days since XXXX" to integers
            if var_name == "TIME":
                variable = self._return_years_integers(variable)
            else:
                # Get the array out from the variable object
                variable = variable[:]

            # Individual post-processing of variable data before returning.
            # Some variables come as a "masked" variable. Just give us the basic data.
            if numpy.ma.is_masked(variable):
                variable = variable.data

            # Convert the mask to a boolean (rather than integer) format.
            if var_name == "MASK":
                variable = numpy.array(variable, dtype=numpy.bool)

            # For some weird reason, these three variables in the 2014ERA dataset are given as (kg m-2 day-1), rather than (kg m-2 year-1).
            # The numbers verify this.  Multiply by 1-year of days to get yearly totals.
            if (dataset_key == "2014ERA") and (var_name in ("EVAP_SUBL", "RAIN", "SNOW")):
                variable = variable * 365.25

            # Convert Kelvin to Celcius
            if var_name == "TEMP":
                variable = variable - 273.15

            # Convert cell area from m2 to km2
            if var_name == "CELLAREA":
                variable = variable / 1e6

            # Basins are a bit weird.  Get the basin integers here.
            if var_name == "BASINS":
                # The fetched variable is the basin masks, (19,602,402)
                maskall = variable
                # Retrieve the basin numbers, (19,)
                basins = self._DICT_datasets[dataset_key][var_name].variables["basins"][:]
                # Create a zeroed-out array
                basin_integers = numpy.zeros(maskall.shape[1:], dtype=numpy.integer)

                # Fill it with each of the basin numbers
                for i in range(basins.size):
                    basin_integers[maskall[i,:,:] == 1.0] = numpy.round(basins[i])

                # This is our new variable: doesn't yet include Flade Isblake
                variable = basin_integers

            # Get rid of empty 1-d dimensions in this data.
            variable.shape = [d for d in variable.shape if d>1]

            # Save the variable into our dictionary, for later use.
            self._DICT_variables[dataset_key][var_name] = variable

        ################################################################################
        # STEP 3, subset by years if asked for.  Should be a 2-tuple of integers
        if (yearspan is not None) and (var_name not in ("MASK","LAT","LON","BASINS","ELEV","CELLAREA")):
            # If we haven't gotten out the variable yet (if above steps were skipped), get it.
            variable = self._DICT_variables[dataset_key][var_name]
            yearstart, yearend = yearspan
            # Subset by yearspan.
            assert (type(yearspan) in (list, tuple)) and (len(yearspan) == 2)
            if var_name == "TIME":
                years = variable
            else:
                years = self.get_variable("TIME",dataset_key=dataset_key,yearspan=None)

            # Make sure these year boundaries are good.
            assert (yearspan[0] >= years[0]) and (yearspan[1] <= years[-1])

            i_start = yearspan[0] - years[0]
            i_end   = yearspan[1] - years[0] + 1

            if var_name == "TIME":
                return variable[i_start:i_end]
            else:
                var_subset = variable[i_start:i_end,:,:]
                if ((i_end - i_start) == 1) and (len(var_subset.shape) == 2):
                    var_subset.shape = (1,var_subset.shape[0], var_subset.shape[1])
                return var_subset

        ################################################################################
        # STEP 4: Return the variable
        return self._DICT_variables[dataset_key][var_name]

        # /_get_variable_subclass()
        ################################################################################

    def _return_years_integers(self, variable_in):
        '''The time is always given in "days since YYYY-MM-DD HH:MM:SS", which is kinda unhelpful.
        We have annual data, so return the time as integer years.'''
        numyears = variable_in.shape[0]
        # Read in "days since YYYY-MM-DD HH:MM:SS" and extract the "YYYY"
        start_year = int(variable_in.units[11:15])
        return numpy.arange(start_year, start_year + numyears, 1.0, dtype=numpy.int)

    def _get_basin_integers(self, dataset):
        '''Uses the Zwally Basins dataset to return an integer array for all the 8 major basins derived in Zwally (??)
        There are 19 basins derived in the file, we use the 8 basic ones.  All pixels off the ice sheet are given 0.'''
        if self.basin_integers is None:
            self.basin_integers = self._get_variable_SUBCLASS("BASINS")
#            dset = self._get_basin_dataset()
#            # Retrieve the basin masks
#            maskall = dset.variables["maskall"][:]
#            # Retrieve the basin numbers
#            basins = dset.variables["basins"][:]
#            # Create a zeroed-out array
#            basin_integers = numpy.zeros(maskall.shape[1:], dtype=numpy.integer)
#
#            # Fill it with each of the basin numbers
#            for i in range(basins.size):
#                basin_integers[maskall[i,:,:] == 1.0] = numpy.round(basins[i])
#
#            self.basin_integers = basin_integers

        return self.basin_integers

    def _get_elevation_data(self, export_picklefile_if_not_created=True):
        '''Return elevation data.  If not already in the variable, return from the picklefile.
        If the picklefile doesn't exist, build the picklefile by interpolating from RACMO elevation data.'''
        # (just doing NEAREST NEIGHBOR for now)
        if (self._DICT_variables[self._LIST_dataset_keys[0]]["ELEV"] is None) \
           or (self._DICT_variables[self._LIST_dataset_keys[0]]["ELEV"].shape != self.get_variable("LAT").shape):

            if os.path.exists(HIRHAM_ELEV_PICKLEFILE):
                print "Reading", os.path.split(HIRHAM_ELEV_PICKLEFILE)[-1]
                f = open(HIRHAM_ELEV_PICKLEFILE,'r')
                d = pickle.load(f)
                f.close()
            else:
                MAR = MAR_Manager()
                d = self.interpolate_between_grids(MAR, "ELEV", interpolation_type="nearest neighbor")

                if export_picklefile_if_not_created:
                    print "Writing", os.path.split(HIRHAM_ELEV_PICKLEFILE)[-1]
                    f = open(HIRHAM_ELEV_PICKLEFILE,'w')
                    pickle.dump(d, f)
                    f.close()
                    print os.path.split(HIRHAM_ELEV_PICKLEFILE)[-1], "written."

            return d

        else:
            return self._DICT_variables[self._LIST_dataset_keys[0]]["ELEV"]


########################################################################################################
#########################################  RACMO  ######################################################
########################################################################################################
class RACMO_Manager(RCM_Manager):

    _LIST_dataset_keys = ["RACMO2.3_ERA-Int_1958-2015","RACMO2.1_HadGEM2_1971-2004","RACMO2.1_HadGEM2_RCP4.5_2005-2098"]
    _LIST_variable_keys = ("LAT", "LON", "ELEV", "MASK", "CELLAREA", "RUNOFF", "SMB", "SNOW", "RAIN", "TIME", "MELT", "SUBL", "TEMP", "BASINS", "TEMPMAGNITUDE")

    def __init__(self):
        # Initialize the parent class
        super(RACMO_Manager, self).__init__()
        # Get filler dictionaires
        fnames_dict, dsets_dict, vars_dict = self._build_empty_filename_and_dataset_dictionaries()
        # Define datasets here
        self._DICT_filenames = fnames_dict # Maintains the filename for each dataset and variable in that dataset_key
        self._DICT_datasets = dsets_dict   # Maintains the NetCDF4 dataset object for each variable in that dataset_key
        self._DICT_variables = vars_dict  # Maintains the actual values for each variable in that dataset_key

        self._LAST_DATASET_KEY_USED = None
        return

    def _return_export_shapefile_folder(self):
        return RACMO_EXPORT_SHAPEFILE_FOLDER

    def _return_short_dataset_key(self, dataset_key):
        spos = re.search('_\d{4}-\d{4}', dataset_key)
        return dataset_key[:spos.start()]

    def _return_datatype_base(self):
        '''Useful for creating filenames with the base moniker in them, used by the parent class.'''
        return "RACMO"

    def _split_dataset_into_decade_yearspans(self, dataset_key):
        '''Take each dataset, and return a list of year-spans that correspond to approximate "decades" within that dataset.'''
        assert dataset_key in self._LIST_dataset_keys
        if dataset_key == "RACMO2.3_ERA-Int_1958-2015":
            return [(1958,1965),
                    (1966,1975),
                    (1976,1985),
                    (1986,1995),
                    (1996,2005),
                    (2006,2015)]
        else:
            raise ValueError("Unknown dataset key '{0}' in RACMO_Manager.".format(dataset_key))

    def _build_empty_filename_and_dataset_dictionaries(self):
        '''Build a dictionary from which to get datasets on demand.
        These won't all be immediately opened, but will be opened as needed as they're requested.'''
        _racmo_filenames_dict = {}
        _racmo_datasets_dict = {}
        _racmo_variables_dict = {}

        # for dataset_key,variable--> filename mappings.
        for dataset_key in self._LIST_dataset_keys:
            _racmo_filenames_dict[dataset_key] = {}
            _racmo_filenames_dict[dataset_key]["TIME"  ] = RACMO_FILENAMES_DICT[dataset_key]["MELT"] # "time" is the same in all datasets, just pick one to pull it from.
            _racmo_filenames_dict[dataset_key]["SUBL"  ] = RACMO_FILENAMES_DICT[dataset_key]["SUBL"]
            _racmo_filenames_dict[dataset_key]["RAIN"  ] = RACMO_FILENAMES_DICT[dataset_key]["RAIN"]
            _racmo_filenames_dict[dataset_key]["SNOW"  ] = RACMO_FILENAMES_DICT[dataset_key]["SNOW"]
            _racmo_filenames_dict[dataset_key]["SMB"   ] = RACMO_FILENAMES_DICT[dataset_key]["SMB"]
            _racmo_filenames_dict[dataset_key]["RUNOFF"] = RACMO_FILENAMES_DICT[dataset_key]["RUNOFF"]
            _racmo_filenames_dict[dataset_key]["MELT"  ] = RACMO_FILENAMES_DICT[dataset_key]["MELT"]
            _racmo_filenames_dict[dataset_key]["TEMP"  ] = RACMO_FILENAMES_DICT[dataset_key]["TEMP_2M"]
            _racmo_filenames_dict[dataset_key]["LAT"   ] = RACMO_FILENAMES_DICT[dataset_key]["GRID"]
            _racmo_filenames_dict[dataset_key]["LON"   ] = RACMO_FILENAMES_DICT[dataset_key]["GRID"]
            _racmo_filenames_dict[dataset_key]["CELLAREA"] = RACMO_FILENAMES_DICT[dataset_key]["GRID"]
            _racmo_filenames_dict[dataset_key]["MASK"  ] = RACMO_FILENAMES_DICT[dataset_key]["GRID"]
            _racmo_filenames_dict[dataset_key]["ELEV"  ] = RACMO_FILENAMES_DICT[dataset_key]["GRID"]
            _racmo_filenames_dict[dataset_key]["BASINS"] = None
            _racmo_filenames_dict[dataset_key]["TEMPMAGNITUDE"] = RACMO_FILENAMES_DICT[dataset_key]["TEMP_2M"]

            # Build an empty dictionary, for performance's sake.  Fill in the None's with datasets as they're called.
            _racmo_datasets_dict[dataset_key] = {}
            _racmo_variables_dict[dataset_key] = {}
            for key in self._LIST_variable_keys:
                _racmo_datasets_dict[dataset_key][key] = None
                _racmo_variables_dict[dataset_key][key] = None

        return _racmo_filenames_dict, _racmo_datasets_dict, _racmo_variables_dict

    def _get_variable_SUBCLASS(self, var_name, dataset_key="RACMO2.3_ERA-Int_1958-2015", yearspan=None, annual=True):
        # Define variable mappings from datasets here.
        var_name = var_name.upper()
        assert var_name in self._LIST_variable_keys

        if dataset_key is None:
            if self._LAST_DATASET_KEY_USED is None:
                dataset_key = self._LIST_dataset_keys[0]
            else:
                dataset_key = self._LAST_DATASET_KEY_USED

        assert dataset_key in self._LIST_dataset_keys

        self._LAST_DATASET_KEY_USED = dataset_key

        ################################################################################
        # STEP 1: If we haven't opened the dataset yet, do it and fill in the rest of _DICT_datasets for that dataset.
        #       The BASIN datasets are different (don't have a netCDF4 dataset)
        if (var_name != "BASINS") and (self._DICT_datasets[dataset_key][var_name] is None):

            dset = netCDF4.Dataset(self._DICT_filenames[dataset_key][var_name])

            if var_name in ("LAT","LON","ELEV","MASK","CELLAREA"):
                # Get the grid dataset, and then fill it in for all the rest of the variables.
                for vname in ("LAT","LON","ELEV","MASK","CELLAREA"):
                    self._DICT_datasets[dataset_key][vname] = dset
            elif var_name in ("MELT","TIME"):
                # In my implementation, the "TIME" variable is derived from the same dataset as the "MELT" variable.
                for vname in ("MELT","TIME"):
                    self._DICT_datasets[dataset_key][vname] = dset
            else:
                # For the rest, just use the one it comes with.
                self._DICT_datasets[dataset_key][var_name] = dset

        ################################################################################
        # STEP 2: If we haven't opened that variable yet, do it and fill in the variable.

        if self._DICT_variables[dataset_key][var_name] is None:
            dset = self._DICT_datasets[dataset_key][var_name]
            if var_name == "LAT":
                var = dset.variables['lat'][:]
            elif var_name == "LON":
                var = dset.variables['lon'][:]
            elif var_name == "ELEV":
                var = dset.variables['topography'][:]
            elif var_name == "CELLAREA":
                var = dset.variables['gridarea'][:]
            elif var_name == "RUNOFF":
                var = dset.variables['runoff'][:]
            elif var_name == "SMB":
                var = dset.variables['smb'][:]
            elif var_name == "SNOW":
                var = dset.variables['snowfall'][:]
            elif var_name == "RAIN":
                # The RACMO "precip" variable includes both rain and snow.  Subtract snow to get rain.
                precip = dset.variables['precip'][:]
                precip.shape = [d for d in precip.shape if d>1]
                var = precip - self._get_variable_SUBCLASS("SNOW", dataset_key = dataset_key, yearspan=None, annual=False)
            elif var_name == "TIME":
                var = dset.variables['time'][:] # "time" is the same in all datasets, just pick one to pull it from.
            elif var_name == "MELT":
                var = dset.variables['snowmelt'][:]
            elif var_name == "SUBL":
                var = dset.variables['subl'][:]
            elif var_name == "TEMP":
                # Convert Kelvin to Celcius
                var = dset.variables['t2m'][:] - 273.15
            elif var_name == "MASK":
                # Combine the main ice sheet and the ice caps into one variable.
                var1 = numpy.array(dset.variables['GrIS_mask'][:], dtype=numpy.bool)
                var2 = numpy.array(dset.variables['GrIS_caps_mask'][:], dtype=numpy.bool)
                var = var1 | var2
            elif var_name == "BASINS":
                var = self._get_basin_integers()

            var.shape = [d for d in var.shape if d > 1]

            self._DICT_variables[dataset_key][var_name] = var
        else:
            var = self._DICT_variables[dataset_key][var_name]

        # Get rid of 1-d variable shape dimensions

        ################################################################################
        # STEP 3, turn monthly totals into annual totals/averages, if asked for.
        # WE WILL NEED TO CHANGE THIS IF USING A DATSET THAT ISN'T EXACTLY 58 YEARS LONG.
        if annual:
            # Get rid of 1-dimensional axes here.
            var.shape = [d for d in var.shape if d>1]
            if var_name in ("RUNOFF","SMB","SNOW","RAIN","MELT","SUBL"):
                var = numpy.sum(numpy.reshape(var, (var.size/(12*var.shape[1]*var.shape[2]),12,var.shape[1], var.shape[2])), axis=1)
            elif var_name == "TEMP":
                var = numpy.mean(numpy.reshape(var, (var.size/(12*var.shape[1]*var.shape[2]),12,var.shape[1], var.shape[2])), axis=1)
            elif var_name == "TIME":
                if dataset_key == "RACMO2.3_ERA-Int_1958-2015":
                    # RACMO 2.3 TIME is "months since 1958-01-15
                    var = numpy.array(numpy.round(numpy.reshape(numpy.round(var), (58,12))[:,0] / 12.0 + 1958.0), dtype=numpy.int)
                else:
                    assert dataset_key in ("RACMO2.1_HadGEM2_1971-2004","RACMO2.1_HadGEM2_RCP4.5_2005-2098")
                    time_object = self._DICT_datasets[dataset_key]['TIME'].variables['time']
                    # RACMO 2.1 TIME is "days since 1950-01-01", using a 360-day/year calendar (odd but true)
                    # It comes with a .dtgstart attribute string which explicitly states the starting date, easy to extrapolate from there.
                    var = numpy.array(numpy.round((numpy.reshape(var, (-1,12))[:,0] - var[0])/360.0) + int(time_object.dtgstart[:4]), dtype=numpy.int)
            else:
                pass # The grid variables don't need reshaping.

        ################################################################################
        # STEP 4, subset by years if asked for.  Should be a 2-tuple of integers
        if (yearspan is not None) and (var_name not in ("MASK","LAT","LON","BASINS","ELEV","CELLAREA")):
            # If we haven't gotten out the variable yet (if above steps were skipped), get it.
            assert (type(yearspan) in (list, tuple)) and (len(yearspan) == 2)
            yearstart, yearend = yearspan
            if var_name == "TIME":
                years = var
            else:
                years = self.get_variable("TIME",dataset_key=dataset_key,yearspan=None,annual=annual)

            if annual:
                # Make sure these year boundaries are good.
                assert (yearstart >= years[0]) and (yearend <= years[-1])

                i_start = yearspan[0] - years[0]
                i_end   = yearspan[1] - years[0] + 1

                if var_name == "TIME":
                    return var[i_start:i_end]
                else:
                    var_subset = var[i_start:i_end,:,:]
                    if ((i_end - i_start) == 1) and (len(var_subset.shape) == 2):
                        var_subset.shape = (1,var_subset.shape[0], var_subset.shape[1])
                    return var_subset

            else: # Not annual.
                # I'm seriously just working on annual dat right now.
                # I'll implement the monthly stuff if/when I need it.
                raise NotImplementedError()

        ################################################################################
        # STEP 5: Return the variable
        return var

        # /_get_variable_subclass()
        ################################################################################

    def _get_basin_integers(self, dataset_key=None, export_picklefile_if_not_created=True):
        # For MAR (picklefile and resample)
        if dataset_key is None:
            dataset_key = self._LIST_dataset_keys[0]

        if (self._DICT_variables[dataset_key]["BASINS"] is None) \
           or (self._DICT_variables[dataset_key]["BASINS"].shape != self.get_variable("LAT").shape):

            if os.path.exists(RACMO_GRID_INTEGERS_PICKLEFILE):
                print "Reading", os.path.split(RACMO_GRID_INTEGERS_PICKLEFILE)[-1]
                f = open(RACMO_GRID_INTEGERS_PICKLEFILE,'r')
                d = pickle.load(f)
                f.close()
            else:
                H5 = HIRHAM5_Manager()
                d = self.interpolate_between_grids(H5, "BASINS", interpolation_type="nearest neighbor")

                if export_picklefile_if_not_created:
                    print "Writing", os.path.split(RACMO_GRID_INTEGERS_PICKLEFILE)[-1]
                    f = open(RACMO_GRID_INTEGERS_PICKLEFILE,'w')
                    pickle.dump(d, f)
                    f.close()
                    print os.path.split(RACMO_GRID_INTEGERS_PICKLEFILE)[-1], "written."

            return d

        else:
            return self._DICT_variables[self._LIST_dataset_keys[0]]["BASINS"]


########################################################################################################
#########################################  MAR  ########################################################
########################################################################################################
class MAR_Manager(RCM_Manager):

    # FOR NOW, this just uses the Montly outputs interpolated at 5 km.  There is a version of those for all the datasets
    # Later I could look into implementing different resolution data in the raw formats.
    _LIST_dataset_keys = MAR_ANNUAL_FILENAMES_DICT.keys()
    _LIST_dataset_keys.sort()
    _LIST_variable_keys = ("LAT", "LON", "ELEV", "CELLAREA", "MASK", "RUNOFF", "SMB", "SNOW", "RAIN", "TIME", "MELT", "SUBL", "TEMP", "BASINS", "TEMPMAGNITUDE")

    def __init__(self):
        # Initialize the parent class
        super(MAR_Manager, self).__init__()
        # Get filler dictionaires
        fnames_dict, vars_dict = self._build_empty_filename_and_dataset_dictionaries()
        # Define datasets here
        self._DICT_filenames = fnames_dict # Maintains the filename for each dataset and variable in that dataset_key
        self._DICT_variables = vars_dict  # Maintains the actual values for each variable in that dataset_key
        self._LAST_DATASET_KEY_USED = None
        return

    def _return_export_shapefile_folder(self):
        return MAR_EXPORT_SHAPEFILE_FOLDER

    def _return_short_dataset_key(self, dataset_key):
        spos = re.search('_\d{4}-\d{4}', dataset_key)
        return dataset_key[:spos.start()]

    def _return_datatype_base(self):
        '''Useful for creating filenames with the base moniker in them, used by the parent class.'''
        return "MAR"

    def _split_datasets_into_decade_yearspans(self, dataset_key):
        '''Take each dataset, and return a list of year-spans that correspond to approximate "decades" within that dataset.'''
        assert dataset_key in self._LIST_dataset_keys

        if dataset_key.find("1950-2005") > -1:
            return [(1950,1965),
                    (1966,1975),
                    (1976,1985),
                    (1986,1995),
                    (1996,2005)]
        elif dataset_key.find("1900-2005") > -1:
            return [(1900,1915),
                    (1916,1925),
                    (1926,1935),
                    (1936,1945),
                    (1946,1955),
                    (1956,1965),
                    (1966,1975),
                    (1976,1985),
                    (1986,1995),
                    (1996,2005)]
        elif dataset_key.find("1979-2014") > -1:
            return [(1979,1994),
                    (1995,2004),
                    (2005,2014)]
        elif dataset_key.find("1948-2015") > -1:
            return [(1948,1965),
                    (1966,1975),
                    (1976,1985),
                    (1986,1995),
                    (1996,2005),
                    (2006,2015)]
        elif dataset_key.find("1900-2010") > -1:
            return [(1900,1910),
                    (1911,1920),
                    (1921,1930),
                    (1931,1940),
                    (1941,1950),
                    (1951,1960),
                    (1961,1970),
                    (1971,1980),
                    (1981,1990),
                    (1991,2000),
                    (2001,2010)]
        elif dataset_key.find("2006-2100") > -1:
            return [(2006,2010),
                    (2011,2020),
                    (2021,2030),
                    (2031,2040),
                    (2041,2050),
                    (2051,2060),
                    (2061,2070),
                    (2071,2080),
                    (2081,2090),
                    (2091,2100)]

        else:
            raise ValueError("Unknown dataset key '{0}' in RACMO_Manager.".format(dataset_key))

    def _build_empty_filename_and_dataset_dictionaries(self):
        '''Build a dictionary from which to get datasets on demand.
        These won't all be immediately opened, but will be opened as needed as they're requested.
        STRUCTURE: self.hirham_datasets {
            ["2010GCM", "2014ERA", "2050GCM45", "2050GCM85", "2100GCM45", "2100GCM85"]
            each of which is a dict with keys of the variables.'''
        _mar_filenames_dict = {}
        _mar_variables_dict = {}

        for dset_key in self._LIST_dataset_keys:
            # In the MAR data, all the raw variables are in one data-set, hence there is not need to index filenames and datasets by variable.
            _mar_filenames_dict[dset_key] = MAR_ANNUAL_FILENAMES_DICT[dset_key]
            _mar_variables_dict[dset_key] = None

        return _mar_filenames_dict, _mar_variables_dict

    def _get_variable_SUBCLASS(self, var_name, dataset_key=None, yearspan=None, annual=None):

        if dataset_key is None:
            dataset_key = self._LAST_DATASET_KEY_USED if (self._LAST_DATASET_KEY_USED != None) else "ERA-int_1979-2014_10km"

        self._LAST_DATASET_KEY_USED = dataset_key

        # Convert all keys to upper-case.
        var_name = var_name.upper()
        # Make sure this variable is in our options.
        assert var_name in self._LIST_variable_keys

        ################################################################################
        # STEP 1: If we haven't opened the picklefile yet, do it and fill in the self._DICT_variables[dataset_key] dictionary.
        if self._DICT_variables[dataset_key] is None:
            self._read_or_build_annual_dataset_picklefile(dataset_key, force_rebuild=False)

        ################################################################################
        # STEP 2, subset by years if asked for.  Should be a 2-tuple of integers
        if (yearspan is not None) and (var_name not in ("MASK","LAT","LON","BASINS","ELEV","CELLAREA")):
            # If we haven't gotten out the variable yet (if above steps were skipped), get it.
            variable = self._DICT_variables[dataset_key][var_name]
            yearstart, yearend = yearspan
            # Subset by yearspan.
            assert (type(yearspan) in (list, tuple)) and (len(yearspan) == 2)
            if var_name == "TIME":
                years = variable
            else:
                years = self.get_variable("TIME",dataset_key=dataset_key,yearspan=None)

            # Make sure these year boundaries are good.
            assert (yearspan[0] >= years[0]) and (yearspan[1] <= years[-1])

            i_start = yearspan[0] - years[0]
            i_end   = yearspan[1] - years[0] + 1

            if var_name == "TIME":
                return variable[i_start:i_end]
            else:
                var_subset = variable[i_start:i_end,:,:]
                if ((i_end - i_start) == 1) and (len(var_subset.shape) == 2):
                    var_subset.shape = (1,var_subset.shape[0], var_subset.shape[1])
                return var_subset

        # If not subsetting timestep, just return the whole variable.
        else:
            return self._DICT_variables[dataset_key][var_name]

        # /_get_variable_subclass()
        ################################################################################

    def _build_all_annual_picklefiles_from_monthly_netCDF4_datasets(self):
        '''Run through ALL the datsets, and build an annually-averaged picklefile from all the monthly datasets.
        This is preferable to using the annual datasets provided, which don't include the "corrected" SMB values and they don't
        include 3m air temperature, etc.  NOTE: in part-tundra, part-ice pixels, this averages ONLY the on-ice portion.
        Must redo this if we want the total of the pixel (land + ice).'''
        for dset_key in self._LIST_dataset_keys:

            print "======", dset_key, "======"
            self._read_or_build_annual_dataset_picklefile(dset_key, force_rebuild=True)

        return

    def _read_or_build_annual_dataset_picklefile(self, dataset_key, force_rebuild=False):
        '''The MAR datasets that only provided monthly data, with no annual subsets, use a picklefile instead.
        The picklefile contains copy of the fully-populated self._DICT_variables dictionary, containing all dataset[var] key/value pairs.

        If the picklefie exists, just read it and return to the calling function.  If it doesn't exist, read all the monthly files, build it and save it to disk.
        If "force_rebuild", delete the old picklefile and rebuild it.'''

        filename = self._DICT_filenames[dataset_key]


        if os.path.exists(filename):
            if force_rebuild:
                # Delete the file here
                print "Deleting old", os.path.split(filename)[-1]
                os.remove(filename)
                pass
            else:
                f = open(filename, 'r')
                print "Reading", os.path.split(filename)[-1]
                d = pickle.load(f)
                f.close()
                self._DICT_variables[dataset_key] = d
                return d

        picklefilename = filename
        foldername = os.path.split(picklefilename)[0]
        # Get all the filenames that include "monthly" and end in ".nc"
        annual_filenames = [os.path.join(foldername,fname) for fname in os.listdir(foldername) if ((fname.find("monthly") != -1) and (fname[-3:] == ".nc"))]
        annual_filenames.sort()
        annual_datasets = [None] * len(annual_filenames)
        for i,fname in enumerate(annual_filenames):
            print "Reading", os.path.split(fname)[-1]
            annual_datasets[i] = netCDF4.Dataset(fname)

        YEARS = numpy.array([int(fname[-7:-3]) for fname in annual_filenames], dtype=numpy.int)
        N_years = len(YEARS)
        # Initialize to a new, empty dictionary.
        self._DICT_variables[dataset_key] = {}
        for var_key in self._LIST_variable_keys:
            self._DICT_variables[dataset_key][var_key] = None

        for i, dataset in enumerate(annual_datasets):
            print YEARS[i]
            for var_name in self._LIST_variable_keys:

                var_key = {"TIME"     : None    ,
                           "SUBL"     :"SU"     ,
                           "RAIN"     :"RF"     ,
                           "SNOW"     :"SF"     ,
                           "RUNOFF"   :"RUcorr" ,
                           "SMB"      :"SMBcorr",
                           "MELT"     :"MEcorr" ,
                           "TEMP"     :"TTcorr" ,
                           "MASK"     :"MSK_MAR",
                           "BASINS"   : None    ,
                           "LAT"      :"LAT"    ,
                           "LON"      :"LON"    ,
                           "CELLAREA" :"AREA"   ,
                           "ELEV"     :"SRF_bam13",
                           "TEMPMAGNITUDE": None,
                           } [var_name]


                if var_name == "BASINS":
                    if i == 0:
                        self._DICT_variables[dataset_key][var_name] = self._get_basin_integers()
                    continue
                elif var_name == "TIME":
                    if i == 0:
                        self._DICT_variables[dataset_key][var_name] = YEARS
                    continue
                else:
                    # Individual post-processing of variable data before returning.
                    # Some variables come as a "masked" variable. Just give us the basic data.
                    variable = dataset.variables[var_key][:]
                    variable = variable[:]
                    if numpy.ma.is_masked(variable):
                        is_masked = True
                        variable = variable.data[:]
                        empty_val = variable[0,0] if (len(variable.shape) == 2) else variable[0,0,0]
                    else:
                        is_masked = False


                if var_name in ("LAT", "LON", "ELEV", "CELLAREA"):
                    if i == 0:
                        self._DICT_variables[dataset_key][var_name] = variable
                    else:
                        assert numpy.all(self._DICT_variables[dataset_key][var_name] == variable)


                elif var_name == "MASK":
                    emptyval = variable[0,0]
                    mask = ((variable != emptyval) & (variable > 1.5))
                    # In MAR, mask out the portions of the dataset on Ellesmere Island, in the far northwest of the grid.
                    lat, lon = self.get_lats_lons()
                    mask[((lat > 78.0) & (lon < -75.0)) | ((lat > 79.0) & (lon < -70.0)) | ((lat > 80.0) & (lon < -68.0)) | ((lat > 81.0) & (lon < -65.0))] = False

                    if i == 0:
                        self._DICT_variables[dataset_key][var_name] = mask
                    else:
                        assert numpy.all(self._DICT_variables[dataset_key][var_name] == mask)

                elif var_name not in ["BASINS", "TIME"]:
                    if i == 0:
                        var_array = numpy.empty(([N_years] + list(variable.shape[1:])), dtype=variable.dtype)
                    else:
                        var_array = self._DICT_variables[dataset_key][var_name]

                    # For temperature, take the mean of the monthly values for the annual value.
                    if var_name == "TEMP":
                        var_array[i,:,:] = numpy.mean(variable, axis=0)
                    # For all the other variables, take the sum of the montly values.
                    else:
                        var_array[i,:,:] = numpy.sum(variable, axis=0)

                        if is_masked:
                            # The sum screws up the empty-values, so flatten the 2-3 axes of the array to create an empty mask, then re-fill back in the empty-values.
                            orig_shape = variable.shape
                            orig_shape_var_array = var_array.shape
                            variable.shape = variable.shape[0], variable.shape[1]*variable.shape[2]
                            var_array.shape = var_array.shape[0], var_array.shape[1]*var_array.shape[2]

                            empty_mask = (variable[0,:] == empty_val)
                            var_array[i,empty_mask] = empty_val

                            variable.shape = orig_shape
                            var_array.shape = orig_shape_var_array

                    self._DICT_variables[dataset_key][var_name] = var_array
                else:
                    raise NotImplementedError(var_name)


            dataset.close()

        for var_name in self._LIST_variable_keys:
            variable = self._DICT_variables[dataset_key][var_name]
            mask = self._DICT_variables[dataset_key]["MASK"]
            if var_name == "BASINS":
                mean = None
            elif var_name == "TIME":
                mean = numpy.mean(variable)
            elif var_name in ("LAT", "LON", "CELLAREA", "ELEV", "MASK"):
                mean = numpy.mean(variable[mask])
            else:
                variable = variable.copy()
                variable.shape = variable.shape[0], variable.shape[1]*variable.shape[2]
                mean = numpy.mean(variable[:,mask.flatten()])

            print '\t', var_name, self._DICT_variables[dataset_key][var_name].shape, mean

        print "Writing", os.path.split(picklefilename)[-1]
        f = open(picklefilename, 'w')
        pickle.dump(self._DICT_variables[dataset_key], f)
        f.close()

        return self._DICT_variables[dataset_key]

    def _get_basin_integers(self, export_picklefile_if_not_created=True):
        # For MAR (picklefile and resample)
        if (self._DICT_variables[self._LIST_dataset_keys[0]] is None) \
           or (self._DICT_variables[self._LIST_dataset_keys[0]]["BASINS"] is None) \
           or (self._DICT_variables[self._LIST_dataset_keys[0]]["BASINS"].shape != self.get_variable("LAT").shape):

            if os.path.exists(MAR_BASIN_INTEGERS_PICKLEFILE):
                print "Reading", os.path.split(MAR_BASIN_INTEGERS_PICKLEFILE)[-1]
                f = open(MAR_BASIN_INTEGERS_PICKLEFILE,'r')
                d = pickle.load(f)
                f.close()
            else:
                H5 = HIRHAM5_Manager()
                d = self.interpolate_between_grids(H5, "BASINS", interpolation_type="nearest neighbor")

                if export_picklefile_if_not_created:
                    print "Writing", os.path.split(MAR_BASIN_INTEGERS_PICKLEFILE)[-1]
                    f = open(MAR_BASIN_INTEGERS_PICKLEFILE,'w')
                    pickle.dump(d, f)
                    f.close()
                    print os.path.split(MAR_BASIN_INTEGERS_PICKLEFILE)[-1], "written."

            return d

        else:
            return self._DICT_variables[self._LIST_dataset_keys[0]]["BASINS"]


########################################################################################################
##################################  STANDALONE FUNCTIONS  ##############################################
########################################################################################################


def UNIT_TEST_loading_data():
    '''Perform a unit test of loading data from the data set, for each of the keys.'''
    # THIS WORKS. LOADING DATA IN HIRHAM5 HAS PASSED.
    print "_________HIRHAM5__________"
    H5 = HIRHAM5_Manager()
    for dset_key in H5._LIST_dataset_keys:
        print
        for var_name in H5._LIST_variable_keys:
            print dset_key + ",", var_name,
            variable = H5.get_variable(var_name, dataset_key=dset_key)
            print variable.shape, variable.dtype,

            yearspan = {'2010GCM': (1995,2003),
                        '2014ERA': (1989,2009),
                        '2050GCM45': (2035,2045),
                        '2050GCM85': (2035,2045),
                        '2100GCM45': (2085,2095),
                        '2100GCM85': (2085,2095),
                        } [dset_key]

            variable = H5.get_variable(var_name, dataset_key=dset_key, yearspan = yearspan)
            print variable.shape

    print "_________RACMO__________"
    RACMO = RACMO_Manager()
    for dset_key in RACMO._LIST_dataset_keys:

        for var_name in RACMO._LIST_variable_keys:
            print dset_key, var_name,
            variable = RACMO.get_variable(var_name, dataset_key = dset_key)
            print variable.shape, variable.dtype,

            yearspan = {"RACMO2.3_ERA-Int_1958-2015": (1963,1992),
                        "RACMO2.1_HadGEM2_1971-2004": (1989,2003),
                        "RACMO2.1_HadGEM2_RCP4.5_2005-2098": (2035,2045),
                        } [dset_key]

            variable = RACMO.get_variable(var_name, dataset_key=dset_key, yearspan = yearspan)

            mask = RACMO.get_variable("MASK", dataset_key=dset_key)
            if var_name == "BASINS":
                mean = "foobar"
            elif len(variable.shape) == 2:
                mean = numpy.mean(variable[mask])
            elif len(variable.shape) == 3:
                mean = numpy.mean(variable[:,mask])
            elif var_name == "TIME":
                mean = numpy.mean(variable)
            else:
                mean = "HUH???????"

            print variable.shape, mean

    print "_________MAR__________"
    MAR = MAR_Manager()
    for dset_key in MAR._LIST_dataset_keys:
        print
        mask = MAR.get_variable("MASK", dataset_key=dset_key)
        for var_name in MAR._LIST_variable_keys:
            print dset_key + ",", var_name,
            variable = MAR.get_variable(var_name, dataset_key=dset_key)
            if var_name == "BASINS":
                mean = "foobar"
            elif len(variable.shape) == 2:
                mean = numpy.mean(variable[mask])
            elif len(variable.shape) == 3:
                mean = numpy.mean(variable[:,mask])
            elif var_name == "TIME":
                mean = numpy.mean(variable)
            else:
                mean = "HUH???????"
            print variable.shape, variable.dtype, mean,

            yearspan = {'CanESM2-histo_1950-2005_25km': (1995,2003),
                        'CanESM2-rcp45_2006-2100_25km': (2045,2055),
                        'CanESM2-rcp85_2006-2100_25km': (2045,2055),
                        'ERA-int_1979-2014_10km'      : (1995,2003),
                        'MIROC5-histo_1900-2005_25km' : (1995,2003),
                        'MIROC5-rcp45_2006-2100_25km' : (2045,2055),
                        'MIROC5-rcp85_2006-2100_25km' : (2045,2055),
                        'NorESM1-histo_1950-2005_25km': (1995,2003),
                        'NorESM1-rcp45_2006-2100_25km': (2045,2055),
                        'NorESM1-rcp85_2006-2100_25km': (2045,2055),
                        'ERA_20C_1900-2010_20km'      : (1995,2003),
                        'NCEPv1_1948-2015_20km'       : (1995,2003),
                        } [dset_key]

            variable = MAR.get_variable(var_name, dataset_key=dset_key, yearspan = yearspan)
            print variable.shape

def quick_look_at_variables():
    for M in (HIRHAM5_Manager(), MAR_Manager(), RACMO_Manager()):
        MASK = M.get_variable("MASK", dataset_key=M._LIST_dataset_keys[0])
        for dataset_key in M._LIST_dataset_keys:
            print
            print M.get_manager_name(), dataset_key
            for var_name in ["LAT","LON","ELEV","MELT","RAIN","SNOW","TEMP"]:
                data = M.get_variable(var_name, dataset_key = dataset_key)
                if var_name in ["MELT","RAIN","SNOW","TEMP"]:
                    data = numpy.mean(data, axis=0)


                print "   ", var_name, numpy.mean(data[MASK])

            print "   ", "EMI", numpy.mean(M.compute_EMI_values(dataset_key=dataset_key)[MASK])

def plot_EMI_of_icebridge_pixels():
    '''Plot the average EMI yearly trend for all IceBridge pixels in each Reanalysis dataset.'''
    MAR, RACMO, HIRHAM = MAR_Manager(), RACMO_Manager(), HIRHAM5_Manager()
    managers = [HIRHAM, MAR, MAR, RACMO]
    dataset_keys = ['2014ERA',
                    'ERA-int_1979-2014_10km',
                    'NCEPv1_1948-2015_20km',
                    "RACMO2.3_ERA-Int_1958-2015"]

    LEGEND_KEYS = ["HIRHAM 5, ERA-Interim",
                   "MAR 3.5.2, ERA-Interim",
                   "MAR 3.5.2, NCEPv1",
                   "RACMO 2.3, ERA-Interim"]

    # Colors are supposed to be in 0-1 range
    COLORS = numpy.array([(0,158,115), #blue-ish green
                          (230,159,0), #orange
                          (0,114,178), #blue, (below: reddish-purple)
                          (204,121,167)]) \
                          / 255.0

    years_lists = [None] * len(dataset_keys)
    EMI_lists = [None] * len(dataset_keys)

    fig, ax = plt.subplots(1,1)

    for i, (manager, dataset_key, legend_key, color) in \
                enumerate(zip(managers, dataset_keys, LEGEND_KEYS, COLORS)):
        years, EMI = manager.compute_EMI_of_IceBridge_pixels(dataset_key)

        # Save these populations for the next plot.
        years_lists[i] = years
        EMI_lists[i] = EMI

        yearmask = (years >= 1970)
        years = years[yearmask]
        EMI_trend = numpy.mean(EMI[yearmask,:], axis=1)

        ax.plot(years, EMI_trend, label=legend_key, color=color)

    ax.set_xlabel("Year")
    ax.set_ylabel("Excess Melt (mm w.e. yr$^{-1}$)")
    ax.legend()

    figname = "EMI_Trends_Reanalysis.tif"

    print "Plotting", figname
    fig.savefig(os.path.join(STATISTICS_OUTPUT_FOLDER,figname), dpi=600)

    plt.close()
    plt.cla()

    fig, ax = plt.subplots(1,1)


def interpolate_HIRHAM_runoff_totals(input_years, input_annual_runoff):
    '''The HIRHAM GCM datasets only have runoff for years 2000-2010, 2040-2050, and 2090-2100.
    This means the cumulative runoff totals reset at zero every period and are invalid.

    Given the years HIRHAM annual runoff is given for, interpolate annual runoff totals linearly,
    And use those totals to calculate cumulative runoff from 1990 through 2100.
    Return both interpolated annual values and cumulative values, 1990 through 2100.'''
    assert (min(input_years) == 2000) and (max(input_years) == 2100) and (len(input_years) == 33)

    output_years = numpy.arange(1990,2101,1, dtype=numpy.int)
    output_annual_runoff = numpy.zeros((len(output_years),),dtype=numpy.float)

    # First, fill in the annual runoff totals with the values from the start, SKIPPING the first
    # year in each interval which is necessarily zero.
    mask_input_2001_2010 = (input_years >= 2001) & (input_years <= 2010)
    mask_input_2041_2050 = (input_years >= 2041) & (input_years <= 2050)
    mask_input_2091_2100 = (input_years >= 2091) & (input_years <= 2100)

    mask_output_2001_2010 = (output_years >= 2001) & (output_years <= 2010)
    mask_output_2041_2050 = (output_years >= 2041) & (output_years <= 2050)
    mask_output_2091_2100 = (output_years >= 2091) & (output_years <= 2100)


    output_annual_runoff[mask_output_2001_2010] = input_annual_runoff[mask_input_2001_2010]
    output_annual_runoff[mask_output_2041_2050] = input_annual_runoff[mask_input_2041_2050]
    output_annual_runoff[mask_output_2091_2100] = input_annual_runoff[mask_input_2091_2100]

    # Since year 2000 & 2001 totals are zero, we will just leave 1990-1999 at zero.
    # Just verify this is true.
    assert numpy.count_nonzero(input_annual_runoff[(input_years <= 2001)]) == 0
    assert numpy.count_nonzero(output_annual_runoff[output_years <= 2001]) == 0

    # Find the mean for the first few years at the start and end of each interpolated period.
    mean_annual_runoff_2008_2010 = numpy.mean(input_annual_runoff[(input_years >= 2008) & (input_years <= 2010)])
    mean_annual_runoff_2041_2043 = numpy.mean(input_annual_runoff[(input_years >= 2041) & (input_years <= 2043)])
    mean_annual_runoff_2048_2050 = numpy.mean(input_annual_runoff[(input_years >= 2048) & (input_years <= 2050)])
    mean_annual_runoff_2091_2093 = numpy.mean(input_annual_runoff[(input_years >= 2091) & (input_years <= 2093)])

    step_change_2011_2040 = (mean_annual_runoff_2041_2043 - mean_annual_runoff_2008_2010) / float(2041-2010)
    mean_annual_interpolated_2011_2040 = numpy.arange(mean_annual_runoff_2008_2010 + step_change_2011_2040, mean_annual_runoff_2041_2043 - 1e-10, step_change_2011_2040)

    step_change_2051_2090 = (mean_annual_runoff_2091_2093 - mean_annual_runoff_2048_2050) / float(2091-2050)
    mean_annual_interpolated_2051_2090 = numpy.arange(mean_annual_runoff_2048_2050 + step_change_2051_2090, mean_annual_runoff_2091_2093 - 1e-10, step_change_2051_2090)

    # Fill in the gaps with the interpolated values.
    output_annual_runoff[(output_years >= 2011) & (output_years <= 2040)] = mean_annual_interpolated_2011_2040[:]
    output_annual_runoff[(output_years >= 2051) & (output_years <= 2090)] = mean_annual_interpolated_2051_2090[:]

    # Sum up the cumulative totals
    output_cumulative_runoff = numpy.cumsum(output_annual_runoff)

    return output_years, output_annual_runoff, output_cumulative_runoff


def plot_ICE_SLABS_decadal_outputs_v3():
    '''Plot the ice slab areas both for RCM outputs and GCM outputs.

    If GCM_greyscale is true, plot the GCM outputs as a greyscale span instead of a colored spaghetti plot.
    When doing this, omit the RACMO one since we only have it for RCP 4.5 and not 8.5'''
    data_array = RCM_Manager()._read_EMI_ICE_SLAB_decadal_outputs_v3()

    # 1: Get the range of REANALYSIS datasets at present-day, 2004-2013
    REANALYSIS_KEYS = ["HIRHAM2014ERA",
                       "MARERA-int"   ,
                       "MARNCEPv1"    ,
                       "RACMORACMO2.3_ERA-Int",
                       ]

    LEGEND_KEYS = ["HIRHAM 5, ERA-Int",
                   "MAR 3.5.2, ERA-Int",
                   "MAR 3.5.2, NCEPv1",
                   "RACMO 2.3, ERA-Int"]

    COLORS = ["green", "red", "blue", "purple"]

    # A handy datatype for the final summary table.
    # There are two numbers for each entry since, for the paer, we're summarizing over two periods for each dataset.
    summary_dt = numpy.dtype([('DATASET'    ,numpy.string_, 21),
                              ('YEAR_START' ,numpy.int, 2),
                              ('YEAR_END'   ,numpy.int, 2),
                              ('AREA_TREND_KM_A',numpy.float, 2),
                              ('AREA_KM2',numpy.float, 2),
                              ('RUNOFF_TREND_GT_A',numpy.float, 2),
                              ('RUNOFF_ACC_GT_A2',numpy.float, 2),
                              ('RUNOFF_TOTAL_GT',numpy.float, 2),
                            ])
    summary_dict = {}

    fig, axes = plt.subplots(3,2,figsize=(5.0,6.5)) #, gridspec_kw = {'wspace':0.34, 'hspace':0.23})
    ax_area = axes[0][0]
    ax_runoff = axes[0][1]


    data13_min = 1e12
    data13_max = -1e12
    # Only in the runoff estimates are accelerations calculated, those fields left blank for the area calculations
    data_string = "DATASET_KEY,YEARSTART,YEAREND,AREA_OR_RUNOFF,SLOPE,INTERCEPT,R_VALUE,P_VALUE,STDERR," + \
                  "ACC_SLOPE,ACC_INT,R_VALUE,P_VALUE,STDERR\n"

    years_RE = numpy.arange(1990,2014,1,dtype=numpy.int)
    areas_RE = numpy.zeros((4,len(years_RE)),dtype=numpy.float)

    for i,(key, legend_key, color) in enumerate(zip(REANALYSIS_KEYS, LEGEND_KEYS, COLORS)):

        # Add an entry for our table summary here.
        if key not in summary_dict.keys():
            summary_dict[key] = numpy.empty([1,], dtype=summary_dt)[0]
            summary_dict[key]['DATASET'] = key

        # Get the dataset rows that match each dataset key.
        rows_mask = numpy.array([(row_value.find(key) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool)
        years = data_array['YEAR_END']
        # Get just the years beyond 1990.
        mask_1990 = ((years >= 1990) & (years <= 2013))
        mask_1990_2000 = ((years >= 1990) & (years <= 2000)) & rows_mask
        mask_2001_2013 = ((years >= 2001) & (years <= 2013)) & rows_mask
        mask_total = rows_mask & mask_1990

        # Mask the years and area totals to begin at 1990
        years_since_1990 = years[mask_total]
        area = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_total]
        cum_runoff = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][mask_total]

        # Save the reanalysis areas for later use
        assert numpy.all( years_RE == years_since_1990 )
        areas_RE[i,:] = area

        data13_min = min(data13_min, area[-1])
        data13_max = max(data13_max, area[-1])

        # Divide by 1000, plot as 10^3 km^2
        ax_area.plot(years_since_1990, area/1000., color=color, label=legend_key, lw=1)
        # Plot runoff
        ax_runoff.plot(years_since_1990, cum_runoff/360.0, color=color, label=legend_key, lw=1)

        # Get regression stats out.
        # YEARS 1990 - 2000
        years_1990_2000 = years[mask_1990_2000]
        area_1990_2000 = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_1990_2000]
        cum_runoff_1990_2000 = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][mask_1990_2000]
        ann_runoff_1990_2000 = data_array['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000'][mask_1990_2000]

        slope_1990_2000, intercept_1990_2000, rvalue_1990_2000, pvalue_1990_2000, stderr_1990_2000 = \
            scipy.stats.linregress(years_1990_2000, area_1990_2000)

        data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(
                                                    key,
                                                    years_1990_2000[0],
                                                    years_1990_2000[-1],
                                                    "AREA",
                                                    slope_1990_2000,
                                                    intercept_1990_2000,
                                                    rvalue_1990_2000,
                                                    pvalue_1990_2000,
                                                    stderr_1990_2000)

        cslope_1990_2000, cintercept_1990_2000, crvalue_1990_2000, cpvalue_1990_2000, cstderr_1990_2000 = \
            scipy.stats.linregress(years_1990_2000, cum_runoff_1990_2000)
        aslope_1990_2000, aintercept_1990_2000, arvalue_1990_2000, apvalue_1990_2000, astderr_1990_2000 = \
            scipy.stats.linregress(years_1990_2000, ann_runoff_1990_2000)

        # Fill in the table summary here.
        summary_dict[key]['YEAR_START'       ][0] = years_1990_2000[ 0]
        summary_dict[key]['YEAR_END'         ][0] = years_1990_2000[-1]
        summary_dict[key]['AREA_TREND_KM_A'  ][0] = slope_1990_2000
        summary_dict[key]['AREA_KM2'         ][0] = area_1990_2000[-1]
        summary_dict[key]['RUNOFF_TREND_GT_A'][0] = cslope_1990_2000
        summary_dict[key]['RUNOFF_ACC_GT_A2' ][0] = aslope_1990_2000
        summary_dict[key]['RUNOFF_TOTAL_GT'  ][0] = cum_runoff_1990_2000[-1]

        data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n".format(
                                                    key,
                                                    years_1990_2000[0],
                                                    years_1990_2000[-1],
                                                    "RUNOFF",
                                                    cslope_1990_2000,
                                                    cintercept_1990_2000,
                                                    crvalue_1990_2000,
                                                    cpvalue_1990_2000,
                                                    cstderr_1990_2000,
                                                    aslope_1990_2000,
                                                    aintercept_1990_2000,
                                                    arvalue_1990_2000,
                                                    apvalue_1990_2000,
                                                    astderr_1990_2000)

        # YEARS 2001-2013
        years_2001_2013 = years[mask_2001_2013]
        area_2001_2013 = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_2001_2013]
        cum_runoff_2001_2013 = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][mask_2001_2013]
        ann_runoff_2001_2013 = data_array['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000'][mask_2001_2013]

        slope_2001_2013, intercept_2001_2013, rvalue_2001_2013, pvalue_2001_2013, stderr_2001_2013 = \
            scipy.stats.linregress(years_2001_2013, area_2001_2013)

        data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format(
                                                    key,
                                                    years_2001_2013[0],
                                                    years_2001_2013[-1],
                                                    "AREA",
                                                    slope_2001_2013,
                                                    intercept_2001_2013,
                                                    rvalue_2001_2013,
                                                    pvalue_2001_2013,
                                                    stderr_2001_2013)

        cslope_2001_2013, cintercept_2001_2013, crvalue_2001_2013, cpvalue_2001_2013, cstderr_2001_2013 = \
            scipy.stats.linregress(years_2001_2013, cum_runoff_2001_2013)
        aslope_2001_2013, aintercept_2001_2013, arvalue_2001_2013, apvalue_2001_2013, astderr_2001_2013 = \
            scipy.stats.linregress(years_2001_2013, ann_runoff_2001_2013)

        data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n".format(
                                                    key,
                                                    years_2001_2013[0],
                                                    years_2001_2013[-1],
                                                    "RUNOFF",
                                                    cslope_2001_2013,
                                                    cintercept_2001_2013,
                                                    crvalue_2001_2013,
                                                    cpvalue_2001_2013,
                                                    cstderr_2001_2013,
                                                    aslope_2001_2013,
                                                    aintercept_2001_2013,
                                                    arvalue_2001_2013,
                                                    apvalue_2001_2013,
                                                    astderr_2001_2013)

        summary_dict[key]['YEAR_START'       ][1] = years_2001_2013[ 0]
        summary_dict[key]['YEAR_END'         ][1] = years_2001_2013[-1]
        summary_dict[key]['AREA_TREND_KM_A'  ][1] = slope_2001_2013
        summary_dict[key]['AREA_KM2'         ][1] = area_2001_2013[-1]
        summary_dict[key]['RUNOFF_TREND_GT_A'][1] = cslope_2001_2013
        summary_dict[key]['RUNOFF_ACC_GT_A2' ][1] = aslope_2001_2013
        summary_dict[key]['RUNOFF_TOTAL_GT'  ][1] = cum_runoff_2001_2013[-1]


    # Put a grey dashed line at zero
    ax_area.axhline(0, color='lightgrey', linestyle="dashed", zorder=-1, linewidth=1)
    ax_runoff.axhline(0, color='lightgrey', linestyle="dashed", zorder=-1, linewidth=1)

    # Remove the top and right spines
    ax_area.spines['right'].set_visible(False)
    ax_area.spines['top'].set_visible(False)
    ax_runoff.spines['right'].set_visible(False)
    ax_runoff.spines['top'].set_visible(False)

    print "2013 areas (km2):", data13_min, data13_max
    ax_area.legend(fontsize="xx-small", loc="upper left")
    ax_area.set_ylabel("Area (10$^3$ km$^2$)")
    ax_runoff.set_ylabel("Runoff (mm SLE)")

    # Get a time-series of each GCM forced dataset, plot that up 'till 2100.
    GCM_KEYS_45 = [("HIRHAM2010GCM","HIRHAM2050GCM45","HIRHAM2100GCM45"),
                "MARNorESM1-rcp45",
                "MARMIROC5-rcp45",
                "MARCanESM2-rcp45",
                "RACMORACMO2.1_HadGEM2_RCP4.5",
                ]

    LEGEND_KEYS_45 = ["HIRHAM 5, ECEarth",
                   "MAR 3.5.2, NorESM1",
                   "MAR 3.5.2, MIROC5",
                   "MAR 3.5.2, CanESM2",
                   "RACMO 2.1, HadGEM2",
                   ]

    COLORS_45 = ["green","goldenrod","blue","red","purple"]

    GCM_KEYS_85 = [("HIRHAM2010GCM","HIRHAM2050GCM85","HIRHAM2100GCM85"),
                "MARNorESM1-rcp85",
                "MARMIROC5-rcp85",
                "MARCanESM2-rcp85",
                ]

    LEGEND_KEYS_85 = ["HIRHAM 5, ECEarth",
                   "MAR 3.5.2, NorESM1",
                   "MAR 3.5.2, MIROC5",
                   "MAR 3.5.2, CanESM2",
                   ]

    COLORS_85 = COLORS_45[0:4]

    axes_area = axes[1][0], axes[1][1]
    axes_runoff = axes[2][0], axes[2][1]

    #########################
    ## Plot GCM ice slab area
    #########################
    for ax_i,ax in enumerate(axes_area):

        # Choose whether to plot the RCP4.5 or RCP8.5 data, depending upon the axis.
        GCM_KEYS    = GCM_KEYS_45    if (ax_i == 0) else GCM_KEYS_85
        LEGEND_KEYS = LEGEND_KEYS_45 if (ax_i == 0) else LEGEND_KEYS_85
        COLORS      = COLORS_45      if (ax_i == 0) else COLORS_85

        years_GCMs = numpy.arange(1990,2101,1,dtype=numpy.int)
        areas_GCMs = numpy.zeros((4, len(years_GCMs)), dtype=numpy.float )

        for i,key in enumerate(GCM_KEYS):

            if key not in summary_dict.keys():
                summary_dict[key] = numpy.empty([1,], dtype=summary_dt)[0]
                summary_dict[key]['DATASET'] = key if type(key) == str else key[-1]

            color = COLORS[i]
            legend_key = LEGEND_KEYS[i]

            # If we're
            if type(key) == str:
                rows_mask = numpy.array([(row_value.find(key) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool)
            else:
                assert type(key) == tuple
                rows_mask = numpy.array([(row_value.find(key[0]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool) | \
                            numpy.array([(row_value.find(key[1]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool) | \
                            numpy.array([(row_value.find(key[2]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool)

            years = data_array['YEAR_END']
            mask_1990 = (years >= (1990 if (legend_key[0:6] != "HIRHAM") else 2000))
            mask_total = rows_mask & mask_1990

            years_since_1990 = years[mask_total]
            area = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_total]

            # For HIRHAM, get years and interpolate values between them.
            if type(key) == tuple:
                HIRHAM_yearspans = ((2000,2010),(2040,2050),(2090,2100))
                old_ye = None
                for y_i,(ys,ye) in enumerate(HIRHAM_yearspans):
                    # Fill in the HIRHAM values for each of these
                    y_mask = ((years_since_1990 >= ys) & (years_since_1990 <= ye))
                    y_mask_GCMs = ((years_GCMs >= ys) & (years_GCMs <= ye))
                    areas_GCMs[i,y_mask_GCMs] = area[y_mask]

                    # Only put the legend for the first line-segment, skip for the rest.
                    if y_i == 0:
                        ax.plot(years_since_1990[y_mask], area[y_mask]/1000., color=color, lw=1, zorder=2, label=legend_key)
                    else:
                        ax.plot(years_since_1990[y_mask], area[y_mask]/1000., color=color, lw=1, zorder=2)

                    if y_i < 2:
                        ax.plot([ye,HIRHAM_yearspans[y_i+1][0]], numpy.array([area[y_mask][-1],area[numpy.where(y_mask)[0][-1]+1]])/1000., color=color, lw=1, zorder=1, linestyle="dotted")

                    # Fill in the gaps between these spans. Zero for the starting gaps, interpolate for the rest.
                    if y_i == 0:
                        areas_GCMs[i,(years_GCMs < 2000)] == 0.0
                    else:
                        years_gap_mask = ((years_GCMs > old_ye) & (years_GCMs < ys))
                        start_area = areas_GCMs[i,(years_GCMs == old_ye)]
                        end_area   = areas_GCMs[i,(years_GCMs == ys)]
                        areas_GCMs[i,years_gap_mask] = (end_area - start_area) * (numpy.array(years_GCMs[years_gap_mask],dtype=numpy.float) - old_ye)/float(ys-old_ye) + start_area

                    old_ye = ye
            else:
                # Skip RACMO dataset. Since the RACMO is the very last key in the RCP4.5 keys, this shouldn't screw up the indexing
                if (key[:5] != "RACMO"):
                    areas_GCMs[i,:] = area
                    ax.plot(years_since_1990, area/1000., color=color,zorder=2, lw=1, label=legend_key)

            # Get regression stats out.
            # YEARS 1990 - 2050
            mask_1990_2050 = ((years >= 1990) & (years <= 2050)) & rows_mask
            years_1990_2050 = years[mask_1990_2050]
            area_1990_2050 = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_1990_2050]
            slope_1990_2050, intercept_1990_2050, rvalue_1990_2050, pvalue_1990_2050, stderr_1990_2050 = \
                scipy.stats.linregress(years_1990_2050, area_1990_2050)

            data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format( \
                                        key if (type(key) == str) else key[-1],
                                        years_1990_2050[0],
                                        years_1990_2050[-1],
                                        "AREA",
                                        slope_1990_2050,
                                        intercept_1990_2050,
                                        rvalue_1990_2050,
                                        pvalue_1990_2050,
                                        stderr_1990_2050)

            # Fill in the table summary here.
            summary_dict[key]['YEAR_START'       ][0] = years_1990_2050[ 0]
            summary_dict[key]['YEAR_END'         ][0] = years_1990_2050[-1]
            summary_dict[key]['AREA_TREND_KM_A'  ][0] = slope_1990_2050
            summary_dict[key]['AREA_KM2'         ][0] = area_1990_2050[-1]

            # YEARS 2050-2100
            mask_2050_2100 = ((years >= (2051 if (type(key) == str) else 2041)) & (years <= 2100)) & rows_mask
            years_2050_2100 = years[mask_2050_2100]
            area_2050_2100 = data_array['ICE_SLAB_CONTINUOUS_AREA_KM2'][mask_2050_2100]
            slope_2050_2100, intercept_2050_2100, rvalue_2050_2100, pvalue_2050_2100, stderr_2050_2100 = \
                scipy.stats.linregress(years_2050_2100, area_2050_2100)

            data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8}\n".format( \
                                        key if (type(key) == str) else key[-1],
                                        years_2050_2100[0],
                                        years_2050_2100[-1],
                                        "AREA",
                                        slope_2050_2100,
                                        intercept_2050_2100,
                                        rvalue_2050_2100,
                                        pvalue_2050_2100,
                                        stderr_2050_2100)

            # Fill in the table summary here.
            summary_dict[key]['YEAR_START'       ][1] = years_2050_2100[ 0]
            summary_dict[key]['YEAR_END'         ][1] = years_2050_2100[-1]
            summary_dict[key]['AREA_TREND_KM_A'  ][1] = slope_2050_2100
            summary_dict[key]['AREA_KM2'         ][1] = area_2050_2100[-1]


        # Shade dark grey between min and max Reanalysis outputs
        ax.fill_between(years_RE, numpy.amin(areas_RE, axis=0)/1000., numpy.amax(areas_RE, axis=0)/1000., color="black",zorder=4)

        # Horizontal grid lines every 200K, light grey.
        for hline_height in (0,200,400,600):
            ax.axhline(hline_height, color='lightgrey', zorder=0, linestyle="dashed", linewidth=1)

        if ax_i == 0:
            ax.legend(fontsize="xx-small", loc="upper left")

        # Set the x-label years to true.
        [label.set_visible(True) for label in ax.get_xticklabels()]
        # Set the x-label years to true.
        [label.set_visible(True) for label in ax.get_yticklabels()]


        ax.set_ylabel("Area (10$^3$ km$^2$)")

        # Remove the top and right spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

    #########################
    ## Plot GCM runoff
    #########################
    for ax_i,ax in enumerate(axes_runoff):
        # Choose whether to plot the RCP4.5 or RCP8.5 data, depending upon the axis.
        GCM_KEYS    = GCM_KEYS_45    if (ax_i == 0) else GCM_KEYS_85
        LEGEND_KEYS = LEGEND_KEYS_45 if (ax_i == 0) else LEGEND_KEYS_85
        COLORS      = COLORS_45      if (ax_i == 0) else COLORS_85

        years_GCMs = numpy.arange(1990,2101,1,dtype=numpy.int)
        runoff_GCMs = numpy.zeros((4, len(years_GCMs)), dtype=numpy.float )

        for i,key in enumerate(GCM_KEYS):

            assert key in summary_dict.keys()

            color = COLORS[i]

            # If we're
            if type(key) == str:
                rows_mask = numpy.array([(row_value.find(key) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool)

                years_all = data_array['YEAR_END']
                mask_1990 = (years_all >= 1990)
                mask_total = rows_mask & mask_1990

                years = years_all[mask_total]
                runoff = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][mask_total]
                annual_runoff = data_array['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000'][mask_total]

            else:
                assert type(key) == tuple
                rows_mask = numpy.array([(row_value.find(key[0]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool) | \
                            numpy.array([(row_value.find(key[1]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool) | \
                            numpy.array([(row_value.find(key[2]) > -1) for row_value in data_array["DATASET"]], dtype=numpy.bool)

                hirham_years = data_array['YEAR_END'][rows_mask]
                hirham_annual_runoff = data_array['ICE_SLAB_RUNOFF_ANNUAL_GT_START_0000'][rows_mask]

                # Get the interpolated runoff values for the HIRHAM plots.
                years, annual_runoff, runoff = interpolate_HIRHAM_runoff_totals(hirham_years, hirham_annual_runoff)
                # Save this variable for later.
                hirham_cum_runoff = runoff
                print "HIRHAM_CUM_RUNOFF, {0}: {1}, {2}: {3}".format(years[years==2050], hirham_cum_runoff[years==2050],
                                                                     years[years==2100], hirham_cum_runoff[years==2100])


            if (not (type(key) == str) and (key[:5] == 'RACMO')):
                # Enter these values into the runoff calculations for below.
                runoff_GCMs[i,:] = runoff

            # Plot the lines.
            if type(key) == tuple:
                # For hirham, use dotted lines for interpolated values.
                ax.plot(years[years <= 2010], runoff[years <= 2010]/360., color=color,zorder=2, lw=1)
                ax.plot(years[(years > 2010) & (years <= 2040)], runoff[(years > 2010) & (years <= 2040)]/360.,
                                                      color=color, zorder=1, lw=1, linestyle="dotted")
                ax.plot(years[(years > 2040) & (years <= 2050)], runoff[(years > 2040) & (years <= 2050)]/360.,
                                                      color=color,zorder=2, lw=1)
                ax.plot(years[(years > 2050) & (years <= 2090)], runoff[(years > 2050) & (years <= 2090)]/360.,
                                                      color=color, zorder=1, lw=1, linestyle="dotted")
                ax.plot(years[(years > 2090) & (years <= 2100)], runoff[(years > 2090) & (years <= 2100)]/360.,
                                                      color=color,zorder=2, lw=1)
            else:
                # Skip RACMO dataset. Since the RACMO is the very last key in the RCP4.5 keys, this shouldn't screw up the indexing
                if (key[:5] != "RACMO"):
                    ax.plot(years, runoff/360., color=color,zorder=2, lw=1)

            # Get regression stats out.
            # YEARS 1990 - 2050
            if type(key) == tuple:
                # For hirham, use only non-interpolated values, 2000-2010 and 2041-2050
                mask_1990_2050 = ((years >= 2000) & (years <= 2010)) | ((years >= 2041) & (years <= 2050))
                years_1990_2050 = years[mask_1990_2050]
                runoff_1990_2050 = hirham_cum_runoff[mask_1990_2050]
                ann_runoff_1990_2050 = annual_runoff[mask_1990_2050]
            else:
                mask_1990_2050 = ((years >= 1990) & (years <= 2050))
                years_1990_2050 = years[mask_1990_2050]
                runoff_1990_2050 = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][rows_mask][mask_1990_2050]
                ann_runoff_1990_2050 = annual_runoff[mask_1990_2050]

            slope_1990_2050, intercept_1990_2050, rvalue_1990_2050, pvalue_1990_2050, stderr_1990_2050 = \
                scipy.stats.linregress(years_1990_2050, runoff_1990_2050)

            aslope_1990_2050, aintercept_1990_2050, arvalue_1990_2050, apvalue_1990_2050, astderr_1990_2050 = \
                scipy.stats.linregress(years_1990_2050, ann_runoff_1990_2050)

            data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n".format( \
                                        key if (type(key) == str) else key[-1],
                                        years_1990_2050[0],
                                        years_1990_2050[-1],
                                        "RUNOFF",
                                        slope_1990_2050,
                                        intercept_1990_2050,
                                        rvalue_1990_2050,
                                        pvalue_1990_2050,
                                        stderr_1990_2050,
                                        aslope_1990_2050,
                                        aintercept_1990_2050,
                                        arvalue_1990_2050,
                                        apvalue_1990_2050,
                                        astderr_1990_2050)

            summary_dict[key]['RUNOFF_TREND_GT_A'][0] = slope_1990_2050
            summary_dict[key]['RUNOFF_ACC_GT_A2' ][0] = aslope_1990_2050
            summary_dict[key]['RUNOFF_TOTAL_GT'  ][0] = runoff_1990_2050[-1]

            # YEARS 2050-2100
            if type(key) == tuple:
                # For hirham, use only non-interpolated values, 2000-2010 and 2041-2050
                mask_2050_2100 = ((years >= 2041) & (years <= 2050)) | ((years >= 2091) & (years <= 2100))
                years_2050_2100 = years[mask_2050_2100]
                runoff_2050_2100 = hirham_cum_runoff[mask_2050_2100]
                ann_runoff_2050_2100 = annual_runoff[mask_2050_2100]
            else:
                mask_2050_2100 = ((years >= 2051) & (years <= 2100))
                years_2050_2100 = years[mask_2050_2100]
                runoff_2050_2100 = data_array['ICE_SLAB_RUNOFF_TOTAL_GT_START_0000'][rows_mask][mask_2050_2100]
                ann_runoff_2050_2100 = annual_runoff[mask_2050_2100]

            slope_2050_2100, intercept_2050_2100, rvalue_2050_2100, pvalue_2050_2100, stderr_2050_2100 = \
                scipy.stats.linregress(years_2050_2100, runoff_2050_2100)

            aslope_2050_2100, aintercept_2050_2100, arvalue_2050_2100, apvalue_2050_2100, astderr_2050_2100 = \
                scipy.stats.linregress(years_2050_2100, ann_runoff_2050_2100)

            data_string = data_string + "{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13}\n".format( \
                                        key if (type(key) == str) else key[-1],
                                        years_2050_2100[0],
                                        years_2050_2100[-1],
                                        "RUNOFF",
                                        slope_2050_2100,
                                        intercept_2050_2100,
                                        rvalue_2050_2100,
                                        pvalue_2050_2100,
                                        stderr_2050_2100,
                                        aslope_2050_2100,
                                        aintercept_2050_2100,
                                        arvalue_2050_2100,
                                        apvalue_2050_2100,
                                        astderr_2050_2100)

            summary_dict[key]['RUNOFF_TREND_GT_A'][1] = slope_2050_2100
            summary_dict[key]['RUNOFF_ACC_GT_A2' ][1] = aslope_2050_2100
            summary_dict[key]['RUNOFF_TOTAL_GT'  ][1] = runoff_2050_2100[-1]


        # Horizontal grid lines every 10 mm SLE, light grey.
        for hline_height in (0,20,40,60):
            ax.axhline(hline_height, color='lightgrey', zorder=0, linestyle="dashed", linewidth=1)

        # Set the x-label years to true.
        [label.set_visible(True) for label in ax.get_xticklabels()]
        # Set the x-label years to true.
        [label.set_visible(True) for label in ax.get_yticklabels()]

        # Remove the top and right spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_ylabel("Runoff (mm SLE)")

    # Instead of explicitly sharing all the axes, just link the ones we need here at the end.
    axes_area[0].set_ylim(axes_area[1].get_ylim())
    axes_runoff[0].set_ylim(axes_runoff[1].get_ylim())

    fig.tight_layout()

    # Save the figure
    figname = "ICE_SLAB_AREAS_COMBINED.png"

    print "Plotting", figname
    fig.savefig(os.path.join(STATISTICS_OUTPUT_FOLDER,figname), dpi=2400)

    # Output the stats table summary data
    f = open(os.path.join(STATISTICS_OUTPUT_FOLDER, 'RUNOFF_STATS_SUMMARY.csv'), 'w')
    f.write('YEAR_START,YEAR_END,AREA_TREND_KM_A,AREA_KM2,RUNOFF_TREND_GT_A,RUNOFF_ACC_GT_A2,RUNOFF_TOTAL_GT\n')

    summary_keys = summary_dict.keys()
    summary_keys.sort()
    for key in summary_keys:
        for i in (0,1):
            f.write("{0},{1:d},{2:d},{3:0.1f},{4:0.3f},{5:0.3f},{6:0.3f},{7:0.2f}\n".format(
                    key if type(key) == str else key[-1],
                    summary_dict[key]['YEAR_START'][i],
                    summary_dict[key]['YEAR_END'][i],
                    summary_dict[key]['AREA_TREND_KM_A'][i],
                    summary_dict[key]['AREA_KM2'][i],
                    summary_dict[key]['RUNOFF_TREND_GT_A'][i],
                    summary_dict[key]['RUNOFF_ACC_GT_A2'][i],
                    summary_dict[key]['RUNOFF_TOTAL_GT'][i]
                  ))
    f.close()
    print 'RUNOFF_STATS_SUMMARY.csv written.'

    data_csv_fname = "ICE_SLAB_TRENDLINES_OUTPUT.csv"
    f = open(os.path.join(STATISTICS_OUTPUT_FOLDER, data_csv_fname), 'w')
    f.write(data_string)
    f.close()
    print "Exported", data_csv_fname


def run_decadal_RUNOFF_outputs_and_save_to_CSV_v3(skip_shapefiles=False):
    csv_fname = EMI_ICE_SLAB_DECADAL_OUTPUTS_CSV_FILE
    f = open(csv_fname,'w')
    for i,M in enumerate([HIRHAM5_Manager(), MAR_Manager(), RACMO_Manager()]):
        M.compute_ice_layers_and_excess_melt_runoff_V3(f,
                                                       forget_header=(False if i==0 else True),
                                                       skip_shapefiles=skip_shapefiles)

    f.close()
    print os.path.split(csv_fname)[-1], "written."
    return



def compute_area_of_ice_slab_polygons():
    shapefile_fname = r'C:/Users/mmacferrin/Dropbox/Research/Papers/MacFerrin et al 2016 - IceBridge Lenses/Data/IceBridge Area Shapefiles/IceBridgeArea_Shape.shp'

    # Create the data source and file.
    driver = ogr.GetDriverByName("ESRI Shapefile")

    data_source = driver.Open(shapefile_fname, 0)
    print "Opened", os.path.split(shapefile_fname)[-1]

    layer = data_source.GetLayer()
    features = [layer.GetFeature(i) for i in range(layer.GetFeatureCount())]

    polygon_areas = numpy.array([feature.geometry().Area() for feature in features])

    print polygon_areas
    poly_sum = numpy.sum(polygon_areas)
    print "SUM", poly_sum/1e6
    return polygon_areas


if __name__ == "__main__":
    run_decadal_RUNOFF_outputs_and_save_to_CSV_v3(skip_shapefiles=True)
