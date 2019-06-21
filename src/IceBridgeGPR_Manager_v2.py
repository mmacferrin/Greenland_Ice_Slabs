# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 12:20:43 2017

@author: mmacferrin
"""
#from FirnCore_Manager import FirnCore_Manager
from InSituGPR_Manager import InSituGPR_Manager, RadarSpeedPicker, InSituGPR_Track
from GPR_FileData import ICEBRIDGE_DATA_FOLDER, \
                         ICEBRIDGE_DATA_H5FILE, \
                         ICEBRIDGE_EXPORT_FOLDER, \
                         ICEBRIDGE_ICELENS_QUICKLOOK_FOLDER, \
                         ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE, \
                         ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE, \
                         ICEBRIDGE_SURFACE_INDICES_PICKLEFILE_FOLDER, \
                         ICEBRIDGE_SURFACE_SLICE_PICKLEFILE_FOLDER, \
                         ICEBRIDGE_SAMPLE_DEPTHS_PICKLEFILE_FOLDER, \
                         ICEBRIDGE_EXCLUSIONS_LAKES_OTHER_FILE, \
                         ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE, \
                         ICEBRIDGE_ROLL_CORRECTED_PICKLEFILE_FOLDER, \
                         ICEBRIDGE_ROLL_CORRECTION_OUTPUT_FILE, \
                         ICEBRIDGE_DEPTH_CORRECTED_PICKLEFILE_FOLDER, \
                         ICEBRIDGE_DEPTH_CORRECTION_OUTPUT_FILE, \
                         ACT13_DOWNSAMPLED_PICKLEFILE, \
                         ACT13_DOWNSAMPLED_COUNT_PICKLEFILE, \
                         ACT13_ICEBRIDGE_SUBSET__MASK_HMASK_VMASK_PICKLEFILE, \
                         ACT13_ICEBRDIGE__LAT_LON_E_N_DEPTHS_PICKLEFILE, \
                         ACT13_SUBSET__LAT_LON_ELEV_E_N_DEPTHS_PICKLEFILE, \
                         ICEBRIDGE_VALIDATION_OUTPUT_CSV_FILE, \
                         ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER, \
                         ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, \
                         ICEBRIDGE_SMOOTHED_ICE_LAYER_SHAPEFILE

import os
import numpy
import matplotlib.pyplot as plt
import tables
import scipy.io
import png
import re
import h5py
import pickle
import osgeo.ogr as ogr
import osgeo.osr as osr
import simplekml

# Define a quick guassian function to scale the cutoff mask above
def _gaussian(x,mu,sigma):
    return numpy.exp(-numpy.power((x-mu)/sigma, 2.)/2.)

class IceBridgeGPR_Manager_v2():
    '''IceBridgeGPR_Manager: Manages processes for all IceBridge Accumulation radar,
    including mapping.
    Uses: FirnCore_Manager, FirnCore_Profile, InSituGPR_Manager, InSituGPR_Track, IceBridgeGPR_Track
    Used by: MasterPlanFunction
    '''
    def __init__(self, h5filename = ICEBRIDGE_DATA_H5FILE, verbose=True):
        self.track_names = None
        self.tracks = None
        self.verbose = verbose
        self.h5filename = h5filename
        self.h5file = tables.open_file(h5filename, 'r')

    def get_2013_IceBridge_reference_track(self):
        '''Retrieve the 2013 IceBridge track that serves as a "reference" against the
        ACT in situ GPR.  This will be used to do all the tuning.'''
        ACT13_FILE_ID = "20130409_01_010_012"
        if self.tracks is None:
            return IceBridgeGPR_Track_v2(self.h5file, ACT13_FILE_ID)
        else:
            return self.tracks[self.track_names.index(ACT13_FILE_ID)]

    def delete_exported_images(self, regex_filter=None):
        '''Utility for cleaning up files.  If regex_filter is None, delete them all.  Otherwise use the filter to find matches.'''
        export_folder = ICEBRIDGE_EXPORT_FOLDER
        fnames_all = [os.path.join(export_folder, f) for f in os.listdir(export_folder)]

        for fname in fnames_all:
            if regex_filter is None or re.search(regex_filter, fname) != None:
                if self.verbose:
                    print "Deleting", os.path.split(fname)[-1]
                os.remove(fname)

        return

    def export_images_to_subfolder(self, regex_filter):
        '''Take all the images/files in a folder that fit a particular regular expression.
        Move them all to a subfolder of that same expression.'''
        base_folder_name = ICEBRIDGE_EXPORT_FOLDER
        file_list = os.listdir(base_folder_name)
        # Ignore folders
        file_list = [f for f in file_list if ((not os.path.isdir(os.path.join(base_folder_name,f))) and (re.search(regex_filter, f) != None))]
        subfolder_path = os.path.join(base_folder_name, regex_filter)
        if not (os.path.exists(subfolder_path)):
            os.mkdir(subfolder_path)

        file_matches_old = [os.path.join(base_folder_name, f) for f in file_list]
        file_matches_new = [os.path.join(subfolder_path, f) for f in file_list]
        for f1, f2 in zip(file_matches_old, file_matches_new):
            print "Moving {0} to {1}".format(os.path.split(f1)[-1], os.path.split(subfolder_path)[-1]),
            # IF the file exists already, overwrite it.
            if os.path.exists(f2):
                print "(Deleting previous)"
                os.remove(f2)
            else:
                print
            # This moves the actual file.
            os.rename(f1, f2)

    def compile_icebridge_tracks_with_ice_lenses(self, quicklook_directory = ICEBRIDGE_ICELENS_QUICKLOOK_FOLDER):
        '''Similar to "::compile_icebridge_tracks()", this will complile a dictionary of IceBridge_Track objects.
        However, instead of using all tracks and all files within the tracks, this will peruse a directory containing only files that
        have been flagged to be part of IceBridge files that may (or do) contain ice lenses in them.

        The directory contains only the quicklook PDF files that visualize the tracks.
        '''
        track_count_N = 0

        # 1) Open the directory, get a list of all the files in there
        # 2) Peruse files, separate out "map" files from "echo" files, use just the echo files.
        all_files = [f for f in os.listdir(quicklook_directory) if f.find("1echo.jpg") != -1]

        # 3) Compile lists of all sequential adjoining files (part of the same sequence) that will create tracks.
        track_list_dict = {}
        current_list = []
        for f in all_files:
            # If we have a new list, OR if the last file listed in the current list has the
            # same FILE_ID, and the file number is one greater than the last one, append it onto the list.
            if len(current_list) == 0 or \
                (f[0:11] == current_list[-1][0:11] and (int(f[12:15]) - int(current_list[-1][12:15])) == 1):

                current_list.append(f)
            else:
                # Convert the track_id substring into a TRACK_ID integer, append it.
                track_id = current_list[0][0:11] + current_list[0][11:15] + current_list[-1][11:15]
                # Stick it into the dictionary
                track_list_dict[track_id] = current_list

                # Create a new "current list" with the new file in it.
                current_list = [f]
                track_count_N += 1

        # Make sure the very last one gets on there.
        track_list_dict[current_list[0][0:11] + current_list[0][11:15] + current_list[-1][11:15]] = current_list
        track_count_N += 1

        if self.verbose:
            print track_count_N, "tracks from", len(all_files), "files.\n"


        # 4) Create a (subset) track from each file.
        self.track_names = track_list_dict.keys()
        self.track_names.sort()

        tracks = [IceBridgeGPR_Track_v2(self.h5file, name, verbose=self.verbose) for name in self.track_names]

        # 5) Save the dictionary.
        self.tracks = tracks

        return tracks

    def compute_number_of_traces_omitted_by_surface_mismatch(self):
        TOTAL_numtraces = 0
        TOTAL_excluded_at_start = 0
        TOTAL_excluded_by_bad_surfacepicks = 0

        for track in self.tracks:
            track_numtraces = track.numtraces()
            track_excluded_at_start = numpy.count_nonzero(
                                        ~track.mask_manager.compute_mask(track.NAME,
                                                                         ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE,
                                                                         track_numtraces)
                                                          )
            track_excluded_by_bad_surfacepicks = numpy.count_nonzero(
                                        ~track.mask_manager.compute_mask(track.NAME,
                                                                         ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE,
                                                                         track_numtraces)
                                                                     )
            TOTAL_numtraces += track_numtraces
            TOTAL_excluded_at_start += track_excluded_at_start
            TOTAL_excluded_by_bad_surfacepicks += track_excluded_by_bad_surfacepicks

            print track.NAME, track_numtraces, track_excluded_at_start, track_excluded_by_bad_surfacepicks

        print TOTAL_numtraces, "traces over", len(self.tracks), "flight lines."
        print "{0} excluded from start, {1:0.5f}%".format(TOTAL_excluded_at_start, float(TOTAL_excluded_at_start)/TOTAL_numtraces*100)
        print "{0} of remaining excluded from poor surface picks, {1:0.5f}%".format(TOTAL_excluded_by_bad_surfacepicks, float(TOTAL_excluded_by_bad_surfacepicks)/(TOTAL_numtraces-TOTAL_excluded_at_start)*100)

        return

    def compute_total_distance(self):

        total_distance = 0
        for t in self.tracks:
            total_distance += numpy.sum(t.compute_distances())

        print "\nTOTAL DISTANCE: {0:f} km".format(total_distance / 1000.0)

        return total_distance

    def plot_correlation_between_roll_and_curvature(self):
        '''Go through all the tracks, get all those that have both roll and curvature, and compute a scatterplot of the correlation between them.'''

        roll_curve_picklefile = r'F:\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\accum\exported\Roll_Curve_Picklefiles\Roll_Curve_TOTAL.pickle'
        if os.path.exists(roll_curve_picklefile):
            f = open(roll_curve_picklefile, 'r')
            all_roll, all_curve = pickle.load(f)
            f.close()
        else:

            total_roll_list = []
            total_curve_list = []
            # Collect the roll and curvature for all tracks which have both.
            for track in self.tracks:
                print track.NAME
                roll = track._return_aircraft_roll()
                if roll is None:
                    print "NO ROLL, MOVING ON."
                    continue
                curvature = track._compute_path_curvature(normalize=True)
                # Every once in awhile one of the curvature calculations give a NaN. Mask it out and continue.
                curve_nan_mask = numpy.isnan(curvature)
                num_nans = numpy.count_nonzero(curve_nan_mask)
                if num_nans > 0:
                    print num_nans, "of", len(curvature), "curvatures are NaN.  Eliminating..."
                    # Mask out the ones that are NaNs
                # The curvature should have 2 less elements than the roll.  Do a sanity check for this.
                assert (len(roll) - 2) == len(curvature)
                # Cut the first and last roll point.
                roll = roll[1:-1]
                # Maks out anywhere the curvatures were nans
                curvature = curvature[~curve_nan_mask]
                roll = roll[~curve_nan_mask]

                total_roll_list.append(roll)
                total_curve_list.append(curvature)
                print

            # Now flatten the lists out.
            N = numpy.sum([r.shape[0] for r in total_roll_list])

            # Create empty vectors to deposit them.
            all_roll  = numpy.empty((N,), dtype=total_roll_list[0].dtype)
            all_curve = numpy.empty((N,), dtype=total_curve_list[0].dtype)
            list_index = 0
            for r,c in zip(total_roll_list, total_curve_list):
                assert r.shape[0] == c.shape[0]
                all_roll[list_index:list_index+r.shape[0]] = r
                all_curve[list_index:list_index+c.shape[0]] = c
                list_index += r.shape[0]

            # Save to picklefile
            f = open(roll_curve_picklefile, 'w')
            pickle.dump((all_roll, all_curve), f)
            f.close()
            print "Exported", os.path.split(roll_curve_picklefile)[-1]

            print N, "points."

        roll_deg = numpy.degrees(all_roll)
        # We need to use the absolute value of the roll here.
        roll_abs = numpy.abs(roll_deg)
        curve_deg = numpy.degrees(all_curve)

        # Compute curve fit.
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(roll_abs, curve_deg)
        print "SLOPE {0}\nINT {1:f}\nR: {2:f}, R^2: {3:f}\nP: {4:f}\nSTD_ERR: {5:f}\n".format(slope,intercept,r_value,r_value**2,p_value,std_err)
        xmin = numpy.min(roll_abs)
        xmax = numpy.max(roll_abs)
        x = numpy.arange(xmin, xmax, (xmax-xmin)/100.0)
        y = slope*x + intercept

        # Plot every 50th point for now.  (See what we have to do after that.)
        plt.scatter(roll_abs[::20], curve_deg[::20], marker='.', color='blue')
        # Plot the linear regression through the points.
        plt.plot(x,y,color='darkred', linewidth=3)
        plt.text(1,0.23,"y = 0.0099x - 0.0013", size="large", weight="heavy")
        plt.text(1.25,0.18,"p $\leq$ 1e-7\nR$^2$ = 0.027\nR$^2$ = 0.85 (Roll$>$5$^\circ$)", size="large")

        plt.xlabel("Aircraft Roll (unsigned) ($^\circ$)")
        plt.ylabel("Flight Path Curvature ($^\circ$)")
        plt.xlim(-0.25, xmax*1.01)
        plt.ylim(-.004, numpy.max(curve_deg[::20])*1.01)
        plt.show()
        plt.savefig(r"C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\Roll_Curve_Correlation.png",dpi=600)
        print "Exported Roll_Curve_Correlation.png"

        # Now get a table here for when the correlation breaks down.
        # Print this stuff in a CSV format for easy copy-paste-analyze
        cutoffs = numpy.arange(0,xmax,0.5)
        print
        print "ROLL_CUTOFF,P_BELOW,R_BELOW,R2_BELOW,P_ABOVE,R_ABOVE,R2_ABOVE"
        for cutoff in cutoffs:
            roll_below_mask = roll_abs < cutoff
            xbelow = roll_abs[roll_below_mask]
            xabove = roll_abs[~roll_below_mask]
            ybelow = curve_deg[roll_below_mask]
            yabove = curve_deg[~roll_below_mask]

            sb, ib, rb, pb, sb = scipy.stats.linregress(xbelow, ybelow)
            sa, ia, ra, pa, sa = scipy.stats.linregress(xabove, yabove)

            print "{0},{1},{2},{3},{4},{5},{6}".format(cutoff, pb, rb, rb**2, pa, ra, ra**2)

        print

    def perform_roll_corrections(self):
        '''Perform the roll corrections for all the tracks, and export the correction parameters to a text file.'''
        fout = open(ICEBRIDGE_ROLL_CORRECTION_OUTPUT_FILE, 'w')
        fout.write("NAME,A_R,A_S,A_p_value,A_r_value,C_T,C_U,C_V,C_p_value,C_r_value,Roll_or_Curve,Max_Roll,A_avg_20m,C_avg_20m,mean_correction_20m_below_5deg,mean_correction_20m_above_5deg\n")

        # Loop through all the tracks.
        for track in self.tracks:
            print track.NAME
            # Perform the roll correction, output picklefiles
            pout = track.perform_roll_correction(export=True)
            fout.write("{0},{1:0.8f},{2:0.8f},{3:0.8f},{4:0.8f},{5:0.8f},{6:0.8f},{7:0.8f},{8:0.8f},{9:0.8f},{10},{11:0.4f},{12:0.8f},{13:0.8f},{14:0.8f},{15:0.8f}\n".format(track.NAME, *pout))

        fout.close()
        print "Exported", os.path.split(ICEBRIDGE_ROLL_CORRECTION_OUTPUT_FILE)[-1]

    def perform_depth_corrections(self, export=True):
        '''Perform depth corrections for all the tracks, and export the correction parameters to a text file.'''
        fout = open(ICEBRIDGE_DEPTH_CORRECTION_OUTPUT_FILE, 'w')
        fout.write("NAME,A,B,C\n")

        for track in self.tracks:
            print track.NAME
            pout = track.perform_depth_correction(export=export)
            fout.write("{0},{1},{2},{3}\n".format(track.NAME, *pout))

        fout.close()
        print "Exported", os.path.split(ICEBRIDGE_DEPTH_CORRECTION_OUTPUT_FILE)[-1]

    def export_ice_layer_lat_lon_distance_thicknesses(self):
        '''Export to a CSV, ice layer lat,lon,& thicknesses.  Omit all zero values.'''
        tracks = self.compile_icebridge_tracks_with_ice_lenses()

        fout = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'w')
        header = "Track_name,Tracenumber,lat,lon,alongtrack_distance_m,20m_ice_content_m\n"
        fout.write(header)

        for track in tracks:
            print track.NAME,
            lats, lons, distances, ice_contents = track.return_ice_layers_lat_lon_distance_thickness(masked=False)
            # The one "really long" track has artifacts in the center that aren't real ice layers.  Filter these out.
            if track.NAME == "20120412_01_095_095":
                ice_contents[0:9000] = 0.0

            assert len(lats) == len(lons) == len(ice_contents)
            tracenums = numpy.arange(len(lats), dtype=numpy.int)

            tracecount = 0
            for lat, lon, tracenum, distance, ice_content in zip(lats, lons, tracenums, distances, ice_contents):
                # Record ONLY traces that have > 1 m ice content in them.  We're not interested in thinner stuff here.
                # If it has > 16 m ice content (80%), we also omit it, just to keep pure ice out of it.
                if 1.0 <= ice_content <= 16.0:
                    line = "{0},{1},{2},{3},{4},{5}\n".format(track.NAME, tracenum, lat, lon, distance, ice_content)
                    fout.write(line)
                    tracecount += 1
            print tracecount, "of", len(lats), "traces."
            print

        fout.close()
        print "Exported", os.path.split(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE)[-1]


    def import_ice_layer_lat_lon_distance_thicknesses(self):
        '''Import from a CSV, return the track_name, tracenum, lat, lon, cumulative_distance, and ice layer thickness at each trace.'''
        print "Reading", os.path.split(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE)[-1]
        fin = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'r')
        lines = [line.strip().split(",") for line in fin.readlines()]
        header = dict([(name,i) for i,name in enumerate(lines[0])])
        lines = lines[1:]
        track_names = numpy.array([items[header["Track_name"]] for items in lines])
        tracenums = numpy.array([int(items[header["Tracenumber"]]) for items in lines])
        lats = numpy.array([float(items[header["lat"]]) for items in lines])
        lons = numpy.array([float(items[header["lon"]]) for items in lines])
        distances = numpy.array([float(items[header["alongtrack_distance_m"]]) for items in lines])
        ice_content_m = numpy.array([float(items[header["20m_ice_content_m"]]) for items in lines])

        return track_names, tracenums, lats, lons, distances, ice_content_m

    def export_smoothed_ice_layer_shapefile(self):
        '''Take the CSV ice layer data and create a shapefile from it, with smoothed data averaged every 2.5 (?) km.'''

        fname = ICEBRIDGE_SMOOTHED_ICE_LAYER_SHAPEFILE

        # Create the data source and file.
        driver = ogr.GetDriverByName("ESRI Shapefile")

        if os.path.exists(fname):
            print "Deleting previous", os.path.split(fname)[-1]
            driver.DeleteDataSource(fname)

        data_source = driver.CreateDataSource(fname)

        # Create a WGS84 geographic spatial reference to go with it.
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)

        # Create the layer
        layer = data_source.CreateLayer("Ice Lenses", srs, geom_type=ogr.wkbLineString)

        # Add the variable field we're interested in
        field_name = ogr.FieldDefn("TrackName", ogr.OFTString)
        field_name.SetWidth(19)
        layer.CreateField(field_name)

        field_mean = ogr.FieldDefn("Mean_Ice_m", ogr.OFTReal)
        layer.CreateField(field_mean)

        field_std = ogr.FieldDefn("Std_Ice_m", ogr.OFTReal)
        layer.CreateField(field_std)

        field_count = ogr.FieldDefn("NumTraces", ogr.OFTInteger)
        layer.CreateField(field_count)

        # Read the CSV data.
        track_names, tracenums, lats, lons, distances, ice_content_m = self.import_ice_layer_lat_lon_distance_thicknesses()

        for t_name in numpy.unique(track_names):
            print t_name,
            # Subset the data for just this track.
            t_mask = (track_names == t_name)
            t_lats = lats[t_mask]
            t_lons = lons[t_mask]
            t_distances = distances[t_mask]
            t_ice_content_m = ice_content_m[t_mask]

            STEP_m = 2500 # Meters per step (2.5 km here)

            step_count = 0
            for d in numpy.arange(0,numpy.max(t_distances), STEP_m):
                step_mask = (d <= t_distances) & (t_distances < (d + STEP_m))
                # SKip if there are no points in this step.
                if not numpy.any(step_mask):
                    continue

                mean_ice_content = numpy.mean(t_ice_content_m[step_mask])
                std_ice_content = numpy.std(t_ice_content_m[step_mask])
                step_numtraces = numpy.count_nonzero(step_mask)
                step_lats = t_lats[step_mask]
                step_lons = t_lons[step_mask]

                # Create the feature:
                feature = ogr.Feature(layer.GetLayerDefn())
                # Set the attributes
                feature.SetField("TrackName", t_name)
                feature.SetField("Mean_Ice_m", mean_ice_content)
                feature.SetField("Std_Ice_m", std_ice_content)
                feature.SetField("NumTraces", step_numtraces)

                wkt_geom_string = "LINESTRING (" + \
                                  ",".join(["{0} {1}".format(lon, lat) for (lat, lon) in zip(step_lats, step_lons)]) \
                                  + ")"

                linestring = ogr.CreateGeometryFromWkt(wkt_geom_string)
                feature.SetGeometry(linestring)
                # Put this feature in the layer
                layer.CreateFeature(feature)
                # Dereference the feature
                step_count += 1
                feature = None

            print step_count, "steps."

        # Now dereference the data source to save it.
        data_source = None
        print "Exported", os.path.split(fname)[-1]
        return

    def export_KML_reference_tracks(self):
        '''There are three (four?) reference tracks we're using to compare with LandSat imagery.
        Load them here, create a KML with those tracks in them, just where ice slabs show up.'''

        csv_file = open(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE, 'r')
        print "Reading", os.path.split(ICEBRDIGE_ICE_LAYER_OUTPUT_CSV_FILE)[-1]
        lines = csv_file.readlines()
        csv_file.close()

        header = dict([(name.strip(), i) for i,name in enumerate(lines[0].split(','))])
        lines = [line.strip().split(',') for line in lines[1:]]

        REFTRACKS_LIST = ['20110502_01_079_083',
                          '20100526_01_077_078',
                          '20130409_01_010_012',
                          '20110418_01_173_175']


        kml = simplekml.Kml()

        output_folder = r'C:\Users\mmacferrin\Dropbox\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\Greenland Ice Lens Transects'
        for track_name in REFTRACKS_LIST:
            print track_name

            # Read the points for that track from the CSV
            lines_subset = [line for line in lines if (line[header["Track_name"]] == track_name)]
            print len(lines_subset)
            lats = [float(line[header["lat"]]) for line in lines_subset]
            lons = [float(line[header["lon"]]) for line in lines_subset]

            print len(lats), len(lons)

            kml.newlinestring(name = track_name, description= "Ice Slabs in Track " + track_name,
                              coords = zip(lons, lats))

        kml_filename = os.path.join(output_folder, "Output_Ice_Slabs.kml")
        kml.save(kml_filename)
        print os.path.split(kml_filename)[-1], "written."


class Mask_Manager():
    '''A class for managing mask files for IceBridgeGPR_Track_v2 instances.
    This is a helper class so that each icebridgeGPR instance doesn't have to open, parse and
    manage all the files themselves.  This becomes a class variable (not an instance variable) of the
    IceBridgeGPR_Track_v2() class.'''

    def __init__(self):
        ###############################################
        ## FILES/VARIABLES FOR MASKING OUT PIXELS
        self.mask_dictionary = None
        self.valid_mask_filenames = [ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE, ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE, ICEBRIDGE_EXCLUSIONS_LAKES_OTHER_FILE]
        self._initialize_mask_dictionary()
        ###############################################

        ###############################################
        ## FILES/VARIABLES FOR SUGGESTING STARTING SURFACE PIXELS
        self.starting_surface_pixel_dictionary = None
        self._initialize_surface_pixel_dictionary()
        ###############################################

    def _initialize_surface_pixel_dictionary(self):
        '''Reads the ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE, and gets a dictionary of
        flight lines whose surface picks were mis-picked by the auto-picker, and uses
        an improved starting location instead.  Most of these will be None.'''

        sdict = {}
        f = open(ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE,'r')
        lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
        f.close()

        for line in lines:
            items = line.split()
            if len(items) == 1:
                sdict[items[0]] = None
            else:
                assert len(items) == 2
                sdict[items[0]] = int(items[1])

        self.starting_surface_pixel_dictionary = sdict

    def return_suggested_starting_pixel(self, flightline_name):
        return self.starting_surface_pixel_dictionary[flightline_name]

    def _initialize_mask_dictionary(self):
        '''Reads the mask files, populates the dictionary.
        The dictionary is double-indexed as such:
        mask_dictionary[mask_filename][IceBridge_track_name] = list-of-masks
            list-of-masks = None if nothing masked.
                          = ((None,##),(##,##),(##,None)) <-- any combination of these for masks.
        '''
        mdict = {}
        # Iterate through all the files
        for i,fname in enumerate(self.valid_mask_filenames):
            # Open the file
            f = open(fname, 'r')
            # Populate the top level of the dictionary
            mdict[fname] = {}

            # Iterate flightlines
            lines = [line.strip() for line in f.readlines() if len(line.strip()) > 0]
            for line in lines:
                items = line.split()
                # The flightline name will always be the first item
                linename = items[0]
                assert len(linename) == 19
                if len(items) == 1:
                    # If nothing appears after the line name, there are no exclusions on this file.
                    mdict[fname][linename] = None
                else:
                    # If there are exclusions listed, add them to our list.
                    gaps = []
                    for j,gaptxt in enumerate(items[1:]):
                        if gaptxt[0] == "-":
                            # -NNN , always at the start
                            assert j==0
                            gaps.append((None,int(gaptxt[1:])))
                        elif gaptxt[-1] == "-":
                            # NNN- , always at the end
                            assert (j+2)==len(items)
                            gaps.append((int(gaptxt[:-1]),None))
                        elif gaptxt.find("-") > -1:
                            # NNN-NNN , somewhere in the middle
                            gap = tuple([int(g) for g in gaptxt.split('-')])
                            assert len(gap) == 2
                            gaps.append(gap)
                        else:
                            gap = (int(gaptxt),)
                            gaps.append(gap)

                    # Add the list of gaps to our dictionary.
                    mdict[fname][linename] = gaps

        self.mask_dictionary = mdict
        return

    def compute_mask(self, flightline_name, mask_filenames, array_length):
        '''Create a boolean mask of values included in this particular mask iteration.
        Valid options for "mask_filenames"
        - "ALL": use them all combined.
        - Any complete filenames in self.valid_mask_filenames
        - A list or tuple of filenames in self.valid_mask_filenames'''

        if type(mask_filenames) == str:
            assert mask_filenames == "ALL" or mask_filenames in self.valid_mask_filenames
        else:
            assert type(mask_filenames) in (list,tuple)
            assert numpy.all([(fname in self.valid_mask_filenames) for fname in mask_filenames])

        boolean_array = numpy.ones((array_length,), dtype=numpy.bool)

        # Iterate over just one file, or all of them.
        if type(mask_filenames) == str:
            if mask_filenames == "ALL":
                filenames_to_iterate = self.valid_mask_filenames
            else:
                filenames_to_iterate = (mask_filenames,)
        elif type(mask_filenames) in (list,tuple):
            filenames_to_iterate = mask_filenames

        # Iterate over each of the files
        for fname in filenames_to_iterate:
            # Find the gaps from that file.
            gaps = self.mask_dictionary[fname][flightline_name]

            if gaps is not None:
                for gap in gaps:
                    # Subtract the gap from the overall array.
                    # Some of these gaps might overlap, that's okay.
                    if len(gap) == 2:
                        boolean_array[gap[0]:(gap[1]+1 if gap[1] is not None else None)] = False
                    else:
                        assert len(gap) == 1
                        boolean_array[gap[0]] = False

        return boolean_array

class IceBridgeGPR_Track_v2():
    '''IceBridgeGPR_Track2: Stores information and processing on a single
       IceBridge Accumulation Radar flight line. Version 2
    Uses: InSituGPR_Manager, FirnCoreManager
    Used by: IceBridgeGPR_Manager
    '''
    # Determine a universal radar speed to use for the IceBridge tracks.
    # Use the density of bubbly refrozen ice from Machguth, et al (2016), of 873+/-25 kg m-3.
    # The ideal Robin's speed factors was determined for the in situ GPR in InSituGPR_Manager.
    # (However, we can revisit this if it proves to need adjusting)
    # Define this as a class variable rather than an object variable, so we only have to define it once.
    radar_speed_m_s = RadarSpeedPicker().modified_Robin_speed(0.873, 0.734)
    mask_manager = Mask_Manager()

    def __init__(self, h5file, flightline_ID, verbose=True):
        ##########################################################
        # Attributes -- convention is variables start with ALL_CAPS words, functions don't.
        self.H5FILE = h5file
        self.NAME = str(flightline_ID)

        # Sanity check on the name:
        # YYYYMMDD_XX_YYY_ZZZ format. XX = Track sub-number, YYY = starting fileID, ZZZ = ending fileID
        assert len(self.NAME) == 19

        self.FNAME_surface_picks_picklefile = None # Define here?
        self.FNAME_roll_corrected_picklefile = os.path.join(ICEBRIDGE_ROLL_CORRECTED_PICKLEFILE_FOLDER, self.NAME + "_ROLL_CORRECTED.pickle")
        self.FNAME_depth_corrected_picklefile = os.path.join(ICEBRIDGE_DEPTH_CORRECTED_PICKLEFILE_FOLDER, self.NAME + "_DEPTH_CORRECTED.pickle")
        self.FNAME_ice_lenses_picklefile = os.path.join(ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER, self.NAME + "_SG1_CUTOFF_-0.45_THRESHOLD_350.pickle")

        self.TRACES = None
        self.SAMPLE_TIMES = None
        self.SAMPLE_DEPTHS = None
        self.TRACES_surface_slice_100m = None
        self.TRACES_roll_corrected = None
        self.TRACES_depth_corrected = None
        self.TRACES_boolean_ice_layers = None

        self.YEAR = None
        self.FILENAMES = None
        self.VERBOSE = verbose
        self.TABLE_coords_table = None
        self.TABLE_file_table = None
        self.LIST_filenames = None
        self.LIST_file_ids = None
        self.LIST_original_surface_picks = None
        self.LIST_original_surface_indices = None
        self.LIST_improved_surface_indices = None
        ##########################################################

    def DO_IT_ALL(self):
        '''From start to finish, process all the data for this track.  The "do it all" function.
        This will rarely (if ever) be called from the start, but is helpful for reproduction as well as
        for keeping track of the processing steps needed to transform the entire file.'''
        self._read_metadata()
        self.compute_surface_picks(export=True)
        self.perform_roll_correction(export=True, max_depth_m=100)
        self.perform_depth_correction(export=True, max_depth_m=100)
        self.identify_ice_lenses(export=True, max_depth_m=20)

    def _read_metadata(self):
        '''Read all the metadata from the self.h5file IceBridge metadata file.  Get the needed data for this set.'''
        if isinstance(self.H5FILE, tables.File):
            f = self.H5FILE
        elif isinstance(self.H5FILE, str):
            f = tables.open_file(self.H5FILE, 'r')
            self.H5FILE = f
        else:
            raise TypeError("Unknown h5 object type: " + str(type(self.H5FILE)))

        flights_table = f.get_node('/Accumulation/Flight_Lines')
        coords_table = f.get_node('/Accumulation/Coordinates')
        file_table = f.get_node('/Accumulation/File_Names')

        # Extract the year of this flight.
        flightline_id = (self.NAME[0:8] + self.NAME[9:11])
        self.YEAR = flights_table.read_where('Flightline_ID == {0}'.format(flightline_id))['year'][0]

        # Get a list of all the subset files we need to read for the data.
        file_table_temp = file_table.read_where('Flightline_ID == {0}'.format(flightline_id))
        filenames_temp = file_table_temp['relative_filepath']
        file_ids_temp = file_table_temp['File_ID']
        # Retreive only files that are in the "flightline_id_start_end" range we specified
        relative_filenum_start = int(self.NAME[-7:-4])
        relative_filenum_end   = int(self.NAME[-3:])
        # Get the files that span between the starting and ending identified file names
        filenames_subset = [(i,fname) for (i,fname) in enumerate(filenames_temp) if relative_filenum_start <= int(fname[-7:-4]) <= relative_filenum_end]
        file_ids_subset = [item[0] for item in filenames_subset]
        filenames_subset = [item[1] for item in filenames_subset]
        # Get the file IDs for the start and end file ID in the table.
        file_ids_start, file_ids_end = file_ids_temp[file_ids_subset[0]], file_ids_temp[file_ids_subset[-1]]

        self.TABLE_file_table = file_table.read_where(\
                          '(Flightline_ID == {0}) & (File_ID >= {1}) & (File_ID <= {2})'.format( \
                            flightline_id, file_ids_start, file_ids_end))

        relative_paths = self.TABLE_file_table['relative_filepath']
        self.FILENAMES = [os.path.join(ICEBRIDGE_DATA_FOLDER, rp) for rp in relative_paths]
        # Ensure that the file names are in alpha-numerical order (they should already be, so this is likely redundant, but still.)
        self.FILENAMES.sort()

        # A subset of the 'Accumulation/Coordinates' table, referring only to
        # records from this particular flight track.
        if len(self.NAME) == 10:
            self.TABLE_coords_table = coords_table.read_where('Flightline_ID == {0}'.format(flightline_id))
        else:
            self.TABLE_coords_table = coords_table.read_where( \
                                    '(Flightline_ID == {0}) & (File_ID >= {1}) & (File_ID <= {2})'.format( \
                                     flightline_id, file_ids_start, file_ids_end))

        return

    def _convert_h5_dataset_to_dict(self, dataset):
        '''The files read as h5py files don't return arrays when you index them.
        For this, we'll convert the h5py "Dataset" into a simple key/value dictionary.'''
        if type(dataset) == dict:
            return

        assert type(dataset) in (h5py.Dataset, h5py.File)

        data_dict = {}

        for key in dataset.keys():
            # Only include "Data" field if it's explicity included.
            # Leave out metadata fields, which are categorized in h5py.Group objects
            if (key != "Data") and (type(dataset[key]) != h5py.Group):
                data_dict[key] = dataset[key].value.flatten()
            elif (key == "Data"):
                # The data field in the 2014 H5PY files are transposed with respect to the rest of the data sets.
                data_dict[key] = dataset[key].value.transpose()

        return data_dict

    def _import_from_picklefile(self, FNAME):
        '''Import a picklefile and return the object within.'''
        f = open(FNAME, "r")
        print "Reading", os.path.split(FNAME)[-1]
        data = pickle.load(f)
        f.close()
        return data

    def _export_to_picklefile(self, data, FNAME):
        '''Export a variable to a picklefile.'''
        f = open(FNAME, "w")
        pickle.dump(data, f)
        f.close()
        print "Exported", os.path.split(FNAME)[-1]
        return

    def _compute_score_accuracy_wrt_reference_track_BOOLEAN(self, boolean_array, reference_track):
        '''Given a particular boolean array and the downsampled BOOLEAN reference track,
        return two scores:
        1) The Type-1 Error rate (false positives)
        2) The Type-2 Error rate (false negatives)'''
        N = float(boolean_array.size)
        type_1_errors = float(numpy.count_nonzero(boolean_array & ~reference_track)) / N
        type_2_errors = float(numpy.count_nonzero(~boolean_array & reference_track)) / N

        return type_1_errors, type_2_errors


    def _compute_score_accuracy_wrt_reference_track(self, boolean_array, GPR_resampled):
        '''Given a particular boolean array and the downsampled reference track,
        return two scores:
        1) The Type-1 Error rate (false positives)
        2) The Type-2 Error rate (false negatives)'''
        type_1_errors = (1.0 - numpy.mean(GPR_resampled[ boolean_array])) if (numpy.count_nonzero(boolean_array)  > 0) else 0.0
        type_2_errors =        numpy.mean(GPR_resampled[~boolean_array])  if (numpy.count_nonzero(~boolean_array) > 0) else 0.0

        return type_1_errors, type_2_errors

    def plot_in_situ_and_reference_track_together(self):
        '''Create a plot of the in-situ data, atop the OIB reference track data.'''
        # Make sure only the reference track is used.  No other track should be calling this.
        # (Ideally I'd make the reference track a sub-class of the general OIB track, but this is easier for now and it works the same.)
        assert self.NAME == "20130409_01_010_012"

        # Get In Situ GPR object
        in_situ_manager = InSituGPR_Manager()
        in_situ_track = in_situ_manager.GPR_tracks[in_situ_manager.track_names.index("ACT_Transect")]
        is_traces, is_cmap = in_situ_track.export_boolean_ice_array(export_image=False)

        # Subset in situ traces to 45 km.
        distances = numpy.arange(is_traces.shape[1], dtype=numpy.float)*1.5 / 1000
        _, is_lons, _ = in_situ_track.extract_lats_lons_elevs()
        depths = in_situ_track.calculate_trace_depths()
        N = numpy.count_nonzero(distances <= 45.0)
        M = numpy.count_nonzero(depths <= 15.0)
        print is_traces.shape
        is_traces = is_traces[:M,:N]
        is_lons = is_lons[:N]
        print is_traces.shape
        print is_cmap

        # Get IceBridge Traces
        traces = self.get_processed_traces(datatype="depth_corrected")
        # Retreive the depths, subset to max_depth_m (20 m)
        depths = self.get_sample_depths(trace_array = traces)
        depth_start_N = numpy.nonzero(depths >= 1.0)[0][0]
        depth_end_N = numpy.nonzero(~(depths <= 26.0))[0][0]
        # Rather than distances, just subset this to the longitudes of the in-situ radar above it.
        _, lons = self.return_coordinates_lat_lon()
        print lons
        lons_mask = (lons >= numpy.min(is_lons)) & (lons <= numpy.max(is_lons))
        lons_start_N = numpy.nonzero(lons_mask)[0][0]
        lons_stop_N = numpy.nonzero(~lons_mask[lons_start_N:])[0][0] + lons_start_N + 1
        print lons_start_N, lons_stop_N
        print traces.shape
        traces = traces[depth_start_N:depth_end_N,lons_start_N:lons_stop_N]
        # Trip top & bottom 1% for display purposes (contrast enhancement)
        low_cutoff, high_cutoff = numpy.percentile(traces, [1,99])
        traces[traces < low_cutoff] = low_cutoff
        traces[traces > high_cutoff] = high_cutoff
        print traces.shape

        # FIGURE PLOTTING
        # 8" wide by 4" tall
        figsize=(8.0,2.5)
        # One figure, three axes.
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, dpi=600, figsize=figsize)
        # Take care of spacing and font sizes
        fig.subplots_adjust(bottom=0.25, hspace=0.08, left=0.07, right=0.98)
        ax1.tick_params(axis="both", which="major", labelsize=8)
        ax2.tick_params(axis="both", which="major", labelsize=8)

        # Plot In Situ Data (ice-lens part at the beginning) in color
        ax1.set_xlim(0,45)
        ax1.set_ylim(20,0)
        ax1.set_ylabel("Depth (m)            ")
        ax1.imshow(is_traces, cmap=is_cmap, interpolation="nearest", aspect="auto", extent=(0,45,20,0))
        ax1.text(0.975, 0.70, "A", transform=ax1.transAxes, fontsize="x-small", weight="bold", bbox=dict(facecolor="white", edgecolor="white"))

        # Subset the IceBridge data to the same edge boundaries, plot in greyscale.
        ax2.set_ylim(20,0)
        ax2.imshow(traces, cmap="gray", interpolation="nearest", aspect="auto", extent=(0,45,20,0))
        # Omit the "zero" y tick label
        ticks = ax2.get_yticks()
        ax2.set_yticks(ticks[1:])
        ax2.set_xlabel("Distance (km)")
        ax2.text(0.975, 0.70, "B", transform=ax2.transAxes, fontsize="x-small", weight="bold", bbox=dict(facecolor="white", edgecolor="white"))

        figname = r'C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\Ice Layer Cross-Validation\GPR_Icelayer_Figure.png'

        fig.savefig(figname, figsize=figsize)
        plt.close()
        print "Saved", os.path.split(figname)[1]

    def plot_validation_data_and_find_minima(self):
        '''Read the CSV produced in ::validate_reference_track_w_in_situ_data() and
        spit it out to some plots for interpolation.'''
        f = open(ICEBRIDGE_VALIDATION_OUTPUT_CSV_FILE,'r')
        lines = f.readlines()
        f.close()

        header=dict([(name.strip(),i) for i,name in enumerate(lines[0].split(","))])
        items = [line.strip().split(',') for line in lines[1:]]
        algorithm_names = numpy.array([line[header["ALGORITHM"]] for line in items])
        cutoff_values = numpy.array([float(line[header["CUTOFF_VALUE"]]) for line in items])
        continuity_ts = numpy.array([int(line[header["CONTINUITY_THRESHOLD"]]) for line in items])
        type_1_errors = numpy.array([float(line[header["TYPE_1_ERROR"]]) for line in items])
        type_2_errors = numpy.array([float(line[header["TYPE_2_ERROR"]]) for line in items])

        # Compute minimum and maximum errors here.
        max_t1 = numpy.max(type_1_errors)
        min_t1 = numpy.min(type_1_errors)
        max_t2 = numpy.max(type_2_errors)
        min_t2 = numpy.min(type_2_errors)
        print "(t1)", min_t1, max_t1,"(t2)", min_t2, max_t2

        for name in ("orig", "SG1", "SG2", "S1", "S2"):
            # Get the name for the plot.
            name_mask = (algorithm_names == name)
            if name == "orig":
                name = "Original"

            # Extract data from the CSV data.
            this_cutoff_values = cutoff_values[name_mask]
            this_continuity_ts = continuity_ts[name_mask]
            # Rescale the errors accordingly (see how this works)
            this_type_1_errors = type_1_errors[name_mask]
            this_type_2_errors = type_2_errors[name_mask]
            # Reshape
            this_cutoff_values.shape = (40,46)
            this_continuity_ts.shape = (40,46)
            this_type_1_errors.shape = (40,46)
            this_type_2_errors.shape = (40,46)

            # Scale the data for the color bars.
            ysteps = this_cutoff_values[:,0]
            xsteps = this_continuity_ts[0,:]
            minx = min(xsteps) - (xsteps[1]-xsteps[0])/2.0
            maxx = max(xsteps) + (xsteps[1]-xsteps[0])/2.0
            # We want the y-axis to go from -1.50 (top) to -0.05 (bottom), backward
            maxy = min(ysteps) - (ysteps[1]-ysteps[0])/2.0
            miny = max(ysteps) + (ysteps[1]-ysteps[0])/2.0

            fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(8,4))

            ax1.set_xlim(minx, maxx)
            ax1.set_ylim(miny, maxy)
            ax1.set_title("Type 1 Errors")
            ax1.set_ylabel("Sensitivity Cutoff (dB)")
            im1 = ax1.imshow(this_type_1_errors, cmap="inferno_r", aspect="auto", interpolation="nearest", vmin=min_t1, vmax=max_t1, extent=(minx, maxx, miny, maxy))
            plt.colorbar(im1, ax=ax1)
            ax1.text(0.065, 0.94, "A", transform=ax1.transAxes, fontsize="medium", weight="bold", bbox=dict(facecolor="white", edgecolor="white"))

            ax2.set_xlim(minx, maxx)
            ax2.set_ylim(miny, maxy)
            ax2.set_title("Type 2 Errors")
            ax2.set_xlabel("Continuity Threshold (pixels)")
            im2 = ax2.imshow(this_type_2_errors, cmap="viridis_r", aspect="auto", interpolation="nearest", vmin=min_t2, vmax=max_t2, extent=(minx, maxx, miny, maxy))
            plt.colorbar(im2, ax=ax2)
            ax2.text(0.065, 0.94, "B", transform=ax2.transAxes, fontsize="medium", weight="bold", bbox=dict(facecolor="white", edgecolor="white"))

            # Compute the overall accuracy as the combination of rescaled Type-1 and Type-2 errors, which we want ot minimize.
            normalized_error_score = this_type_1_errors + this_type_2_errors


            ax3.set_xlim(minx, maxx)
            ax3.set_ylim(miny, maxy)
            ax3.set_title("Type 1 + Type 2")
            im3 = ax3.imshow(normalized_error_score, cmap="gnuplot2_r", aspect="auto", interpolation="nearest",
                             vmin=numpy.min(normalized_error_score), vmax=numpy.max(normalized_error_score),
                             extent=(minx, maxx, miny, maxy))
            plt.colorbar(im3, ax=ax3)
            ax3.text(0.065, 0.94, "C", transform=ax3.transAxes, fontsize="medium", weight="bold", bbox=dict(facecolor="white", edgecolor="white"))

            mini = numpy.unravel_index(numpy.argmin(normalized_error_score), normalized_error_score.shape)

            ax3.scatter([this_continuity_ts[mini]], [this_cutoff_values[mini]], c="black", marker="+")

            print "Min Score: {0} (Type1: {1}, Type2: {2}), Cutoff {3}, Cont. Threshold {4}".format(normalized_error_score[mini],
                          this_type_1_errors[mini], this_type_2_errors[mini], this_cutoff_values[mini], this_continuity_ts[mini])

            plt.tight_layout()
            plt.show()

            figname = os.path.join(ICEBRIDGE_EXPORT_FOLDER, "{0}_VALIDATION_IMAGE_{1}.png".format(self.NAME, name))
            plt.savefig(figname, dpi=600)

            print "Exported", os.path.split(figname)[-1]
            plt.cla()
            plt.close()

    def validate_reference_track_w_in_situ_data(self, export_images=True):
        '''Here is where we try out different algorithsm and cutoff values, see what works best.'''
        # Check to make sure we're only doing this with the reference track.
        assert self.NAME == "20130409_01_010_012"

        # import regridded data:
        GPR_resampled = self._import_from_picklefile(ACT13_DOWNSAMPLED_PICKLEFILE)
        IceBridge_total_mask, IceBridge_hmask, IceBridge_vmask \
                      = self._import_from_picklefile(ACT13_ICEBRIDGE_SUBSET__MASK_HMASK_VMASK_PICKLEFILE)
        IceBridge_lat, IceBridge_lon, IceBridge_E, IceBridge_N, IceBridge_depths \
                      = self._import_from_picklefile(ACT13_ICEBRDIGE__LAT_LON_E_N_DEPTHS_PICKLEFILE)

        # Create a "perfect" track that identifies all pixels w/ >= 50% ice lenses
        # from GPR data as 1, all pixels <50% ice lens data as 0.  See what the "best case" score is there.
        PERFECT_TRACK = (GPR_resampled >= 0.50)
        ALL_ZEROS = numpy.zeros(GPR_resampled.shape, dtype=numpy.bool)
        ALL_ONES  = numpy.ones(GPR_resampled.shape, dtype=numpy.bool)

        print "ZEROS:  ", self._compute_score_accuracy_wrt_reference_track(ALL_ZEROS, GPR_resampled), "(all false negatives)"
        print "ONES:   ", self._compute_score_accuracy_wrt_reference_track(ALL_ONES , GPR_resampled), "(all false positives)"
        print "PERFECT:", self._compute_score_accuracy_wrt_reference_track(PERFECT_TRACK, GPR_resampled), "mean (false +, false -)"

        # Get the icebridge traces to work with:
        orig_traces = self.get_processed_traces(datatype="depth_corrected")[:len(IceBridge_depths),:]
        orig_traces = orig_traces[:,IceBridge_hmask]

        # Create ice lens array with various metrics:
        cutoff_values = numpy.arange(-1.5,0.5,0.05) # 30 increments from -1.50 to +0.45
        continuity_thresholds = numpy.arange(0,451,10) # from 0 pixels to 300, in increments of 10

        # Create empty dictionaries to store the results.
        # Each key will be an (algorithm_name, cutoff, continuity) tuple, with the (Type-1 errors, Type-2 errors) tuple as the result.
        results_dict = {}

        for cutoff_value in cutoff_values:
            print "\n~~~~~~~~~~~~~~~CUTOFF {0:0.2f}~~~~~~~~~~~~~~~~~~~~".format(cutoff_value)
            algorithm_names = ("orig", "SG1", "SG2", "S1", "S2")
            boolean_arrays = [None] * len(algorithm_names)

            # Get all pixels below that cutoff value
            boolean_arrays[0] = (orig_traces <= cutoff_value)

            # Compute the modifications here:
            # Shrink one and then grow one, then back to original size: SG1
            boolean_arrays[1] = self._boolean_shrink_and_grow(boolean_arrays[0], N=1)
            # Shrink two then grow two, then back to original size: SG2
            boolean_arrays[2] = self._boolean_shrink_and_grow(boolean_arrays[0], N=2)
            # Shrink one then back to original size: S1
            boolean_arrays[3] = self._boolean_grow_by_1(self._boolean_shrink_by_1(boolean_arrays[0], N=1), N=1)
            # Shrink two then back to original size: S2
            boolean_arrays[4] = self._boolean_grow_by_1(self._boolean_shrink_by_1(boolean_arrays[0], N=2), N=2)

            # Calculate continuity (copy/paste from v1 code) dictionary
            # Group_ID_arrays contains an MxN array of numbers for the group id's each ice lens belongs to
            group_id_arrays = [None] * len(algorithm_names)
            # Group_size_dicts contains a dictionary containing the (id:size) size of each ice lens group
            group_size_dicts = [None] * len(algorithm_names)

            for i in xrange(len(algorithm_names)):
                group_id_arrays[i], group_size_dicts[i] = self._caluculate_icelens_connectedness(boolean_arrays[i])

            boolean_image_blank = numpy.zeros(boolean_arrays[0].shape, dtype=numpy.bool)

            for continuity_t in continuity_thresholds:
                for name, boolean_array, group_id_array, group_size_dict in zip(algorithm_names,
                                                                                boolean_arrays,
                                                                                group_id_arrays,
                                                                                group_size_dicts):
                    # Use continuity dictionary to filter results and create final images.  Try this out for awhile.
                    ice_lenses_above_cutoff_size = boolean_image_blank.copy()
                    # Set each group of pixels in that category to True, only for groups larger than the cutoff
                    for group_ID in [ID for (ID,size) in group_size_dict.items() if size >= continuity_t]:
                        ice_lenses_above_cutoff_size[group_id_array == group_ID] = True

                    if export_images:
                        fname = "_{0}_CUTOFF_{1:0.2f}_THRESHOLD_{2:0>3d}.png".format(name, cutoff_value, continuity_t)
                        self.export_boolean_image(~ice_lenses_above_cutoff_size, fname)

                    results_dict[(name, cutoff_value, continuity_t)] = self._compute_score_accuracy_wrt_reference_track_BOOLEAN(ice_lenses_above_cutoff_size, PERFECT_TRACK)


        # Let's spit out a CSV of this stuff, should be pretty easy.
        print "Exporting results to", os.path.split(ICEBRIDGE_VALIDATION_OUTPUT_CSV_FILE)[-1]
        fout = open(ICEBRIDGE_VALIDATION_OUTPUT_CSV_FILE, 'w')
        fout.write("ALGORITHM,CUTOFF_VALUE,CONTINUITY_THRESHOLD,TYPE_1_ERROR,TYPE_2_ERROR\n")
        for name in algorithm_names:
            for cutoff_value in cutoff_values:
                for continuity_t in continuity_thresholds:
                    fout.write("{0},{1:f},{2:d},{3:f},{4:f}\n".format(name, cutoff_value, continuity_t,
                                                               *results_dict[(name, cutoff_value, continuity_t)]))

        return

    def _caluculate_icelens_connectedness(self, boolean_image):
        '''A parent function that iterates over an image until all pixels have found which "group" they belong to.
        Return an int array of group_ID numbers (zero are empty pixels), and a dictionary of (ID:size) pairs.'''
        group_id_array = numpy.zeros(boolean_image.shape, dtype=numpy.int)
        visited_mask_empty = numpy.zeros(boolean_image.shape, dtype=numpy.bool)
        # Visited mask cumulative -- a boolean array of all the pixels we've visited.  Starts out empty, should match boolean_image in the end
        visited_mask_cumulative = visited_mask_empty.copy()
        # Keeps track of how many pixels are in each group.
        group_size_dict = {}

        # Keep iterating while there are still pixels we haven't visited.
        still_to_visit = boolean_image
        current_group_id = 1
        while numpy.count_nonzero(still_to_visit) > 0:
            # Grab the next non-zero pixel location in still_to_visit
            next_location = numpy.unravel_index(numpy.argmax(still_to_visit), boolean_image.shape)
            # Make a copy of the empty mask for this recursive trip.
            this_visited_mask = visited_mask_empty.copy()
            # Recurse!
            self._icelens_connectedness_iterator_subfunction(boolean_image, this_visited_mask, next_location)
            # Mark all pixels in the big map with this group's id
            group_id_array[this_visited_mask] = current_group_id
            # Save the size of this group to the dictionary
            this_group_size = numpy.count_nonzero(this_visited_mask)
            group_size_dict[current_group_id] = this_group_size
            # Add these pixels to the cumulative "visited" pixels
            visited_mask_cumulative = visited_mask_cumulative | this_visited_mask
            # Subtract them from the pixels still to visit.
            still_to_visit = (boolean_image & (~visited_mask_cumulative))
            # Add one to the current running group ID
            current_group_id += 1

        return group_id_array, group_size_dict

    def _icelens_connectedness_iterator_subfunction(self, boolean_image, visited_mask, pixel_coords):
        '''An iterative function for finding all connected pixels in a region of the image.'''
        # If THIS pixel is not an ice layer OR this pixel has already been visited, return zero.
        pixel_coords_list = [pixel_coords]
        visit_directions = [0] # 0=right, 1=down, 2=left, 3=up

        fright = lambda y,x:(y,x+1)
        fdown  = lambda y,x:(y+1,x)
        fleft  = lambda y,x:(y,x-1)
        fup    = lambda y,x:(y-1,x)
        direction_dict = {0:fright, 1:fdown, 2:fleft, 3:fup}

        in_bounds = lambda y,x: (x>=0 and y>=0 and x<boolean_image.shape[1] and y<boolean_image.shape[0])

        while len(pixel_coords_list) > 0:
            direction = visit_directions[-1]
            y,x = pixel_coords_list[-1]
            visited_mask[y,x] = True

            if (0 <= direction <= 3):
                next_loc = direction_dict[direction](y,x)
                # We'll look the next direction at THIS pixel
                visit_directions[-1] += 1
                # If it's in bounds, still in the ice lens and not yet visited, go there next.
                if in_bounds(*next_loc) and boolean_image[next_loc] and (not visited_mask[next_loc]):
                    pixel_coords_list.append(next_loc)
                    visit_directions.append(0)
            elif (direction == 4):
                # Done!  Pop back to the last pixel.
                pixel_coords_list.pop()
                visit_directions.pop()
            else: # Shouldn't get here.
                assert False

        return visited_mask

    def export_downsampled_InSituGPR_track_to_picklefiles(self):
        if self.NAME != "20130409_01_010_012":
            raise ValueError("IceBridge Track must be the reference track '20130409_01_010_012', {0} doesn't fit.".format(self.NAME))

        # Collect In Situ GPR datasets
        GPR_m = InSituGPR_Manager()
        GPR = GPR_m.GPR_tracks[GPR_m.track_names.index("ACT_Transect")]
        assert isinstance(GPR, InSituGPR_Track)
        GPR_boolean = GPR.calculate_boolean_icelens_array(skip_calculations=True)
        GPR_lat, GPR_lon, GPR_elev = GPR.extract_lats_lons_elevs()
        GPR_E, GPR_N = self.return_coordiates_polarstereo(lats=GPR_lat, lons=GPR_lon)
        GPR_depths = GPR.calculate_trace_depths()
        # Collect IceBridge GPR metadata (wait on the boolean arrays)
        IceBridge_lat, IceBridge_lon = self.return_coordinates_lat_lon()
        IceBridge_E, IceBridge_N = self.return_coordiates_polarstereo()
        IceBridge_depths = self.get_sample_depths().flatten()
        IceBridge_depths = IceBridge_depths - IceBridge_depths[0]

        # The InSitu GPR extends to the right of the Icebridge AR data, so trim the GPR on the right.
        # This mask selects all the In Situ GPR points to the *left* of the right-most IceBridge point.
        GPR_hmask = (GPR_E <= numpy.max(IceBridge_E))
        # Trim the IceBridge on the left, since it extends to the left of the in situ data.
        # This mask selects (True) all the IceBridge points to the *right* of the left-most InSitu GPR point.
        IceBridge_hmask = (IceBridge_E >= numpy.min(GPR_E))
        # The IceBridge data goes deeper, create a vertical mask to just use the top ~18 m slice of it.
        #   start at 1 cm down and the IceBridge depths start at 0 m?
        # Is there an ideal "nudge" offset we need to consider for this?...
        IceBridge_vmask = (IceBridge_depths <= numpy.max(GPR_depths)).flatten()

        total_h = IceBridge_hmask.copy()
        total_h.shape = 1,total_h.shape[0]
        total_v = IceBridge_vmask.copy()
        total_v.shape = total_v.shape[0], 1
        IceBridge_total_mask = total_h & total_v
        assert IceBridge_total_mask.shape == (IceBridge_vmask.size, IceBridge_hmask.size)

        ###################################
        # Trim all the data
        ###################################
        GPR_boolean      = self._subset_array(GPR_boolean     , mask=GPR_hmask      )
        GPR_lat          = self._subset_array(GPR_lat         , mask=GPR_hmask      )
        GPR_lon          = self._subset_array(GPR_lon         , mask=GPR_hmask      )
        GPR_elev         = self._subset_array(GPR_elev        , mask=GPR_hmask      )
        GPR_E            = self._subset_array(GPR_E           , mask=GPR_hmask      )
        GPR_N            = self._subset_array(GPR_N           , mask=GPR_hmask      )
        IceBridge_E      = self._subset_array(IceBridge_E     , mask=IceBridge_hmask)
        IceBridge_N      = self._subset_array(IceBridge_N     , mask=IceBridge_hmask)
        IceBridge_lat    = self._subset_array(IceBridge_lat   , mask=IceBridge_hmask)
        IceBridge_lon    = self._subset_array(IceBridge_lon   , mask=IceBridge_hmask)
        IceBridge_depths = IceBridge_depths[IceBridge_vmask]

        '''Take the ACT13 GPR track, and resample (downsample) it into the same
        grid resolution as the adjacent IceBridge reference track. Return a boolean (?) array
        of resampled points, or perhaps a floating-point array of fraction of that grid-cell occupied by ice in the GPR track.'''
        resampled = numpy.empty((len(IceBridge_depths), len(IceBridge_E)), dtype=numpy.float)
        count = numpy.empty(resampled.shape, dtype=numpy.long)
        # Create spatial "boundaries" for each adjacent IceBridge pixel, both along-track and
        # depthwise. This helps us bin the in situ points into each closest IceBridge point.
        horizontal_boundaries = (IceBridge_E[1:] + IceBridge_E[:-1]) / 2.0
        # Add extra values on the ends to account for points on the end.
        min_leftside =  (2 * horizontal_boundaries[0]) - horizontal_boundaries[1]
        max_rightside = (2 * horizontal_boundaries[-1]) - horizontal_boundaries[-2]
        horizontal_boundaries = numpy.append(numpy.append([min_leftside], IceBridge_E),max_rightside)
        # The icebridge depths start at zero, so just go down from there and add a point to the end.
        vertical_boundaries = numpy.append(IceBridge_depths, [(2*IceBridge_depths[-1]) - IceBridge_depths[-2]])

        # Not entirely sure how to vectorize this, just iterate for this time.
        for i in xrange(resampled.shape[0]):
            print i+1, "of", resampled.shape[0]

            for j in xrange(resampled.shape[1]):
                hmask = (horizontal_boundaries[j] <= GPR_E) & (GPR_E < horizontal_boundaries[j+1])
                vmask = (vertical_boundaries[i] <= GPR_depths) & (GPR_depths < vertical_boundaries[i+1])
                subset = GPR_boolean[vmask,:][:,hmask]
                # Save the fraction of these points that are True (contain ice lenses)
                resampled[i,j] = float(numpy.count_nonzero(subset)) / float(subset.size)
                count[i,j] = subset.size

        # Export boolean images.
        self.export_image(1.0-resampled, "_ACT13_downsampled")
        self.export_boolean_image(~GPR_boolean, "_ACT13")

        # EXPORT TO PICKLEFILES
        self._export_to_picklefile(resampled, ACT13_DOWNSAMPLED_PICKLEFILE)
        self._export_to_picklefile(count, ACT13_DOWNSAMPLED_COUNT_PICKLEFILE)
        self._export_to_picklefile([IceBridge_total_mask, IceBridge_hmask, IceBridge_vmask],
                                   ACT13_ICEBRIDGE_SUBSET__MASK_HMASK_VMASK_PICKLEFILE)
        self._export_to_picklefile([IceBridge_lat, IceBridge_lon, IceBridge_E, IceBridge_N, IceBridge_depths],
                                   ACT13_ICEBRDIGE__LAT_LON_E_N_DEPTHS_PICKLEFILE)
        self._export_to_picklefile([GPR_lat, GPR_lon, GPR_elev, GPR_E, GPR_N, GPR_depths],
                                   ACT13_SUBSET__LAT_LON_ELEV_E_N_DEPTHS_PICKLEFILE)

        return resampled

    def get_trace_array(self):
        '''Return a MxN array for all the traces along that flightline, probably spanning multiple files.
        If "save_to_attrs" is set, save the traces to the Track object attributes.
        This is more efficient if using more than once, but far less space efficient if only using once.
        '''
        # If we haven't ingested the metadata yet, do so.
        if self.FILENAMES is None:
            self._read_metadata()

        if self.TRACES is not None:
            return self.TRACES

        if self.VERBOSE:
            print "Getting flightline traces", self.NAME

        table_idx = 0
        time_col = None
        if self.VERBOSE:
            print "Reading {0} files".format(len(self.FILENAMES)),

        for i,fname in enumerate(self.FILENAMES):
            if self.VERBOSE:
                print ".",

            # Most the data is in MATLAB format, however the 2014 files or not.
            # Use the h5py library to open the 2014 files and convert them into a
            # dictionary with consistent fields as the original MATLAB objects have.
            try:
                fdata = scipy.io.loadmat(fname)
            except:
                mdata = h5py.File(fname, 'r')
                fdata = self._convert_h5_dataset_to_dict(mdata)

            file_traces = fdata['Data']

            if time_col is None:
                time_col = fdata['Time']
            else:
                time_col_temp = fdata['Time']
                # Make sure all our time columns are identical.
                assert numpy.all(numpy.equal(time_col, time_col_temp))
                time_col = time_col_temp

            # Get the dimensions of the data, how many samples per trace
            if i==0:
                traces = numpy.empty((file_traces.shape[0], self.TABLE_coords_table.shape[0]), dtype=file_traces.dtype)

            # Save the file traces into the array
            traces[:, table_idx:(table_idx + file_traces.shape[1])] = file_traces
            # Increment the table index
            table_idx += file_traces.shape[1]

        if self.VERBOSE:
            print

        # Sanity check.  We should have completely filled the array here.
        assert table_idx == self.TABLE_coords_table.shape[0]

        self.TRACES = traces

        self.SAMPLE_TIMES = time_col
        self.compute_sample_depths()

        return traces

    def compute_sample_depths(self):
        '''Using the self.SAMPLE_TIMES field, compute self.sample_depths.'''
        if self.SAMPLE_TIMES is None:
            self.get_trace_array()

        self.SAMPLE_DEPTHS = self.radar_speed_m_s * self.SAMPLE_TIMES / 2.0
        return

    def return_coordinates_lat_lon(self):
        '''Return the lat/lon coordinates of every trace in the flight line. Useful for subsetting or mapping.'''
        # If we haven't ingested the metadata yet, do so.
        if self.TABLE_coords_table is None:
            self._read_metadata()

        tracenums = self.TABLE_coords_table["Flight_tracenum"]
        # Make sure all the trace numbers are in numerical order, separated by 1
        assert numpy.all((tracenums[1:] - tracenums[0:-1]) == 1)

        lats = self.TABLE_coords_table["Latitude"]
        lons = self.TABLE_coords_table["Longitude"]

        return lats, lons

    def return_coordiates_polarstereo(self, lats=None, lons=None):
        '''Return the coordinates in North Polar Stereo projection, in (eastings, northings)'''
        if (lats is None) or (lons is None):
            lats, lons = self.return_coordinates_lat_lon()
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

        return eastings, northings

    def compute_distances(self):
        '''Compute the distance (in km) of the traces in the file.'''
        eastings, northings = self.return_coordiates_polarstereo()

        # C = sqrt(A^2  + B^2)
        distances = numpy.power(numpy.power((eastings[1:] - eastings[:-1]),2) + numpy.power((northings[1:] - northings[:-1]),2), 0.5)

        if self.VERBOSE:
            print self.NAME, "Avg: {0: 5.2f}, Total: {1: 9.2f}".format(numpy.mean(distances), numpy.sum(distances))
        return distances


    def numtraces(self):
        '''Return the TOTAL number of traces in this file.  To get masked subsets,
        read the subsetted files and compute masks with the Mask_Manager object.'''
        if self.TABLE_coords_table is None:
            self._read_metadata()

        return len(self.TABLE_coords_table)

    def _read_original_surface_picks(self):
        '''Get the default surface picks from the IceBridge radar.  We will use these to improve it.'''
        if self.LIST_original_surface_picks is not None:
            return self.LIST_original_surface_picks

        if self.TABLE_coords_table is None:
            self._read_metadata()

        self.LIST_original_surface_picks = self.TABLE_coords_table['Surface']
        return self.LIST_original_surface_picks

    def _compute_original_surface_indices(self):
        surface_picks = self._read_original_surface_picks()

        # If we won't have the sample times down, get it.
        if self.SAMPLE_TIMES is None:
            self.get_trace_array()

        sample_time = self.SAMPLE_TIMES.copy()
        surface_picks = surface_picks.copy()
        # Create 2D arrays for projections
        surface_picks.shape = (surface_picks.shape[0], 1)
        sample_time.shape   = (1, sample_time.shape[0])
        match_outputs = (surface_picks == sample_time)

        # Make sure every surface pick actually had an index in the array.
        if not numpy.all(numpy.any(match_outputs, axis=1)):
            # Some of the surface picks don't seem to actually line up with time sample values.
            # We can fix this by simply getting the "closest" pick to it.
            mismatched_mask = ~numpy.any(match_outputs, axis=1)
            if self.VERBOSE:
                print "Surface mismatches.  Correcting {0} values.".format(numpy.sum(mismatched_mask))
            closest_matches_i = numpy.argmin(sample_time - surface_picks[mismatched_mask], axis=1)
            # Assign these "closest matches" as the matches themselves
            match_outputs[mismatched_mask, closest_matches_i] = True

        # Get the indices of the closest matches along that axis.
        output_array = numpy.where(match_outputs)[1]

        self.LIST_original_surface_indices = output_array
        return output_array


    def compute_surface_picks(self, export=True):
        # 1) Read original traces, subset them to get any foobars off the barfoos.
        # 2) Get masks from the SURFACE_PICK_GUIDE file

        # ONLY GO THROUGH ALL THESE COMPUTATIONS IF WE WANT TO COMPUTE AND EXPORT THIS STUFF.
        # OTHERWISE JUST READ FROM THE PICKLEFILES.
        picklefile_name = self.NAME + "_SURFACE.pickle"
        picklefile_path = os.path.join(ICEBRIDGE_SURFACE_INDICES_PICKLEFILE_FOLDER, picklefile_name)

        if not export and os.path.exists(picklefile_path):
            f = open(picklefile_path, 'r')
            improved_indices_expanded = pickle.load(f)
            f.close()
            return improved_indices_expanded

        surface_maskname = ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE
        traces = self._subset_array(self.get_trace_array(), mask=surface_maskname)
        # SHOULD I TAKE LOG HERE?  I THINK SO, LET'S TRY IT FIRST.
        traces = numpy.log10(traces)
        # Get the original indicies to use as a starter
        original_indices = self._subset_array(self._compute_original_surface_indices(), mask=surface_maskname)

        # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in original file.)
        # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
        MASK_RADIUS = 50
        vertical_span_mask = numpy.empty([MASK_RADIUS*2,], dtype=numpy.float)
        vertical_span_mask[:MASK_RADIUS] = -1.0
        vertical_span_mask[MASK_RADIUS:] = +3.0

        vertical_span_mask = vertical_span_mask * _gaussian(numpy.arange(vertical_span_mask.shape[0]),mu=(MASK_RADIUS-5),sigma=(float(MASK_RADIUS)/3.0))

        # Expand the shape to handle array broadcasting below
        vertical_span_mask.shape = vertical_span_mask.shape[0], 1

        # This is the vertical window size of the extent of the search.  Should be bigger than any jump from one surface pixel to the next.
        MASK_SEARCH_RADIUS = 150

        improved_indices = numpy.empty(original_indices.shape, dtype=original_indices.dtype)

        # Start at the left suggested vertical pixel starting point
        # IF THEY GOT IT WRONG, use the hand-picked "suggested surface pick" in the ICEBRIDGE_SURFACE_PICK_SUGGESTIONS_FILE instead.
        alternate_suggested_pick = self.mask_manager.return_suggested_starting_pixel(self.NAME)
        if alternate_suggested_pick is None:
            last_best_index = original_indices[0]
        else:
            last_best_index = alternate_suggested_pick

        # A template graph to use, just have to add in the center vertical index at each point and go from there.
        search_indices_template = numpy.sum(numpy.indices((vertical_span_mask.shape[0], 2*MASK_SEARCH_RADIUS)),axis=0) - MASK_SEARCH_RADIUS - MASK_RADIUS
        for i in xrange(traces.shape[1]):
            # Create an array of indices spanning the top-to-bottom of the MASK_SEARCH_RADIUS, and fanning out MASK_RADIUS above and below that point.
            search_indices = search_indices_template + last_best_index
            # Handle overflow indices if below zero or above max (shouldn't generally happen)... just assign to the top or bottom pixel
            search_indices[search_indices < 0] = 0
            search_indices[search_indices >= traces.shape[0]] = traces.shape[0]-1

            bestfit_sum = numpy.sum(traces[:,i][search_indices] * vertical_span_mask, axis=0)

            assert bestfit_sum.shape[0] == 2*MASK_SEARCH_RADIUS

            # Get the best fit (with the highest value from the transformation fit)
            last_best_index = search_indices[MASK_RADIUS,numpy.argmax(bestfit_sum)]
            improved_indices[i] = last_best_index

        # Erase most the little "jump" artifacts in the surface picker.
        improved_indices = self._get_rid_of_false_surface_jumps(improved_indices)
        # Must re-expand the surface indices to account for masked values (filled w/ nan)
        improved_indices_expanded = self._refill_array(improved_indices, surface_maskname)

        if export:
            radar_slice = self._return_radar_slice_given_surface(traces, improved_indices, meters_cutoff_above=5, meters_cutoff_below=10)
            idx_above, idx_below = self._radar_slice_indices_above_and_below(meters_cutoff_above=5, meters_cutoff_below=10)
            mean_returns = numpy.mean(radar_slice, axis=1)
            plt.axhline(y=0,linestyle="--",color="black")
            plt.plot(mean_returns, -numpy.arange(-idx_above,idx_below,1)*numpy.mean(self.SAMPLE_DEPTHS[1:] - self.SAMPLE_DEPTHS[:-1]))
            plt.title(self.NAME)
            plt.xlabel("Mean Strength (dB)")
            plt.ylabel("Depth (m)")
            plt.show()
            fname = self.NAME + "_MEAN_HORIZONTAL_RETURNS_5m_30m.png"
            plt.savefig(os.path.join(ICEBRIDGE_EXPORT_FOLDER,fname),dpi=600)
            print "Saved", fname
            plt.cla()
            plt.close()

            # Plot the span of the improved indices
            plt.plot(-improved_indices)
            plt.title(self.NAME)
            plt.show()
            figname = os.path.join(ICEBRIDGE_EXPORT_FOLDER, self.NAME + "__BESTFIT_V1_PLOT.png")
            plt.savefig(figname, dpi=400)
            print "Exported", os.path.split(figname)[-1]
            plt.cla()
            plt.close()

            # Export the images
            # 0-30 m depth, improves traces
            radar_slice = self._return_radar_slice_given_surface(traces, improved_indices, meters_cutoff_above=0, meters_cutoff_below=30)
            radar_slice_expanded = self._refill_array(radar_slice, surface_maskname)
            self.export_image(radar_slice_expanded, image_label="_0m_30m_BESTFIT_V1")

            # -5-30 m depth, improved traces NO DOTTED LINE
            radar_slice = self._return_radar_slice_given_surface(traces, improved_indices, meters_cutoff_above=5, meters_cutoff_below=30)
            radar_slice_expanded = self._refill_array(radar_slice, surface_maskname)
            self.export_image(radar_slice_expanded, image_label="_5m_30m_BESTFIT_V1")

            # -5-30 m depth, improved traces WITH DOTTED LINE
            idx_above, idx_below = self._radar_slice_indices_above_and_below(meters_cutoff_above=5, meters_cutoff_below=30)
            # Add a dotted line to the surface so we can see it well on the figure
            dotted_line_indices = numpy.arange(radar_slice.shape[1])
            dotted_line_mask = numpy.mod((dotted_line_indices / 3),3) == 1
            radar_slice[idx_above,dotted_line_mask] = numpy.nan
            radar_slice_expanded = self._refill_array(radar_slice, surface_maskname)
            self.export_image(radar_slice_expanded, image_label="_5m_30m_BESTFIT_DOTTED_LINE")

            # -5-30 m depth, original traces
            radar_slice = self._return_radar_slice_given_surface(traces, original_indices, meters_cutoff_above=5, meters_cutoff_below=30)
            radar_slice_expanded = self._refill_array(radar_slice, surface_maskname)
            self.export_image(radar_slice_expanded, image_label="_5m_30m_ORIGINAL")

            # Save output to a picklefile.
            f = open(picklefile_path, 'w')
            pickle.dump(improved_indices_expanded, f)
            f.close()
            print "Exported", picklefile_name

        print

        return improved_indices_expanded

    def _radar_slice_indices_above_and_below(self, meters_cutoff_above, meters_cutoff_below):
        delta_distance = numpy.mean(self.SAMPLE_DEPTHS[1:] - self.SAMPLE_DEPTHS[:-1])
        idx_above = int(numpy.round(float(meters_cutoff_above) / delta_distance))
        # Add one to the index below to include that last pixel when array-slicing
        idx_below = int(numpy.round(float(meters_cutoff_below) / delta_distance)) + 1
        return idx_above, idx_below

    def _get_rid_of_false_surface_jumps(self, surface_indices):
        '''Some of the 2011 files especially, have strong echos that are errantly being picked up as the surface.  Find these big "jumps", and get rid of them.  Use the suggested surface instead.'''
        improved_surface = surface_indices.copy()

        jumps = improved_surface[1:] - improved_surface[:-1]
        # Substitute any large jumps with brightest pixel in a window of original surface.  Do this until large jumps either go away or have all been corrected to original surface.
        for i in xrange(len(jumps)):

            # Slope windowsize = number of pixels we use to average the previous slope.
            slope_windowsize = 10
            if i < slope_windowsize:
                continue
            mean_slope = numpy.mean(numpy.array(jumps[i-slope_windowsize:i], dtype=numpy.float))

            # Find the difference of this slope from the last five stops
            difference_from_mean_slope = jumps[i] - mean_slope
            # Ignore if it's jumped less than 3 from the mean recent slope, or less than 50% greater than the mean slope at this time.
            if (difference_from_mean_slope < 5) or (difference_from_mean_slope < (1.5*mean_slope)):
                continue

            # tune settings
            jump_lookahead = 20 # Number of pixels to look ahead to see if we can find a matching down-jump
            if i+jump_lookahead > len(jumps):
                jump_lookahead = len(jumps) - i

            # This is how close the surface on the "other side" of the jump must be to the original slope to be considered for it.
            jump_magnitude_threshold = 1.10

            # See if we can find a point in the near future that would approximate the current slope.
            slopes_ahead = numpy.cumsum(jumps[i:i+jump_lookahead]) / numpy.arange(1,jump_lookahead+1)
            opposite_match = numpy.argmax(slopes_ahead <= (mean_slope * jump_magnitude_threshold))

            if opposite_match > 0:
                # We found a match, onward!
                opposite_match_index = i + opposite_match
                for j in range(i+1,opposite_match_index+1):
                    improved_surface[j] = numpy.round(improved_surface[i] + float(improved_surface[opposite_match_index+1] - improved_surface[i])*(j-i)/(opposite_match_index+1-i))

                # now recompute jumps
                jumps = improved_surface[1:] - improved_surface[:-1]
                continue

            # IF THE ABOVE DIDN'T WORK, TRY THE 'JUMP' TECHNIQUE, SEEING WHETHER AN ANOMALOUS 'JUMP' IS COUNTERBALANCED BY AN
            # OPPOSITE AND (APPROXIMATELY) EQUAL JUMP IN THE OPPOSITE DIRECTION.
            # Don't worry about any trends less than 12 pixels.  Hills do that.
            jump = jumps[i]
            if abs(jump) < 5:
                continue

            # tune settings
            jump_lookahead = 50 # Number of pixels to look ahead to see if we can find a matching down-jump
            jump_magnitude_threshold = 0.50 # What fraction of the original jump the new jump has to be (in the opposite direction) to qualify.

            # see if we can find a jump in the near-future that crosses this threshold in the other direction.  If so, we've found our counter-part
            if jump < 0:
                opposite_jump_index = numpy.argmax((jumps[i:i+jump_lookahead]) > (-jump*jump_magnitude_threshold))
            elif jump > 0:
                opposite_jump_index = numpy.argmax((jumps[i:i+jump_lookahead]) < (-jump*jump_magnitude_threshold))

            if opposite_jump_index > 0:
                opposite_jump_index += i
            else: # If we didn't find a partner opposite offset, skip and move along.
                continue

            # Linearly interpolate, get to the closest pixel
            try:
                for j in range(i+1,opposite_jump_index+1):
                    improved_surface[j] = numpy.round(improved_surface[i] + float(improved_surface[opposite_jump_index+1] - improved_surface[i])*(j-i)/(opposite_jump_index+1-i))
            except IndexError:
                print "i", i, "j", j, "opposite_jump_index", opposite_jump_index, improved_surface.shape, jumps.shape
                # Break the program here.
                100/0

            # now recompute jumps
            jumps = improved_surface[1:] - improved_surface[:-1]
            continue
        return improved_surface

    def get_radar_slice_100m(self):
        '''A quick helper function for just getting the top 100m surface slice.'''
        if self.TRACES_surface_slice_100m is not None:
            return self.TRACES_surface_slice_100m

        slice_from_file = self.read_radar_slice_from_picklefile()

        if slice_from_file is not None:
            return slice_from_file
        else:
            self.save_radar_slice_to_picklefile(meters_cutoff_above=0, meters_cutoff_below=100)
            return self.TRACES_surface_slice_100m

    def read_radar_slice_from_picklefile(self):
        '''Read the picklefile that contains the surface slice, and return it.
        If the file doesn't exist, return None.'''
        fname = self.NAME + "_SURFACE_SLICE_100M.pickle"
        pathname = os.path.join(ICEBRIDGE_SURFACE_SLICE_PICKLEFILE_FOLDER, fname)

        if os.path.exists(pathname):
            f = open(pathname, 'r')
            radar_slice = pickle.load(f)
            f.close()
            self.TRACES_surface_slice_100m = radar_slice
            return radar_slice
        else:
            return None

    def save_radar_slice_to_picklefile(self, meters_cutoff_above=0, meters_cutoff_below=100):
        '''Take the surface slice and save it to a picklefile.'''
        fname = self.NAME + "_SURFACE_SLICE_100M.pickle"
        pathname = os.path.join(ICEBRIDGE_SURFACE_SLICE_PICKLEFILE_FOLDER, fname)

        # Get the masks that tell us where we need to mask stuff out.
        maskfilenames = [ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE,ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE, ICEBRIDGE_EXCLUSIONS_LAKES_OTHER_FILE]
        # Get the surface indices and traces, mask them out according to masks above.
        surface_indices = self.compute_surface_picks(export=False)
        surface_indices = self._subset_array(surface_indices,mask=maskfilenames)
        traces = self._subset_array(self.get_trace_array(),mask=maskfilenames)

        # Get our slice
        radar_slice = self._return_radar_slice_given_surface(traces,
                                                             surface_indices,
                                                             meters_cutoff_above=meters_cutoff_above,
                                                             meters_cutoff_below=meters_cutoff_below)

        radar_slice_expanded = self._refill_array(radar_slice, maskfilenames)

        # Save it to the picklefile.
        f = open(pathname, 'w')
        pickle.dump(radar_slice_expanded,f)
        f.close()
        print "Exported", fname

        self.TRACES_surface_slice_100m = radar_slice_expanded

        return

    def _return_radar_slice_given_surface(self, traces,
                                                surface_indices,
                                                meters_cutoff_above=0,
                                                meters_cutoff_below=30):
        '''From this radar track, return a "slice" of the image above and below the surface by
        (meters_cutoff_above, meters_cutoff_below), respectively.

        Return value:
            A ((idx_below+idx_above), numtraces]-sized array of trace sample values.
        '''
        idx_above, idx_below = self._radar_slice_indices_above_and_below(meters_cutoff_above, meters_cutoff_below)

        output_traces = numpy.empty((idx_above + idx_below, traces.shape[1]), dtype=traces.dtype)

        for i,s in enumerate(surface_indices):
            try:
                output_traces[:,i] = traces[(s-idx_above):(s+idx_below), i]
            except ValueError:
                # If the surf_i is too close to one end of the array or the other, it extends beyond the edge of the array and breaks.
                if s < idx_above:
                    start, end = None, idx_above+idx_below
                elif s > (traces.shape[0] - idx_below):
                    start, end = traces.shape[0] - (idx_above + idx_below), None
                else:
                    # SHouldn't get here.
                    print i, s, traces.shape
                    assert False
                output_traces[:,i] = traces[start:end, i]

        return output_traces

    def perform_roll_correction(self, export=True, max_depth_m=100):
        '''Defines that A,C of the best-fit roll/curvature correction for the flight,
        applies that correction, exports the images, and saves the results of the slice to a picklefile.'''

        ####################################
        ## 1) Read surface slice... potentially from files created above in ::compute_surface_picks().
        trace_masks = [ICEBRIDGE_EXCLUSIONS_SURFACE_PICK_FILE,ICEBRIDGE_EXCLUSIONS_SURFACE_MISMATCH_FILE, ICEBRIDGE_EXCLUSIONS_LAKES_OTHER_FILE]
        # Read the raw traces
        traces = self.get_radar_slice_100m()
        # Mask out erroneous points, if any
        traces = self._subset_array(traces, trace_masks)
        # Log transform
        traces = numpy.log10(traces)

        #######################################
        # 2) Extract airplane roll
        roll = self._return_aircraft_roll()

        #   2a) If roll isn't available (2012), extract airplane path curvature
        if roll is None:
            # Get path curvature
            curvature = self._compute_path_curvature()
            # The curvature omits the first and last point, must fill them.
            curvature_expanded = numpy.empty([len(curvature)+2], dtype=curvature.dtype)
            curvature_expanded[1:-1] = curvature
            curvature_expanded[0] = numpy.nan
            curvature_expanded[-1] = numpy.nan

            curvature = self._subset_array(curvature_expanded, trace_masks)
        else:
            roll = numpy.abs(self._subset_array(self._return_aircraft_roll(), trace_masks))
            curvature = None

        ##############################
        # 3) Compute parabolic function to correct IceBridge data... test whether it's significant?
        if self.SAMPLE_DEPTHS is None:
            self.get_sample_depths()

        depths = (self.SAMPLE_DEPTHS - self.SAMPLE_DEPTHS[0]).flatten()
        depths = depths[0:traces.shape[0]]

        degrees_to_plot = numpy.degrees(roll if roll is not None else curvature)

        # Mask out any NaN values... these would just be in the
        data_mask = ~numpy.isnan(degrees_to_plot)
        traces_to_plot = traces[:,data_mask]
        degrees_to_plot = degrees_to_plot[data_mask]

        # Toyed with using specific weights, but found them unncessary.
        # Everything weighted equally here.
        weights = numpy.ones(degrees_to_plot.shape)

        # Our quadratic function for fitting the curve.
        def func(x,a,c):
            return a*x**2 + c

        A = numpy.empty(depths.shape, dtype=numpy.float64)
        C = numpy.empty(depths.shape, dtype=numpy.float64)

        # Compute the best A,C quadratic fit for each depth (i). We will then fit a curve to this.
        for i in range(len(depths)):
            # Perform a weighted linear regression. We really want a quadratic regression but with B==0, so just AX2 + C
            popt, pcov = scipy.optimize.curve_fit(func,
                                                  degrees_to_plot,
                                                  traces_to_plot[i,:],#mean_signal_strength,
                                                  p0 = (-0.03,numpy.mean(traces_to_plot[i,:])),
                                                  sigma = weights,
                                                  check_finite = True,
                                                  bounds = ((-numpy.inf, -numpy.inf),(0,0))
                                                  )

            A[i],C[i] = popt

        # Fit exponential curves to the A lines with depth.
        def exfuncA(x,R,S):
            return R * numpy.exp(S * x)

        try:
            Apopt, Apcov = scipy.optimize.curve_fit(exfuncA,
                                                depths,
                                                A,
                                                check_finite=True,
                                                sigma = numpy.ones(depths.shape),
                                                bounds=((-numpy.inf, -numpy.inf),(0,0)),
                                                max_nfev=1000000
                                                )
        except RuntimeError:
            Apopt, Apcov = scipy.optimize.curve_fit(exfuncA,
                                                depths,
                                                A,
                                                check_finite=True,
                                                sigma = numpy.ones(depths.shape),
                                                bounds=((-numpy.inf, -numpy.inf),(0,0)),
                                                method="dogbox",
                                                maxfev=1000000
                                                )

        A_R, A_S = Apopt
        A_computed = exfuncA(depths, A_R, A_S)

        # Calculate p-values and r-values for the function we just fit above, usling linregress
        A_slope,A_intercept,A_r_value,A_p_value,A_std_err = scipy.stats.linregress(A,A_computed)
        print "A p-value:", A_p_value
        print "A R-squared:", A_r_value**2

        def exfuncC(x,T,U,V):
            return T * numpy.exp(U * x) + V

        try:
            Cpopt, _ = scipy.optimize.curve_fit(exfuncC,
                                                depths,
                                                C,
                                                check_finite=True,
                                                sigma = numpy.ones(depths.shape),
                                                bounds=((0,-numpy.inf,-numpy.inf),(numpy.inf,0,0)),
                                                max_nfev=1000000
                                                )
        except RuntimeError:
            Cpopt, _ = scipy.optimize.curve_fit(exfuncC,
                                                depths,
                                                C,
                                                check_finite=True,
                                                sigma = numpy.ones(depths.shape),
                                                bounds=((0,-numpy.inf,-numpy.inf),(numpy.inf,0,0)),
                                                method="dogbox",
                                                maxfev=1000000
                                                )

        C_T, C_U, C_V = Cpopt
        C_computed = exfuncC(depths, C_T, C_U, C_V)


        # Calculate p-values and r-values for the function we just fit above, usling linregress
        C_slope,C_intercept,C_r_value,C_p_value,C_std_err = scipy.stats.linregress(C,C_computed)
        print "C p-value:", C_p_value
        print "C R-squared:", C_r_value**2

        ########################################
        # 5) Correct data for roll

        # Let's conflate roll/curvature and just one use array here:
        if roll is None:
            # Convert curvature to roll... this is just a scaling factor, determined by the line fit correlation btw roll & curvature derived before.
            roll = curvature
            # Fill in the first & the last values with closest matches
            roll[0] = roll[1]
            roll[-1] = roll[-2]
            # Make sure we don't have any NaN values at this point

        # There may still be an occasional point where the roll is NaN in the middle of the array.
        # For this exercise, just use the roll next to it.
        if numpy.any(numpy.isnan(roll)):
            roll_is_nan_indices = numpy.where(numpy.isnan(roll))[0]
            roll[roll_is_nan_indices] = (roll[roll_is_nan_indices-1] + roll[roll_is_nan_indices+1]) / 2
        assert not numpy.any(numpy.isnan(roll))


        # Expand a dimension for array broadcasting
        # convert to degrees
        roll = numpy.degrees(roll)
        # add dimension for array broadcasting
        roll.shape = (1,roll.shape[0])
        # add dimension for array broadcasting
        A_computed.shape = (A_computed.shape[0], 1)
        C_computed.shape = (C_computed.shape[0], 1)

        # T_side / T_nadir = (A*R^2 + C) / C
        # Therefore to convert, T_nadir = T_side * C / (A*R^2 + C)
        traces_roll_corrected = traces * (C_computed / (A_computed * numpy.power(roll, 2.0) + C_computed))

        # Get stats for the top 20 meters
        pixelcount_20m = numpy.count_nonzero((self.SAMPLE_DEPTHS - self.SAMPLE_DEPTHS[0]) < 20.0)
        mean_signal_strength_20m = numpy.mean(traces[0:pixelcount_20m,:], axis=0)
        mean_signal_strength_20m_corrected = numpy.mean(traces_roll_corrected[0:pixelcount_20m,:], axis=0)
        A_avg_20m = numpy.mean(A_computed[0:pixelcount_20m])
        C_avg_20m = numpy.mean(C_computed[0:pixelcount_20m])

        Points_below_5deg_roll = roll[roll <= 5.0]
        Points_above_5deg_roll = roll[roll > 5.0]

        mean_correction_20m_below_5deg = numpy.mean(C_avg_20m / ((A_avg_20m * Points_below_5deg_roll**2) + C_avg_20m))
        if len(Points_above_5deg_roll) == 0:
            mean_correction_20m_above_5deg = numpy.nan
        else:
            mean_correction_20m_above_5deg = numpy.mean(C_avg_20m / (A_avg_20m * Points_above_5deg_roll**2 + C_avg_20m))

        # 4) Export Images
        if export:
            fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(5,5))
            ax1.scatter(depths, A, color="tomato", marker=".")
            ax1.plot(depths, A_computed.flatten(), color="darkred", linewidth=2)
            ax1.set_ylabel("A")
            ax1.set_title(self.NAME)
            ax1.set_xlabel("Depth (m)")

            xlim = ax1.get_xlim()
            ylim = ax1.get_ylim()
            text_xpos = xlim[0] + (xlim[1] - xlim[0])*0.400
            text_ypos = ylim[0] + (ylim[1] - ylim[0])*0.070
            text_str = "$A= {0:0.4f}".format(A_R) + r"\cdot e^{" + "{0:0.4f}".format(A_S) + r"\cdot depth}$"
            ax1.text(x=text_xpos,y=text_ypos,s=text_str, bbox=dict(facecolor="white",alpha=0.80, edgecolor="white"))
            A_xpos = xlim[0] + (xlim[1] - xlim[0])*0.025
            A_ypos = ylim[0] + (ylim[1] - ylim[0])*0.870
            ax1.text(x=A_xpos,y=A_ypos,s="A", bbox=dict(facecolor="white",alpha=1.0, edgecolor="black"))


            ax2.scatter(depths, C, color="lightgreen", marker=".")
            ax2.plot(depths, C_computed.flatten(), color="darkgreen", linewidth=2)
            ax2.set_ylabel("C")
            ax2.set_xlabel("Depth (m)")

            xlim = ax2.get_xlim()
            ylim = ax2.get_ylim()
            text_xpos = xlim[0] + (xlim[1] - xlim[0])*0.400
            text_ypos = ylim[0] + (ylim[1] - ylim[0])*0.850
            text_str = "$C= {0:0.2f}".format(C_T) + r"\cdot e^{" + "{0:0.4f}".format(C_U) + r"\cdot depth}" + "{0:0.2f}$".format(C_V)
            ax2.text(x=text_xpos,y=text_ypos,s=text_str, bbox=dict(facecolor="white",alpha=0.80, edgecolor="white"))
            B_xpos = xlim[0] + (xlim[1] - xlim[0])*0.025
            B_ypos = ylim[0] + (ylim[1] - ylim[0])*0.870
            ax2.text(x=B_xpos,y=B_ypos,s="B", bbox=dict(facecolor="white",alpha=1.0, edgecolor="black"))

            plt.tight_layout()

            fname = self.NAME + "_PLOT_A_C_Curves.png"
            plt.savefig(os.path.join(ICEBRIDGE_EXPORT_FOLDER, fname), dpi=600)
            print "Exported", fname
            plt.cla()
            plt.close()

            #########################
            # Plot roll correction info
            #########################

            # Axis 1, line plot of roll/curvature.
            fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, gridspec_kw=dict(height_ratios=(1.0,1.0,2.0)), figsize=(5,5))
            ax1.plot(degrees_to_plot, color="b")
            ax1.set_ylabel(("Curvature" if curvature is not None else "Roll") + " ($^\circ$)")
            ax1.set_xlabel("Trace Number")
            ax1.set_title(self.NAME)
            xlim = ax1.get_xlim()
            ylim = ax1.get_ylim()
            A_xpos = xlim[0] + (xlim[1] - xlim[0])*0.950
            A_ypos = ylim[0] + (ylim[1] - ylim[0])*0.600
            ax1.text(x=A_xpos,y=A_ypos,s="A", bbox=dict(facecolor="white",alpha=1.0, edgecolor="black"))

            # Many y-tick labels are too big, shrink them down.
            yticks = ax1.yaxis.get_major_ticks()
            yfontsize_dict = {6:9,7:8,8:7,9:6,10:5,11:5}
            if len(yticks) > 5:
                for tick in yticks:
                    tick.label.set_fontsize(yfontsize_dict[len(yticks)])

            # Some x-tick labels are too big, shrink them down.
            xticks = ax1.xaxis.get_major_ticks()
            xfontsize_dict = {7:12,8:11,9:10,10:10,11:10,12:9}
            if len(xticks) > 6:
                for tick in xticks:
                    tick.label.set_fontsize(xfontsize_dict[len(xticks)])

            # Axis 2, mean signal strength along the
            ax2.plot(mean_signal_strength_20m, color="darkred")
            ax2.set_ylabel("GPR $\Omega$ (dB)")
            ax2.set_xlabel("Trace Number")

            xlim = ax2.get_xlim()
            ylim = ax2.get_ylim()
            B_xpos = xlim[0] + (xlim[1] - xlim[0])*0.950
            B_ypos = ylim[0] + (ylim[1] - ylim[0])*0.600
            ax2.text(x=B_xpos,y=B_ypos,s="B", bbox=dict(facecolor="white",alpha=1.0, edgecolor="black"))

            # Many y-tick labels are too big, shrink them down.
            yticks = ax2.yaxis.get_major_ticks()
            yfontsize_dict = {6:9,7:8,8:7,9:6,10:5,11:5}
            if len(yticks) > 5:
                for tick in yticks:
                    tick.label.set_fontsize(yfontsize_dict[len(yticks)])

            # Some x-tick labels are too big, shrink them down.
            xticks = ax2.xaxis.get_major_ticks()
            xfontsize_dict = {7:12,8:11,9:10,10:10,11:10,12:9}
            if len(xticks) > 6:
                for tick in xticks:
                    tick.label.set_fontsize(xfontsize_dict[len(xticks)])

            degree_range = numpy.arange(0,numpy.max(degrees_to_plot),0.01)
            # Axis 3, scatter plot with curve fit.
            ax3.scatter(degrees_to_plot, mean_signal_strength_20m_corrected[data_mask], marker=".", color="lightpink")
            ax3.scatter(degrees_to_plot, mean_signal_strength_20m[data_mask], marker=".", color="blue")
            C_20m_broadcast = numpy.ones(degree_range.shape) * C_avg_20m
            ax3.plot(degree_range, C_20m_broadcast, linewidth=2, color="red", linestyle="--")
            ax3.plot(degree_range, A_avg_20m*(degree_range**2)+C_avg_20m,linewidth=2,color="darkblue")

            xlim = ax3.get_xlim()
            ylim = ax3.get_ylim()
            text_xpos = xlim[0] + (xlim[1] - xlim[0])*0.025
            text_ypos = ylim[0] + (ylim[1] - ylim[0])*0.070
            ax3.text(x=text_xpos,y=text_ypos,s="$\Omega(\Theta)={0:0.4f}\Theta^2{1:+0.2f}$".format(A_avg_20m,C_avg_20m), bbox=dict(facecolor="white",alpha=0.80, edgecolor="white"))
            ax3.set_ylabel("GPR $\Omega$ (dB)")
            ax3.set_xlabel(("Curvature" if curvature is not None else "Roll") + " ($^\circ$)")

            xlim = ax3.get_xlim()
            ylim = ax3.get_ylim()
            C_xpos = xlim[0] + (xlim[1] - xlim[0])*0.950
            C_ypos = ylim[0] + (ylim[1] - ylim[0])*0.800
            ax3.text(x=C_xpos,y=C_ypos,s="C", bbox=dict(facecolor="white",alpha=1.0, edgecolor="black"))

            plt.tight_layout()
#            plt.show()
            fname = self.NAME + "_PLOT_ROLLCORRECT_20m.png"
            plt.savefig(os.path.join(ICEBRIDGE_EXPORT_FOLDER, fname), dpi=600)
            print "Exported", fname
            plt.cla()
            plt.close()
#            foobar

            pixelcount_30m = len(depths[depths <= 30.0])
            # 7a) Export uncorrected to image.  Get only the top 30 m
            traces_inflated = self._refill_array(traces, trace_masks)
            self.export_image(traces_inflated[:pixelcount_30m,:], image_label = "_ROLLCORRECT__BEFORE")

            # 7) Export roll-corrected to image.
            traces_roll_corrected_inflated = self._refill_array(traces_roll_corrected, trace_masks)
            self.export_image(traces_roll_corrected_inflated[:pixelcount_30m,:], image_label = "_ROLLCORRECT_AFTER")

            #######################################
            # 6) Export roll-corrected traces to picklefile.
            self.TRACES_roll_corrected = traces_roll_corrected_inflated

            picklefile_name = self.FNAME_roll_corrected_picklefile
            f = open(picklefile_name, 'w')
            pickle.dump(traces_roll_corrected_inflated, f)
            f.close()
            print "Exported", os.path.split(picklefile_name)[-1]

        print

        #################
        # 8) Return parameters
        return (A_R,
                A_S,
                A_p_value,
                A_r_value,
                C_T,
                C_U,
                C_V,
                C_p_value,
                C_r_value,
                "R" if curvature is None else "C",
                numpy.max(degrees_to_plot),
                A_avg_20m,
                C_avg_20m,
                mean_correction_20m_below_5deg,
                mean_correction_20m_above_5deg)

    def _return_aircraft_roll(self):
        '''Gets the roll of the aircraft at each trace, if available in the file.
        If unavailable, return None.'''
        if self.TABLE_coords_table is None:
            self._read_metadata()
        roll = self.TABLE_coords_table['Roll']
#        print "ROLL: min {0}, max {1}, mean {2}".format(numpy.min(roll), numpy.max(roll), numpy.mean(roll))

        if numpy.all(numpy.isnan(roll)):
            return None
        else:
            return roll

    def _compute_path_curvature(self, normalize=False):
        '''According to the X,Y locations of the plane for each trace, compute what
        the curvature of the aircraft is from the point immediately preceding each
        trace and the point immediately after.

        This helps compute a function between the roll of the airplane (computed
        and provided in all years but 2012) and the curvature of the path, which
        can be computed for all years.

        This cannot be computed for the first and last trace in the path... therefore this will return N-2
        curvatures (in radians), omitting n=0 and n=N-1 points at the start and end.'''
        ## 1. Get latitude and longitude
        lat, lon = self.return_coordinates_lat_lon()
        points = zip(lon, lat)

        ## 2. Convert all points to Polar Stereo grid projection
        gcsSR = osr.SpatialReference()
        gcsSR.SetWellKnownGeogCS("WGS84") # WGS84 geographical coordinates
        npsSR = osr.SpatialReference()
        npsSR.ImportFromEPSG(3413) # NSIDC North Pole Stereographic
        GCS_TO_NPS = osr.CoordinateTransformation(gcsSR, npsSR)

        # Transform to North Polar Stereo
        nps_points = numpy.array(GCS_TO_NPS.TransformPoints(points))
        # Cut off the zero "height" dimension
        nps_points = nps_points[:,0:2]

        ## 3. Compute rays (length N-1) from one point to the next
        rays = nps_points[1:,:] - nps_points[:-1,:]

        ## 4. Compute dot products (length N-2) between each ray.  This gives the dot product at all points but the last
        dot_products = numpy.sum(rays[1:,:] * rays[:-1,:], axis=1)

        ## 5. Compute curvatures for each point:
        ray_lengths = numpy.sqrt(numpy.sum(numpy.power(rays,2), axis=1))

        ##      Theta(n) = arccos( (ray_n .dot. ray_n+1) / (len(ray_n)*len(ray_n+1)) )
        thetas = numpy.arccos( dot_products / (ray_lengths[1:] * ray_lengths[:-1]) )
        nan_mask = numpy.isnan(thetas)

        print "CURVE:", "min", numpy.min(thetas[~nan_mask]), "max", numpy.max(thetas[~nan_mask]), "mean", numpy.mean(thetas[~nan_mask])

        ## 6. Return curvatures
        # If "normalize" is set, we normalize the curvatures for the path length at each point (using the mean of the ray before and after)
        if normalize:
            # Returns "radians per meter"
            return thetas / numpy.mean([ray_lengths[:-1],ray_lengths[1:]],axis=0)
        else:
            # Returns "radians"
            return thetas


    def perform_depth_correction(self, export=True, max_depth_m = 100):
        # Get our traces and the trace depths
        traces_all = self.get_processed_traces(datatype="roll_corrected")
        # Subset traces to mask out all NaN values (previously masked)
        mask = self._compute_boolean_mask(traces=traces_all, mask=None)
        traces = self._subset_array(traces_all, mask=None)
        depths = self.get_sample_depths(trace_array = traces)

        # Use array broadcasting here.
        depths_expanded = numpy.empty(traces.shape, dtype=depths.dtype)
        # Use array broadcasting to copy the depths into all the trace values
        depths.shape = depths.shape[0],1
        depths_expanded[:] = depths
        depths.shape = depths.shape[0]

        assert traces.shape == depths_expanded.shape

        # 1) Get the exponential curve fit
        def exfunc(y,A,B,C):
            return A * numpy.exp(B * y) + C

        popt, pcov = scipy.optimize.curve_fit(exfunc, depths_expanded.flatten(), traces.flatten(),
                                              bounds=((-numpy.inf, -numpy.inf, -numpy.inf),
                                                      ( numpy.inf,          0,          0)),
                                              max_nfev=1000000
                                             )

        A,B,C = popt
        print popt

        if export:
            # Correct the traces and normalize them.
            # Original function is Z = A * e^(By) + C
            # Inverse function to normalize AND get rid of heteroscedasticitiy is 0 = ((Z - C)/A * e^(-By) - 1.0) * e^(By)
            traces_norm = ((traces - C) / A * numpy.exp(-B * depths_expanded) - 1.0) * numpy.exp(B * depths_expanded)
            # Then divide by the standard deviation of the traces to have them normalized for variance
            # All traces  for all tracks will have a MEAN of zero and a STDDEV of 1
            traces_norm = traces_norm / (numpy.std(traces_norm))


            ###################################################
            ## Depth-correction and normalization PLOT
            ###################################################
            # We don't need to plot all the traces, just a subset (100,000 will do)
            if traces.size > 100000:
                # Subset to only plot 100000 (?) of the points
                gap = int(traces.size / 100000)
                traces_subset = traces.flatten()[::gap]
                # Contract the variability of the points to have ~ the same variability as the original points, for display only
                norm_subset = (traces_norm.flatten()/4.0)[::gap]
                depths_subset = depths_expanded.flatten()[::gap]
            else:
                traces_subset = traces.flatten()
                norm_subset = (traces_norm/4.0).flatten()
                depths_subset = depths_expanded.flatten()

            curve_fit_y = exfunc(depths, *popt)

            # 2) Make a plot, save it.
            fig = plt.figure(figsize=(5,3))
            # Plot the corrected points below, in pink/red
            plt.plot(depths_subset, norm_subset, "o", ms=1, color="salmon", fillstyle="full", mec="salmon")
            plt.axhline(y=0,color="red",ls="--",label="corrected")

            # Plot the original points atop, in blue
            plt.plot(depths_subset, traces_subset, "o", ms=1, color="lightblue", fillstyle="full", mec="lightblue")
            plt.plot(depths, curve_fit_y, color="blue",label="uncorrected")

            ax = fig.axes[0]

            equation_string = "$\Omega(y) = {0:0.3f} ".format(A) + "\cdot e^{" + "{0:0.5f}\cdot y".format(B) + "}" + "{0:0.3f}$".format(C)
            plt.text(0.04,0.10,equation_string,
                     horizontalalignment="left",
                     verticalalignment="center",
                     transform=ax.transAxes)

            # Plot legend
            handles, labels = ax.get_legend_handles_labels()
            # Even thought we plotted the corrected first, put the uncorrected first in the legend
            handles = handles[::-1]
            labels = labels[::-1]
            ax.legend(handles, labels, loc="upper right", fontsize="x-small", markerscale=0.70)

            # Title and axis labels
            plt.title(self.NAME)
            plt.xlabel("Depth $y$ (m)")
            plt.ylabel("GPR $\Omega$ (dB)")

            plt.tight_layout()
            figname = os.path.join(ICEBRIDGE_EXPORT_FOLDER, self.NAME + "_DEPTH_CURVE_PLOT.png")
            plt.savefig(figname, dpi=600)
            print "Exported", os.path.split(figname)[1]
            plt.cla()
            plt.close()

            ######################################
            ## Export picklefile
            ######################################
            traces_norm_inflated = self._refill_array(traces_norm, mask)

            f = open(self.FNAME_depth_corrected_picklefile, 'w')
            pickle.dump(traces_norm_inflated, f)
            f.close()
            print "Exported", os.path.split(self.FNAME_depth_corrected_picklefile)[-1]

            # Save to object
            self.TRACES_depth_corrected = traces_norm_inflated

            ######################################
            ## Export corrected image
            ######################################
            cutoff_30m = depths[(depths <= 30.0)].size
            traces_export = traces_norm_inflated[:cutoff_30m, :]
            self.export_image(traces_export,"_XDEPTHCORRECT_AFTER")

        # 3) Return depth-correction parameters
        print
        return popt

    def identify_ice_lenses(self, export=True, max_depth_m=20):
        '''From the ACT13 track validation performed in validate_reference_track_w_in_situ_data()
        and plot_validation_data_and_find_minima(), create ice lens images from each algorithm.'''
        # Read the traces and depths
        traces = self.get_processed_traces(datatype="depth_corrected")
        # Grab the mask and subset the traces to get rid of NaNs
        mask = self._compute_boolean_mask(traces, mask=None)
        traces = self._subset_array(traces, mask=None)
        # Retreive the depths, subset to max_depth_m (20 m)
        depths = self.get_sample_depths(trace_array = traces)
        depth_N = numpy.count_nonzero(depths <= max_depth_m)
        traces = traces[:depth_N,:]

        # We identified the minimum signal-cutoff and continuity-threshold values for each algorithm.  They
        # gave very close results.  Produce one image from each Algorithm and we will evaluate which one worked best on all the datasets.
        # These sets of cutoffs were produced from examination done in validate_reference_track_w_in_situ_data(),
        # and plot_validation_data_and_find_minima()
        ALGORITHMS = ("orig","SG1","SG1")
        CUTOFFS = (-0.45, -0.45, -0.45)
        THRESHOLDS = (0, 0, 350)

        for algorithm, cutoff, continuity_threshold in zip(ALGORITHMS, CUTOFFS, THRESHOLDS):
            # Apply the cutoff.
            boolean_traces = (traces <= cutoff)

            # Get rid of truth values in the top 2 pixels.  These are most often artifacts and should not be included
            boolean_traces[:2,:] = False

            # Apply the filter.
            if algorithm == "orig":
                pass
            elif algorithm == "SG1":
                boolean_traces = self._boolean_shrink_and_grow(boolean_traces, N=1)
            elif algorithm == "SG2":
                boolean_traces = self._boolean_shrink_and_grow(boolean_traces, N=2)
            elif algorithm == "S1":
                boolean_traces = self._boolean_grow_by_1(self._boolean_shrink_by_1(boolean_traces, N=1), N=1)
            elif algorithm == "S2":
                boolean_traces = self._boolean_grow_by_1(self._boolean_shrink_by_1(boolean_traces, N=2), N=2)

            # Perform the continuity thresholding.
            group_id_array, group_size_dict = self._caluculate_icelens_connectedness(boolean_traces)
            ice_lenses_above_cutoff_size = numpy.zeros(boolean_traces.shape, dtype=numpy.bool)

            # Set each group of pixels in that category to True, only for groups larger than the cutoff
            for group_ID in [ID for (ID,size) in group_size_dict.items() if size >= continuity_threshold]:
                ice_lenses_above_cutoff_size[group_id_array == group_ID] = True

            if export:
                traces_refilled = self._refill_array(ice_lenses_above_cutoff_size, mask)

                fname_ext = "_{0}_CUTOFF_{1:0.2f}_THRESHOLD_{2:0>3d}.png".format(algorithm, cutoff, continuity_threshold)
                self.export_boolean_image(~traces_refilled, fname_ext)

                picklefile_fname = os.path.splitext(os.path.join(ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER, self.NAME + fname_ext))[0] + ".pickle"
                self._export_to_picklefile(traces_refilled, picklefile_fname)

        if export:
            mask_picklefile_fname = os.path.join(ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER, self.NAME + "_mask.pickle")
            self._export_to_picklefile(mask, mask_picklefile_fname)

        print

        return

    def _boolean_shrink_by_1(self, orig, N=1):
        # The & operator will "shrink" the True values of the array.
        # If that pixel or any adjacent pixel (L,R,U,D) is not true, it will make that pixel not true.
        for _ in range(N):
            new = orig.copy()
            new[ :  , :-1] = new[ :  , :-1] & orig[ :  ,1:  ] # SUBSET_LEFT
            new[ :  ,1:  ] = new[ :  ,1:  ] & orig[ :  , :-1] # SUBSET_RIGHT
            new[ :-1, :  ] = new[ :-1, :  ] & orig[1:  , :  ] # SUBSET_UP
            new[1:  , :  ] = new[1:  , :  ] & orig[ :-1, :  ] # SUBSET_DOWN
            orig = new

        return new

    def _boolean_grow_by_1(self, orig, N=1):
        # The | operator will "grow" the True values of the array.
        # If that pixel or any adjacent pixel (L,R,U,D) is true, it will make that pixel True.
        for _ in range(N):
            new = orig.copy()
            new[ :  , :-1] = new[ :  , :-1] | orig[ :  ,1:  ] # SUBSET_LEFT
            new[ :  ,1:  ] = new[ :  ,1:  ] | orig[ :  , :-1] # SUBSET_RIGHT
            new[ :-1, :  ] = new[ :-1, :  ] | orig[1:  , :  ] # SUBSET_UP
            new[1:  , :  ] = new[1:  , :  ] | orig[ :-1, :  ] # SUBSET_DOWN
            orig = new

        return new

    def _boolean_shrink_and_grow(self, boolean_array, N=1):
        '''Take a boolean T/F array (assuming True == ice lenses) and use a
        "shrink and grow" method to get rid of noise.  Shrink, then grow the pixels
        by N steps to get rid of small errant strands and noise.'''
        array = boolean_array
        # SHRINK N TIMES
        for _ in range(N):
            array = self._boolean_shrink_by_1(array)
        for _ in range(N*2):
            array = self._boolean_grow_by_1(array)
        for _ in range(N):
            array = self._boolean_shrink_by_1(array)

        return array

    def _export_to_8bit_array(self, array):
        '''In order to export a function to a PNG image, use this funciton to
        export to an 8 bit unsigned integer array of scaled values.'''

        output_array = numpy.zeros(array.shape, dtype=numpy.uint8)
        excluded_mask = numpy.isnan(array)

        range_min = 0
        range_max = 2**8 - 1
        # Get the data minimum and maximum while cutting off 0.5% of outliers
        nonzero_values = array[~excluded_mask]
        data_cutoff_min = numpy.percentile(nonzero_values,  0.5)
        data_cutoff_max = numpy.percentile(nonzero_values, 99.5)

        export_array_rescaled = (array - data_cutoff_min) / (data_cutoff_max - data_cutoff_min) * range_max
        # Round to integer values
        export_array_rescaled_int = numpy.rint(export_array_rescaled)
        # Saturate at top & bottom
        export_array_rescaled_int[export_array_rescaled_int < range_min] = range_min
        export_array_rescaled_int[export_array_rescaled_int > range_max] = range_max
        # Set all numpy.nan values to zero
        export_array_rescaled_int[excluded_mask] = range_min
        # plug into the integer array (conversion from larger to smaller integers)
        output_array[:,:] = export_array_rescaled_int[:,:]

        return output_array

    def export_boolean_image(self, array, image_label=""):
        '''Create a black-and-white boolean image of this track.'''
        outfilename = self.NAME + ("_" if (len(image_label)>0 and image_label[0] != "_") else "") + image_label + (".png" if ((len(image_label) < 4) or (len(image_label) > 4 and image_label[-4:] != '.png')) else "")
        outfilepath = os.path.join(ICEBRIDGE_EXPORT_FOLDER, outfilename)
        png_file = png.from_array(array, mode="L;1")
        png_file.save(outfilepath)

        if self.VERBOSE:
            print "Exported", outfilename

        return


    def export_image(self, array, image_label=""):
        '''Create an black-and-white output image of this track.'''
        export_integers = self._export_to_8bit_array(array)

        outfilename = self.NAME + ("_" if (len(image_label)>0 and image_label[0] != "_") else "") + image_label + (".png" if ((len(image_label) < 4) or (len(image_label) > 4 and image_label[-4:] != '.png')) else "")
        outfilepath = os.path.join(ICEBRIDGE_EXPORT_FOLDER, outfilename)
        png_file = png.from_array(export_integers, mode="L")
        png_file.save(outfilepath)

        if self.VERBOSE:
            print "Exported", outfilename

        return

    def _compute_boolean_mask(self, traces=None, mask="ALL"):
        '''Reads a mask guide file, returns a boolean array of T/F values determining
        which points should be included in the dataset and which points omitted.'''
        if type(mask) in (str, list, tuple):
            return self.mask_manager.compute_mask(self.NAME, mask_filenames=mask, array_length=self.numtraces())

        if traces is not None:
            return ~numpy.any(numpy.isnan(traces), axis=0)

        raise ValueError("Must provide either traces or valid mask to _compute_boolean_mask()")

    def _subset_array(self, array, mask="ALL"):
        '''Take an inpute array and subset it by the make file name.  Return the subset array.'''
        if mask is None:
            mask = self._compute_boolean_mask(traces=array, mask=None)
        elif type( mask ) in (list,tuple,str):
            mask = self._compute_boolean_mask(traces=array, mask=mask)

        if len(array.shape) == 1:
            return array[mask]
        elif len(array.shape) == 2:
            return array[:,mask]

    def _refill_array(self, array, mask):
        '''This does the converse of _subset_array().  If an array has been subset,
        this fills back in the values that were orgininally omitted, with Nan if floating-point, else 0'''
        if type(mask) in (str,tuple,list):
            mask = self._compute_boolean_mask(mask=mask)

        # There was no subset, just return the original array.
        if numpy.all(mask):
            return array

        newshape = list(array.shape)
        newshape[-1] = len(mask)
        array_expanded = numpy.empty(newshape, dtype=array.dtype)

        if len(array.shape) == 1:
            array_expanded[:] = numpy.nan if array.dtype not in (bool,int) else 0
            array_expanded[mask] = array
        else:
            assert len(array.shape) == 2
            array_expanded[:,:] = numpy.nan if array.dtype not in (bool,int) else 0
            array_expanded[:,mask] = array

        return array_expanded

    def get_roll_corrected_traces(self):
        '''Get the roll corrected traces from the picklefile.
        If it doesn't exist, compute them to the picklefile and recursively call again.'''
        # If we've already read these traces, just return them.
        if self.TRACES_roll_corrected is not None:
            return self.TRACES_roll_corrected

        if os.path.exists(self.FNAME_roll_corrected_picklefile):
            print "Reading", os.path.split(self.FNAME_roll_corrected_picklefile)[-1]
            f = open(self.FNAME_roll_corrected_picklefile, 'r')
            self.TRACES_roll_corrected = pickle.load(f)
            f.close()
            return self.TRACES_roll_corrected

        # Create the picklefile and then return the traces
        self.perform_roll_correction(export=True)
        assert self.TRACES_roll_corrected is not None
        return self.TRACES_roll_corrected

    def get_depth_corrected_traces(self):
        '''Get the depth corrected traces from the picklefile.
        If it doesn't exist, compute them to the picklefile and recursively call again.'''
        # If we already have these traces, just return them.
        if self.TRACES_depth_corrected is not None:
            return self.TRACES_depth_corrected

        if os.path.exists(self.FNAME_depth_corrected_picklefile):
            print "Reading", os.path.split(self.FNAME_depth_corrected_picklefile)[-1]
            f = open(self.FNAME_depth_corrected_picklefile, 'r')
            self.TRACES_depth_corrected = pickle.load(f)
            f.close()
            return self.TRACES_depth_corrected

        # Create the picklefile then return the traces
        self.perform_depth_correction(export=True)
        assert self.TRACES_depth_corrected is not None
        return self.TRACES_depth_corrected

    def get_boolean_ice_traces(self):
        '''Return the boolean (T/F) traces.'''
        if self.TRACES_boolean_ice_layers is not None:
            return self.TRACES_boolean_ice_layers

        if os.path.exists(self.FNAME_ice_lenses_picklefile):
            fname = self.FNAME_ice_lenses_picklefile
            print "Reading", os.path.split(fname)[-1]
            f = open(fname, 'r')
            self.TRACES_boolean_ice_layers = pickle.load(f)
            f.close()
            return self.TRACES_boolean_ice_layers

        self.identify_ice_lenses(export=True)
        return self.TRACES_boolean_ice_layers

    def get_boolean_ice_mask(self):
        '''Return the linear trace ice mask output along with the boolean files.  Good for filtering everything.'''
        fname = os.path.join(ICEBRIDGE_BOOLEAN_RESULTS_PICKLEFILE_FOLDER, self.NAME + "_mask.pickle")
        print "Reading", os.path.split(fname)[-1]
        f = open(fname, 'r')
        mask = pickle.load(f)
        f.close()
        if self.NAME == "20120412_01_095_095":
            mask[0:9000] = False
        return mask

    def get_processed_traces(self, datatype="original"):
        assert datatype.lower() in ("original", "surface_slice", "roll_corrected", "depth_corrected", "boolean_ice_layers")

        if datatype == "original":
            return self.get_trace_array()

        elif datatype == "surface_slice":
            return self.get_radar_slice_100m()

        elif datatype == "roll_corrected":
            return self.get_roll_corrected_traces()

        elif datatype == "depth_corrected":
            return self.get_depth_corrected_traces()

        elif datatype == "boolean_ice_layers":
            return self.get_boolean_ice_traces()

    def get_sample_depths(self, trace_array = None):
        '''Either read them from the picklefile, or get the trace array and derive them.
        If "trace_array" is provided, only return the top M sample depths in that MxN array.'''
        if self.SAMPLE_DEPTHS is not None:
            if trace_array is not None and self.SAMPLE_DEPTHS.size != trace_array.shape[0]:
                return self.SAMPLE_DEPTHS[:trace_array.shape[0]]
            else:
                return self.SAMPLE_DEPTHS

        fname = self.NAME + "_SAMPLE_DEPTHS.pickle"
        pathname = os.path.join(ICEBRIDGE_SAMPLE_DEPTHS_PICKLEFILE_FOLDER, fname)

        if os.path.exists(pathname):
            f = open(pathname, 'r')
            self.SAMPLE_DEPTHS = pickle.load(f)
            f.close()
        else:
            self.get_trace_array()

        if trace_array is None:
            return self.SAMPLE_DEPTHS
        else:
            return (self.SAMPLE_DEPTHS - self.SAMPLE_DEPTHS[0]).flatten()[0:trace_array.shape[0]]


    def export_sample_depths_picklefile(self):
        '''Get the sample depths from reading the files, then export just the depths to an array.
        This saves a lot of file I/O later on.'''
        fname = self.NAME + "_SAMPLE_DEPTHS.pickle"
        pathname = os.path.join(ICEBRIDGE_SAMPLE_DEPTHS_PICKLEFILE_FOLDER, fname)

        if self.SAMPLE_DEPTHS is None:
            self.get_trace_array()

        assert self.SAMPLE_DEPTHS is not None

        f = open(pathname, 'w')
        pickle.dump(self.SAMPLE_DEPTHS, f)
        f.close()
        print "Exported", fname
        return

    def return_ice_layers_lat_lon_distance_thickness(self, masked=False):
        '''Once we have boolean ice layers calculated, return the latitude,
        longitude, elevation, and ice thickness for each trace.  If masked=True,
        return them masked out.  If masked=False, don't bother masking them.'''
        lats, lons = self.return_coordinates_lat_lon()
        boolean_traces = self.get_processed_traces(datatype="boolean_ice_layers")

        depths = self.get_sample_depths(trace_array = boolean_traces)
        depth_delta_m = numpy.mean(depths[1:] - depths[:-1])
        distances = numpy.cumsum(self.compute_distances())
        distances = numpy.append([0.0], distances)
        # Number of pixels times the thickness of each pixel
        ice_content_m = numpy.sum(boolean_traces, axis=0) * depth_delta_m

        if masked:
            mask = self.get_boolean_ice_mask()
            lats = lats[mask]
            lons = lons[mask]
            distances = distances[mask]
            ice_content_m = ice_content_m[mask]

        return lats, lons, distances, ice_content_m

def plot_surface_picking_mask_curve():
    # 3) Perform surface pick crawling threshold behavior mask (assume a step-change analysis [goes from weak->strong at surface], and continuity of surface in original file.)
    # Create a step-change mask to optimze where the returns transition from "dark" to "bright"
    MASK_RADIUS = 50
    vertical_span_mask = numpy.empty([MASK_RADIUS*2,], dtype=numpy.float)
    vertical_span_mask[:MASK_RADIUS] = -1.0
    vertical_span_mask[MASK_RADIUS:] = +3.0

    vertical_span_mask = vertical_span_mask * _gaussian(numpy.arange(vertical_span_mask.shape[0])-MASK_RADIUS,mu=0,sigma=(float(MASK_RADIUS)/3.0))

    # Expand the shape to handle array broadcasting below
    vertical_span_mask.shape = vertical_span_mask.shape[0], 1

    # A routine to plot the shape of the mask, for the supplemental section. #only need to do this once, then comment out.
    plt.axvline(x=0,linestyle="--", color="black")
    plt.plot(vertical_span_mask, numpy.arange(vertical_span_mask.shape[0])-MASK_RADIUS)
    plt.xlabel(r"Mask value $\kappa$(y)")
    plt.ylabel("Vertical pixel y")
    plt.xlim(-1.5,3.5)
    plt.ylim(50,-50)
    plt.show()
    fname = "_Surface_Identifier_Mask.png"
    plt.savefig(os.path.join(ICEBRIDGE_EXPORT_FOLDER, fname), dpi=600)
    print "Saved", fname


if __name__ == "__main__":

    ib = IceBridgeGPR_Manager_v2()
    ib.export_KML_reference_tracks()

    ib.export_ice_layer_lat_lon_distance_thicknesses()
    ib.export_smoothed_ice_layer_shapefile()

    for track in ib.tracks:
        track.DO_IT_ALL()