# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 13:37:22 2016

@author: mmacferrin

InSituGPR_Manager.py -- Handles code for processing in-situ GPR data and detecting
ice lenses.
"""
# Native Python libaries
import os
import pickle
import re

# Installed libraries
import numpy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import scipy.stats
import scipy.ndimage.filters
import osgeo.ogr as ogr
import osgeo.osr as osr

# Project files
from FirnCore_Manager import FirnCore_Manager, FirnCore_Profile
from GPR_FileData import RESAMPLED_COORDINATE_FILES,            \
                         RESAMPLED_GPR_LOGVARIANCE_PICKLEFILES, \
                         GPR_DETRENDED_PICKLEFILES,             \
                         DATA_FIGURES_OUTPUT_FOLDER,            \
                         FIGSHARE_BASE,                         \
                         EXPORT_DIR

class RadarSpeedPicker:
    '''RadarSpeedPicker: An object for determining the speed of radar signals through
    firn and ice.  Simple methods were proving a bit untenable, the ranges uncertain.
    This will (hopefully) make decisions a bit more justifiable.'''

    def __init__(self):
        '''Not much to initialize here. Just so basic variable definitions.'''
        self.ICE_DENSITY = 917.00 # kg/m3
        self.C = 299792458        # m/s, speed of light in a vaccuum.
        self.N_PURE_ICE = 1.3078  # Index of refraction of light through pure ice.
        self.N_AIR = 1.000293     # Index of refraction of light in air.


    def speed_of_light_in_ice(self, density_kg_m3, method="standard", verbose_warnings=False):
        '''Given an ice density (between 0 and 917 kg/m3), compute the approximate
        speed of light through that medium given an index of refraction of 1.3078 and
        speed in pure ice of .2294 m/ns.  (It's just a scaling between this.).
        Return value in m/s.

        'method' == 'standard' (for "pure" ice) or 'Robin' (for glacier ice)'''
        try:
            assert 0 <= density_kg_m3 <= self.ICE_DENSITY
        except AssertionError:
            if density_kg_m3 > self.ICE_DENSITY:
                if verbose_warnings:
                    print "Warning:", density_kg_m3, " kg/m3 > 917.00"
            else:
                # if we have a density less than zero, go ahead and crash.  Should never happen.
                print "DENSITIES CANNOT BE LESS THAN ZERO:", density_kg_m3
                assert False

        if method.lower() =="standard": # Assumes completely "pure" ice with zero impurities.
            # This is not a good assumption for propagation in glacial ice.
            ni = self.N_PURE_ICE   # index of refraction in pure ice, may have to increase this for impurities.
            na = self.N_AIR        # index of refraction of air (at 1 atm, 0*C)
            vi = self.C/(na+(ni-na)*(density_kg_m3/self.ICE_DENSITY))

        elif method.lower() == "robin": # Robin, 1969 and Robin, 1975, referred to in Kovacs, et al. (1995), equation 5.
            # Robin puts index of refraction of solid glacial ice at 1.78
            vi = self.C/(1+0.851*float(density_kg_m3)/1000.0)

        elif method.lower() == "robin_adjusted": # Robin seems to be giving us depths just a bit too shallow (slow) to match to the cores.  Speed up the signal just a bit here.
            vi = self.C/(1+0.700*float(density_kg_m3)/1000.0)

        else:
            raise ValueError("Unknown method for ice velocity (accept 'standard' and 'Robin'): {0}".format(method))

        return vi

    def modified_Robin_speed(self, density_kg_m3, coefficient):
        '''Using the forumula outlined in Kovacs (1995) Eq. 5, derived from Robin (1975 & 69),
        determine the speed of light using a given coefficient.  This is used by
        RadarSpeedPicker::determine_speed_by_correlation() to do a correlation study.'''
        return self.C / (1.0 + (coefficient*density_kg_m3/1000.0))


class InSituGPR_Manager():
    '''InSituGPR_Manager: Manages processes for all in situ GPR tracks and files.
    Uses: FirnCore_Manager, InSituGPR_Track
    Used by: IceBridgeGPR_Manager
    '''
    def __init__(self, GPR_FILES=RESAMPLED_GPR_LOGVARIANCE_PICKLEFILES, COORD_FILES=RESAMPLED_COORDINATE_FILES, DETRENDED_FILES=GPR_DETRENDED_PICKLEFILES, verbose=True):
        # Test whether GPR and COORDINATE files align correctly with each other
        for c,f,d in zip(COORD_FILES,GPR_FILES,DETRENDED_FILES):
            assert os.path.split(c)[-1][0:11] == os.path.split(f)[-1][0:11] == os.path.split(d)[-1][0:11]
        self.GPR_LOGVARIANCE_FILES = GPR_FILES
        self.GPR_COORDINATE_FILES = COORD_FILES
        self.GPR_DETRENDED_FILES = DETRENDED_FILES
        self.GPR_tracks = [InSituGPR_Track(f,c,d,verbose=verbose) for (f,c,d) in zip(self.GPR_LOGVARIANCE_FILES, self.GPR_COORDINATE_FILES, self.GPR_DETRENDED_FILES)]
        self.track_names = [track.name() for track in self.GPR_tracks]
        self.verbose=verbose

    def get_GPR_Track_object_from_name(self, name):
        if name in self.track_names:
            return self.GPR_tracks[self.track_names.index(name)]
        else:
            return None

    def merge_GPR_transects(self):
        '''This was originally done in the "Firn Compaction/Code/Python/GPR_Data.py" file,
        by the "GPR_Merger" object.  May leave it there, or port it over here if needed.'''

        print "GPR Transects already merged in GPR_Data.py::GPR_Merger object.  Not repeated here (yet)."
        return

    def create_GPR_variance_images(self):
        '''This was originally done in the "Firn Compaction/Code/Python/GPR_Data.py" file,
        by the "GPR_DataFile" object.  May leave it there, or port it over here if needed.
        After doing those calculations, each GPR data file is saved as a 3x13 window variance image,
        in the MALAGS/TraceImages/resampled directory, in a [NAME]_resampled_logvariance.pickle file.'''

        print "GPR variance images already cleated in GPR_Data.py::GPR_DataFile object.  Not repeated here (yet)."
        return

    def _retrieve_closest_traces(self, firn_cores):
        '''An internal utility function.  Takes the 10 closest traces from each firn core
        to each GPR track, and sorts the track distances by the closest traces.
        Closest_traces_dict is defined in InSituGPR_Manager::calculate_icelens_cutoff().
        The values contained in closest_traces_dict are contained in InSituGPR_Track::closest_traces().

        Returns a dictionary: {firncore_name, (track_names, trace_indices, trace_distances)}
        key:    firncore_name: The name of the firn core to which each group of traces is identified.
        values:   track_names:     (1) track names for each of of N closest traces.
                  track_indices:   (2) indices in each track for each of N closest traces.
                  trace_distances: (3) distance (in m) of that closest trace.
                                      Output will be sorted by this value, with
                                      lowest distance (closest trace) first.
        '''
        NUM_TRACES_PER_CORE=10

        # Get locations of each firn core
        firncore_names = [firncore.name() for firncore in firn_cores]
        lat_lon_elevs = numpy.array([firncore.lat_lon_elev() for firncore in firn_cores], dtype=numpy.float64)
        lats = lat_lon_elevs[:,0]
        lons = lat_lon_elevs[:,1]
        # Retrieve the 10 closest traces from each GPR track to each firn core.
        closest_traces = dict([(track.name(), track.closest_traces(lats, lons, topN=NUM_TRACES_PER_CORE)) for track in self.GPR_tracks])


        track_names_dict = dict([[core_name, [None]*(len(self.GPR_tracks)*NUM_TRACES_PER_CORE)] for core_name in firncore_names])
        track_indices_dict = dict([[core_name, numpy.empty([len(self.GPR_tracks)*NUM_TRACES_PER_CORE], dtype=numpy.uint)] for core_name in firncore_names])
        track_distances_dict = dict([[core_name, numpy.empty([len(self.GPR_tracks)*NUM_TRACES_PER_CORE], dtype=numpy.float)] for core_name in firncore_names])

        # Compile all the track_names, track_indices and track_distances into arrays keyed by core_names
        for i,core_name in enumerate(firncore_names):
            for j,track_name in enumerate(closest_traces.keys()):
                trace_indices, trace_distances = closest_traces[track_name]
                j_start = j*NUM_TRACES_PER_CORE
                j_end = (j+1)*NUM_TRACES_PER_CORE
                track_names_dict[core_name][j_start:j_end] = [track_name]*NUM_TRACES_PER_CORE
                track_indices_dict[core_name][j_start:j_end] = trace_indices[i,:]
                track_distances_dict[core_name][j_start:j_end] = trace_distances[i,:]

        output_dictionary = {}

        # Iterate through all the firn cores, and compile the NUM_TRACES_PER_CORE closest traces to that core
        for core_name in firncore_names:
            top_N_track_idx = numpy.argsort(track_distances_dict[core_name])[0:NUM_TRACES_PER_CORE]
            top_N_track_names = [track_names_dict[core_name][i] for i in top_N_track_idx]
            top_N_track_indices = track_indices_dict[core_name][top_N_track_idx]
            top_N_track_distances = track_distances_dict[core_name][top_N_track_idx]

            output_dictionary[core_name] = ((top_N_track_names, top_N_track_indices, top_N_track_distances))

        return output_dictionary

    def determine_best_radar_coefficient_by_correlation(self, firncore_manager, plot_figures=True, skip_calculations=False):
        '''Given a firn core and GPR signals, determine the best "speed of light" signal
        through the core that best correlates density with GPR backscatter variation.
        Since thick ice lenses are most dense, and correspond with lowest local variation,
        there should be a strong negative correlation in the data to see what the best
        speed to use is. Use only the top 10 meters

        Two tasks in this function:
        1) If "plot_figures", create plots of coefficient/correlation to see what gives the best (most negative)
        correlations.
        2) Return the coefficient that calculates the best coefficient for an
        index of refraction for the propagation of an 800 MHz radar signal through firn.
        '''
        # This simulation takes a while.  We've already done it once, if all we want
        # is the number, just return what we found through prior simulations.
        # We can run the function any time with skip_calculations=False to verify this result.
        if skip_calculations:
            return 0.734

        assert isinstance(firncore_manager, FirnCore_Manager)

        # We will ONLY be using three cores from 2013: Cores 1, 2 and 3 that
        # correspond with thick enough ice lenses to make out in the GPR signal
        firn_cores = [firncore_manager.return_core_profile(core) for core in ("core_1_2013","core_2_2013","core_3_2013")]
        # A sanity check to make sure we have the right thing back
        for firncore in firn_cores:
            assert isinstance(firncore, FirnCore_Profile)

        # Retrieve the 10 closest traces to each firn core is a dictionary.
        # Output format outlined in InSituGPR_Manager::_retrieve_closest_traces()
        firncore_traces_dict = self._retrieve_closest_traces(firn_cores)

        # INITIAL ANALYSES DETERMINE THAT THE ROBINS COEFFICIENTS PROBABLY LIE BETWEEN
        # THESE TWO VALUES.  GIVE A RANGE TO GUESS WHAT THE BEST COEFFICIENT IS TO
        # PROVIDE THE MOST "OPTIMAL" SPEED OF LIGHT FOR THIS PURPOSE.
        robins_coefficients = numpy.arange(.600,.900,0.001)

        # This is likely not the most efficient way of doing this, but for this purpose we will loop
        # through all the various possible coefficients and speeds from those. This will give us
        list_of_gpr_sample_depth_dicts = [None] * len(robins_coefficients)

        # Fetch the actual trace/sample values from the GPR files.
        gpr_traces_samples = self.extract_trace_values_from_files(firncore_traces_dict, filetype="detrended")

        # Loop through all the possible robin's coefficients, getting lists of gpr/depth profiles.
        for i,N in enumerate(robins_coefficients):

            # Get TWT along each core segment, useful for determining depth and thickness of nearby GPR layers
            gpr_core_twt_dict = dict([(core.name(), self.calculate_GPR_times_in_firn_core(core, max_depth_m=25, robins_coeff=N)) for core in firn_cores])

            # Just make sure we're not getting any misfits here.
            for key in gpr_core_twt_dict.keys():
                assert type(gpr_core_twt_dict[key]) != type(None)

            # Get the depths for each GPR sample, for traces near each core.
            # Do this ONLY for cores that have density profiles, which have been figured out in gpr_core_depths_dict
            gpr_sample_depths_dict = dict([(core.name(), self.calculate_GPR_sample_depths(core, gpr_core_twt_dict[core.name()])) \
                                            for core in firn_cores if core.name() in gpr_core_twt_dict.keys()])
            list_of_gpr_sample_depth_dicts[i] = gpr_sample_depths_dict

        mean_correlations_all = numpy.empty([len(firn_cores), robins_coefficients.shape[0]], dtype=numpy.float)

        # Then, Loop through all the firn cores and compute the correlations with firn core densities.
        for c,core in enumerate(firn_cores):
            #1. For each robin's coefficient:
            #   2. get densities of each core, and signal strength at each density (must sample gpr_depths to get the signal strength)
            #   3. create a Spearman correlation coefficient for that travel time for that core.
            #   4. Put that into an array item for that core, of length len(robins_coefficients)

            if self.verbose and plot_figures:
                print core.name(), "+++++++++++++++++++"

            track_names, track_indices, track_distances = firncore_traces_dict[core.name()]
            firncore_densities = core.densities(correct_zeros=True)
            firncore_depths = core.depths()

            spearman_correlation_rhos = numpy.empty([robins_coefficients.shape[0], len(track_names)], dtype=robins_coefficients.dtype)
            spearman_correlation_pvals = numpy.empty([robins_coefficients.shape[0], len(track_names)], dtype=robins_coefficients.dtype)

            for i in range(len(robins_coefficients)):
                gpr_sample_depths = list_of_gpr_sample_depth_dicts[i][core.name()]
                gpr_traces = gpr_traces_samples[core.name()]

                for j in range(len(track_names)):

                    trace = gpr_traces[:,j]
                    trace_values_at_densities = numpy.empty(firncore_densities.shape)
                    for k in range(len(firncore_densities)):
                        # Find the GPR sample index closest to the core depth.
                        closest_depth_idx = numpy.argmin(numpy.abs(gpr_sample_depths - (firncore_depths[k]/100.0)))
                        # Extract the GPR sample at that depth
                        trace_values_at_densities[k] = trace[closest_depth_idx]

                    # Calculate the Spearman Correlation Coefficient.  A slightly weaker
                    # measure than typical Pearson, but doesn't rely upon variables to both be
                    # normally distributed.
                    # We're only interested in looking at the TOP 10 METERS (use indices only up to 1000)
                    spearman_correlation_rhos[i,j], spearman_correlation_pvals[i,j] = \
                        scipy.stats.spearmanr(firncore_densities[:1000], trace_values_at_densities[:1000])

            if plot_figures:
                plt.figure(figsize=(5,3), dpi=150)
                for j in range(spearman_correlation_rhos.shape[1]):
                    plt.plot(robins_coefficients, spearman_correlation_rhos[:,j], color="lightgrey")
                mean_correlations = numpy.mean(spearman_correlation_rhos, axis=1)
                mean_correlations_all[c,:] = mean_correlations
                plt.plot(robins_coefficients, mean_correlations, color="darkgreen")
                plt.ylabel(r"Spearman Correlation Coefficient $\rho$", fontsize=10)
                plt.xlabel(r"GPR Velocity Coefficient $\lambda$")
                plt.xlim(0.6,0.9)
                plt.title(core.name())
                plt.tight_layout()

                if self.verbose:
                    min_correlation_idx = numpy.argmin(mean_correlations)
                    print "Strongest Correlation {0} at factor {1}".format( \
                        mean_correlations[min_correlation_idx], robins_coefficients[min_correlation_idx] )

                figure_filename = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, r"GPR_Speed_correlation_coefficients\{0}_correlations.png".format(core.name()))
                if self.verbose:
                    print "Saving", os.path.split(figure_filename)[-1]

                plt.savefig(figure_filename)

                plt.cla()
                plt.close()

            if self.verbose: print

        mean_correlations_total = numpy.mean(mean_correlations_all, axis=0)
        std_correlations_total = numpy.std(mean_correlations_all, axis=0)
        std_top = mean_correlations_total + std_correlations_total
        std_bottom = mean_correlations_total - std_correlations_total

        if plot_figures:
            plt.figure(figsize=(5,3))
            plt.plot()
            plt.fill_between(robins_coefficients, std_bottom, std_top, color="lightgrey")
            plt.plot(robins_coefficients, mean_correlations_total, color="darkgreen")
            plt.ylabel(r"Spearman Correlation Coefficient $\rho$", fontsize=10)
            plt.xlabel(r"GPR Velocity Coefficient $\lambda$")
            plt.xlim(0.6,0.9)
            plt.tight_layout()

            figure_filename = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, r"GPR_Speed_correlation_coefficients\TOTAL_correlations.png")
            if self.verbose:
                print "Saving", os.path.split(figure_filename)[-1]

            plt.savefig(figure_filename, dpi=600)

            plt.cla()
            plt.close()

        min_correlation_idx = numpy.argmin(mean_correlations_total)
        if self.verbose:
            print "Strongest Correlation {0} at factor {1}".format( \
                mean_correlations_total[min_correlation_idx], robins_coefficients[min_correlation_idx] )

        # CALCULATIONS DETERMINED THE BEST/CLOSEST CORRELATION IS USING A FACTOR OF 0.734.
        # THIS IS REGARDLESS OF WHETHER CORE 3 IS USED OR NOT. VERY STRONG RESULT.
        return robins_coefficients[min_correlation_idx]


    def calculate_GPR_times_in_firn_core(self, firn_core, max_depth_m=25, robins_coeff=0.734):
        '''Take a firn core profile and calculate the TWT totals along each cm
        of the core from density profiles.  If the max_depth is greater than the depth
        of the core, use the last meter of the core for an "average density" and
        extend it to depth from there.  Return a cm-resolution array of TWT depths in seconds.

        If the firn_core doesn't include density values, return None.'''
        # So far this seems to calculate speeds that are too high, giving us a >23-24 m depth for our GPR traces.
        # Plus, it's inconsistent with other literature.  Do we need to scale this to give the greatest correlation between
        # densities and depths?  I.e. given different scaling factors, choose one that gives the strongest correlation between
        # GPR sample logvariance values and density?
        # Let's see what the data says and go from there.
        #
        # Answer: From Kovacs (1995) we found the references to Robin (1969 & 1975) that
        # give a simple equation for light propagation in ice that are more in line with
        # our observations.  It's obvious my assumptions using "pure ice" values are not
        # valid in glacial ice with impurities in it, which slow down the propagation of light
        # significantly.  It also gives a far wider range of values for propagation in ice
        # versus light snow, which makes this type of analysis much more important.
        # Use the Robin formulas for now and see how they do.
        assert isinstance(firn_core, FirnCore_Profile)

        densities = firn_core.densities(correct_zeros=True)
        if densities is None:
            return None

        # Get an output array of the correct distance.
        if max_depth_m != None and max_depth_m > len(densities)/100.0:
            output_twt = numpy.empty((max_depth_m*100),dtype=numpy.float)
        else:
            output_twt = numpy.empty(densities.shape, dtype=numpy.float)

        # If the max_depth_m > depth of the core densities, extend the core densities.
        #  Use the mean of the last 2 m of core section to assign densities (not perfect, but good enough for this purpose)
        if len(output_twt) > len(densities):
            densities_extended = numpy.empty(output_twt.shape, dtype=densities.dtype)
            densities_extended[:len(densities)] = densities[:]
            densities_extended[len(densities):] = numpy.mean(densities[-200:])
            densities = densities_extended

        # Should be the same length now.
        assert len(densities) == len(output_twt)

        speed_object = RadarSpeedPicker()

        # Find out the speed the signal will travel through each layer.
        depth_speeds = numpy.array([speed_object.modified_Robin_speed(p, robins_coeff) for p in densities])

        # Compute the two-way travel time in each 1-cm segment (0.01 m * 2 for two-way-travel)
        # time = distance / velocity
        for i in range(len(depth_speeds)):
            output_twt[i] = (0 if i==0 else output_twt[i-1]) + \
                ((0.01*2) / depth_speeds[i])

        return output_twt

    def calculate_GPR_sample_depths(self, firn_core, core_TWTs):
        '''Given an array of TWT from a firn core--output by
        InSituGPR_Manager::calculate_GPR_distances_from_firn_core()--determine the depth
        in m of GPR traces along the vertical axis of the GPR traces near that core.  Return
        a 2024-element array of depths in meters.'''
        # The time interval per sample in the original recording was 0.1 ns, or 0.1e-9 s
        time_per_sample_s = 0.1e-9
        numsamples = 2024
        GPR_times = numpy.arange(1,numsamples+1,dtype=numpy.float64)*time_per_sample_s
        depths_interpolated = numpy.empty(GPR_times.shape, dtype=numpy.float64)

        # Each core segment is 1 cm.  Calculate the array of core depths as a step function.
        core_depths = numpy.arange(1,len(core_TWTs)+1,dtype=numpy.float64)*0.01

        assert isinstance(firn_core, FirnCore_Profile)
        # Our GPR core Two-way-travel times should be deeper than the actual GPR depths.
        # If they aren't, we need to go back and makes sure our GPR_core_TWTs are caculated to deeper ranges
        assert core_TWTs[-1] > GPR_times[-1]

        for i,sample_t in enumerate(GPR_times):
            # Find the first spot in the core that the TWT is greater than the GPR sample time
            first_idx = numpy.where(core_TWTs > sample_t)[0][0]
            # should be deeper than the first centimeter of the core
            # (we don't have cm-resolution on the GPR samples)
            assert first_idx > 0
            assert core_TWTs[first_idx] >= sample_t
            assert core_TWTs[first_idx-1] <= sample_t

            # Find the ratio of the time between this layer and the layer above it that this sample was taken.
            time_ratio = (sample_t - core_TWTs[first_idx-1]) / (core_TWTs[first_idx] - core_TWTs[first_idx-1])
            # Use this ratio to interpolate between the distances (we're talking about fractions of a centimeter here)
            depths_interpolated[i] = time_ratio * (core_depths[first_idx] - core_depths[first_idx-1]) + core_depths[first_idx-1]

            assert core_depths[first_idx-1] <= depths_interpolated[i] <= core_depths[first_idx]

        return depths_interpolated

    def extract_trace_values_from_files(self, firncore_traces_dict, filetype="detrended"):
        '''Given a dictionary of firncore_GPR_traces returned by InSituGPR_Manager::_retrieve_closest_traces(),
        open up the GPR files and pull out the actual traces.  Return a dictionary of NxS values, as such:

        input: firncore_traces_dict, a dictionary returned by InSituGPR_Manager::_retrieve_closest_traces().

        return value: dictionary of GPR trace sample values.
            keys:  Firncore names, same keys as firncore_traces_dict
            values: SxN array of GPR sample values from logvariance images.
                S: Number of sample values in the GPR, 2024 in this case.
                N: "topN" trace numbers, as defined in InSituGPR_Track::closest_traces()
            The order of samples along the N axis will be the same order as defined in the input arrays.

        If filetype=="both", make it a dictionary of two dictionaries.
            output_traces["logvariance"] = one dict with values from the original logvariance GPR files
            output_traces["detrended"] = one dict with the value from the detrended GPR files
        '''
        assert filetype in ("logvariance", "detrended", "both")

        assert type(firncore_traces_dict) == dict

        output_traces = dict()
        if filetype=="both":
            output_traces["logvariance"] = dict()
            output_traces["detrended"] = dict()

        for key in firncore_traces_dict.keys():

            for ftype in (["logvariance","detrended"] if filetype == "both" else [filetype]):
                # Extract the "topN" closest traces information,
                # from the function InSituGPR_Manager._retrieve_closest_traces()
                GPR_names, trace_indices, trace_distances = firncore_traces_dict[key]
                GPR_Tracks = [self.get_GPR_Track_object_from_name(name) for name in GPR_names]

                # Create an empty SxN array to hold trace values from all extracted traces.
                trace_array = numpy.empty((2024,len(trace_indices)), dtype=numpy.float32)

                for i in range(len(GPR_names)):
                    trace_array[:,i] = GPR_Tracks[i].get_individual_trace(trace_indices[i], filetype=ftype)

                if filetype == "both":
                    # Put the traces into the output dictionary
                    output_traces[ftype][key] = trace_array
                else:
                    output_traces[key] = trace_array

        # Return to sender.
        return output_traces


    def plot_GPR_trace_histograms_horizontal(self, firn_cores, firncore_traces_dict, gpr_sample_depths_dict, gpr_traces_samples):
        '''A diagnostic tool. For the given firn core(s), plot a histogram of the GPR sample values
        in the traces near the firn core.
        '''

        output_dir = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, "core_trace_histograms")

        for core in firn_cores:
            assert isinstance(core, FirnCore_Profile)
            core_name = core.name()

            # Extract the needed data
            gpr_depths = gpr_sample_depths_dict[core_name]
            track_names, track_indices, trace_distances = firncore_traces_dict[core_name]

            gpr_samples_log = gpr_traces_samples["logvariance"][core_name]
            gpr_samples_detrended = gpr_traces_samples["detrended"][core_name]

            # Create the figure
            outfile_tif = os.path.join(output_dir, core.name_year_first() + ".tif")
            outfile_png = os.path.join(output_dir, core.name_year_first() + ".png")
            fig, ax = plt.subplots(figsize=(6,3))

            # Set the axis labels
            ax.set_xlabel("Depth (m)")
            ax.set_ylabel("GPR Variance (dB)")

            # Plot the 10 closest traces, one at a time.
            for i in range(len(track_names)):
                ax.plot(gpr_depths, gpr_samples_log[:,i], color="grey", alpha=0.25, linewidth=0.5)
                ax.plot(gpr_depths, gpr_samples_detrended[:,i], color="red", alpha=0.25, linewidth=0.5)

            # Plot trend line through this and put up the equation.
            # A bit of cheating here, just going through the first one.
            gpr_depths_expanded = numpy.empty((gpr_depths.shape[0], gpr_samples_log.shape[1]), dtype=gpr_depths.dtype)
            # Broadcast the depths to all the GPR traces
            gpr_depths_copy = gpr_depths.copy()
            gpr_depths_copy.shape = gpr_depths.shape[0], 1
            gpr_depths_expanded[:,:] = gpr_depths_copy

            p0_log, p1_log = numpy.polyfit(gpr_depths_expanded.flatten(), gpr_samples_log.flatten(), deg=1)
            p0_det, p1_det = numpy.polyfit(gpr_depths_expanded.flatten(), gpr_samples_detrended.flatten(), deg=1)

            ax.plot(gpr_depths, gpr_depths * p0_log + p1_log, color="grey")
            ax.plot(gpr_depths, gpr_depths * p0_det + p1_det, color="red")
            ax.text(0.1, 0.15, "y = {0:0.4f}x + {1:0.2f}".format(p0_log, p1_log), ha="left", transform=ax.transAxes)

            # Set overall figure properties
            fig.suptitle("Traces  adjacent  to  " + core.name_year_first(), size=15)
            plt.tight_layout()
            fig.subplots_adjust(top=0.9)

            if self.verbose:
                print "Writing", os.path.split(outfile_tif)[-1]
            plt.savefig(outfile_tif, dpi=600)

            if self.verbose:
                print "Writing", os.path.split(outfile_png)[-1]
            plt.savefig(outfile_png, dpi=100)

            plt.cla()
            plt.close()


    def plot_GPR_trace_histograms(self, firn_cores, firncore_traces_dict, gpr_sample_depths_dict, gpr_traces_samples, filetype="both"):
        '''A diagnostic tool only.  For each firn core, plot a histogram of the GPR sample values
        in the traces near tha firn core.  This helps determine what range of cutoff
        values to use for different cores, and perhaps shed some light on what depth to use as well.'''
        output_dir = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, "core_trace_histograms")

        for core in firn_cores:
            assert isinstance(core, FirnCore_Profile)
            core_name = core.name()

            # Extract the needed data
            gpr_depths = gpr_sample_depths_dict[core_name]
            track_names, track_indices, trace_distances = firncore_traces_dict[core_name]
            if filetype=="both":
                gpr_samples_log = gpr_traces_samples["logvariance"][core_name]
                gpr_samples_detrended = gpr_traces_samples["detrended"][core_name]
            else:
                gpr_samples = gpr_traces_samples[core_name]

            core_lenses = core.create_ice_lens_array()
            core_densities = core.densities(correct_zeros=True)
            core_depths = (numpy.arange(0,len(core_densities))+1)*0.01

            # Set up the figure.
            outfile = os.path.join(output_dir, core.name_year_first() + ".png")
            fig = plt.figure(figsize=(8,6), dpi=150)
            gs = gridspec.GridSpec(1,2, width_ratios=(1,3))

            # Create the first axis, density & stratigraphy.
            ax1 = plt.subplot(gs[0])
            ax1.set_ylim(numpy.ceil(max(core_depths[-1], gpr_depths[-1])),0)
            ax1.set_xlim(0,1000)

            # Plot the ice lenses in the plot, in light blue
            if len(core_lenses.shape) == 2 and core_lenses.shape[1] == 2:
                for i in range(core_lenses.shape[0]):
                    ax1.axhspan(ymin=core_lenses[i,0]/100.0, ymax=core_lenses[i,1]/100.0, color="lightblue")

            # Plot the core density as a black-line step function
            ax1.step(core_densities, core_depths, color="black")

            # grey-out bottom
            ax1.axhspan(ymin=core_depths[-1]+0.01,ymax=ax1.get_ylim()[0],color="lightgrey")

            # Set axis and label settings.
            plt.setp(ax1, title="Density")
            ax1.set_xlabel("Density (kg m$^{-3}$)")
            ax1.set_ylabel("Depth (m)")
            plt.setp(ax1.get_xticklabels(), rotation="vertical", fontsize=10)

            # Create the second plot axis, GPR trace values.
            ax2 = plt.subplot(gs[1], sharey = ax1)
            plt.setp(ax2, title="GPR Traces")
            plt.setp(ax2.get_yticklabels(), text='')

            # Plot the 10 closest traces, one at a time.
            for i in range(len(track_names)):
                # Offset by +3 for each sample, to separate the lines.
                offset = i*4
                if filetype=="both":
                    ax2.plot(gpr_samples_log[:,i] + offset, gpr_depths, color="lightgrey")
                    ax2.plot(gpr_samples_detrended[:,i] + offset, gpr_depths)
                else:
                    ax2.plot(gpr_samples[:,i] + offset, gpr_depths)
                ax2.axvline(x=offset, linestyle=":", color="black")
            ax2.set_xlim(0,offset+8)

            # Set overall figure properties
            fig.suptitle(core.name_year_first(), size=15)
            plt.tight_layout()
            fig.subplots_adjust(top=0.9)
            if self.verbose:
                print "Writing", os.path.split(outfile)[-1]
            plt.savefig(outfile)

            plt.cla()
            plt.close()

        return

    def determine_sample_cutoff_values(self, gpr_traces, gpr_depths, max_depth_m=10.0, skip_calculations=False):
        '''Determine a reasonable range of candidate
        GPR cutoffs to use for plotting sample values and determining ice lenses.  This may
        be a hand-chosen range based upon the plots from InSituGPR_Manager::plot_GPR_trace_histograms()

        Only compute ranges based on the samples found in the top "max_depth_m" meters of firn.
        If "skip_calculations", simply return a range value rather than opening up all the radar
        files and analyzing the trace arrays.  Can only do this after we've already done
        the initial analysis.

        "gpr_traces" is assumed to be a dictionary of trace values as returned from
        InSituGPR_Manager::extract_trace_values_from_files()
        "gpr_depths" is assumed to be a dictionary of GPR depth values for each set
        of traces, as returned by InSituGPR_Manager::calculate_GPR_sample_depths() for
        each core and turned into a dictionary in InSituGPR_Manager::calculate_icelens_cutoff().

        Returns a numpy array of GPR logvariance values to use as simulations for finding
        Type-1 and Type-2 errors in the GPR analysis.
        '''
        if skip_calculations:
            return numpy.arange(3.60, 7.1001, 0.01)

        assert type(gpr_traces == dict) and type(gpr_depths) == dict and len(gpr_traces.keys()) == len(gpr_depths.keys())

        min_maxes = numpy.empty((len(gpr_traces), 2), dtype=numpy.float)
        for i,core_name in enumerate(gpr_traces.keys()):
            traces = gpr_traces[core_name]
            depths = gpr_depths[core_name]

            # Use only the traces that are shallower than the max_depth_m
            max_depth_idx = numpy.nonzero(depths > max_depth_m)[0][0]

            # Calculate the minimum and maximum value of all the traces.
            traces_subset = traces[:max_depth_idx, :]
            min_maxes[i,:] = [numpy.min(traces_subset), numpy.max(traces_subset)]
            if self.verbose:
                print core_name, min_maxes[i,:]

        # Calculate the overall min and max from all the cores.
        min_max_all = [numpy.min(min_maxes[:,0]), numpy.max(min_maxes[:,1])]
        if self.verbose:
            print "ALL", min_max_all

        # In our calculations, these values turn out to be [3.6123, 7.0807].
        # Rounded values are [3.60, 7.10]
        min_rounded = numpy.floor(min_max_all[0]*10.0)/10.0
        max_rounded = numpy.ceil(min_max_all[1]*10.0)/10.0

        step=0.01
        # Add the tiny bit to the "max_rounded" to make sure the top value is included in the range.
        return numpy.arange(min_rounded, max_rounded+(step/100.0), step)


    def plot_all_GPR_type_1_2_errors(self, firn_cores,
                                           firncore_traces_dict,
                                           gpr_sample_depths_dict,
                                           gpr_traces_samples,
                                           sample_cutoff_values_array,
                                           ice_thickness_cutoff=0.50,
                                           max_depth_m=10.0):
        '''Same functionality as ::plot_GPR_type_1_2_errors(), but just output to
        one single plot, not to 9 different plots.'''

        assert type(firncore_traces_dict) == type(gpr_sample_depths_dict) == type(gpr_traces_samples) == dict
        assert type(firn_cores) == list

        fig, axes = plt.subplots(nrows=9, ncols=1, sharex=True, sharey=False, squeeze=True, figsize=(5,9), gridspec_kw={'hspace':0})

        for i,core_profile in enumerate(firn_cores):
            assert isinstance(firn_cores[0], FirnCore_Profile)
            core_name = core_profile.name()

            # Get the core_depths, convert from cm to meters.
            core_depths = core_profile.depths() / 100.0
            # truncate core_depths to top "max_depth_m" depths
            core_depths = core_depths[core_depths <= max_depth_m]

            core_ice_lens_array = core_profile.create_ice_lens_array()
            core_ice_lens_thicknesses = core_profile.ice_lens_thicknesses() / 100.0

            gpr_depths = gpr_sample_depths_dict[core_name]
            # This isn't entirely necessary, but to save computational space, truncate
            # GPR depths to go only slightly deeper than core_depths
            gpr_depths = gpr_depths[gpr_depths <= (max_depth_m * 1.05)]

            # Get the array of boolean (T/F) values of ice content from core.
            # If it's already been calculated earlier, skip re-doing the calculations
            # and return what's already there.
            core_boolean_ice_array = core_profile.create_boolean_ice_array(skip_calculations=True)[:core_depths.shape[0]]

            # Get equivalent trace values aligned with each core depth.
            # USE numpy array projection here to speed calculations (although using more memory, which is fine here).
            # Convert core_depths to a vertical array and gpr_depths to a horizontal array
            core_depths.shape = core_depths.shape[0], 1
            gpr_depths.shape  = 1, gpr_depths.shape[0]

            # This calculates the index of the "closest trace" of the GPR sample to the
            # equivalent core depth, for each core depth.
            closest_depth_gpr_indices = numpy.argmin(numpy.abs(gpr_depths - core_depths), axis=1)
            # This array should have the same number of indices as the core depths
            assert closest_depth_gpr_indices.shape[0] == core_depths.shape[0]

            # Get the actual GPR trace values for each core_depth
            gpr_samples = gpr_traces_samples[core_name]
            gpr_samples_aligned_w_core = gpr_samples[closest_depth_gpr_indices, :]

            # Now use array projection to calculate the boolean values above & below the cutoff.
            # Create an empty 3rd dimension to the gpr_samples, this will be the sample_cutoff_values_array
            gpr_samples_aligned_w_core.shape = gpr_samples_aligned_w_core.shape + (1,)
            # Create a 1st & 2nd empty dimnesion to the sample_cutoff_values_array
            sample_cutoff_values_array.shape = (1,1) + sample_cutoff_values_array.shape

            boolean_gpr_array = (gpr_samples_aligned_w_core <= sample_cutoff_values_array)

            # Re-flatten sample_cutoff_values_array (remove extra dimensions)
            sample_cutoff_values_array = sample_cutoff_values_array.flatten()

            core_boolean_ice_copy = core_boolean_ice_array.copy()

            # Now filter out any ice lenses that are thinner than the cutoff.
            for j, thickness in enumerate(core_ice_lens_thicknesses):
                if thickness < ice_thickness_cutoff:
                    min_idx, max_idx = core_ice_lens_array[j,:]
                    if min_idx > len(core_boolean_ice_copy):
                        # We've gone off the end of the core, just break it off here.
                        break

                    assert numpy.all(core_boolean_ice_copy[min_idx:max_idx])
                    core_boolean_ice_copy[min_idx:max_idx] = False

            # Re-project core_boolean_ice_array to directly compare with this one using array projection.
            core_boolean_ice_copy.shape = core_boolean_ice_copy.shape + (1,1)

            # Calculate false positives (Type-1) and false negatives (Type-2)
            # Shape is [core_depths, "topN" GPR samples, sample_cutoff_values_array]
            # As of this running (2016.10.04), it's (1000, 10, 351), or 3.51 million elements
            # When summed over all the depths, shape goes to (10, 351)
            type_1_errors = numpy.sum((boolean_gpr_array & ~core_boolean_ice_copy), axis=0)
            type_2_errors = numpy.sum((~boolean_gpr_array & core_boolean_ice_copy), axis=0)
            # Divide by the number of layers, to give a normalized proportion.
            type_1_errors_norm = type_1_errors / float(core_boolean_ice_copy.shape[0])
            type_2_errors_norm = type_2_errors / float(core_boolean_ice_copy.shape[0])

            # Get the subplot axis to work with here.
            ax = axes[i]
            ax.axvline(x=5.0, color="black", linestyle=":")

            # Plot 10 lines, type-1 and type-2 against sample_cutoff_values
            for j in range(type_1_errors.shape[0]):
                ax.plot(sample_cutoff_values_array, type_1_errors_norm[j,:], color="lightpink")
                ax.plot(sample_cutoff_values_array, type_2_errors_norm[j,:], color="lightblue")

            # Plot the mean values
            mean_type_1_errors = numpy.mean(type_1_errors_norm, axis=0)
            mean_type_2_errors = numpy.mean(type_2_errors_norm, axis=0)
            ax.plot(sample_cutoff_values_array, mean_type_1_errors, color="magenta"   , linewidth=2, label="Type 1 Errors")
            ax.plot(sample_cutoff_values_array, mean_type_2_errors, color="mediumblue", linewidth=2, label="Type 2 Errors")
            # Find the index where the cutoff is 5.0.  It's not performing a
            # strict "==" search correctly, so we'll go with a range that should just include that single 5.0 value.
            index_50 = numpy.argmax((sample_cutoff_values_array>=4.9999)&(sample_cutoff_values_array<=5.0001))
            print core_name, "Type1 mean @ 5.0: {0}, Type2 mean @ 5.0: {1}".format(mean_type_1_errors[index_50], mean_type_2_errors[index_50])

            # include the plot name
            ax.text(5.25, 0.87, core_name, fontsize=11,
                                           verticalalignment="center",
                                           horizontalalignment="center",
                                           alpha=1.0,
                                           backgroundcolor="white")
            # include sub-figure labels, A-I, with a black box around each.
            ax.text(7.00, 0.87, ['A','B','C','D','E','F','G','H','I'][i],
                                           fontsize=11,
                                           verticalalignment="center",
                                           horizontalalignment="center",
                                           alpha=1.0,
                                           backgroundcolor="white",
                                           bbox={"pad":5,"linestyle":"solid", "edgecolor":"black","linewidth":1.25})

            if i==2:
                ax.legend(fontsize=9, loc="upper left")

            if i==4:
                ax.set_ylabel("Error Rate (fraction)")


        plt.xlim(numpy.min(sample_cutoff_values_array), numpy.max(sample_cutoff_values_array))
        # Get rid of the top "1.0" label in all but the top axis.
        for i,ax in enumerate(axes):
            ax.set_ylim(-0.05, 1.05)
            if i==0:
                continue
            locs = ax.get_yticks()
            labels = ax.get_yticklabels()
            label_texts = [label.get_text() for label in labels]
            label_texts = [('' if str(loc)=='1.0' else str(loc)) for loc in locs]
            ax.set_yticklabels(label_texts)
        plt.xlabel("GPR Logvariance Cutoff")
        plt.tight_layout()

        figure_name = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, r"Type-1-2 GPR Errors\Type12Errors_ALL_{0:0.2f}.png".format(ice_thickness_cutoff))
        plt.savefig(figure_name, dpi=600)

        if self.verbose:
            print "Plotted", os.path.split(figure_name)[-1]

        # Take-away conclusions:  In ICE cores (1-3), anything below cutoff of 5.0 minimizes
        # Type-1 Errors to near-zero.  Type-2 Errors remain but are reduced.  Virtually zero
        # Type-2 errors exist in higher cores.
        # So here's the plan.
        # 1) Pick 5.0 as a cutoff.  Apply to all transects.
        # 2) Plot ice lenses below 5.0 and take a look at what we have.
        # 3) Compare boolean arrays to cores again, look at Type-1 and Type-2 errors of THICK LENSES ONLY.

        return

    def plot_GPR_type_1_2_errors(self, firn_cores,
                                       firncore_traces_dict,
                                       gpr_sample_depths_dict,
                                       gpr_traces_samples,
                                       sample_cutoff_values_array,
                                       ice_thickness_cutoff_array=None,
                                       max_depth_m=10.0):
        '''Create paper-ready plots of Type-1 and Type-2 errors for determining the best
        GPR cutoff.  If we can simply do it here, we can return the value.  Otherwise save that
        for ::determine_GPR_cutoff_value().

        '''
        assert type(firncore_traces_dict) == type(gpr_sample_depths_dict) == type(gpr_traces_samples) == dict
        assert type(firn_cores) == list

        for core_profile in firn_cores:
            assert isinstance(firn_cores[0], FirnCore_Profile)
            core_name = core_profile.name()

            # Get the core_depths, convert from cm to meters.
            core_depths = core_profile.depths() / 100.0
            # truncate core_depths to top "max_depth_m" depths
            core_depths = core_depths[core_depths <= max_depth_m]

            core_ice_lens_array = core_profile.create_ice_lens_array()
            core_ice_lens_thicknesses = core_profile.ice_lens_thicknesses() / 100.0
            if ice_thickness_cutoff_array is None:
                ice_thickness_cutoff_array = [0.0,]

            gpr_depths = gpr_sample_depths_dict[core_name]
            # This isn't entirely necessary, but to save computational space, truncate
            # GPR depths to go only slightly deeper than core_depths
            gpr_depths = gpr_depths[gpr_depths <= (max_depth_m * 1.05)]

            # Get the array of boolean (T/F) values of ice content from core.
            # If it's already been calculated earlier, skip re-doing the calculations
            # and return what's already there.
            core_boolean_ice_array = core_profile.create_boolean_ice_array(skip_calculations=True)[:core_depths.shape[0]]

            # Get equivalent trace values aligned with each core depth.
            # USE numpy array projection here to speed calculations (although using more memory, which is fine here).
            # Convert core_depths to a vertical array and gpr_depths to a horizontal array
            core_depths.shape = core_depths.shape[0], 1
            gpr_depths.shape  = 1, gpr_depths.shape[0]

            # This calculates the index of the "closest trace" of the GPR sample to the
            # equivalent core depth, for each core depth.
            closest_depth_gpr_indices = numpy.argmin(numpy.abs(gpr_depths - core_depths), axis=1)
            # This array should have the same number of indices as the core depths
            assert closest_depth_gpr_indices.shape[0] == core_depths.shape[0]

            # Get the actual GPR trace values for each core_depth
            gpr_samples = gpr_traces_samples[core_name]
            gpr_samples_aligned_w_core = gpr_samples[closest_depth_gpr_indices, :]

            # Now use array projection to calculate the boolean values above & below the cutoff.
            # Create an empty 3rd dimension to the gpr_samples, this will be the sample_cutoff_values_array
            gpr_samples_aligned_w_core.shape = gpr_samples_aligned_w_core.shape + (1,)
            # Create a 1st & 2nd empty dimnesion to the sample_cutoff_values_array
            sample_cutoff_values_array.shape = (1,1) + sample_cutoff_values_array.shape

            boolean_gpr_array = (gpr_samples_aligned_w_core <= sample_cutoff_values_array)

            # Re-flatten sample_cutoff_values_array (remove extra dimensions)
            sample_cutoff_values_array = sample_cutoff_values_array.flatten()

            for ice_thickness_minimum in ice_thickness_cutoff_array:
                core_boolean_ice_copy = core_boolean_ice_array.copy()

                # Now filter out any ice lenses that are thinner than the cutoff.
                for i, thickness in enumerate(core_ice_lens_thicknesses):
                    if thickness < ice_thickness_minimum:
                        min_idx, max_idx = core_ice_lens_array[i,:]
                        if min_idx > len(core_boolean_ice_copy):
                            # We've gone off the end of the core, just break it off here.
                            break

                        assert numpy.all(core_boolean_ice_copy[min_idx:max_idx])
                        core_boolean_ice_copy[min_idx:max_idx] = False

                # Re-project core_boolean_ice_array to directly compare with this one using array projection.
                core_boolean_ice_copy.shape = core_boolean_ice_copy.shape + (1,1)

                # Calculate false positives (Type-1) and false negatives (Type-2)
                # Shape is [core_depths, "topN" GPR samples, sample_cutoff_values_array]
                # As of this running (2016.10.04), it's (1000, 10, 351), or 3.51 million elements
                # When summed over all the depths, shape goes to (10, 351)
                type_1_errors = numpy.sum((boolean_gpr_array & ~core_boolean_ice_copy), axis=0)
                type_2_errors = numpy.sum((~boolean_gpr_array & core_boolean_ice_copy), axis=0)
                # Divide by the number of layers, to give a normalized proportion.
                type_1_errors_norm = type_1_errors / float(core_boolean_ice_copy.shape[0])
                type_2_errors_norm = type_2_errors / float(core_boolean_ice_copy.shape[0])

                # Plot name, in our folder.
                figure_name = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, r"Type-1-2 GPR Errors\Type12Errors_{0}_{1:0.2f}.png".format(core_name, ice_thickness_minimum))
                plt.figure(figsize=(5,3), dpi=150)

                # Plot 10 lines, type-1 and type-2 against sample_cutoff_values
                for i in range(type_1_errors.shape[0]):
                    plt.plot(sample_cutoff_values_array, type_1_errors_norm[i,:], color="lightpink")
                    plt.plot(sample_cutoff_values_array, type_2_errors_norm[i,:], color="lightblue")

                plt.plot(sample_cutoff_values_array, numpy.mean(type_1_errors_norm, axis=0), color="magenta", linewidth=2, label="Type 1 Errors")
                plt.plot(sample_cutoff_values_array, numpy.mean(type_2_errors_norm, axis=0), color="mediumblue", linewidth=2, label="Type 2 Errors")

                plt.xlim(numpy.min(sample_cutoff_values_array), numpy.max(sample_cutoff_values_array))
                plt.ylim(0.0, 1.0)
                # Legend, font size 9, place in upper-left (position 2)
                plt.legend(fontsize=9, loc="upper left")
                plt.xlabel("GPR Logvariance Cutoff")
                plt.ylabel("Error Rate")
                plt.title(core_name)

                plt.tight_layout()
                plt.savefig(figure_name)
                if self.verbose:
                    print "Plotted", os.path.split(figure_name)[-1]
                plt.cla()
                plt.close()

        # Take-away conclusions:  In ICE cores (1-3), anything below cutoff of 5.0 minimizes
        # Type-1 Errors to near-zero.  Type-2 Errors remain but are reduced.  Virtually zero
        # Type-2 errors exist in higher cores.
        # So here's the plan.
        # 1) Pick 5.0 as a cutoff.  Apply to all transects.
        # 2) Plot ice lenses below 5.0 and take a look at what we have.
        # 3) If needed: "Shrink" lenses by x-pixels to recude small noise, then expand again (only to adjacent pixels) to keep layers.
        # 4) Compare boolean arrays to cores again, look at Type-1 and Type-2 errors of THICK LENSES ONLY.

        return

    def determine_GPR_cutoff_value(self):
        '''Calculate and return an appropriate cutoff value that maximizes correlation between
        ice lenses in cores and GPR values in traces.  This is key to the whole algorithm.'''
        # THIS FUNCTION IS ESSENTIALLY EMPTY.  From InSituGPR_Manager::plot_GPR_type_1_2_errors(), we
        # determined that the initial cutoff value would be 5.0.  Use this.
        return 5.0

    def create_gpr_boolean_ice_layers(self):
        '''An extension function of InSituGPR_Manager::calculate_ice_lenses().  In the function
        InSituGPR_Manager::plot_GPR_type_1_2_errors(), we realized that a cutoff of roughly 5.0
        reduced Type-2 errors in ice-laden traces will keeping Type-1 errors to a minimum.
        Use this cutoff to maximize the ice lenses we can find.  This will take a bit of
        playing but is basically outlined below in comments (copied from
        InSituGPR_Manager::plot_GPR_type_1_2_errors()).

        Goal: Compile a boolean_ice_array in each GPR track.
        '''
        cutoff = self.determine_GPR_cutoff_value()

        for track in self.GPR_tracks:
            track.calculate_boolean_icelens_array(starting_cutoff=cutoff)


    def calculate_ice_lenses(self, firncore_manager, make_all_plots=False):
        assert isinstance(firncore_manager, FirnCore_Manager)

        ## For this analysis, we only need 2013 firn cores (ignore the other years for this purpose).
        firn_cores = firncore_manager.return_core_profiles_from_year(2013)

        for firncore in firn_cores:
            # A sanity check to make sure we have the right thing back
            assert isinstance(firncore, FirnCore_Profile)

        # Retrieve the 10 closest traces to each firn core is a dictionary.
        # Output format outlined in InSituGPR_Manager::_retrieve_closest_traces()
        firncore_traces_dict = self._retrieve_closest_traces(firn_cores)

        # Get TWT along each core segment, useful for determining depth and thickness of nearby GPR layers
        gpr_core_twt_dict = dict([(core.name(), self.calculate_GPR_times_in_firn_core(core, max_depth_m=25)) for core in firn_cores])

        # Filter out cores without density measurements from our measurements (namely all core_7(b-f)_2013 cores)
        for key in gpr_core_twt_dict.keys():
            if gpr_core_twt_dict[key] is None:
                del gpr_core_twt_dict[key]
                del firncore_traces_dict[key]
                del firn_cores[[f.name() for f in firn_cores].index(key)]

        # Get the depths for each GPR sample, for traces near each core.
        # Do this ONLY for cores that have density profiles, which have been figured out in gpr_core_depths_dict
        gpr_sample_depths_dict = dict([(core.name(), self.calculate_GPR_sample_depths(core, gpr_core_twt_dict[core.name()])) \
                                        for core in firn_cores if core.name() in gpr_core_twt_dict.keys()])


        # Fetch the actual trace/sample values from the GPR files.
        gpr_traces_samples = self.extract_trace_values_from_files(firncore_traces_dict, filetype="both")

        # For testing purposes only... looking at the GPR data to determine the range of
        # cutoff values to use.  Comment this out when we're done with it.
        if make_all_plots:
            self.plot_GPR_trace_histograms_horizontal(firn_cores, firncore_traces_dict, gpr_sample_depths_dict, gpr_traces_samples)

        # Get an array of possible cutoff values to use.
        sample_cutoff_values_array = self.determine_sample_cutoff_values(gpr_traces_samples, gpr_sample_depths_dict, max_depth_m=10.0, skip_calculations=True)

        if make_all_plots:
            self.plot_all_GPR_type_1_2_errors(firn_cores,
                                              firncore_traces_dict,
                                              gpr_sample_depths_dict,
                                              gpr_traces_samples,
                                              sample_cutoff_values_array,
                                              ice_thickness_cutoff = 0.5,
                                              max_depth_m = 10.0)


        return

    def apply_vertical_variance_corrections(self):
        '''Many of the post-processed logvariance images still contain a vertical
        depth-dependent component in the signal.  Use linear-trend-line analysis to
        remove the vertical trends in the signal and give us a "detrended" signal
        instead.  This will increase noise in the vertical direction but should reduce
        bulk false positives in the signal.

        Take "_resampled_logvariance" images and convert to "resampled_logvariance_detrended" images.
        '''
        # Just call the function to process the image. After this we can always use
        # track.GPR_DETRENDED_FILENAME (which should no longer be zero) to
        for track in self.GPR_tracks:
            track.apply_vertical_variance_correction()
        return

    def export_radar_images(self, filetype="logvariance"):
        # Create PNG output images, mapped to color tables, for each of the track files.
        for track in self.GPR_tracks:
            track.export_image(filetype=filetype)

    def print_minmaxes(self):
        for t in self.GPR_tracks:
            t.print_minmax()

    def export_boolean_icelens_images(self):
        '''Create boolean ice arrays and export PNG images of those arrays.'''
        for t in self.GPR_tracks:
            t.export_boolean_ice_array()

    def create_csv_files_of_ice_content(self, max_depth_m=10.0):
        '''Create a CSV file of ice content in the top "max_depth_m" of firn
        for each radar file.'''
        for t in self.GPR_tracks:
            t.csv_of_ice_content(max_depth_m)

class InSituGPR_Track():
    '''InSituGPR_Track: Stores information and manages data for each in situ
        GPR track.
    Uses: FirnCore_Manager
    Used by: InSituGPR_Manager
    '''
    def __init__(self, GPR_FILENAME, COORD_FILENAME, GPR_DETRENDED_FILENAME, verbose=True):
        self.GPR_FILENAME = GPR_FILENAME
        self.COORD_FILENAME = COORD_FILENAME
        self.GPR_DETRENDED_FILENAME = GPR_DETRENDED_FILENAME
        self.utm_converter = None
        self.traces_logvariance = None
        self.traces_detrended = None
        self.boolean_ice_array = None
        self.verbose = verbose

    def name(self):
        match_filename = re.match("(DAT_\d{4}_A1)|(ACT_Transect)|(KANU)", os.path.split(self.GPR_FILENAME)[-1])
        match_coordfilename = re.match("(DAT_\d{4}_A1)|(ACT_Transect)|(KANU)", os.path.split(self.COORD_FILENAME)[-1])

        assert (match_filename != None) and (match_filename.group() == match_coordfilename.group())
        return match_filename.group()

    def _open(self, filetype="logvariance"):
        '''Open the GPR file, save the traces to self.traces'''
        filetype = filetype.lower()
        assert filetype in ("logvariance", "detrended")

        if (self.traces_logvariance if filetype == "logvariance" else self.traces_detrended) is None:
            fname = (self.GPR_FILENAME if filetype == "logvariance" else self.GPR_DETRENDED_FILENAME)
            if self.verbose:
                print "Opening", os.path.split(fname)[-1]
            f = open(fname, 'r')
            traces = pickle.load(f)

            if filetype=="logvariance":
                self.traces_logvariance = traces
            elif filetype=="detrended":
                self.traces_detrended = traces
            f.close()
        else:
            traces = (self.traces_logvariance if filetype == "logvariance" else self.traces_detrended)

        return traces

    def close(self):
        '''An option to save memory, close the file.  Can always be reopened later.'''
        if self.traces != None:
            del self.traces
        self.traces = None
        return

    def get_individual_trace(self, index, filetype="logvariance"):
        '''Return an individual trace from the file.
        Filetypes can be "logvariance" or "detrended", dependin on the file you want.'''
        # Get the traces.  Open the file if yet unopened.
        assert filetype.lower() in ("logvariance", "detrended")
        traces = self._open(filetype=filetype)
        return traces[:,index]

    def _set_utm_converter(self, utmZone=23, north=1):
        '''Sets and returns an osr.CoordinateTransformation object to convert WGS84 lat/lon
        coordinates into WGS84 UTM coordinates.  All coordinates must be in the
        same UTM zone.  If the object has already been created, just return it.'''
        if self.utm_converter is None:
            utmSR = osr.SpatialReference()
            utmSR.SetWellKnownGeogCS("WGS84")
            utmSR.SetUTM(utmZone,int(north))
            gcsSR = utmSR.CloneGeogCS()
            self.utm_converter = osr.CoordinateTransformation(gcsSR, utmSR)

        return self.utm_converter

    def export_image(self, filetype="logvariance"):
        '''Apply a color map and create a visual export of the image to a PNG file.
        Save in the EXPORT_DIR directory.
        filetype can be "logvariance" or "detrended".'''
        assert filetype.lower() in ("logvariance", "detrended")

        filebase = os.path.split(os.path.splitext((self.GPR_FILENAME if filetype == "logvariance" else self.GPR_DETRENDED_FILENAME))[0])[-1]
        outfilename = os.path.join(EXPORT_DIR, filebase + ".png")

        traces = self._open(filetype=filetype).copy()

#        # Unforunately PyPlot doesn't like to save any image more than 32768 pixels wide
        if traces.shape[1] > 32768:
            traces = traces[:,:32768]

        if self.verbose:
            print "Min:", numpy.min(traces), "Max:", numpy.max(traces)
        # IMG_0129: min 2.9657, max 6.7695
        # TRANSECT: min 2.8206, max 6.9996

        cmap = self.create_radar_colormap()
        traces, minval, maxval = self.adjust_minmax_traces_to_colormap(traces)

        if self.verbose:
            print "Saving", os.path.split(outfilename)[-1],

        plt.imsave(fname=outfilename, arr=traces, cmap=cmap)

        if self.verbose:
            print "Saved."

        plt.close()
        return

    def create_radar_colormap(self, old_minmax=None, new_minmax=None):
        '''Based on using the log-variance of the radar moving windows,
        this function defines and returns a colormap that will span between the
        min and max values given.  This provides a consistent color mapping for
        separate independent portions/files of the radar transect.

        NOTE: it is necessary to "adjust" the outliers of each dataset to conform
        to this new min and max (call "adjust_minmax_traces_to_colormap" below).
        Although the grand majority of each file's points lie approximately
        within this span, the actual outlier points in each log-variance array
        differ slightly, and if unadjusted would change the ranges displayed
        on each plot.  With the exception of 1-10 outlier point values (in over
        a million data points per file), the rest of the data will remain unchanged,
        and will then be plotted on a consistent color scale.

        Return value: 1) a matplotlib colormap object'''

        # I will use the custom "Rainbow" color map that has gone over well in previous iterations of this plot
        colortable = [(0,"blue"), (.42, "blue"), (.55, "cyan"), (0.62, "green"), (0.685, "yellow"), (0.75, "red"), (1,"red")]

        # If we'd like to "clip" the color table to only show values above and below thresholds (rather than the entire color table),
        # We'll clip and rescale the color map here.

        if old_minmax != None and new_minmax != None:
            oldmin, oldmax = old_minmax
            newmin, newmax = new_minmax
            newmin_fraction = (float(newmin) - oldmin) / (float(oldmax) - oldmin)
            newmax_fraction = (float(newmax) - oldmin) / (float(oldmax) - oldmin)
            expansion_ratio = 1.0/(newmax_fraction - newmin_fraction)

            for i in range(len(colortable)):
                v,c = colortable[i]
                v_new = (v-newmin_fraction)*expansion_ratio
                v_new = max(0.0, v_new)
                v_new = min(1.0, v_new)
                colortable[i] = (v_new,c)

        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("radar", colortable, N=1000)
        return cmap

    def adjust_minmax_traces_to_colormap(self, traces, defined_MIN = 2.90, defined_MAX = 7.0):
        '''Adjust the outlier points in a set of traces to have the min/max values
        provided to the function.  The justification for this (for the purpuse of display)
        is outlined in the "create_radar_colormap" object above.

        Return values:
            1) numpy array of traces with min/max points modified to fit the colormap range.
            2) the minimum of the data used (defined_MIN)
            3) the maximum of the data used (defined_MAX)'''
        # Copy the traces array
        traces_new = traces.copy()
        # Adjust minimum value in the array
        traces_new[traces_new <= max(defined_MIN, numpy.min(traces_new))] = defined_MIN
        # Adjust maximum value in the array
        traces_new[traces_new >= min(defined_MAX, numpy.max(traces_new))] = defined_MAX

        return traces_new, defined_MIN, defined_MAX

    def apply_vertical_variance_correction(self, overwrite=True, create_diagnostic_plots=False, plot_all_values=False):
        '''For a single "resampled_logvariance" image, apply a vertical detrending of the data.
        Many of the post-processed logvariance images still contain a vertical
        depth-dependent component in the signal.  Use linear-trend-line analysis to
        remove the vertical trends in the signal and give us a "detrended" signal
        instead.  This will increase noise in the vertical direction but should reduce
        bulk false positives in the signal.

        Take "_resampled_logvariance" images and convert to "resampled_logvariance_detrended" images.
        '''
        if self.verbose:
            print "Detrending", self.name()

        # If we've already created this file and don't wish to overwrite,
        #  just return the file name and fuggedaboudit.
        if os.path.exists(self.GPR_DETRENDED_FILENAME) and overwrite==False:
            return self.GPR_DETRENDED_FILENAME

        # The original "logvariance" traces.
        log_traces = self._open(filetype="logvariance")
        if numpy.any(numpy.isnan(log_traces)):
            raise ValueError(os.path.split(self.GPR_FILENAME)[-1] + " contains NaN values.")

        # Vertically detrended traces.  Make an empty copy of the log_traces array (no need to copy data)
        det_traces = numpy.empty(log_traces.shape, dtype=log_traces.dtype)

        # Make sure the axes are in the order I think they are.
        assert log_traces.shape[0] == 2024 and len(log_traces.shape) == 2

        # Create vectors for linear regression outputs
        v_slopes = numpy.empty([log_traces.shape[1]], dtype=numpy.float64)
        v_intercepts = numpy.empty([log_traces.shape[1]], dtype=numpy.float64)
        v_r_values = numpy.empty([log_traces.shape[1]], dtype=numpy.float64)
        v_p_values = numpy.empty([log_traces.shape[1]], dtype=numpy.float64)
        v_stderrs = numpy.empty([log_traces.shape[1]], dtype=numpy.float64)

        #x-axis is simply a step array of sample numbers
        sample_numbers = numpy.arange(0, log_traces.shape[0], step=1, dtype=numpy.float)
        trace_numbers  = numpy.arange(0, log_traces.shape[1], step=1, dtype=numpy.float)

        # could *probably* vectorize this, but a for() loop is working for now.
        for i in range(log_traces.shape[1]):
            v_slopes[i], v_intercepts[i], v_r_values[i], v_p_values[i], v_stderrs[i] = \
            scipy.stats.linregress(sample_numbers, log_traces[:,i])

        # Create gaussian filter of slopes
        # Sigma == about every 1 kilometer
        N = 1000/1.5
        v_slopes_gaussian = scipy.ndimage.gaussian_filter(v_slopes, sigma=N, order=0, mode="mirror")
        # Create a version where all positive values are equal to zero.
        # We don't want to weaken deep signals just because there's big ice near the surface.
        v_slopes_gaussian_neg_only = v_slopes_gaussian.copy()
        slopes_mask = v_slopes_gaussian > 0
        v_slopes_gaussian_neg_only[slopes_mask] = 0.0

        # Create the same smoothed signal for intercept, to normalize all the data at zero.
        v_intercepts_gaussian = scipy.ndimage.gaussian_filter(v_intercepts, sigma=N, order=0, mode="mirror")
        # In traces where we're omitting the slope correction, just use the mean of the vertical traces as a bulk offset.
        v_intercepts_neg_only = v_intercepts.copy()
        v_intercepts_neg_only[slopes_mask] = numpy.mean(log_traces, axis=0)[slopes_mask]
        v_intercepts_neg_only_gaussian = scipy.ndimage.gaussian_filter(v_intercepts_neg_only, sigma=N, order=0, mode="mirror")

        # Create some diagnostic output plots here, only if asked for.
        if create_diagnostic_plots: # and self.name() != "DAT_0154_A1":
            outfolder = os.path.join(DATA_FIGURES_OUTPUT_FOLDER, "GPR_detrending_figures")
            plotfilename = os.path.join(outfolder, self.name() + "_radar_trendlines_{0}.tif".format("all" if plot_all_values else "SlopeInt"))
            if plot_all_values == True:
                fig = plt.figure(figsize=(5.,12.),dpi=300,facecolor='w')
                gs = gridspec.GridSpec(5,1, height_ratios=(1,1,1,1,1))
            else:
                fig = plt.figure(figsize=(5.,5.),dpi=300,facecolor='w')
                gs = gridspec.GridSpec(2,1, height_ratios=(1,1))

            # Create Axis-1, slopes.
            ax1 = plt.subplot(gs[0])
            slope_factor = 10e-5
            ax1.plot(trace_numbers, v_slopes / slope_factor, color="lightgrey")
            ax1.plot(trace_numbers, v_slopes_gaussian / slope_factor, color="lightblue")
            ax1.plot(trace_numbers, v_slopes_gaussian_neg_only / slope_factor, color="darkblue")
            # Set Axis-1 and label settings.
            ax1.set_ylabel(r"Slope ($\times 10^{-5}$)")
            ax1.xaxis.set_ticklabels([])
            plt.tick_params(axis="y", labelsize=9)

            # Create Axis-2, intercepts
            if plot_all_values:
                ax2 = plt.subplot(gs[1], sharex=ax1)
            else:
                ax2 = plt.subplot(gs[1])
            ax2.plot(trace_numbers, v_intercepts, color="lightgrey")
            ax2.plot(trace_numbers, v_intercepts_gaussian, color="pink")
            ax2.plot(trace_numbers, v_intercepts_neg_only_gaussian, color="red")
            # Set Axis-1 and label settings.
            ax2.set_ylabel("Intercept")
            plt.tick_params(axis="y", labelsize=9)

            if plot_all_values:
                # Create Axis-3, r values
                ax3 = plt.subplot(gs[2], sharex=ax1)
                ax3.plot(trace_numbers, v_r_values, color="green")
                # Set Axis-1 and label settings.
                ax3.set_ylabel("R value")
                plt.tick_params(axis="y", labelsize=9)

                # Not using P-values, omit here.
                # Create Axis-4, p values
                # Extremely small p-values seem to be breaking the plot.
                # All of the values are extremely tiny, so just set to zero for this particular plot.
                v_p_values_to_plot = v_p_values.copy()
                tiny_mask = (v_p_values_to_plot < 1e-20)
                v_p_values_to_plot[tiny_mask] = 0.0

                ax4 = plt.subplot(gs[3], sharex=ax1)
                ax4.plot(trace_numbers, v_p_values_to_plot, color="purple")
                # Set Axis-1 and label settings.
                ax4.set_ylabel("P value")
                plt.tick_params(axis="y", labelsize=9)

                # Create Axis-5, Standard Errors
                # Don't share the x-axis on this one, we want the tickvalue labels to show here.
                ax5 = plt.subplot(gs[4])
                ax5.plot(trace_numbers, v_stderrs, color="darkred")
                # Set Axis-1 and label settings.
                ax5.set_ylabel("Std. Error")
                plt.tick_params(axis="y", labelsize=9)

            # Adjust x-axis.
            (ax5 if plot_all_values else ax2).set_xlabel("Trace Number")
            plt.setp((ax5 if plot_all_values else ax2).get_xticklabels(), rotation=45, size=9)

            # Set overall figure properties
            fig.suptitle(self.name(), size=13)
            fig.set_size_inches(5.0,5.0)
            plt.tight_layout()
            fig.subplots_adjust(top=0.94)

            # Save the file.
            if self.verbose:
                print "Plotting", os.path.split(plotfilename)[-1]

            # Save at 300 dpi
            fig.savefig(plotfilename,dpi=300)

            # Clear the plotting axes and close the plot.
            plt.cla()
            plt.close()
            #
            # / if create_diagnostic_plots
            #

        # Now take the log traces and subtracts off gaussian-smoothed slopes for
        # each trace.  This should provide a "detrended" image to work with, which we can then plot.
        # To plot, we will need to borrow code from GPR_Data.py where we apply a color mask.
        if False: # FOR NOW, just skip this. Delete this line later.
            # Detrend each trace according to the calculated slopes
            for i in range(det_traces.shape[1]):
                if v_slopes_gaussian_neg_only[i] == 0.0:
                    det_traces[:,i] = log_traces[:,i]
                else:
                    det_traces[:,i] = log_traces[:,i] - (sample_numbers * v_slopes_gaussian_neg_only[i])

            if self.verbose:
                print "Writing", os.path.split(self.GPR_DETRENDED_FILENAME)[-1]

            f = open(self.GPR_DETRENDED_FILENAME, 'w')
            pickle.dump(det_traces, f)
            f.close()

        return self.GPR_DETRENDED_FILENAME

    def _read_coord_lines(self, strip_comments=True, check_resampled=True):
        '''Opens the coordinates file, returns lines.
        If "strip_comments", then strip out any blank lines or lines with comments.
        if "check_resampled", check to ensure that the file is the "resampled" 1.5-meter/trace version of the file.'''
        if check_resampled:
            path, fname = os.path.split(self.COORD_FILENAME)
            if fname.lower().find("_resampled") == -1 and os.path.split(path)[-1].find("resampled") == -1:
                raise ValueError("Must read 'resampled' coordinate files only: {0}".format(fname))

        f = open(self.COORD_FILENAME, 'r')
        lines = f.readlines()
        f.close()

        # Strip out all blank lines, and empty whitespace
        lines = [line.strip() for line in lines if len(line.strip()) > 0]

        # Strip comments if "strip_comments" set
        if strip_comments:
            lines = [line for line in lines if line[0] != "'"]

        return lines

    def extract_lats_lons_elevs(self, correct_lons=True):
        '''Return the latitude and longitude of each line in the file, each as an 1xN-length numpy array.
        If "correct_lons", turn westerly longitudes into negative numbers.'''
        lines = self._read_coord_lines(strip_comments=True, check_resampled=True)

        items = [line.split('\t') for line in lines]
        lats = numpy.array([float(line_items[3]) for line_items in items], dtype=numpy.float64)
        lons = numpy.array([(-float(line_items[5]) if (line_items[6] == "W" and correct_lons) else float(line_items[5])) for line_items in items], dtype=numpy.float64)
        elevs = numpy.array([float(line_items[7]) for line_items in items], dtype=numpy.float64)

        return lats, lons, elevs

    def closest_traces(self, lats, lons, topN=1):
        '''Given a latitude & longitude point (or list), find the closest GPR
        trace (or traces) to that point.
        Returns: (indices, distances)
            - indices: an MxN unsigned integer array,
                       giving the indices of the traces closest to that coordinate.
            - distances: an MxN float array,
                         giving the distances of each trace closest to that coordinate.
                         Traces are ordered by shortest distance first, in identical
                         order as the indices array
            M: the number of lat/lon points inputed.  Will be 1 if only one lat/lon value is entered.
            N: the "topN" value input.  Will be 1 if only the single closest value is wanted.
        '''
        # Convert lat, lon to utm_e, utm_n
        GCS2UTM = self._set_utm_converter()
        # Ensure the inputs are numpy arrays, convert to 1-element arrays if only single values
        if type(lats) in (float, int):
            lats = numpy.array([lats])
        if type(lons) in (float, int):
            lons = numpy.array([lons])

        # Make sure the longitude is negative (Western hemisphere)
        lons = (lons*-1 if numpy.all(lons > 0.0) else lons)

        points = zip(lons, lats)
        utm_points = numpy.array(GCS2UTM.TransformPoints(points))
        core_utm_es, core_utm_ns = utm_points.T[0,:], utm_points.T[1,:]

        # Get the lat/lon of all GPR traces in this file (from the coordinates file)
        # Turn westerly longitudes negative ("correct_lons")
        gpr_lats, gpr_lons, _ = self.extract_lats_lons_elevs(correct_lons=True)

        gpr_points = zip(gpr_lons, gpr_lats)
        gpr_utm_points = numpy.array(GCS2UTM.TransformPoints(gpr_points))
        gpr_utm_es, gpr_utm_ns = gpr_utm_points.T[0,:], gpr_utm_points.T[1,:]

        output_indices = numpy.empty((len(core_utm_es), topN), dtype=numpy.uint32)
        output_distances = numpy.empty((len(core_utm_es), topN), dtype=numpy.float)

        # rescale shapes of arrays to allow broadcasting along each axis.
        core_utm_es.shape = len(core_utm_es), 1
        core_utm_ns.shape = len(core_utm_ns), 1
        gpr_utm_es.shape  = 1, len(gpr_utm_es)
        gpr_utm_ns.shape  = 1, len(gpr_utm_ns)

        # For each core... if we can parallelize this loop, we'll try
        dist_squared = (gpr_utm_es - core_utm_es)**2 + (gpr_utm_ns - core_utm_ns)**2
        if topN == 1:
            ## NOTE: THIS FUNCTION HAS NOT YET BEEN TESTED FOR topN == 1
            # IF we're just looking for the single top value, let's make this quick.
            output_indices[:,0] = numpy.argmin(dist_squared, axis=1)
            output_indices.shape = [output_indices.shape[0]]
            output_distances[:,0] = numpy.sqrt(dist_squared[numpy.arange(0,len(core_utm_es),1,dtype=numpy.uint), output_indices])
            output_distances.shape = [output_distances.shape[0]]
        else:
            # Use numpy.partition to make this quicker than sorting ALL the distances.

            # create indices along the axis-1
            distance_partition_i1 = numpy.argpartition(dist_squared, topN, axis=1)[:,0:topN]
            # artificially create indices along axis-0 (should be easy but isn't yet)
            distance_partition_i0 = numpy.zeros(output_indices.shape, dtype=numpy.uint)
            i0_vector = numpy.arange(0,len(core_utm_es),1,dtype=numpy.uint)
            i0_vector.shape = i0_vector.shape[0], 1
            distance_partition_i0[:] = i0_vector[:]
            # Now copy the values into the distance_array
            distance_partition = numpy.sqrt(dist_squared[distance_partition_i0, distance_partition_i1])

            # Now, sort the (much smaller) list of distance values
            output_indices_i = numpy.argsort(distance_partition, axis=1)
            output_indices = distance_partition_i1[distance_partition_i0, output_indices_i]
            output_distances = distance_partition[distance_partition_i0, output_indices_i]

        return output_indices, output_distances

    def create_boolean_radar_colormap(self):
        '''When picking out boolean ice layers, color any pixels identified (True)
        in "bool_mask" as blue or purple.  Plot the rest of the image in greyscale.'''
        # Get a copy of the traces array, so we can adjust it.

        colortable = [(0.0, "blue"), (0.5/7.0, "blue"), (1.0/7.0,'#dddddd'), (3.0/7.0, '#dddddd'), (1.0,'#555555')]

        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("radar", colortable, N=1000)
        return cmap

    def export_boolean_ice_array(self, filetype="detrended", export_image=True):
        if self.boolean_ice_array is None:
            self.calculate_boolean_icelens_array()

        # Get a copy of traces array so we don't muck up the original.
        traces = self._open(filetype=filetype).copy()

        if self.boolean_ice_array is not None:
            bool_array = self.boolean_ice_array
        else:
            bool_array = self.calculate_boolean_icelens_array()

        if export_image:
            # Unforunately PyPlot doesn't like to save any image more than 32768 pixels wide
            if traces.shape[1] > 32768:
                traces = traces[:,:32768]
                bool_array = bool_array[:,:32768]

        minval,maxval = numpy.min(traces), numpy.max(traces)

        assert minval > 0.0
        # Set any pixel that's identified as ice to 0.0.  That will be our cutoff.
        traces[bool_array] = 0.0

        # Adjust the maximum to a known number so that all the values are scaled appropriately.
        if maxval > 7.0:
            # If some traces are above 7.0, scale to 7.0
            traces[traces > 7.0] = 7.0
        elif maxval < 7.0:
            # If all traces are below 7.0, set so that at least one pixel (the maximum one) is 7.0.
            # This is NOT to affect the data itself, just for display purposes so all
            # radar images are displayed on the same color map.
            traces[numpy.unravel_index(numpy.argmax(traces),traces.shape)] = 7.0
        maxval = 7.0

        if export_image:
            filebase = os.path.split(os.path.splitext((self.GPR_FILENAME if filetype == "logvariance" else self.GPR_DETRENDED_FILENAME))[0])[-1]
            outfilename = os.path.join(EXPORT_DIR, filebase + "_ice_lenses.png")

        cmap = self.create_boolean_radar_colormap()

        if export_image:
            if self.verbose:
                print "Saving", os.path.split(outfilename)[-1],

            plt.imsave(fname=outfilename, arr=traces, cmap=cmap)

            if self.verbose:
                print "Saved."

            plt.close()

        return traces, cmap

    def shrink_and_grow_boolean_array(self, bool_array, pixels_shrink_by=10, pixels_grow_by=10):
        '''In order to filter out noise in the radar data, shrink the radar boolean array
        by "pixels_shring_by" pixels, which will eliminate errant or stray "bits" of ice lenses.
        Then, grow the remaining ice lenses back by "pixels_grow_by".  The two numbers can
        be different if we wish to adjust how thick the lenses should be detected.

        Return the modified array.'''

        b = bool_array
        ydim, xdim = b.shape
        ymax = ydim-1
        xmax = xdim-1

        # Shrink the pixels
        if self.verbose:
            print "Shrinking {0} pixels.".format(pixels_shrink_by)
        for _ in range(pixels_shrink_by):
            b_temp = b.copy()
            # Top left pixel, turn false if it's already false or an adjacent pixel is false.  Else true.
            b_temp[0,0] = b[0,0] & b[0,1] & b[1,0]
            # Top right pixel, turn false if it's already false or an adjacent pixel is false.  Else true.
            b_temp[0,xmax] = b[0,xmax] & b[0,xmax-1] & b[1,xmax]
            # Bottom left pixel, turn false if it's already false or an adjacent pixel is false.  Else true.
            b_temp[ymax,0] = b[ymax,0] & b[ymax-1,0] & b[ymax,1]
            # Bottom right pixel, turn false if it's already false or an adjacent pixel is false.  Else true.
            b_temp[ymax,xmax] = b[ymax,xmax] & b[ymax,xmax-1] & b[ymax-1,xmax]
            # Top row
            b_temp[0,1:xmax] = b[0,1:xmax] & b[1,1:xmax] & b[0,0:xmax-1] & b[0,2:xmax+1]
            # Bottom row
            b_temp[ymax,1:xmax] = b[ymax,1:xmax] & b[ymax-1,1:xmax] & b[ymax,0:xmax-1] & b[ymax,2:xmax+1]
            # Left column
            b_temp[1:ymax,0] = b[1:ymax,0] & b[1:ymax,1] & b[0:ymax-1,0] & b[2:ymax+1,0]
            # Right column
            b_temp[1:ymax,xmax] = b[1:ymax,xmax] & b[1:ymax,xmax-1] & b[0:ymax-1,xmax] & b[2:ymax+1,xmax]
            # All the center pixels, must compare with self and 4 adjacent pixels
            b_temp[1:ymax,1:xmax] = b[1:ymax,1:xmax] & b[0:ymax-1,1:xmax] & b[2:ymax+1,1:xmax] & b[1:ymax,0:xmax-1] & b[1:ymax,2:xmax+1]
            b = b_temp

        # Grow the pixels
        if self.verbose:
            print "Growing {0} pixels.".format(pixels_grow_by)
        for _ in range(pixels_grow_by):
            b_temp = b.copy()
            # Top left pixel, turn true if it's already true or an adjacent pixel is true.  Else False.
            b_temp[0,0] = b[0,0] | b[0,1] | b[1,0]
            # Top right pixel, turn true if it's already true or an adjacent pixel is true.  Else False.
            b_temp[0,xmax] = b[0,xmax] | b[0,xmax-1] | b[1,xmax]
            # Bottom left pixel, turn true if it's already true or an adjacent pixel is true.  Else False.
            b_temp[ymax,0] = b[ymax,0] | b[ymax-1,0] | b[ymax,1]
            # Bottom right pixel, turn true if it's already true or an adjacent pixel is true.  Else False.
            b_temp[ymax,xmax] = b[ymax,xmax] | b[ymax,xmax-1] | b[ymax-1,xmax]
            # Top row
            b_temp[0,1:xmax] = b[0,1:xmax] | b[1,1:xmax] | b[0,0:xmax-1] | b[0,2:xmax+1]
            # Bottom row
            b_temp[ymax,1:xmax] = b[ymax,1:xmax] | b[ymax-1,1:xmax] | b[ymax,0:xmax-1] | b[ymax,2:xmax+1]
            # Left column
            b_temp[1:ymax,0] = b[1:ymax,0] | b[1:ymax,1] | b[0:ymax-1,0] | b[2:ymax+1,0]
            # Right column
            b_temp[1:ymax,xmax] = b[1:ymax,xmax] | b[1:ymax,xmax-1] | b[0:ymax-1,xmax] | b[2:ymax+1,xmax]
            # All the center pixels, must compare with self and 4 adjacent pixels
            b_temp[1:ymax,1:xmax] = b[1:ymax,1:xmax] | b[0:ymax-1,1:xmax] | b[2:ymax+1,1:xmax] | b[1:ymax,0:xmax-1] | b[1:ymax,2:xmax+1]
            b = b_temp

        return b

    def calculate_boolean_icelens_array(self, starting_cutoff=5.0, shrink_and_grow=False, skip_calculations=True):
        '''Using values calculated earlier in InSituGPR_Manager::calculate_ice_lenses(),
        we start with a cutoff value and use that the create a boolean ice mask.
        Then (if requested), we filter that mask and see if we can filter out noise while expanding signal.
        See how well this goes.'''
        # If we've already done these calculations and saved the results, no need to redo them.
        if skip_calculations and self.boolean_ice_array is not None:
            return self.boolean_ice_array

        traces = self._open(filetype="detrended")
        # First round, just take traces less than the cutoff
        boolean_array = (traces <= starting_cutoff)

        # This works (sort of), but not great.  Will leave out for further analysis, try other approaches.
        # Only perform this operation if "shrink_and_grow" is set (default to False)
        if shrink_and_grow:
            boolean_array = self.shrink_and_grow_boolean_array(boolean_array, pixels_shrink_by=10, pixels_grow_by=15)

        self.boolean_ice_array = boolean_array

        return boolean_array

    def calculate_trace_depths(self, density_kg_m3=873.0):
        #  Use the Robin's
        # coefficient of 0.734 previously calculated above for time/depth conversion.
        speed_picker = RadarSpeedPicker()
        speed_m_s = speed_picker.modified_Robin_speed(density_kg_m3, 0.734)

        # Trace samples are taken every 0.1 nano-second.
        trace_twt_increment = 0.10e-9

        # Calculate the approximate depth of each GPR sample
        trace_depths = numpy.arange(1,2025,1)*trace_twt_increment * speed_m_s / 2.0
        return trace_depths


    def calculate_ice_content_along_track(self, max_depth_m, density_kg_m3=873.0):
        '''Return an Nx4 numpy array of [lat,lon,elev,ice_content_percent] for all
        N traces of the radar track, within the top "max_depth_m" traces.

        Use a firn density of bubbly refrozen ice here (from Machguth, et al. 2016) to
        calculate depths and thicknesses.  This will under-estimate depth somewhat
        where there is firn, but will give more accurate results of actual ice thickness
        within those layers.'''
        boolean_array = self.calculate_boolean_icelens_array()

        # Calculate depths of traces down to "max_depth_m".
        trace_depths = self.calculate_trace_depths(density_kg_m3)
        # Cut off depths > max_depth_m
        trace_depths = trace_depths[trace_depths < max_depth_m]
        # Cut off bottom of boolean array that's too deep
        boolean_array = boolean_array[:len(trace_depths), :]

        # Calculate the % of ice in each trace of the array.
        ice_pct_array = numpy.sum(boolean_array, axis=0) / float(len(trace_depths)) * 100.0

        # Create the output array
        lat_lon_elev_pct = numpy.empty((len(ice_pct_array), 4), dtype=numpy.float64)
        # Get latitudes, longitudes, elevations
        lats, lons, elevs = self.extract_lats_lons_elevs()

        # Fill in array.
        lat_lon_elev_pct[:,0] = lats[:]
        lat_lon_elev_pct[:,1] = lons[:]
        lat_lon_elev_pct[:,2] = elevs[:]
        lat_lon_elev_pct[:,3] = ice_pct_array[:]

        return lat_lon_elev_pct

    def csv_of_ice_content(self, max_depth_m):
        '''For this track, spit out a CSV file of ice content in the top "max_depth_m" meters.
        Output the file into a CSV file in a folder.'''
        filename = os.path.join(FIGSHARE_BASE, "IceBridge_Accumulation_Radar", "txt", r"{0}_{1:0.1f}.csv".format(self.name(), max_depth_m))

        lat_lon_elev_pct = self.calculate_ice_content_along_track(max_depth_m=max_depth_m)
        csv_lines = ["{0:f},{1:f},{2:f},{3:f}".format(lat_lon_elev_pct[i,0], \
                                                      lat_lon_elev_pct[i,1], \
                                                      lat_lon_elev_pct[i,2], \
                                                      lat_lon_elev_pct[i,3]) \
                     for i in range(lat_lon_elev_pct.shape[0])]
        header_line = "latitude,longitude,elevation,pct_ice_in_top_{0:0.1f}_m\n".format(max_depth_m)

        f = open(filename, 'w')
        if self.verbose:
            print "Writing", os.path.split(filename)[-1]
        f.write(header_line + "\n".join(csv_lines))
        f.close()

    def print_minmax(self, filetype="detrended"):
        traces = self._open(filetype=filetype)
        print self.name(), numpy.min(traces), numpy.max(traces)

def plot_ice_content_vs_distance():
    '''Along the ACT-13 transect, plot the ice content from this lenses (m) in the
    top 10 meters of firn against the distance (from the ice border?) and the elevation.

    We already have CSV files output by the InSituGPR_Manager::create_csv_files_of_ice_content() function,
    use the output from those to plot ice content vs. elevtion for the ACT-transect.  Should be pretty quick.'''

    csv_name = r'C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\Ice Content CSVs\ACT_Transect_5.0.csv'
    fig_name = os.path.join(os.path.split(csv_name)[0], "ACT_Ice_Content_5.0_distance.tif")

    f = open(csv_name, 'r')
    lines = f.readlines()
    f.close()

    header = dict([(colname, i) for i, colname in enumerate(lines[0].strip().split(','))])
    lines = lines[1:]

    items = [line.strip().split(',') for line in lines]
    distance_km = numpy.arange(0,len(items)*1.5,1.5)/1000.0
    pct_ice = numpy.array([float(item[header["pct_ice_in_top_5.0_m"]]) for item in items])

    # Smooth each with a Gaussian filter.
    distance_smoothed = scipy.ndimage.gaussian_filter(distance_km, 1500/1.5)
    pct_ice_smoothed = scipy.ndimage.gaussian_filter(pct_ice, 1500/1.5)

    fig = plt.figure(dpi=600)
    fig.set_size_inches(12,1.65)
    plt.plot(distance_km, pct_ice * 5.0 / 100.0, color="lightgrey")
    plt.plot(distance_smoothed, pct_ice_smoothed * 5.0 / 100.0, color="blue", lw=1.5)
    plt.xlabel("Distance (km)")
    plt.ylabel("Ice (m)")
    plt.xlim(0,107)
    plt.ylim(0,4.5)
    plt.yticks(numpy.arange(0.,5.,1.0))
    plt.tight_layout()
    print "Saving", os.path.split(fig_name)[-1]
    plt.savefig(fig_name, dpi=600)
    plt.cla()
    plt.close()


def plot_ice_content_vs_elevation():
    '''Along the ACT-13 transect, plot the ice content from this lenses (m) in the
    top 10 meters of firn against the distance (from the ice border?) and the elevation.

    We already have CSV files output by the InSituGPR_Manager::create_csv_files_of_ice_content() function,
    use the output from those to plot ice content vs. elevtion for the ACT-transect.  Should be pretty quick.'''

    csv_name = r'C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\Ice Content CSVs\ACT_Transect_10.0.csv'
    fig_name = os.path.join(os.path.split(csv_name)[0], "ACT_Ice_Content_10.0.png")

    f = open(csv_name, 'r')
    lines = f.readlines()
    f.close()

    header = dict([(colname, i) for i, colname in enumerate(lines[0].strip().split(','))])
    lines = lines[1:]

    items = [line.strip().split(',') for line in lines]
    elevations = numpy.array([float(item[header["elevation"]]) for item in items])
    pct_ice = numpy.array([float(item[header["pct_ice_in_top_10.0_m"]]) for item in items])

    # Smooth each with a Gaussian filter.
    elev_smoothed = scipy.ndimage.gaussian_filter(elevations, 1500/1.5)
    pct_ice_smoothed = scipy.ndimage.gaussian_filter(pct_ice, 1500/1.5)

    plt.figure(figsize=(5,2.5), dpi=150)
    plt.plot(elevations, pct_ice * 10.0 / 100.0, color="lightgrey")
    plt.plot(elev_smoothed, pct_ice_smoothed * 10.0 / 100.0, color="blue", lw=1.5)
    plt.xlabel("Elevation (m a.s.l)")
    plt.ylabel("Ice content (m)")
    plt.tight_layout()
    plt.xlim(1700, 2260)
    plt.ylim(0,9)
    print "Saving", os.path.split(fig_name)[-1]
    plt.savefig(fig_name, dpi=150)
    plt.cla()
    plt.close()

def output_shapefile_first_45km():
    '''Opens one of the CSV files, outputs a simple linestring shapefile with the first
    45 km of transect points. Just uses every 100th point (300 total) because we don't
    need that much detail. The shapefile is useful for Figure 1.1 of the paper.'''
    csv_name = r'C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\Ice Content CSVs\ACT_Transect_10.0.csv'
    shapefile_name = r'C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Data\ACT_Transect_Shapefiles\ACT_Transect_45km.shp'

    f = open(csv_name, 'r')
    lines = f.readlines()
    f.close()

    header = dict([(colname, i) for i, colname in enumerate(lines[0].strip().split(','))])
    lines = lines[1:]
    items = [line.strip().split(',') for line in lines]

    lats = numpy.array([float(item[header["latitude"]]) for item in items])
    lons = numpy.array([float(item[header["longitude"]]) for item in items])
    distance_km = numpy.arange(0,len(items)*1.5,1.5)/1000.0

    # Subset to first 45 kilometers
    mask_45km = distance_km <= 45.0
    lats = lats[mask_45km]
    lons = lons[mask_45km]
    distance_km = distance_km[mask_45km]
    # now subset to every 100th point
    mask_100th = numpy.arange(0,len(lats),100,dtype=numpy.long)
    # If it doesn't include the "last" point, include it
    if mask_100th[-1] < (len(lats)-1):
        mask_100th = numpy.append(mask_100th, [len(lats)-1])
    lats = lats[mask_100th]
    lons = lons[mask_100th]
    distance_km = distance_km[mask_100th]

    # Create the Shapefile data source and file.
    driver = ogr.GetDriverByName("ESRI Shapefile")

    # If the shapefile already exists, delete it.
    if os.path.exists(shapefile_name):
        print "Deleting previous", os.path.split(shapefile_name)[-1]
        driver.DeleteDataSource(shapefile_name)

    # Create the shapefile
    data_source = driver.CreateDataSource(shapefile_name)

    # Create a WGS84 geographic spatial reference to go with it.
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)

    # Create the layer
    layer = data_source.CreateLayer("ACT13", srs, geom_type=ogr.wkbLineString)

    # Create the line feature:
    feature = ogr.Feature(layer.GetLayerDefn())

    wkt_geom_string = "LINESTRING (" + \
                      ",".join(["{0} {1}".format(lon, lat) for (lat, lon) in zip(lats, lons)]) \
                      + ")"

    linestring = ogr.CreateGeometryFromWkt(wkt_geom_string)
    feature.SetGeometry(linestring)
    # Put this feature in the layer
    layer.CreateFeature(feature)
    # Dereference the feature and data source to save the shapefile.
    feature = None
    data_source = None
    print os.path.split(shapefile_name)[-1], "written."

    return

if __name__ == "__main__":

    cores = FirnCore_Manager(verbose=True)

    gpr = InSituGPR_Manager(verbose=True)

    print gpr.track_names
    gpr.determine_best_radar_coefficient_by_correlation(cores, plot_figures=True)
    gpr.calculate_ice_lenses(cores, make_all_plots=True)
    gpr.export_boolean_icelens_images()

    gpr.apply_vertical_variance_corrections()
    gpr.export_radar_images(filetype="logvariance")
    gpr.export_radar_images(filetype="detrended")

    gpr.create_csv_files_of_ice_content(max_depth_m=20.0)

    output_shapefile_first_45km()
