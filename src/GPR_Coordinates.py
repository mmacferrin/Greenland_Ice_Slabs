# -*- coding: utf-8 -*-
"""
Created on Fri Oct 31 21:32:01 2014

GPR_Coordinates.py
* Works with post-processed GPR Coordinates
* Outputs interpolated lines with gaps filled in
* Maintains coordinates both in lat/lon and UTM
* Finds closest traces to each core from each GPR profile.

@author: Mike MacFerrin
"""

from GPR_FileData import RESAMPLED_COORDINATE_FILES, \
                         RESAMPLED_COR_DIR,

import numpy as np
from osgeo import osr
import os
import matplotlib.pyplot as plt
import datetime

class GPR_Coords:
    def __init__(self, GPR_CORFILE):
        self.filename = GPR_CORFILE
        self.datetimes = None
        self.lats = None
        self.lons = None
        self.utm_n = None
        self.utm_e = None
        self.elevs = None
        self.utm_converter = None
        # If the file is already fully interpreted with no gaps, this simply reads the data
        self.interpolate_lines()

    def set_utm_converter(self, utmZone=23, north=1):
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

    def interpolate_lines(self):
        '''Reads the file, fills in all the arrays and interpolates all missing values,
        usually not more than 1 in any given gap.'''
        f = open(self.filename, 'r')
        lines = f.readlines()
        f.close()

        # Strip off potential empty lines at end
        while len(lines[-1].strip()) == 0:
            lines = lines[:-1]
        # Get size of array by getting the last trace-num (may underestimate, which is okay)
        try:
            N = int(lines[-1].split()[0])
        except ValueError:
            N = int(lines[-1].split(',')[0])

        datetimes = [None] * N
        lats =  np.zeros([N], dtype=np.float64)
        lons =  np.zeros([N], dtype=np.float64)
        utm_n = np.zeros([N], dtype=np.float64)
        utm_e = np.zeros([N], dtype=np.float64)
        elevs = np.zeros([N], dtype=np.float64)

        last_i = 0
        for line in lines:
            if line.lstrip()[0] == "'": # filter out comments
                continue
            items = line.split()

            assert len(items) >= 10
            this_i = int(items[0])-1
            datestr, timestr = items[1], items[2]
            try:
                datetimes[this_i] = datetime.datetime.strptime(datestr + " " + timestr, "%Y-%m-%d %H:%M:%S")
            except ValueError:
                datetimes[this_i] = datetime.datetime.strptime(datestr + " " + timestr[:-3], "%Y-%m-%d %H:%M:%S")
                datetimes[this_i] = datetimes[this_i] + datetime.timedelta(milliseconds=int(float(timestr[-3:])*1e3))
            lats[this_i] = float(items[3])
            lons[this_i] = float(items[5])
            elevs[this_i] = float(items[7])

            # Interpolate over any missing values, assume constant speed
            for i in range(last_i+1, this_i):
                lats[i] = (lats[this_i] - lats[last_i]) * float(i-last_i)/float(this_i-last_i) + lats[last_i]
                lons[i] = (lons[this_i] - lons[last_i]) * float(i-last_i)/float(this_i-last_i) + lons[last_i]
                elevs[i] = (elevs[this_i] - elevs[last_i]) * float(i-last_i)/float(this_i-last_i) + elevs[last_i]
                datetimes[i] = ((datetimes[this_i] - datetimes[last_i])*(i-last_i))/(this_i-last_i) + datetimes[last_i]

            last_i = this_i

        # Assure that no values were skipped
        assert np.count_nonzero(lats) == np.size(lats)

        # Convert to UTM zone 23N
        GCS2UTM = self.set_utm_converter()

        # Must convert positive longitudes (*W) to negative
        utm_points = GCS2UTM.TransformPoints(zip(-lons, lats))
        for i in range(N):
            utm_e[i] = utm_points[i][0]
            utm_n[i] = utm_points[i][1]

        self.datetimes = datetimes
        self.lats = lats
        self.lons = lons
        self.utm_e = utm_e
        self.utm_n = utm_n
        self.elevs = elevs
        return

    def distances(self):
        '''Return an N-1 length list of distances between interpolated coordinate points.'''
        return np.sqrt((self.utm_e[1:] - self.utm_e[:-1])**2 + (self.utm_n[1:] - self.utm_n[:-1])**2)

    def plot_trace_distances(self):

        distances = self.distances()
        print "----", os.path.split(self.filename)[-1], "----"
        print "mean:", np.mean(distances)
        print "min:", np.min(distances), "  max:", np.max(distances)

        threshold = 5.0
        count_over_threshold = np.count_nonzero(distances > threshold)
        if count_over_threshold > 0:
            print count_over_threshold, "traces above", threshold, "m/trace."
            print np.where(distances > threshold)[0] + 1
        print


        plt.subplot(2,1,1)
        plt.hist(distances, bins=40)
        plt.title(os.path.split(self.filename)[1])

        plt.subplot(2,1,2)
        plt.plot(range(1,len(distances)+1), distances)
        plt.xlabel("tracenum")
        plt.ylabel("distance per trace (m)")

        plt.show()
        plt.close()

    def closest_traces(self, latlon, numtraces=2):
        '''Given a lat/lon tuple, return the "numtraces" closest traces, as well
        as distances in meters.  Traces will be returned as a list of tuples, each
        tuple containing (gpr_filename, tracenum, lat, lon, utm_n, utm_e, elev, distance_in_m).'''
        if (self.lats is None) or (self.utm_e is None) or (self.utm_n is None):
            self.interpolate_lines()

        # Convert lat, lon to utm_e, utm_n
        GCS2UTM = self.set_utm_converter()
        lat, lon = latlon
        # Make sure the longitude is negative (Western hemisphere)
        lon = (-lon if lon > 0 else lon)
        point_utm_e, point_utm_n, _ = GCS2UTM.TransformPoint(lon, lat)

        # Find distances to point from coordinates:
        distances = np.sqrt((self.utm_e - point_utm_e)**2 + (self.utm_n - point_utm_n)**2)
        sorted_indices = np.argsort(distances)[:numtraces]

        tracenums = sorted_indices + 1
        lats = self.lats[sorted_indices]
        lons = self.lons[sorted_indices]
        utm_ns = self.utm_n[sorted_indices]
        utm_es = self.utm_e[sorted_indices]
        elevs = self.elevs[sorted_indices]
        distances = distances[sorted_indices]

        return_list = [None] * numtraces
        for i in range(numtraces):
            return_list[i] = (self.filename, tracenums[i], lats[i], lons[i],
                              utm_ns[i], utm_es[i], elevs[i], distances[i])

        return return_list


    def resample(self, interval=1.5, outcorfile=None, verbose=True):
        # Perform a nearest-neighbor resampling at "interval" spans (in m).
        # Output a new .cor file corresponding to this new list of resampled traces.
        # The new corfile is written (if specified), but the points in this object
        # remain the same.  To create a new object, simply create a new GPR_Coords object
        # with the new .cor file as its base.
        #
        # Return a numpy int-array of tracenumbers (from the original file) to put
        # into a new image file.
        #
        # If "verbose", print out how many samples were skipped and how many were repeated.
        if verbose:
            print os.path.split(self.filename)[-1]

        old_distances = self.distances()
        cum_distances = np.empty((len(old_distances)+1,), dtype=np.double)
        cum_distances[0] = 0.0
        for i in range(len(old_distances)):
            cum_distances[i+1] = cum_distances[i] + old_distances[i]

        total_distance = np.sum(old_distances)
        numtraces = int(np.floor(total_distance / interval))

        new_lats  = np.empty((numtraces,), dtype=np.double)
        new_lons  = np.empty((numtraces,), dtype=np.double)
        new_elevs = np.empty((numtraces,), dtype=np.double)
        new_distances = np.arange(0.0, interval*numtraces, interval, dtype=np.double)
        new_timestamps = [None] * numtraces
        old_tracenums = np.empty((numtraces,), dtype=np.int)

        # "old_j" refers to the counting index of the old traces
        # "i" refers to the index into he new traces
        old_j = 0
        for i in range(len(new_lats)):
            distance = new_distances[i]
            while cum_distances[old_j] < distance:
                old_j += 1

            if cum_distances[old_j] == distance:
                new_lats[i] = self.lats[old_j]
                new_lons[i] = self.lons[old_j]
                new_elevs[i] = self.elevs[old_j]
                new_timestamps[i] = self.datetimes[old_j]
                old_tracenums[i] = (old_j + 1)
            elif cum_distances[old_j] > distance:
                ratio_along_line = (distance - cum_distances[old_j-1])/(cum_distances[old_j] - cum_distances[old_j-1])
                assert (0 < ratio_along_line < 1)
                new_lats[i] = (self.lats[old_j] - self.lats[old_j-1])*ratio_along_line + self.lats[old_j-1]
                new_lons[i] = (self.lons[old_j] - self.lons[old_j-1])*ratio_along_line + self.lons[old_j-1]
                new_elevs[i] = (self.elevs[old_j] - self.elevs[old_j-1])*ratio_along_line + self.elevs[old_j-1]
                # Interpolate timestamp here, using ms.
                old_timedelta = self.datetimes[old_j] - self.datetimes[old_j-1]
                # Have to use integer math, will calculate an approximate fraction here.
                new_timedelta = old_timedelta * int(ratio_along_line * 1e6) / int(1e6)
                new_timestamps[i] = self.datetimes[old_j-1] + new_timedelta

                old_tracenums[i] = old_j if (ratio_along_line >= 0.5) else (old_j+1)

            else:
                assert False

        if verbose:
            # Check for skipped or repeated numbers here.
            iters = old_tracenums[1:] - old_tracenums[:-1]
            num_repeated = np.count_nonzero(iters == 0)
            num_skipped = np.sum(iters[iters > 1]) - np.count_nonzero(iters > 1)
            print num_skipped, "traces skipped. ", num_repeated, "traces repeated before fixing."

        # The synchronization artifacts with the GPS leave a number of traces repeated immediately adjacent
        # to other traces that are skipped, which shouldn't happen.  This will go through all the
        # old_tracenums and adjust individual ones so a greater number of unique traces are included without
        # unnecessarily skipping adjacent ones.
        shifted_forward = 0
        shifted_backward = 0
        for i in range(1, len(old_tracenums)-1):
            if (old_tracenums[i] == old_tracenums[i-1]) and (old_tracenums[i+1]-old_tracenums[i] >= 2):
                old_tracenums[i] = old_tracenums[i] + 1
                shifted_forward += 1
            elif (old_tracenums[i] - old_tracenums[i-1] >= 2) and (old_tracenums[i+1] == old_tracenums[i]):
                old_tracenums[i] = old_tracenums[i] - 1
                shifted_backward += 1
        if verbose:
            print shifted_forward, "points shifted forward,", shifted_backward, "points shifted backward."
            iters = old_tracenums[1:] - old_tracenums[:-1]
            num_repeated = np.count_nonzero(iters == 0)
            num_skipped = np.sum(iters[iters > 1]) - np.count_nonzero(iters > 1)
            print num_skipped, "traces skipped. ", num_repeated, "traces repeated after fixing."

        if outcorfile != None:
            if verbose: print "Writing resampled\\" + os.path.split(outcorfile)[-1], '\n'
            f = open(outcorfile, 'w')
            # Print out a new .COR file.
            outstr = "{0:d}\t{1:04d}-{2:02d}-{3:02d}\t{4:02d}:{5:02d}:{6:02d}.{7:02d}\t{8:14.11f}\tN\t{9:14.11f}\tW\t{10:8.3f}\tM\t1\n"

            for i in range(len(new_lats)):
                dt = new_timestamps[i]
                f.write(outstr.format(i+1, dt.year, dt.month, dt.day, dt.hour, \
                                      dt.minute, dt.second, (dt.microsecond)/10000, \
                                      new_lats[i], new_lons[i], new_elevs[i]))

            f.close()

        return old_tracenums


def plot_orig_and_resampled_utm_coords(old_coords, new_coords):
    '''Given two coordinate files, one original and one resampled, plot a few of them
    to visually check how the resampling went.'''
    assert isinstance(old_coords, GPR_Coords)
    assert isinstance(new_coords, GPR_Coords)
    assert len(old_coords.utm_e) > 50 and len(new_coords.utm_e) > 50

    plt.plot(old_coords.utm_e[:1000], old_coords.utm_n[:1000], 'b-')
    plt.scatter(old_coords.utm_e[:1000], old_coords.utm_n[:1000], c='b', marker='o')
    plt.scatter(new_coords.utm_e[:1000], new_coords.utm_n[:1000], c='r', marker="^")
    plt.title(os.path.split(old_coords.filename)[-1])
    plt.show()


if __name__ == "__main__":

    for f in RESAMPLED_COORDINATE_FILES:

        g = GPR_Coords(f)
        g = GPR_Coords(f)
        new_f = os.path.join(RESAMPLED_COR_DIR, os.path.split(f)[-1])
        g.resample(interval=1.5,outcorfile=new_f,verbose=True)
