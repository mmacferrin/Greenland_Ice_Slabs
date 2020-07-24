# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 14:52:50 2016

@author: mmacferrin

IceBridge Flight Line Database -- for locating, reading, subsetting and ingesting OIB data

This code is primarily responsible for creating an hdf5 metadata database for all the 
OIB Accumulation Radar flightlines. This allows code to open one database to get all the 
metadata (lat,lon,elev,roll,pitch,heading,etc) for each trace without having to
traverse the directory structure of the files and open each individual data file
just to pull out metadata on the traces. It takes a few minutes to build the database,
but saves significant processing time when doing the IceBridge data processing.

There is also a member function for outputting a shapefile (.shp) from the
IceBridge metadata, useful for GIS visualizations.

The primary function in this code is BUILD_DATABASE(), toward the bottom.

"""

import numpy
import os
import re
import scipy.io
import h5py
import tables as T
from osgeo import ogr, osr
# import matplotlib.pyplot as plt

# Stand-in variables for directories.
#BASEDIR = r'C:\Users\mmacferrin\Desktop'
BASEDIR = r'F:\Research\DATA\ATM-IceBridge\IceBridge - AccumRadar\accum'
DB_FILE = os.path.join(BASEDIR, 'IceBridgeDB.h5')
CC_LINE_SHAPEFILE = r'F:\Research\DATA\Camp Century, Greenland\CC_Icebridge_shapefile'

class IceBridgeFlightLineDescriptor(T.IsDescription):
    Flightline_ID       = T.UInt32Col(pos=0)
    name                = T.StringCol(itemsize=20, pos=1)
    year                = T.UInt16Col(pos=2)
    month               = T.UInt8Col(pos=3)
    day                 = T.UInt8Col(pos=4)

class IceBridgeFileTableDescriptor(T.IsDescription):
    File_ID             = T.UInt32Col(pos=0)
    Flightline_ID       = T.UInt32Col(pos=1) # Links to "FlightLine" table, see above
    relative_filepath   = T.StringCol(itemsize=256, pos=2)

class IceBridgeCoordinatesDescriptor(T.IsDescription):
    Flightline_ID       = T.UInt32Col(pos=0) # Links to "FlightLine" table, see above
    File_ID             = T.UInt32Col(pos=1) # Links to "File" table, see above
    Flight_tracenum     = T.UInt32Col(pos=2) # Tells which unique tracenumber in the entire flightline
    File_tracenum       = T.UInt32Col(pos=3) # Tells which unique tracenumber in the file
    Latitude            = T.Float64Col(pos=4)
    Longitude           = T.Float64Col(pos=5)
    Elevation           = T.Float64Col(pos=6)
    Roll                = T.Float64Col(pos=7)
    Pitch               = T.Float64Col(pos=8)
    Heading             = T.Float64Col(pos=9)
    GPS_time            = T.Float64Col(pos=10)
    Surface             = T.Float64Col(pos=11)
    OnGreenlandIce      = T.BoolCol(pos=12) # A boolean flag to tell whether it's on the Greenland ice sheet
    OnAntarcticIce      = T.BoolCol(pos=13) # A boolean flag to tell whether it's on the Antarctic ice sheet

class IceBridgeRadarDB:
    PyTableGroupName = "Accumulation"
    PyTableFlightLineTableName = "Flight_Lines"
    PyTableFileTableName = "File_Names"
    PyTableCoordinatesTableName = "Coordinates"

    def __init__(self, dbfile = DB_FILE, datadir = BASEDIR):
        self.filename          = dbfile
        self.basedir           = datadir
        self.h5file            = None
        self.flightline_table  = None
        self.file_table        = None
        self.coordinates_table = None

    def __del__(self):
        if self.h5file != None:
            self.h5file.close()

    def build_database(self):
        '''
        # Read all the data from the self.radar_basedir.
        #   - Create a db5 Pytables database
        #   - Peruse the appropriate subdirectories and pull out all radar locations.
        #   - Read those into the database tables.'''
        # This assigns variables self.h5file (object), self.flightline_table,
        #  self.file_table and self.coordinates_table
        self._create_new_database(self.filename)

        # Fill the datdabase tables
        self._populate_database(self.basedir)

        # Once the table is created and populated, put indices in place to optimize read times.
        self._create_table_indices()

    def open_database(self, status='r+'):
        '''Open the database with read permissions, populate the table object variables.'''
        if self.h5file != None and self.h5file.isopen:
            return
        if os.path.exists(self.filename):
            self.h5file = T.openFile(self.filename, status, title="IceBridge")
            self.flightline_table = self.h5file.getNode("/{0}/{1}".format(self.PyTableGroupName, self.PyTableFlightLineTableName))
            self.file_table = self.h5file.getNode("/{0}/{1}".format(self.PyTableGroupName, self.PyTableFileTableName))
            self.coordinates_table = self.h5file.getNode("/{0}/{1}".format(self.PyTableGroupName, self.PyTableCoordinatesTableName))
        else:
            raise IOError("File \"{0}\" does not exist.".format(self.filename))

        return

    def _create_new_database(self, fname):
        ''' Call this from build_database().
        Creates a new, empty database and the tables, returns the database object
        with write permissions, ready for data entry.'''

        self.h5file = T.openFile(fname, 'w', title="IceBridge")
        print('Writing file "{0}"'.format(fname))
        assert self.h5file.isopen

        # Create the base "group" of tables
        group = self.h5file.createGroup("/", self.PyTableGroupName, title="Accumulation")

        # Create the "Flightline" table
        self.flightline_table = self.h5file.createTable(group,
                                                  self.PyTableFlightLineTableName,
                                                  IceBridgeFlightLineDescriptor)
        self.flightline_table.attrs.TABLE_DESCRIPTION = "Describes each unique flight line in the OIB data.  Flight lines may be described by multiple files."

        # Create the "File" table
        self.file_table = self.h5file.createTable(group,
                                             self.PyTableFileTableName,
                                             IceBridgeFileTableDescriptor)
        self.file_table.attrs.TABLE_DESCRIPTION = "Describes the location of each IceBridge file, and the flight line it belongs to."

        # Create the "Coordinates" table
        self.coordinates_table = self.h5file.createTable(group,
                                             self.PyTableCoordinatesTableName,
                                             IceBridgeCoordinatesDescriptor)
        self.coordinates_table.attrs.TABLE_DESCRIPTION = "Lists XYZ coordinates of each trace, and which file and flight line it belongs to.  All trace numbers should be present in order in the table."

        return

    def _populate_database(self, basedir):
        '''First, traverse through all the subdirectories and accumulate files & flight lines
        Then, enter each file and flight line into each respective table.
        Then, enter each trace coordinate into the coordinates table.'''

        flightline_dict = self._recurse_directories(basedir)
        self._populate_flightline_table(flightline_dict)
        self._populate_file_table(flightline_dict)
        self._populate_coordinate_table(flightline_dict)

        return

    def _recurse_directories(self, basedir):
        '''Recursively traverse through all the sub-directories.
        Build a dictionary with key/value: (flightline: (python list of filenames)).
        Filenames will all be relative paths to the basedir.'''

        print("Recursing directories...")

        flightline_dict ={}

        Nfiles = 0
        file_list = []

        # First level directory, contains folders named YYYY_Greenland_P3 or YYYY_Antarctica_P3
        folder_list_1 = [os.path.join(basedir,f) for f in os.listdir(basedir) if \
                        (os.path.isdir(os.path.join(basedir,f)) and re.search(r'\d{4}_((Greenland)|(Antarctica))_P3', f) != None)]
        for year_folder in folder_list_1:
            # Look in "CSARP_qlook" folder
            year_folder_qlook = os.path.join(year_folder, "CSARP_qlook")
            # Look for all folders in there with a YYYYMMDD_FF numerical format.
            folder_list_2 = [os.path.join(year_folder_qlook, f) for f in os.listdir(year_folder_qlook) if \
                             (os.path.isdir(os.path.join(year_folder_qlook, f)) and re.search(r'\d{8}_\d{2}',f) != None)]

            # Look through each of those directories
            for flightline_folder in folder_list_2:
                # Collect all *.mat files in that directory
                datafile_list = [os.path.join(flightline_folder, f) for f in os.listdir(flightline_folder) if f[-4:] == ".mat"]
                Nfiles += len(datafile_list)
                file_list.extend(datafile_list)

#                N_total_files = len(os.listdir(flightline_folder))
#                N_left_out = N_total_files - len(datafile_list)
#
#                print os.path.split(flightline_folder)[-1], len(datafile_list), "files,", Nfiles, "total,", N_left_out, "left out."

        for fname in file_list:
            flightline = self._get_flightline_text_from_filename(os.path.split(fname)[-1])
            if flightline in list(flightline_dict.keys()):
                flightline_dict[flightline].append(fname)
            else:
                flightline_dict[flightline] = [fname]

#        print flightline_dict.keys(),
#        print len(flightline_dict[flightline_dict.keys()[0]])

        return flightline_dict

    def _compute_ice_flags(self):
        '''Compute two vectors with boolean flags for each coordinate, whether or not it's on the Greenland ice sheet
        and on the Antarctic ice sheet (uphill of the MOA grounding line), respectively.
        Return two vectors of boolean T/F flags for Greenland and Antartctica, respectively.'''
        if self.h5file is None:
            self.open_database()

        # Retrieve IceBridge latitudes & longitudes
        ib_lats = self.coordinates_table.cols.Latitude[:]
        ib_lons = self.coordinates_table.cols.Longitude[:]

        ############################################################################
        # Retrieve GREENLAND lat/lons
        gr_txt_file = r'F:\Research\DATA\ASTER GDEM Greenland Misc\Shapefiles\Eispolygon_GEO2.txt'
        f = open(gr_txt_file)
        lines = [line.split() for line in f.readlines() if line.strip() != '']
        f.close()
        # Convert all longitudes to positive-only.  See what this does (may need to offset and come back later).
        latlons = numpy.array([(float(item[0]), float(item[1])) for item in lines])
        gr_lats = latlons[:,0]
        gr_lons = latlons[:,1]
        # Seperate out polygons for peripheral ice caps.
        gr_polys = []
        temp_i = 0
        while temp_i < latlons.shape[0]:
            #Each polygon makes a loop, ending with the same coordinate upon which it begins.  Use this.
            lat_i, lon_i = gr_lats[temp_i], gr_lons[temp_i]
            match_i = numpy.where((gr_lats[temp_i:] == lat_i) & (gr_lons[temp_i:] == lon_i))[0]
            assert len(match_i) >= 2
            gr_polys.append( (gr_lats[temp_i : temp_i + match_i[1]+1], gr_lons[temp_i : temp_i + match_i[1]+1]) )
            temp_i = temp_i + match_i[1] + 1
        ############################################################################

        # CAN WE DRAMATICALLY SPEED THIS UP BY USING A BOUNDING BOX (especialy for Antarctica?)?

        # Initialize flags to zero.
        empty_flags = numpy.zeros(ib_lats.shape, dtype=numpy.bool)
        gr_poly_flags = [None] * len(gr_polys)

        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Define a reusable function inside this function to perform the polygon analysis.
        def inside_polygon_func(coordlats, coordlons, polylats, polylons):
            '''Return a boolean array of length coordlats that describes whether
            each point is inside the polygon.  Uses vertical lines going up, so
            it works on the Antarctica polygon we have, even if it's not technically a polygon.'''
            # First, let's use a bounding box to subset our data and make this quicker by doing computations on far fewer points.
            bbox_lats_min = numpy.min(polylats)
            # Must handle specially for Antarctica, put the minimum latitude at -91.
            # If it's a negative minimum latitude, we're talking about Antarctica
            if bbox_lats_min < -70:
                bbox_lats_min = -91.0

            bbox_lats_max = numpy.max(polylats)
            bbox_lons_min = numpy.min(polylons)
            bbox_lons_max = numpy.max(polylons)

            bbox_mask = (coordlats <= bbox_lats_max) & (coordlats >= bbox_lats_min) & (coordlons <= bbox_lons_max) & (coordlons >= bbox_lons_min)
            sub_coordlons = coordlons[bbox_mask].copy()
            sub_coordlats = coordlats[bbox_mask].copy()

            progress_i = 0

            flags = empty_flags.copy()

            if sub_coordlons.size == 0:
                return flags
            sub_flags = flags[bbox_mask].copy()

            for i in range(polylats.size-1):
                i_flags_1 = (polylons[i] > sub_coordlons) != (polylons[i+1] > sub_coordlons)

                # Subset again by the first set of conditions (in between the two adjacent longitudes)
                # This subset of points consists only of those that met the first condition.
                sub_coordlons_2 = sub_coordlons[i_flags_1]
                sub_coordlats_2 = sub_coordlats[i_flags_1]

                if sub_coordlats_2.size > 0:
                    # Test the second set of conditions, just with this reduced dataset.
                    i_flags_2 = (sub_coordlats_2 < ((polylats[i+1] - polylats[i]) * (sub_coordlons_2 - polylons[i]) / (polylons[i+1] - polylons[i]) + polylats[i]))

                    i_flags_1[i_flags_1] = i_flags_2

                    # if it crosses that line, switch the flags.
                    sub_flags[i_flags_1] = ~sub_flags[i_flags_1]

                # Output progress dots.
                progress_i += 1

                if (progress_i % 50) == 0:
                    print('.', end=' ')
                    if (progress_i % 1000) == 0:
                        print(progress_i)
                    if progress_i == (plats.size - 2):
                        print()

            flags[bbox_mask] = sub_flags

            return flags
        ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # We could vectorize all this, but NxM array would be too huge for memory.  Must iterate instead.
        for pi,(plats, plons) in enumerate(gr_polys):
            print("Poly", pi+1, ",", plats.size, "points.")

            gr_poly_flags[pi] = inside_polygon_func(ib_lats, ib_lons, plats, plons)

            print(numpy.sum(gr_poly_flags[pi]), "inside polygon.")

        # Compile into a single array with or operators.
        gr_flags = empty_flags.copy()
        for poly_flags in gr_poly_flags:
            gr_flags = gr_flags | poly_flags

        print("TOTAL GREENLAND")
        print(gr_flags.size, numpy.sum(gr_flags), float(numpy.sum(gr_flags)) / gr_flags.size * 100.0, "%")
        print()

        ############################################################################
        # Retrieve ANTARCTICA lat/lons
        an_sav_file = r'F:\Research\DATA\Boundaries\Antarctica\moa_groundingline.sav'
        sav = scipy.io.readsav(an_sav_file)
        lats = sav['olat']
        lons = sav['olon']
        # Find the index where it flips over the longitude cutoff... make that the start/end of the line.
        lon_gap_i = numpy.argmin(lons[1:] - lons[:-1]) + 1
        an_lons = numpy.append(lons[lon_gap_i:], lons[:lon_gap_i])
        an_lats = numpy.append(lats[lon_gap_i:], lats[:lon_gap_i])
        # For completeness, add the final point from each (disjointed) end
        # onto the other end, at the longitude cutoffs of +/- 180*
        # This will ensure there's not a "gap".  We will look for intersections looking northward, not east/west.
        an_lons = numpy.append([an_lons[-1] - 360.0], lons)
        an_lats = numpy.append([an_lats[-1]], an_lats)
        an_lons = numpy.append(an_lons, [an_lons[1] + 360.0])
        an_lats = numpy.append(an_lats, [an_lats[1]])
        # Now we have a single W->E line against which to run our analysis w/ S->N vectors.
        ############################################################################

        # initialize flags to zero
        print("processing Antarctica... this'll take awhile.")
        an_flags = inside_polygon_func(ib_lats, ib_lons, an_lats, an_lons)
        print("TOTAL ANTARCTICA")
        print(an_flags.size, numpy.sum(an_flags), float(numpy.sum(an_flags)) / an_flags.size * 100.0, "%")
        print()

        return gr_flags, an_flags

    def _get_flightline_text_from_filename(self, fname):
        match = re.search(r'(?<=Data_)(img_01_)?\d{8}_\d{2}(?=_\d{3}\.mat$)', fname)
        assert match != None
        match_str = match.group(0)
        if match_str[0:7] == "img_01_":
            return match_str[7:]
        else:
            return match_str

    def _get_flightline_number_from_flightline_text(self, flightline_str):
        return int(flightline_str[0:8] + flightline_str[-2:])

    def _populate_flightline_table(self, flightline_dict):
        '''Using the flightline dictionary, enter all the keys into the flightline table.'''
        print("Populating Flightline Table")

        i=0
        for flightline in list(flightline_dict.keys()):
            i+=1
            row = self.flightline_table.row
            row['Flightline_ID'] = self._get_flightline_number_from_flightline_text(flightline)
            row['name']  = flightline
            row['year']  = int(flightline[0:4])
            row['month'] = int(flightline[4:6])
            row['day']   = int(flightline[6:8])
            row.append()

        self.flightline_table.flush()
        print("\t", i, "rows written.")
        return

    def _populate_file_table(self, flightline_dict):
        '''Using the flightline dictionary, enter all the files in each value set
        into the file database.'''
        # This is an artificial iterator for the File_ID field, perhaps look
        # at ways to use an actual iterator with enforced uniqueness?
        print("Populating File Table")
        file_iterator = 0

        flightlines = list(flightline_dict.keys())
        flightlines.sort()

        for flightline in flightlines:
            flightline_id = self._get_flightline_number_from_flightline_text(flightline)

            for filename in flightline_dict[flightline]:
                # Generate a relative filename to the base directory.
                Nbasechars = len(self.basedir)
                assert filename[0:Nbasechars] == self.basedir
                relative_filepath = filename[Nbasechars+1:]

                row =self.file_table.row
                row['File_ID'] = file_iterator
                row['Flightline_ID'] = flightline_id
                row['relative_filepath'] = relative_filepath
                row.append()
                file_iterator += 1

        self.file_table.flush()
        print("\t", file_iterator, "rows written.")
        return

    def _convert_h5_dataset_to_dict(self, dataset, include_Data=False):
        '''The files read as h5py files don't return arrays when you index them.
        For this, we'll convert the h5py "Dataset" into a simple key/value dictionary.

        The "include_Data" flag tells whether or not to copy over the large "Data" field.
        Keeping that out when it's unneeded speeds this up quite a bit.'''
        if type(dataset) == dict:
            return

        assert type(dataset) in (h5py.Dataset, h5py.File)

        data_dict = {}

        for key in list(dataset.keys()):
            # Only include "Data" field if it's explicity included.
            # Leave out metadata fields, which are categorized in h5py.Group objects
            if (include_Data or key != "Data") and (type(dataset[key]) != h5py.Group):
                data_dict[key] = dataset[key].value

        return data_dict

    def _populate_coordinate_table(self, flightline_dict):
        '''Using the flightline dictionary, enter all the coordinates from each file
        into the coordinates table.'''
        print("Populating Coordinates Table")
        file_ids = self.file_table.col('File_ID')
        flightline_ids = self.file_table.col('Flightline_ID')
        file_paths = self.file_table.col('relative_filepath')
        trace_counter = 0
        file_counter = 0
        line_file_counter = 0

        last_line_id = None
        last_line_tracenum = 0

        for fname, fid, line_id in zip(file_paths, file_ids, flightline_ids):

            fpath = os.path.join(self.basedir, fname)
            try:
                data = scipy.io.loadmat(fpath)
                datafile = None
            except:
#                print "File in H5 (not Matlab):", fname
                data = h5py.File(fpath, 'r')
                datafile = data
                data = self._convert_h5_dataset_to_dict(datafile)

            # Keep track of files and trace numbers in each flight line.
            # This ASSUMES the files are kept in the same order in which they were read,
            # i.e. that flightlines have their files listed consecutively.
            if line_id != last_line_id:
                if last_line_id != None:
                    print(":", last_line_tracenum, "traces.", line_file_counter, "files.")
                print(line_id, end=' ')
                last_line_tracenum = 0
                last_line_id = line_id
                line_file_counter = 0

            for i in range(len(data['Latitude'].flatten())):
                row = self.coordinates_table.row

                row['Flightline_ID']   = line_id
                row['File_ID']         = fid
                row['Flight_tracenum'] = last_line_tracenum + i # Tells which unique tracenumber in the entire flightline
                row['File_tracenum']   = i # Tells which unique tracenumber in the file
                row['Latitude']        = data['Latitude'].flatten()[i]
                row['Longitude']       = data['Longitude'].flatten()[i]
                row['Elevation']       = data['Elevation'].flatten()[i]
                try:
                    row['Roll']            = data['Roll'].flatten()[i]
                    row['Pitch']           = data['Pitch'].flatten()[i]
                    row['Heading']         = data['Heading'].flatten()[i]
                except KeyError:
                    # Some files don't have Roll, Pitch or Heading parameters.  If not, fill with empty value.
                    row['Roll']    = numpy.nan
                    row['Pitch']   = numpy.nan
                    row['Heading'] = numpy.nan

                row['GPS_time']        = data['GPS_time'].flatten()[i]
                row['Surface']         = data['Surface'].flatten()[i]
                row.append()

            if datafile != None:
                datafile.close()

            last_line_tracenum += (i + 1)
            trace_counter += i
            file_counter += 1
            line_file_counter += 1

        print(":", last_line_tracenum, "traces.", line_file_counter, "files.")

        # Save what we've done so far.
        self.coordinates_table.flush()

        # Get the greenland and antarctica ice flags prepared.
        gr_flags, an_flags = self._compute_ice_flags()
        self.coordinates_table.modify_column(column=gr_flags, colname="OnGreenlandIce")
        self.coordinates_table.modify_column(column=an_flags, colname="OnAntarcticIce")

        # Save again.
        self.coordinates_table.flush()

        print("\t", file_counter, "files read.")
        print("\t", trace_counter, "rows written.")
        return

    def _create_table_indices(self):
        '''Creates indices to optimize future table queries by file names and/or flightlines.'''
        # See "Indexed Searches" at http://www.pytables.org/usersguide/optimization.html
        print("Indexing file_table, Flightline_ID")
        try:
            self.file_table.cols.Flightline_ID.create_index()
        except ValueError:
            self.file_table.cols.Flightline_ID.reindex()

        print("Indexing file_table, File_ID")
        try:
            self.file_table.cols.File_ID.create_index()
        except ValueError:
            self.file_table.cols.File_ID.reindex()

        print("Indexing flightline_table, Flightline_ID")
        try:
            self.flightline_table.cols.Flightline_ID.create_index()
        except ValueError:
            self.flightline_table.cols.Flightline_ID.reindex()

        print("Indexing coordinates_table, Flightline_ID")
        try:
            self.coordinates_table.cols.Flightline_ID.create_index()
        except ValueError:
            self.coordinates_table.cols.Flightline_ID.reindex()

        print("Indexing coordinates_table, File_ID")
        try:
            self.coordinates_table.cols.File_ID.create_index()
        except ValueError:
            self.coordinates_table.cols.File_ID.reindex()

        print("Indexing coordinates_table, Latitude")
        try:
            self.coordinates_table.cols.Latitude.create_index()
        except ValueError:
            self.coordinates_table.cols.Latitude.reindex()

        print("Indexing coordinates_table, Longitude")
        try:
            self.coordinates_table.cols.Longitude.create_index()
        except ValueError:
            self.coordinates_table.cols.Longitude.reindex()

        self.file_table.flush()
        self.flightline_table.flush()
        self.coordinates_table.flush()
        return

    def subset_by_bounding_box(self, LatLatLonLon=(-100,+100,-370,+370)):
        '''Given a set of two latitudes and two longitudes (definding a bounding box),
        return all the coordinate rows that fit withing that box.  Return them as a dictionary
        of key-value pairs.

        Default value simply creates a global subset with all the points.'''
        lat1, lat2, lon1, lon2 = LatLatLonLon

        # Ensure lat2 > lat1, Lon2 > Lon1
        lat1, lat2 = (lat1, lat2) if (lat2 > lat1) else (lat2, lat1)
        lon1, lon2 = (lon1, lon2) if (lon2 > lon1) else (lon2, lon1)

        ret_array = self.coordinates_table.read_where('''(Latitude  >= lat1) & \
                                                         (Latitude  <= lat2) & \
                                                         (Longitude >= lon1) & \
                                                         (Longitude <= lon2)''')

        return ret_array

    def create_line_shapefile_from_subset_array(self, subset_array, shapefile_out):
        '''Given an array (returned by "table.read_where()") of a subset of the data,
        create a shapefile of lines, each shape object defined by the flightline number.'''

        # Create a shapefile Datasource
        shp_driver = ogr.GetDriverByName("ESRI Shapefile")
        data_source = shp_driver.CreateDataSource(shapefile_out)

        # This for creating a geographic lat/lon shapefile
#        # create the spatial reference, WGS84
#        srs = osr.SpatialReference()
#        srs.ImportFromEPSG(4326)

        # This for creating a UTM Zone 20N Shapefile
        utmSR = osr.SpatialReference()
        utmSR.SetWellKnownGeogCS("WGS84")
        utmSR.SetUTM(20,1) # Zone 20, 1 for North
        gcsSR = utmSR.CloneGeogCS()
        utm_converter = osr.CoordinateTransformation(gcsSR, utmSR)

        # Create a layer for the lines
        # create the layer
        layer = data_source.CreateLayer("CC_IceBridge_Lines", utmSR, ogr.wkbLineString)

        # Add the fields we're interested in
        layer.CreateField(ogr.FieldDefn("Flightline", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("Trace_beg", ogr.OFTInteger))
        layer.CreateField(ogr.FieldDefn("Trace_end", ogr.OFTInteger))

        # Organize our subset_array into a list of wktLineStrings
        # Organize into a set of unique flightlines, ordered into a list
        flightlines = list(set(subset_array['Flightline_ID']))
        flightlines.sort()

        for line_id in flightlines:
            # Get the indices of all the points in that flightline.
            indices = numpy.where(subset_array["Flightline_ID"] == line_id)[0]
            # Get the trace numbers.  They *should* be sorted already but let's not assume.
            tracenums = subset_array["Flight_tracenum"][indices]
            trace_argsort = numpy.argsort(tracenums)
            # Sort both the list of indices, and the tracenumbers, by tracenumer
            indices = indices[trace_argsort]
            tracenums = tracenums[trace_argsort]


            # See whether this tracenum is all consecutive.
#            tracenum_diffs = tracenums[1:] - tracenums[:-1] - 1
            if False:
                '''Note: There was code here to break the flightlines up into
                'segments' and write out those segments into separate features.
                
                It is currently disabled. If you want to re-enable it, you'll need
                to uncomment out the code below, and toy with it a while (and
                comment-out the "if False: pass" code-block here) to get it to work.'''
                pass

#            if len(numpy.nonzero(tracenum_diffs)[0]) > 0:
#                print "disjoint segment"
#                # There is more than one segment, must split them up.
#                disjoint_indices = numpy.where(tracenum_diffs != 0)[0] + 1
#                print disjoint_indices
#                print tracenums[disjoint_indices-2], tracenums[disjoint_indices-1], tracenums[disjoint_indices], tracenums[disjoint_indices+1], tracenums[disjoint_indices+2]
#                file_ids = subset_array["File_ID"]
#                print file_ids[indices][disjoint_indices-1], file_ids[indices][disjoint_indices], file_ids[indices][disjoint_indices+1]
#                foobar
#                segments = [(0,disjoint_indices[0])]
#                for i,index in enumerate(disjoint_indices):
#                    if i==len(disjoint_indices)-1:
#                        segments.append((index,len(tracenums)))
#                    else:
#                        segments.append((index,disjoint_indices[i+1]))

            else:
                segments = [(0,len(tracenums))]


            # Create a line feature for each line segment (even if only one)
            for segment in segments:
                lats = subset_array["Latitude"][indices][segment[0]:segment[1]]
                lons = subset_array["Longitude"][indices][segment[0]:segment[1]]

                # For UTM Points
                lonlat_points = [(lons[i], lats[i]) for i in range(len(lats))]
                utm_points = utm_converter.TransformPoints(lonlat_points)

                # Create geometry well-known-text
                wkt = "LINESTRING( "
                for i in range(len(lats)):
#                    wkt = wkt + str(lons[i]) + " " + str(lats[i]) + ", "
                    wkt = wkt + str(utm_points[i][0]) + " " + str(utm_points[i][1]) + ", "
                wkt = wkt[:-2] + ")"
                line = ogr.CreateGeometryFromWkt(wkt)

                # Create the feature
                feature = ogr.Feature(layer.GetLayerDefn())
                # Set Feature attributes
                feature.SetField("Flightline", int(line_id))
                feature.SetField("Trace_beg", float(tracenums[segment[0]]))
                feature.SetField("Trace_end", float(tracenums[segment[1]-1]))
                feature.SetGeometry(line)

                # Insert the feature into the layer
                layer.CreateFeature(feature)
                feature.Destroy()

        data_source.Destroy()


    def subset_by_shapefile(self, shapefile_clip, shapefile_out = None, return_fields=['File_ID','File_Tracenum','Flightline_ID','Flight_Tracenum']):
        '''With a single-value polygon shapefile, subset all the points that fit
        within that shapefile.  Return a shapefile (in the same coordinate
        projection as the shapefile provided), with "return_fields" designating
        the columns of the coordinates_table that will be listed as attributes in the shapefile.'''

        # TODO: input logic for opening the shapefile and outputting coordinates.
        # In order to do this, we must be able to read in different coordinate systems and
        # reproject the native lat/lon coordinates to perform clipping
        # For now, handle Lon/Lat and UTM E/N coordinate systems.

        return

# First global function, build the database from scratch.
def BUILD_DATABASE():
    db = IceBridgeRadarDB(dbfile=DB_FILE, datadir=BASEDIR)
    db.build_database()
    db._create_table_indices()

#    db.open_database(status="r+")
    gr_ice, an_ice = db._compute_ice_flags()
        # Get the greenland and antarctica ice flags prepared.
    db.coordinates_table.modify_column(column=gr_ice, colname="OnGreenlandIce")
    db.coordinates_table.modify_column(column=an_ice, colname="OnAntarcticIce")
    db.coordinates_table.flush()

def CREATE_CAMP_CENTURY_SHAPEFILE():
    db = IceBridgeRadarDB(dbfile=DB_FILE, datadir=BASEDIR)
    db.open_database()
    print("Subsetting data...")
    subset = db.subset_by_bounding_box((77.167, 77.195,-61.22525,-60.985))
    db.create_line_shapefile_from_subset_array(subset, CC_LINE_SHAPEFILE)
    return db


if __name__ == "__main__":

    BUILD_DATABASE()
    
#    db = IceBridgeRadarDB()
#    db.open_database(status='r')

#    lats = db.coordinates_table.cols.Latitude[:]
#    lons = db.coordinates_table.cols.Longitude[:]
#    grice = db.coordinates_table.cols.OnGreenlandIce[:]
#    anice = db.coordinates_table.cols.OnAntarcticIce[:]
#
#    print lats.shape, lats.dtype
#    print lons.shape, lons.dtype
#    print grice.shape, grice.dtype
#
#    plt.scatter(lons[::100][~(anice[::100])], lats[::100][~(anice[::100])], color="red", marker=".")
#    plt.scatter(lons[::100][anice[::100]], lats[::100][anice[::100]], color="blue", marker=".")
#    plt.xlim(-180, 180)
#    plt.ylim(-90, -55)
#    plt.show()

#    db = CREATE_CAMP_CENTURY_SHAPEFILE()