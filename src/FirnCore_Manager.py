# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 15:27:22 2016

@author: mmacferrin

FirnCore_Manager.py
    - Creates an object for each core. -- use / adapt previous code
    - Creates an object for the in situ GPR segments -- use previous code
    - Creates variance images for in situ GPR, finds cutoff that maximises core agreement (with thick ice lenses)
    - Creates an object for IceBridge Accumulation GPR profiles.
        - Use that object to find detection algorithm for IceBridge GPR, both against cores & in situ GPR.
    - Run IceBridge picker algorithm over all of Greenland Ice Sheet.  Hoorah.  Create map of Greenland ice lens GPRs lines.
"""
# Python 2.7.9 included libraries
import os
import datetime

# External open-source libraries
# xlrd version 1.0.0, for reading and parsing Microsoft Excel (.xls or .xlsx) files, rather than exporthing them all to CSV.
import xlrd
import matplotlib.pyplot as plt
import numpy

from GPR_FileData import CORE_XLSX_FILE, \
                         CORE_FOLDER


class FirnCore_Manager():
    '''FirnCoreManager: Stores information and manages processes on all firn cores.
    Uses: FirnCore_Profile()
    Used by: InSituGRP_Manager, IceBridgeGPR_Manager
    '''
    def __init__(self, core_fname=CORE_XLSX_FILE, verbose = True):
        # Open Excel file
        self.verbose = verbose

        if self.verbose:
            print "Opening", os.path.split(core_fname)[-1]

        self.workbook = xlrd.open_workbook(core_fname)
        self.cover_sheet = self.workbook.sheet_by_index(0)
        # All sheet names exclusing the cover sheet, exclude "core_template_2"
        self.sheet_names = [self.workbook.sheet_by_index(i).name for i in range(1,self.workbook.nsheets) if self.workbook.sheet_by_index(i).name != "core_template_2"]


    def _test_sheet_headers(self, xlsx_workbook):
        '''A quick test to make sure the worksheet headers (apart from the cover page) are all the same.
        If so, return.
        If not, throw a ValueError

        ## NOT ALL THE SHEETS AGREE ON HEADER VALUES.  2015 ONWARD INCLUDES A 'CIRCUMFERENCE' TAG THAT
        THE PREVIOUS SHEETS DON'T.  ALL THE OTHER HEADER COLUMNS ARE THE SAME THOUGH, WITH THE SAME MEANING.
        WE WILL SIMPLY INDEX THE COLUMNS BY NAME USING A DICTIONARY.

        This function will be omitted from the main processing.
        '''
        # Skip sheet 0 (cover sheet)
        for i in range(1,xlsx_workbook.nsheets):
            if i==1:
                # Returns a python list of cells
                base_header = xlsx_workbook.sheet_by_index(i).row(0)
            else:
                header = xlsx_workbook.sheet_by_index(i).row(0)
                for j,cell in enumerate(header):
                    if base_header[j].value != header[j].value:
                        raise ValueError("Sheet {0}, column {1} ('{2}'): Sheet headers do not agree: ('{3}')".format(\
                        xlsx_workbook.sheet_by_index(i).name, j, header[j].value, base_header[j].value))
        return

    def return_core_profile(self, name):
        return FirnCore_Profile(self.workbook.sheet_by_name(name), title_sheet=self.cover_sheet, verbose = self.verbose)

    def return_core_profiles_from_year(self, year):
        assert type(year) in (int, float, str)
        # Type conversion just assures that it'll end up an integer string, no matter which input type is given
        year_str = str(int(numpy.round(float(year))))

        # Extract years with "year_str" string, create Firn Core profiles for eah.
        profile_names = [self.return_core_profile(name) for name in self.sheet_names if name.find(year_str) > -1]
        return profile_names

    def plot_core_profiles(self, depth_cm):
        '''Plots all the core profiles with density, stratigraphy.'''
        for name in self.sheet_names:
            sheet = self.workbook.sheet_by_name(name)
            profile = FirnCore_Profile(sheet, self.cover_sheet)
            profile.plot_ice_lenses(depth_cm=depth_cm)
        return


    # Outputs statistics of core profiles for Table X in paper.
    def output_core_statistics(self, filename="Core_Statistics.csv", depth_cm=None):
        '''Outputs core statistics to a CSV file, to be used in Table X in Karen Alleys Scatterometry paper.

        Also used for getting the core locations and  basic stats for GPR comparision'''
        outstring = ''
        for i,name in enumerate(self.sheet_names):
            sheet = self.workbook.sheet_by_name(name)
            profile = FirnCore_Profile(sheet, self.cover_sheet)
            outstring = outstring + profile.core_stats_CSV(depth_cm=depth_cm,
                                    return_with_header=(True if i==0 else False))

        if filename != None:
            f = open(os.path.join(CORE_FOLDER, filename), 'w')
            f.write(outstring)
            f.close()

            if self.verbose:
                print filename, "written."

        return outstring


class FirnCore_Profile():
    '''FirnCore_Profile: Stores information and processing for a single firn core.
    Uses: None
    Used by: FirnCore_Manager
    '''
    def __init__(self, xlsx_sheet, title_sheet=None, verbose=True):
        self.sheet = xlsx_sheet
        self.title_sheet = title_sheet
        self.title_sheet_header = self._cover_sheet_header()
        self.header = self._create_header()
        self.verbose = verbose
        if verbose: print "Core", self.name()

        # Boolean 1-d array of True/False for each cm of the core, whether it's ice
        self.boolean_ice_array = None
        # 2xN array of top/bottom depths of each "continuous" ice lens
        self.ice_lens_array = None

    def name(self):
        return self.sheet.name

    def name_year_first(self):
        '''Return the core name, with the name rearranged so that the year comes first.'''
        name = self.sheet.name
        year = name[-4:]
        try:
            int(year)
        except ValueError:
            print "Warning: something wrong with putting year first,", name, "-->", year

        return year + "_" + name[:-5]

    def lat_lon_elev(self):
        '''Return the latitude, longitude, elevation of the core.'''
        row_idx = self.title_sheet.col_values(self.title_sheet_header["core"]).index(self.name())
        return (self.title_sheet.cell_value(row_idx, self.title_sheet_header["N1"]),
                self.title_sheet.cell_value(row_idx, self.title_sheet_header["W1"]),
                self.title_sheet.cell_value(row_idx, self.title_sheet_header["Z1"])
                )

    def date_drilled(self, datemode):
        row_idx = self.title_sheet.col_values(self.title_sheet_header["core"]).index(self.name())
        col_idx = self.title_sheet_header["date cored"]
        date_tuple = xlrd.xldate_as_tuple(self.title_sheet.cell_value(row_idx, col_idx), datemode)
        return datetime.datetime(*date_tuple)

    def depths(self):
        '''Return the depths in cm of the core.'''
        depths = self.sheet.col_values(self.header["depth"])
        assert depths[0].lower() == u'depth'
        depths = depths[1:]
        if type(depths[-1]) == str and depths[-1].lower() == "end":
            depths = depths[:-1]

        return numpy.array(depths)

    def densities(self, correct_zeros=False):
        '''Return the density profile of the core.
        If "correct_zeros" is True, then fill in any zero values.  At the surface, use 340 kg/m3.
        At deeper spans of the core, interpolate linearly between the non-zero values above and below the gap.'''
        if "density" not in self.header.keys():
            return None
        density = self.sheet.col_values(self.header["density"])
        assert density[0].lower() == "density"
        density = numpy.array(density[1:])
        # Assure we're getting an array of all float values (occasionally returns strings)
        if density.dtype != numpy.float64:
            density_float = numpy.empty(density.shape, dtype=numpy.float64)
            # If the cells contain some non-float values, such a blank lines or text,
            # insert "empty" 0.0 values
            for i,d in enumerate(density):
                try:
                    density_float[i] = float(d)
                except ValueError:
                    density_float[i] = 0.0
            density = density_float

        if not correct_zeros:
            return density

        while True:
            # Find the first index (after the placeholder) where the density is zero
            zeros_idx = numpy.where(density == 0.0)[0]
            # If found none, exit. We're done.
            if len(zeros_idx) == 0:
                break
            # If we find one or more, retrieve the first index.
            else:
                zeros_idx = zeros_idx[0]

            # From the start, find the first index where the density is nonzero
            nonzeros_idx = numpy.nonzero(density[zeros_idx:])[0]
            if len(nonzeros_idx) == 0:
                # If the zero-span is at the end, simply extend the density of the last
                # known value out to the end.
                density[zeros_idx:] == density[zeros_idx-1]
                assert len(numpy.where(density==0.0)[0]) == 0
                break
            else:
                nonzeros_idx = nonzeros_idx[0] + zeros_idx

            # If the span is at the beginning, set snow to a bulk value of 340.00 kg/m3
            if zeros_idx == 0:
                density[zeros_idx:nonzeros_idx] = 340.0
            else:
                # It's a span in the middle.  In this case, linearly interpolate from the
                # last known value to the next known value
                # NOTE: This could infinitely loop if one of the interpolated values was zero,
                # but since all our densities must be positive, this shouldn't happen.
                prev_density = density[zeros_idx-1]
                next_density = density[nonzeros_idx]
                numvalues = nonzeros_idx - zeros_idx
                # Linear interpolation between prev_density and next_density
                density[zeros_idx:nonzeros_idx] = prev_density + \
                    (next_density - prev_density)*(numpy.arange(1,numvalues+1,dtype=numpy.float))/float(numvalues)

        return density

    def _create_header(self):
        '''Create a header-index-lookup for the sheet, using the top row of the worksheet.
        Header will be a dictionary of [name]:[index] pairs'''
        header = {}
        header_values = [cell.value for cell in self.sheet.row(0)]
        for i,value in enumerate(header_values):
            if value.strip() != '':
                header[value] = i

        return header

    def _cover_sheet_header(self):
        title_header_values = [cell.value for cell in self.title_sheet.row(0) if cell.value.strip() != '']
        title_header = dict(zip(title_header_values, range(len(title_header_values))))
        return title_header

    def create_boolean_ice_array(self, skip_calculations = False):
        # If "skip_calculations" is True and we've already calculated it, just return what's already there.
        if skip_calculations and type(self.boolean_ice_array) != type(None):
            return self.boolean_ice_array

        depths = self.sheet.col_values(self.header["depth"])
        type_values = self.sheet.col_values(self.header["type"])
        assert depths[0] == u'depth'
        assert type_values[0] == u'type'
        # After verifying the correct column, cut off the header cell
        depths = depths[1:]
        type_values = type_values[1:]

        depths = numpy.array(depths)
        bool_array = numpy.array([(cell.lower().find('ice') > -1) for cell in type_values])
        assert depths.shape == bool_array.shape

        self.boolean_ice_array = bool_array

        return bool_array

    def create_ice_lens_array(self):
        '''Returns a 2xN array of where ice lenses are contained in the core segment.
        Each row is a lens, with the top index and the bottom index of the array,
        inclusive and non-inclusive, respectively.  This means that a 2-cm thick array
        starting at layer 68 (zero-indexed) will have [68,70] in that row of the array.
        This allows for easy calculation of thicknesses (bottom - top).'''
        if self.boolean_ice_array is None:
            self.create_boolean_ice_array()
        assert self.boolean_ice_array is not None

        # If this is already done, don't bother re-doing, just return it as is.
        if self.ice_lens_array is not None:
            return self.ice_lens_array

        ice_lens_list = []

        ice_array = self.boolean_ice_array
        temp_idx = 0
        # Iterate over the core, finding the spans where the boolean_ice_array is true.
        while temp_idx < len(ice_array):
            # Find the first index (after the placeholder) where the boolean array is "True"
            start_idx = numpy.nonzero(ice_array[temp_idx:])[0]
            # If found none, exit. We're done.
            if len(start_idx) == 0:
                break
            # If we find one, add it to the placeholder to get the index of the whole array.
            else:
                start_idx = start_idx[0] + temp_idx

            # From the start, find the first index where the boolean array is "False"
            end_idx = numpy.nonzero(~ice_array[start_idx:])[0]
            # If not found, make it the end of the array (+1, for easy thickness computations)
            if len(end_idx) == 0:
                end_idx = len(ice_array)
            # If found, set it there.
            else:
                end_idx = end_idx[0] + start_idx

            # The lens is the start_idx, end_idx, marking the top & bottom of the ice layer.
            ice_lens_list.append((start_idx, end_idx))

            # Reset the temporary placeholder to search further down the core.
            temp_idx = end_idx

        ice_lens_array = numpy.array(ice_lens_list)
        self.ice_lens_array = ice_lens_array

        return ice_lens_array

    def ice_lens_thicknesses(self):
        '''Return the thicknesses of the lenses calculated in "create_ice_lens_array()".'''
        if self.ice_lens_array is None:
            self.create_ice_lens_array()
        return self.ice_lens_array[:,1] - self.ice_lens_array[:,0]

    def ice_content_percent(self, depth_cutoff_cm=None, use_partials=True):
        '''Return the percentage of ice in the firn core, as logged.
        If depth_cutoff is used, return only % from the top to that depth in the core.
        The "use_partials" parameter determines whether we use the "percentages" of ice
        in each layer as outlined in the workbook.  If False, each layer with any ice counts
        as all-ice.  If True, a "50" value in the "type_perc" column will mean only 50% of that
        layer is countes as ice-content.  Blank cells in "type_perc" are treated as 100%.'''
        boolean_ice_array = self.create_boolean_ice_array()
        if boolean_ice_array == None:
            return 0

        if use_partials:
            # First filter out only the layers that contain ice, using boolean_ice_array
            ice_layer_perc_values = [None] * int(numpy.sum(boolean_ice_array))
            cell_values = self.sheet.col_values(self.header["type_perc"])[1:]
            perc_values_i = 0
            for i in range(len(cell_values)):
                if self.boolean_ice_array[i] == True:
                    ice_layer_perc_values[perc_values_i] = cell_values[i]
                    perc_values_i += 1
                else:
                    continue
            # Then multiply by percentages and sum it all up.
            total_ice_thickness = numpy.sum([value/100.0 if type(value) == float else 1.0 for value in ice_layer_perc_values])
            return total_ice_thickness / float(len(boolean_ice_array))
        else:
            return numpy.sum(boolean_ice_array) / float(len(boolean_ice_array))

        # Should never get here.
        assert False

    def thickest_layer_top_N_cm(self, max_depth_cm=None):
        '''Find the thickness of the thickest ice layer in the top N centimeters of core.'''
        ice_lens_array = self.create_ice_lens_array()
        if ice_lens_array is None or len(ice_lens_array) == 0:
            return 0

        ice_lens_tops = ice_lens_array[:,0]
        thicknesses = self.ice_lens_thicknesses()

        if max_depth_cm != None:
            # Concatenate ice_lens_array to max thickness (omit layers deeper than that)
            ice_lens_tops = ice_lens_tops[ice_lens_tops <= max_depth_cm]
            thicknesses = thicknesses[:len(ice_lens_tops)]

        return numpy.max(thicknesses)

    def depth_of_first_layer_of_N_thickness(self, thickness, max_depth_cm=None):
        '''Look for layers of a given thickness in the core in the top "max_depth_cm",
        and return the depth of the top of the first layer that meets that thickness, in cm.
        Return -1 if no such layer exists of that thickness.'''
        # We just need the top of the ice lens depths, not the bottom
        ice_lens_array = self.create_ice_lens_array()
        if ice_lens_array is None or len(ice_lens_array) == 0:
            return -1

        ice_lens_tops = ice_lens_array[:,0]
        thicknesses = self.ice_lens_thicknesses()
        depths = self.depths()

        if max_depth_cm != None:
            # Concatenate ice_lens_array to max thickness (omit layers deeper than that)
            ice_lens_tops = ice_lens_tops[ice_lens_tops <= max_depth_cm]
            thicknesses = thicknesses[:len(ice_lens_tops)]

        lens_indices = numpy.where(thicknesses >= thickness)[0]
        if len(lens_indices) == 0:
            return -1
        else:
            return depths[ice_lens_tops[lens_indices[0]]]


    def plot_ice_lenses(self, depth_cm = None, ax=None):
        '''Create a plot with the top "depth_cm" of firn's ice lenses included in it.'''
        figure_name = os.path.join(CORE_FOLDER, self.name_year_first() + ".png")

        # Gather the depths & densities, omit the blank or zero values
        depths = self.depths()
        densities = self.densities()
        if densities != None:
            density_mask = [False if (D=="" or D==0.0) else True for D in densities]
            depths_good = numpy.array([depths[i] for i in range(len(density_mask)) if density_mask[i]], dtype=numpy.float)
            densities_good = numpy.array([densities[i] for i in range(len(density_mask)) if density_mask[i]], dtype=numpy.float)

        if ax is None:
            fig, ax = plt.subplots(figsize=(1.75,5.0), dpi=150)

        if self.verbose:
            print os.path.split(figure_name)[-1]

        ice_lens_array = self.create_ice_lens_array()
        # Skip plotting ice lenses for any core in which they're not available or don't exist.
        if len(ice_lens_array.shape) == 2 and ice_lens_array.shape[1] == 2:
            for i in range(ice_lens_array.shape[0]):
                ax.axhspan(ymin=ice_lens_array[i,0]/100.0, ymax=ice_lens_array[i,1]/100.0, color="lightblue")

        if densities != None:
            ax.step(densities_good, depths_good/100.0, color="black")

        ax.set_xlim(0,1000)
        ax.set_ylim((2000 if depth_cm is None else depth_cm)/100.0, -0.15)
        ax.axhspan(ymin=-1,ymax=0, color="lightgrey")
        ax.axhline(y=0, color="black")

        if ax is None:
            ax.set_xlabel("Density (kg m$^{-3}$)")
            ax.set_ylabel("Depth (m)")

        plt.setp(ax.get_xticklabels(), rotation="vertical", fontsize=10)

        if ax is None:
            ax.title(self.name())

            plt.tight_layout()

            plt.savefig(figure_name)
            plt.cla()
            plt.close()

        return

    def core_stats_CSV(self, depth_cm = None, return_with_header = False):
        '''Output the statistics of the core in the following format:
        name,nearest code location,N1,W1,Z1,depth,corer,ice fraction,thickest_ice,depth ice >5cm,depth ice >10cm,depth ice >20cm,depth ice >50cm'''
        output_header = "name,nearest code location,N1,W1,Z1,depth,corer,ice fraction,thickest_ice,depth ice >5cm,depth ice >10cm,depth ice >20cm,depth ice >50cm\n"
        title_header = self.title_sheet_header

        name = self.name()
        row_idx = self.title_sheet.col_values(title_header["core"]).index(name)
        row = self.title_sheet.row(row_idx)

        nearest_code_location = row[title_header["nearest code location"]].value
        N1 = row[title_header["N1"]].value
        if N1 == '': N1 = 0
        W1 = row[title_header["W1"]].value
        if W1 == '': W1 = 0
        Z1 = row[title_header["Z1"]].value
        if Z1 == '': Z1 = 0
        depth = row[title_header["depth"]].value
        corer = row[title_header["corer"]].value
        pct_ice = self.ice_content_percent(depth_cutoff_cm=depth_cm)
        thickest_ice = self.thickest_layer_top_N_cm(max_depth_cm=depth_cm)
        depth_5cm = self.depth_of_first_layer_of_N_thickness(5,max_depth_cm=depth_cm)
        depth_10cm = self.depth_of_first_layer_of_N_thickness(10,max_depth_cm=depth_cm)
        depth_20cm = self.depth_of_first_layer_of_N_thickness(20,max_depth_cm=depth_cm)
        depth_50cm = self.depth_of_first_layer_of_N_thickness(50,max_depth_cm=depth_cm)

        try:
            outstr = "{0:s},{1:s},{2:f},{3:f},{4:f},{5:f},{6:s},{7:f},{8:f},{9:f},{10:f},{11:f},{12:f}\n".format( \
            name,nearest_code_location,N1,W1,Z1,depth,corer,float(pct_ice),float(thickest_ice),float(depth_5cm),float(depth_10cm),float(depth_20cm),float(depth_50cm))
        except ValueError:
            print name,nearest_code_location,N1,W1,Z1,depth,corer,pct_ice,thickest_ice,depth_5cm,depth_10cm,depth_20cm,depth_50cm
            print "name", type(name)
            print "nearest_code_location", type(nearest_code_location)
            print "N1", type(N1)
            print "W1", type(W1)
            print "Z1", type(Z1)
            print "depth", type(depth)
            print "corer", type(corer)
            print "pct_ice", type(pct_ice)
            print "thickest_ice", type(thickest_ice)
            print "depth_5cm", type(depth_5cm)
            print "depth_10cm", type(depth_10cm)
            print "depth_20cm", type(depth_20cm)
            print "depth_50cm", type(depth_50cm)
            # This is a critical mistake, if we get here go ahead and print diagnostic info and crash.
            assert False

        if return_with_header:
            return output_header + outstr
        else:
            return outstr

