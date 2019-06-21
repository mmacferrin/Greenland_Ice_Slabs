# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 12:31:34 2013

@author: macferrin
"""
import numpy
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.image as mpimage
#from scipy import ndimage
#from scipy import misc
#from FileData import GPR_ASCII_GRD_FILES,
#                     GPR_IMG_FILES,            \
#                     IMG_DIR,                  \
#                     RESAMPLED_DIR,            \
#                     KANU_TRANSECT_COLLECTION, \
#                     NON_TRANSECT_COR_FILES,   \
#                     NON_TRANSECT_GPR_FILES
#                     GPR_COR_FILES_PROCESSED,  \
#                     GPR_MERGED,               \
#                     COR_MERGED,               \
#                     GPR_MERGED_RESAMPLED,     \
#                     COR_MERGED_RESAMPLED,     \
#                     GPR_IMG_FILES_RESAMPLED,  \
#                     GPR_COR_FILES_RESAMPLED,  \
#                     CORE_FILENAMES,           \

# SOME FILENAMES WERE NOT UPDATED WITH FileData and are currently breaking the script.
# Just define empty names here:

GPR_ASCII_GRD_FILES = None
GPR_IMG_FILES = None
GPR_COR_FILES_RESAMPLED = None
CORE_FILENAMES = None
FirnCoreProfiler = None
GPR_IMG_FILES_RESAMPLED = None
GPR_COR_FILES_PROCESSED = None
COR_MERGED = None
GPR_MERGED = None
COR_MERGED_RESAMPLED = None
GPR_MERGED_RESAMPLED = None

from FileData import EXPORT_DIR,               \
                     MAIN_TRANSECT_COLLECTION, \
                     KANU_GPR_MERGED_RESAMPLED, \
                     KANU_COR_MERGED_RESAMPLED, \
                     GPR_DETRENDED_PICKLEFILES, \
                     TRANSECT_GPR_MERGED_RESAMPLED, \
                     TRANSECT_COR_MERGED_RESAMPLED


import mahotas # Basic image saving/processing. -- No longer supports floating point save??
from PIL import Image
import tifffile as tiff
#import PythonMagick as PM
from scipy import stats, optimize
#import FirnCoreProfiler
from GPR_Coordinates import GPR_Coords
import pickle
import warnings
warnings.simplefilter('ignore', Image.DecompressionBombWarning)

def speed_of_light_in_ice(density_kg_m3):
    '''Given an ice density (between 0 and 917 kg/m3), compute the approximate
    speed of light through that medium given an index of refraction of 1.3078 and
    speed in pure ice of .2294 m/ns.  (It's just a scaling between this.).
    Return value in m/s.'''
    ice_density = 917.00
    assert 0 <= density_kg_m3 <= ice_density
    c = 299792458
    ni = 1.31
    na = 1.0
    vi = c/(na+(ni-na)*(density_kg_m3/ice_density))
    return vi


class GPR_Collection:
    '''GPR_Collection does batch processing on individual GPR_Datafile objects.'''
    def __init__(self, grd_files=GPR_ASCII_GRD_FILES, img_files=GPR_IMG_FILES):
        self.grd_files = grd_files
        self.img_files = img_files

    def import_all_images(self):
        for g in self.grd_files:
            gpr = GPR_DataFile(g)
            del gpr

    def subtract_mean_all_images(self):
        for f in self.img_files:
            gpr = GPR_DataFile(f)
            gpr.subtract_horizontal_mean_transform()

    def export_all_images(self):
        for f in self.img_files:
            gpr = GPR_DataFile(f)
            gpr.export_all_images()

    def display_all_histograms(self, fignum=2):
        assert len(self.img_files) == 20

        gs = gridspec.GridSpec(4,5)

        for i,f in enumerate(self.img_files):
            ax = plt.subplot(gs[i])
            gpr = GPR_DataFile(f)
            gpr.display_histogram(fignum, show=False)
            ax.set_yticklabels([])
            del gpr

        plt.show()

    def display_all_std_curves(self, fignum=2):
        assert len(self.img_files) == 20

        gs = gridspec.GridSpec(4,5)

        for i,f in enumerate(self.img_files):
            ax = plt.subplot(gs[i])
            gpr = GPR_DataFile(f)
            gpr.display_mean_std_curves(fignum, show=False)
            ax.set_yticklabels([])
            del gpr

        plt.show()
        plt.close()

class GPR_IceLensIdentifier:
    '''A class for combining core data with gpr data to identify ice lenses.'''
    def __init__(self, GPR_FILES = GPR_IMG_FILES_RESAMPLED, \
                       GPR_CORFILES = GPR_COR_FILES_RESAMPLED, \
                       COREFILES = CORE_FILENAMES):
        self.cores = [FirnCoreProfiler.CoreProfile(f) for f in CORE_FILENAMES]
        self.gpr_files = GPR_FILES
        self.gpr_corfiles = GPR_CORFILES
        self.closest_traces = {}

    def find_closest_traces(self, numtraces = 2):
        '''For each core, find the GPR traces that are closest to it (omit Core 9).
        Save the output in a dictionary {} with the CORE_FILENAME as keys, and a
        "numtraces"-length list of tuples as items.  The tuples are defined in the
        GPR_Coordinates.GPR_Coords.closest_traces() function.'''
        cor_objects = [GPR_Coords(f) for f in self.gpr_corfiles]
        for core in self.cores:
            # Skip core 9
            if core.filename.find("core_9_2013") != -1:
                continue

            trace_tuples_all = []
            latlon = (core.latitude, core.longitude)
            # Collect the closest traces from each gpr file
            for gprcor_object in cor_objects:
                trace_tuples_all.extend(gprcor_object.closest_traces(latlon, numtraces = numtraces))

            # Find the N (numtraces) closest traces from all the files.
            distances = numpy.array([tup[-1] for tup in trace_tuples_all], dtype=numpy.double)
            sorted_indices = numpy.argsort(distances)

            self.closest_traces[core.filename] = [trace_tuples_all[i] for i in sorted_indices[:numtraces]]

#        for key in self.closest_traces.keys():
#            print key
#            for trace in self.closest_traces[key]:
#                print trace
#            print

        return self.closest_traces

    def best_fit_minimization(self):
        '''This is the major function which will evaluate the cutoff variables to get the
        best "fit" of variables for identifying ice lenses.  Namely, I need to find
        - (a) Attentuation coefficient (alpha) for solid ice lenses
        - (b)     ""           ""         ""   for porous firn
        - (c) Roughness cutoff (variance threshold) for identifying "ice"

        Evaulation will be based upon two minimization factors, in the order of importance:
        1. False positives - layers identified as "ice" that are actually firn
                - It is important to get rid of these.
        2. False negatives - layers identified as "firn" that are actually ice in the cores.
                - There will be some of these, especially with thin layers.

        Should seek to minimize false positives (near-zero) while also minimizing false-negatives.
        Let's look at the false-positive relationship first.  How to minimize those?

        APPROACH:  Get alpha_firn first from firn lenses, and a window size! (to eliminate false positives)
        Might be able to get roughness cutoff then too.  Then, apply to ice/firn mixes to get alpha_ice.
        '''
        pass

class GPR_Merger:
    '''A class for merging 2 different GPR files into one.'''
    def __init__(self, TRANSECT_COLLECTION=MAIN_TRANSECT_COLLECTION, resampled=True):
        self.master_traces = numpy.array([],dtype=numpy.int16)
        self.master_traces.shape = (2024,0)
        self.master_corlines = []
        self.transect_indices =  [t[0] for t in TRANSECT_COLLECTION]
        self.transect_switches = [t[2] for t in TRANSECT_COLLECTION]
        self.transect_merges   = [t[3] for t in TRANSECT_COLLECTION]

        self.default_GPR_file = GPR_MERGED_RESAMPLED if resampled else GPR_MERGED
        self.default_COR_file = COR_MERGED_RESAMPLED if resampled else COR_MERGED
        self.GPR_file = None
        self.COR_file = None

        self.is_resampled = resampled
        return

    def merge_all(self, GPR_outfile=None, COR_outfile=None):
        '''COR will be given a ".mod.cor" extension, to indicate that the
                COR_outfile may contain comments relating to the original files.
                All comments are noted with a ' tag at the beginning of the line.
           SWITCH2:  True if second coordinate file needs to be switched.  False if not.
           MERGETYPE: 0:  Don't crop either files.  Just append them fully.
                      1:  Keep all the first, crop the second, according to whichever
                          point's coordinates most closely follow the last point of
                          the first GPR.  Since this transect goes all westward,
                          the points will be sorted by ascending order of longitude.
                      2:  Keep all the second, crop the first.  Same rules as above.'''
        if GPR_outfile is None:
            GPR_outfile = self.default_GPR_file
        if COR_outfile is None:
            COR_outfile = self.default_COR_file

        self.GPR_file = GPR_outfile
        self.COR_file = COR_outfile

        # Loop through each transect file
        GPR_FILES = GPR_IMG_FILES_RESAMPLED if self.is_resampled else GPR_IMG_FILES
        for i in range(len(self.transect_indices)):
            gpr = GPR_DataFile(GPR_FILES[self.transect_indices[i]])
            merge_type = self.transect_merges[i]
            switched = self.transect_switches[i]

            print
            print "============================================================="
            print os.path.split(gpr.filename)[1]

            if gpr.traces is None:
                gpr.open_image()
            traces = gpr.traces
            corlines = gpr.cor_lines()

            # Keep the comments separate from the data lines, for now.
            corcomments = [self._corfile_comment_line(gpr.filename)]

            # trace_start_i -- the beginning trace index
            # trace_end_i -- the end trace index
            # cor_start_i -- the beginning corline INDEX
            # cor_end_i -- the end corline INDEX

            # Handle next lines, and each switching case (elegant way to do that?)
            trace_increment = -1 if switched else 1

            if merge_type == 0:
                if switched:
                    trace_start_i = traces.shape[1]-1
                    cor_start_i = len(corlines)-1
                else:
                    trace_start_i = 0
                    cor_start_i = 0

            elif merge_type == 1:
                assert len(self.master_corlines) > 0
                # ensure last line isn't a comment
                assert self.master_corlines[-1][0] != "'"

                if switched:
                    cor_start_i = len(corlines)-1
                else:
                    cor_start_i = 0

                # Fetch the longitude of the last master trace:
                items = self.master_corlines[-1].split()
                assert(items[6] in ("W","E"))
                last_longitude = float(items[5]) * (-1 if items[6]=="W" else 1)

                for i,corline in zip(range(cor_start_i,-1 if switched else len(corlines),trace_increment),
                                     corlines[cor_start_i:None:trace_increment]):
                    items = corline.split()
                    assert(items[6] in ("W","E"))
                    this_longitude = float(items[5]) * (-1 if items[6]=="W" else 1)
                    if this_longitude > last_longitude:
                        cor_start_i = i
                        # Get the beginning trace number from the corfile line.
                        # The trace numbers are 1-indexed, so subtract one to get the array index
                        trace_start_i = int(items[0])-1
                        break

            elif merge_type == 2:
                if switched:
                    raise ValueError("Merge 'case 2, switched' defined, but unhandled yet.")
                else:
                    cor_start_i = 0

                # Fetch the longitude of the last master trace:
                items = self.master_corlines[-1].split()
                assert(items[6] in ("W","E"))
                last_longitude = float(items[5]) * (-1 if items[6]=="W" else 1)

                for i,corline in zip(range(cor_start_i,-1 if switched else len(corlines),trace_increment),
                                     corlines[cor_start_i:None:trace_increment]):
                    items = corline.split()
                    assert(items[6] in ("W","E"))
                    this_longitude = float(items[5]) * (-1 if items[6]=="W" else 1)
                    if this_longitude > last_longitude:
                        cor_start_i = i
                        # Get the beginning trace number from the corfile line.
                        # The trace numbers are 1-indexed, so subtract one to get the array index
                        trace_start_i = int(items[0])-1
                        break


            else:
                raise ValueError("Unknown merge case: {0}".format(merge_type))

            # Modify corlines to reflect new tracefile
            corlines = corlines[cor_start_i:None:trace_increment]
            traces = traces[:,trace_start_i:None:trace_increment]

            corlines = self._modify_corlines(corlines,
                                             self.master_traces.shape[1],
                                             traces.shape[1],
                                             int(corlines[0].split()[0]),
                                             switched)

            # Add the number of traces to corlines comments
            corcomments.append("' {0} traces\n".format(traces.shape[1]))

            # Append corcomments, corlines to master_corlines
            self.master_corlines.extend(corcomments)
            self.master_corlines.extend(corlines)

            # Append traces to master_traces
            self.master_traces = numpy.append(self.master_traces, traces, axis=1)

            print len(corlines), "new corlines."
            print traces.shape[1], "new traces."

            for comment in corcomments:
                print comment,
            print corlines[0],
            print corlines[-1]
            print len(self.master_corlines), "total corlines."
            print self.master_traces.shape[1], "total traces."

        print "Saving", os.path.split(GPR_outfile)[1]
        mahotas.imsave(GPR_outfile, self.master_traces)
        print "Saving", os.path.split(COR_outfile)[1]
        f = open(COR_outfile, 'w')
        f.write("".join(self.master_corlines))
        f.close()


    def _modify_corlines(self, corlines, master_traces_len, traces_len, trace_start_num, switched):
        new_corlines = []
        for line in corlines:
            items = line.split()
            ti = int(items[0])
            other_items = items[1:]
            new_ti = master_traces_len + 1 + ((ti-trace_start_num) if (not switched) else (trace_start_num-ti))
            new_items = [str(new_ti)] + other_items + [items[0]]
            new_corlines.append("\t".join(new_items) + "\n")
        return new_corlines

    def _corfile_comment_line(self, filepath):
        '''Returns a string to add as a comment to a corfile, indicating the file
        that the path came from.'''
        return "'===== " + os.path.splitext(os.path.split(filepath)[1])[0] + " =====\n"

    def save_merged_fig(self):
        '''Save the merged figure in 3 rows.'''
        if self.master_traces is None:
            self.merge_all()




class GPR_DataFile:
    '''GPR_DataFile is a class performing operations on a single GPR file.'''

    def __init__(self, filename):
        self.filename = filename
        # Get the resampled corfile if needed.  Else the original.
        if filename == TRANSECT_GPR_MERGED_RESAMPLED:
            self.corfile = TRANSECT_COR_MERGED_RESAMPLED
        elif filename == KANU_GPR_MERGED_RESAMPLED or (filename.find("KANU") >= 0):
            self.corfile = KANU_COR_MERGED_RESAMPLED
        elif filename.find("_resampled") == -1:
            corfile_i = GPR_IMG_FILES.index(self.filename)
            self.corfile = GPR_COR_FILES_PROCESSED[corfile_i]
        else:
            corfile_i = GPR_IMG_FILES_RESAMPLED.index(self.filename)
            self.corfile = GPR_COR_FILES_RESAMPLED[corfile_i]

        self.traces = None
#        self.traces = self.open_image(try_picklefile=True)

    def cor_lines(self):
        corfile = open(self.corfile, 'r')
        lines = [l for l in corfile.readlines() if l.strip() != '']
        corfile.close()
        return lines

    def show_image(self, fname=None, subset_range=None, outfile=None):
#        colormap = [(0,"white"), (.27, "white"), (0.7, "black"), (1,"black")]
#        cmap = colors.LinearSegmentedColormap.from_list("radar", colormap, N=1000)

        self.open_image(fname=fname, try_picklefile=True)

        cmap = self.create_radar_colormap()
        traces, minval, maxval = self.adjust_minmax_traces_to_colormap(self.traces)

        if not (subset_range is None):
            traces = traces[:,subset_range[0]:subset_range[1]]

        plt.imshow(traces, cmap=cmap)
        plt.colorbar()

        if outfile:
            plt.savefig(outfile, dpi=800)

        plt.show()

    def plot_transect_image(self, fname=None):
        traces = self.open_image(fname=fname, try_picklefile=True)
        # Cut off at the first 28000 traces for ACT_transect
        if traces.shape[1] > 28000:
            traces = traces[:,:28000]

        traces, minval, maxval = self.adjust_minmax_traces_to_colormap(traces)

        colormin = 4.55
        colormax = 6.05

        cmap = self.create_radar_colormap(old_minmax=(minval, maxval), new_minmax=(colormin,colormax))

        plt.figure(figsize=(13,3))
        plt.imshow(numpy.clip(traces,colormin,colormax), cmap=cmap,
                   extent=[0,45,20.5,0], aspect=0.20, vmin=colormin, vmax=colormax)
        plt.colorbar(fraction=0.05, pad=0.03, extend='both', label="log$_{10}$ (variance)")
        plt.xlabel("Distance (km)")
        plt.ylabel("Depth (m)")
        plt.tight_layout()
        plt.show()
        for fmt in ("svg", "eps", "tif", "jpg"):
            fname = r"C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\GPR_Main_2013_Transect.{0}".format(fmt)
            print os.path.split(fname)[-1]
            plt.savefig(fname, dpi=600)

    def plot_KANU_image(self, fname=None):
        #TODO: STILL FIX THE IMAGE
        traces = self.open_image(fname=fname, try_picklefile=True)

        traces, minval, maxval = self.adjust_minmax_traces_to_colormap(traces)

        colormin = 4.55
        colormax = 6.05

        length_pixels = traces.shape[1]

        cmap = self.create_radar_colormap(old_minmax=(minval, maxval), new_minmax=(colormin,colormax))

        plt.figure(figsize=(10,3))
        plt.imshow(numpy.clip(traces,colormin,colormax), cmap=cmap,
                   extent=[0,int(length_pixels*1.5/1000.0),20.5,0], aspect=0.10, vmin=colormin, vmax=colormax)
        plt.colorbar(fraction=0.05, pad=0.05, extend='both', label="log$_{10}$ (variance)")
        plt.xlabel("Distance (km)")
        plt.ylabel("Depth (m)")
        plt.tight_layout()
        plt.show()
        for fmt in ("svg", "eps", "tif", "jpg"):
            fname = r"C:\Users\mmacferrin\Dropbox\Research\Papers\MacFerrin et al 2016 - IceBridge Lenses\Figures\KAN_U_Transect.{0}".format(fmt)
            print os.path.split(fname)[-1]
            plt.savefig(fname, dpi=600)

    def open_image(self, fname=None, try_picklefile=False):
        if fname is None:
            fname = self.filename

        if fname == self.filename and not (self.traces is None):
            return self.traces

        if try_picklefile:
            if os.path.splitext(fname)[-1] == ".pickle":
                picklename = fname
            else:
                picklename = os.path.splitext(fname)[0] + ".pickle"

            if os.path.exists(picklename):
                print "Reading", os.path.split(picklename)[-1],
                self.traces = self.read_traces_from_picklefile(picklename)
                print "... traces read."
                return self.traces

        print "Reading", os.path.split(fname)[-1],

        try:
            img = Image.open(fname, mode='r')
            traces = numpy.array(img.getdata())
            del img
            assert ((traces.shape[0] / 2024) * 2024) == traces.shape[0]
            traces.shape = (2024, traces.shape[0]/2024)
        except ValueError:
            # VERY slow, but works correctly for 16-bit signed integer arrays
            traces = tiff.imread(fname)

        print "... traces read."

        if try_picklefile:
            print "Saving", os.path.split(picklename)[-1]
            self.save_traces_as_picklefile(traces, picklename)

        self.traces = traces
        return traces

    def create_variance_image(self, width=3, skip_process=False):

        if not skip_process:
            traces = self.open_image(fname=None, try_picklefile=True)
            assert not (traces is None)

            # Apply gain to traces here.  They still need it!!! (even if it weakens sensitivity at depth)
            traces = self.apply_gain(70,.693/2000.0)

            X,Y = traces.shape
            variances = numpy.zeros(traces.shape, dtype=numpy.float32)
            counts = numpy.zeros(traces.shape, dtype=numpy.uint8)
            means = numpy.zeros(traces.shape, dtype=numpy.float32)

            # Box width, height
            bX = width
            bY = 13 # 13 pixels encompasses an entire radar wavelength in the vertical direction

            # Half-box sizes
            hX = int(bX/2)
            hY = int(bY/2)

            # Compute neighborhood sums in first loop (convert to means immediately after)
            for i in range(-hX,hX+1):
                for j in range(-hY,hY+1):
                    xmin, xmax = max(0,i),min(X,X+i)
                    ymin, ymax = max(0,j),min(Y,Y+j)
                    txmin, txmax = max(0,-i),min(X,X-i)
                    tymin, tymax = max(0,-j),min(Y,Y-j)

                    counts[xmin:xmax,ymin:ymax] += 1
                    # Compute the sum here, then divide by counts after the loop
                    means[xmin:xmax,ymin:ymax] += traces[txmin:txmax,tymin:tymax]

            means = means / counts
            diffmeanssqr = (traces - means)**2

            # Compute
            for i in range(-hX,hX+1):
                for j in range(-hY,hY+1):
                    xmin, xmax = max(0,i),min(X,X+i)
                    ymin, ymax = max(0,j),min(Y,Y+j)
                    txmin, txmax = max(0,-i),min(X,X-i)
                    tymin, tymax = max(0,-j),min(Y,Y-j)

                    variances[xmin:xmax,ymin:ymax] += diffmeanssqr[txmin:txmax,tymin:tymax]

            variances = (variances / counts)
            assert numpy.count_nonzero(variances < 0.0) == 0

            if numpy.count_nonzero(variances == 0.0) > 0:
                print "Correcting", numpy.count_nonzero(variances==0.0), "variances of 0.0."
                variances[variances == 0.0] = numpy.min(variances[variances>0.0])

            log_var = numpy.log10(variances)
    #        sqr_var = numpy.sqrt(variances)

    #        plt.hist(log_var.flatten())
    #        plt.hist(sqr_var.flatten())
    #        plt.show()
    #        foobar

        fname = self.filename
        base, ext = os.path.splitext(fname)
        extender = "_logvariance"
        fname = base + extender + ext

        if not skip_process:
            print "Saving", os.path.split(fname)[1]
            img = Image.fromarray(log_var)
            img.save(fname)

            picklename = os.path.splitext(fname)[0] + ".pickle"
            print "Saving", os.path.split(picklename)[-1]
            self.save_traces_as_picklefile(log_var, picklename)

            self.export_single_image(fname)

        return fname

    def save_traces_as_picklefile(self, traces, fname=None):
        if fname is None:
            fname = os.path.splitext(self.filename)[0] + ".pickle"

        if os.path.splitext(fname)[-1] != ".pickle":
            fname = os.path.splitext(fname)[0] + ".pickle"

        f = open(fname, 'w')
        pickle.dump(traces, f)
        f.close()

        return fname

    def read_traces_from_picklefile(self, fname=None):
        if fname is None:
            fname = os.path.splitext(self.filename)[0] + ".pickle"

        f = open(fname, 'r')
        self.traces = pickle.load(f)
        f.close()
        return self.traces

    def export_single_image(self, infilename = None, extender = None):
        if infilename is None:
            infilename = self.filename

        if extender != None:
            base, ext = os.path.splitext(infilename)
            fname = base + extender + ext
        else:
            fname = infilename
        base, filename = os.path.split(fname)
        outfile = os.path.splitext(os.path.join(EXPORT_DIR, filename))[0] + ".png"

        traces = self.open_image(infilename, try_picklefile=True)


        # Unforunately PyPlot doesn't like to save any image more than 32768 pixels wide
#        if traces.shape[1] > 32768:
#            traces = traces[:,:32768]
        TRACES_CUTOFF_I =28000
        traces = traces[:,:TRACES_CUTOFF_I]

#        plt.hist(traces.flatten())
#        plt.show()
#        foobar

        # aspect, as w/h
#        aspect = float(traces.shape[1]) / float(traces.shape[0])
#
#        fig = plt.figure(figsize=(4*aspect + 1, 5))
#        gs = gridspec.GridSpec(1,1,top=.99, bottom=0.04, left=0.05, right=0.99)
#        ax = plt.subplot(gs[0])
#        kwargs = dict(vmin=0, vmax=numpy.max(traces)*0.15)
#        plt.imshow(traces, **kwargs)


        base, filename = os.path.split(outfile)
        base, dirname = os.path.split(base)

        print "Min:", numpy.min(traces), "Max:", numpy.max(traces)
        # IMG_0129: min 2.9657, max 6.7695
        # TRANSECT: min 2.8206, max 6.9996

        cmap = self.create_radar_colormap()
        traces, minval, maxval = self.adjust_minmax_traces_to_colormap(traces)

        print "Saving", dirname + "/" + filename
        plt.imsave(fname=outfile, arr=traces, cmap=cmap)

        print "Saved."

#        plt.savefig(outfile)
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
#        colortable = [(0,"blue"), (.3, "blue"), (.55, "cyan"), (0.60, "green"), (0.70, "yellow"), (0.8, "red"), (1,"red")]
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

        cmap = colors.LinearSegmentedColormap.from_list("radar", colortable, N=1000)
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
        traces_new = numpy.array(traces)
        # Adjust minimum value in the array
        traces_new[traces_new <= max(defined_MIN, numpy.min(traces_new))] = defined_MIN
        # Adjust maximum value in the array
        traces_new[traces_new >= min(defined_MAX, numpy.max(traces_new))] = defined_MAX

        return traces_new, defined_MIN, defined_MAX

    def export_merged_image(self, merged_image = GPR_MERGED):
        traces = mahotas.imread(merged_image)
        print traces.shape
        std = numpy.std(traces)
        # Get where the file breaks are.  Subtract one to give file indices instead of 1-indexed tracenumbers.
        breaks = numpy.array(self.get_tracenum_breaks_from_corfile())-1
        print std
        aspect = (float(traces.shape[1]) / (float(traces.shape[0])*4)) * 0.35
        figsz = (aspect*8, 8)
        print figsz
        plt.figure(figsize=figsz)
        gs = plt.GridSpec(3,1)
        for i in range(3):
            plt.subplot(gs[i])
            tr = traces[:,traces.shape[1]*i/3:traces.shape[1]*(i+1)/3]
            print i,tr.shape
            plt.imshow(tr, vmin=-(std*3), vmax=std*3)
            for br in breaks:
                if (traces.shape[1]*i/3) <= br <= (traces.shape[1]*(i+1)/3):
                    plt.axvline(x=br-(traces.shape[1]*i/3), color="b", linestyle="-")
    #    print "Showing..."
    #    plt.show()
        print "Saving..."
        plt.savefig(os.path.join(EXPORT_DIR,"ACT_Transect.png"))
        plt.close()

    def get_tracenum_breaks_from_corfile(self,merged_corfile=COR_MERGED):
        '''Run through the merged corfile, get the tracenumbers (1-indexed in the file)
        where each new file begins.  Return a list of these tracenumbers, omitting the first one.'''
        f = open(merged_corfile,'r')
        lines = f.readlines()
        f.close()

        breaks = []
        reached_comment = False
        for line in lines[2:]:
            if line[0] == "'":
                reached_comment = True
                continue
            assert line[0] != "'"
            if not reached_comment:
                continue

            breaks.append(int(line.split()[0]))
            reached_comment = False
        return breaks

    def subtract_horizontal_mean_transform(self):
        '''This helps "destripe" the image from horizontal stripes that still exist
        after the REFLEXW de-wow filter is applied.  It is important to get rid of these
        stripes in order not to have them affect our measurements.  Problem is, this
        "smooths out" sections of the data that would otherwise be noisy enough to NOT
        be ice lenses, such as the surface snow.  What combination of filters needs to be
        applied to detect both?  However, for now, keep it a simple mean subtraction.

        This function takes the base image, subtracts the mean, and saves it as a secondary processed image.'''
        self.open_image()

        trace_means = numpy.mean(self.traces,axis=1)
        traces_mean_corrected = numpy.copy(self.traces)
        for i in range(self.traces.shape[0]):
            traces_mean_corrected[i,:] = self.traces[i,:] - trace_means[i]

        print "Saving", os.path.split(self.img_filename_2)[1]
        mahotas.imsave(self.img_filename_2, traces_mean_corrected)

    def show_image_w_meanplot(self):
        self.open_image()
        gs = gridspec.GridSpec(2,2, width_ratios=(2,1))

        trace_means = numpy.mean(self.traces,axis=1)
        trace_std = numpy.std(trace_means)
        traces_mean_corrected = numpy.copy(self.traces)
        print self.traces.shape
        for i in range(self.traces.shape[0]):
            traces_mean_corrected[i,:] = self.traces[i,:] - trace_means[i]

        ax1 = plt.subplot(gs[0,0])
        ax1
        plt.imshow(self.traces, cmap="gray", vmin=numpy.min(self.traces)*0.4, vmax=numpy.max(self.traces)*0.4)
        ax2 = plt.subplot(gs[1,0])
        ax2
        plt.imshow(traces_mean_corrected, cmap="gray", vmin=numpy.min(self.traces)*0.4, vmax=numpy.max(self.traces)*0.4)
        ax3 = plt.subplot(gs[:,1])

        plt.plot(trace_means, numpy.arange(self.traces.shape[0]))
        plt.axvline(x=trace_std)
        plt.axvline(x=-trace_std)
        ax3.set_ylim((numpy.max(self.traces.shape[0]),0))
        plt.show()


    def display_histogram(self, fignum=2, show=True):
        imgname = self.__get_img_step_filename(fignum)
        plt.hist(mahotas.imread(imgname).flatten(), bins=50)
        plt.title(os.path.split(imgname)[1])
        if show:
            plt.show()


    def display_histogram_split(self, fignum=2, show=True):
        imgname = self.__get_img_step_filename(fignum)
        traces = mahotas.imread(imgname)
        gs = gridspec.GridSpec(2,1)
        plt.subplot(gs[0])
        traces_top = traces[0:traces.shape[0]/2,:]
        plt.hist(traces_top.flatten(),bins=100)
        plt.title(os.path.split(imgname)[1] + " TOP")
        std_top = numpy.std(traces_top)
        plt.axvline(x=-std_top, color='r')
        plt.axvline(x=std_top, color='r')
        plt.draw()
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        plt.text((xmin + (xmax-xmin)*.25),(ymin + (ymax-ymin)*.75), "{0:0.1f}".format(std_top))

        plt.subplot(gs[1])
        traces_bottom = traces[traces.shape[0]/2:traces.shape[0],:]
        plt.hist(traces_bottom.flatten(),bins=100)
        plt.title(os.path.split(imgname)[1] + " BOTTOM")
        std_bottom = numpy.std(traces_bottom)
        plt.axvline(x=-std_bottom, color='r')
        plt.axvline(x=std_bottom, color='r')
        plt.draw()
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        plt.text((xmin + (xmax-xmin)*.25),(ymin + (ymax-ymin)*.75), "std = {0:0.1f}".format(std_bottom))

        if show:
            plt.show()

    def display_plot_mean_std(self, fignum=2, show=True):
        imgname = self.__get_img_step_filename(fignum)
        traces = mahotas.imread(imgname)
        if float(traces.shape[0])/float(traces.shape[1]) >= 0.7:
            gs = gridspec.GridSpec(1,2)
        else:
            gs = gridspec.GridSpec(2,1)
        plt.subplot(gs[0])
        std = numpy.std(traces)
        plt.imshow(traces, vmin=-(std*3), vmax=std*3)
        plt.title(os.path.split(imgname)[1])

        depths = numpy.arange(traces.shape[0])

        ax2 = plt.subplot(gs[1])
        means = numpy.mean(traces,axis=1)
        stds = numpy.std(traces,axis=1)
        plt.plot(means, depths, color="b")
        plt.plot(stds, depths, color="r")

        # Find the linear regression through that line.
        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(depths, stds)
        print "Slope:", slope
        print "Intercept:", intercept
        print "r_value:", r_value
        print "p_value:", p_value
        print "std_err:", std_err

        plt.plot(slope*depths + intercept, depths, "g-")

        ax2.set_ylim((numpy.max(traces.shape[0]),0))

        if show:
            plt.show()


    def display_mean_std_curves(self, fignum=2, show=True):
        if self.filename != GPR_MERGED:
            imgname = self.__get_img_step_filename(fignum)
        else:
            imgname = self.filename
        traces = mahotas.imread(imgname)

        shortname = os.path.split(imgname)[1]
        print shortname

        plt.title(shortname)

        depths = numpy.arange(traces.shape[0])
        means = numpy.mean(traces,axis=1)
        stds = numpy.std(traces,axis=1)
        plt.plot(means, depths, color="b")
        plt.plot(stds, depths, color="r")

        # Find the linear regression through that line.
        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(depths, stds)

        plt.plot(slope*depths + intercept, depths, "g-")

        # Find the exponential regression through that line.
        def expfunc(x, a, b, c):
            return a * numpy.exp(-b * x) + c

        popt, pcov = optimize.curve_fit(expfunc, depths, stds, p0=(800,0.001,0), maxfev=10000)
        print zip(["a","b","c"],popt)
        plt.plot(expfunc(depths,*popt), depths, 'm-')

        plt.ylim((numpy.max(traces.shape[0]),0))

        if show:
            plt.show()

    def apply_gain(self, start_depth, gain_factor, traces=None):

        if traces is None:
            traces = self.open_image()

        # Convert to 32-bit floating point to handle overflow (must cap numbers when converting back to 16-bit integers.)
        tracesf = numpy.array(traces, dtype=numpy.float32)

        for i in range(traces.shape[0]):
            if i <= start_depth:
                continue
            tracesf[i,:] = tracesf[i,:]*numpy.exp((i-start_depth)*gain_factor)

        # Handle overflow and underflow
#        assert traces.dtype == numpy.int16
        tracesf[tracesf > (2**15-1)] = 2**15-1
        tracesf[tracesf < (-2**15)]  = -2**15

        self.traces = numpy.array(tracesf, dtype=numpy.int16)
        return self.traces

    def display_sample_curves(self, numsamples=10, show=True):
        self.open_image()
        print self.traces.shape

    def resample(self, sample_distance_m = 1.5):
        '''Resample the image at every 1.5 meters using a nearest-neighbor approach.
        This filters out "standing still" portions of each traverse, and makes the
        standard deviation convolutions span consistent distances.  Save the file
        when finished.'''
        infilename  = self.img_filename
        outfilename = self.img_filename_resampled
        corobject = GPR_Coords(self.corfilename)

        # First, concentrate on the points themselves, from the corfile.
        trace_resamples_i = corobject.resample(verbose=False)

        # Read data
        traces = self.img_data(fignum=1) # Use the de-striped (but not stddev) image

        print os.path.split(infilename)[-1], traces.shape, "-->", len(trace_resamples_i)
        resampled_traces = numpy.empty([traces.shape[0], len(trace_resamples_i)], dtype=traces.dtype)

        # i refers to the new (resampled) trace index
        # j refers to the old (original) trace index
        for i,old_tracenum in enumerate(trace_resamples_i):
            j = old_tracenum - 1
            resampled_traces[:,i] = traces[:,j]

        print "Writing", os.path.split(outfilename)[-1]
        mahotas.imsave(outfilename, resampled_traces)

        return

def plot_vertical_variance_in_firn(gpr):
    '''In order to adjust for attenuation in firn, we will attempt to plot graphs
    of our known "firn only" gpr files, plotting downward variance of the image
    as a function of depth.  We assume that variance *should* be equal at all
    depths in an all-firn column if gain is set correctly.  We will fit an
    exponential curve to the data (almost all gain functions are exponential) and
    adjust our gain for firn from there, isolating and eliminating one variable
    needed in the calculations for identifying ice lenses (noted above).'''
    assert isinstance(gpr, GPR_DataFile)
    # TODO: Plot vertical variance here!!!!

def MAIN_plot_2013_Transects():
    '''A "main" program for plotting the 2013 GPR transect.
    Created 2015.07.21 '''
#    gprm_transect = GPR_DataFile(GPR_MERGED_RESAMPLED)
#    fname = gprm_transect.create_variance_image(skip_process=True)
#    gprm_transect.plot_transect_image(fname)

    gprm_KANU = GPR_DataFile(KANU_GPR_MERGED_RESAMPLED)
    fname = gprm_KANU.create_variance_image(skip_process=True)
    gprm_KANU.plot_KANU_image(fname)

#    ## Create plots for non-transect images
#    for line_GPR in NON_TRANSECT_GPR_FILES:
#        gpr = GPR_DataFile(line_GPR)
#        gpr.create_variance_image()


def MAIN_plot_cores_next_to_transects(corenum=1):
    '''Main program for plotting core profiles next to GPR transects in two vertical stripes.'''

    #1: Make a 2-panel figure
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(2.8,8),
                                  gridspec_kw={'wspace'      :0.1,
                                               'width_ratios':(1.5,1),
                                               'left':0.21,
                                               'right':1.0-0.21}
                                  )

    core_filename = FirnCoreProfiler.CORE_FILENAMES[corenum-1]

    # TODO: Replace this temporary exported image with a "real" matched trace at KAN_U
    KANU_EXPORTED = r'C:\Users\mmacferrin\Documents\Greenland2013 Data\MALAGS\TraceImages\resampled\exported\DAT_0129_A1_resampled_logvariance.png'
    TRANSECT_EXPORTED = r'C:\Users\mmacferrin\Documents\Greenland2013 Data\MALAGS\TraceImages\resampled\exported\ACT_Transect_resampled_logvariance.png'
    img_filename = KANU_EXPORTED if corenum in (1,2) else TRANSECT_EXPORTED

    gpr_filename = KANU_GPR_MERGED_RESAMPLED if corenum in (1,2) else TRANSECT_GPR_MERGED_RESAMPLED

    #1.5: Retreive core info
    cp = FirnCoreProfiler.CoreProfile(core_filename)
    print cp.filename
    #3: Fill right panel in with ice lens stratigraphy (same width)
    cp.plot_stratigraphy_only(ax=ax2)

    # TODO: Get "N" GPR traces adjacent to core.
    gpr = GPR_DataFile(gpr_filename)
    lat, lon = cp.latitude, cp.longitude
    gprc = GPR_Coords(gpr.corfile)
    closest_trace = gprc.closest_traces((lat,lon), numtraces=1)[0][1]
    # TODO: Replace this temporary exported image with a "real" matched trace at KAN_U
    if corenum in (1,2):
        closest_trace = 175

    print closest_trace
#    gpr.open_image()
    N = 300
    trace_distance = N*1.5
    img = mpimage.imread(img_filename)
    print img.shape
    img_traces = img[:, (closest_trace - N/2) : (closest_trace + N/2), :]

    print img_traces.shape
    ax1.imshow(img_traces, aspect="auto", extent=(0,trace_distance,20,0))
    ax1.plot((trace_distance/2.0, trace_distance/2.0),(0,20), color="black", lw=2.0)

    # TODO: Fill left panel in with GPR profile (N pixels wide)

    # TODO: Adjust depths?

    # Set axis limits and labels
    ax1.set_xlim(0,trace_distance)
    ax1.set_ylim(20,0)
    ax1.set_ylabel("Radar Depth (m)")
    ax1.set_xlabel("Distance (m)")
    ax1.set_xticks([0,200,400])
    ax1.yaxis.tick_left()
    ax2.set_xticks([])
    ax2.set_ylabel("Core Depth (m)")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    plt.show()

if __name__ == "__main__":


    KANU_fname = [f for f in GPR_DETRENDED_PICKLEFILES if f.find("KANU_") >= 0][0]
    gpr = GPR_DataFile(KANU_fname)

#    MAIN_plot_2013_Transects()
#    pass

#    MAIN_plot_cores_next_to_transects(1)
#    MAIN_plot_2013_Transect()

#    ## 2 lines: Merging KAN_U Transect
#    gprm = GPR_Merger(TRANSECT_COLLECTION = KANU_TRANSECT_COLLECTION, resampled=True)
#    gprm.merge_all(GPR_outfile=KANU_MERGED_RESAMPLED,COR_outfile=KANU_COR_MERGED_RESAMPLED)

    ## Now, open and plot KAN_U Transect
#    gpr = GPR_DataFile(KANU_MERGED_RESAMPLED)
#    gpr.create_variance_image()


#    f = [f for f in GPR_IMG_FILES_RESAMPLED if f.find("0129") != -1][0]
#    gpr = GPR_DataFile(f)
#
#    fname = gpr.create_variance_image(skip_process=True)
#    gpr.plot_KANU_image(fname)

#    gprm = GPR_DataFile(GPR_MERGED_RESAMPLED)
#
#    fname = gprm.create_variance_image(skip_process=True)
#    gprm.plot_transect_image(fname)

#    gpri = GPR_IceLensIdentifier()
#    ct = gpri.find_closest_traces()
#    for core in ct.keys():
#        print "===", os.path.split(core)[1], "==="
#        for trace in ct[core]:
#            print trace
#        print

#    for f in GPR_IMG_FILES:
#        gpr = GPR_DataFile(f)
#        gpr.resample()
#        gpr.display_sample_curves()

#    gpr = GPR_DataFile(GPR_MERGED)
#    gpr.apply_gain(70,.693/2000.0)
#    gpr.create_std_images()

#    gpr.export_merged_image()
#    gpr.display_mean_std_curves()

#    gpr = GPR_DataFile(GPR_IMG_FILES[0])
#    gpr.create_std_image()
##    gpr.display_histogram()
#    gpr.show_image_w_meanplot()
#    for f in GPR_IMG_FILES:
#        gpr = GPR_DataFile(f)
##        gpr.export_single_image(1)
#        gpr.create_std_images()
#        #gpr.export_all_images()
#        del gpr

#    gprc = GPR_Collection()
#    gprc.merge_whole_transect()
#    gprc.display_all_std_curves()
#    gprc.import_all_images()
#    gprc.display_all_histograms()

#    gprm = GPR_Merger()
#    gprm.merge_all()