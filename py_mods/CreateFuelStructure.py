import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
import struct


class CreateFuelLoadfiles:

    def __init__(self):
        self.n_frames = 0
        self.raw_files = ''
        self.corr_files = ''
        self.master_file = ''



    def create_fuel_loadfiles(self,filelist):

        '''
        Creating a data structure for storage of information about the observation to be
        referenced easily by the pipeline.  Basically trying to translate the Flame DRP to python.
        '''

        #print("input is: %s"%(dir(input)))



        self.n_frames = len(filelist)
        self.corr_files = ['%s_corr'%(os.path.basename(fname)) for fname in filelist]
        self.raw_files = filelist

        return self


class CreateUtilStructure:

    #util substructure for referencing files and their locations.
    def __init__(self):
        self.science = ''
        self.arc = ''
        self.slitflat = ''
        self.intermediate_dir = ''


    def create_util_structure(self,input):

        #Make the directory withing LMask*/ to put the newly corrected frames after each step
        if not os.path.exists('%s/%s'%(input.slit_mask,'corr_frames')): os.mkdir('%s/%s'%(input.slit_mask,'corr_frames'))

        science = CreateFuelLoadfiles().create_fuel_loadfiles(input.science_filelist)
        print('%s science frames read'%(science.n_frames))

        calib_arc = CreateFuelLoadfiles().create_fuel_loadfiles(input.arc_filelist)
        calib_arc.master_file = '%s/%s'%(input.slit_mask,'master_arc.fits')
        calib_slitflat = CreateFuelLoadfiles().create_fuel_loadfiles(input.flat_filelist)
        calib_slitflat.master_file = '%s/%s'%(input.slit_mask,'master_slitflat.fits')

        self.science = science
        self.arc = calib_arc
        self.slitflat = calib_slitflat
        self.intermediate_dir = input.slit_mask

        return self
'''
class CreateSettingStructure:

        def __init__(self):



            settings = { $
    sky_emission_filename: data_dir + 'sky_emission_model_nir.dat', $
    linelist_sky_filename: '', $
    linelist_arcs_filename: '', $
    star_x_range: [1000, 1200], $
    star_y_window: 0, $
    clean_individual_frames: 0, $
    badpix_useflat: 1, $
    badpix_usedark: 1, $
    badpix_sigma: 7.0, $
    badpix_flatcorrection : 0.20, $
    flatfield_data: 1, $
    darksub_data: 1, $
    trace_slit_with_emlines : 0, $
    trace_slit_xmargin: 20, $
    trace_slit_ymargin: 12, $
    trace_slit_polydegree: 2, $
    trace_longslit: 0, $
    trim_slit: [0, 0], $
    spatial_resampling: 1.0, $
    roughwavecal_R : [500, 1000, 3000], $
    roughwavecal_smooth_window: 20, $
    roughwavecal_split: 0, $
    findlines_stack_rows: 0, $
    findlines_poly_degree: 5, $
    findlines_Nmin_lines: 6, $
    findlines_linefit_window: 6.0, $       ; in units of expected linewidth
    wavesolution_order_x: 3, $
    wavesolution_order_y: 2, $
    wavecal_sigmaclip: 3.0, $
    shift_arcs: 1, $
    shift_arcs_Nmin_lines: 1, $
    illumination_correction : 1, $
    skysub: 1, $
    skysub_plot: 0, $
    skysub_plot_range: [0.4, 0.6], $
    skysub_bspline_oversample: 1.0, $
    skysub_reject_fraction: 0.10, $
    skysub_reject_loops: 3, $
    skysub_reject_window: 2.0, $
    interpolation_method: 'NaturalNeighbor', $
    frame_weights: 'None', $
    combine_sigma_clip : 2.0, $
    combine_min_framefrac : combine_min_framefrac, $
    extract_optimal : 1, $
    extract_gaussian_profile : 1, $
    stop_on_error:1 $
   }

'''

class CreateFuelStructure:

    def __init__(self):
        self.input = None
        self.settings = None
        self.util = None
        self.instrument = None
        self.slits = None


    def create_fuel_structure(self,input):

        self.input = input
        self.util = CreateUtilStructure().create_util_structure(input)

        return self
