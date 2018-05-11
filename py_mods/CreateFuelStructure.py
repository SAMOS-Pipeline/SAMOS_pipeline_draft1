import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
import struct


class CreateFuelStructure:

    def __init__(self):
        self.input = None
        self.setings = None
        self.util = None
        self.instrument = None
        self.slits = None

    

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
