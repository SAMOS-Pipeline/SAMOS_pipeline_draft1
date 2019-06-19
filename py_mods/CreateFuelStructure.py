from __future__ import print_function
import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
import struct
from CreateSlitStructure import CreateSlitStructure


class CreateFuelLoadfiles:

    def __init__(self):
        self.n_frames = 0
        self.raw_files = ''
        self.corr_files = ''
        self.transf_files = ''
        self.master_file = ''



    def create_fuel_loadfiles(self,filelist):

        '''
        Creating a data structure for storage of information about the observation to be
        referenced easily by the pipeline.  Basically trying to translate the Flame DRP to python.
        '''

        #print("input is: %s"%(dir(input)))



        self.n_frames = len(filelist)
        #self.corr_files =
        self.raw_files = filelist

        return self


class CreateUtilStructure:

    #util substructure for referencing files and their locations.
    def __init__(self):
        self.science = ''
        self.arc = ''
        self.slitflat = ''
        self.field = ''
        self.intermediate_dir = ''


    def create_util_structure(self,input):


        #Make the directory withing LMask*/ to put the newly corrected frames after each step
        #if not os.path.exists('%s/%s'%(input.slit_mask,'corr_frames')): os.mkdir('%s/%s'%(input.slit_mask,'corr_frames'))

        science = CreateFuelLoadfiles().create_fuel_loadfiles(input.science_filelist)
        print('%s science frames read'%(science.n_frames))

        calib_arc = CreateFuelLoadfiles().create_fuel_loadfiles(input.arc_filelist)
        calib_arc.master_file = '%s/%s'%(input.slit_mask,'master_arc.fits')
        calib_slitflat = CreateFuelLoadfiles().create_fuel_loadfiles(input.flat_filelist)
        calib_slitflat.master_file = '%s/%s'%(input.slit_mask,'master_slitflat.fits')
        field = CreateFuelLoadfiles().create_fuel_loadfiles(input.field_filelist)

        self.science = science
        self.arc = calib_arc
        self.slitflat = calib_slitflat
        self.field = field
        self.intermediate_dir = '%s/%s'%(input.slit_mask,'intermediates')
        if not os.path.exists(self.intermediate_dir): os.mkdir(self.intermediate_dir)
        #intermediate_dir is the main directory where
        #fuel.util.(science,arc,slitflat).corr_files live


        return self

class CreateIdentifyStructure:

        def __init__(self):
            self.lgfile = '' #logfile written to by PyRAF
            self.db = '' #database dir for identify output
            self.linelist_ext = 'linelist_full' #full linelist for autoidentify
            self.linelist = 'linelist' #concatenated linelist for aidpars

        def create_id(self):

            if os.path.exists('logfile'): os.system('rm logfile')

            os.system('touch logfile')

            if not os.path.exists('database'): os.mkdir('database')

            self.lgfile = os.path.join(os.getcwd(),'logfile')
            self.db = os.path.join(os.getcwd(),'database')

            return self

class CreateFuelStructure:

    def __init__(self):
        self.input = None
        self.util = None
        self.instrument = None
        self.slits = []
        self.identify = None
        self.fuelsave_dir = None


    def create_fuel_structure(self,input,dump_dir):
        if not os.path.exists(dump_dir):os.mkdir(dump_dir)
        
        self.fuelsave_dir = dump_dir
        self.input = input
        self.util = CreateUtilStructure().create_util_structure(input)
        self.identify = CreateIdentifyStructure().create_id()
        #self.slits = CreateSlitStructure().create_slit_structure(input)


        return self
