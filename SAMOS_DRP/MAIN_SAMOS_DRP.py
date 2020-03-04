from __future__ import print_function
import sys
import os
import json
from astropy.io import fits
from astropy import units
from CreateInput import CreateInput
from DoOverscan import Overscan,Bias
import glob
import pandas as pd
import ccdproc
from ccdproc import CCDData

#from SAMOSHelpers import IsItHere

class MAIN_SAMOS_DRP:

    def __init__(self,obsid,workdir,raw_data_dir):

        self.obsid = obsid
        self.working_dir = workdir
        self.raw_data_dir = raw_data_dir
        self.bias_filelist = None
        self.flat_filelist = None
        self.comp_filelist = None
        self.science_filelist = None
        self.field_filelist = None
        self.grating_list = None
        self.use_bias_frames = False
        self.buckets = []
        self.full_data_bucket = None
        self.header_keys = ['naxis','date','slit','date-obs','obstype','object',
                    'exptime','obsra','obsdec','grating','cam_targ','grt_targ',
                    'filter','filter2','gain','rdnoise','roi','wavmode']



    def __call__(self):


        DRPin = CreateInput(self.obsid,self.working_dir,self.raw_data_dir)
        grouped_data = DRPin.group_spec_obs()
        print(DRPin.raw_data_dir)
        #create input information for the observations you want to process.
        #obsid is the date of observations in the format YYYYMMDD
        self.raw_data_dir = DRPin.raw_data_dir
        self.flat_filelist = DRPin.flat_filelist
        self.comp_filelist = DRPin.comp_filelist
        self.science_filelist = DRPin.science_filelist
        self.field_filelist = DRPin.field_filelist
        self.bias_filelist = DRPin.bias_filelist
        self.grating_list = DRPin.grating_list


        if not os.path.exists('%s/%s_processing'%(self.working_dir,self.obsid)):
            os.mkdir('%s/%s_processing'%(self.working_dir,self.obsid))

        self.processing_dir = '%s/%s_processing'%(self.working_dir,self.obsid)

        self.buckets = DRPin.data_bucket

    def overscan_and_trim(self):

        inlamps = self.comp_filelist
        inflats = self.flat_filelist
        intargs = self.science_filelist
        infields = self.field_filelist
        biases = self.bias_filelist
        corr_dir = self.processing_dir


        #make the master bias frame so it can be subtracted from the data
        if self.use_bias_frames==True:
            mb = Bias(biases,corr_dir)
        else:
            mb = None

        #how the outfiles are named depends on what the basefile name of the data is.
        #for now the syntax is based on the Goodman MOS test data that I have.
        outlamps = ["%s/b_%s" % (corr_dir,os.path.basename(i)[11:]) for i in inlamps]
        outflats = ["%s/b_%s" % (corr_dir,os.path.basename(i)[11:]) for i in inflats]
        outtargs = ["%s/b_%s" % (corr_dir,os.path.basename(i)[11:]) for i in intargs]
        outfield = ["%s/b_%s" % (corr_dir,os.path.basename(i)[11:]) for i in infields]

        [Overscan(i,o,mb) for i,o in zip(inlamps,outlamps)]
        [Overscan(i,o,mb) for i,o in zip(inflats,outflats)]
        [Overscan(i,o,mb) for i,o in zip(intargs,outtargs)]

        db = open('%s'%self.inputRAW.db_file,"a")
        db.write("\n## Trimmed/Overscan Subtracted\n")
        db.write("# OSLamps\n")
        [db.write("%s\n"%(olamp)) for olamp in outlamps]
        db.write("# OSFlats\n")
        [db.write("%s\n"%(oflat)) for oflat in outflats]
        db.write("# OSTargs\n")
        [db.write("%s\n"%(otarg)) for otarg in outtargs]
        db.close()

        #print("saving this step as OverscanTrim")

        return self
