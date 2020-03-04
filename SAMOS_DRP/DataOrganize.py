from __future__ import print_function
import sys
import os
import json
#print('Main working directory is: %s'%(sys.path[0]))
from SAMOSHelpers import check_header_notes
from astropy.io import fits
from astropy.coordinates import SkyCoord
import glob
import pandas as pd
import numpy as np
import ccdproc
from ccdproc import CCDData, ImageFileCollection


class DataOrganize(object):

    self.working_dir = workdir
    self.obsid = obsid
    self.bias_filelist = None
    self.flat_filelist = None
    self.comp_filelist = None
    self.science_filelist = None
    self.field_filelist = None
    self.raw_data_dir = ''
    self.groups = None
    self.header_keys = ['naxis','date','slit','date-obs','obstype','object',
                'exptime','obsra','obsdec','grating','cam_targ','grt_targ',
                'filter','filter2','gain','rdnoise','roi','wavmode']
