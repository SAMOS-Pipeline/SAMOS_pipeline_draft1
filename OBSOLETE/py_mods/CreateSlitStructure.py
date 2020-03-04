from __future__ import print_function
import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
import struct


class CreateSlitStructure:
    #creating empty slit structure for fuel structure.

    def __init__(self):
        self.number = 0
        self.obj = None
        self.science = None
        self.arc = None
        self.flats = None
        self.skip = 0 #skip some number of slits
        self.position_angle = np.nan
        self.x_mm = 0.0
        self.y_mm = 0.0
        self.obj_ra  = ''
        self.obj_dec = ''
        self.approx_bottom = 0.0
        self.approx_top = 0.0
        self.approx_target = 0.0
        self.width_arcsec = 0.0
        self.width_mm = np.nan
        self.alen_mm = 0.0
        self.blen_mm = 0.0
        self.angle = np.nan
        self.approx_R = 0.0
        self.cent_lambda = 0.0
        self.range_lambda0 = [0.0,0.0]
        #self.range_delta_lambda = [0.0,0.0]


    def create_slit_structure(self,maskSMF):

        smf = open(maskSMF,'r')
        for line in smf.readlines():
            if 'WLIMIT' in line:
                self.range_lambda0 = [float(line.split(' ')[1].strip()),
                float(line.split(' ')[2].strip())]
            elif 'WAVELENGTH' in line:
                self.cent_lambda = float(line.split(' ')[1].strip())

        return self
