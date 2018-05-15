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
        self.slit_name = ''
        self.skip = 0 #skip some number of slits
        self.position_angle = np.nan
        self.approx_bottom = 0.0
        self.approx_top = 0.0
        self.approx_target = 0.0
        self.width_arcsec = np.nan
        self.approx_R = 0.0
        self.range_lambda0 = [0.0,0.0]
        self.range_delta_lambda = [0.0,0.0]


    def create_slit_structure(self,fuel):

        smf = fuel.input.mask_SMF
