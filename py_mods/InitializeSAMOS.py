import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
from CreateFuelStructure import CreateFuelStructure
from CreateInput import CreateInput

class InitializeSAMOSInstrument:

    def __init__(self):
        self.instrument_name = ''
        self.grating = ''
        self.grating_order = ''
        self.central_wavelength: ''
        self.pixel_scale: ''
        self.filter: ''
        self.readnoise: ''
        self.gain: ''
        self.resolution_slit1arcsec: 0.0
        self.linearity_correction: ''
        self.trim_edges: 4.0
        self.default_badpixel_mask: ''
        self.default_dark: ''
        self.default_pixelflat: ''
        self.default_illumflat: ''
        self.default_arc: ''

    def initialize_SAMOS_instrument(self,science_frame):

        fname = fits.open(science_frame)
        hdr = fname[0].header
        date = hdr['DATE-OBS']
        self.instrument_name = hdr['INSTRUME']
        self.grating = hdr['GRISM']
        data_dir = os.path.split(science_frame)[0]
        obsdef = open('%s/L3VPH_all.obsdef'%(data_dir),'r')
        for line in obsdef.readline():
            if 'GR_ORDER' in line:
                self.grating_order = line.split(' ')[1]

        if self.grating == 'VPH-Red':
            self.central_wavelength = 8000 #angstroms
        elif self.grating == 'VPH-Blue':
            self.central_wavelength = 5000
        elif self.grating == 'VPH-All':
            self.central_wavelength = 7100

        self.pixel_scale = hdr['SCALE'] #arcsec/pixel
        self.filter = hdr['FILTER']
        self.readnoise = hdr['ENOISE']
        speed = hdr['SPEED']
        if speed == 'Slow':
            if fname[-2:]=='c1':
                self.gain = 0.15 #e-/ADU
            elif fname[-2:]=='c2':
                self.gain = 0.13
        elif speed == 'Fast':
            if fname[-2:]=='c1':
                self.gain = 1.67
            elif fname[-2:]=='c2':
                self.gain = 1.45
        elif speed == 'Turbo':
            if fname[-2:]=='c1':
                self.gain = 3.04
            elif fname[-2:]=='c2':
                self.gain = 2.63
        if int(date[:4])>=2005:
            self.default_badpixel_mask = '%s/SITE2'%(data_dir)
        else:
            self.default_badpixel_mask = '%s/SITE'%(data_dir)

        return self


def initialize_SAMOS(datedir,mask):

    input_structure = CreateInput().create_input(datedir,mask)

    fuel = CreateFuelStructure().create_fuel_structure(input_structure)

    instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(input_structure.science_filelist[0])

    fuel.instrument = instrument
    #instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(fuel.input.science_filelist)
    return fuel
