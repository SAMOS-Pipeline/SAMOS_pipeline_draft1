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
        self.date = ''
        self.instrument = ''
        self.mask = ''
        self.dewoff = ''
        self.camera = ''
        self.mode = ''
        self.ns = ''
        self.dewar = ''
        self.grism = ''
        self.gr_order = ''
        self.gr_angle = ''
        self.d_alignrot = ''
        self.central_wavelength= ''
        self.pixel_scale= ''
        self.filter= ''
        self.readnoise= ''
        self.gain= ''
        self.celestial_ra = ''
        self.celestial_dec = ''
        self.resolution_slit1arcsec= 0.0
        self.linearity_correction= ''
        self.trim_edges= 4.0
        self.default_badpixel_mask= ''
        self.default_dark= ''
        self.default_pixelflat= ''
        self.default_illumflat= ''
        self.default_arc= ''


    def initialize_SAMOS_instrument(self,science_frame,mask):

        fname = fits.open(science_frame)
        hdr = fname[0].header
        self.date = hdr['DATE-OBS']
        self.instrument = hdr['INSTRUME']
        self.grism = hdr['GRISM']
        data_dir = os.path.split(science_frame)[0]
        obsdef = open('%s/L3VPH_all.obsdef'%(data_dir),'r')
        for line in obsdef.readline():
            if 'MODE' in line:
                self.mode = line.split(' ')[1]
            elif 'DEWAR' in line:
                self.dewar = line.split(' ')[1]
            elif 'DEWOFF' in line:
                self.dewoff = line.split(' ')[1]
            elif 'GR_ORDER' in line:
                self.gr_order = line.split(' ')[1]
            elif 'GR_ANGLE' in line:
                self.gr_angle = line.split(' ')[1]
            elif 'D_ALIGNROT' in line:
                self.d_alignrot = line.split(' ')[1]

        if self.grism == 'VPH-Red':
            self.central_wavelength = 8000 #angstroms
        elif self.grism == 'VPH-Blue':
            self.central_wavelength = 5000
        elif self.grism == 'VPH-All':
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
        if int(self.date[:4])>=2005:
            self.default_badpixel_mask = '%s/SITE2'%(data_dir)
        else:
            self.default_badpixel_mask = '%s/SITE'%(data_dir)

        for line in open(mask,'r'):
            if "POSITION" in line:
                print('found celestial coords')
                print(line)
                print(line.split(' ')[2])
                self.celestial_ra = line.split(' ')[2].strip()
                self.celestial_dec = line.split(' ')[3].strip()


        return self


def initialize_SAMOS_slits(input,instrument):

    mask_file = open(input.mask_SMF, 'r')
    for line in mask_file:
        pass

    slit_obj = np.genfromtxt(input.mask_SMF,skip_header=13,skip_footer=9,
                    usecols=1,unpack=True,dtype=str)
    x_mm,y_mm = np.genfromtxt(input.mask_SMF,skip_header=13,skip_footer=9,
                    usecols=(4,5),unpack=True,dtype=float)
    width_mm = np.genfromtxt(input.mask_SMF,skip_header=13,skip_footer=9,
                    usecols=6,unpack=True,dtype=float)
    len_mm_neg,len_mm_pos = np.genfromtxt(input.mask_SMF,skip_header=13,skip_footer=9,
                    usecols=(7,8),unpack=True,dtype=float)

    width_arcsec = width_mm





def initialize_SAMOS(datedir,mask):

    input_structure = CreateInput().create_input(datedir,mask)

    fuel = CreateFuelStructure().create_fuel_structure(input_structure)

    instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(input_structure.science_filelist[0],input_structure.mask_SMF)

    fuel.instrument = instrument
    #instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(fuel.input.science_filelist)



    return fuel
