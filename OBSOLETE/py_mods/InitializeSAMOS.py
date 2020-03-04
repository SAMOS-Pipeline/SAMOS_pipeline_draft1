from __future__ import print_function
import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
from CreateFuelStructure import CreateFuelStructure
from CreateInput import CreateInput
from CreateSlitStructure import CreateSlitStructure
from SaveFuel import save_fuel_step


class InitializeSAMOSInstrument:

    def __init__(self):
        self.date = ''
        self.instrument = ''
        self.mask = ''
        self.epoch = 0.0
        self.temp = 0.0
        self.dref = 0.0
        self.hangle = 0.0
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
        self.resolution = 0.0
        self.resolution_slit1arcsec= 0.0
        self.linearity_correction= [0.0,1.0,4e-6] #copied from FLAME, no LDSS3 linearity info
        self.linear_disp = 0.0
        self.scale_arc_per_mm = 0.0
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
            self.linear_disp = 1.175 #ang/pix
            self.resolution = 1810 #0.75" slit
        elif self.grism == 'VPH-Blue':
            self.central_wavelength = 5000
            self.linear_disp = 0.682
            self.resolution = 1900
        elif self.grism == 'VPH-All':
            self.central_wavelength = 7100
            self.linear_disp = 1.890
            self.resolution = 860


        self.pixel_scale = hdr['SCALE'] #arcsec/pixel
        self.scale_arc_per_mm = 12.6
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
                self.celestial_ra = line.split(' ')[2].strip()
                self.celestial_dec = line.split(' ')[3].strip()
            elif "TEMPERATURE" in line:
                self.temp = float(line.split()[1])
            elif "EPOCH" in line:
                self.epoch = float(line.split()[1])
            elif "DREF" in line:
                self.dref = float(line.split()[1])
            elif "HANGLE" in line:
                self.hangle = float(line.split()[1])

        obsdef.close()
        return self


def initialize_SAMOS_slits(input_mask):

    slits = []
    this_slit = CreateSlitStructure().create_slit_structure(input_mask)
    """
    mask_file = open(input_mask, 'r')
    slit_num = 1
    for line in mask_file:
        line = line.split()
        if "SLIT" in line:

            this_slit.number = slit_num
            this_slit.obj = line[1]
            this_slit.obj_ra,this_slit.obj_dec = line[2],line[3]
            this_slit.x_mm,this_slit.y_mm = float(line[4]),float(line[5])
            this_slit.alen_mm,this_slit.blen_mm = float(line[7]),float(line[8])
            this_slit.width_mm = float(line[6])
            this_slit.angle = float(line[9])

            slit_num+=1

            slits.append(this_slit)
    """
    return np.asarray(slits)


def initialize_SAMOS(datedir,mask):

    input_structure = CreateInput().create_input(datedir,mask)

    data_dump_dir = '%s/stored_%s'%(input_structure.working_dir,mask)

    fuel = CreateFuelStructure().create_fuel_structure(input_structure,data_dump_dir)

    slit_init = initialize_SAMOS_slits(input_structure.mask_SMF)

    instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(input_structure.science_filelist[0],input_structure.mask_SMF)

    fuel.instrument = instrument
    #instrument = InitializeSAMOSInstrument().initialize_SAMOS_instrument(fuel.input.science_filelist)
    #fuel.slits = slit_init

    save_fuel_step(fuel,'Initialize')

    return fuel
