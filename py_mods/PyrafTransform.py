from __future__ import print_function
import glob
from astropy.io import fits
from SAMOSHelpers import *
import sys
import os
from pyraf import iraf
import iraf_steps as irf_stp
import iraf_params as irf_par
import shutil


def pyraf_trnsf(fuel):

    dbdir = fuel.identify.db
    intargs = fuel.util.science.corr_files
    transf_dir = '%s/%s'%(fuel.input.slit_mask,'transformed_targs')
    if not os.path.exists(transf_dir): os.mkdir(transf_dir)
    wvcal_targs = []
    for fc in glob.glob('%s/fcslit*'%(dbdir)):
        slitn = os.path.basename(fc).split('.')[0][2:]
        for targ in intargs:
            targslitn,targnm = os.path.basename(targ).split('_')
            if targslitn==slitn:
                print(targslitn,slitn)
                fcnm = os.path.basename(fc)[2:]
                if os.path.exists('current_targ.fits'): os.system('rm current_targ.fits')
                os.symlink(targ,'current_targ.fits')

                outtargnm = '%s.t%s'%(targslitn,targnm)
                irf_stp.transform('current_targ.fits',outtargnm,fcnm)

                wvcal_targ = '%s/%s'%(transf_dir,outtargnm)
                os.rename(outtargnm,wvcal_targ)
                wvcal_targs.append(wvcal_targ)

                os.rename('current_targ.fits','previous_targ.fits')
                os.unlink('%s/%s'%(os.getcwd(),'previous_targ.fits'))

    fuel.util.science.corr_files = wvcal_targs

    return fuel
