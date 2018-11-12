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

"""

This module is meant to use PyRAF to identify the features in the comparison lamp.
If you have activated the correct SAMOS working environment for the current branch of the
pipeline, this step should work.
If not, the user should activate the correct environment by typing
' conda env create -f EnvSAMOS.yml ' into their terminal while in the working directory.

"""

def pyraf_identify(fuel):

    irf_stp.initialize_iraf()

    arc_list = fuel.util.arc.corr_files
    arc_fits = [fits.open(arc) for arc in arc_list]
    arc_file_names = [os.path.split(arc)[1][:-9] for arc in arc_list]
    #get a general name for the arc files that strips the 81**.fits from the file names.
    obj_names = [obj[0].header["OBJECT"] for obj in arc_fits]
    fcnames = ["%s%s"%(arc_slit.replace('_','.'),nm[6:].replace(' ','.')) for arc_slit,nm \
                in zip(arc_file_names,obj_names)]

    identified_arcs = []
    for arc in range(len(arc_list)):
        arcdir,arcfile = os.path.split(arc_list[arc])

        if arcfile in identified_arcs:
            continue
        if os.path.exists('current_arc.fits'): os.system('rm current_arc.fits')
        os.symlink(arc_list[arc],'current_arc.fits')
        #os.chdir('%s/intermediates/slit_cutouts/slit_%s'% \
        #        (fuel.input.slit_mask,arc+1))
                #go to slit directory for the arc.

        #if not os.path.exists('linelist'): shutil.copyfile('%s/linelist'%(os.getcwd()),'linelist')
        #if not os.path.exists('database'): os.mkdir('database')
        #makes database directory for iraf to use for dispersion solutions.
        irf_par.set_aidpars_calibration(iraf.aidpars)
        irf_par.set_autoidentify_calibration(iraf.autoidentify)
        irf_par.set_reidentify_calibration(iraf.reidentify)
        irf_par.set_fitcoords_calibration(iraf.fitcoords)
        print(arcfile)

        iraf.autoidentify.interactive = 'NO'
        auto_solution = irf_stp.auto_identify_arcs('current_arc.fits',fcnames[arc])
        if auto_solution == True:
            irf_stp.reident('current_arc.fits',fcnames[arc])
        elif auto_solution == False:
            pass

        os.rename('current_arc.fits','previous_arc.fits')
        os.unlink('%s/%s'%(os.getcwd(),'previous_arc.fits'))
        identified_arcs.append(arcfile)


    return fuel
