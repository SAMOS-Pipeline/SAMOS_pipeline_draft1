from __future__ import print_function
import sys
import os
#sys.path.insert(0,os.path.split(os.getcwd())[0])
from pyraf import iraf
import iraf_params as irf_par


def initialize_iraf():
    #initialize IRAF and tasks needed.
    #give user the option to edit task parameters??
    iraf.noao(_doprint=0)
    #iraf.imred(_doprint=0)
    #iraf.ccdred(_doprint=0)
    #iraf.specred(_doprint=0)
    #iraf.onedspec(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.apextract(_doprint=0)
    iraf.longslit(_doprint=0)
    return

def auto_identify_arcs(input_arc,fcname,logfile='logfile'):
    """
    uses the iraf task noao.twodspec.longslit.identify to get a wavelength solution using
    the comparison arclamp images for each slit after they've been cut out.  Called by the
    PyrafIdentify module.

    irf_par.set_aidpars_calibration(iraf.aidpars)
    irf_par.set_autoidentify_calibration(iraf.autoidentify)
    irf_par.set_reidentify_calibration(iraf.reidentify)
    irf_par.set_fitcoords_calibration(iraf.fitcoords)
    """

    iraf.autoidentify(images=input_arc,crval=8000,cdelt=1.175)

    check_fit = os.popen("tail -n 3 logfile").read().strip().split('\n')
    #print(check_fit)
    #print(check_fit[-1])
    if 'No solution found' in check_fit[-1]:
        print("No solution found in %s"%(input_arc))
        iraf.aidpars.crpix = 685
        print(iraf.aidpars.crpix)
        #iraf.autoidentify(images=input_arc)
        return False
    else:

        return True



def reident(input_arc,fcname,logfile='logfile'):

    iraf.reidentify(reference=input_arc,images=input_arc)
    iraf.fitcoords(images=input_arc[:-5],fitname=fcname)


    return

def transform(input_targ,output_targ,fcname):
    irf_par.set_transform(iraf.transform)

    iraf.transform(input=input_targ,output=output_targ,fitnames=fcname)
    return
