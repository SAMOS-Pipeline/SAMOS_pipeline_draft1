from __future__ import print_function
import glob
from astropy.io import fits
import numpy as np
import pickle
import os
import sys

def save_fuel_step(this_fuel,fuel_fname):

    """
    This function takes an fuel instance created or updated from one of the data reduction steps
    and dumps it into a file for future access.  This will useful for if the user wants to stop
    the pipeline and pick back up instead of restarting the reduction process from scratch.
    """
    dump_dir = this_fuel.fuelsave_dir
    filename = '%s/%s_dump.txt'%(dump_dir,fuel_fname)
    outfile = open(filename, 'w')
    try:
        pickle.dump(this_fuel,outfile)
        print('WRITTEN: this reduction step %s/%s '%(os.path.basename(dump_dir),os.path.split(filename)[1]))

    except:
        print('WARNING: could not write this fuel step to file %s/%s'%(os.path.basename(dump_dir),os.path.split(filename)[1]))

    outfile.close()

    return

#def refuel_data(dbfile):
