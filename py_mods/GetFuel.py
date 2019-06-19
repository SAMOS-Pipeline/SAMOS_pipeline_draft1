from __future__ import print_function
import glob
from astropy.io import fits
import numpy as np
import pickle
import os
import sys
from InitializeSAMOS import initialize_SAMOS


def get_fuel_step(fuelsave_dir,which_fuel):

    fname = '%s/%s.txt'%(fuelsave_dir,which_fuel)
    infile = open(fname, 'r')

    try:
        fuel = pickle.load(infile)
        print('RETRIEVED: this reduction step %s/%s '%(os.path.basename(fuelsave_dir),os.path.basename(fname)))
    except:
        print("Could not load %s"%(fname))

    return fuel
