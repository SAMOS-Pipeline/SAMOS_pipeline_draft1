import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
import struct

def create_fuel_structure(input):

    '''
    Creating a data structure for storage of information about the observation to be
    referenced easily by the pipeline.  Basically trying to translate the Flame DRP to python.
    '''

    print("input is: %s"%(dir(input)))
