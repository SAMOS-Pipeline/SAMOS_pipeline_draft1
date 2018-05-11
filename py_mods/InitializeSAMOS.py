import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
from CreateFuelStructure import create_fuel_structure


def initialize_SAMOS(input):

    fuel = create_fuel_structure(input)
