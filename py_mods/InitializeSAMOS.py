import sys
import os
import json
from SAMOSHelpers import *
from astropy.io import fits
import glob
from CreateFuelStructure import CreateFuelStructure
from CreateInput import CreateInput


def initialize_SAMOS(datedir,mask):

    input_structure = CreateInput().create_input(datedir,mask)

    fuel = CreateFuelStructure().create_fuel_structure(input_structure)

    return fuel
