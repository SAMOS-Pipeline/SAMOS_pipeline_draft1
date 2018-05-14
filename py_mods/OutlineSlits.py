#!/usr/bin/env python

import glob
from astropy.io import fits
from FlatNorm import *
from SAMOSHelpers import *
from SlitID import *
from SlitCutout import *


def outline_slits(fuel):
    inlamps = fuel.util.arc.corr_files
    inflats = fuel.util.slitflat.corr_files
    intargs = fuel.util.science.corr_files

    if not inlamps or not inflats or not intargs:
        MuyMalo("Problem with what files are available!")

    db = fuel.input.db_file #sys.argv[1]
    IsItHere(db)
    mask = os.path.split(fuel.input.slit_mask)[1] #db[:-3]
    corr_dir = fuel.util.intermediate_dir

    masterFlat = fuel.util.slitflat.master_file
    slit_positions = fuel.input.slit_position_file

    #I will want to import a file with approximate slit top edges to read into get_edges().

    x_edges_main,y_edges_main = get_edges(masterFlat,slit_positions,fuel.input.slit_mask)

    print("Working on cutting out slits....")

    cutout_slit_arrs = cutout_slit(intargs[0],x_edges_main,y_edges_main)
