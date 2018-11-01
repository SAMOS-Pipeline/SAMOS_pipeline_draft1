#!/usr/bin/env python
from __future__ import print_function
import glob
from astropy.io import fits
from FlatNorm import *
from SAMOSHelpers import *
from SlitID import get_edges
from SlitCutout import cutout_slit


def outline_slits(fuel):
    inlamps = fuel.util.arc.corr_files
    masterFlat = fuel.util.slitflat.master_file
    intargs = fuel.util.science.corr_files

    if not inlamps or not masterFlat or not intargs:
        MuyMalo("Problem with what files are available!")

    db = fuel.input.db_file #sys.argv[1]
    IsItHere(db)
    mask = os.path.split(fuel.input.slit_mask)[1] #db[:-3]
    corr_dir = fuel.util.intermediate_dir
    slit_positions = fuel.input.slit_position_file

    #I will want to import a file with approximate slit top edges to read into get_edges().

    #obj_id,x_edges_main,y_edges_main = get_edges(masterFlat,slit_positions,fuel.input.slit_mask)
    x_edges_main,y_edges_main = get_edges(masterFlat,slit_positions,fuel.input.slit_mask)

    print("Working on cutting out slits....")

    targs_cutout_paths = []
    for targ in range(len(intargs)):
        slit_cutout_paths = cutout_slit(intargs[targ],x_edges_main,y_edges_main)
        for slit in range(len(slit_cutout_paths)):
            targs_cutout_paths.append(slit_cutout_paths[slit])

    targs_cutout_paths = np.asarray(targs_cutout_paths)

    fuel.util.science.corr_files = targs_cutout_paths

    flat_cutout_paths = []
    slit_cutout_paths = cutout_slit(masterFlat,x_edges_main,y_edges_main)
    for slit in range(len(slit_cutout_paths)):
        flat_cutout_paths.append(slit_cutout_paths[slit])

    flat_cutout_paths = np.asarray(flat_cutout_paths)

    fuel.util.slitflat.corr_files = flat_cutout_paths

    """
    for flat in range(len(inflats)):
        slit_cutout_paths = cutout_slit(inflats[flat],x_edges_main,y_edges_main)
        for slit in range(len(slit_cutout_paths)):
            flats_cutout_paths.append(slit_cutout_paths[slit])

    flats_cutout_paths = np.asarray(flats_cutout_paths)

    fuel.util.slitflat.corr_files = flats_cutout_paths
    """
    arcs_cutout_paths = []
    for arc in range(len(inlamps)):
        slit_cutout_paths = cutout_slit(inlamps[arc],x_edges_main,y_edges_main)
        for slit in range(len(slit_cutout_paths)):
            arcs_cutout_paths.append(slit_cutout_paths[slit])

    arcs_cutout_paths = np.asarray(arcs_cutout_paths)

    fuel.util.arc.corr_files = arcs_cutout_paths

    return fuel
