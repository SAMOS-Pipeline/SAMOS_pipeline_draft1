#!/usr/bin/env python
from __future__ import print_function
import glob
from astropy.io import fits
from FlatNorm import *
from SAMOSHelpers import *
from StartFrom import save_fuel_step


def norm_div_flats(fuel):

    inlamps = fuel.util.arc.corr_files
    inflats = fuel.util.slitflat.corr_files
    intargs = fuel.util.science.corr_files

    if not inlamps or not inflats or not intargs:
        MuyMalo("Problem with what files are available!")

    db = fuel.input.db_file #sys.argv[1]
    print(db)
    IsItHere(db)
    mask = os.path.split(fuel.input.slit_mask)[1] #db[:-3]
    corr_dir = fuel.util.intermediate_dir

    [IsItHere(i) for i in inflats]

    print("Created list of flats in %s\n" %(mask))

    master_flat = "%s/fb%s%s" %(corr_dir,mask,'master_flat.fits')

    fuel.util.slitflat.master_file = master_flat

    [MkFlatNorm(inflats,master_flat)]


    targflat_fielded_path = ["%s/%s" % (corr_dir,'f'+os.path.basename(i)) for i in intargs]
    lampflat_fielded_path = ["%s/%s" % (corr_dir,'f'+os.path.basename(i)) for i in inlamps]

    if not os.path.exists(os.path.split(targflat_fielded_path[0])[0]):  os.mkdir(os.path.split(targflat_fielded_path[0])[0])


    [NormFlatDiv(i,master_flat,j) for i,j in zip(intargs,targflat_fielded_path)]
    [NormFlatDiv(i,master_flat,j) for i,j in zip(inlamps,lampflat_fielded_path)]

    fuel.util.science.corr_files = targflat_fielded_path
    fuel.util.arc.corr_files = lampflat_fielded_path


    db = open("%s"%(fuel.input.db_file),"a")
    db.write("\n")
    db.write("## Flat Fielded\n")
    db.write("# FFLamps\n")
    [db.write("%s\n"%(olamp)) for olamp in lampflat_fielded_path]
    db.write("# FFTargs\n")
    [db.write("%s\n"%(otarg)) for otarg in targflat_fielded_path]
    db.write("# MasterFlat\n")
    db.write("%s" %(master_flat))

    db.close()

    save_fuel_step(fuel,'NormFlatDiv')

    return fuel
