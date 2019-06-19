from __future__ import print_function
import glob
from SAMOSHelpers import *
from Overscan import *
from SaveFuel import save_fuel_step

def overscan_and_trim(fuel):

    inlamps = fuel.util.arc.raw_files
    inflats = fuel.util.slitflat.raw_files
    intargs = fuel.util.science.raw_files
    infield = fuel.util.field.raw_files

    if not inlamps or not inflats or not intargs:
        MuyMalo("Problem with what files are available!")

    db = fuel.input.db_file #sys.argv[1]
    IsItHere(db)
    mask = os.path.split(fuel.input.slit_mask)[1] #db[:-3]

    corr_dir = fuel.util.intermediate_dir



    outlamps = ["%s/b%s%s.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in inlamps]
    outflats = ["%s/b%s%s.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in inflats]
    outtargs = ["%s/b%s%s.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in intargs]
    outfield = ["%s/b%sfield%s.fits" % (corr_dir,mask,os.path.basename(i)[7:-5]) for i in infield]

    # IF YOU THINK YOU MIGHT SPAWN JOBS IN THE BACKGROUND OR WHATEVER
    #infiles = inlamps+inflats+intargs
    #outfiles = outlamps+outflats+outtargs
    #commands = ["python Overscan.py %s %s" % (i,o) for i,o in zip(infiles,outfiles)]
    #list(map(os.system,commands))

   # IF YOU LIKE YOUR SCRIPTS TO GO AND DO EVERYTHING
    [Overscan(i,o) for i,o in zip(inlamps,outlamps)]
    [Overscan(i,o) for i,o in zip(inflats,outflats)]
    [Overscan(i,o) for i,o in zip(intargs,outtargs)]
    [FieldTrim(i,o) for i,o in zip(infield,outfield)]


    stitchlamps = [stitch(o) for o in outlamps]
    stitchflats = [stitch(o) for o in outflats]
    stitchtargs = [stitch(o) for o in outtargs]
    stitchfield = [stitch(o) for o in outfield]

    print(stitchlamps)

    fuel.util.arc.corr_files = list(set(stitchlamps))
    fuel.util.science.corr_files = list(set(stitchtargs))
    fuel.util.slitflat.corr_files = list(set(stitchflats))
    fuel.util.field.corr_files = list(set(stitchfield))



    db = open('%s'%fuel.input.db_file,"a")
    db.write("\n## Trimmed/Overscan Subtracted\n")
    db.write("# OSLamps\n")
    [db.write("%s\n"%(olamp)) for olamp in outlamps]
    db.write("# OSFlats\n")
    [db.write("%s\n"%(oflat)) for oflat in outflats]
    db.write("# OSTargs\n")
    [db.write("%s\n"%(otarg)) for otarg in outtargs]
    db.close()



    save_fuel_step(fuel,'OverscanTrim')

    return fuel
