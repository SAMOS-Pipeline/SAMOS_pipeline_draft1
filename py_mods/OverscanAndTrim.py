import glob
from SAMOSHelpers import *
from Overscan import *


def overscan_and_trim(fuel):

    inlamps = fuel.util.arc.raw_files
    inflats = fuel.util.slitflat.raw_files
    intargs = fuel.util.science.raw_files
    infield = fuel.input.field_filelist

    if not inlamps or not inflats or not intargs:
        MuyMalo("Problem with what files are available!")

    db = fuel.input.db_file #sys.argv[1]
    IsItHere(db)
    mask = os.path.split(fuel.input.slit_mask)[1] #db[:-3]
    corr_dir = fuel.util.intermediate_dir


    outlamps = ["%s/%s%s_corr.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in inlamps]
    outflats = ["%s/%s%s_corr.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in inflats]
    outtargs = ["%s/%s%s_corr.fits" % (corr_dir,mask,os.path.basename(i)[3:-5]) for i in intargs]
    outfield = ["%s/%sfield%s_corr.fits" % (corr_dir,mask,os.path.basename(i)[7:-5]) for i in infield]

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

    fuel.util.arc.corr_files = outlamps
    fuel.util.science.corr_files = outtargs
    fuel.util.slitflat.corr_files = outflats
    fuel.util.field.corr_files = outfield



    db = open('%s'%fuel.input.db_file,"a")
    db.write("\n## Trimmed/Overscan Subtracted\n")
    db.write("# OSLamps\n")
    [db.write("%s\n"%(olamp)) for olamp in outlamps]
    db.write("# OSFlats\n")
    [db.write("%s\n"%(oflat)) for oflat in outflats]
    db.write("# OSTargs\n")
    [db.write("%s\n"%(otarg)) for otarg in outtargs]
    db.close()

    return fuel
