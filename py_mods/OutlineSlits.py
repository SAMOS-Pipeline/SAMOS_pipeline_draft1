#!/usr/bin/env python
from __future__ import print_function
import glob
from astropy.io import fits
from FlatNorm import *
from SAMOSHelpers import *
from SlitID import get_edges
from SlitCutout import cutout_slit
from CreateSlitStructure import CreateSlitStructure
from SaveFuel import save_fuel_step

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
            if slit_cutout_paths[slit] in targs_cutout_paths:
                continue
            targs_cutout_paths.append(slit_cutout_paths[slit])

    targs_cutout_paths = np.asarray(sorted(targs_cutout_paths,key=lambda x: \
                        int(os.path.split(x)[1].split('_')[0][4:])))

    fuel.util.science.corr_files = targs_cutout_paths


    flat_cutout_paths = []
    slit_cutout_paths = cutout_slit(masterFlat,x_edges_main,y_edges_main)
    for slit in range(len(slit_cutout_paths)):
        flat_cutout_paths.append(slit_cutout_paths[slit])

    flat_cutout_paths = np.asarray(flat_cutout_paths)

    fuel.util.slitflat.corr_files = flat_cutout_paths


    arcs_cutout_paths = []
    for arc in range(len(inlamps)):
        slit_cutout_paths = cutout_slit(inlamps[arc],x_edges_main,y_edges_main)
        for slit in range(len(slit_cutout_paths)):
            if slit_cutout_paths[slit] in arcs_cutout_paths:
                continue
            arcs_cutout_paths.append(slit_cutout_paths[slit])

    arcs_cutout_paths = np.asarray(sorted(arcs_cutout_paths,key=lambda x: \
                        int(os.path.split(x)[1].split('_')[0][4:])))

    fuel.util.arc.corr_files = arcs_cutout_paths

    update_slits(fuel,fuel.util.science.corr_files,fuel.input.mask_SMF,slit_positions)
    print('slits have been identified and excised. Ready for PyRAF identify.')

    save_fuel_step(fuel,'OutlineSlits')

    return fuel

def update_slits(fuel,cutout_targs,maskf,slit_pos):


    mask_file = open(maskf, 'r')
    y_pos_pix = np.genfromtxt(slit_pos)
    slit_num = 1
    nums = []
    xs,ys,objs,ras,decs = [],[],[],[],[]
    for line in mask_file:
        line = line.split()
        if "SLIT" in line:

            number = slit_num
            obj = line[1]
            obj_ra,obj_dec = line[2],line[3]
            x_mm,y_mm = float(line[4]),float(line[5])
            alen_mm,blen_mm = float(line[7]),float(line[8])
            width_mm = float(line[6])
            angle = float(line[9])

            xs.append(x_mm)
            ys.append(y_mm)
            objs.append(obj)
            ras.append(obj_ra)
            decs.append(obj_dec)
    xs,ys,objs,ras,decs = zip(*sorted(zip(xs,ys,objs,ras,decs),reverse=True))
    fuel_updated_slits = []
    for slit in range(int(len(cutout_targs)/2.)):
        this_slit = CreateSlitStructure().create_slit_structure(maskf)
        #this_slit = fuel.slits[slit]
        this_slit.number = slit_num
        this_slit.obj = objs[slit]
        this_slit.obj_ra = ras[slit]
        this_slit.obj_dec = decs[slit]
        this_slit.x_mm = xs[slit]
        this_slit.y_mm = ys[slit]
        fuel_updated_slits.append(this_slit)
        slit_num +=1

    fuel.slits = fuel_updated_slits

    return
