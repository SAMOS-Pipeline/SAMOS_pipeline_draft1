from __future__ import print_function
import sys
import os
import json
sys.path.insert(0,os.path.split(os.getcwd())[0])
print('Main working directory is: %s'%(sys.path[0]))
from SAMOSHelpers import *
from astropy.io import fits
import glob

class CreateInput:

    def __init__(self):
        self.working_dir = os.path.split(os.getcwd())[0]
        self.db_file = ''
        self.science_filelist = ''
        self.arc_filelist = ''
        self.flat_filelist = ''
        self.field_filelist = ''
        self.slit_position_file = ''
        self.slit_mask = ''
        self.mask_SMF = ''


    def create_input(self,datedir,mask):

        #working_dir = os.path.split(os.getcwd())[0]

        instdir  = "%s/LDSS3"%(self.working_dir)

        datadir = "%s/%s" % (instdir,datedir)

        IsItHere(datadir)

        print("Scanning %s\n" %(datadir))
        images = glob.glob("%s/ccd????c?.fits"%(datadir))
        images.sort()

        fitslist = [fits.open(image) for image in images]
        fitslist, images = zip(*[(f,i) for f,i in zip(fitslist,images) if f[0].header["aperture"]==mask]) # IS THE MASK IN?

        fieldlist,fieldimages = zip(*[(f,i) for f,i in zip(fitslist,images) if f[0].header["grism"]=="Open"])
        fieldlist,fieldimages = fieldlist[-2:],fieldimages[-2:]
        print("Found field image:")
        print(list(map(os.path.basename,fieldimages)))
        fitslist, images = zip(*[(f,i) for f,i in zip(fitslist,images) if f[0].header["grism"]!="Open"])  # IS A DISPERSER IN?


        flatlist, flats = zip(*[(f,i) for f,i in zip(fitslist,images) if "flat" in f[0].header["object"].lower()])
        print("Found flats:")
        print(list(map(os.path.basename,flats)))

        lamplist, lamps = zip(*[(f,i) for f,i in zip(fitslist,images) if len([s for s in ["he","ne","ar"] if s in f[0].header["object"].lower()])>=2])
        print("Found lamps:")
        print(list(map(os.path.basename,lamps)))

        targlist, targs = zip(*[(f,i) for f,i in zip(fitslist,images) if i not in flats and i not in lamps])
        print("Found science:")
        print(list(map(os.path.basename,targs)))

        print("With exposure times:")
        exptimes = np.array([f[0].header["exptime"] for f in targlist])
        print(exptimes)

        targlist, targs = zip(*[(f,i) for f,i,t in zip(fitslist,images,exptimes) if t>0.5*exptimes.max()])
        print("Keeping science:")
        print(list(map(os.path.basename,targs)))


        self.science_filelist = targs
        self.arc_filelist = lamps
        self.flat_filelist = flats
        self.field_filelist = fieldimages
        self.slit_position_file = '%s/helper_files/%s_ycoords.txt'%(self.working_dir,mask)
        self.slit_mask = mask
        self.mask_SMF = '%s/helper_files/%s.SMF'%(self.working_dir,mask)

        db = open("%s/%s.db"%(self.working_dir,mask),"w")
        db.write("# lamp\n")
        [db.write("%s\n"%(lamp)) for lamp in self.arc_filelist]
        db.write("# flat\n")
        [db.write("%s\n"%(flat)) for flat in self.flat_filelist]
        db.write("# targ\n")
        [db.write("%s\n"%(targ)) for targ in self.science_filelist]
        db.write("# field\n")
        [db.write("%s\n"%(field)) for field in self.field_filelist]
        db.close()

        if not os.path.exists('%s/%s'%(self.working_dir,self.slit_mask)): os.mkdir('%s/%s'%(self.working_dir,self.slit_mask))

        self.slit_mask = '%s/%s'%(self.working_dir,self.slit_mask)
        self.db_file = "%s/%s.db"%(self.working_dir,mask)
        return self
