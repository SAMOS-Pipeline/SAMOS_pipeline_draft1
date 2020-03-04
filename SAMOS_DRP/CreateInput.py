from __future__ import print_function
import sys
import os
import json
#print('Main working directory is: %s'%(sys.path[0]))
from SAMOSHelpers import check_header_notes
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as units
import glob
import pandas as pd
import numpy as np
import ccdproc
from ccdproc import CCDData, ImageFileCollection
from DataBucket import DataBucket


class CreateInput(object):

    """
    Based on goodman_pipeline `classify_spectroscopic_data`,
    which uses :class: `~ccdproc.ImageFileCollection`.
    This function is meant to group the spectroscopic observations
    based on the grism used.
    There will be five grisms (2 lo-res, 3 hi-res) available with SAMOS,
    but only 4 can be mounted at a given time.
    Two low-res: 400-600nm blue, 702.5 l/mm, R = 3137
                 600-950nm red, 416.9 l/mm, R = 2791
    Three high-res: H-alpha, 1250 l/mm, R = 8111
                    H-beta, 1775 l/mm, R=9601
                    *CaIII region, 1198.5 l/mm, R=10141
    *not in the default grism configuration so it needs to be specified in
    the instrument setup report for the telescope operator.
    (more info on spectroscopic setup in the SAMOS opto-mechanical paper!)

    Current version of this is for test data from Goodman.

    FITS header key words used (from the Goodman pipeline):

        - ``slit``
        - ``radeg``
        - ``decdeg``
        - ``grating``
        - ``cam_targ``
        - ``grt_targ``
        - ``filter``
        - ``filter2``
        - ``gain``
        - ``rdnoise``
    """

    def __init__(self,obsid,workdir,rawdir):
        """
        obsid corrosponds to the night on which the observations were taken
        this obsid will eventually point to a table in a SQL database that
        will have information about the telescope settings that night
        """

        self.working_dir = workdir
        self.obsid = obsid
        self.header_keys = ['naxis','date','slit','date-obs','obstype','object',
                    'exptime','obsra','obsdec','grating','cam_targ','grt_targ',
                    'filter','filter2','gain','rdnoise','roi','wavmode']


        print("Scanning %s\n" %(rawdir))
        rawimages = glob.glob("%s/psg_1403??_??????_???.fits.fz"%(rawdir))
        rawimages.sort()
        #the ccdproc class ImageFileCollection cannot read compressed files
        #so here is where I make a new raw data dir with fits images for the
        #main data
        uncomp_dir = '%s/DATA_RAW/UNCOMP'%(self.working_dir)
        if not os.path.exists(uncomp_dir):
            os.mkdir(uncomp_dir)

        images = glob.glob('%s/psg_1403??_??????_???.fits'%(uncomp_dir))
        if len(images)==0:
            for raw in rawimages:
                #test files are compressed, which can't be read by ImageFileCollection so
                #I have to save new versions of the fits files with only the relevant hdu
                uncomp_img_name = os.path.split(raw)[1][:-3]
                uncomp_path = '%s/%s'%(uncomp_dir,uncomp_img_name)
                ccd = CCDData.read(raw,unit=units.adu,hdu=1,hdu_mask=[1,0])
                copy = ccd.copy()

                copy.data = copy.data[0]
                copy.write(uncomp_path,overwrite=True)
                uncomp_imgs.append(uncomp_img_name)

        #images = ['%s'%(uname) for uname in uncomp_imgs]
        fitslist = []
        for fits in images:
            fitsccd = CCDData.read(fits)
            fitslist.append(fitsccd)

        gratings = [f.header["GRATING"] for f in fitslist]
        print("Gratings used:")
        print(set(list(gratings)))

        flatlist,flats = zip(*[(f,i) for f,i in zip(fitslist,images) if f.header["OBSTYPE"].strip(' ')=="FLAT"])
        print("Found flats:")
        print(list(map(os.path.basename,flats)))

        biaslist,biases = zip(*[(f,i) for f,i in zip(fitslist,images) if f.header["OBSTYPE"].strip(' ')=="BIAS"])
        print("Found bias frames:")
        print(list(map(os.path.basename,biases)))

        complist,comps = zip(*[(f,i) for f,i in zip(fitslist,images) if f.header["OBSTYPE"].strip(' ')=="COMP"])
        print("Found comparison lamps:")
        print(list(map(os.path.basename,comps)))

        targlist,targs = zip(*[(f,i) for f,i in zip(fitslist,images) if np.logical_and(f.header["OBSTYPE"].strip(' ')=="OBJECT",
                                                                                    check_header_notes(f.header)==False)])
        print("Found science:")
        print(list(map(os.path.basename,targs)))

        fieldlist,fields = zip(*[(f,i) for f,i in zip(fitslist,images) if np.logical_and(f.header["OBSTYPE"].strip(' ')=="OBJECT",
                                                                                      "field" in f.header["OBJECT"])])
        print("Found field images:")
        print(list(map(os.path.basename,fields)))

        self.flat_filelist = flats
        self.comp_filelist = comps
        self.science_filelist = targs
        self.field_filelist = fields
        self.bias_filelist = biases
        self.grating_list = gratings
        self.raw_data_dir = uncomp_dir
        self.full_data_bucket = np.concatenate(([flats],[comps],[targs],[biases]))

        self.db_file = "%s/%s.db"%(self.working_dir,obsid)
        db = open(self.db_file,"w")
        db.write("## Initialization\n")
        db.write("# lamp\n")
        [db.write("%s\n"%(lamp)) for lamp in self.comp_filelist]
        db.write("# flat\n")
        [db.write("%s\n"%(flat)) for flat in self.flat_filelist]
        db.write("# targ\n")
        [db.write("%s\n"%(targ)) for targ in self.science_filelist]
        db.write("# field\n")
        [db.write("%s\n"%(field)) for field in self.field_filelist]
        db.write("# bias\n")
        [db.write("%s\n"%(bias)) for bias in self.field_filelist]
        db.close()


    def group_spec_obs(self,images,keywords):

        """
        Based on goodman_pipeline `classify_spectroscopic_data`,
        which uses :class: `~ccdproc.ImageFileCollection`.
        This function is meant to group the spectroscopic observations
        based on the grism used.
        There will be five grisms (2 lo-res, 3 hi-res) available with SAMOS,
        but only 4 can be mounted at a given time.
        Two low-res: 400-600nm blue, 702.5 l/mm, R = 3137
                     600-950nm red, 416.9 l/mm, R = 2791
        Three high-res: H-alpha, 1250 l/mm, R = 8111
                        H-beta, 1775 l/mm, R=9601
                        *CaIII region, 1198.5 l/mm, R=10141
        *not in the default grism configuration so it needs to be specified in
        the instrument setup report for the telescope operator.
        (more info on spectroscopic setup in the SAMOS opto-mechanical paper!)

        Current version of this is for test data from Goodman.

        FITS header key words used (from the Goodman pipeline):

            - ``slit``
            - ``radeg``
            - ``decdeg``
            - ``grating``
            - ``cam_targ``
            - ``grt_targ``
            - ``filter``
            - ``filter2``
            - ``gain``
            - ``rdnoise``
        """

        grp_ifc = ccdproc.ImageFileCollection(filenames=images,keywords=keywords)

        self.org_fits = grp_ifc.summary.to_pandas() #to a pandas data frame!

        self.org_fits['radeg'] = ''
        self.org_fits['decdeg'] = ''
        #convert hh:mm:ss,dd:mm:ss to degrees.
        for i in self.org_fits.index.tolist():
            SC = SkyCoord(self.org_fits.obsra.iloc[i], self.org_fits.obsdec.iloc[i],unit=(units.hourangle,units.deg))
            radeg,decdeg = SC.ra,SC.dec

            self.org_fits.iloc[i, self.org_fits.columns.get_loc('radeg')] = '{:.2f}'.format(radeg)

            self.org_fits.iloc[i, self.org_fits.columns.get_loc('decdeg')] = '{:.2f}'.format(decdeg)

        #organize files by the instrument configuration and targets defined from the headers
        readout_configs = self.org_fits.groupby(['roi','gain','rdnoise']).size().reset_index().rename(columns={0: 'count'})

        sub_buckets = []

        for i in readout_configs.index:
            spec_group = DF_grp_ifc[((DF_grp_ifc['slit'] == readout_configs.iloc[i]['slit']) &
                               (DF_grp_ifc['radeg'] == readout_configs.iloc[i]['radeg']) &
                               (DF_grp_ifc['decdeg'] == readout_configs.iloc[i]['decdeg']) &
                               (DF_grp_ifc['grating'] == readout_configs.iloc[i]['grating']) &
                               (DF_grp_ifc['cam_targ'] == readout_configs.iloc[i]['cam_targ']) &
                               (DF_grp_ifc['grt_targ'] == readout_configs.iloc[i]['grt_targ']) &
                               (DF_grp_ifc['filter'] == readout_configs.iloc[i]['filter']) &
                               (DF_grp_ifc['filter2'] == readout_configs.iloc[i]['filter2']) &
                               (DF_grp_ifc['gain'] == readout_configs.iloc[i]['gain']) &
                               (DF_grp_ifc['rdnoise'] == readout_configs.iloc[i]['rdnoise']))]

            group_obstype = spec_group.obstype.unique()
            data_bucket.config_bucket(roi=readout_configs.iloc[i]['roi'],
                                      gain=readout_configs.iloc[i]['gain'],
                                      rdnoise=readout_configs.iloc[i]['rdnoise'])

            if 'COMP' in group_obstype and len(group_obstype) == 1:
                data_bucket.add_comp(spec_group)

            if 'OBJECT' in group_obstype and len(group_obstype) == 1:
                data_bucket.add_science(spec_group)


            sub_buckets.append(data_bucket)

        self.buckets = sub_buckets
        return self
