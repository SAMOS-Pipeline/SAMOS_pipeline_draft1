from __future__ import print_function
import sys
import os
import json
from astropy.io import fits
from astropy import units
import glob
import pandas as pd
import numpy as np
import ccdproc
from ccdproc import CCDData
import logging
import argparse
from astropy.io.fits.verify import VerifyError
from ccdproc import ImageFileCollection
from .drp_mods import cutout_slit,read_fits,write_fits


class SlitBuckets:
    """
    This class is to organize the slit images from a
    comp/targ bucket and pair each slit with its
    respective comparison lamp and object image
    __init__ arg takes an image bucket that has
    already been processed by ImageProcessor
    """

    def __init__(self,imgproc_night_bucket):

        self.log = logging.getLogger(__name__)
        self.raw_data_dir = imgproc_night_bucket.raw_data_dir
        self.processing_dir = imgproc_night_bucket.processing_dir
        self.mask_proc_dir = None
        self.gain = imgproc_night_bucket.gain
        self.rdnoise = imgproc_night_bucket.rdnoise
        self.ccdsum = imgproc_night_bucket.ccdsum
        self.empty_bucket = True
        self.slit_num = None
        self.slit_targs = None
        self.slit_comps = None
        self.slit_mask = None
        self.slit_edges_x = imgproc_night_bucket.slit_edges_x
        self.slit_edges_y = imgproc_night_bucket.slit_edges_y
        self.spec_buckets = imgproc_night_bucket.spec_buckets
        self.in_prefix = imgproc_night_bucket.out_prefix
        self.combined_bucket = imgproc_night_bucket.combined_bucket
        self.out_prefix_list = None
        #each new slit ccd file will have sNNN_prefix_filename.fits, where N
        #is the slit number
    def __call__(self):


        self.log.debug("Working on slit cutouts.")

        if self.combined_bucket is not None:
            #combined_bucket is a single dataframe, not a list of dataframes
            image_list = self.combined_bucket['file'].tolist()
            sample_img = os.path.join(self.processing_dir,
                                    np.random.choice(image_list))
        else:
            image_list = self.spec_buckets[0]['file'].tolist()

            sample_img = os.path.join(self.raw_data_dir,
                                    np.random.choice(image_list))

        sample_ccd = read_fits(sample_img)
        self.slit_mask =  sample_ccd.header['SLIT']
        self.mask_proc_dir = os.path.join(self.processing_dir,"{:s}".format(self.slit_mask))
        if not os.path.exists(self.mask_proc_dir):
            self.log.debug("making slitmask directory {:s}".format(self.mask_proc_dir))
            os.mkdir(self.mask_proc_dir)
        elif len(os.listdir(self.mask_proc_dir))>1:
            self.log.debug("cleaning slitmask directory of previous files {:s}".format(self.mask_proc_dir))
            [os.remove(os.path.join(self.mask_proc_dir,i)) for i in os.listdir(self.mask_proc_dir)]

        for fname in image_list:
            if self.combined_bucket is not None:
                ccdpath = os.path.join(self.processing_dir,fname)

            else:
                ccdpath = os.path.join(self.processing_dir,"{:s}{:s}".format(self.in_prefix,fname))

            ccd = read_fits(ccdpath)

            ccd_slits,trimsecs = cutout_slit(ccd,self.slit_edges_x,self.slit_edges_y)
            self.slit_num = 1
            for nslit in range(len(ccd_slits)):
                slit = ccd_slits[nslit]
                if self.combined_bucket is not None:
                    out_prefix = "s{:s}_".format(str(self.slit_num))
                else:
                    out_prefix = "s{:s}_".format(str(self.slit_num))+out_prefix
                slit_ccd = CCDData(slit,unit='adu')
                new_slit_fname = os.path.join(self.mask_proc_dir,"{:s}{:s}".format(out_prefix,fname))

                slit_ccd.header = ccd.header.copy()
                slit_ccd.header.set('SLITNUM',
                                    value="{:s}".format(str(self.slit_num).zfill(2)),
                                    comment='')
                slit_ccd.write(new_slit_fname,overwrite=True)
                self.slit_num+=1

                if self.out_prefix_list is None:
                    self.out_prefix_list = [out_prefix]
                else:
                    self.out_prefix_list.append(out_prefix)

                if slit_ccd.header['OBSTYPE']=='OBJECT':
                    self.log.debug("appending slit target")
                    self.add_slit_targ(new_slit_fname)

                elif slit_ccd.header['OBSTYPE']=='COMP':
                    self.log.debug("appending slit comparison lamp")
                    self.add_slit_comp(new_slit_fname)



    def add_slit_comp(self,slit_comp):

        if self.slit_comps is None:
            self.slit_comps = [slit_comp]
        else:
            self.slit_comps.append(slit_comp)

    def add_slit_targ(self,slit_targ):

        if self.slit_targs is None:
            self.slit_targs = [slit_targ]
        else:
            self.slit_targs.append(slit_targ)
