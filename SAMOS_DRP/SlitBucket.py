from __future__ import print_function
import sys
import os
import json
from astropy.io import fits
from astropy import units
import glob
import pandas as pd
import ccdproc
from ccdproc import CCDData
import logging
import argparse
from astropy.io.fits.verify import VerifyError
from ccdproc import ImageFileCollection
from drp_mods import cutout_slit,read_fits,write_fits


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
        self.out_prefix_list = None
        #each new slit ccd file will have sNNN_prefix_filename.fits, where N
        #is the slit number
    def __call__(self):


        self.log.debug("Working on slit cutouts.")
        image_list = self.spec_buckets[0]['file'].tolist()
        sample_image = os.path.join(self.raw_data_dir,
                                    np.random.choice(image_list))
        sample_ccd = read_fits(sample_image)
        self.slit_mask =  sample_ccd.header['SLIT']
        self.mask_proc_dir = os.path.join(self.processing_dir,"{:s}".format(self.slit_mask))
        if not os.path.exists(self.mask_proc_dir):
            log.debug("making slitmask directory {:s}".format(self.mask_proc_dir))
            os.makedir(self.mask_proc_dir)

        for fname in image_list:

            ccdpath = os.path.join(self.mask_proc_dir,"{:s}{:s}.fits".format(self.in_prefix,fname))
            ccd = read_fits(ccdpath)
            ccd_slits,trimsecs = cutout_slit(ccd,self.slit_edges_x,self.slit_edges_y)
            slit_num = 1
            for slit in ccd_slits:
                out_prefix = "s{:s}_".format(str(slit_num))+self.in_prefix
                slit_ccd = CCDData(slit,unit='adu')
                new_slit_fname = os.path.join(self.mask_proc_dir,"{:s}{:s}.fits".format(out_prefix,fname))

                slit_ccd.write(new_slit_fname,overwrite=True)
                slit_num+=1
                if self.out_prefix_list is None:
                    self.out_prefix_list = [out_prefix]
                else:
                    self.out_prefix_list.append(out_prefix)



    #def add_slit_cutout(self):
