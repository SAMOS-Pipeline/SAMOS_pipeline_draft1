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



class BucketClassifier:
    """
    Based on Goodman Pipeline data_classifier.
    Classifies the data being presented to the pipeline.
    Data classifier is intended to define the camera that is being used and the
    technique in use. This will be used later to make important decisions
    regarding the process to be used.
    """

    def __init__(self):
        """Initialization method for the DataClassifier class
        The general arguments of the program are parsed and become part of the
        class attributes. The rest of attributes are initialized as None.
        """
        self.log = logging.getLogger(__name__)
        self.raw_path = None
        self.nights_dict = None
        self.instrument = None
        self.image_collection = None
        self.objects_collection = None

    def __repr__(self):
        """String representation of the information contained."""
        return str("Raw Path: {:s}".format(self.raw_path))

    def __call__(self, raw_path):
        """Call method for the DataClassifier class
        This method call specific method that define all the attributes of the
        class. The intention is to define the instrument and technique in use.
        Args:
            raw_path (str): Full Path to raw data
        """
        self.raw_path = raw_path

        # define the ImageFileCollection instance right away.

        try:
            ifc = ImageFileCollection(self.raw_path)

        except VerifyError as error:  # pragma: no cover
            self.log.error("Raised VerifyError: {:}".format(error))
            self.log.critical("Some keywords are not FITS compliant. Trying "
                              "to fix the headers.")

            fix_keywords(path=self.raw_path)

            self.log.info("Headers have been fixed, please rerun the pipeline!")
            sys.exit()

        self.image_collection = ifc.summary.to_pandas()

        self.objects_collection = self.image_collection[
            self.image_collection.obstype != 'BIAS']

        self.nights_dict = {}
        self.log.debug('Raw path: {:s}'.format(self.raw_path))

        self._get_instrument()
        if self.instrument is not None:
            self.log.info('Instrument: {:s} Camera'.format(self.instrument))
        else:
            self.log.critical("Unable to determine which camera was used.")
            self.log.info("Make sure you only have 'Blue' or 'Red' camera data "
                          "only, not both.")
            sys.exit()

        self._get_obs_technique()
        if self.technique is not None:
            self.log.info('Observing Technique: {:s}'.format(self.technique))
        # else:
        #     self.log.critical("Unable to determine observing technique used.")
        #     sys.exit()

        if self.instrument is not None and self.technique is not None:

            # folder name is used as key for the dictionary
            night = os.path.basename(self.raw_path)

            self.nights_dict[night] = {'full_path': self.raw_path,
                                       'instrument': self.instrument,
                                       'technique': self.technique}
        else:
            self.log.error('Failed to determine Instrument or Technique '
                           'for the night: {:s}'.format(self.raw_path))

    def _get_instrument(self):
        """Identify Goodman's Camera
        The header keyword of the camera is `INSTCONF`.
        Notes:
            This methods no longer offers backwards compatibility.
        """

        instconf = self.objects_collection.instconf.unique()

        if len(instconf) > 1:
            for _inst in instconf:
                self.log.debug("INSTCONF = {:s} is present.".format(_inst))
            self.log.warning("Camera changes are forbidden during the night")
        elif len(instconf) == 1:
            self.instrument = instconf[0]
            self.log.debug("Detected {:s} camera.".format(self.instrument))
        # else:
        #     self.log.error("Impossible to determine which camera was used.")

    def _get_obs_technique(self):
        """Identify if the data is Imaging or Spectroscopy
        For imaging data the keyword `WAVMODE` is `Imaging` therefore the logic
        here is: If there is only one value for `WAVMODE` and it is `Imaging`
        then the technique is `Imaging`. If `Imaging` is in the result along
        with other then it will assume the technique is Spectroscopy and will
        ignore all the Imaging data. If none of the conditions above are met it
        will assume the technique is Spectroscopy.
        The result is stored as an attribute of the class.
        """

        wavmodes = self.objects_collection.wavmode.unique()

        if len(wavmodes) == 1 and wavmodes[0] == 'Imaging':
                self.technique = 'Imaging'

        elif 'Imaging' in wavmodes and len(wavmodes) > 1:
                self.log.error('There seems to be Imaging and Spectroscopic '
                               'data. I will assume the Imaging data are '
                               'acquisition images therefore they will be '
                               'ignored.')
                self.log.info("If you really have Imaging data, please process "
                              "them in a separated folder.")
                self.technique = 'Spectroscopy'
        else:
                self.technique = 'Spectroscopy'
        # inform the results, no need to return
        self.log.info('Detected {:s} Data from {:s} '
                      'Camera'.format(self.technique, self.instrument))


if __name__ == '__main__':
    passnique = None


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
