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

log = logging.getLogger(__name__)
class DataBucket:

    """
    Like the Goodman pipeline's NightDataContainer class, this
    is a list of :class:`~pandas.DataFrame` objects.
    """

    def __init__(self,obsid,data_path):

        self.obsid = obsid
        self.data_path = data_path
        self.gain = None
        self.rdnoise = None
        self.ccdsum = None
        self.empty_bucket = True
        self.bias_buckets = None
        self.flat_buckets = None
        self.comp_buckets = None
        self.targ_buckets = None
        self.spec_buckets = None
        self.slit_buckets = None
        self.fits_slit_locs = None
        self.sunset_time = None
        self.sunrise_time = None
        self.evening_twilight = None
        self.morning_twilight = None
        self.header_keys = ['naxis','date','slit','date-obs','obstype','object',
                    'exptime','obsra','obsdec','grating','cam_targ','grt_targ',
                    'filter','filter2','gain','rdnoise','ccdsum','wavmode']


    def __repr__(self):
        """Produces a nice summary of the information contained"""
        if self.empty_bucket:
            return str("Empty Data Container")
        else:

            class_info = str("{:s}\n"
                             "SAMOS Observation ID: {:s}\n"
                             "Data Directory: {:s}\n"
                             "Working Directory: {:s}".format(str(self.__class__),
                                                      self.obsid,
                                                      os.path.split(self.data_path)[1],
                                                      os.path.split(os.getcwd())[1]))

            if all([self.gain, self.rdnoise, self.ccdsum]):
                class_info += str("\nGain: {:.2f}\n"
                                  "Readout Noise: {:.2f}\n"
                                  "ccdsum: {:s}".format(self.gain,
                                                     self.rdnoise,
                                                     self.ccdsum))
            class_info += str("\nIs Empty: {:s}\n".format(str(self.empty_bucket)))

            BUCKET_info = "\nBucket Grouping Information\n"
            BUCKET_info += "\n(small buckets inside bigger bucekts!)\n"
            BUCKET_info += "BIAS Group:\n"
            BUCKET_info += self._get_group_repr(self.bias_buckets)
            BUCKET_info += "FLATs Group:\n"
            BUCKET_info += self._get_group_repr(self.flat_buckets)
            BUCKET_info += "COMP Group:\n"
            BUCKET_info += self._get_group_repr(self.comp_buckets)
            BUCKET_info += "SPEC Group\n"
            BUCKET_info += self._get_group_repr(self.spec_buckets)

            class_info += BUCKET_info
            return class_info

    @staticmethod
    def _get_group_repr(BUCKET):
        """Converts the filenames in each bucket group to string
        This class has a __repr__ method and in this method the file names
        contained in the different groups gets formatted as a string for
        displaying in a readable way.
        """
        BUCKET_str = ""
        if BUCKET is not None:
            for i in range(len(BUCKET)):
                if len(BUCKET) == 1:
                    BUCKET_str += "Files in BUCKET\n"
                else:
                    BUCKET_str += "Files in BUCKET {:d}\n".format(i + 1)
                for _file in BUCKET[i]['file']:
                    BUCKET_str += "  {:s}\n".format(_file)
            return BUCKET_str
        else:

            return "  Group is Empty\n"

    def add_bias_bucket(self,bias_bucket):

        if self.bias_buckets is None:
            self.bias_buckets = [bias_bucket]
        else:
            self.bias_buckets.append(bias_bucket)

        if self.bias_buckets is not None:
            self.empty_bucket = False


    def add_flat_bucket(self,flat_bucket):

        if self.flat_buckets is None:
            self.flat_buckets = [flat_bucket]
        else:
            self.flat_buckets.append(flat_bucket)

        if self.flat_buckets is not None:
            self.empty_bucket = False
        print("added flat bucket %s"%(flat_bucket))

    def add_spec_bucket(self,spec_bucket):
        """add a data bucket with science frames and comparison
        lamps coupled together."""
        if self.spec_buckets is None:
            self.spec_buckets = [spec_bucket]
        else:
            self.spec_buckets.append(spec_bucket)
        if self.spec_buckets is not None:
            self.empty_bucket = False
        comp_bucket = spec_bucket[spec_bucket.obstype == 'COMP']
        self.add_comp_bucket(comp_bucket=comp_bucket)


    def add_comp_bucket(self,comp_bucket):

        if self.comp_buckets is None:
            self.comp_buckets = [comp_bucket]
        else:
            self.comp_buckets.append(comp_bucket)

        if self.comp_buckets is not None:
            self.empty_bucket = False

    def add_target_bucket(self,science_bucket):

        if self.targ_buckets is None:
            self.targ_buckets = [science_bucket]
        else:
            self.targ_buckets.append(science_bucket)

        if self.targ_buckets is not None:
            self.empty_bucket = False

    def add_combined_img_bucket(self,combined_img):
        """like the other buckets, combined_img is a
        list of pandas.DataFrame objects"""
        if self.combined_imgs_bucket is None:
            self.combined_imgs_bucket = [combined_img]
        else:
            self.combined_imgs_bucket.append(combined_img)

        if self.combined_imgs_bucket is not None:
            self.empty_bucket = False


    def set_sun_times(self, sun_set, sun_rise):
        """Sets values for sunset and sunrise
        Args:
            sun_set (str): Sun set time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
            sun_rise (str):Sun rise time in the format 'YYYY-MM-DDTHH:MM:SS.SS'
        """

        self.sun_set_time = sun_set
        self.sun_rise_time = sun_rise

    def set_twilight_times(self, evening, morning):
        """Sets values for evening and morning twilight
        Args:
            evening (str): Evening twilight time in the format
                'YYYY-MM-DDTHH:MM:SS.SS'
            morning (str): Morning twilight time in the format
                'YYYY-MM-DDTHH:MM:SS.SS'
        """

        self.evening_twilight = evening
        self.morning_twilight = morning

    def set_readout(self, gain, rdnoise, ccdsum):
        """Set Gain, Read noise and ccdsum.
        Args:
            gain (float): Gain from header
            rdnoise (float): Read noise from header.
            ccdsum (str): ccdsum from header.
        """
        self.gain = gain
        self.rdnoise = rdnoise
        self.ccdsum = ccdsum
