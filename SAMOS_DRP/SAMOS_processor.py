from __future__ import print_function
import ccdproc
import datetime
import glob
import logging
import numpy as np
import matplotlib.pyplot as plt
import random
import re
import os

from astropy import units as u
from astropy.io import fits
from ccdproc import CCDData

class ImageProcessor(object):
    """Image processing class
    This class contains methods for performing CCD image reduction for
    spectroscopy and imaging.
    """

    def __init__(self, args, data_container):
        """Initialization method for ImageProcessor class
        Args:
            args (Namespace): argparse instance
            data_container (DataFrame): Contains relevant information of the
            night and the data itself.
        """
        # TODO (simon): Check how inheritance could be used here.

        self.gain = data_container.gain
        self.rdnoise = data_container.rdnoise
        self.ccdsum = data_container.roi
        self.bias_buckets = data_container.bias_buckets
        self.flat_buckets = data_container.flat_buckets
        self.comp_buckets = data_container.comp_buckets
        self.targ_buckets = data_container.targ_buckets
        self.spec_buckets = data_container.spec_buckets
        self.sun_set = data_container.sun_set_time
        self.sun_rise = data_container.sun_rise_time
        self.morning_twilight = data_container.morning_twilight
        self.evening_twilight = data_container.evening_twilight
        self.pixel_scale = 0.15 * u.arcsec
        self.queue = None
        self.trim_section = None
        self.overscan_region = None
        self.master_bias = None
        self.master_bias_name = None
        self.out_prefix = None


        for data_bucket in [self.bias_bucket,
                            self.flat_bucket,
                            self.targ_buckets,
                            self.spec_buckets]

            if data_bucket is not None:
