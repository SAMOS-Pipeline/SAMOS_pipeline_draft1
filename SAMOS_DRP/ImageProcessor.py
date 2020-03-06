# -*- coding: utf8 -*-
import sys
import os
import json
from astropy.io import fits
from astropy import units
import glob
import pandas as pd
import numpy as np
import ccdproc
from ccdproc import CCDData,ImageFileCollection
import logging
import argparse
from .DataBucket import DataBucket
from .SAMOSHelpers import save_bucket_status,MuyMalo
from .SAMOS_NIGHT import SAMOSNight
from .drp_mods import classify_spectroscopic_data
from .drp_mods import (define_trim_section,get_overscan_region,
                     image_trim,image_overscan,read_fits,write_fits,
                     create_master_bias,create_master_flats,name_master_flats,
                     normalize_master_flat,get_best_flat,combine_data,
                     call_cosmic_rejection,identify_slits)
import SAMOS_DRP.Spectroscopy.wcs
import SAMOS_DRP.Spectroscopy.wavelength
from SAMOS_DRP.Spectroscopy.wavelength import WavelengthCalibration

import warnings
warnings.filterwarnings('ignore')

log = logging.getLogger(__name__)
class ImageProcessor:

    """Methods for CCD image reduction.  Based on Goodman Pipeline ImageProcessor.
    """

    def __init__(self,night_bucket,ignore_bias=True,
                ignore_flats=False):
    #args are [ignore_bias,ignore_flats,raw data,intermediate data]

        self.log = logging.getLogger(__name__)
        self.ignore_bias = ignore_bias
        self.ignore_flats = ignore_flats
        print(night_bucket.raw_data_dir)
        self.raw_data_dir = night_bucket.raw_data_dir
        self.processing_dir = night_bucket.processing_dir
        self.gain = night_bucket.FULL_bucket.gain
        self.rdnoise = night_bucket.FULL_bucket.rdnoise
        self.ccdsum = night_bucket.FULL_bucket.ccdsum
        self.data_buckets = night_bucket.data_buckets
        self.bias_buckets = night_bucket.data_buckets.bias_buckets
        self.flat_buckets = night_bucket.data_buckets.flat_buckets
        self.targ_buckets = night_bucket.data_buckets.targ_buckets
        self.spec_buckets = night_bucket.data_buckets.spec_buckets
        self.slit_buckets = night_bucket.data_buckets.slit_buckets
        self.combined_bucket = None
        self.slit_edges_x = None
        self.slit_edges_y = None
        self.sun_set = night_bucket.data_buckets.sunset_time
        self.sun_rise = night_bucket.data_buckets.sunrise_time
        self.morning_twilight = night_bucket.data_buckets.morning_twilight
        self.evening_twilight = night_bucket.data_buckets.evening_twilight
        self.pixel_scale = 0.15 * units.arcsec
        self.clean_cosmic = False
        self.queue = None
        self.trim_section = None
        self.overscan_region = None
        self.master_bias = None
        self.master_bias_name = None
        self.out_prefix = None
        self.combine = True #keeping this here in case I want to use combine
                             #later.

    def __call__(self):

        """Call method for the main DRP class
        This method call the higher level functions in order to do the
        spectroscopic data reduction.
        Args:
            args (object): argparse.Namespace instance that contains all the
                arguments.
        """

        self.log.debug("Processing ccd images.")

        # data_bucket instance of NightDataContainer defined in core


        for bucket in [self.bias_buckets,self.flat_buckets,
                          self.spec_buckets]:
            if bucket is not None:

                image_list = bucket[0]['file'].tolist()
                sample_image = os.path.join(self.raw_data_dir,
                                            np.random.choice(image_list))

                self.trim_section = define_trim_section(image=sample_image)


                self.overscan_region = get_overscan_region(image=sample_image)


                log.info("trimsec: %s,\n osec: %s"%(self.trim_section,
                                                    self.overscan_region))

                for sub_bucket in bucket:
                    #print(sub_bucket[['file','obstype','gain','rdnoise']])
                    sbucket_obstype = sub_bucket.obstype.unique()
                    if len(sbucket_obstype) == 1 and sbucket_obstype[0] == 'BIAS':

                        bias_files = sub_bucket.file.tolist()
                        self.master_bias, self.master_bias_name = \
                            create_master_bias(
                                bias_files=bias_files,
                                raw_data=self.raw_data_dir,
                                reduced_data=self.processing_dir)


                    elif len(sbucket_obstype) == 1 and sbucket_obstype[0] == 'FLAT':

                        flat_files = sub_bucket.file.tolist()
                        sample_header = fits.getheader(os.path.join(
                            self.raw_data_dir, flat_files[0]))
                        master_flat_name = name_master_flats(
                            header=sample_header,
                            reduced_data=self.processing_dir,
                            sun_set=self.sun_set,
                            sun_rise=self.sun_rise,
                            evening_twilight=self.evening_twilight,
                            morning_twilight=self.morning_twilight)

                        create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.raw_data_dir,
                            reduced_data=self.processing_dir,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            ignore_bias=self.ignore_bias)

                    else:
                        log.debug('Process Data Group')
                        self.process_spectroscopy_science(sub_bucket)


    def process_spectroscopy_science(self, science_bucket, save_all=True):
        """Process Spectroscopy science images.
        This function handles the full image reduction process for science
        files. if save_all is set to True, all intermediate steps are saved.
        Args:
            science_bucket (object): :class:`~pandas.DataFrame` instance that contains a
                list of science images that where observed at the same pointing
                and time. It also contains a set of selected keywords from the
                image's header.
            save_all (bool): If True the pipeline will save all the intermadiate
                files such as after overscan correction or bias corrected and
                etc.
        """
        # TODO (simon): The code here is too crowded.
        # TODO cont. Create other functions as necessary

        target_name = ''
        slit_trim = None
        master_bias = None
        master_flat = None
        master_flat_name = None
        obstype = science_bucket.obstype.unique()
        # print(obstype)
        if 'OBJECT' in obstype or 'COMP' in obstype:
            object_comp_group = science_bucket[
                (science_bucket.obstype == 'OBJECT') |
                (science_bucket.obstype == 'COMP')]

            print(object_comp_group)

            if 'OBJECT' in obstype:
                target_name = science_bucket.object[science_bucket.obstype ==
                                                   'OBJECT'].unique()[0]

                log.info('Processing Science Target: '
                              '{:s}'.format(target_name))
            else:
                # TODO (simon): This does not make any sense
                log.info('Processing Comparison Lamp: '
                              '{:s}'.format(target_name))

            if 'FLAT' in obstype and not self.ignore_flats:
                flat_sub_group = science_bucket[science_bucket.obstype == 'FLAT']
                flat_files = flat_sub_group.file.tolist()
                sample_header = fits.getheader(os.path.join(
                    self.raw_data_dir, flat_files[0]))
                master_flat_name = name_master_flats(
                    header=sample_header,
                    reduced_data=self.processing_dir,
                    sun_set=self.sun_set,
                    sun_rise=self.sun_rise,
                    evening_twilight=self.evening_twilight,
                    morning_twilight=self.morning_twilight)

                master_flat, master_flat_name = \
                    create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.raw_data_dir,
                            reduced_data=self.processing_dir,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            ignore_bias=self.ignore_bias)
            elif self.ignore_flats:
                log.warning('Ignoring creation of Master Flat by request.')
                master_flat = None
                master_flat_name = None
            else:
                log.info('Attempting to find a suitable Master Flat')
                object_list = object_comp_group.file.tolist()

                # grab a random image from the list
                random_image = np.random.choice(object_list)

                # define random image full path
                random_image_full = os.path.join(self.raw_data_dir,
                                                 random_image)

                # read the random chosen file
                ccd = CCDData.read(random_image_full)

                if not self.ignore_flats:
                    # define the master flat name
                    master_flat_name = name_master_flats(
                        header=ccd.header,
                        reduced_data=self.processing_dir,
                        sun_set=self.sun_set,
                        sun_rise=self.sun_rise,
                        evening_twilight=self.evening_twilight,
                        morning_twilight=self.morning_twilight,
                        get=True)

                    # load the best flat based on the name previously defined
                    master_flat, master_flat_name = \
                        get_best_flat(flat_name=master_flat_name,
                                      path=self.processing_dir)
                    if (master_flat is None) and (master_flat_name is None):
                        log.critical('Failed to obtain master flat')
                    try:

                        log.debug('Attempting to find slit trim section')
                        slit_reference_file = os.path.join(\
                                        os.path.split(self.working_dir)[0],'slit_refs')
                        self.slit_edges_x,self.slit_edges_y = \
                                    identify_slits(master_flat=master_flat,
                                    slit_reference_file=slit_reference_file)
                    except AttributeError:
                        log.critical("Master flat inexistent, can't find slit trim "
                                      "section, exiting")
                        MuyMalo("something went wrong, no slit edges found"
                                      "section")
                    #slit_flat_arrs,slit_trimsec = cutout_slit(master_flat,x_slit_edges,y_slit_edges)

            norm_master_flat = None
            norm_master_flat_name = None
            all_object_images = []
            all_comp_images = []
            all_slit_images = []
            for science_image in object_comp_group.file.tolist():
                self.out_prefix = ''

                # define image full path
                image_full_path = os.path.join(self.raw_data_dir,
                                               science_image)

                # load image
                ccd = read_fits(image_full_path)

                # apply overscan

                ccd = image_overscan(ccd, overscan_region=self.overscan_region)
                self.out_prefix += 'o_'

                if save_all:
                    full_path = os.path.join(self.processing_dir,
                                             self.out_prefix + science_image)

                    # ccd.write(full_path, clobber=True)
                    write_fits(ccd=ccd, full_path=full_path)

                ccd = image_trim(ccd=ccd,
                                 trim_section=self.trim_section,
                                 trim_type='trimsec')

                self.out_prefix = 't' + self.out_prefix

                if save_all:
                    full_path = os.path.join(
                        self.processing_dir,
                        self.out_prefix + science_image)

                    # ccd.write(full_path, clobber=True)
                    write_fits(ccd=ccd, full_path=full_path)

                if not self.ignore_bias:
                    # TODO (simon): Add check that bias is compatible
                    print("master bias shape: {}\n current ccd shape: {}".format(
                                               master_bias.shape,ccd.shape))

                    master_bias = image_trim(ccd=master_bias,
                                        trim_section=self.trim_section,
                                        trim_type='trimsec')
                    ccd = ccdproc.subtract_bias(ccd=ccd,
                                                master=master_bias,
                                                add_keyword=False)

                    self.out_prefix = 'z' + self.out_prefix

                    ccd.header['SRP_BIAS'] = (
                        os.path.basename(self.master_bias_name),
                        'Master bias image')

                    if save_all:
                        full_path = os.path.join(
                            self.processing_dir,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)
                else:
                    log.warning('Ignoring bias correction by request.')

                # Do flat correction
                if master_flat is None or master_flat_name is None:
                    log.warning('The file {:s} will not be '
                                     'flatfielded'.format(science_image))
                elif self.ignore_flats:
                    log.warning('Ignoring flatfielding by request.')
                else:

                    if norm_master_flat is None:

                        norm_master_flat, norm_master_flat_name = \
                            normalize_master_flat(
                                master=master_flat,
                                name=master_flat_name,
                                method='simple')


                    ccd = ccdproc.flat_correct(ccd=ccd,
                                               flat=norm_master_flat,
                                               add_keyword=False)

                    self.out_prefix = 'f' + self.out_prefix

                    ccd.header['SRP_FLAT'] = (
                        os.path.basename(norm_master_flat_name),
                        'Master flat image')

                    # propagate master flat normalization method
                    ccd.header['SRP_NORM'] = norm_master_flat.header['SRP_NORM']

                    if save_all:
                        full_path = os.path.join(
                            self.processing_dir,
                            self.out_prefix + science_image)

                        # ccd.write(full_path, clobber=True)
                        write_fits(ccd=ccd, full_path=full_path)



                ccd, prefix = call_cosmic_rejection(
                    ccd=ccd,
                    image_name=science_image,
                    out_prefix=self.out_prefix,
                    red_path=self.processing_dir,
                    save=False)
                    #other vars are set to Goodman Pipeline
                    #defaults, I will add option things later

                if ccd is not None:
                    if ccd.header['OBSTYPE'] == 'OBJECT':
                        log.debug("Appending OBJECT image for combination")
                        all_object_images.append(ccd)
                    elif ccd.header['OBSTYPE'] == 'COMP':
                        log.debug("Appending COMP image for combination")
                        all_comp_images.append(ccd)
                    else:
                        log.error("Unknown OBSTYPE = {:s}"
                                       "".format(ccd.header['OBSTYPE']))
                else:
                    log.warning("Cosmic ray rejection returned a None.")
            combined_bucket = []
            if self.combine:
                log.warning("Combination of data is experimental.")
                if len(all_object_images) > 1:

                    log.info("Combining {:d} OBJECT images"
                                  "".format(len(all_object_images)))

                    object_group = object_comp_group[
                        object_comp_group.obstype == "OBJECT"]


                    #print(object_group, len(all_object_images))

                    combined_image,combined_full_path = \
                                combine_data(all_object_images,
                                 dest_path=self.processing_dir,
                                 prefix=self.out_prefix,
                                 save=True)
                    combined_object_row = object_group.iloc[0].copy()
                    combined_object_row['file'] = os.path.split(combined_full_path)[1]
                    combined_bucket.append(pd.DataFrame(
                                        data=[combined_object_row.values.flatten()],
                                        columns=combined_object_row.index))
                elif len(all_object_images) == 1:
                    # write_fits(all_object_images[0])
                    pass

                else:

                    log.error("No OBJECT images to combine")

                if len(all_comp_images) > 1:
                    log.info("Combining {:d} COMP images"
                                  "".format(len(all_comp_images)))
                    comp_group = object_comp_group[
                         object_comp_group.obstype == "COMP"]
                    combined_image,combined_full_path = \
                                 combine_data(all_comp_images,
                                 dest_path=self.processing_dir,
                                 prefix=self.out_prefix,
                                 save=True)
                    combined_comp_row = comp_group.iloc[0].copy()
                    combined_comp_row['file'] = os.path.split(combined_full_path)[1]
                    combined_bucket.append(pd.DataFrame(
                                            data=[combined_comp_row.values.flatten()],
                                            columns=combined_comp_row.index))
                else:
                    log.error("No COMP images to combine")
            else:
                log.debug("Combine is disabled (default)")

            if len(combined_bucket)>1:
                self.combined_bucket = pd.concat(combined_bucket)
        elif 'FLAT' in obstype:
            self.queue.append(science_bucket)
            log.warning('Only flats found in this group')
            flat_sub_group = science_bucket[science_bucket.obstype == 'FLAT']
            # TODO (simon): Find out if these variables are useful or not
            flat_files = flat_sub_group.file.tolist()

            sample_header = fits.getheader(os.path.join(self.raw_data_dir,
                                                        flat_files[0]))
            master_flat_name = name_master_flats(
                header=sample_header,
                reduced_data=self.processing_dir,
                sun_set=self.sun_set,
                sun_rise=self.sun_rise,
                evening_twilight=self.evening_twilight,
                morning_twilight=self.morning_twilight)

            create_master_flats(
                            flat_files=flat_files,
                            raw_data=self.raw_data_dir,
                            reduced_data=self.processing_dir,
                            overscan_region=self.overscan_region,
                            trim_section=self.trim_section,
                            master_bias_name=self.master_bias_name,
                            new_master_flat_name=master_flat_name,
                            saturation=self.saturation_threshold)
        else:
            log.error('There is no valid datatype in this group')
