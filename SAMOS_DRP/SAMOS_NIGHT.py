# -*- coding: utf8 -*-
from __future__ import print_function
import sys
import os
import json
from astropy.io import fits
from astropy import units
from astropy.coordinates import SkyCoord
import glob
import pandas as pd
import ccdproc
from ccdproc import CCDData,ImageFileCollection
import logging
import argparse
from .DataBucket import DataBucket
import warnings
from .SAMOSHelpers import check_header_notes
warnings.filterwarnings('ignore')




class SAMOSNight:

    """Defines and initializes the variables for processing the data
    from a specific night with SAMOS.
    """

    def __init__(self,obsid,raw_data_dir,proc_dir,LOG_FILENAME,altworkdir=False,
                ignore_bias=False,ignore_flats=False):

        self.obsid = obsid

        self.raw_data_dir = raw_data_dir
        self.processing_dir = proc_dir #where pipeline products are stored
        self.work_dir = os.getcwd() #where the user is running the pipeline
                                    #usually one directory up from proc_dir
        self.ignore_bias = ignore_bias
        self.ignore_flats = ignore_flats
        self.FULL_bucket = None
        self.all_data_types = None
        self.data_buckets = DataBucket(obsid,raw_data_dir)
        self.header_keys = ['naxis',
                         'date',
                         'slit',
                         'date-obs',
                         'obstype',
                         'object',
                         'exptime',
                         'obsra',
                         'obsdec',
                         'grating',
                         'cam_targ',
                         'grt_targ',
                         'filter',
                         'filter2',
                         'gain',
                         'rdnoise',
                         'ccdsum',
                         'wavmode']



        if os.path.exists(LOG_FILENAME):
            os.remove(LOG_FILENAME)
        logging.basicConfig(filename=LOG_FILENAME,level=logging.DEBUG)
        #is this how this works or am I being redundant??
        self.log = logging.getLogger(__name__)
        if not os.path.exists(proc_dir):
            os.mkdir(proc_dir)
        elif os.listdir(proc_dir):
            self.log.debug("starting with a clean processing directory for %s"%(proc_dir))
            [os.remove(i) for i in glob.glob(os.path.join(proc_dir,"*.fits"))]


    def __call__(self):

        FULL_bucket = ImageFileCollection(self.raw_data_dir,self.header_keys)
        self.FULL_bucket = FULL_bucket.summary.to_pandas()

        print("the pipeline is sorting your buckets")

        self.FULL_bucket['radeg'] = ''
        self.FULL_bucket['decdeg'] = ''
        for i in self.FULL_bucket.index.tolist():
            SCs = SkyCoord(self.FULL_bucket.obsra.iloc[i],
                           self.FULL_bucket.obsdec.iloc[i],
                           unit=(units.hourangle,units.deg))
            radeg,decdeg = SCs.ra,SCs.dec

            self.FULL_bucket.iloc[
                i, self.FULL_bucket.columns.get_loc('radeg')] = \
                '{:.2f}'.format(radeg)

            self.FULL_bucket.iloc[
                i, self.FULL_bucket.columns.get_loc('decdeg')] = \
                '{:.2f}'.format(decdeg)
            # now we can compare using degrees

        readout_configurations = self.FULL_bucket.groupby(
            ['gain',
             'rdnoise',
             'ccdsum']).size().reset_index().rename(columns={0: 'count'})

        data_buckets_list = []
        for i in readout_configurations.index:
            self.log.info("Organizing data for this configuration: "
                           "Gain: {:.2f}, Noise: {:.2f}, ccdsum: {:s}"
                           "".format(readout_configurations.iloc[i]['gain'],
                                     readout_configurations.iloc[i]['rdnoise'],
                                     readout_configurations.iloc[i]['ccdsum']))
            if not self.data_buckets.empty_bucket:
                self.log.debug("Reset data container thank u")
                self.data_buckets = DataBucket(obsid=self.obsid,
                    data_path=self.raw_data_dir)

            self.data_buckets.set_readout(
                gain=readout_configurations.iloc[i]['gain'],
                rdnoise=readout_configurations.iloc[i]['rdnoise'],
                ccdsum=readout_configurations.iloc[i]['ccdsum'])

            configuration_buckets = self.FULL_bucket[
                ((self.FULL_bucket['gain'] ==
                  readout_configurations.iloc[i]['gain']) &
                 (self.FULL_bucket['rdnoise'] ==
                  readout_configurations.iloc[i]['rdnoise']) &
                 (self.FULL_bucket['ccdsum'] ==
                  readout_configurations.iloc[i]['ccdsum']))]

            self.all_data_types = configuration_buckets.obstype.unique()
            print(isinstance(self.data_buckets,DataBucket))
            print(vars(self.data_buckets))
            self.spectroscopy_buckets(FULL_bucket=configuration_buckets,
                                        data_buckets=self.data_buckets)


            if self.data_buckets.empty_bucket:
                self.log.warning("The following files will not be processed:")
                for _file in sub_collection['file'].tolist():
                    self.log.warning("{:s}".format(_file))
                self.log.debug('data_buckets is empty')
            else:
                self.log.info('Found valid data, appending to data container '
                              'list')
                data_buckets_list.append(self.data_buckets)

        if len(data_buckets_list) == 0:
            return [None]
        elif not all(data_buckets_list):
            self.log.warning("It is possible that there is no valid data.")
            return [None]
        else:
            return data_buckets_list




    def spectroscopy_buckets(self,FULL_bucket,data_buckets):
        """ Modified from the Goodman Pipeline! :)
        Organizes data for spectroscopy
        This method identifies all combinations of nine **key** keywords that
        can set apart different objects with their respective calibration data
        or not. The keywords used are:
          - ``GAIN``
          - ``RDNOISE``
          - ``GRATING``
          - ``FILTER2``
          - ``CAM_TARG``
          - ``GRT_TARG``
          - ``SLIT``
          - ``OBSRA``
          - ``OBSDEC``
        This method populates the `data_buckets` class attribute for
        the insance of this :class:`SAMOS_NIGHT.SAMOSNight.
        A data group is an instance of a :class:`pd.DataFrame`.
        """

        assert isinstance(FULL_bucket, pd.DataFrame)
        assert isinstance(data_buckets, DataBucket)

                # obtain a list of timestamps of observing time
        # this will only be used for naming flats
        dateobs_list = FULL_bucket['date-obs'].tolist()

        # get times for twilights, sunset an sunrise
        #afternoon_twilight, morning_twilight, sun_set, sun_rise = \
        #    get_twilight_time(date_obs=dateobs_list)

        # set times in data container
        #data_buckets.set_sun_times(sun_set,
        #                             sun_rise)
        #data_buckets.set_twilight_times(afternoon_twilight,
        #                                  morning_twilight)

        # process bias
        bias_buckets = FULL_bucket[FULL_bucket.obstype == 'BIAS']

        if not self.ignore_bias:
            if len(bias_buckets) == 0:
                self.log.critical('There are no BIAS images for this '
                                  'configuration. Use --ignore-bias '
                                  'to proceed without BIAS.')
                # sys.exit('CRITICAL ERROR: BIAS not Found.')
                return False
            else:
                bias_conf = bias_buckets.groupby(
                    ['gain',
                     'rdnoise',
                     'radeg',
                     'decdeg']).size().reset_index().rename(
                    columns={0: 'count'})

                # bias_conf
                for i in bias_conf.index:
                    bias_group = bias_buckets[(
                        (bias_buckets['gain'] == bias_conf.iloc[i]['gain']) &
                        (bias_buckets['rdnoise'] == bias_conf.iloc[i][
                            'rdnoise']) &
                        (bias_buckets['radeg'] == bias_conf.iloc[i][
                            'radeg']) &
                        (bias_buckets['decdeg'] == bias_conf.iloc[i][
                            'decdeg']))]

                    data_buckets.add_bias_bucket(bias_bucket=bias_group)
        else:
            self.log.warning('Ignoring BIAS by request.')

        if 'FLAT' not in FULL_bucket.obstype.unique() and \
                not self.ignore_flats:
            self.log.critical('There is no FLAT images. Use --ignore-flats to '
                              'continue without FLATs.')
            sys.exit('CRITICAL ERROR: FLAT not Found.')
        elif self.ignore_flats:
            self.log.warning('Ignoring FLAT images on request.')
            bucket_collection = FULL_bucket[
                ((FULL_bucket.obstype != 'BIAS') &
                 (FULL_bucket.obstype != 'FLAT'))]
        else:
            # process non-bias i.e. flats and object ... and comp
            bucket_collection = FULL_bucket[FULL_bucket.obstype != 'BIAS']

        configurations = bucket_collection.groupby(['slit',
                          'radeg',
                          'decdeg',
                          'grating',
                          'cam_targ',
                          'grt_targ',
                          'filter',
                          'filter2',
                          'gain',
                          'rdnoise']).size().reset_index().rename(columns={0: 'count'})

        for i in configurations.index:

            data_bucket_group = bucket_collection[(
                (bucket_collection['slit'] == configurations.iloc[i]['slit']) &
                (bucket_collection['radeg'] == configurations.iloc[i]['radeg']) &
                (bucket_collection['decdeg'] == configurations.iloc[i]['decdeg']) &
                (bucket_collection['grating'] == configurations.iloc[i]['grating']) &
                (bucket_collection['cam_targ'] == configurations.iloc[i]['cam_targ']) &
                (bucket_collection['grt_targ'] == configurations.iloc[i]['grt_targ']) &
                (bucket_collection['filter2'] == configurations.iloc[i]['filter2']) &
                (bucket_collection['gain'] == configurations.iloc[i]['gain']) &
                (bucket_collection['rdnoise'] == configurations.iloc[i]['rdnoise']))]

            group_obstype = data_bucket_group.obstype.unique()
            if 'FLAT' in group_obstype and len(group_obstype) == 1:
                data_buckets.add_flat_bucket(data_bucket_group)
            else:
                data_buckets.add_spec_bucket(data_bucket_group)
                # Comparison lamps are processed as science data.
                if 'FLAT' in group_obstype:
                    print('flat in spec group %s'%(data_bucket_group))
                    # grab flats and put them in the flats group as well
                    object_flat_group = data_bucket_group[
                        data_bucket_group['obstype'] == 'FLAT']
                    data_buckets.add_flat_bucket(object_flat_group)

        return data_buckets


class GenerateDcrParFile:
    """
    from Goodman Pipeline.
    Creates dcr.par file based on lookup table
    `dcr` parameters depend heavily on binning, this class generates a file
    using the default format. The lookup table considers camera and binning.
    """
    _format = [
        "THRESH  = {:.1f} // Threshold (in STDDEV)",
        "XRAD    = {:d}   // x-radius of the box (size = 2 * radius)",
        "YRAD    = {:d}   // y-radius of the box (size = 2 * radius)",
        "NPASS   = {:d}   // Maximum number of cleaning passes",
        "DIAXIS  = {:d}   // Dispersion axis: 0 - no dispersion, 1 - X, 2 - Y",
        "LRAD    = {:d}   // Lower radius of region for replacement statistics",
        "URAD    = {:d}   // Upper radius of region for replacement statistics",
        "GRAD    = {:d}   // Growing radius",
        "VERBOSE = {:d}   // Verbose level [0,1,2]",
        "END"]

    _columns = ['parameter',
                'red-1',
                'red-2',
                'red-3',
                'blue-1',
                'blue-2',
                'blue-3']

    _lookup = [
        ['thresh', 3.0, 4.0, 3.0, 3.0, 3.0, 3.0],
        ['xrad', 9, 7, 9, 8, 9, 9],
        ['yrad', 9, 9, 9, 8, 9, 9],
        ['npass', 5, 5, 5, 5, 5, 5],
        ['diaxis', 0, 0, 0, 0, 0, 0],
        ['lrad', 1, 1, 1, 1, 1, 1],
        ['urad', 3, 3, 3, 3, 3, 3],
        ['grad', 1, 0, 1, 1, 1, 1],
        ['verbose', 1, 1, 1, 1, 1, 1]
    ]

    def __init__(self, par_file_name='dcr.par'):
        """
        Args:
            par_file_name:
        """
        self._file_name = par_file_name
        self._df = pd.DataFrame(self._lookup, columns=self._columns)
        self._binning = "{:s}-{:s}"
        self._data_format = "\n".join(self._format)

    def __call__(self, instrument='Red', binning='1', path='default'):
        """
        Args:
            instrument (str): Instrument from INSTCONF keyword
            binning (str): Serial (dispersion) Binning from the header.
            path (str): Directory where to save the file.
        """
        assert any([instrument == option for option in ['Red', 'Blue']])
        b = self._binning.format(instrument.lower(), binning)
        self._data_format = self._data_format.format(
            self._df[b][self._df.parameter == 'thresh'].values[0],
            int(self._df[b][self._df.parameter == 'xrad'].values[0]),
            int(self._df[b][self._df.parameter == 'yrad'].values[0]),
            int(self._df[b][self._df.parameter == 'npass'].values[0]),
            int(self._df[b][self._df.parameter == 'diaxis'].values[0]),
            int(self._df[b][self._df.parameter == 'lrad'].values[0]),
            int(self._df[b][self._df.parameter == 'urad'].values[0]),
            int(self._df[b][self._df.parameter == 'grad'].values[0]),
            int(self._df[b][self._df.parameter == 'verbose'].values[0]))
        self._create_file(path=path)

    def _create_file(self, path):
        """Creates `dcr.par` file
        Args:
            path (str): Path to where to save the `dcr.par` file.
        """
        if os.path.isdir(path):
            full_path = os.path.join(path, self._file_name)
        else:
            full_path = os.path.join(os.getcwd(), self._file_name)

        with open(full_path, 'w') as dcr_par:
            dcr_par.write(self._data_format)
