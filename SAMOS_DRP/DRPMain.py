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
from ccdproc import CCDData
import warnings
warnings.filterwarnings('ignore')
import logging
import argparse
import logging
from importlib import reload
from .drp_mods import identify_slits,cutout_slit,\
                      image_trim,write_fits

from .ImageProcessor import ImageProcessor
from .SAMOSHelpers import load_bucket_status,save_bucket_status
from .SlitBuckets import SlitBuckets

def get_args(arguments=None):
    """Handles the argparse library and returns the arguments
    The list of arguments can be found with running ``redspec -h`` or
    ``redspec --help``.
    Notes:
        The full list of arguments are not listed here as the may evolve in
        which case is impossible to keep this up to date.
    Returns:
        An object that contains all the variables parsed through the argument
        system
    """
    # getLogger without __name__ so that we get the root logger.
    log = logging.getLogger()

    parser = argparse.ArgumentParser(
        description="Extracts goodman spectra and does automatic wavelength "
                    "calibration.\nPipeline Version: {:s}".format(__version__))

    parser.add_argument('--data-path',
                        action='store',
                        default=os.getcwd(),
                        type=str,
                        metavar='<Source Path>',
                        dest='source',
                        help='Path for location of raw data. Default <./>')

    parser.add_argument('--proc-path',
                        action='store',
                        default=os.getcwd(),
                        type=str,
                        metavar='<Destination Path>',
                        dest='destination',
                        help='Path for destination of processed data. Default '
                             '<./>')

    parser.add_argument('--search-pattern',
                        action='store',
                        default='cfzsto',
                        type=str,
                        metavar='<Search Pattern>',
                        dest='pattern',
                        help="Pattern for matching the goodman's reduced data.")

    parser.add_argument('--output-prefix',
                        action='store',
                        default='w',
                        metavar='<Out Prefix>',
                        dest='output_prefix',
                        help="Prefix to add to calibrated spectrum.")

    parser.add_argument('--extraction',
                        action='store',
                        default='fractional',
                        type=str,
                        metavar='<Extraction Type>',
                        dest='extraction_type',
                        choices=['fractional', 'optimal'],
                        help="Method to perform extraction. 'fractional' or "
                             "'optimal'. Only fractional pixel extraction is "
                             "implemented. Default 'fractional'.")

    parser.add_argument('--fit-targets-with',
                        action='store',
                        default='moffat',
                        type=str,
                        dest='target_fit_model',
                        choices=['moffat', 'gaussian'],
                        help="Model to fit peaks found on spatial profile "
                             "while searching for spectroscopic targets.")

    parser.add_argument('--reference-files',
                        action='store',
                        default='data/ref_comp/',
                        metavar='<Reference Dir>',
                        dest='reference_dir',
                        help="Directory of Reference files location")

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug_mode',
                        help="Debugging Mode")

    parser.add_argument('--debug-plot',
                        action='store_true',
                        dest='debug_with_plots',
                        help="Debugging show debugging plots")

    parser.add_argument('--max-targets',
                        action='store',
                        dest='max_n_targets',
                        metavar='<max targets>',
                        type=int,
                        default=3,
                        help="Maximum number of targets to be found in a "
                             "single image. Default 3")

    parser.add_argument('--background-threshold',
                        action='store',
                        dest='background_threshold',
                        type=int,
                        default=3,
                        help="Multiplier for background level used to "
                             "discriminate usable targets. Default 3 times "
                             "background level")

    parser.add_argument('--save-plots',
                        action='store_true',
                        dest='save_plots',
                        help="Save all plots in a directory")

    # parser.add_argument('--combine',
    #                     action='store_true',
    #                     dest='combine',
    #                     help="Combine compatible data")

    parser.add_argument('--plot-results',
                        action='store_true',
                        dest='plot_results',
                        help="Show wavelength calibrated spectrum at the end.")

    parser.add_argument('--version',
                        action='store_true',
                        dest='show_version',
                        help="Show current version of the Goodman Pipeline")

    args = parser.parse_args(args=arguments)

    try:
        ref_full_path = os.path.join(
            os.path.dirname(sys.modules['goodman.pipeline'].__file__),
            args.reference_dir)
    except KeyError as error:
        log.debug("KeyError {:s}".format(str(error)))
        ref_full_path = os.path.join(
            os.path.dirname(sys.modules['goodman_pipeline'].__file__),
            args.reference_dir)
    if not os.path.isdir(ref_full_path):
        log.info("Reference files directory doesn't exist.")
        try:
            os.path.os.makedirs(ref_full_path)
            log.info('Reference Files Directory is: %s', ref_full_path)
            args.reference_dir = ref_full_path
        except OSError as err:
            log.error(err)
    else:
        args.reference_dir = ref_full_path

    if not os.path.isabs(args.source):
        args.source = os.path.join(os.getcwd(), args.source)
    if not os.path.isdir(args.source):
        log.error("Source Directory {:s} doesn't exist.".format(args.source))
        if 'test' not in parser.prog:
            parser.print_help()
        parser.exit(0, "Leaving the Program.")

    if not os.path.isabs(args.destination):
        args.destination = os.path.join(os.getcwd(), args.destination)

    if not os.path.isdir(args.destination):
        log.error("Destination folder doesn't exist.")
        try:
            os.path.os.makedirs(args.destination)
            log.info('Destination folder created: %s', args.destination)
        except OSError as err:
            log.error(err)
            parser.print_help()
            parser.exit(0, "Leaving the Program.")
    return args

class DRPMain:

    """Main pipeline reduction class.
    """

    def __init__(self,raw_data_dir,altworkdir=False):


        self.log = logging.getLogger(__name__)
        self.args = None
        self.wavelength_solution_obj = None
        self.wavelength_calibration = None
        self.reference = None

        if not altworkdir:
            proc_dir = os.path.join(os.getcwd(),'%s_processing'%(self.obsid))
        else:
            proc_dir = os.path.join(altworkdir,'%s_processing'%(self.obsid))
        if not os.path.exists(proc_dir):
            os.mkdir(proc_dir)
        self.processing_dir = proc_dir

    def __call__(self):

        """Call method for the main DRP class
        This method call the higher level functions in order to do the
        spectroscopic data reduction.
        Args:
            args (object): argparse.Namespace instance that contains all the
                arguments.
        """

        self.log.debug("Initializing reference data locator.")
        #self.reference = ReferenceData(reference_dir=self.args.reference_dir)
        # data_container instance of NightDataContainer defined in core
        self.log.debug("Calling data classification procedure.")
        data_container = classify_spectroscopic_data(
            path=self.args.source,
            search_pattern=self.args.pattern)

        if data_container.is_empty:
            self.log.debug("Received empty data container.")
            sys.exit("Unable to find or classify data.")
        else:
            self.log.debug("Received non-empty data container.")

        self.log.debug("Calling _run method for MainApp")
        self._run(data_container=data_container,
                  extraction_type=self.args.extraction_type,
                  target_fit_model=self.args.target_fit_model,
                  background_threshold=self.args.background_threshold)

        self.log.info("END")
