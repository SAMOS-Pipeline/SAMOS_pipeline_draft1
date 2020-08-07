import os
import ccdproc
from ccdproc import CCDData
import numpy as np
import logging
import warnings
warnings.filterwarnings('ignore')
from SAMOS_DRP.SAMOS_mods import add_wcs_keys,read_fits
import SAMOS_DRP.Spectroscopy.wcs
import SAMOS_DRP.Spectroscopy.wavelength
from SAMOS_DRP.Spectroscopy.wavelength import WavelengthCalibration

class WaveCalBuckets:
    """
    This class performs wavelength calibration for
    each slitobj/complamp pair.
    __init__ arg requires the slit bucket instance
    as input.
    """

    def __init__(self,slit_buckets):

        self.log = logging.getLogger(__name__)
        self.raw_data_dir = slit_buckets.raw_data_dir
        self.processing_dir = slit_buckets.processing_dir
        self.work_dir = slit_buckets.work_dir
        self.gain = slit_buckets.gain
        self.rdnoise = slit_buckets.rdnoise
        self.ccdsum = slit_buckets.ccdsum
        self.empty_bucket = True
        self.slit_num = None
        self.slit_targs = slit_buckets.slit_targs
        self.slit_comps = slit_buckets.slit_comps
        self.spec_buckets = slit_buckets.spec_buckets
        self.wavecal_spec_buckets = None
        self.in_prefix_list = slit_buckets.out_prefix_list
        self.out_prefix_list = None
        self.reference_lamps = os.path.join(os.path.split(self.work_dir)[0],'comp_refs')
        print("reference lamp location: %s"%(self.reference_lamps))
        #sNNN_prefix_filename.fits, where N
        #is the slit number in 000 format
    def __call__(self):

        #compute wavelength solution From
        #comparison lamps and apply to science
        #image for each slit.  resulting
        #WaveCal instances are saved to a list

        wavecals = []

        for slit in range(len(self.slit_targs)):

            Tccd = read_fits(self.slit_targs[slit])
            Cccd = read_fits(self.slit_comps[slit])
            wcs_starg = add_wcs_keys(Tccd)
            wcs_scomp = add_wcs_keys(Cccd)
            wavecal = WavelengthCalibration()

            targ1D = CCDData(np.median(Tccd.data,axis=0),unit='adu')
            targ1D.header = wcs_starg.header.copy()
            comp1D = CCDData(np.median(Cccd.data,axis=0),unit='adu')
            comp1D.header = wcs_scomp.header.copy()

            this_wavecal = wavecal(ccd=targ1D,
                                  comp_list=[comp1D],
                                  save_data_to=self.processing_dir,
                                  reference_data=self.reference_lamps,
                                  object_number=1,
                                  plot_results=True,
                                  save_plots=True,
                                  plots=False)
            wavecals.append(this_wavecal)

        self.wavecal_spec_buckets = wavecals
