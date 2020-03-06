from SAMOS_DRP.drp_mods import add_wcs_keys,read_fits
import SAMOS_DRP.Spectroscopy.wcs
import SAMOS_DRP.Spectroscopy.wavelength
from SAMOS_DRP.Spectroscopy.wavelength import WavelengthCalibration

import warnings
warnings.filterwarnings('ignore')

def do_wavecal_all_targs(self):
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
        this_wavecal = wavecal(ccd=targ1D,comp_list=[comp1D],save_data_to=self.processing_dir,reference_data='comp_refs',
                                  object_number=1,plot_results=True,save_plots=True,plots=False)
        wavecals.append(this_wavecal)

    self.wavecal_slits = wavecals
