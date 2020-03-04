
import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
import pandas as pd
import ccdproc
from ccdproc import CCDData,ImageFileCollection
import logging

log = logging.getLogger(__name__)

def MuyMalo(x):
    print(x)
    os._exit(1)

def IsItHere(x):
    if not os.path.exists(x):
       print("Error:  %s not found!"%(x))
       MuyMalo("Quitting")

def MakeThumbnail(input,jdir,s1=-10,s2=10):
    f = fits.open(input)
    d = f[0].data.astype("f")
    m = np.median(d)
    s = 1.49*np.median(np.abs(d-m))
    z1 = m-s1*s
    z2 = m+s1*s
    pixels = np.clip(255*(d-z1)/(z2-z1),0,255).astype("b")
    pixmap = P.frombytes("L",(pixels.shape[1],pixels.shape[0]),pixels.ravel().tostring())
    image = P.merge("L",[pixmap])
    jpeg = "%s/%s" % (jdir,os.path.basename(input.replace(".fits",".jpg")))
    image.save(jpeg,"jpeg")


def check_header_notes(header):
    """
    This function checks the notes of the header to sort alignment/field images
    from science frames.  For now, this is only necessary for the Goodman MOS
    test data being used.
    """
    non_sci_flag_words = ["field", "align", "focus", "offset","testing"]
    acq_comment = "Triggered Acquisition"
    header_words = np.asarray(header["NOTES"].split(' '))
    if header["NOTES"]=='':
        return True
    elif any(word in non_sci_flag_words for word in header_words) or \
            (header['COMMENT']==acq_comment):
        #print(header_words)
        return True
    else:
        return False

def is_file_saturated(ccd, threshold=20.):
    """Detects a saturated image
    It counts the number of pixels above the saturation level, then finds
    which percentage they represents and if it is above the threshold it
    will return True. The percentage threshold can be set using the command
    line argument ``--saturation``.
    Args:
        ccd (CCDData): Image to be tested for saturation
        threshold (float): Percentage of saturated pixels allowed. Default 1.
    Returns:
        True for saturated and False for non-saturated
    """

    pixels_above_saturation = np.count_nonzero(
        ccd.data[np.where(
            ccd.data > get_saturation_value(
                ccd=ccd))])

    total_pixels = np.count_nonzero(ccd.data)

    saturated_percent = (pixels_above_saturation * 100.) / total_pixels

    if saturated_percent >= float(threshold):
        log.warning(
            "The current image has more than {:.2f} percent "
            "of pixels above saturation level".format(float(threshold)))
        return True
    else:
        return False

def get_saturation_value(ccd):

    columns = ['camera',
                   'read_rate',
                   'analog_attn',
                   'gain',
                   'read_noise',
                   'half_full_well',
                   'saturates_before']

    saturation_table = [['Blue', 50, 0, 0.25, 3.33, 279600, True],
                        ['Blue', 50, 2, 0.47, 3.35, 148723, True],
                        ['Blue', 50, 3, 0.91, 3.41, 76813, True],
                        ['Blue', 100, 0, 0.56, 3.69, 124821, True],
                        ['Blue', 100, 2, 1.06, 3.72, 65943, True],
                        ['Blue', 100, 3, 2.06, 3.99, 33932, False],
                        ['Blue', 200, 0, 1.4, 4.74, 49928, False],
                        ['Blue', 200, 2, 2.67, 5.12, 26179, False],
                        ['Blue', 400, 0, 5.67, 8.62, 12328, False],
                        ['Red', 100, 3, 1.54, 3.45, 66558, True],
                        ['Red', 100, 2, 3.48, 5.88, 29454, False],
                        ['Red', 344, 3, 1.48, 3.89, 69257, True],
                        ['Red', 344, 0, 3.87, 7.05, 26486, False],
                        ['Red', 750, 2, 1.47, 5.27, 69728, True],
                        ['Red', 750, 2, 1.45, 5.27, 69728, True],
                        ['Red', 750, 0, 3.77, 8.99, 27188, False],
                        ['Red', 750, 0, 3.78, 8.99, 27188, False]]

    _sdf = pd.DataFrame(saturation_table,
                                 columns=columns)

    hfw = _sdf.half_full_well[
            (_sdf.camera == 'Red') & #test data camera is red
            (_sdf.gain == ccd.header['GAIN']) &
            (_sdf.read_noise == ccd.header['RDNOISE'])]

    if hfw.empty:
        log.critical('Unable to obtain saturation level')
        saturation = None
        return None
    else:
        saturation = float(hfw.to_string(index=False))
        log.debug("Set saturation level as {:.0f}".format(
            saturation))
        return saturation

def save_bucket_status(this_bucket,bucket_fname):

    """
    This function takes an fuel instance created or updated from one of the data reduction steps
    and dumps it into a file for future access.  This will useful for if the user wants to stop
    the pipeline and pick back up instead of restarting the reduction process from scratch.
    """
    dump_dir = this_bucket.processing_dir
    filename = '%s/%s.txt'%(dump_dir,bucket_fname)
    outfile = open(filename, 'w')
    try:
        pickle.dump(this_bucket,outfile)
        log.info('WRITTEN: this reduction step %s/%s '%(os.path.basename(dump_dir),os.path.split(filename)[1]))

    except:
        log.warning('WARNING: could not write this fuel step to file %s/%s'%(os.path.basename(dump_dir),os.path.split(filename)[1]))

    outfile.close()

    return


def load_bucket_status(bucketsave_dir,which_bucket):

    fname = '%s/%s.txt'%(bucketsave_dir,which_bucket)
    infile = open(fname, 'r')

    try:
        fuel = pickle.load(infile)
        log.info('RETRIEVED: this reduction step %s/%s '%(os.path.basename(fuelsave_dir),os.path.basename(fname)))
    except:
        log.warning("Could not load %s"%(fname))

    return fuel
