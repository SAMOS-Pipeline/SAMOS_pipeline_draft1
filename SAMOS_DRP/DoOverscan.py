from __future__ import print_function
import glob
import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
import astroscrappy
from astropy import units
from SAMOSHelpers import MakeThumbnail
import ccdproc
from ccdproc import CCDData, ImageFileCollection

# Must read the data array
# Must convert to float32
# Must analyze bias level (as a function of row and/or column) and subtract
# Trim the array
# Write to output
def ParseSec(r):
    (x0,x1),(y0,y1) = map(lambda x: map(int,x.split(":")), r[1:-1].split(","))
    return (x0,x1),(y0,y1)


def Bias(bias_frames,corr_dir):


    master_list = []
    for bfu in bias_frames:
        bccd = CCDData.read(bfu,unit='count')
        (xb0,xb1),(yb0,yb1) = np.subtract(ParseSec(ccd.header["TRIMSEC"]),1)
        obccd = ccdproc.subtract_overscan(ccd=bccd,median=True,overscan=bccd[:,xb0:xb1],overscan_axis=1)

        tobccd = ccdproc.trim_image(ccd=obccd,fits_section=obccd.header['TRIMSEC'])
        master_list.append(tobccd)

    master_bais = ccdproc.combine(master_list,method='median')

    #for n in range(len(bias_frames)):
    #    master_bias.header.add_history('{}'.format(bias_frames[n]))

    master_bias_name = "master_bias_{}.fits".format(
        "x".join(master_bais.header['CCDSUM'].split()))

    master_bias_path = os.path.join(corr_dir,master_bias_name)
    master_bais.write(master_bias_path,overwrite=True)
    print("Master bias frame %s created and put in processing dir"%master_bias_name)
    return master_bias_path


def Overscan(input,output,bias):

    ccd = CCDData.read(input,unit='count')
    serial_binning,parallel_binning = [int(x) for x in ccd.header['CCDSUM'].split()]
    #will incorporate binning later when we have SAMOS data.
    print("Working on %s" %(input))
    pixel_cols = ccd.header['NAXIS1']
    pixel_rows = ccd.header['NAXIS2']

    print(ccd.header['TRIMSEC'])

   #by inspection, overscan region is 16 pixels in the dispersion direction.
    #overscan and data regions both span the full y-range of pixels.
    (xb0,xb1),(yb0,yb1) = np.subtract(ParseSec(ccd.header["TRIMSEC"]),1)
    ccd_os = ccdproc.subtract_overscan(ccd=ccd,median=True,overscan=ccd[:,xb0:xb1],overscan_axis=1)
    ccd_os_trim = ccdproc.trim_image(ccd=ccd_os,fits_section=ccd.header['TRIMSEC'])
    ccd = ccd_os_trim.copy()
    if bias is not None:
        mbias = CCDData.read(bias,unit='adu')
        ccd_os_trim_bias = ccdproc.subtract_bias(ccd=ccd_os_trim,master=mbias)
        ccd = ccd_os_trim_bias.copy()

    #clean cosmic rays here bc I need to get rid of them somehow
    #mask,ccd.data = astroscrappy.detect_cosmics(ccd.data,cleantype='medmask')

    ccd.header['DRP_COSM'] = ('LACosmic', "Cosmic rejection method")

    print("Writing %s" % (output))

    ccd.header.add_history("trimmed and overscan and bias subtracted")
    ccd.write(output,overwrite=True)
    print("%s written." % (output))

    p = os.path.dirname(output)
    pj = "%s/jpeg" % (p)
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)

if __name__ == "__main__":
   input,output = sys.argv
   Overscan(input,output)
