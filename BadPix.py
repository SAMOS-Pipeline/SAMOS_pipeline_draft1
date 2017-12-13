import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
from SAMOSHelpers import *
import astroscrappy


# First, find bad pixel columns.
# Create bad pixel mask.
# Fill the pixels with the median of the surrounding good pixels.

def fixbad(input,output):

    data,header = fits.getdata(input,header=True)
    
    mask = np.ma.make_mask(data, copy=True,shrink=True,dtype=np.bool)

    mask[:,:] = False
    
    crmask, dataCR = astroscrappy.detect_cosmics(data,inmask=mask,cleantype='medmask')
    
    hdu = fits.PrimaryHDU(dataCR.astype("f"))
    hdu.header = header.copy()
    hdu.header.add_history('Cosmic Ray corrected.')
    hdu.writeto(output,overwrite=True)
    
    #print(os.path.split(output)[0])
    
    pj = "%s/jpeg" % (os.path.split(output)[0])
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)