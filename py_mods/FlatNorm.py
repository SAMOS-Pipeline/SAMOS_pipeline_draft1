import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
from SAMOSHelpers import *

# Combine the dome flats into  a single master flat.
# Fit the average spectrum of the master flat.
# Divide the master flat field frame by the fitted spectrum.
# Divide each data frame by the normalized flat field frame master_flat.fits.


def MkFlatNorm(input,output):


    # Create list of data corresponding to each flat.  Looks like [data_flat1, data_flat2, data_flat3, .....]
    # Scale each image by its median, so final list is [data_flat1/median(data_flat1), data_flat2/median(data_flat2),...]

    flats_data = np.array([fits.open(i)[0].data for i in input])
    scaled_flats_data = np.array([data/np.median(data) for data in flats_data])



    # Median combine and normalize by the mean to measure pixel to pixel variation in sensitivity of detector.
    med_comb_flat = np.median(scaled_flats_data,axis=0)
    master_flat = med_comb_flat/np.mean(med_comb_flat)


    hdu = fits.PrimaryHDU(master_flat.astype("f"))
    hdu.header = fits.open(input[0])[0].header.copy()
    hdu.header.add_history('Combined and normalized flat field')
    hdu.writeto(output,overwrite=True)

    pj = "%s/jpeg" % (os.path.split(output)[0])
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)
    print('Normalized master flat %s has been created.\n  Ready to be divided out of bias subtracted data frames.'%(output))


def NormFlatDiv(input,mflat,output):

    data,header = fits.getdata(input,header=True)
    mflat_data,mheader = fits.getdata(mflat,header=True)
    out_data = data/mflat_data
    #print(out_data.shape)
    #print(input)


    hdu = fits.PrimaryHDU(out_data.astype("f"))
    hdu.header = header.copy()
    hdu.header.add_history('Flat Fielded')
    hdu.writeto(output,overwrite=True)

    pj = "%s/jpeg" % (os.path.split(output)[0])
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)

    print("Flat fielded data has been created and stored in %s" %(output))

    print("Shape of each data frame is: ", (out_data.shape) )
