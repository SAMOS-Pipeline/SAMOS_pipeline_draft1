from __future__ import print_function
import glob
from SAMOSHelpers import *
import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
import astroscrappy

# Must read the data array
# Must convert to float32
# Must analyze bias level (as a function of row and/or column) and subtract
# Trim the array
# Write to output

def ParseSec(r):
    (x0,x1),(y0,y1) = map(lambda x: map(int,x.split(":")), r[1:-1].split(","))
    return (x0,x1),(y0,y1)

def FieldTrim(input,output):
    """
    This function trims the mask field image to the proper size.  I have to use this to
    create the slit_edges text file for now, since the slit positions aren't in the FITS headers.
    """
    f = fits.open(input)
    h = f[0].header
    d = f[0].data.astype("f")
    theta = 180.
    rot_matrix = np.matrix([[int(np.cos(theta)),int(-np.sin(theta))],
                [int(np.sin(theta)),int(np.cos(theta))]])

    print("Working on %s" %(input))
    (xb0,xb1),(yb0,yb1) = np.subtract(ParseSec(h["biassec"]),1)
    (xd0,xd1),(yd0,yd1) = np.subtract(ParseSec(h["datasec"]),1)

    data = d[1300:2800]
    if 'c1' in input:
        data = np.rot90(data,2)
    elif 'c2' in input:
        data = np.fliplr(np.rot90(data,2))


    print("Writing %s" % (output))
    hdu = fits.PrimaryHDU(data.astype("f"))
    hdu.header = h.copy()
    hdu.header["bitpix"] = -32
    hdu.header["NAXIS1"] = data.shape[1]
    hdu.header["NAXIS2"] = data.shape[0]
    [hdu.header.remove(key) for key in ["bscale","bzero"]]
    hdu.header.add_history("trimmed")
    hdu.writeto(output,overwrite=True)
    print("%s written." % (output))

    p = os.path.dirname(output)
    pj = "%s/jpeg" % (p)
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)

def Overscan(input,output):
    f = fits.open(input)
    h = f[0].header
    d = f[0].data.astype("f")
    theta = 180.
    rot_matrix = np.matrix([[int(np.cos(theta)),int(-np.sin(theta))],
                [int(np.sin(theta)),int(np.cos(theta))]])

    print("Working on %s" %(input))


    (xb0,xb1),(yb0,yb1) = np.subtract(ParseSec(h["biassec"]),1)
    (xd0,xd1),(yd0,yd1) = np.subtract(ParseSec(h["datasec"]),1)
    print(input)


    overscan = d[yb0:yb1+1,xb0:xb1+1]

#    o = np.add.reduce(np.transpose(overscan),0)/overscan.shape[1] # MEAN
    o = np.median(overscan,axis=1) # MEDIAN
    data = d[yd0:yd1+1,xd0:xd1+1] - o[::,np.newaxis]

# KLUGE FOR PRACTICING WITH THE LDSS3 NOD-SHUFFLE DATA FROM TOM

    data = data[1300:2800]
    if 'c1' in input:
        data = np.rot90(data,2)
    elif 'c2' in input:
        data = np.fliplr(np.rot90(data,2))
#clean cosmic rays here bc I need to get rid of them somehow
    mask = np.ma.make_mask(data, copy=True,shrink=True,dtype=np.bool)
    mask[:,:] = False
    crmask, dataCR = astroscrappy.detect_cosmics(data,inmask=mask,cleantype='medmask')


    print("Writing %s" % (output))
    hdu = fits.PrimaryHDU(dataCR.astype("f"))
    hdu.header = h.copy()
    hdu.header["bitpix"] = -32
    hdu.header["NAXIS1"] = dataCR.shape[1]
    hdu.header["NAXIS2"] = dataCR.shape[0]
    [hdu.header.remove(key) for key in ["bscale","bzero"]]
    hdu.header.add_history("overscan subtracted")
    hdu.writeto(output,overwrite=True)
    print("%s written." % (output))

    p = os.path.dirname(output)
    pj = "%s/jpeg" % (p)
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)


def stitch(input):
    #take the bias corrected chip images and combine them.

    basename = os.path.basename(input)[:-7]
    outfile = input[:-7]+'.fits'
    if "c2" in input:
        fc2 = fits.open(input)
        hc2 = fc2[0].header
        dc2 = fc2[0].data.astype("f")

        fc1 = fits.open(input[:-7]+'c1.fits')
        hc1 = fc1[0].header
        dc1 = fc1[0].data.astype("f")
        dfull = np.concatenate((dc2,dc1),axis=1)
        dfull = np.fliplr(dfull)
#now clean some bad columns from the stitched data
        mask = np.ma.make_mask(dfull, copy=True,shrink=True,dtype=np.bool)
        mask[:,:] = False
        mask[:,1602:1607] = True
        mask[564:613,1548:1553] = True
        mdfull = np.ma.masked_array(dfull,mask=mask,fill_value=np.nan)
        s = 2
        for row in range(mdfull.shape[0]):
            for col in range(mdfull.shape[1]):
                if np.isnan(mdfull[row,col]):
                    y1 = row-s
                    y2 = row+s+1
                    x1 = col-s
                    x2 = col+s+1
                    if y1<0.0:
                        y1 = 0.0
                    if y2>mdfull.shape[0]:
                        y2 = mdfull.shape[0]
                    if x1<0.0:
                        x1 = 0.0
                    if x2>mdfull.shape[1]:
                        x2 = mdfull.shape[1]
                    dfull[row,col] = np.med(mdfull[y1:y2,x1:x2],axis=0)
        
        hdu = fits.PrimaryHDU(dfull.astype("f"))
        hdu.header = hc2.copy()
        hdu.header["NAXIS1"] = dfull.shape[1]
        hdu.header["NAXIS2"] = dfull.shape[0]
        hdu.header.add_history("stitched")
        hdu.writeto(outfile,overwrite=True)
    else:
        pass


    return outfile


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

if __name__ == "__main__":
   input,output = sys.argv[1:3]
   Overscan(input,output)
