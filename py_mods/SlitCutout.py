from __future__ import print_function
import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P
from SAMOSHelpers import *
import pyregion
import matplotlib.pyplot as plt
import regions
import matplotlib.cm as cm
from PIL import Image as P
from SlitID import get_edges
from numpy.polynomial.polynomial import polyval


'''
This is the main module responsible for extracting the slits.
Each slit will be put into its own FITS file for further data reduction.
'''

def cutout_slit(input,x_edges,y_edges):

    #x_edges, y_edges = get_edges(input,'LMask2_ycoords_c1.txt','LMask2')

    margin = 5 #extra space to include in cutout_bin
    slit_height = 18
    input_fpath = os.path.split(input)
    print(input_fpath)
    data,header = fits.getdata(input,header=True)
    sz = data.shape
    N_pixel_x = sz[1]
    N_pixel_y = sz[0]
    print("number of cols: %f"%N_pixel_x)
    print("number of rows: %f"%N_pixel_y)
    x_axis = np.linspace(0,N_pixel_x-1,num=N_pixel_x,dtype=int)
    y_axis = np.linspace(0,N_pixel_y-1,num=N_pixel_y,dtype=int)
    pixel_x = np.zeros((N_pixel_y,N_pixel_x),dtype=int)
    for row in range(N_pixel_y):
        for col in range(N_pixel_x):
            pixel_x[row,col] = col

    pixel_y = np.zeros((N_pixel_y,N_pixel_x),dtype=int)
    for row in range(N_pixel_y):
        for col in range(N_pixel_x):
            pixel_y[row,col] = row


    cutout_slit_arrays_full = []
    new_y_edges = np.empty_like(y_edges)

    new_top_ys = np.linspace(50,data.shape[0],num=y_edges.shape[0],dtype=int)
    for row in range(y_edges.shape[0]):
        new_y_edges[row].fill(new_top_ys[row])
    #print(new_y_edges)

    data = data[:,20:]
    for slit in range(len(y_edges)):
        im = np.copy(data)
        im[:,:] = np.nan
        top_y = y_edges[slit]
        bottom_y = top_y-slit_height


        #in the new image array, I want to shift the columns of included data
        #so they all align in the same range of rows. That way the slit cutouts
        #won't be all slanted.

        new_top_y = new_y_edges[slit]
        new_bot_y = new_top_y-slit_height




        for col in range(im.shape[1]-20):
            col = int(col)
            im[int(new_bot_y[col]):int(new_top_y[col]),col] = \
                   data[int(bottom_y[col+20]):int(top_y[col+20]),col+20]

        cutout_slit_arrays_full.append(im)



    individual_slit_cutouts = []


    slit1 = cutout_slit_arrays_full[-1]
    #print(slit1.shape)
    xi = 0
    xf = slit1.shape[1]-1

    slit_rects = []
    slit_cutout_paths = []
    slit_num = 1
    for slit in cutout_slit_arrays_full:

        in_slit = []
        for row in range(slit.shape[0]):
            if np.isnan(slit[row]).all()==True:
                pass
            else:
                in_slit.append(row)

        yi = in_slit[0]
        yf = in_slit[-1]
        #print(yi,yf)
        slit_rects.append(slit[yi:yf,xi:xf])
        hdu = fits.PrimaryHDU(slit[yi:yf,xi:xf].astype("f"))
        hdu.header = header.copy()
        outname = "slit%s_%s"%(str(slit_num),str(input_fpath[-1]))
        #print(outname)

        slit_fits_path = "%s/slit_cutouts"%input_fpath[:-1]
        #print(slit_fits_path)
        if not os.path.exists(slit_fits_path): os.mkdir(slit_fits_path)

        this_slit_dir = "%s/slit_%s/"%(slit_fits_path,str(slit_num))
        if not os.path.exists(this_slit_dir): os.mkdir(this_slit_dir)

        this_slit = this_slit_dir+outname
        hdu.writeto(this_slit,overwrite=True)


        slit_cutout_paths.append(this_slit)
        slit_num += 1


    return np.asarray(slit_cutout_paths)
