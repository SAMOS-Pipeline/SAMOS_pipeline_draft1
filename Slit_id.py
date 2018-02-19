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

"""

This module is responsible for identifying the edges of the slits.  It creates a region file that is plotted 
on top of the master_flat to show the regions picked out by the function.  Based on Flame pipeline.
This still needs work for identifying slits that are close to each other.

"""

input = 'LMask1/LMask1master_flat.fits'
data,header = fits.getdata(input,header=True)

def get_edges(data_in,approx_edge=1414,cutout_size=30,binsize=20,starting_pixel=0):
    """
    approx_edge is approximate y-pixel for top of the slit.
    cutout size is number of y rows of pixels to check for variation.
    binsize is number of x columns to include each step.
    algorithms adapted for python based on Flame data reduction pipeline written in IDL (Belli, Contursi, and Davies (2017))
    """
    
    sz = data.shape
    N_pixel_x = sz[1] #number of x pixels
    x_edge = []
    y_edge = []

    previous_ycoord = approx_edge


    while starting_pixel < N_pixel_x:

        end_pixel = np.asarray([starting_pixel + binsize-1,N_pixel_x-1]).min()
        cutout_bin = data[np.int(previous_ycoord) - np.int(cutout_size/2): np.int(previous_ycoord) + np.int(cutout_size/2),\
                        np.int(starting_pixel) : np.int(end_pixel)]
                                          
        profile = np.median(cutout_bin, axis=1)
    
        derivative = np.roll(profile,1)-profile
        derivative[0] = 0
        derivative[-1] = 0
        peak,peak_location = derivative.max(),np.argmax(derivative)               
                        
                        
        peak_location += (previous_ycoord - (cutout_size/2))                  
                        
        x_edge.append(int(np.round(0.5*(starting_pixel + end_pixel),0)))
        y_edge.append(int(peak_location))    
    
        previous_ycoord = peak_location
    
        starting_pixel += binsize
    


    nan_arr = np.empty(N_pixel_x)
    nan_arr[:] = np.nan
    y_edge_full = nan_arr.copy()
    for i in range(len(y_edge_full)):
        for j in range(len(x_edge)):
            if i==x_edge[j]:
                y_edge_full[x_edge[j]]=y_edge[j]
            else:
                pass
                
    return x_edge,y_edge


"""
data,header = fits.getdata('LMask1/flat_fielded/fnLMask18150c1.fits',header=True)

N_pix_x = data.shape[1]

num_slits = 1
space_top = 84.0 #empty space between top slit and edge of detector in pixels
space_between_slits = 2.0
bar_height = 22.0
index_first_bar = 1
index_last_bar = 2

slit_top = 1500 - space_top - bar_height * (index_first_bar-1.0) - 0.5*space_between_slits
slit_bottom = 1500 - space_top - bar_height*(index_last_bar) + 0.5*space_between_slits

target_position = 0.5*(slit_bottom+slit_top)
"""



region_filename = 'LMask1/slits.reg'
hdu = fits.open('LMask1/LMask1master_flat.fits')

xt,yt = get_edges(hdu[0].data)
xb,yb = get_edges(hdu[0].data,approx_edge=1396,cutout_size=4,binsize=10)

reg_string_top = """
# Region file format: DS9 version 4.1  
# Filename: LMask1/LMask1master_flat.fits
global color=green width=1 font=helvetica 10 normal select=1 highlite=1 fixed=0 edit=1 move=1 delete=1 include=1 source=1
image; line """+"(" + str(xt[0]) +','+str(yt[0])+','+str(xt[int(len(xt)/2)])+ \
    ','+str(yt[int(len(yt)/2)])+")"+""" # color=cyan background
image; line """+"(" + str(xt[int(len(xt)/2)]) +','+str(yt[int(len(yt)/2)])+ \
    ','+str(xt[-1])+','+str(yt[-1])+")"+""" # color=cyan background 
image; line """+"(" + str(xb[0]) +','+str(yb[0])+','+str(xb[int(len(xb)/2)])+ \
    ','+str(yb[int(len(yb)/2)])+")"+""" # color=pink background
image; line  """+"(" + str(xb[int(len(xb)/2)]) +','+str(yb[int(len(yb)/2)])+ \
    ','+str(xb[-1])+','+str(yb[-1])+")"+""" # color=pink background  
    """


r = regions.write_ds9(reg_string_top,region_filename)

fig = plt.figure()
ax = plt.subplot(111)
ax.imshow(hdu[0].data, cmap=cm.gray,origin="lower")

r2 = pyregion.parse(reg_string_top).as_imagecoord(header=hdu[0].header)


patch_list, artist_list = r2.get_mpl_patches_texts()


for p in patch_list:
    ax.add_patch(p)
for t in artist_list:
    ax.add_artist(t)


ax.imshow(hdu[0].data, cmap=cm.gray,origin="lower")


plt.show()




'''
hdu = fits.PrimaryHDU(out_data.astype("f"))
hdu.header = header.copy()
hdu.header.add_history('Flat Fielded')
hdu.writeto(output,overwrite=True)
'''








