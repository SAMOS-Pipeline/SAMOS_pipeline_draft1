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



def get_edges(input,slit_txt,cutout_size=30,binsize=20):
    """
    approx_edge is approximate y-pixel for top of the slit.
    cutout size is number of y rows of pixels to check for variation.
    binsize is number of x columns to include each step.
    algorithms adapted for python based on Flame data reduction pipeline written in IDL (Belli, Contursi, and Davies (2017))
    """
    print(slit_txt)
    
    
    # read in file of approximate slit edges
    approx_edges = np.genfromtxt(slit_txt)
    #print(approx_edges)
    data,header = fits.getdata(input,header=True)
    sz = data.shape
    N_pixel_x = sz[1] #number of x pixels
   #print(sz)
    #print(data[2249-15:2249+15,0:19])
    x_edges_main = []
    y_edges_main = []

    for approx_edge in approx_edges:
    
        starting_pixel = 0    
        x_edge = []
        y_edge = []
        previous_ycoord = approx_edge
        #print(previous_ycoord)


        while starting_pixel < N_pixel_x:
            

            end_pixel = np.asarray([starting_pixel + binsize-1,N_pixel_x-1]).min()
            cutout_bin = data[np.int(previous_ycoord) - np.int(cutout_size/2): np.int(previous_ycoord) + np.int(cutout_size/2),\
                        np.int(starting_pixel) : np.int(end_pixel)]
            #print(cutout_bin)
                                          
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
    
        x_edges_main.append(x_edge)
        y_edges_main.append(y_edge)
        
    #print(x_edges_main)
        
    #fill in all other values with nans.  Not sure why yet but leaving it here for now.
    nan_arr = np.empty(N_pixel_x)
    nan_arr[:] = np.nan
    y_edge_full = nan_arr.copy()
    for i in range(len(y_edge_full)):
        for j in range(len(x_edge)):
            if i==x_edge[j]:
                y_edge_full[x_edge[j]]=y_edge[j]
            else:
                pass
                


    #hdu = fits.open('LMask1/flat_fielded/fnLMask18150c1.fits')


    str_tops = []
    str_bottoms = []
        
    for i in range(len(x_edges_main)):
    
        
        xt,yt = x_edges_main[i],y_edges_main[i]
        xti,yti = int(xt[0]),int(yt[0])
        xtf,ytf = int(xt[-1]),int(yt[-1])
        xbi,ybi = xti,yti-21
        xbf,ybf = xtf,ytf-21
        str_top ='''line('''+str(xti)+","+str(yti)+","+str(xtf)+","+str(ytf)+''') # line=0 0 color=cyan\n'''
        str_bottom = '''line('''+str(xbi)+","+str(ybi)+","+str(xbf)+","+str(ybf)+''') # line=0 0 color=yellow\n'''
        
        str_tops.append(str_top)
        str_bottoms.append(str_bottom)
           
        
        
    region_txt = open('LMask2/reg_txt.reg','w')
    region_txt.write( """# Region file format: DS9 version 4.1"""+"\n"+"""# File name: LMask2/LMask2master_flat.fits"""+"\n"+\
    """global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 """+\
    """fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical \n""")
    
    for i in range(len(str_tops)):
        region_txt.write(str(str_tops[i])+str(str_bottoms[i]))
    
    
    region_txt.close()
    
    reg_string = open("LMask2/reg_txt.reg","r").read()

    #r = regions.write_ds9(reg_string,region_filename)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.imshow(data, cmap=cm.gray,origin="lower")

    r2 = pyregion.parse(reg_string).as_imagecoord(header=header)


    patch_list, artist_list = r2.get_mpl_patches_texts()


    for p in patch_list:
        ax.add_patch(p)
    for t in artist_list:
        ax.add_artist(t)

    logscale_data = np.log(1000*data + 1)/np.log(1000)
    ax.imshow(logscale_data, cmap=cm.gray,origin="lower")


    plt.show()
    
    
    return



'''
hdu = fits.PrimaryHDU(out_data.astype("f"))
hdu.header = header.copy()
hdu.header.add_history('One slit extracted')
hdu.writeto('LMask1/single_slit.fits',overwrite=True)
'''


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


