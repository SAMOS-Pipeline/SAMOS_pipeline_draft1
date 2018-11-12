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

"""

This module is responsible for identifying the edges of the slits.  It creates a region file that is plotted
on top of the master_flat to show the regions picked out by the function.  Based on Flame pipeline.
This still needs work for identifying slits that are close to each other.

"""



def get_edges(input,slit_txt,mask,cutout_size=10,binsize=20):
    """
    approx_edge is approximate y-pixel for top of the slit.
    cutout size is number of y rows of pixels to check for variation.
    binsize is number of x columns to include each step.
    algorithms adapted for python based on Flame data reduction pipeline written in IDL (Belli, Contursi, and Davies (2017))
    """
    #print(slit_txt)


    # read in file of approximate slit edges
    approx_edges = np.genfromtxt(slit_txt,unpack=True,usecols=0)

    #print(approx_edges)
    data,header = fits.getdata(input,header=True)
    sz = data.shape
    N_pixel_x = sz[1] #number of x pixels
    N_pixel_y = sz[0]

   #print(sz)
    #print(data[2249-15:2249+15,0:19])
    x_edges_main = []
    y_edges_main = []

    for approx_edge in approx_edges:
        #try to put starting pixel in the slit a little more to remove overlapping cutout_slit
        approx_edge = approx_edge-3
        starting_pixel = int(N_pixel_x/2.)
        x_edge = []
        y_edge = []
        previous_ycoord = approx_edge
        #print(previous_ycoord)


        while starting_pixel < N_pixel_x:

            end_pixel = np.asarray([starting_pixel + binsize-1,N_pixel_x-1]).min()
            cutout_bin = data[np.int(previous_ycoord-cutout_size/2.):
                        np.int(previous_ycoord+cutout_size/2.),
                            starting_pixel:end_pixel]

            profile = np.median(cutout_bin, axis=1)

            derivative = np.roll(profile,1)-profile
            derivative[0] = 0
            derivative[-1] = 0
            peak,peak_location = derivative.max(),np.argmax(derivative)

            #print("previous ycoord - cutout:previous ycoord =[%s:%s]"\
            #%(np.int(previous_ycoord) - np.int(cutout_size),np.int(previous_ycoord)))
            #print("starting pixel %s"%(starting_pixel))
            #print("end pixel %s"%(end_pixel))
            #print("peak location: %s"%(peak_location))


            peak_location += (previous_ycoord - (cutout_size/2.))
            if peak_location < 0:
                peak_location = 0
            if peak_location > N_pixel_y-1:
                peak_location = N_pixel_y-1

            #print("peak location updated: %s"%(peak_location))

            x_edge.append(int(np.round(0.5*(starting_pixel + end_pixel),0)))
            y_edge.append(int(peak_location))
            #print(previous_ycoord)

            previous_ycoord = peak_location

            starting_pixel += binsize

        starting_pixel = int(N_pixel_x/2.)
        previous_ycoord = approx_edge


        while starting_pixel > 0:

            end_pixel = np.asarray([starting_pixel-binsize-1,15]).max()
            cutout_bin = data[int(previous_ycoord-cutout_size/2.):
                            int(previous_ycoord+cutout_size/2.),
                            end_pixel:starting_pixel]

            profile = np.median(cutout_bin, axis=1)
            derivative = np.roll(profile,1)-profile
            derivative[0] = 0
            derivative[-1] = 0
            peak,peak_location = derivative.max(),np.argmax(derivative)

            peak_location += (previous_ycoord - (cutout_size/2.))
            if peak_location < 0:
                peak_location = 0
            if peak_location > N_pixel_y-1:
                peak_location = N_pixel_y-1


            x_edge.insert(0,int(0.5*(starting_pixel+end_pixel)))
            y_edge.insert(0,int(peak_location))

            previous_ycoord = peak_location
            starting_pixel -= binsize


        #x_edges_main.insert(0,x_edge)
        #y_edges_main.insert(0,y_edge)
        x_edge_full = np.linspace(0,N_pixel_x-1,num=N_pixel_x,dtype=int)
        nan_arr = np.empty(N_pixel_x)
        nan_arr[:] = np.nan
        y_edge_full = nan_arr.copy()
        for i in range(len(y_edge_full)):
            for j in range(len(x_edge)-1):
                if i==x_edge[j]:
                    y_edge_full[0:x_edge[0]]=y_edge[0]
                    y_edge_full[x_edge[-1]:]=y_edge[-1]
                    y_edge_full[x_edge[j]:x_edge[j+1]]=y_edge[j]
                else:
                    pass
        x_edges_main.append(x_edge_full)
        y_edges_main.append(y_edge_full)



    #hdu = fits.open('LMask1/flat_fielded/fnLMask18150c1.fits')

    #print(data[x_edges_main,y_edges_main])

    x_edges_main, y_edges_main = np.asarray(x_edges_main),np.asarray(y_edges_main)
    #print(x_edges_main.shape,y_edges_main.shape)
    #print(y_edges_main)
    return x_edges_main,y_edges_main

"""
    slope,yint = np.polyfit(x_edges_main[0],y_edges_main[0],1)
    print(slope,yint)


    y_top = np.dot(np.asarray((slope*x) + yint).astype(int),np.ones(len(y_edges_main[0])))
    y_bottom = np.dot(np.asarray((slope*x)+ yint).astype(int)-21,np.ones(len(y_edges_main[0])))




    fig = plt.figure()
    ax = plt.subplot(111)
    #ax.imshow(data, cmap=cm.gray,origin="lower")
    logscale_data = np.log(1000*data + 1)/np.log(1000)
    #ax.imshow(logscale_data, cmap=cm.gray,origin="lower")
    ax.plot(x,y)

    print(x,y)

    ax.imshow(logscale_data[y[-1]-30:y[0]+30,x[0]-2:x[-1]+2],cmap=cm.gray,origin="lower")



    plt.show()

    #print(x_edges_main.shape)
    #print(y_edges_main.shape)



    polygons = []
    for i in range(x_edges_main.shape[0]):

        numx = x_edges_main.shape[1]
        xt,yt = x_edges_main[i],y_edges_main[i]
        xti,yti = int(xt[0]),int(yt[0])
        xtmid,ytmid = int(xt[int(numx/2.)]),int(yt[int(numx/2.)])
        xtf,ytf = int(xt[-1]),int(yt[-1])
        xbi,ybi = xti,yti-21
        xbmid,ybmid = xtmid,ytmid-21
        xbf,ybf = xtf,ytf-21
        box = '''polygon('''+str(xti)+","+str(yti)+","+str(xtf)+","+str(ytf)+","+str(xbf)+","+str(ybf)+","+str(xbi)+","+str(ybi)+''') \
        # line=0 0 color=blue\n'''


        polygons.append(box)



    region_txt = open('%s/reg_txt.reg'%mask,'w')
    if 'LMask1' in mask:
        region_txt.write( '''# Region file format: DS9 version 4.1'''+"\n"+'''# File name: LMask1/LMask1master_flat.fits'''+"\n"+\
        '''global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '''+\
        '''fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical \n''')
        print(mask)
    elif 'LMask2' in mask:
        region_txt.write( '''# Region file format: DS9 version 4.1'''+"\n"+'''# File name: LMask2/LMask2master_flat.fits'''+"\n"+\
        '''global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 '''+\
        '''fixed=0 edit=1 move=1 delete=1 include=1 source=1\n physical \n''')

    for i in range(len(polygons)):
        #region_txt.write(str(str_tops[i])+str(str_bottoms[i]))
        region_txt.write(str(polygons[i]))


    region_txt.close()

    reg_string = open("%s/reg_txt.reg"%mask,"r").read()


    r = regions.write_ds9(reg_string,region_filename)
"""
