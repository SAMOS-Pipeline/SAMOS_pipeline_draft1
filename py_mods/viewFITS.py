from __future__ import print_function
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import os

class viewFITS:

    def __init__(self):
        self.info = None
        self.header = None
        self.data = None
        self.fig = None


    def view_FITS(self,fname):

        hdu_list = fits.open(fname)
        self.info = hdu_list.info()
        self.header = hdu_list[0].header
        self.data = hdu_list[0].data
        title = os.path.split(fname)[1]

        self.fig, ax = plt.subplots(figsize=(10,10))
        img = ax.pcolormesh(self.data,cmap='gray',norm=colors.LogNorm());
        ax.set_xlabel('dispersion direction',fontsize=20)
        ax.set_ylabel('spatial direction',fontsize=20)
        ax.set_title(title,fontsize=20)
        cb = self.fig.colorbar(img,ax=ax,orientation='vertical',fraction=0.054, pad=0.04,format=ticker.ScalarFormatter())
        tick_locator = ticker.MaxNLocator(nbins=6)
        cb.locator = tick_locator
        cb.update_ticks()
        cb.set_label('log pixel value',fontsize=15)

        plt.close()
        hdu_list.close()

        return self
