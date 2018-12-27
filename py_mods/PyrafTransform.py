from __future__ import print_function
import glob
from astropy.io import fits
from SAMOSHelpers import *
import sys
import os
from pyraf import iraf
import iraf_steps as irf_stp
import iraf_params as irf_par
import shutil
from matplotlib import pyplot as plt
import matplotlib.cm as cm
from astropy import wcs
from numina.array import combine
from numina.array import combine_shape
from numina.array import subarray_match


def plot_transf_imgs(targs):

        #targs = [i for i in targs if '28158' in os.path.basename(i).split('.')[1]]
        filepaths = []
        names = []
        root_dir = os.path.split(targs[0])[0]
        for i in targs:
            fn = os.path.basename(i)
            if len(fn.split('.')[0])<6:
                fn = fn[:4]+'0'+fn[4:]
            if fn not in names:
                names.append(fn)
        for slit in sorted(names):
            filepaths.append(root_dir+'/'+slit)
        #print(filepaths)

        lambdas = []
        cdelts = []
        first_col_positions = []
        cdata = []
        slit_data = []
        for targ in targs:
            hdulist = fits.open(targ)
            cdata.append(hdulist)
            data,header = fits.getdata(targ,header=True)
            data = data[:][:2010]

            for row in range(data.shape[0]):
                #data[row] = data[row][:2010]
                for col in range(data.shape[1]):
                    pix = data[row][col]
                    if pix<0.1:
                        #print(np.where(data[row]>0))
                        pix_replace = np.median(data[row][np.where(data[row]>=0.1)])
                        #data[row][col] = pix_replace
                        pix = pix_replace

                    data[row][col] = pix

            #logscale_data = np.log(1000*data + 1)/np.log(1000)
            #slit_data.append(logscale_data)
            slit_data.append(data)


            """
            lambda1 = header['CRVAL1']
            cdelt = np.round(header['CDELT1'],2)
            cdelts.append(cdelt)
            lambda2 = lambda1+cdelt*data.shape[1]
            #first_col_pos = first_col_positions.append()#(7467.9773-lambda_min)*cdelt)
            wavelengths = np.arange(lambda1,lambda2,cdelt)#np.linspace(lambda1,lambda2,data.shape[1])
            first_col_positions.append(int((7467.9773-7222.1523)*cdelt))
            lambdas.append(wavelengths)
            #ax.imshow(logscale_data, cmap=cm.gray,origin="lower")
            #plt.show(fig)

            #ax.imshow(logscale_data, cmap=cm.gray,origin="lower")
            """
        cdata = np.asarray(cdata)
        slit_data = np.asarray(slit_data)
        baseheader = header
        basedata = cdata[0][0].data #fits.getdata(targs[0],header=False)
        print(slit_data.shape)
        baseshape = cdata[0][0].data.shape
        subpixshape = cdata[0][0].data.shape
        basearr = np.asarray([baseshape],dtype='int')
        print(basearr,basearr.shape)
        refpix = np.divide(basearr, 2).astype('float')
        print(refpix.shape)
        offsets_xy = np.zeros((len(targs),refpix.shape[1]))
        with fits.open(targs[0]) as hdulist:
             wcsh = wcs.WCS(hdulist[0].header)
             skyref = wcsh.wcs_pix2world(refpix, 1)

        for idx, frame in enumerate(targs[1:]):
            with fits.open(frame) as hdulist:
                wcsh = wcs.WCS(hdulist[0].header)
                pixval = wcsh.wcs_world2pix(skyref,1)
                offsets_xy[idx+1] = -(pixval[0]-refpix[0])

        offsets_fc = offsets_xy[:,::-1]
        print('offsets_fc=',offsets_fc)
        offsets_fc_t = np.round(offsets_fc).astype('int')
        offsets_fc = offsets_xy[:, ::-1]
        print('offsets_fc_t= %s'%(str(offsets_fc_t)))
        finalshape, offsetsp = combine_shape(subpixshape, offsets_fc_t)
        print('final shape = %s, \n offset shape = %s'%(str(finalshape),str(offsetsp)))
        rhduls, regions = resize_hdulists(cdata, subpixshape, offsetsp, finalshape)
        regions = np.asarray(regions)
        print(regions.shape)
        method = combine.mean
        rect_arr = []
        for d in rhduls:
            rect_arr.append(d[0].data)
        rect_arr = np.asarray(rect_arr)
        #rect_arr = np.asarray([d[0].data for d in rhduls], dtype='float32') #method([d[0].data for d in rhduls], dtype='float32')


        arrmap = np.zeros(finalshape, dtype='int')

        result = []
        #masks = np.asarray([(arrmap[region] > 0) for region in regions])
        #print(masks)

        #data2 = method([d[0].data for d in cdata], masks=masks, dtype='float32')
        for spec in range(slit_data.shape[0]):
            arrmap = np.zeros(finalshape, dtype='int')
            mask = regions[spec][1]
            start = mask.start
            end = mask.stop
            print(start)
            print(mask)
            #for row in range(slit_data.shape[1]):

            arrmap[:,start:end] = slit_data[spec,:,:]
            result.extend(arrmap)
            #mask_start = regions[spec][1][0]
            #mask_end = regions[spec][1][1]
            #print(mask_start,mask_end)

        result = np.asarray(result)
        print(result.shape)
        """
        slit_data = np.asarray(slit_data)
        lambdas = np.asarray(lambdas)
        lambda_min = lambdas.min()
        lambda_max = lambdas.max()
        cdelts = np.asarray(cdelts)
        #slit_len = slit_data.shape[1]
        #ax.imshow(logscale_data,cmap=cm.gray,origin="lower")


        #full_lambda_row_labels = np.asarray([int(i) for i in full_lambda_row])
        #rect_row = np.full(full_lambda_row.shape,np.nan)
        rect_arr = []
        print(cdelts)
        full_lambda_row = np.arange(lambda_min,lambda_max,cdelt)
        first_col_pos = (7467.9773-lambda_min)*cdelt
        for img in range(lambdas.shape[0]):

            rect_mask = np.where(np.logical_and(full_lambda_row>=lambdas[img].min(),full_lambda_row<=lambdas[img].max()))
            #print(rect_mask)
            #print(len(lambdas[img]))
            for row in range(slit_data[img].shape[0]):
                #print(slit_data[img].shape)
                rect_row = np.full(full_lambda_row.shape,np.nan)
                diff = len(rect_row[rect_mask])-len(slit_data[img][row])

                #print(len(rect_row[rect_mask]))

                rect_row[rect_mask] = slit_data[img][row]#[:-diff]

                rect_arr.append(rect_row)


        fig = plt.figure(figsize=(15,15))
        ax = plt.subplot(111)
        ax.imshow(np.asarray(rect_arr),cmap=cm.gray,origin='lower')
        xticks = ax.get_xticks()
        full_lambda_row_labels = np.full(len(xticks),lambda_min)+(xticks*cdelt)
        full_lambda_row_labels = np.asarray([int(i) for i in full_lambda_row])
        ax.set_xticklabels(full_lambda_row_labels)
        ax.set_ylabel('pixel position in y-direction')
        ax.set_xlabel('wavelength (angstroms)')
        plt.show(fig)
        rect_arr = np.asarray(rect_arr)
        print(rect_arr.shape)
        print(rect_arr)
        """
        #ax.imshow(rect_arr,cmap=cm.gray,origin="lower")
        #plt.show()
        hdu = fits.PrimaryHDU(result.astype("f"), header=baseheader)
        #hdu = fits.PrimaryHDU(rect_arr.astype("f"))
        hdu.header = header.copy()
        hdu.writeto('combine_test.fits',overwrite=True)

def resize_hdulists(hdulists, shape, offsetsp, finalshape, window=None):
         rhdulist = []
         regions = []
         for hdulist, rel_offset in zip(hdulists, offsetsp):
             region, _ = subarray_match(finalshape, rel_offset, shape)
             rframe = resize_hdul(hdulist, finalshape, region)
             rhdulist.append(rframe)
             regions.append(region)
         return rhdulist, regions

def resize_hdul(hdul, newshape, region, extensions=None, window=None,
                     scale=1, fill=0.0, conserve=True):
     from numina.frame import resize_hdu

     if extensions is None:
         extensions = [0]

     nhdul = [None] * len(hdul)
     for ext, hdu in enumerate(hdul):
         if ext in extensions:
              nhdul[ext] = resize_hdu(hdu, newshape,
                                     region, fill=fill,
                                     window=window,
                                     scale=scale,
                                     conserve=conserve)
         else:
             nhdul[ext] = hdu
     return fits.HDUList(nhdul)

def pyraf_trnsf(fuel):

    dbdir = fuel.identify.db
    intargs = fuel.util.science.corr_files
    transf_dir = '%s/%s'%(fuel.input.slit_mask,'transformed_targs')
    if not os.path.exists(transf_dir): os.mkdir(transf_dir)
    wvcal_targs = []
    for fc in glob.glob('%s/fcslit*'%(dbdir)):
        slitn = os.path.basename(fc).split('.')[0][2:]
        for targ in intargs:
            targslitn,targnm = os.path.basename(targ).split('_')
            if targslitn==slitn:
                print(targslitn,slitn)
                fcnm = os.path.basename(fc)[2:]
                if os.path.exists('current_targ.fits'): os.system('rm current_targ.fits')
                os.symlink(targ,'current_targ.fits')

                outtargnm = '%s.t%s'%(targslitn,targnm)
                irf_stp.transform('current_targ.fits',outtargnm,fcnm)

                wvcal_targ = '%s/%s'%(transf_dir,outtargnm)
                os.rename(outtargnm,wvcal_targ)
                wvcal_targs.append(wvcal_targ)

                os.rename('current_targ.fits','previous_targ.fits')
                os.unlink('%s/%s'%(os.getcwd(),'previous_targ.fits'))

    fuel.util.science.corr_files = wvcal_targs

    return fuel
