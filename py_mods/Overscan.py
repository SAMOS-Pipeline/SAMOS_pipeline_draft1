
import glob
from SAMOSHelpers import *

# Must read the data array •
# Must convert to float32 •
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
    print("Working on %s" %(input))
    (xb0,xb1),(yb0,yb1) = np.subtract(ParseSec(h["biassec"]),1)
    (xd0,xd1),(yd0,yd1) = np.subtract(ParseSec(h["datasec"]),1)

    data = d[1300:2800]

    print("Writing %s" % (output))
    hdu = fits.PrimaryHDU(data.astype("f"))
    hdu.header = h.copy()
    hdu.header["bitpix"] = -32
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

    print("Writing %s" % (output))
    hdu = fits.PrimaryHDU(data.astype("f"))
    hdu.header = h.copy()
    hdu.header["bitpix"] = -32
    [hdu.header.remove(key) for key in ["bscale","bzero"]]
    hdu.header.add_history("overscan subtracted")
    hdu.writeto(output,overwrite=True)
    print("%s written." % (output))

    p = os.path.dirname(output)
    pj = "%s/jpeg" % (p)
    if not os.path.exists(pj): os.mkdir(pj)
    MakeThumbnail(output,pj)


if __name__ == "__main__":
   input,output = sys.argv[1:3]
   Overscan(input,output)
