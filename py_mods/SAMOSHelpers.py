
import sys
import os
import numpy as np
from astropy.io import fits
from PIL import Image as P

def MuyMalo(x):
    print(x)
    os._exit(1)

def IsItHere(x):
    if not os.path.exists(x):
       print("Error:  %s not found!"%(x))
       MuyMalo("Quitting")

def MakeThumbnail(input,jdir,s1=-10,s2=10):
    f = fits.open(input)
    d = f[0].data.astype("f")
    m = np.median(d)
    s = 1.49*np.median(np.abs(d-m))
    z1 = m-s1*s
    z2 = m+s1*s
    pixels = np.clip(255*(d-z1)/(z2-z1),0,255).astype("b")
    pixmap = P.frombytes("L",(pixels.shape[1],pixels.shape[0]),pixels.ravel().tostring())
    image = P.merge("L",[pixmap])
    jpeg = "%s/%s" % (jdir,os.path.basename(input.replace(".fits",".jpg")))
    image.save(jpeg,"jpeg")


