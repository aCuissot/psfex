import cosmoHammer
from numpy.linalg import svd
from math import sqrt
import numpy as np
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
from scipy.stats.stats import kurtosis
import cv2

filename = 'your_image_file.fits'
image = fits.getdata(filename)
# TODO: create S the matrix containing all the star patches nb_stars x 15*15 
u, singleValues = svd(S)
singleValues = singleValues[:12] 

radius = np.ones((15,15))
for x in range(15):
    for y in range(15):
        radius[x][y] = sqrt((x-7.5)**2 + (y-7.5)**2)

print(radius)
beta1 = 1.233
beta2 = 3.419
gamma = 0.387
alpha1=1/sqrt(2**(1/(beta1-1))-1)
alpha2=1/sqrt(2**(1/(beta2-1))-1)

def moffat1(r):
    return (1+(r/alpha1)**2)**(-beta1)

def moffat2(r):
    return (1+(r/alpha2)**2)**(-beta2)

def moffat():
    I1 = np.copy(radius)
    I2 = np.copy(radius)
    m1vec = np.vectorize(moffat1(r))
    m2vec = np.vectorize(moffat2(r))
    I1 = m1vec(I1)
    I2 = m2vec(I2)
    return I1 + I2



def profiles(image):
    ypix, xpix = np.where(image==1)
    x = np.take(image, ypix[0], axis=0)
    y = np.take(image, xpix[0], axis=1)

    return x, y #these are the horizontal and vertical profiles through the star's centroid

def interpolate_width(axis):
    half_max = 1/2
    x = np.linspace(0, len(axis), len(axis))

    # Do the interpolation
    spline = UnivariateSpline(x, axis-half_max, s=0)
    r1, r2 = spline.roots()

    return r2-r1 #this is the FWHM along the specified axis

horizontal, vertical = profiles(image)
fwhm_x = interpolate_width(horizontal)
fwhm_y = interpolate_width(vertical)
F = (fwhm_x + fwhm_y)/2
model = moffat()

"""SIZE"""

r1 = F * sqrt(2**(1/(beta1-1))-1)/(2*sqrt(2**(1/beta1)-1))
r2 = F * sqrt(2**(1/(beta2-1))-1)/(2*sqrt(2**(1/beta2)-1))

"""ELLIPTICITY FLEXION KURTOSIS"""

for x in range(15):
    for y in range(15):
        teta1, teta2 = x, y
        teta1*=r1
        teta1*=r2
        e1, e2
        e=sqrt(e1**2+e2**2)
        if e>0:
            A = 1/sqrt(2) * np.array([[(e2/abs(e2))*sqrt((1+e)*(1+e1/e)), -sqrt((1-e)*(1-e1/e))][sqrt((1+e)*(1-e1/e)), (e2/abs(e2))*sqrt((1-e)*(1+e1/e))]])
        else:
            A = np.ones((2,2))
        f1, f2, g1, g2
        D1 = -1/2 * (np.array([[3*f1, f2], [f2, f1]]) +  np.array([[g1, g2],[g2, -g1]]))
        D2 = -1/2 * (np.array([[f2, f1], [f1, 3*f2]]) +  np.array([[g2, -g1],[-g1, -g2]]))

E = kurtosis(image)        

"""CENTROID"""

gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

ret,thresh = cv2.threshold(gray_image,127,255,0)

M = cv2.moments(thresh)

cX = int(M["m10"] / M["m00"])
cY = int(M["m01"] / M["m00"])

"""
Applying all this transforation to the profil should give the postion of the photon on the simulated image.
"""
