import numpy as np
from define import *
from pyfftw import *

firsttimeflag = 0

def fft_init(nthreads):
    if (firsttimeflag!=0):
        firsttimeflag = 1

    return



def  fft_end(nthreads):

    if (firsttimeflag!=0):
        firsttimeflag = 0
        fftwf_cleanup()

    return



def fft_conv(data1, fdata2, width, height):
    npix = width*height
    npix2 = (((width>>1) + 1)<< 1) * height

    plan = fftwf_plan_dft_r2c_2d(height, width, data1, fdata1, FFTW_ESTIMATE|FFTW_DESTROY_INPUT)

    fftwf_execute(plan)

    fftwf_destroy_plan(plan)

    fac = 1.0/npix
    fdata1p = fdata1
    fdata2p = fdata2
    index = 0
    for i in range(npix2/2):
        real = fdata1p[index] *fdata2p[index] - fdata1p[index+1]*fdata2p[index+1]
        imag = fdata1p[index+1]*fdata2p[index] + fdata1p[index]*fdata2p[index+1]
        fdata1p[index] = fac*real
        index+=1
        fdata1p[index] = fac*imag
        index+=1

    plan = fftwf_plan_dft_c2r_2d(height, width, fdata1,
            data1, FFTW_ESTIMATE|FFTW_DESTROY_INPUT)

    fftwf_execute(plan)

    fftwf_destroy_plan(plan);

    return



def fft_rtf(data, width, height):

    npix2 = (((width>>1) + 1)<< 1) * height

    plan = fftwf_plan_dft_r2c_2d(height, width, data,
            fdata, FFTW_ESTIMATE)

    fftwf_execute(plan)

    fftwf_destroy_plan(plan)

    return fdata


def fft_ctf(data, width, height, sign):

    plan = fftwf_plan_dft_2d(height, width, data,
            data, sign, FFTW_ESTIMATE)

    fftwf_execute(plan);

    fftwf_destroy_plan(plan)

    return



def fft_shift(data, width, height):
    npix = width*height
    datat = data
    xc = width/2
    yc = height/2
    index =0
    for y in range(height):
        y2 = y-yc
        if (y2<0):
            y2 += height
        
        y2 *= width
        for x in range(width):
            x2 = x-xc
            if (x2<0):
                x2 += width
            
            temp[x2+y2] = datat[index]
            index+=1

    tempt = temp
    for i in range(npix):
        data[i] = tempt[i]
    
    return
