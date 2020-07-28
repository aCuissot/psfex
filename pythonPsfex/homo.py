import numpy as np
from define import *
from vignet import *

def psf_homo(psf, filename, homopsf_params, homobasis_number, homobasis_scale, ext, next):

    npix = psf.size[0]*psf.size[1]
    poly = psf.poly
    ndim = poly.ndim
    coeff = poly.basis
    ncoeff = poly.ncoeff

    nbasis = psf_pshapelet(basis, psf.size[0],psf.size[1], homobasis_number, sqrt(homobasis_number+1.0)*homobasis_scale)
    nfree = nbasis*ncoeff

    moffat.xc[0] = (psf.size[0]/2)
    moffat.xc[1] = (psf.size[1]/2)
    moffat.amplitude = 1.0
    moffat.fwhm_min = moffat.fwhm_max = homopsf_params[0]
    moffat.theta = 0.0
    moffat.beta = homopsf_params[1]
    moffat.nsubpix = 1
    psf_moffat(psf, moffat)
    target = None
    psf.loc = target

    bigsize[0] = psf.size[0]*2
    bigsize[1] = psf.size[1]*2
    nbigpix = bigsize[0]*bigsize[1]

    fft_init(prefs.nthreads)
    f1 = 0
    for i in range(nbasis):
        vignet_copy(basis[i*npix], psf.size[0],psf.size[1], bigbasis, bigsize[0],bigsize[1], 0,0, VIGNET_CPY)
        fft_shift(bigbasis, bigsize[0], bigsize[1])
        fbigbasis = fft_rtf(bigbasis, bigsize[0], bigsize[1])
        for c in range(ncoeff):
            f1+=1
            vignet_copy(psf.comp[c*npix], psf.size[0],psf.size[1], bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY)
            fft_conv(bigconv, fbigbasis, bigsize[0], bigsize[1])
            vignet_copy(bigconv, bigsize[0],bigsize[1], basisc[f1*npix], psf.size[0],psf.size[1], 0,0, VIGNET_CPY)
            basis2 = basisc
            b2_index = 0
            for f2 in range(f1+1):
                basis1 = basisc
                b1_index = f1*npix
                dval = 0.0
                for p in range(npix):
                    dval += basis1[b1_index]*basis2[b2_index]
                    b2_index+=1
                    b1_index+=1
   
                cross[f1+f2*nfree] = cross[f2+f1*nfree] = dval
            
            b1_index = f1*npix
            basis2 = target
            dval = 0.0
            for p in range(npix):
                dval += basis1[b1_index]*basis2[p]
                b1_index+=1
            
            tcross[f1] = dval
            
    fft_end(prefs.nthreads)

    nt = HOMO_NSNAP**ndim
   
    dstep = 1.0/HOMO_NSNAP
    dstart = (1.0-dstep)/2.0
    for d in range(ndim):
        dpos[d] = -dstart
    
    for n in range(nt):
        poly_func(poly, dpos)
        amatt = amat
        amat_index = 0
        bmatt = bmat
        bmat_index
        for j1 in range(nbasis):
            for c1 in range(ncoeff):
                for j2 in range(nbasis):
                    for c2 in range(ncoeff):
                        dval = 0.0
                        for c3 in range(ncoeff):
                            f1 = j1*ncoeff + c3
                            for c4 in range(ncoeff):
                                dval += coeff[c3]*coeff[c4]*cross[f1*nfree+j2*ncoeff + c4]
                        amatt[amat_index] += coeff[c1]*coeff[c2]*dval
                        amat_index+=1

                dval = 0.0
                for c3 in range(ncoeff):
                    dval += coeff[c3]*tcross[j1*ncoeff+c3]
                
                bmatt[bmat_index] += coeff[c1]*dval
                bmat_index+=1
    
        for d in range(ndim):
            if (dpos[d]<dstart-0.01) :
                dpos[d] += dstep
                break
            else:
                dpos[d] = -dstart

    bmatt = bmat
    bmat_index = 0
    for j in range(nbasis):
        kernelt = kernel
        kernel_index = 0
        for c in range(ncoeff):
            basis1 = basis
            dval = bmatt[bmat_index]
            bmat_index +=1 
            for p in range(npix):
                kernelt[kernel_index]+= dval*(basis1[j*npix + p])
                kernel_index+=1

    vignet_copy(kernel, psf.size[0],psf.size[1],
                bigconv, bigsize[0],bigsize[1], 0,0, VIGNET_CPY)
    fft_conv(bigconv, fbigpsf, bigsize[0], bigsize[1]);
    vignet_copy(bigconv, bigsize[0],bigsize[1],
                kernorm, psf.size[0],psf.size[1], 0,0, VIGNET_CPY);
    
    dval = 0.0
    kernelt = kernel
    for p in range(npix):
        dval += kernelt[p]
    
    if (dval>0.0):
        a = 1.0/dval
        kernelt = kernel
        for p in range(npix*ncoeff):
            kernelt[p] *= a

    psf.homo_kernel = kernel
    psf.homopsf_params[0] = homopsf_params[0]
    psf.homopsf_params[1] = homopsf_params[1]
    psf.homobasis_number = homobasis_number
    psf_savehomo(psf, filename, ext, next)

    return



def psf_savehomo(psf, filename, ext, next):
    if (not ext):
        cat = new_cat(1)
        init_cat(cat)
        cat.filename = filename
        if (open_cat(cat, WRITE_ONLY) != RETURN_OK):
            error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename)
        
        if (next>1):
            save_tab(cat, cat.tab)
        
    poly = psf.poly
    tab = new_tab("HOMO_DATA")
    addkeywordto_head(tab, "POLNAXIS", "Number of context parameters")
    fitswrite(tab.headbuf, "POLNAXIS", poly.ndim, H_INT, T_LONG)
    for i in range(poly.ndim):
        str = "POLGRP%1d" % (i+1)
        addkeywordto_head(tab, str, "Polynom group for this context parameter")
        temp = poly.group[i]+1
        fitswrite(tab.headbuf, str, temp, H_INT, T_LONG)
        str = "POLNAME%1d" % (i+1)
        addkeywordto_head(tab, str, "Name of this context parameter");
        fitswrite(tab.headbuf, str, psf.contextname[i], H_STRING, T_STRING)
        str = "POLZERO%1d" % (i+1)
        addkeywordto_head(tab, str, "Offset value for this context parameter")
        fitswrite(tab.headbuf, str, psf.contextoffset[i], H_EXPO, T_DOUBLE)
        str = "POLSCAL%1d" % (i+1)
        addkeywordto_head(tab, str, "Scale value for this context parameter")
        fitswrite(tab.headbuf, str, psf.contextscale[i], H_EXPO, T_DOUBLE)

    addkeywordto_head(tab, "POLNGRP", "Number of context groups")
    fitswrite(tab.headbuf, "POLNGRP", poly.ngroup, H_INT, T_LONG)
    for i in range(poly.ngroup):
        str = "POLDEG%1d" % (i+1)
        addkeywordto_head(tab, str, "Polynom degree for this context group")
        fitswrite(tab.headbuf, str, poly.degree[i], H_INT, T_LONG)
    
    addkeywordto_head(tab, "PSF_FWHM", "FWHM of target PSF")
    fitswrite(tab.headbuf, "PSF_FWHM", psf.homopsf_params[0], H_FLOAT,T_DOUBLE)
    addkeywordto_head(tab, "PSF_BETA", "Moffat Beta of target PSF")
    fitswrite(tab.headbuf, "PSF_BETA", psf.homopsf_params[1],
            H_FLOAT,T_DOUBLE)
    addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data")
    fitswrite(tab.headbuf, "PSF_SAMP", psf.pixstep, H_FLOAT, T_FLOAT)
    tab.bitpix = BP_FLOAT
    tab.bytepix = t_size[T_FLOAT]
    if (poly.ncoeff>1):
        tab.naxis = 3
        tab.naxisn[0] = psf.size[0]
        tab.naxisn[1] = psf.size[1]
        tab.naxisn[2] = poly.ncoeff
        tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]*tab.naxisn[2]
    else:
        tab.naxis = 2
        tab.naxisn[0] = psf.size[0]
        tab.naxisn[1] = psf.size[1]
        tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
    
    tab.bodybuf = psf.homo_kernel
    if (next == 1):
        prim_head(tab)
    
    fitswrite(tab.headbuf, "XTENSION", "IMAGE   ", H_STRING, T_STRING)

    save_tab(cat, tab)
    tab.bodybuf = None

    if (ext==next-1):
        free_cat(cat , 1)

    return