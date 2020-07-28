import numpy as np
from define import *

def pca_onsnaps(psfs, ncat, npc):

    ndim = psfs[0].poly.ndim
    w = psfs[0].size[0]
    h = psfs[0].size[1]
    npix = w*h
    nt = 1
    for d in range(ndim):
        nt *= PCA_NSNAP
    
    dstep = 1.0/PCA_NSNAP
    dstart = (1.0-dstep)/2.0

    for c in range(ncat):
        psf = psfs[c]
        for d in range(ndim):
            dpos[d] = -dstart
        
        for n in range(nt):
            psf_build(psf, dpos)
            pix1 = pix = psf.loc
            covmatt = covmat
            covmatt_index = 0
            str = "Setting-up the PCA covariance matrix (%.0f%%)..." % (100.0*(n/nt+c)/ncat)
            for j in range(npix):
                pix2 = pix;
                dval = pix1[j]
                for i in range(npix):
                    covmatt[covmatt_index] += dval*pix2[i]
                    covmatt_index+=1
                
            for d in range(ndim):
                if (dpos[d]<dstart-0.01):
                    dpos[d] += dstep
                    break
                else: 
                    dpos[d] = -dstart
          
    for p in range(npc):
        str = "Computing Principal Component vector #%d..." % (p)
        pca_findpc(covmat, basis[p*npix], npix);

    return basis


def pca_oncomps(psfs, next, ncat, npc):
    
    ndim = psfs[0].poly.ndim
    npix = psfs[0].size[0]*psfs[0].size[1]
    nt = 1
    for d in range(ndim):
        nt *= PCA_NSNAP
    
    npixt = npix*nt*next
    dstep = 1.0/PCA_NSNAP
    dstart = (1.0-dstep)/2.0

    compt = comp
    compt_index=  0
    for c in range(ncat):
        str = "Setting-up the PCA covariance matrix (%.0f%%)..." % (50.0*c/ncat)
        for e in range(next):
            psf = psfs[c*next+e]
            for d in range(ndim):
                dpos[d] = -dstart
            
            for n in range(nt):
                psf_build(psf, dpos)
                pix = psf.loc
                for i in range(npix):
                    compt[compt_index] = pix[i]
                    compt_index+=1
                
                for d in range(ndim):
                    if (dpos[d]<dstart-0.01):
                        dpos[d] += dstep
                        break
                    else :
                        dpos[d] = -dstart
                    
    for i in range(npixt):
        dval = 0
        compt = comp
        compt_index = i
        for c in range(ncat):
            compt_index += npixt
            dval += compt[compt_index]
        
        dval /= ncat
        compt_index = i
        for c in range(ncat):
            compt_index += npixt
            compt[compt_index] -= dval
        
    covmatt = covmat
    covmatt_index = 0
    for c1 in range(ncat):
        str = "Setting-up the PCA covariance matrix (%.0f%%)..." % (50.0+50.0*(c1)/ncat)
        for c2 in range(ncat):
            comp1 = comp + c1*npixt
            comp2 = comp + c2*npixt
            dval = 0.0
            for i in range(npixt):
                dval += comp1[i]*comp2[i]
            
            covmatt[covmatt_index] = dval;
            covmatt_index+=1
        
    for p in range(npc):
        str = "Computing Principal Component vector #%d..." % (p)
        pca_findpc(covmat, vector, ncat)
        for c in range(ncat):
            pc[c*npc+p] = vector[c]
        
    return pc


def pca_findpc(covmat, vec, nmat):
    t=tmat
    for i in range(nmat):
        t[i] = 1.0

    dtnorm = 1.0
    for n in range(PCA_NITER) :
        if not dtnorm>PCA_CONVEPS:
            break
        xnorm = 0.0
        c=covmat
        c_index = 0
        x=xmat
        for j in range(nmat):
            xval=0.0
            t=tmat
            for i in range(nmat):
                xval += c[c_index]*t[i]
                c_index+=1
            
            xnorm += xval*xval
            x[j] = xval
        
        xnorm = 1.0/sqrt(xnorm)
        dtnorm = 0.0
        t=tmat
        x=xmat
        for i in range(nmat):
            dtval = t[i]
            t[i] = x[i]*xnorm
            dtval -= t[i]
            dtnorm += dtval*dtval
        
        dtnorm = sqrt(dtnorm)
    
    lambda_var = 0.0
    c=covmat
    c_index = 0
    x=tmat
    for j in range(nmat):
        xval=0.0
        t=tmat
        for i in range(nmat):
            xval += c[c_index]*t[i]
            c_index+=1
        
        lambda_var += xval*x[j]
    
    c=covmat
    c_index=0
    x=tmat
    for j in range(nmat):
        xval=x[j]*lambda_var
        t=tmat
        for i in range(nmat):
            c[c_index] -= xval*t[i]
            c_index+=1
    v=vec
    t=tmat
    for i in range(nmat):
        v[i] = float(t[i])

    return lambda_var