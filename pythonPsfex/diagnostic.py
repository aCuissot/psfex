import numpy as np
from define import *

lm_opts = [1.0e-2, 1.0e-12, 1.0e-12, 1.0e-12, 1.0e-4]

def psf_wcsdiagnostic(psf, wcs):
    raw = np.zeros(NAXIS, dtype=np.float64)

    if (not psf.samples_accepted or not wcs or wcs.naxis<2):
        psf.moffat_fwhm_wcs = psf.moffat_fwhm
        psf.moffat_fwhm_wcs_min = psf.moffat_fwhm_min
        psf.moffat_fwhm_wcs_max = psf.moffat_fwhm_max
        psf.pfmoffat_fwhm_wcs = psf.pfmoffat_fwhm
        psf.pfmoffat_fwhm_wcs_min = psf.pfmoffat_fwhm_min
        psf.pfmoffat_fwhm_wcs_max = psf.pfmoffat_fwhm_max
        psf.pixscale_wcs = psf.pixscale_wcs_min = psf.pixscale_wcs_max = 0.0
        return
    
    naxis = wcs.naxis
    ncx = ncy = nt = 1
    if psf.nsnap>1:
        nsnap2 = psf.nsnap 
    else:
        nsnap2 = 2
    for n in range(psf.poly.ndim):
        nt *= nsnap2
        if (psf.cx>=0 and n<psf.cx):
            ncx *= nsnap2
        if (psf.cy>=0 and n<psf.cy):
            ncy *= nsnap2
        
    for i in range(naxis):
        raw[i] = wcs.naxisn[i]/2.0 + 0.5
    
    xstep = wcs.naxisn[0] / (nsnap2-1)
    ystep = wcs.naxisn[1] / (nsnap2-1)
    raw[1] = 0.5
    fwhmmin = pffwhmmin = pixscalemin = BIG
    fwhmmax = pffwhmmax = pixscalemax = -BIG
    mfwhm = mpffwhm = mpixscale = 0.0
    for j in range(nsnap2):
        raw[0] = 0.5
        for i in range(nsnap2):
            pixscale = sqrt(wcs_scale(wcs, raw))
            dval = dpfval = 0.0
            nfwhm = 0
            for n in range(nt):
                if ((psf.cx<0 or (n/ncx)%nsnap2 == i) and (psf.cy<0 or (n/ncy)%nsnap2 == j)):
                    if psf.nsnap>1:
                        n2=n
                    else:
                        n2=0
                    dval += pixscale* sqrt(psf.moffat[n2].fwhm_min*psf.moffat[n2].fwhm_max)
                    dpfval += pixscale * sqrt(psf.pfmoffat[n2].fwhm_min*psf.pfmoffat[n2].fwhm_max)
                    nfwhm+=1
                
            
            fwhm = dval / nfwhm
            mfwhm += fwhm
            if (fwhm < fwhmmin) :
                fwhmmin = fwhm
            if (fwhm > fwhmmax):
                fwhmmax = fwhm

            pffwhm = dpfval / nfwhm
            mpffwhm += pffwhm
            if (pffwhm < pffwhmmin):
                pffwhmmin = pffwhm
            if (pffwhm > pffwhmmax):
                pffwhmmax = pffwhm

            mpixscale += pixscale
            if (pixscale < pixscalemin):
                pixscalemin = pixscale

            if (pixscale > pixscalemax):
                pixscalemax = pixscale

            raw[0] += xstep
        raw[1] += ystep
    

    psf.moffat_fwhm_wcs = mfwhm/nsnap2/nsnap2*DEG/ARCSEC
    psf.moffat_fwhm_wcs_min = fwhmmin*DEG/ARCSEC
    psf.moffat_fwhm_wcs_max = fwhmmax*DEG/ARCSEC

    psf.pfmoffat_fwhm_wcs = mpffwhm/nsnap2/nsnap2*DEG/ARCSEC
    psf.pfmoffat_fwhm_wcs_min = pffwhmmin*DEG/ARCSEC
    psf.pfmoffat_fwhm_wcs_max = pffwhmmax*DEG/ARCSEC

    psf.pixscale_wcs = mpixscale/nsnap2/nsnap2*DEG/ARCSEC
    psf.pixscale_wcs_min = pixscalemin*DEG/ARCSEC
    psf.pixscale_wcs_max = pixscalemax*DEG/ARCSEC
    return


def psf_diagnostic(psf):
    dpos = np.zeros(POLY_MAXDIM, dtype=np.float64)
    nmed = 0
    npc = psf.poly.ndim
    nt = 1
    for i in range(npc):
        nt *= psf.nsna
        nmed += nmed*psf.nsnap + (psf.nsnap-1)/2
    
    psf.nmed = nmed

    w = psf.size[0]
    h = psf.size[1]
    dstep = 1.0/psf.nsnap
    dstart = (1.0-dstep)/2.0

    memset(dpos, 0, POLY_MAXDIM*sizeof(double))
    for i in rangde(npc):
        dpos[i] = -dstart

    fwhm = psf.fwhm / psf.pixstep
    moffat_parammin[0] = 1.0/(psf.fwhm*psf.fwhm)/10.0
    moffat_parammax[0] = 1.0/(psf.fwhm*psf.fwhm)*10.0
    moffat_parammin[1] = 0.0
    moffat_parammax[1] = w - 1.0
    moffat_parammin[2] = 0.0
    moffat_parammax[2] = h - 1.0
    moffat_parammin[3] = fwhm/3.0
    moffat_parammax[3] = fwhm*3.0
    moffat_parammin[4] = fwhm/3.0
    moffat_parammax[4] = fwhm*3.0
    moffat_parammin[5] = moffat_parammax[5] = 90.0
    moffat_parammin[6] = PSF_BETAMIN
    moffat_parammax[6] = 10.0

    for n in range(nt):
        psf_compdiag(psf, moffat[n], dpos, 1)
        psf_compdiag(psf, pfmoffat[n], dpos, PSF_NSUBPIX)

        for i in range(npc):
            if (dpos[i]<dstart-0.01):
                dpos[i] += dstep
                break
            else:
                dpos[i] = -dstart

    psf.pfmoffat_fwhm_min = psf.pfmoffat_ellipticity_min = psf.pfmoffat_ellipticity1_min = psf.pfmoffat_ellipticity2_min = psf.pfmoffat_beta_min= psf.pfmoffat_residuals_min= psf.moffat_fwhm_min = psf.moffat_ellipticity_min= psf.moffat_ellipticity1_min = psf.moffat_ellipticity2_min= psf.moffat_beta_min= psf.moffat_residuals_min = psf.sym_residuals_min= psf.noiseqarea_min = BIG
    psf.pfmoffat_fwhm_max = psf.pfmoffat_ellipticity_max = psf.pfmoffat_ellipticity1_max  = psf.pfmoffat_ellipticity2_max  = psf.pfmoffat_beta_max = psf.pfmoffat_residuals_max = psf.moffat_fwhm_max = psf.moffat_ellipticity_max = psf.moffat_ellipticity1_max = psf.moffat_ellipticity2_max = psf.moffat_beta_max  = psf.moffat_residuals_max = psf.sym_residuals_max= psf.noiseqarea_max = -BIG

    for n in range(nt):
        fwhm=0.5*moffat[n].fwhm_min+moffat[n].fwhm_max
        if (fwhm < psf.moffat_fwhm_min):
            psf.moffat_fwhm_min = fwhm
        if (fwhm > psf.moffat_fwhm_max):
            psf.moffat_fwhm_max = fwhm
        
        ellip = moffat[n].fwhm_max + moffat[n].fwhm_min
        if (ellip > 0.0) :
            ellip = (moffat[n].fwhm_max - moffat[n].fwhm_min) / ellip

        if (ellip < psf.moffat_ellipticity_min):
            psf.moffat_ellipticity_min = ellip
        
        if (ellip > psf.moffat_ellipticity_max):
            psf.moffat_ellipticity_max = ellip
        

        ellip1 = ellip*cosf(2.0*moffat[n].theta*DEG)
        ellip2 = ellip*sinf(2.0*moffat[n].theta*DEG)

        if (ellip1 < psf.moffat_ellipticity1_min):
            psf.moffat_ellipticity1_min = ellip1
        
        if (ellip1 > psf.moffat_ellipticity1_max):
            psf.moffat_ellipticity1_max = ellip1
        
        if (ellip2 < psf.moffat_ellipticity2_min):
            psf.moffat_ellipticity2_min = ellip2
        
        if (ellip2 > psf.moffat_ellipticity2_max):
            psf.moffat_ellipticity2_max = ellip2
        

        if (moffat[n].beta < psf.moffat_beta_min):
            psf.moffat_beta_min = moffat[n].beta
        
        if (moffat[n].beta > psf.moffat_beta_max):
            psf.moffat_beta_max = moffat[n].beta
        

        if (moffat[n].residuals < psf.moffat_residuals_min):
            psf.moffat_residuals_min = moffat[n].residuals
        
        if (moffat[n].residuals > psf.moffat_residuals_max):
            psf.moffat_residuals_max = moffat[n].residuals


        if (moffat[n].symresiduals < psf.sym_residuals_min):
            psf.sym_residuals_min = moffat[n].symresiduals
        
        if (moffat[n].symresiduals > psf.sym_residuals_max):
            psf.sym_residuals_max = moffat[n].symresiduals
        
        if (moffat[n].noiseqarea < psf.noiseqarea_min):
            psf.noiseqarea_min = moffat[n].noiseqarea
        
        if (moffat[n].noiseqarea > psf.noiseqarea_max):
            psf.noiseqarea_max = moffat[n].noiseqarea
        
        fwhm=0.5*(pfmoffat[n].fwhm_min+pfmoffat[n].fwhm_max)
        if (fwhm< psf.pfmoffat_fwhm_min):
            psf.pfmoffat_fwhm_min = fwhm
        
        if (fwhm > psf.pfmoffat_fwhm_max):
            psf.pfmoffat_fwhm_max = fwhm
        
        ellip = pfmoffat[n].fwhm_max + pfmoffat[n].fwhm_min
        if (ellip > 0.0) :
            ellip = (pfmoffat[n].fwhm_max - pfmoffat[n].fwhm_min) / ellip
        

        ellip1 = ellip*cosf(2.0*pfmoffat[n].theta*DEG)
        ellip2 = ellip*sinf(2.0*pfmoffat[n].theta*DEG)

        if (ellip < psf.pfmoffat_ellipticity_min):
            psf.pfmoffat_ellipticity_min = ellip
        
        if (ellip > psf.pfmoffat_ellipticity_max):
            psf.pfmoffat_ellipticity_max = ellip
        
        if (ellip1 < psf.pfmoffat_ellipticity1_min):
            psf.pfmoffat_ellipticity1_min = ellip1
        
        if (ellip1 > psf.pfmoffat_ellipticity1_max):
            psf.pfmoffat_ellipticity1_max = ellip1
        
        if (ellip2 < psf.pfmoffat_ellipticity2_min):
            psf.pfmoffat_ellipticity2_min = ellip2
        
        if (ellip2 > psf.pfmoffat_ellipticity2_max):
            psf.pfmoffat_ellipticity2_max = ellip2
        

        if (pfmoffat[n].beta < psf.pfmoffat_beta_min):
            psf.pfmoffat_beta_min = pfmoffat[n].beta
        
        if (pfmoffat[n].beta > psf.pfmoffat_beta_max):
            psf.pfmoffat_beta_max = pfmoffat[n].beta
        

        if (pfmoffat[n].residuals < psf.pfmoffat_residuals_min):
            psf.pfmoffat_residuals_min = pfmoffat[n].residuals
        
        if (pfmoffat[n].residuals > psf.pfmoffat_residuals_max):
            psf.pfmoffat_residuals_max = pfmoffat[n].residuals

    psf.moffat_fwhm = 0.5*(moffat[nmed].fwhm_min + moffat[nmed].fwhm_max)

    ellip = moffat[nmed].fwhm_max + moffat[nmed].fwhm_min
    if ellip > 0.0:
        psf.moffat_ellipticity = (moffat[nmed].fwhm_max - moffat[nmed].fwhm_min) / ellip 
    else:
        psf.moffat_ellipticity = 0.0
    psf.moffat_ellipticity1 = psf.moffat_ellipticity * cosf(2.0*moffat[nmed].theta*DEG)
    psf.moffat_ellipticity2 = psf.moffat_ellipticity * sinf(2.0*moffat[nmed].theta*DEG)

    psf.moffat_beta = moffat[nmed].beta
    psf.moffat_residuals = moffat[nmed].residuals
    psf.sym_residuals = moffat[nmed].symresiduals
    psf.noiseqarea = moffat[nmed].noiseqarea

    psf.pfmoffat_fwhm = 0.5*(pfmoffat[nmed].fwhm_min + pfmoffat[nmed].fwhm_max)

    ellip = pfmoffat[nmed].fwhm_max + pfmoffat[nmed].fwhm_min
    
    if (ellip > 0.0):
         psf.pfmoffat_ellipticity = pfmoffat[nmed].fwhm_max - pfmoffat[nmed].fwhm_min / ellip
    else:
         psf.pfmoffat_ellipticity = 0.0
   
    psf.pfmoffat_ellipticity1 = psf.pfmoffat_ellipticity*cosf(2.0*pfmoffat[nmed].theta*DEG)
    psf.pfmoffat_ellipticity2 = psf.pfmoffat_ellipticity*sinf(2.0*pfmoffat[nmed].theta*DEG)

    psf.pfmoffat_beta = pfmoffat[nmed].beta
    psf.pfmoffat_residuals = pfmoffat[nmed].residuals

    psf.moffat = moffat
    psf.pfmoffat = pfmoffat

    return



def psf_compdiag(psf0, moffat, dpos, oversamp):
    
    dparam = np.zeros(PSF_DIAGNPARAM, dtype=np.float64)
    param = np.zeros(PSF_DIAGNPARAM, dtype=np.float32)
    
    psf = psf0
    psf.nsubpix = oversamp
    npc = psf.poly.ndim
    w = psf.size[0]
    h = psf.size[1]
    m = w*h

    if (psf.samples_accepted) :
        psf_build(psf, dpos)
        fwhm = psf.fwhm / psf.pixstep
        param[0] = 1.0/(psf.fwhm*psf.fwhm)
        param[1] = (w-1)/2.0
        param[2] = (h-1)/2.0
        param[3] = fwhm
        param[4] = fwhm
        param[5] = 0.0
        param[6] = 3.0
        psf_boundtounbound(param, dparam)
        niter = dlevmar_dif(psf_diagresi, dparam, NULL, PSF_DIAGNPARAM, m, PSF_DIAGMAXITER, lm_opts, NULL, NULL, NULL, psf)
        psf_unboundtobound(dparam,param)

    moffat.nsubpix = psf.nsubpix
    moffat.amplitude = param[0]/(psf.pixstep*psf.pixstep)
    moffat.xc[0] = param[1]
    moffat.xc[1] = param[2]
    if (param[3] > param[4]):
        moffat.fwhm_max = param[3]*psf.pixstep
        moffat.fwhm_min = param[4]*psf.pixstep
        moffat.theta = (fmod(param[5]+360.0, 180.0))
    else:
        moffat.fwhm_max = param[4]*psf.pixstep
        moffat.fwhm_min = param[3]*psf.pixstep
        moffat.theta = (fmod(param[5]+450.0, 180.0))
    
    if (moffat.theta > 90.0):
        moffat.theta -= 180.0

    moffat.beta = param[6]

    for i in range(npc):
        moffat.context[i] = dpos[i]*psf.contextscale[i]+psf.contextoffset[i]

    moffat.residuals = psf_normresi(param, psf)
    if (oversamp==1):
        moffat.symresiduals = float(psf_symresi(psf))
        moffat.noiseqarea = float(psf_noiseqarea(psf))

    return


def psf_diagresi(dparam, fvec, m, n, adata):
    par = np.zeros(PSF_DIAGNPARAM, dtype=np.float32)

    psf = adata
    nsubpix = psf.nsubpix
    psf_unboundtobound(dparam, par)
    w = psf.size[0]
    h = psf.size[1]
    ct = cosf(par[5]*PI/180.0)
    st = sinf(par[5]*PI/180.0)
    if par[6]>PSF_BETAMIN:
        fac = 4.0*(powf(2.0, 1.0/par[6])-1.0)
    else:
        fac = 4.0*(powf(2.0, 1.0/PSF_BETAMIN) - 1.0)
    if par[3]>PSF_FWHMMIN:
        inva2 = fac/(par[3]*par[3])
    else:
        inva2 = fac/(PSF_FWHMMIN*PSF_FWHMMIN)
    if par[3]>PSF_FWHMMIN:
        invb2 = fac/(par[4]*par[4])
    else:
        invb2 = fac/(PSF_FWHMMIN*PSF_FWHMMIN)
    invb2 = fac/(PSF_FWHMMIN*PSF_FWHMMIN)
    cxx = inva2*ct*ct + invb2*st*st
    cyy = inva2*st*st + invb2*ct*ct
    cxy = 2.0*ct*st*(inva2 - invb2)
    a = par[0] / (nsubpix*nsubpix)
    beta = -par[6]
    dxstep = psf.pixsize[0]/(nsubpix*psf.pixstep)
    dystep = psf.pixsize[1]/(nsubpix*psf.pixstep)
    dy0 = -par[2] - 0.5*(nsubpix - 1.0)*dystep
    loc = psf.loc
    fvect = fvec
    for i in range(w*h):
        fvect[i] = -loc[i]
    
    for yd in range(nsubpix): 
        dy0+=dystep
        dx0 = -par[1] - 0.5*(psf.nsubpix - 1.0)*dxstep
        for xd in range(nsubpix):
            dx0+=dxstep
            fvect = fvec
            fvect_index = 0
            dy = dy0
            for y in range (h):
                dy+=1.0
                dy2 = dy*dy
                dx = dx0
                for x in range (w):
                    dx+=1.0
                    fvect[fvect_index] += (a*powf(1.0+cxx*dx*dx+cyy*dy2+cxy*dx*dy, beta))
                    fvect_index += 1
    
    psf_boundtounbound(par, dparam)
    return


def psf_normresi(par, psf):

    w = psf.size[0]
    h = psf.size[1]
    ct = cosf(par[5]*PI/180.0)
    st = sinf(par[5]*PI/180.0)
    nsubpix = psf.nsubpix
    fvec = np.zeros(w*h, dtype=np.float32)
    if par[6]>PSF_BETAMIN:
        fac = 4.0*(pow(2.0, 1.0/par[6]) - 1.0)
    else:
        fac = 4.0*(pow(2.0, 1.0/PSF_BETAMIN) - 1.0)

    if par[3]>PSF_FWHMMIN:
        inva2 = fac/(par[3]*par[3])
    else:
        inva2 = fac/(PSF_FWHMMIN*PSF_FWHMMIN)

    if par[4]>PSF_FWHMMIN:
        invb2 = fac/(par[4]*par[4])
    else:
        invb2 = fac/(PSF_FWHMMIN*PSF_FWHMMIN)

    cxx = inva2*ct*ct + invb2*st*st
    cyy = inva2*st*st + invb2*ct*ct
    cxy = 2.0*ct*st*(inva2 - invb2)
    a = par[0] / (nsubpix*nsubpix)
    beta = -par[6]
    dxstep = psf.pixsize[0]/(nsubpix*psf.pixstep)
    dystep = psf.pixsize[1]/(nsubpix*psf.pixstep)
    dy0 = -par[2] - 0.5*(nsubpix - 1.0)*dystep
    fvect = fvec
    for i in range(w*h):
        fvect[i] = 0.0
    
    for yd in range(nsubpix):
        dy0+=dystep
        dx0 = -par[1] - 0.5*(nsubpix - 1.0)*dxstep;
        for xd in range(nsubpix):
            dx0+=dxstep
            fvect = fvec
            fvect_index = 0
            dy = dy0;
            for y in range(h):
                dy+=1.0
                dy2 = dy*dy;
                dx = dx0;
                for x in range(w):
                    dx+=1.0
                    fvect[fvect_index] += a*powf(1.0+cxx*dx*dx + cyy*dy2 + cxy*dx*dy, beta)
                    fvect_index+=1

    resi = norm = 0.0
    fvect = fvec
    loc = psf.loc
    for i in range(w*h):
        val = loc[i]+fvect[i];
        resi += val*fabsf(fvect[i] - loc[i])
        norm += val*val

    if norm>0:
        return 2.0 * resi / norm
    return 1.0


def psf_symresi(psf):
    loc = psf.loc
    resi = norm = 0.0
    s = psf.size[0]*psf.size[1]
    for i in range(s):
        val = loc[i]
        valsym = loc[s-i]
        valmean = val + valsym
        norm += valmean*valmean
        resi += valmean * fabs((val - valsym))

    if norm>0:
        return 2.0 * resi / norm
    return 1.0


def   psf_noiseqarea(psf):

    loc = psf.loc
    sum = sum2 = 0.0
    for i in range(psf.size[0]*psf.size[1]):
        val = loc[i]
        sum += val
        sum2 += val*val
    
    if sum2 > 0.0:
        return psf.pixstep*psf.pixstep*sum*sum/sum2
    return 0.0


def psf_moffat(psf, moffat):
    w = psf.size[0]
    h = psf.size[1]
    xc = moffat.xc[0]
    yc = moffat.xc[1]
    nsubpix = moffat.nsubpix
    a = moffat.amplitude*psf.pixstep*psf.pixstep / (nsubpix*nsubpix)
    beta = -moffat.beta
    ct = cos(moffat.theta*PI/180.0)
    st = sin(moffat.theta*PI/180.0)
    fac = 4*(pow(2, -1.0/beta) - 1.0)
    inva2 = fac/(moffat.fwhm_max*moffat.fwhm_max)*psf.pixstep*psf.pixstep
    invb2 = fac/(moffat.fwhm_min*moffat.fwhm_min)*psf.pixstep*psf.pixstep
    cxx = inva2*ct*ct + invb2*st*st
    cyy = inva2*st*st + invb2*ct*ct
    cxy = 2.0*ct*st*(inva2 - invb2)
    dxstep = psf.pixsize[0]/(nsubpix*psf.pixstep)
    dystep = psf.pixsize[1]/(nsubpix*psf.pixstep)
    dy0 = -yc - 0.5*(nsubpix - 1.0)*dystep
    loc = psf.loc
    for i in range(w*h):
        loc[i] = 0.0
    
    for yd in range(nsubpix):
        dy0+=dystep
        dx0 = -xc - 0.5*(nsubpix - 1.0)*dxstep
        for xd in range(nsubpix):
            dx0+=dxstep
            loc = psf.loc
            loc_index = 0
            dy = dy0
            for y in range(h):
                dy+=1.0
                dy2 = dy*dy
                dx = dx0
                for x in range(w):
                    dx+=1.0
                    loc[loc_index] += a*powf(1.0+cxx*dx*dx + cyy*dy2 + cxy*dx*dy, beta)
                    loc_index +=1
                
    return

def psf_boundtounbound(param, dparam):
    for p in range(PSF_DIAGNPARAM):
        if (moffat_parammin[p]!=moffat_parammax[p]):
            num = param[p] - moffat_parammin[p]
            den = moffat_parammax[p] - param[p]
            if num>1e-50:
                if den>1e-50:
                    dparam[p] = log(num/den)
                else:
                    dparam[p] = 50.0
            else:
                dparam[p] = -50.0    
                    
        elif (moffat_parammax[p] > 0.0 or moffat_parammax[p] < 0.0):
            dparam[p] = param[p] / moffat_parammax[p]

    return


def psf_unboundtobound(dparam, param):

    for p in range(PSF_DIAGNPARAM):
        if moffat_parammin[p]!=moffat_parammax[p]:
            if dparam[p]>50.0:
                param[p] = (moffat_parammax[p] - moffat_parammin[p])/ (1.0 + exp(-50.0))+ moffat_parammin[p]
            if dparam[p]<-50.0:
                param[p] = (moffat_parammax[p] - moffat_parammin[p])/ (1.0 + exp(50.0))+ moffat_parammin[p]
            else:
                param[p] = (moffat_parammax[p] - moffat_parammin[p])/ (1.0 + exp(-(dparam[p])))+ moffat_parammin[p]
        else:
                param[p] = dparam[p]*moffat_parammax[p]
