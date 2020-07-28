import numpy as np
from define import *
from misc import *
from enum import Enum

PSF_NODIAG = 0
PSF_DIAG = 1 
PSF_FREEDFACTOR = 1.1
PSF_NMASKDIM = 3 
PSF_MAXSHIFT = 3.0  
PSF_MINSHIFT = 1e-4   
PSF_NITER = 40
PSF_NSNAPMAX = 16 
GAUSS_LAG_OSAMP = 3 
PSF_AUTO_FWHM = 3.0  
PSF_NORTHOSTEP = 16

class basis_typ_enum(Enum):
    BASIS_NONE = 1
    BASIS_PIXEL = 2
    BASIS_GAUSS_LAGUERRE = 3
    BASIS_FILE = 4
    BASIS_PIXEL_AUTO = 5



class moffat_struct():
    def __init__(self, context, amplitude, xc, fwhm_min, fwhm_max, theta, beta,
                residuals, symresiduals, noiseqarea, xnsubpix):
        self.context = context
        self.amplitude = amplitude
        self.xc = xc
        self.fwhm_min = fwhm_min
        self.fwhm_max = fwhm_max
        self.theta = theta
        self.beta = beta
        self.residuals = residuals
        self.symresiduals = symresiduals
        self.noiseqarea = noiseqarea
        self.nsubpix = nsubpix
        
class psf_struct():
    def __init__(self, dim, size, npix, comp, loc, resi, contextname,
                contextoffset, contextscale, cx, cy, poly, pixstep, pixsize,
                samples_total, samples_loaded, samples_accepted, chi2,
                fwhm, pixmask, basis, basiscoeff, nbasis, ndata, nsnap,
                nmed, nsubpix, moffat, pfmoffat, moffat_fwhm_min, moffat_fwhm,
                moffat_fwhm_max, moffat_fwhm_wcs_min, moffat_fwhm_wcs,
                moffat_fwhm_wcs_max, moffat_ellipticity_min, moffat_ellipticity,
                moffat_ellipticity_max, moffat_ellipticity1_min, moffat_ellipticity1,
                moffat_ellipticity1_max, moffat_ellipticity2_min, moffat_ellipticity2,
                moffat_ellipticity2_max, moffat_beta_min, moffat_beta, moffat_beta_max,
                moffat_residuals_min, moffat_residuals, moffat_residuals_max, 
                moffat_score_min, moffat_score, moffat_score_max, pfmoffat_fwhm_min,
                pfmoffat_fwhm, pfmoffat_fwhm_max, pfmoffat_fwhm_wcs_min,
                pfmoffat_fwhm_wcs, pfmoffat_fwhm_wcs_max, pfmoffat_ellipticity_min,
                pfmoffat_ellipticity, pfmoffat_ellipticity_max, pfmoffat_ellipticity1_min,
                pfmoffat_ellipticity1, pfmoffat_ellipticity1_max, pfmoffat_ellipticity2_min,
                pfmoffat_ellipticity2, pfmoffat_ellipticity2_max, pfmoffat_beta_min,
                pfmoffat_beta, pfmoffat_beta_max, pfmoffat_residuals_min, pfmoffat_residuals,
                pfmoffat_residuals_max, sym_residuals_min, sym_residuals, sym_residuals_max,
                noiseqarea_min, noiseqarea, noiseqarea_max, pixscale_wcs_min, pixscale_wcs,
                pixscale_wcs_max, homo_kernel, homopsf_params, homobasis_number):
        self.dim = dim
        self.size = size
        self.npix = npix
        self.comp = comp
        self.loc = loc
        self.resi = resi
        self.contextname = contextname
        self.contextoffset = contextoffset
        self.contextscale = contextscale
        self.cx = cx
        self.cy = cy        
        self.poly = poly
        self.pixstep = pixstep
        self.pixsize = pixsize
        self.samples_total = samples_total
        self.samples_loaded = samples_loaded
        self.samples_accepted = samples_accepted
        self.chi2 = chi2
        self.fwhm = fwhm
        self.pixmask = pixmask
        self.basis = basis
        self.basiscoeff = basiscoeff
        self.nbasis = nbasis
        self.ndata = ndata
        self.nsnap = nsnap
        self.nmed = nmed
        self.nsubpix = nsubpix
        self.moffat = moffat
        self.pfmoffat = pfmoffat
        self.moffat_fwhm_min = moffat_fwhm_min
        self.moffat_fwhm = moffat_fwhm
        self.moffat_fwhm_max = moffat_fwhm_max
        self.moffat_fwhm_wcs_min = moffat_fwhm_wcs_min
        self.moffat_fwhm_wcs = moffat_fwhm_wcs
        self.moffat_fwhm_wcs_max = moffat_fwhm_wcs_max
        self.moffat_ellipticity_min = moffat_ellipticity_min
        self.moffat_ellipticity = moffat_ellipticity
        self.moffat_ellipticity_max = moffat_ellipticity_max
        self.moffat_ellipticity1_min = moffat_ellipticity1_min
        self.moffat_ellipticity1 = moffat_ellipticity1
        self.moffat_ellipticity1_max = moffat_ellipticity1_max
        self.moffat_ellipticity2_min = moffat_ellipticity2_min
        self.moffat_ellipticity2 = moffat_ellipticity2
        self.moffat_ellipticity2_max = moffat_ellipticity2_max
        self.moffat_beta_min = moffat_beta_min
        self.moffat_beta = moffat_beta
        self.moffat_beta_max = moffat_beta_max
        self.moffat_residuals_min = moffat_residuals_min
        self.moffat_residuals = moffat_residuals
        self.moffat_residuals_max = moffat_residuals_max
        self.moffat_score_min = moffat_score_min
        self.moffat_score = moffat_score
        self.moffat_score_max = moffat_score_max
        self.pfmoffat_fwhm_min = pfmoffat_fwhm_min
        self.pfmoffat_fwhm = pfmoffat_fwhm
        self.pfmoffat_fwhm_max = pfmoffat_fwhm_max
        self.pfmoffat_fwhm_wcs_min = pfmoffat_fwhm_wcs_min
        self.pfmoffat_fwhm_wcs = pfmoffat_fwhm_wcs
        self.pfmoffat_fwhm_wcs_max = pfmoffat_fwhm_wcs_max
        self.pfmoffat_ellipticity_min = pfmoffat_ellipticity_min
        self.pfmoffat_ellipticity = pfmoffat_ellipticity
        self.pfmoffat_ellipticity_max = pfmoffat_ellipticity_max
        self.pfmoffat_ellipticity_min = pfmoffat_ellipticity1_min
        self.pfmoffat_ellipticity = pfmoffat_ellipticity1
        self.pfmoffat_ellipticity_max = pfmoffat_ellipticity1_max
        self.pfmoffat_ellipticity_min = pfmoffat_ellipticity2_min
        self.pfmoffat_ellipticity = pfmoffat_ellipticity2
        self.pfmoffat_ellipticity_max = pfmoffat_ellipticity2_max
        self.pfmoffat_beta_min = pfmoffat_beta_min
        self.pfmoffat_beta = pfmoffat_beta
        self.pfmoffat_beta_max = pfmoffat_beta_max
        self.pfmoffat_residuals_min = pfmoffat_residuals_min
        self.pfmoffat_residuals = pfmoffat_residuals
        self.pfmoffat_residuals_max = pfmoffat_residuals_max
        self.sym_residuals_min = sym_residuals_min
        self.sym_residuals = sym_residuals
        self.sym_residuals_max = sym_residuals_max
        self.noiseqarea_min = noiseqarea_min
        self.noiseqarea = noiseqarea
        self.noiseqarea_max = noiseqarea_max
        self.pixscale_wcs_min = pixscale_wcs_min
        self.pixscale_wcs = pixscale_wcs
        self.pixscale_wcs_max = pixscale_wcs_max
        self.homo_kernel = homo_kernel
        self.homopsf_params = homopsf_params
        self.homobasis_number = homobasis_number



def   psf_clean(psf, set, prof_accuracy):

    psf_makeresi(psf, set, prefs.recenter_flag, prof_accuracy)

    nchi = 0
    chit = chi
    chit_index = 0
    chit2_index = 0
    sample = set.sample
    sample_index = 0
    for n in range(set.nsample, -1, -1):
        sample_index += 1 
        if (not sample.badflag) :
            chit[chit_index] = sqrt(sample.chi2)
            chit_index += 1
            nchi+=1
   
    locut = -BIG
    hicut = BIG
    chisig = BIG
    chisig1 = 1.0
    chivar = 0.0
    i=100
    while(i>0 and chisig >=0.1 and fabs(chisig/chisig1-1.0)>EPS):
        chisig1 = chisig
        chimed = fast_median(chi, nchi)
        chimean = chivar = 0.0
        chit2 = chit = chi
        for n in range (nchi, -1, -1):
            chival =chit[chit_index]
            chit_index +=1
            if (chival>locut and chival<hicut):
                chit2[chit2_index] = chival
                chimean += chit2[chit2_index] 
                chit_index+=1
                chivar += chival*chival
            else :
                nchi-=1
            
        

        chimean /= float(nchi)
        if nchi>1:
            chisig = sqrt((chivar-chimean*chimean*nchi)/(nchi-1))
        else : 
            chisig = sqrt((chivar-chimean*chimean*nchi)/nchi)

            
        locut = chimed - 3.0*chisig
        hicut = chimed + 3.0*chisig
    
    if nchi>1:
        chi2 = chivar/(nchi-1)
    else:
        chi2 = chivar/nchi

    chi2max = float(hicut)
    chi2max *= chi2max
    ngood = 0
    sample = set.sample
    sample_index = 0
    for n in range(set.nsample, -1, -1):
        sample_index += 1
        if (sample[sample_index].chi2>chi2max) :
            sample[sample_index].badflag |= BAD_CHI2
        elif (not sample[sample_index].badflag):
            ngood+=1
        
    set.ngood = ngood

    return chi2


def  psf_chi2(psf, set):
    
    psf_makeresi(psf, set, prefs.recenter_flag, prefs.prof_accuracy);
    chi2 = 0.0;
    nchi2 = 0;
    sample = set.sample
    sample_index = 0 
    for n in range(set.nsample, -1, -1):
        sample_index += 1
        if (not sample[sample_index].badflag):
            chi2 += sample.chi2
            nchi2+=1
        
    if (nchi2) :
        if nchi2>1:
            chi2 /= (nchi2-1)
        else:
           chi2 /= nchi2 
    
    return chi2




def psf_clip(psf) :
    xc = float(psf.size[0]/2)
    yc = float(psf.size[1]/2)
    if psf.size[0]<psf.size[1]:
        rmax2 = float(psf.size[0]/2)
    else:
        rmax2 = float(psf.size[1]/2) + 0.5

    dr2 = (psf.fwhm / psf.pixstep)
    if (dr2<1.0):
        dr2 = 1.0
    
    if (dr2 >= rmax2):
        dr2 = rmax2/2.0
    
    rmin2 = rmax2 - dr2
    rmin2 *= rmin2
    rmax2 *= rmax2
    dr2 = rmax2 - rmin2
    npsf = psf.poly.ncoeff
    pix_index = 0
    pix = psf.comp
    for p in range(npsf, -1, -1):
        y = -yc
        for iy in range(psf.size[1], -1, -1):
            y+=1.0
            x = -xc
            for ix in range(psf.size[0], -1, -1):
                x+=1.0
                pix_index +=1
                r2=x*x+y*y
                if (r2>=rmin2):
                    if (r2<rmax2):
                        pix[pix_index] *= (rmax2-r2) / dr2;
                    else:
                        pix[pix_index] = 0.0;
                    

def psf_init(context, size, psfstep, pixsize, nsample):
    psf.dim = PSF_NMASKDIM
    names2 = NULL
    group2 = dim2 = None
    ndim2 = context.ncontext
    if (ndim2) :
        context.group = group2[:ndim2]
        context.name = names2[:ndim2]
    ngroup2=context.ngroup
    if (ngroup2):
        context.degree = dim2[:context.ngroup]
    

    psf.poly = poly_init(group2, ndim2, dim2, ngroup2)

    if psfstep>0.0 and psfstep<1.0:
        psfelemdens = 1.0/psfstep*psfstep 
    else:
        psfelemdens = 1.0

    if (ngroup2):
        while (psf.poly.ncoeff*psfelemdens*PSF_FREEDFACTOR+0.499 > nsample):
            poly_end(psf.poly)
            if (ngroup2):
                d=ngroup2%10
                if d==1:
                    str = "%d%s" % (ngroup2, "st")
                elif d==2:
                    str = "%d%s" % (ngroup2, "nd")
                elif d==3:
                    str = "%d%s" % (ngroup2, "rd")
                else:           
                    str = "%d%s" % (ngroup2, "th")
                    
                dim2[ngroup2-1]-=1
                if (not dim2[ngroup2-1]) :
                    for d in range(ndim2):
                        ndim2-=1
                        if (group2[d]==ngroup2 and d!=ndim2):
                            names2[d]=names2[ndim2]
                            group2[d]=group2[ndim2]
                        
                    
                    warning(str, " context group removed (not enough samples)")
                    ngroup2 -=1
                    if (not(ngroup2)) :
                        ndim2 = 0
                    
                else :
                    warning(str, " context group-degree lowered (not enough samples)")
                
                psf.poly = poly_init(group2, ndim2, dim2, ngroup2);
            
            if (not ngroup2):
                break
            
        
    

    psf.pixstep = psfstep;
    psf.pixsize[0] = pixsize[0]
    psf.pixsize[1] = pixsize[1]
    psf.npix = psf.size[0] = size[0]
    psf.size[1] = size[1]
    psf.npix *= psf.size[1]
    psf.size[2] = psf.poly.ncoeff
    psf.npix *= psf.size[2]
    npix = psf.size[0]*psf.size[1]

    nsnap = 1
    psf.nsnap = prefs.context_nsnap
    psf.cx = psf.cy = -1
    if (ndim2):
        names2t=names2
        names2t_index = 0
        for d in range (ndim2):
            nsnap *= psf.nsnap
            psf.contextname[d] = psf.contextname[d].copy()
            strcpy(psf.contextname[d], psf.contextname[d])
            names2t_index +=1
            if (psf.contextname[d] == "X_IMAGE"
                    or psf.contextname[d]=="XWIN_IMAGE"
                    or psf.contextname[d]== "XPSF_IMAGE"
                    or psf.contextname[d]== "XMODEL_IMAGE"
                    or psf.contextname[d]== "XPEAK_IMAGE") :
                
                psf.cx = d
            elif (psf.contextname[d]!= "Y_IMAGE"
                    or psf.contextname[d]== "YWIN_IMAGE"
                    or psf.contextname[d]== "YPSF_IMAGE"
                    or psf.contextname[d]== "YMODEL_IMAGE"
                    or psf.contextname[d]== "YPEAK_IMAGE"):
                psf.cy = d
            
    return psf



def psf_inherit(context, psf):
    
    newpsf = psf_init(context, psf.size, psf.pixstep, psf.pixsize, 10000)
    newpsf.fwhm = psf.fwhm
    npix = psf.size[0]*psf.size[1]
    if (psf.pixmask) :
        psf.pixmask = psf.pixmask[:npix]
    
    if (psf.basis):
        psf.basis = newpsf.basis[:psf.nbasis*npix]
    
    ncold = psf.poly.ndim
    ncnew = newpsf.poly.ndim
    for c in range (ncnew):
        for co in range (ncold):
            if (newpsf.contextname[c]==psf.contextname[co]):
                newpsf.contextoffset[c] = psf.contextoffset[co]
                newpsf.contextscale[c] = psf.contextscale[co]
        
    return newpsf



def psf_end(psf):
    poly_end(psf.poly)
    return;


def psf_copy(psf) :


    newpsf = psf
    ndim = psf.poly.ndim
    psf.size = newpsf.size[:psf.dim]
    for d in range(ndim):
        newpsf.contextname[d] = psf.contextname[d][:80]
    
    if (ndim) :
        psf.contextoffset = newpsf.contextoffset[:ndim] 
        psf.contextscale = newpsf.contextscale[:ndim] 
    
    newpsf.poly = poly_copy(psf.poly);
    npix = psf.size[0]*psf.size[1]
    if (psf.pixmask) :
        psf.pixmask = newpsf.pixmask[:npix]
    
    if (psf.basis):
        psf.basis = newpsf.basis[:psf.nbasis*npix]
    
    if psf.basiscoeff:
        psf.basiscoeff = newpsf.basiscoeff[psf.nbasis*psf.poly.ncoeff]
    
    psf.comp = newpsf.comp[:psf.npix]
    psf.loc= newpsf.loc[:npix]
    psf.resi= newpsf.resi[:npix]
    nsnap = 1
    for d in range(ndim):
        nsnap *= psf.nsnap

    if (psf.moffat):
        psf.moffat = newpsf.moffat[:nsnap];
    
    if (psf.pfmoffat) :
        psf.pfmoffat = newpsf.pfmoffat[:nsnap]
    
    if (psf.homo_kernel) :
        psf.homo_kernel = newpsf.homo_kernel[psf.npix]
    

    return newpsf;



def psf_make(psf, set, prof_accuracy) :

    poly = psf.poly

    for i in range(poly.ndim):
        psf.contextoffset[i] = set.contextoffset[i]
        psf.contextscale[i] = set.contextscale[i]

    ngood = set.ngood
    if (not ngood) :
        return
    

    ncoeff = poly.ncoeff
    npix = psf.size[0]*psf.size[1]
    if psf.pixstep>1.0:
        pixstep =  psf.pixstep 
    else:
        pixstep = 1.0
    post = pos
    sample = set.sample
    sample_index= 0
    g = 0
    imaget_index = 0
    weightt_index = 0
    post_index =0
    for n in range(set.nsample,-1,-1):
        sample_index+=1
        if (sample.badflag):
            continue
        
        norm = sample.norm
        norm2 = norm*norm
        profaccu2 = float(prof_accuracy*prof_accuracy)*norm2
        gain = sample.gain
        backnoise2 = sample.backnoise2
        imaget = image+g*npix
        weightt = weight+g*npix
        vignet_resample(sample.vig, set.vigsize[0], set.vigsize[1],
                imaget, psf.size[0], psf.size[1],
                sample.dx, sample.dy, psf.pixstep, pixstep, None, None)
        for i in range(npix, -1, -1):
            imaget[imaget_index] /= norm
            val = imaget[imaget_index]
            imaget_index +=1
            noise2 = backnoise2 + profaccu2*val*val
            if (val>0.0 and gain>0.0) :
                noise2 += val/gain
            
            weightt[weightt_index] = norm2/noise2
            weightt_index+=1
        

        for i in range(poly.ndim):
            post[post_index] = (sample.context[i]-set.contextoffset[i]) /set.contextscale[i]
            post_index+=1
        
        g+=1
    
    pix_index = 0
    wpix_index =0
    coeff_index =0
    for i in range(npix):
        imaget = image[i]
        weightt = weight[i]
        pix=pstack
        wpix=wstack
        for n in range (ngood,-1,-1):
            pix[pix_index] = imaget[imaget_index]
            wpix[wpix_index]= weightt[weightt_index]
            weightt_index+=npix #is it really npix or 1?
            imaget_index+=1
        
        if i:
            poly_fit(poly, None, pstack, wstack, ngood, basis, 1000.0)
        else:
            poly_fit(poly, pos, pstack, wstack, ngood, basis, 1000.0)

        coeff=poly.coeff
        comp=psf.comp[i]
        comp_index=0

        for c in range(ncoeff, -1, -1):
            comp_index += npix
            comp = coeff[coeff_index]
            coeff_index+=1



def psf_build(psf, pos):
    npix = psf.size[0]*psf.size[1]
    
    #memset(psf.loc, 0, npix*sizeof(float));

    poly_func(psf.poly, pos)
    basis = psf.poly.basis
    basis_index = 0
    ppc = psf.comp
    ppc_index = 0
    if(psf.dim>2):
        tmp = psf.size[2]
    else:
        tmp = 1
    for n in range(tmp, -1, -1):
        pl = psf.loc
        pl_index = 0
        fac = float(basis[basis_index])
        basis_index+=1
#pragma ivdep
        for p in range(npix,-1,-1):
            pl[pl_index] +=  fac*ppc[ppc_index]
            pl_index += 1
            ppc_index +=1
    return


def psf_makeresi(psf, set, centflag, prof_accuracy):


    accuflag = (prof_accuracy > 1.0/BIG)
    vigstep = 1/psf.pixstep
    npix = set.vigsize[0]*set.vigsize[1]
    ndim = psf.poly.ndim

    if (centflag):
        cw=ch=int(2*set.fwhm+1.0)
        if (cw>set.vigsize[0]):
            cw=set.vigsize[0]
        
        if (ch>set.vigsize[1]) :
            ch=set.vigsize[1]
        
        ncpix = cw*ch
        hcw = (double)(cw/2)
        hch = (double)(ch/2)
        cvigxt = cvigx
        cvigyt = cvigy
        cvigt_index = 0
        for iy in range(ch):
            yb = iy-hch;
            for ix in range(ch):
                cvigxt[cvigt_index]= ix-hcw
                cvigyt[cvigt_index] = yb
                cvigt_index+=1
    else:
        cvigx = cvigy = None    
        cbasis = cdata = cvigw = None   
        cw = ch = ncpix = 0;            
   
    radmin2 = PSF_MINSHIFT*PSF_MINSHIFT
    radmax2 = PSF_MAXSHIFT*PSF_MAXSHIFT
    okflag = nchi2 = 0
    mse = 0.0
    ngood = 0
    sample = set.sample
    sample_index = 0
    for n in range(set.nsample, -1, -1):
        sample_index+=1
        if (sample.badflag) :
            continue
        
        ngood+=1
        for i in range(ndim):
            pos[i] = (sample.context[i]-set.contextoffset[i])/set.contextscale[i]
        
        psf_build(psf, pos)

        dx = sample.dx
        dy = sample.dy

        if (centflag) :
            vignet_copy(sample.vig, set.vigsize[0], set.vigsize[1], cdata, cw,ch, 0,0, VIGNET_CPY);
            vignet_copy(sample.vigweight, set.vigsize[0], set.vigsize[1], cvigw, cw,ch, 0,0, VIGNET_CPY);
            cdatat=cdata
            cvigwt=cvigw
            cdatat_index = 0
            cvigwt_index = 0
            for i in range(ncpix, -1,-1):
                cdatat[cdatat_index] *= cvigwt[cvigwt_index]
                cdatat_index +=1
                cvigwt_index +=1
            
            for j in range(PSF_NITER):
                vignet_resample(psf.loc, psf.size[0], psf.size[1],
                        cbasis, cw,ch, -dx*vigstep, -dy*vigstep, vigstep, 1.0,
                        sample.vigdgeox, sample.vigdgeoy)

                #memset(amat, 0, 9*sizeof(double));
                bmat[0] = bmat[1] = bmat[2] = mx2=my2=mxy = 0.0
                cvigxt=cvigx
                cvigxt_index = 0
                cvigyt=cvigy
                cvigyt_index = 0
                cvigwt=cvigw
                cvigwt_index = 0
                cbasist=cbasis
                cbasist_index = 0 
                cdatat=cdata
                cdatat_index = 0
                for i in range(ncpix,-1,-1):
                    dval = cbasist[cbasist_index]
                    cbasist_index+=1
                    dwval = dval*cdatat[cdatat_index]
                    cdatat_index +=1
                    bmat[0] += dwval
                    dvalx = cvigxt[cvigxt_index] - dx
                    cvigxt_index+=1
                    bmat[1] += dwval*dvalx
                    dvaly = cvigyt[cvigyt_index] - dy
                    cvigyt_index+=1
                    bmat[2] += dwval*dvaly
                    mx2 += dval*dvalx*dvalx
                    my2 += dval*dvaly*dvaly
                    mxy += dval*dvalx*dvaly
                    amatt=amat
                    amatt_index = 0
                    dval *= dval*cvigwt[cvigwt_index]
                    cvigwt_index+=1
                    amatt[amatt_index] += dval;
                    amatt_index+=1
                    amatt[amatt_index] += dval*dvalx;
                    amatt_index+=1
                    amatt[amatt_index] += dval*dvaly;
                    amatt_index+=2
                    amatt[amatt_index] += dval*dvalx*dvalx;
                    amatt_index+=1
                    amatt[amatt_index] += dval*dvalx*dvaly;
                    amatt[amatt_index] += dval*dvaly*dvaly;
                    amatt_index+=3

                ddx = (bmat[1]*mx2 + bmat[2]*mxy) / bmat[0]
                dx += 0.5*ddx
                ddy = (bmat[2]*my2 + bmat[1]*mxy) / bmat[0]
                dy += 0.5*ddy
                if (ddx*ddx+ddy*ddy < radmin2) :
                    okflag = 1
                    break
                elif (dx*dx+dy*dy > radmax2):
                    break
                
            
            if (okflag) :
                sample.dx = dx
                sample.dy = dy
            
        vignet_resample(psf.loc, psf.size[0], psf.size[1],
                sample.vigresi, set.vigsize[0], set.vigsize[1],
                -dx*vigstep, -dy*vigstep, vigstep, 1.0,
                sample.vigdgeox, sample.vigdgeoy)
        xi2 = xyi = 0.0
        cvigwt=sample.vigweight
        cvigwt_index=0
        cbasist=sample.vigresi
        cbasist_index = 0
        cdatat=sample.vig
        cdatat_index = 0
        for i in range(npix, -1, -1):
            dwval = cvigwt[cvigwt_index]
            cvigwt_index+=1
            dval = cbasist[cbasist_index]
            cbasist_index+=1
            xi2 += dwval*dval*dval;
            xyi += dwval*dval*cdatat[cdatat_index]
            cdatat_index+=1
        if xi2>0.0:
            norm = xyi/xi2
        else:
            norm = sample.norm  

        chi2 = mse = resival = resinorm = 0.0
        dresit = dresi
        psf_extraccu2 = prof_accuracy*prof_accuracy*norm*norm
        xc = (set.vigsize[0]/2)+sample.dx
        yc = (set.vigsize[1]/2)+sample.dy
        y = -yc
        if psf.size[0]<psf.size[1]:
            
            rmax2 = psf.pixstep*(psf.size[0]/2) 
        else:
            rmax2 = psf.pixstep*(psf.size[1]/2)
        rmax2 *= rmax2
        nchi2 = 0
        vig = sample.vig
        vig_index = 0
        vigw = sample.vigweight
        vigw_index = 0
        vigresi=sample.vigresi
        vigchi = sample.vigchi
        vig_index = 0
        vigchi_index = 0
        vigresi_index = 0
        dresit_index = 0
        for iy in range(set.vigsize[1], -1, -1):
            y+=1.0
            x = -xc
#pragma ivdep
            for ix in range(set.vigsize[0], -1, -1):
                x+=1.0
                vig_index+=1
                vigresi_index += 1
                vigchi_index += 1
                dresit_index+=1
                vigchi[vigchi_index] = 0
                wval=vigw[vigw_index]
                vigw_index+=1
                if (wval>0.0):
                    if (accuflag):
                        wval = 1.0/(1.0 / wval + psf_extraccu2* vigresi[vigresi_index]**2);
                    
                    vigresi[vigresi_index] = fval = (vig[vig_index]-vigresi[vigresi_index]*norm)
                    if (x*x+y*y<rmax2):
                        mse += fval*fval
                        nchi2+=1
                        igchi[vigchi_index]=wval*fval*fval
                        chi2 += vigchi[vigchi_index]
                        dresit[dresit_index] += fval
                        sval = vig[vig_index]+vigresi[vigresi_index]*norm
                        resival += sval*fabsf(fval)
                        resinorm += sval*sval
                    
        if nchi2> 1:
            sample.chi2 =  chi2/(nchi2-1)
        else:
            sample.chi2 = chi2
        if resinorm > 0.0:
            sample.modresi = 2.0*resival/resinorm
        else:
            sample.modresi = resival
    
    mse = sqrt(mse/ngood/nchi2)
    if (ngood > 1):
        nm1 = (ngood - 1)
    else:
        nm1 = 1.0
    dresit=dresi
    dresit_index = 0
    fresit=fresi
    fresit_index = 0
    for i in range (npix, -1, -1):
        fresit[fresit_index] = sqrt(dresit[dresit_index]/nm1)
        dresit_index+=1
        fresit_index+=1
    

    vignet_resample(fresi, set.vigsize[0], set.vigsize[1],
            psf.resi, psf.size[0], psf.size[1], 0.0,0.0, psf.pixstep, 1.0,
            None, None)

    



def psf_refine(psf, set):

    if (not set.ngood or not psf.basis) :
        return RETURN_ERROR

    npix = psf.size[0]*psf.size[1]
    nvpix = set.vigsize[0]*set.vigsize[1]
    vigstep = 1/psf.pixstep

    npsf = psf.nbasis
    if psf.ndata :
        ndata = psf.ndata
    else:
        ndata = set.vigsize[0]*set.vigsize[1]+1
    poly = psf.poly
    ncontext = set.ncontext
    ncoeff = poly.ncoeff
    nsample = set.nsample
    nunknown = ncoeff*npsf

    matoffset =nunknown-ncoeff
    sample=set.sample
    for n in range(nsample):
        if (sample[n].badflag):
            continue
        str = "Processing sample #%d" % (n+1)
        dx = -sample[n].dx*vigstep
        dy = -sample[n].dy*vigstep
        norm = sample[n].norm

        for i in range(ncontext):
            pos[i] = (sample[n].context[i]-set.contextoffset[i])/set.contextscale[i];
        
        psf_build(psf, pos)

        basis = poly_ortho(poly, poly.basis, poly.orthobasis)
        basist=basis
        basist_index = 0
        coeffmatt=coeffmat
        coeffmatt_index = 0
        for l in range(ncoeff):
            dval=basist[basist_index]
            basist_index+=1
            basist2=basis
            for i in range(ncoeff):
                coeffmatt[coeffmatt_index] = dval*basist2[i]
                coeffmatt_index += 1
         
        sigvigt=sigvig
        wvig=sample.vigweight
        for i in range(nvpix):
            sigvigt[i] = sqrt(wvig[i])

        desmatt = desmat
        desindext = desindex
        index = 0
        if (psf.pixmask):
            vignet_resample(psf.loc, psf.size[0], psf.size[1], vig, set.vigsize[0], set.vigsize[1], dx, dy, vigstep, 1.0, sample[n].vigdgeox, sample[n].vigdgeoy)
            vigt=vig
            vigt2=sample[n].vig
            for i in range(nvpix):
                vigt[i] = float(vigt2[i] - vigt[i]*norm);
            
        else:
            vigt=vig
            vigt2=sample[n].vig
            for i in range(nvpix):
                vigt[i] = float(vigt2[i])
         
        for i in range(npsf):
            vignet_resample(psf.basis[i*npix], psf.size[0], psf.size[1],
                    vecvig, set.vigsize[0],set.vigsize[1], dx,dy, vigstep, 1.0,
                    sample[n].vigdgeox, sample[n].vigdgeoy)
            vecvigt=vecvig
            sigvigt=sigvig
            desmatt2=desmatt
            desindext2=desindext
            index2 = index
            jo=0
            for j in range(nvpix):
                dval = vecvigt[j] * sigvigt[j]
                if (abs(dval) > (1/BIG)):
                    desmatt2[index2] = norm*dval
                    desindext2[index2] = (j-jo)
                    index2 += 1
                    jo = j

            desindext2[index2] = 0

            index += ndata
        
        vigt=vig
        sigvigt=sigvig
        bmatt=bmat
        for j in range(nvpix):
            bmatt[j] = vigt[j] * sigvigt[j]

        betamatt = betamat
        betamatt_index = 0
        desmat0=desmat
        desindex0=desindex
        des0_index = 0
        for k in range(npsf):
            desmat02=desmat0
            desindex02=desindex0
            des02_index = des0_index
            for j in range(k, npsf):
                dval = 0.0
                desmatt=desmat0
                desmatt2=desmat02
                desindext=desindex0
                dest_index = des0_index
                desindext2=desindex02
                dest2_index = des02_index
                dindex=desindext[dest_index]-desindext2[dest2_index]
                while (desindext[dest_index] and desindext2[dest2_index]):
                    while (desindext[dest_index] and dindex<0):
                        dest_index+=1
                        dindex+=(desindext[dest_index])
                    
                    while (desindext2[dest2_index] and dindex>0):
                        dest2_index+=1
                        dindex-=desindext2[dest2_index]
                    
                    while (desindext[dest_index] and not dindex):
                        dval += desmatt[dest_index]*desmatt2[dest2_index]
                        dest_index+=1
                        dest2_index+=1
                        dindex = desindext[dest_index]-desindext2[dest2_index]
                    
                if (abs(dval) > (1/BIG)):
                    alphamatt = alphamat
                    alphamat_index = (j+k*npsf*ncoeff)*ncoeff
                    coeffmatt=coeffmat
                    coeffmatt_index = 0
                    for l in range(ncoeff):
                        for i in range(ncoeff):
                            alphamatt[alphamat_index] += dval*coeffmatt[coeffmatt_index]
                            coeffmatt_index+=1
                        
                        alphamat_index += matoffset
                 
                des02_index +=ndata
            
            dval = 0.0
            desmatt=desmat0
            desindext=desindex0
            dest_index = des0_index
            bmatt=bmat-1
            bmatt_index = -1
            while (desindext[dest_index]):
                bmatt_index+=desindext[dest_index]
                
                dval += desmatt[dest_index]*bmatt[bmatt_index]
                dest_index+=1
            basist=basis
            basist_index = 0
            for i in range(ncoeff):
                betamatt[betamatt_index] += dval*(basist[basist_index])
                betamatt_index+=1
                basist_index+=1
            
            des0_index += data

    if (psf.pixmask):
        tikfac= 0.01
        tikfac = 1.0/(tikfac*tikfac)
        for i in range(nunknown):
            alphamat[i+nunknown*i] += tikfac

    bcoeff = None
    if (psf.basiscoeff):
        psf.basiscoeff = None
    
    if (not psf.pixmask):
        bcoeff = psf.basiscoeff
    bcoeff_index = 0
    for j in range(npsf):
        ppix = psf.comp
        ppix_index= 0
        betamatt = poly_deortho(poly, betamat + j*ncoeff, betamat2)
        for c in range(ncoeff):
            vec = psf.basis[j*npix];
            dval = betamatt[c]
#pragma ivdep
            for i in range(npix):
                ppix[ppix_index] += dval*vec[i]
                ppix_index+=1
            if (psf.basiscoeff):
                bcoeff[bcoeff_index] = float(dval)
                bcoeff_index+=1

    return RETURN_OK


def  psf_orthopoly(psf, set):
    
    poly = psf.poly
    ncoeff = poly.ncoeff
    ndim = poly.ndim

    ndata = set.ngood
    norm = -1.0/sqrt(ndata)

    sample=set.sample
    for n in range(ndata):
        if (sample[n].badflag):
            continue
        for i in range(ndim):
            pos[i] = (sample.context[i]-set.contextoffset[i])/set.contextscale[i]
        
        poly_func(poly, pos)
        basis = poly.basis
        datat = data + n
        n+=1
#pragma ivdep
        for c in range (ncoeff):
            datat[c] = basis[c]*norm
        
    poly_initortho(poly, data, ndata)

    


def  psf_makebasis(psf, set, basis_type,  nvec):
    npix = psf.size[0]*psf.size[1]

    if basis_type== BASIS_NONE:
        break
    elif basis_type== BASIS_PIXEL:
        npsf = nvec*nvec
        if (npsf>npix):
            npsf=npix
        
 
        psforder = None
        psf.comp =  psforder[:npix] 
        psfordert=psforder
        for i in range(npix):
            psfordert[i] = abs(psfordert[i]);
        
        fqmedian(psforder, npix)

        psfthresh = psforder[npix-npsf]

        psfmask = psf.pixmask
        npsf = 0
        irad = int((set.vigsize[1]-1)/(2*psf.pixstep))
        iymin = psf.size[1]/2 - irad
        iymax = psf.size[1]/2 + irad
        irad = int((set.vigsize[0]-1)/(2*psf.pixstep))
        ixmin = psf.size[0]/2 - irad
        ixmax = psf.size[0]/2 + irad
        ppix=psf.comp
        xc = (psf.size[0]/2)
        yc = (psf.size[1]/2)
        y = -yc
        if psf.size[0]<psf.size[1]:
            rmax2 = (psf.size[0]/2)**2
        else:
            rmax2 =  ((psf.size[1]/2)+0.5)**2
                    
        for iy in range(psf.size[1]):
            y+=1.0
            i = iy*psf.size[0]
            x = -xc
            if (iy>=iymin and iy<=iymax):
                for ix in range(psf.size[0]):
                    x+=1.0
                    i+=1                    
                    if (abs(ppix[i])>=psfthresh and ix>=ixmin and ix<=ixmax and x*x+y*y<rmax2):
                        npsf+=1
                        psfmask[i] = 1

        psf.nbasis = npsf
        basis = psf.basis
        basis_index = 0
        psfmask = psf.pixmask
        for i in range(npix):
            if psfmask[i]:
                basis[basis_index] = 1.0
                basis += 1
 
        psf.ndata = (1+int(INTERPW*psf.pixstep))*(1+int(INTERPW*psf.pixstep))+1
        break

    elif basis_type == BASIS_GAUSS_LAGUERRE:
        psf.nbasis = psf_pshapelet(psf.basis, psf.size[0],psf.size[1],
                nvec, sqrt(nvec+1.0)*prefs.basis_scale)
        break
    elif basis_type ==BASIS_FILE:
        psf.nbasis = psf_readbasis(psf, prefs.basis_name, 0)
        break
    else:
        error(EXIT_FAILURE, "*Internal Error*: unknown PSF vector basis in ",
                "psf_makebasis()")
    


def  psf_laguerre(x, p, q):
    dq = q - 1.0
    if (p==0) :
        return 1.0
    if (p==1):
        return (2.0 - x + dq)
    
    l = 0.0
    lpm2 = 1.0
    lpm1 = 2.0 - x + dq
    dn = 2.0
    for n in range(p-1):
        dn+=1.0
        l = (2.0+(dq-x)/dn)*lpm1 - (1.0+dq/dn)*lpm2
        lpm2 = lpm1
        lpm1 = l
    
    return l



def psf_pshapelet(basis, w, h, nmax, beta):

    kmax = (nmax+1)*(nmax+2)/2

    invbeta2 = 1.0/(beta*beta)
    ostep = 1.0/(GAUSS_LAG_OSAMP)
    ostep2 = ostep*ostep
    odx = 0.5*(ostep - 1.0)
    xc =(w/2)
    yc = (h/2)
    if xc < yc:
        rmax2 = xc
    else:
        rmax2 = yc    
    rmax2 *= rmax2*invbeta2

    fr2t = fr2
    fexpr2t = fexpr2
    fthetat = ftheta
    f_index = 0 
    y = odx - yc
    for iy in range(h) :
        y+=1.0
        x = odx - xc
        for ix in range(w):
            x+=1.0
            y1 = y
            for idy in range(GAUSS_LAG_OSAMP):
                y1+=ostep
                x1 = x
                for idx in range(GAUSS_LAG_OSAMP):
                    x1+=ostep
                    fr2t[f_index] = r2 = (x1*x1+y1*y1)*invbeta2
                    fexpr2t[f_index] = exp(-r2/2.0)
                    fthetat[f_index] = atan2(y1,x1)
                    f_index +=1 
                
    basist = basis[0]
    basist_index = 0
    k=1
    for n in range (nmax):
        for m in range (n%2, n, 2):
            str = "Generating basis vector #%d/%d" % (k, kmax)
            k+=1
            dm = m;
            hnmm = (n-m)/2
            fac = 1.0
            for p in range((n+m)/2, hnmm-1, -1) :
                if (p):
                    fac *= p
            
            fac = sqrt(1.0/(PI*fac))/beta
            if (hnmm%2):
                fac = -fac
            
            fr2t = fr2
            fr2t_index = 0
            fexpr2t = fexpr2
            fthetat = ftheta
            for i in range (w*h):
                val = 0.0
                for j in range(GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP):
                    fr2t_index+=1
                    val += fac*pow(fr2t[fr2t_index], dm/2.0)*psf_laguerre(fr2t[fr2t_index], hnmm, m)
                    fexpr2t[fr2t_index-1]*cos(dm*(fthetat[fr2t_index-1]))
                
                basist[basist_index] = val*ostep2
                basist_index+=1
            
            if (m!=0):
                fr2t = fr2
                fr2t_index = 0

                fexpr2t = fexpr2
                fthetat = ftheta
                for i in range(w*h):
                    val = 0.0
                    for j in range(GAUSS_LAG_OSAMP*GAUSS_LAG_OSAMP):
                        fr2t_index+=1
                    
                        val += fac*pow(fr2t[fr2t_index], dm/2.0)*psf_laguerre(fr2t[fr2t_index], hnmm, m)
                        fexpr2t[fr2t_index]*sin(dm*(fthetat[fr2t_index]))
                    
                    basist[basist_index] = val*ostep2
                    basist_index += 1
                k+=1

    return kmax


def  psf_readbasis(psf, filename, ext):
    cat = read_cat(filename)
    if (not cat):
        error(EXIT_FAILURE, "*Error*: No such catalog: ", filename)
    
    tab = cat.tab
    ntabp1 = cat.ntab+1
    firstab = None
    extp1 = ext+1
    next=0
    while(ntabp1 and next<extp1):
        ntabp1-=1
        tab = tab.nexttab
        if (tab.naxis>=2):
            if (not next):
                firstab = tab
            
            next+=1
   
    if (not ntabp1) :
        if (not next) :
            error(EXIT_FAILURE, "No image data in ", filename)
        
        if (next>extp1):
            warning("Not enough extensions, using only 1st datacube of ",
                    filename)
        
    tab = tab.prevtab
    npixin = tab.naxisn[0]*tab.naxisn[1]
    npixout = psf.size[0]*psf.size[1]
    ncomp = tab.tabsize/tab.bytepix/npixin
    QFSEEK(tab.cat.file, tab.bodypos, SEEK_SET, tab.cat.filename)
    for n in range(ncomp):
        read_body(tab, pixin, npixin)
        vignet_copy(pixin, tab.naxisn[0], tab.naxisn[1],
                psf.basis[n*npixout], psf.size[0], psf.size[1], 0, 0,
                VIGNET_CPY)

    return ncomp