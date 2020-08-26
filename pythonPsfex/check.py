import numpy as np
from define import *
from enum import Enum
from astropy.io.fits import *
from diagnostic import *
from field import *
from astropy.wcs import *
from prefs import *
from sample import *
from vignet import *

class checkenum(Enum):
    PSF_NONE = 1
    PSF_BASIS = 2
    PSF_CHI = 3
    PSF_PROTO = 4
    PSF_RESIDUALS = 5
    PSF_RESIDUALS_GRID = 6
    PSF_SAMPLES = 7
    PSF_SAMPLES_GRID = 8
    PSF_SNAPSHOTS = 9
    PSF_SNAPSHOTS_IMRES = 10
    PSF_WEIGHTS = 11
    PSF_MOFFAT = 12
    PSF_SUBMOFFAT = 13
    PSF_SUBSYM = 14
    PSF_GREAT = 15
    

def check_write(field, set, checkname, checktype, ext, next, cubeflag):
    dpos = np.zeros(POLY_MAXDIM, dtype=np.float64)
    if (not ext):
        cat = new_cat(1)
        init_cat(cat)
        cat.filename= checkname
        pstr = find(cat.filename, '.')
        if (pstr>-1):
            pstr = cat.filename+len(cat.filename)
        
        filename= field.rcatname
        pstr2 = find(filename, '.')
        if (pstf2>-1) :
            pstr2 = filename+len(filename)
        
        pstr2[-1] = '\0'
        pstr = "_%s.fits" % (filename)
        if (open_cat(cat, WRITE_ONLY) != RETURN_OK):
            error(EXIT_FAILURE, "*Error*: cannot open for writing ", cat.filename)
        
        if (next>1):
            addkeywordto_head(cat.tab, "NEXTEND ", "Number of extensions")
            fitswrite(cat.tab.headbuf, "NEXTEND", next, H_INT, T_LONG)
            save_tab(cat, cat.tab)
        
        field.ccat[checktype] = cat
    else:
        cat = field.ccat[checktype]

    str = "chip%02d" % (ext+1)

    psf = field.psf[ext]
    tab = new_tab(str)
    head = tab.headbuf
    tab.bitpix =  BP_FLOAT
    tab.bytepix = t_size[T_FLOAT]

    if checktype == PSF_BASIS:
        if (cubeflag) :
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = psf.size[0]
            tab.naxisn[1] = psf.size[1]
            tab.naxisn[2] = psf.nbasis
            npix = tab.naxisn[0]*tab.naxisn[1]*tab.naxisn[2]
            tab.tabsize = tab.bytepix*npix
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            fpix = psf.basis
            for i in range(npix):
                pix[i] = fpix[i]
            
        else :
            nw = int(sqrt(psf.nbasis))
            nw = ((nw-1)/10+1)*10
            nh = (psf.nbasis-1)/nw + 1
            w = psf.size[0]
            if psf.dim>1:
                h = psf.size[1] 
            else:
                h=1
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            fpix = psf.basis
            fpix_index = 0
            for n in range(psf.nbasis):
                pix = pix0 + ((n%nw) + (n/nw)*nw*h)*w
                pix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        fpix_index += 1
                        pix_index += 1
                    pix += step

    elif checktype ==  PSF_CHI:
         
        if (cubeflag):
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = set.vigsize[0]
            tab.naxisn[1] = set.vigsize[1]
            if set.ngood:
                tab.naxisn[2] =  set.ngood
            else:
                tab.naxisn[2] =  1
            
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            sample = set.sample
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                fpix = sample[n].vigchi
                for i in range(npix):
                    pix[pix_index] = fpix[i]
                

        else:
            nw = int(sqrt(set.ngood))
            if nw:
                nw = ((nw-1)/10+1)*10
            else:
                nw =1
            nh = (set.ngood-1)/nw + 1
            w = set.vigsize[0]
            if set.vigdim>1:
                h =  set.vigsize[1]
            else:
                h=1
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            sample = set.sample
            g = 0
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
        
                pix = pix0 
                pix_index = ((g%nw) + (g/nw)*nw*h)*w;
                fpix = sample.vigchi;
                fpix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        pix_index += 1
                        fpix_index += 1
                    
                    pix_index += step
                g+=1
            
        break

    elif checktype ==  PSF_PROTO:
        if (cubeflag):
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = psf.size[0]
            tab.naxisn[1] = psf.size[1]
            tab.naxisn[2] = psf.size[2]
            npix = tab.naxisn[0]*tab.naxisn[1]*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.tabsize = tab.bytepix*npix
            tab.bodybuf = pix0
            pix = pix0
            fpix = psf.comp
            for i in range(npix):
                pix[i] = fpix[i]
            
        else:
            npc = psf.size[2]
            nw = min(npc, 10)
            nh = (npc-1)/nw + 1
            w = psf.size[0]
            h = psf.size[1]
            step = (nw-1)*w
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            dpost=dpos
            for i in range(psf.poly.ndim):
                dpost[i] = 0.5
            
            poly_func(psf.poly, dpos)
            dpost = psf.poly.basis
            fpix = psf.comp
            fpix_index = 0
            for n in range(pc):
                val = dpost[n]
                pix = pix0 
                pix_index = ((n%nw) + (n/nw)*nw*h)*w
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]*val
                        pix_index += 1
                        fpix_index +=1
                    
                    pix_index += step
        
        break

    elif checktype ==  PSF_RESIDUALS:
        if (cubeflag):
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = set.vigsize[0]
            tab.naxisn[1] = set.vigsize[1]
            if set.ngood:
                tab.naxisn[2] = set.ngood
            else:
                tab.naxisn[2] = 1
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            sample = set.sample
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                fpix = sample[n].vigresi
                for i in range(npix):
                    pix[pix_index] = fpix[i]
                    pix_index += 1

        else:
            nw = int(sqrt(set.ngood))
            if nw :
                nw = ((nw-1)/10+1)*10 
            else:
                nw = 1
            nh = (set.ngood-1)/nw + 1
            w = set.vigsize[0]
            if set.vigdim>1:
                h = set.vigsize[1]
            else:
                h=1
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            sample = set.sample
            g = 0
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                pix = pix0 
                pix_index = ((g%nw) + (g/nw)*nw*h)*w
                fpix = sample.vigresi;
                fpix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        pix_index+=1
                        fpix_index+=1
                    
                    pix_index += step
                
                g+=1
        break

    elif checktype ==  PSF_RESIDUALS_GRID:
        npc = 2
        npos = psf.nsnap
        nt = psf.nsnap**npc
        
        dstepx = prefs.context_nsnap / field.wcs[ext].naxisn[0]
        dstepy = prefs.context_nsnap / field.wcs[ext].naxisn[1]
        sample = set.sample
        for n in range(set.nsample):
            if (sample[n].badflag):
                continue
            
            ipos = int(dstepx * (sample[n].x+0.5)) ;
            if (ipos<0) :
                ipos = 0
            elif (ipos>=npos):
                ipos = npos-1
            
            inpos = ipos
            ipos = int(dstepy * (sample[n].y+0.5))
            if (ipos<0):
                ipos = 0
            elif (ipos>=npos):
                ipos = npos-1
            
            inpos += ipos*npos
            osample = gridsample[inpos]
            if (not osample or osample.norm < sample.norm):
                gridsample[inpos] = sample[n]
            
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
        w = set.vigsize[0]
        h = set.vigsize[1]
        if (cubeflag):
            if npc>2:
                nh = psf.nsnap 
                np = nt/(nw*nh)
            else:
                nh = nt/nw
                np = 1
            tab.naxis = 4
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = w
            tab.naxisn[1] = h
            tab.naxisn[2] = nw
            tab.naxisn[3] = nh
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]*tab.naxisn[3]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            for n in range(nt):
                sample = gridsample[n]
                if sample:
                    pix[pix_index: pix_index+npix]= sample.vigresi
                    
                pix_index += npix
            
        else :
            nh = nt/nw
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            step = (nw-1)*w
            for n in range(nt):
                sample = gridsample[n]
                if sample:
                    pix = pix0 
                    pix_index = ((n%nw) + (n/nw)*nw*h)*w
                    fpix = sample.vigresi
                    fpix_index = 0
                    for y in range(h):
                        for x in range(w):
                            pix[pix_index] = fpix[fpix_index]
                            pix_index+=1
                            fpix_index+=1
                        
                        pix_index += step
                   
        break

    elif checktype == PSF_SAMPLES:
        if (cubeflag):
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = set.vigsize[0]
            tab.naxisn[1] = set.vigsize[1]
            if set.ngood:
                tab.naxisn[2] = set.ngood
            else:
                tab.naxisn[2] = 1
                
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            sample = set.sample
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                fpix = sample[n].vig;
                for i in range(npix):
                    pix[pix_index] = fpix[i]
                    pix_index += 1
                
        else:
            nw = int(sqrt(set.ngood))
            if nw:
                nw = ((nw-1)/10+1)*10
            else:
                nw = 1
    
            nh = (set.ngood-1)/nw + 1
            w = set.vigsize[0]
            if set.vigdim>1:
                h = set.vigsize[1]
            else:
                h=1
                
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            sample = set.sample
            g = 0
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                pix = pix0 
                pix_index = ((g%nw) + (g/nw)*nw*h)*w
                fpix = sample.vig
                fpix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        pix_index+=1
                        fpix_index+=1
                    
                    pix_index += step
                
                g+=1
            
        break

    elif checktype == PSF_SAMPLES_GRID:
        npc = 2
        npos = psf.nsnap
        nt = psf.nsnap**npc;
        
        dstepx = prefs.context_nsnap / field.wcs[ext].naxisn[0]
        dstepy = prefs.context_nsnap / field.wcs[ext].naxisn[1]
        sample = set.sample
        for n in range(set.nsample):
            if (sample[n].badflag):
                continue
            ipos = int(dstepx * (sample.x+0.5))
            if (ipos<0):
                ipos = 0
            elif (ipos>=npos):
                ipos = npos-1
            
            inpos = ipos
            ipos = int(dstepy * (sample.y+0.5))
            if (ipos<0):
                ipos = 0
            elif (ipos>=npos):
                ipos = npos-1
            
            inpos += ipos*npos
            osample = gridsample[inpos]
            if (not osample or osample.norm < sample.norm):
                gridsample[inpos] = sample[n]

        if npc:
            nw = psf.nsnap
        else:
            nw = 1
        w = set.vigsize[0]
        h = set.vigsize[1]
        if (cubeflag):
            if npc > 2:
                nh = psf.nsnap 
                np = nt/(nw*nh)
            else:
                nh = nt/nw
                np = 1
            tab.naxis = 4
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = w
            tab.naxisn[1] = h
            tab.naxisn[2] = nw
            tab.naxisn[3] = nh
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]*tab.naxisn[3]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            for n in range(nt):
                sample = gridsample[n]
                if sample :
                    pix[pix_index:pix_index+npix] = sample.vig
                
                pix_index += npix
        else:
            nh = nt/nw
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            step = (nw-1)*w
            for n in range(nt):
                sample = gridsample[n]
                if sample:
                    pix = pix0 
                    pix_index = ((n%nw) + (n/nw)*nw*h)*w
                    fpix = sample.vig
                    fpix_index= 0
                    for y in range(h):
                        for x in range(w):
                            pix[pix_index] = fpix[fpix_index]
                            pix_index+=1
                            fpix_index+=1
                        
                        pix_index += step
                    
        break

    elif checktype ==  PSF_WEIGHTS:
        if (cubeflag) :
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = set.vigsize[0]
            tab.naxisn[1] = set.vigsize[1]
            if set.ngood:
                tab.naxisn[2] = set.ngood
            else:
                tab.naxisn[2] = 1
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            sample = set.sample
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                fpix = sample[n].vigweight
                for i in range(npix):
                    pix[pix_index] = fpix[i]
                    pix_index+=1
                
        else:
            nw = int(sqrt(set.ngood))
            if nw:
                nw = ((nw-1)/10+1)*10 
            else:
                nw=1
            nh = (set.ngood-1)/nw + 1
            w = set.vigsize[0]
            if set.vigdim>1:
                h = set.vigsize[1]
            else:
                h=1
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            sample = set.sample
            g = 0
            for n in range(set.nsample):
                if (sample[n].badflag):
                    continue
                
                pix = pix0 
                pix_index = ((g%nw) + (g/nw)*nw*h)*w
                fpix = sample[n].vigweight
                fpix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        pix_index+=1
                        fpix_index+=1
                    pix_index += step
                g+=1
        break

    elif checktype ==  PSF_SNAPSHOTS:
        npc = psf.poly.ndim
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
        nt = psf.nsnap**npc;
    
    
        w = psf.size[0]
        if  psf.dim>1:
            h = psf.size[1]
        else:
            h=1
        dstep = 1.0/prefs.context_nsnap
        dstart = (1.0-dstep)/2.0
        for i in range(npc):
            dpos[i] = -dstart
        
        if (cubeflag):
            if npc>2:
                nh = psf.nsnap 
                np = nt/(nw*nh)
            else:
                nh = nt/nw
                np = 1
            tab.naxis = 4
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = w
            tab.naxisn[1] = h
            tab.naxisn[2] = nw
            tab.naxisn[3] = nh
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]*tab.naxisn[3]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            for n in range(nt):
                psf_build(psf, dpos)
                fpix = psf.loc
                for i in range(npix):
                    pix[pix_index] = fpix[i]
                    pix_index+=1
                
                for i in range(npc):
                    if (dpos[i]<dstart-0.01):
                        dpos[i] += dstep
                        break
                    else:
                        dpos[i] = -dstart
        else :
            nh = nt/nw
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            step = (nw-1)*w
            for n in range(nt):
                psf_build(psf, dpos)
                pix = pix0 
                pix_index = ((n%nw) + (n/nw)*nw*h)*w
                fpix = psf.loc
                fpix_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = fpix[fpix_index]
                        pix_index +=1
                        fpix_index += 1
                    
                    pix_index += step
                
                for i in range(npc):
                    if (dpos[i]<dstart-0.01):
                        dpos[i] += dstep
                        break
                    else:
                        dpos[i] = -dstart
                    
        break

    elif checktype ==  PSF_SNAPSHOTS_IMRES:
        npc = psf.poly.ndim
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
            nt = psf.nsnap**npc
        w = set.vigsize[0]
        if set.vigdim>1:
            h = set.vigsize[1]
        else:
            h=1
        vig0 = np.zeros(w*h, dtype=np.float32)
        dstep = 1.0/prefs.context_nsnap
        dstart = (1.0-dstep)/2.0
        for i in range(npc):
            dpos[i] = -dstart
        
        if (cubeflag):
            if npc>2:
                nh = psf.nsnap 
            else:
                nh = nt/nw;
            tab.naxis = 4
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = w
            tab.naxisn[1] = h
            tab.naxisn[2] = nw
            tab.naxisn[3] = nh
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]*tab.naxisn[3]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            pix_index = 0
            for n in range (nt):
                psf_build(psf, dpos)
                vignet_resample(psf.loc, psf.size[0], psf.size[1],
                        vig0, set.vigsize[0], set.vigsize[1], 0.0, 0.0,
                        1.0/psf.pixstep, 1.0, None, None)
                vig = vig0
                for i in range(npix):
                    pix[pix_index] = vig[i]
                    pix_index+=1
                
                for i in range(npc):
                    if (dpos[i]<dstart-0.01):
                        dpos[i] += dstep
                        break
                    else: 
                        dpos[i] = -dstart
                    
        else:
            nh = nt/nw
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            step = (nw-1)*w
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)

            tab.bodybuf = pix0
            for n in range(nt):
                psf_build(psf, dpos)
                vignet_resample(psf.loc, psf.size[0], psf.size[1],
                        vig0, set.vigsize[0], set.vigsize[1], 0.0, 0.0,
                        1.0/psf.pixstep, 1.0, None, None)
                vig = vig0
                pix = pix0 
                pix_index = ((n%nw) + (n/nw)*nw*h)*w
                vig_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = vig[vig_index]
                        pix_index += 1
                        vig_index += 1
                    pix_index += step
                
                for i in range(npc):
                    if (dpos[i]<dstart-0.01):
                        dpos[i] += dstep
                        break
                    else :
                        dpos[i] = -dstart
                    
        break

    elif checktype ==  PSF_MOFFAT:
        npc = psf.poly.ndim
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
        nt = psf.nsnap ** npc
        
        nh = nt/nw
        w = psf.size[0]
        if psf.dim>1:
            h =psf.size[1]
        else:
            h=1
        tab.naxisn[0] = nw*w
        tab.naxisn[1] = nh*h
        step = (nw-1)*w
        tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
        pix0 = np.zeros(tab.tabsize, dtype=np.float32)

        tab.bodybuf = pix0
        dstep = 1.0/psf.nsnap
        dstart = (1.0-dstep)/2.0
        for i in range(npc):
            dpos[i] = -dstart
        
        for n in range(nt):
            psf_moffat(psf, psf.pfmoffat[n])
            pix = pix0 
            pix_index = ((n%nw) + (n/nw)*nw*h)*w
            fpix = psf.loc
            fpix_index = 0
            for y in range(h):
                for x in range(w):
                    pix[pix_index] = fpix[fpix_index]
                    pix_index+=1
                    fpix_index+=1
                pix_index+=step
            
            for i in range(npc):
                if (dpos[i]<dstart-0.01):
                    dpos[i] += dstep
                    break
                else:
                    dpos[i] = -dstart
                
        break

    elif checktype ==  PSF_SUBMOFFAT:
        npc = psf.poly.ndim
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
       
        nt =psf.nsnap**npc
        nh = nt/nw
        w = psf.size[0]
        if psf.dom > 1:
            h = psf.size[1]
        else:
            h=1
        tab.naxisn[0] = nw*w
        tab.naxisn[1] = nh*h
        step = (nw-1)*w
        tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
        pix0 = np.zeros(tab.tabsize, dtype=np.float32)
        tab.bodybuf = pix0
        dstep = 1.0/psf.nsnap
        dstart = (1.0-dstep)/2.0
        for i in range(npc):
            dpos[i] = -dstart
        
        for n in range(nt):
            psf_build(psf, dpos)
            pix = pix0 
            pix_index = ((n%nw) + (n/nw)*nw*h)*w
            fpix = psf.loc
            fpix_index = 0
            for y in range(h):
                for x in range(w):
                    pix[pix_index] = fpix[fpix_index]
                    pix_index+=1
                    fpix_index+=1
                pix_index+=step
        
            psf_moffat(psf, psf.pfmoffat[n])
            pix = pix0 
            pix_index = ((n%nw) + (n/nw)*nw*h)*w
            fpix = psf.loc
            fpix_index = 0
            for y in range(h):
                for x in range(w):
                    pix[pix_index] -= fpix[fpix_index]
                    pix_index+=1
                    fpix_index+=1
                
                pix_index+=step
            
            for i in range(npc):
                if (dpos[i]<dstart-0.01):
                    dpos[i] += dstep
                    break
                else:
                    dpos[i] = -dstart
                
        break
    elif checktype ==  PSF_SUBSYM:
        npc = psf.poly.ndim;
        if npc:
            nw = psf.nsnap
        else:
            nw = 1
        nt = psf.nsnap**npc
        
        nh = nt/nw
        w = psf.size[0]
        if psf.dim>1:
            h=psf.size[1]
        else:
            h=1
        tab.naxisn[0] = nw*w
        tab.naxisn[1] = nh*h
        step = (nw-1)*w
        tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
        pix0 = np.zeros(tab.tabsize, dtype=np.float32)

        tab.bodybuf = pix0
        dstep = 1.0/psf.nsnap
        dstart = (1.0-dstep)/2.0
        for i in range(npc):
            dpos[i] = -dstart
        
        for n in range(nt):
            psf_build(psf, dpos)
            pix = pix0 
            pix_index = ((n%nw) + (n/nw)*nw*h)*w
            fpix = psf.loc
            fpix_index = 0
            fpixsym_index = w*h
            for y in range(h):
                for x in range(w):
                    fpixsym_index -= 1
                    pix[pix_index] = fpix[fpix_index]-fpixsym[fpixsym_index]
                    fpix_index += 1
                    pix_index+=1
                pix_index+=step
            
            for i in range(npc):
                if (dpos[i]<dstart-0.01):
                    dpos[i] += dstep
                    break
                else:
                    dpos[i] = -dstart
                
        break
    elif checktype ==  PSF_GREAT:
        NFPRINTF(OUTPUT, "Reading list of input positions...")
        listfilename = field.catname.copy()
        pstr = rfind(listfilename, '.')
        if (not pstr) :
            pstr = listfilename+len(listfilename)
        
        pstr = "%s" % (".dat")
        listfile = open(listfilename, "r")
        if not listfile:
            error(EXIT_FAILURE, "*Error*: position list not found: ", listfilename)
        
        if (psf.cx>=0):
            scalex = psf.contextscale[psf.cx]
            invcscalex = 1.0/scalex
            coffx = psf.contextoffset[psf.cx]
        else:
            invcscalex = scalex = 1.0
            coffx = 0.0
        
        if (psf.cy>=0):
            scaley = psf.contextscale[psf.cy]
            invcscaley = 1.0/scaley
            coffy = psf.contextoffset[psf.cy]
        else:
            invcscaley = scaley = 1.0
            coffy = 0.0

        list = None
        size = 0
        i=0
        str = filelist.readline()
        while(str):
            if (not str or find("#\t\n", str)):
                continue
            
            if not i:
                ispoon = 1000;
                spoonsize = ispoon*2
                size = spoonsize
                list = np.zeros(size)
            if not(i%1000):
                str2 = "Reading input list... (%d objects)" % (i)
                NFPRINTF(OUTPUT, str2)
                size += spoonsize
                list = np.zeros(size)
            list[2*i] = (float(str)-coffx)*invcscalex
            list[2*i+1] = (float(ptr)-coffy)*invcscaley
            i+=1
            str = filelist.readline()
        
        listfile.close()
        nlist = i
        if (not nlist):
            warning("No valid positions found in ", listfilename)
        
        if (cubeflag):
            tab.naxis = 3
            np.resize(tab.naxisn, naxis)
            tab.naxisn[0] = 48
            tab.naxisn[1] = 48
            if nlist:
                tab.naxisn[2] = nlist
            else:
                tab.naxisn[2] = 1
            npix = tab.naxisn[0]*tab.naxisn[1]
            tab.tabsize = tab.bytepix*npix*tab.naxisn[2]
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            tab.bodybuf = pix0
            pix = pix0
            for n in range(nlist):
                psf_build(psf, list[2*n]);
                list[2*n] = list[2*n]*scalex + coffx
                list[2*n+1] = list[2*n+1]*scaley + coffy
                vignet_resample(psf.loc, psf.size[0], psf.size[1],
                        pix, 48, 48,
                        (floor(list[2*n]+0.49999) - list[2*n] + 0.5) / psf.pixstep,
                        (floor(list[2*n+1]+0.49999) - list[2*n+1] + 0.5) / psf.pixstep,
                        1.0/psf.pixstep, 1.0, None, None)
                pix += 48*48
            
        else:
            nw = int(sqrt(nlist))
            if nw:
                nw = ((nw-1)/10+1)*10
            else:
                nw = 1
            nh = (nlist-1)/nw + 1
            w = 48
            h = 48
            tab.naxisn[0] = nw*w
            tab.naxisn[1] = nh*h
            tab.tabsize = tab.bytepix*tab.naxisn[0]*tab.naxisn[1]
            step = (nw-1)*w
            pix0 = np.zeros(tab.tabsize, dtype=np.float32)
            vig0 = np.zeros(w*h, dtype=np.float32)
            tab.bodybuf = pix0
            for n in range(nlist):
                psf_build(psf, list[2*n])
                list[2*n] = list[2*n]*scalex + coffx
                list[2*n+1] = list[2*n+1]*scaley + coffy
                vignet_resample(psf.loc, psf.size[0], psf.size[1],
                        vig0, 48, 48,
                        (floor(list[2*n]+0.49999) - list[2*n] + 0.5) / psf.pixstep,
                        (floor(list[2*n+1]+0.49999) - list[2*n+1] + 0.5) / psf.pixstep,
                        1.0/psf.pixstep, 1.0, None, None)
                vig = vig0
                pix = pix0 
                pix_index = ((n%nw) + (n/nw)*nw*h)*w
                vig_index = 0
                for y in range(h):
                    for x in range(w):
                        pix[pix_index] = vig[vig_index]
                        pix_index+=1
                        vig_index+=1
                    pix_index+=step
                
        break
    else:
        error(EXIT_FAILURE, "*Internal Error*: Yet unavailable CHECKIMAGE type",
                "")

    fitsread(set.head, "NAXIS1  ", ival1, H_INT, T_LONG)
    fitsread(set.head, "NAXIS2  ", ival2, H_INT, T_LONG)
    if (fitsread(set.head, "CTYPE1  ", str, H_STRING, T_STRING)==RETURN_OK):
        addkeywordto_head(tab,"CTYPE1  ","WCS axis type")
        fitswrite(tab.headbuf, "CTYPE1  ", str,  H_STRING, T_STRING)
        if (fitsread(set.head, "CTYPE2  ", str, H_STRING, T_STRING)==RETURN_OK):
            addkeywordto_head(tab,"CTYPE2  ","WCS axis type")
            fitswrite(tab.headbuf, "CTYPE2  ", str,  H_STRING, T_STRING)
        if (tab.naxis>2):
            addkeywordto_head(tab,"CTYPE3  ","WCS axis type");
            fitswrite(tab.headbuf, "CTYPE3  ", " ",  H_STRING, T_STRING)
        if (tab.naxis>3) :
            addkeywordto_head(tab,"CTYPE4  ","WCS axis type")
            fitswrite(tab.headbuf, "CTYPE4  ", " ",  H_STRING, T_STRING)
    if (fitsread(set.head, "CRVAL1  ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "CRVAL1  ",
                "WCS coordinates of the reference pixel")
        fitswrite(tab.headbuf, "CRVAL1  ", dval1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CRVAL2  ", dval2, H_EXPO, T_DOUBLE)==RETURN_OK):
            addkeywordto_head(tab, "CRVAL2  ",
                    "WCS coordinates of the reference pixel")
            fitswrite(tab.headbuf, "CRVAL2  ", dval2,  H_EXPO, T_DOUBLE)
        if (tab.naxis>2):
            dval1 = 1.0;
            addkeywordto_head(tab,"CRVAL3  ",
                    "WCS coordinates of the reference pixel")
            fitswrite(tab.headbuf, "CRVAL3  ", dval1,  H_EXPO, T_DOUBLE)
        if (tab.naxis>3):
            dval1 = 1.0
            addkeywordto_head(tab,"CRVAL4  ",
                    "WCS coordinates of the reference pixel")
            fitswrite(tab.headbuf, "CRVAL4  ", dval1,  H_EXPO, T_DOUBLE)
    dval1=ival1/tab.naxisn[0]
    dval2=ival2/tab.naxisn[1]
    scalefac = min(dval1, dval2)
    if (fitsread(set.head, "CD1_1   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        dval1 *= scalefac
        addkeywordto_head(tab, "CD1_1   ", "WCS transformation matrix")
        fitswrite(tab.headbuf, "CD1_1   ", dval1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CD1_2   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
            dval1 *= scalefac
            addkeywordto_head(tab, "CD1_2   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD1_2   ", val1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CD2_1   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
            dval1 *= scalefac
            addkeywordto_head(tab, "CD2_1   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD2_1   ", dval1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CD2_2   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
            dval1 *= scalefac
            addkeywordto_head(tab, "CD2_2   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD2_2   ", dval1,  H_EXPO, T_DOUBLE)

        if (tab.naxis>2):
            dval1 = 0.0
            addkeywordto_head(tab, "CD1_3   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD1_3   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD2_3   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD2_3   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD3_1   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD3_1   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD3_2   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD3_2   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 1.0
            addkeywordto_head(tab, "CD3_3   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD3_3   ", dval1,  H_EXPO, T_DOUBLE)
        
        if (tab.naxis>3):
            dval1 = 0.0
            addkeywordto_head(tab, "CD1_4   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD1_4   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD2_4   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD2_4   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD3_4   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD3_4   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD4_1   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD4_1   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD4_2   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD4_2   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 0.0
            addkeywordto_head(tab, "CD4_3   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD4_3   ", dval1,  H_EXPO, T_DOUBLE)
            dval1 = 1.0
            addkeywordto_head(tab, "CD4_4   ", "WCS transformation matrix")
            fitswrite(tab.headbuf, "CD4_4   ", dval1,  H_EXPO, T_DOUBLE)
        
    if (fitsread(set.head, "CDELT1  ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        dval1 *= scalefac
        addkeywordto_head(tab, "CDELT1  ", "WCS pixel scale")
        fitswrite(tab.headbuf, "CDELT1  ", dval1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CDELT2  ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
            dval1 *= scalefac
            addkeywordto_head(tab, "CDELT2  ", "WCS pixel scale")
            fitswrite(tab.headbuf, "CDELT2  ", dval1,  H_EXPO, T_DOUBLE)
        if (tab.naxis>2):
            dval1 = 1.0
            addkeywordto_head(tab,"CDELT3  ","WCS pixel scale")
            fitswrite(tab.headbuf, "CDELT3  ", dval1,  H_EXPO, T_DOUBLE)
        if (tab.naxis>3):
            dval1 = 1.0
            addkeywordto_head(tab,"CDELT4  ","WCS pixel scale")
            fitswrite(tab.headbuf, "CDELT4  ", dval1,  H_EXPO, T_DOUBLE)
        
    if (fitsread(set.head, "CRPIX1  ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "CRPIX1  ",
                "pixel coordinates of the reference pixel")
        dval1 = (dval1 - (ival1+1)/2)/scalefac + (tab.naxisn[0]+1)/2.0
        fitswrite(tab.headbuf, "CRPIX1  ", dval1,  H_EXPO, T_DOUBLE)
        if (fitsread(set.head, "CRPIX2  ", dval2, H_EXPO, T_DOUBLE)==RETURN_OK):
            addkeywordto_head(tab, "CRPIX2  ",
                    "pixel coordinates of the reference pixel")
            dval2 = (dval2 - (ival2+1)/2)/scalefac + (tab.naxisn[1]+1)/2.0
            fitswrite(tab.headbuf, "CRPIX2  ", dval2,  H_EXPO, T_DOUBLE)
        if (tab.naxis>2):
            dval1 = 1.0
            addkeywordto_head(tab,"CRPIX3  ",
                    "pixel coordinates of the reference pixel")
            fitswrite(tab.headbuf, "CRPIX3  ", dval1,  H_EXPO, T_DOUBLE)
        
        if (tab.naxis>3):
            dval1 = 1.0;
            addkeywordto_head(tab,"CRPIX4  ",
                    "pixel coordinates of the reference pixel")
            fitswrite(tab.headbuf, "CRPIX4  ", dval1,  H_EXPO, T_DOUBLE)
    
    if (fitsread(set.head, "EPOCH   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "EPOCH   ", "Epoch of observation")
        fitswrite(tab.headbuf, "EPOCH   ", dval1,  H_EXPO, T_DOUBLE)
    if (fitsread(set.head, "EQUINOX ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "EQUINOX ", "Equinox")
        fitswrite(tab.headbuf, "EQUINOX ", dval1,  H_EXPO, T_DOUBLE)
    if (fitsread(set.head, "RADECSYS", str, H_STRING, T_STRING)==RETURN_OK):
        addkeywordto_head(tab, "RADECSYS", "Coordinate system")
        fitswrite(tab.headbuf, "RADECSYS", str,  H_STRING, T_STRING)
    if (fitsread(set.head, "LONGPOLE", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "LONGPOLE", "Longitude of pole")
        fitswrite(tab.headbuf, "LONGPOLE", dval1,  H_EXPO, T_DOUBLE)
    if (fitsread(set.head, "LATPOLE ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "LATPOLE ", "Latitude of pole")
        fitswrite(tab.headbuf, "LATPOLE ", dval1,  H_EXPO, T_DOUBLE)
    if (fitsread(set.head, "EPOCH   ", dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
        addkeywordto_head(tab, "EPOCH   ", "Epoch of observation")
        fitswrite(tab.headbuf, "EPOCH   ", dval1,  H_EXPO, T_DOUBLE)
    if (fitsfind(set.head, "PV?_????") != RETURN_ERROR):
        for l in range(2):
            for j in range(100):
                str = "PV%d_%d " % (l+1, j)
                if (fitsread(set.head, str, dval1, H_EXPO, T_DOUBLE)==RETURN_OK):
                    addkeywordto_head(tab, str, "Distortion parameter")
                    fitswrite(tab.headbuf, str, dval1,  H_EXPO, T_DOUBLE)

    addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data")
    val = 0.0
    if psf.samples_accepted:
        fitswrite(tab.headbuf, "PSF_SAMP", psf.pixstep, H_FLOAT, T_FLOAT)
    else:
        fitswrite(tab.headbuf, "PSF_SAMP", val, H_FLOAT, T_FLOAT)

    if (next == 1):
        prim_head(tab)
    
    fitswrite(tab.headbuf, "XTENSION", "IMAGE   ", H_STRING, T_STRING)

    save_tab(cat, tab)
    free_tab(tab)
    if (ext==next-1):
        free_cat(cat, 1)

    return
