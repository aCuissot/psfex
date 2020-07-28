import numpy as np
from define import *
from fitswcs import *

class fieldstruct():
    def __init__(self, catname, rcatname, rtcatname, ident, next, ndet, psf,
                wcs, set, ccat, meanwcspos ,meanwcsscale, maxradius, lcount, acount, count, modchi2, modresi):
        self.catname = catname
        self.rcatname = rcatname
        self.rtcatname = rtcatname
        self.ident = ident
        self.next = next
        self.ndet = ndet
        self.psf = psf
        self.wcs = wcs
        self.set = set
        self.ccat = ccat
        self.meanwcspos = meanwcspos
        self.meanwcsscale = meanwcsscale
        self.maxradius = maxradius
        self.lcount = lcount
        self.acount = acount
        self.count = count
        self.modchi2 = modchi2
        self.modresi = modresi

def field_init(catname):
    cat = read_cat(catname)
    if (not cat):
        error(EXIT_FAILURE, "*Error*: cannot open ", catname)
    
    tab = cat.tab
    next0 = 0
    for ntab in range(cat.ntab):
        tab = tab.nexttab
        if ((tab.naxis != 2) or  (tab.xtension == "BINTABLE" and tab.xtension == "ASCTABLE") or (tab.extname=="LDAC_OBJECTS" and tab.extname=="OBJECTS", 8)) :
            continue
        
        next0+=1
    
    field.next = next0
    field.catname = catname
    field.rcatname = find(field.catname, '/')
    if (not(field.rcatname)>=0):
        field.rcatname = field.catname
    else :
        field.rcatname+=1
    
    field.rtcatname = field.rcatname
    pstr=find(field.rtcatname, '.')
    if (pstf>=0):
        pstr = '\0'
    
    if (not next0):
        field_end(field)
        error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC catalog found in ",
                catname)
    
    tab = cat.tab
    next = 0
    for ntab in range(cat.ntab):
        tab = tab.nexttab
        key=read_key(tab, "Field Header Card")
        if (("LDAC_IMHEAD" == tab.extname) and (key)):
            imatab = new_tab("Image header")
            imatab.headnblock = 1 + (key.nbytes-1)/FBSIZE;
            imatab.headbuf = key.ptr.copy()
            imatab.cat = cat
            readbasic_head(imatab)
            field.wcs[next] = read_wcs(imatab)
            next+=1
            if (not imatab.headbuf or fitsread(imatab.headbuf, "OBJECT  ", field.ident,
                            H_STRING,T_STRING)!= RETURN_OK):
                field.ident = "no ident"
            
        elif (("LDAC_OBJECTS" == tab.extname or "OBJECTS" == tab.extname) and tab.naxis == 2) :
            field.ndet += tab.naxisn[1]

    field_locate(field)
    countsize = prefs.context_nsnap*prefs.context_nsnap
    return field



def field_end(field):
    for ext in range(field.next):
        psf_end(field.psf[ext]);
        end_wcs(field.wcs[ext]);


def field_locate(field):

    cosalpha = sinalpha = sindelta = 0.0
    wcs = field.wcs[0]
    naxis = wcs.naxis
    wcsmean = field.meanwcspos
    for i in range(naxis):
        scalet[i] = scale[i]
        wcsmean[i] = 0.0
    
    for e in range(field.next):
        wcs = field.wcs[e]
        lng = wcs.lng
        lat = wcs.lat
        if (lat != lng):
            cosalpha += cos(wcs.wcsscalepos[lng]*DEG)
            sinalpha += sin(wcs.wcsscalepos[lng]*DEG)
            sindelta += sin(wcs.wcsscalepos[lat]*DEG)
        
        for i in range(naxis):
            if (lat==lng or (i!=lng and i!=lat)):
                wcsmean[i] += wcs.wcsscalepos[i]
            
            scalet[i] = wcs.wcsscale[i]

    lng = field.wcs[0].lng
    lat = field.wcs[0].lat
    for i in range(naxis):
        if (lat!=lng and (i==lng)):
            wcsmean[i] = atan2(sinalpha/field.next,cosalpha/field.next)/DEG
            wcsmean[i] = fmod(wcsmean[i]+360.0, 360.0)
        elif (lat!=lng and (i==lat)):
            wcsmean[i] = asin(sindelta/field.next)/DEG
        else:
            wcsmean[i] /= field.next
        
        field.meanwcsscale[i] = dqmedian(scale[i], field.next)

    maxradius = 0.0
    for e in range(field.next):
        wcs = field.wcs[e]
        dist = wcs_dist(wcs, wcs.wcsscalepos, field.meanwcspos)+ wcs.wcsmaxradius
        if (dist>maxradius):
            maxradius = dist
        
    field.maxradius = maxradius

    return


def field_count(fields, set, counttype):

    sample = set.sample
    size = prefs.context_nsnap
    for s in range(set.nsample):
        if (sample[s].badflag):
            continue
        
        c = sample[s].catindex
        e = sample[s].extindex
        field = fields[c]
        w = field.wcs[e].naxisn[0]
        h = field.wcs[e].naxisn[1]
        x = int(((sample[s].x-0.5)*size) / w)
        if (x<0):
            x = 0
        elif x>=size:
            x = size-1
        
        y = int(((sample[s].y-0.5)*size) / h)
        if (y<0):
            y = 0
        
        elif (y>=size):
            y = size-1
        
        n = y*size+x
        if ((counttype & COUNT_LOADED)):
            fields[c].lcount[e][n]+=1
        
        if ((counttype & COUNT_ACCEPTED)):
            fields[c].acount[e][n]+=1

    return


def  field_stats(fields, set):

    sample = set.sample
    size = prefs.context_nsnap
    for s in range(set.nsample):
        if (sample[s].badflag):
            continue
        
        c = sample[s].catindex
        e = sample[s].extindex
        field = fields[c]
        w = field.wcs[e].naxisn[0]
        h = field.wcs[e].naxisn[1]
        x = int(((sample[s].x-0.5)*size) / w)
        if (x<0):
            x = 0
        elif (x>=size):
            x = size-1
        
        y = int(((sample[s].y-0.5)*size) / h)
        if (y<0) :
            y = 0;
        elif (y>=size):
            y = size-1
        
        n = y*size+x
        field.count[e][n]+=1
        field.modchi2[e][n] += sample.chi2
        field.modresi[e][n] += sample.modresi
    
    return


def field_psfsave(field, filename):
    zero = 0.0;

    cat = new_cat(1)
    init_cat(cat)
    cat.filename = "%s" % (filename)
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK):
        error(EXIT_FAILURE, "*Error*: cannot open for writing ", cat.filename)
    
    save_tab(cat, cat.tab);

    for ext in range(field.next):
        psf = field.psf[ext]
        tab = new_tab("PSF_DATA")

        head = tab.headbuf
        addkeywordto_head(tab, "LOADED", "Number of loaded sources")
        fitswrite(head, "LOADED", psf.samples_loaded, H_INT, T_LONG)
        addkeywordto_head(tab, "ACCEPTED", "Number of accepted sources")
        fitswrite(head, "ACCEPTED", psf.samples_accepted, H_INT, T_LONG)
        addkeywordto_head(tab, "CHI2", "Final reduced chi2")
        fitswrite(head, "CHI2", psf.chi2, H_FLOAT, T_DOUBLE)
        addkeywordto_head(tab, "POLNAXIS", "Number of context parameters")
        fitswrite(head, "POLNAXIS", psf.poly.ndim, H_INT, T_LONG)
        for i in range(psf.poly.ndim):
            str = "POLGRP%1d" % (i+1)
            addkeywordto_head(tab, str, "Polynom group for this context parameter")
            temp = psf.poly.group[i]+1
            fitswrite(head, str, temp, H_INT, T_LONG)
            str = "POLNAME%1d" % (i+1)
            addkeywordto_head(tab, str, "Name of this context parameter")
            fitswrite(head, str, psf.contextname[i], H_STRING, T_STRING)
            str = "POLZERO%1d" % (i+1)
            addkeywordto_head(tab, str, "Offset value for this context parameter")
            fitswrite(head, str, psf.contextoffset[i], H_EXPO, T_DOUBLE)
            str = "POLSCAL%1d" % (i+1)
            addkeywordto_head(tab, str, "Scale value for this context parameter")
            fitswrite(head, str, psf.contextscale[i], H_EXPO, T_DOUBLE)

        addkeywordto_head(tab, "POLNGRP", "Number of context groups")
        fitswrite(head, "POLNGRP", psf.poly.ngroup, H_INT, T_LONG)
        for i in range(psf.poly.ngroup):
            str = "POLDEG%1d" % (i+1)
            addkeywordto_head(tab, str, "Polynom degree for this context group")
            fitswrite(head, str, psf.poly.degree[i], H_INT, T_LONG)
    
        addkeywordto_head(tab, "PSF_FWHM", "PSF FWHM in image pixels")
        if psf.samples_accepted:
            fitswrite(head, "PSF_FWHM", psf.fwhm, H_FLOAT, T_FLOAT)
        else:
            fitswrite(head, "PSF_FWHM", zero, H_FLOAT, T_FLOAT)
       
        addkeywordto_head(tab, "PSF_SAMP", "Sampling step of the PSF data in"
                " image pixels")
        
        if psf.samples_accepted:
            fitswrite(head, "PSF_SAMP", psf.pixstep,H_FLOAT, T_FLOAT)
        else:
            fitswrite(head, "PSF_SAMP", zero, H_FLOAT, T_FLOAT)
        
        addkeywordto_head(tab, "PSFNAXIS", "Dimensionality of the PSF data")
        fitswrite(head, "PSFNAXIS", psf.dim, H_INT, T_LONG)
        for i in range(psf.dim):
            str = "PSFAXIS%1d" % (i+1)
            addkeywordto_head(tab, str, "Number of element along this axis")
            fitswrite(head, str, psf.size[i], H_INT, T_LONG)
        key = new_key("PSF_MASK")
        key.naxis = psf.dim
        for i in range(psf.dim):
            key.naxisn[i] = psf.size[i]
        
        key.comment = "Tabulated PSF data"
        key.htype = H_FLOAT
        key.ttype = T_FLOAT
        key.nbytes = psf.npix*t_size[T_FLOAT]
        key.nobj = 1
        key.ptr = psf.comp
        add_key(key, tab, 0)

        if (psf.basiscoeff) :
            key = new_key("PSF_COEFFS")
            key.naxis = psf.dim - 1
            key.naxisn[0] = psf.nbasis
            if (key.naxis>1) :
                key.naxisn[1] = psf.size[2]
            
            key.comment = "PSF basis vector coefficients"
            key.htype = H_FLOAT
            key.ttype = T_FLOAT
            key.nbytes = psf.nbasis*psf.size[2]*t_size[T_FLOAT]
            key.nobj = 1
            key.ptr = psf.basiscoeff
            add_key(key, tab, 0)
        
        save_tab(cat, tab)
        blank_keys(tab)
        free_tab(tab)
    free_cat(cat , 1)

    return
