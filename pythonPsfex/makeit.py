import numpy as np
from define import *
from time import time
from datetime import date, datetime
thetime, thetime2 = 0.0, 0.0
from psf import *
from pca import *
from homo import *
from field import *
from xml import *
from prefs import *
from astropy.io.fits import *
from check import *
from catout import *
from sample import *
from context import *
from diagnostic import *


def  makeit():

    error_installfunc(write_error)

    incatnames = prefs.incat_name
    ncat = prefs.ncat

    thetime = time()
    today = date.today()
    now = datetime.now()
    prefs.sdate_start = today.strftime("%Y-%m-%d")
    prefs.stime_start = now.strftime("%H:%M:%S")

    NFPRINTF(OUTPUT, "")
    if prefs.nthreads>1:
        QPRINTF(OUTPUT,
                "----- %s %s started on %s at %s with %d thread%s\n\n",
                BANNER,
                MYVERSION,
                prefs.sdate_start,
                prefs.stime_start,
                prefs.nthreads,
                "s")
    else:
        QPRINTF(OUTPUT,
                "----- %s %s started on %s at %s with %d thread%s\n\n",
                BANNER,
                MYVERSION,
                prefs.sdate_start,
                prefs.stime_start,
                prefs.nthreads,
                "")


    if ( not ncat):
        today = date.today()
        now = datetime.now()
        prefs.sdate_start = today.strftime("%Y-%m-%d")
        prefs.stime_start = now.strftime("%H:%M:%S")

        if (prefs.xml_flag):
            init_xml(0)
            write_xml(prefs.xml_name)
            end_xml()
        
        return

    NFPRINTF(OUTPUT, "")
    QPRINTF(OUTPUT, "----- %d input catalogs:\n", ncat)
    for c in range(ncat):
        fields[c] = field_init(incatnames[c])
        if fields[c].ndet>1:
            if fields[c].next>1:
                QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
                        fields[c].rcatname, fields[c].ident,
                        fields[c].next, "s",
                        fields[c].ndet, "s")
            else:
                QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
                        fields[c].rcatname, fields[c].ident,
                        fields[c].next, "",
                        fields[c].ndet, "s")
        else:
            if fields[c].next>1:
                QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
                        fields[c].rcatname, fields[c].ident,
                        fields[c].next, "s",
                        fields[c].ndet, "")
            else:
                QPRINTF(OUTPUT, "%-20.20s:  \"%-16.16s\"  %3d extension%s %7d detection%s\n",
                        fields[c].rcatname, fields[c].ident,
                        fields[c].next, "",
                        fields[c].ndet, "")

                                
    
    QPRINTF(OUTPUT, "\n")

    next = fields[0].next

    if (prefs.xml_flag) :
        init_xml(ncat)
    

    psfstep = prefs.psf_step
    psfsteps = None
    nbasis = 0
    psfbasis = None
    psfbasiss = None

    NFPRINTF(OUTPUT, "Initializing contexts..."    )
    context = context_init(prefs.context_name, prefs.context_group,
            prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
            CONTEXT_REMOVEHIDDEN)
    if context.npc :
        fullcontext = context_init(prefs.context_name, prefs.context_group,
                    prefs.ncontext_group, prefs.group_deg, prefs.ngroup_deg,
                    CONTEXT_KEEPHIDDEN)    
    else:
        fullcontext = context

    if (context.npc and ncat<2):
        warning("Hidden dependencies cannot be derived from",
                " a single catalog");
    elif (context.npc and prefs.stability_type == STABILITY_EXPOSURE):
        warning("Hidden dependencies have no effect",
                " in STABILITY_TYPE EXPOSURE mode")
    
    if (not prefs.psf_step):
        NFPRINTF(OUTPUT, "Computing optimum PSF sampling steps...")
        if (prefs.newbasis_type==NEWBASIS_PCACOMMON
                or (prefs.stability_type == STABILITY_SEQUENCE
                        and prefs.psf_mef_type == PSF_MEF_COMMON)):
            set = load_samples(incatnames, 0, ncat, ALL_EXTENSIONS, next, context)
            psfstep = (float)((set.fwhm/2.35)*0.5)
            end_set(set)
        
        elif (prefs.newbasis_type == NEWBASIS_PCAINDEPENDENT
                or context.npc
                or (prefs.stability_type == STABILITY_SEQUENCE
                        and prefs.psf_mef_type == PSF_MEF_INDEPENDENT)):
            psfsteps = np.zeros(next, dtype=np.float32)
            for ext in range(next):
                set = load_samples(incatnames, 0, ncat, ext, next, context)
                if psfstep:
                    psfsteps[ext] = float(psfstep)
                else:
                    psfsteps[ext] = float((set.fwhm/2.35)*0.5)

                end_set(set)
            
    if (prefs.newbasis_type==NEWBASIS_PCACOMMON):
        for ext in range(next):
            for c in range(ncat):
                str = "Computing new PCA image basis from %s..." %(fields[c].rtcatname)
                NFPRINTF(OUTPUT, str)
                set = load_samples(incatnames, c, 1, ext, next, context)
                step = psfstep
                cpsf[c+ext*ncat] = make_psf(set, psfstep, NULL, 0, context)
                end_set(set)
            
        nbasis = prefs.newbasis_number
        psfbasis = pca_onsnaps(cpsf, ncat*next, nbasis)
        for i in range(ncat*next):
            psf_end(cpsf[i])

    elif (prefs.newbasis_type == NEWBASIS_PCAINDEPENDENT):
        nbasis = prefs.newbasis_number
        psfbasiss = np.zeros(next, dtype=np.float32)
        for ext in range(next):
            if (psfsteps):
                step = psfsteps[ext]
            else:
                step = psfstep
            
            for c in range(ncat):
                str = "Computing new PCA image basis from %s..." % (fields[c].rtcatname)
                NFPRINTF(OUTPUT, str)
                set = load_samples(incatnames, c, 1, ext, next, context)
                cpsf[c] = make_psf(set, step, NULL, 0, context)
                end_set(set)
            
            psfbasiss[ext] = pca_onsnaps(cpsf, ncat, nbasis)
            for c in range(ncat):
                psf_end(cpsf[c])

    if (context.npc and prefs.hidden_mef_type == HIDDEN_MEF_COMMON):
        for c in range(ncat):
            str = "Computing hidden dependency parameter(s) from %s..." % (fields[c].rtcatname)
            NFPRINTF(OUTPUT, str)
            for ext in range(next):
                set = load_samples(incatnames, c, 1, ext, next, context)
                if (psfsteps):
                    step = psfsteps[ext]
                else:
                    step = psfstep
                
                if psfbasiss:
                    basis = psfbasiss[ext]
                else:
                    basis = psfbasis
                cpsf[ext+c*next] = make_psf(set, step, basis, nbasis, context)
                end_set(set)
        
    
        fullcontext.pc = pca_oncomps(cpsf, next, ncat, context.npc);
        for c in range(ncat*next):
            psf_end(cpsf[c])
        
    if (prefs.psf_mef_type == PSF_MEF_COMMON):
        if (prefs.stability_type == STABILITY_SEQUENCE):
            set = load_samples(incatnames, 0, ncat, ALL_EXTENSIONS, next, context)
            step = psfstep
            basis = psfbasis
            field_count(fields, set, COUNT_LOADED)
            psf = make_psf(set, step, basis, nbasis, context)
            field_count(fields, set, COUNT_ACCEPTED)
            end_set(set)
            NFPRINTF(OUTPUT, "Computing final PSF model...")
            context_apply(context, psf, fields, ALL_EXTENSIONS, 0, ncat)
            psf_end(psf)
        else:
            for c in range(ncat):
                str = "Computing final PSF model from %s..." % (fields[c].rtcatname)
                NFPRINTF(OUTPUT, str)
                set = load_samples(incatnames, c, 1, ALL_EXTENSIONS, next, context)
                if (psfstep):
                    step = psfstep
                else:
                    step = float((set.fwhm/2.35)*0.5)
                
                basis = psfbasis
                field_count(fields, set, COUNT_LOADED)
                psf = make_psf(set, step, basis, nbasis, fullcontext)
                field_count(fields, set, COUNT_ACCEPTED)
                end_set(set)
                context_apply(fullcontext, psf, fields, ALL_EXTENSIONS, c, 1)
                psf_end(psf)
        
    else :
        for ext in range(next):
            if psfbasiss:
                basis = psfbasiss[ext]
            else:
                basis = psfbasis

            if (context.npc and prefs.hidden_mef_type == HIDDEN_MEF_INDEPENDENT):

                if (psfsteps) :
                    step = psfsteps[ext]
                else :
                    step = psfstep
                
                for c in range(ncat):
                    if (next>1):
                        str = "Computing hidden dependency parameter(s) from %s[%d/%d]..." % (fields[c].rtcatname, ext+1, next)
                    else :
                        str, "Computing hidden dependency parameter(s) from %s..." % (fields[c].rtcatname)
                    
                    NFPRINTF(OUTPUT, str)
                    set = load_samples(incatnames, c, 1, ext, next, context)
                    field_count(fields, set, COUNT_LOADED)
                    cpsf[c] = make_psf(set, step, basis, nbasis, context)
                    field_count(fields, set, COUNT_ACCEPTED)
                    end_set(set)
                
                fullcontext.pc = pca_oncomps(cpsf, 1, ncat, context.npc)
                for c in range(ncat):
                    psf_end(cpsf[c])

            if (prefs.stability_type == STABILITY_SEQUENCE):
                if (next>1) :
                    str = "Computing final PSF model for extension [%d/%d]..." % (ext+1, next)
                    NFPRINTF(OUTPUT, str)
                else:
                    NFPRINTF(OUTPUT, "Computing final PSF model...")
                
                set = load_samples(incatnames, 0, ncat, ext, next, fullcontext)
                if (psfstep):
                    step = psfstep
                elif (psfsteps):
                    step = psfsteps[ext]
                else :
                    step = float((set.fwhm/2.35)*0.5)
                
                field_count(fields, set, COUNT_LOADED)
                psf = make_psf(set, step, basis, nbasis, fullcontext)
                field_count(fields, set, COUNT_ACCEPTED)
                end_set(set)
                context_apply(fullcontext, psf, fields, ext, 0, ncat)
                psf_end(psf)
            else:
                for c in range(ncat):
                    if (next>1):
                        str = "Reading data from %s[%d/%d]..." %  (fields[c].rtcatname, ext+1, next)
                    else:
                        str = "Reading data from %s..." % (fields[c].rtcatname)
                    
                    NFPRINTF(OUTPUT, str)
                    set = load_samples(incatnames, c, 1, ext, next, context)
                    if (psfstep):
                        step = psfstep
                    elif (psfsteps):
                        step = psfsteps[ext]
                    else :
                        step = float((set.fwhm/2.35)*0.5)
                    
                    if (next>1):
                        str = "Computing final PSF model for %s[%d/%d]..." % (fields[c].rtcatname, ext+1, next)
                    else :
                        str = "Computing final PSF model for %s..." % (fields[c].rtcatname)
                    
                    NFPRINTF(OUTPUT, str)
                    field_count(fields, set, COUNT_LOADED)
                    psf = make_psf(set, step, basis, nbasis, context)
                    field_count(fields, set, COUNT_ACCEPTED)
                    end_set(set)
                    context_apply(context, psf, fields, ext, c, 1)
                    psf_end(psf)
                    

    QIPRINTF(OUTPUT,
            "   filename      [ext] accepted/total samp. chi2/dof FWHM ellip."
            " resi. asym.")

    if (prefs.outcat_type != CAT_NONE):
        outcat = init_outcat(prefs.outcat_name, context.ncontext)
    
    for c in range(ncat):
        field = fields[c]
        set = None
        if next >1 : 
            bigset = init_set(context)
        else:
            bigset = None

        for ext in range(next):
            psf = field.psf[ext]
            wcs = field.wcs[ext]
            if (next>1):
                str = "Computing diagnostics for %s[%d/%d]..." % (field.rtcatname, ext+1, next)
            else :
                str = "Computing diagnostics for %s..." % (field.rtcatname)
            
            NFPRINTF(OUTPUT, str)
            set = load_samples(incatnames, c, 1, ext, next, context)
            psf.samples_loaded = set.ngood
            if (set.nsample>1):
                psf_clean(psf, set, prefs.prof_accuracy)
                if set.nsample:
                    psf.chi2 = psf_chi2(psf, set)
                else:
                    psf.chi2 = 0.0
            
            psf.samples_accepted = set.ngood
            psf_diagnostic(psf)
            psf_wcsdiagnostic(psf, wcs)
            nmed = psf.nmed
            field_stats(fields, set)
            if (next>1):
                str = "[%d/%d]" % (ext+1, next)
            else:
                str[0] = '\0'
            
            if ext == 0:
                QPRINTF(OUTPUT, "%-17.17s%-7.7s %5d/%-5d %6.2f %6.2f %6.2f  %4.2f"
                    " %5.2f %5.2f\n",
                    field.rtcatname,
                    str,
                    psf.samples_accepted, psf.samples_loaded,
                    psf.pixstep,
                    psf.chi2,
                    psf.moffat_fwhm,
                    psf.moffat_ellipticity,
                    psf.pfmoffat_residuals,
                    psf.sym_residuals)
            else:
                QPRINTF(OUTPUT, "%-17.17s%-7.7s %5d/%-5d %6.2f %6.2f %6.2f  %4.2f"
                    " %5.2f %5.2f\n",
                    "",
                    str,
                    psf.samples_accepted, psf.samples_loaded,
                    psf.pixstep,
                    psf.chi2,
                    psf.moffat_fwhm,
                    psf.moffat_ellipticity,
                    psf.pfmoffat_residuals,
                    psf.sym_residuals)
                        
            for i in range(prefs.ncheck_type):
                if (prefs.check_type[i]):
                    str = "Saving CHECK-image #%d..." % (i+1)
                    NFPRINTF(OUTPUT, str)
                    check_write(field, set, prefs.check_name[i], prefs.check_type[i], ext, next, prefs.check_cubeflag)

            if (prefs.outcat_type != CAT_NONE):
                write_outcat(outcat, set)
            
            if (bigset) :
                add_set(bigset, set);
                end_set(set);
        
        if bigset:
            end_set(bigset)
        else:
            end_set(set)    

    if (prefs.outcat_type != CAT_NONE):
        end_outcat(outcat)
    
    for c in range(ncat):
        str = "Saving PSF model and metadata for %s..." % (fields[c].rtcatname)
        NFPRINTF(OUTPUT, str)
        if (prefs.psf_dir):
            pstr = find(incatnames[c], '/')
            if (pstr>-1):
                pstr+=1
            else:
                pstr = incatnames[c]
            
            str = "%s/%s" % (prefs.psf_dir, pstr)
        else:
            str = incatnames[c]
        
        pstr = find(str, '.')
        if (pstr >-1):
            pstr = str+len(str)
        
        pstr = "%s" % (prefs.psf_suffix)
        field_psfsave(fields[c], str)
        if (prefs.homobasis_type != HOMOBASIS_NONE):
            for ext in range(next):
                if (next>1):
                    str = "Computing homogenisation kernel for %s[%d/%d]..." % (fields[c].rtcatname, ext+1, next)
                else:
                    str =  "Computing homogenisation kernel for %s..." % (fields[c].rtcatname)
                
                NFPRINTF(OUTPUT, str)
                if (prefs.homokernel_dir):
                    pstr = find(incatnames[c], '/')
                    if (pstr>-1):
                        pstr+=1
                    else:
                        pstr = incatnames[c]
                    
                    str = "%s/%s" % (prefs.homokernel_dir, pstr)
                else:
                    str = incatnames[c]
                pstr = find(str, '.')
                if (pstr > -1):
                    pstr = str+len(str)
                
                pstr = "%s" % (prefs.homokernel_suffix)
                psf_homo(fields[c].psf[ext], str, prefs.homopsf_params, prefs.homobasis_number, prefs.homobasis_scale, ext, next);

        if (prefs.xml_flag):
            update_xml(fields[c])

    thetime2 = time()
    today = date.today()
    now = datetime.now()
    prefs.sdate_start = today.strftime("%Y-%m-%d")
    prefs.stime_start = now.strftime("%H:%M:%S")

    prefs.time_diff = thetime2 - thetime

    if (prefs.xml_flag) :
        NFPRINTF(OUTPUT, "Writing XML file...")
        write_xml(prefs.xml_name)
        end_xml()
    
    for c in range(ncat):
        field_end(fields[c])
    
    if (context.npc):
        context_end(fullcontext)
    
    return



def make_psf(set, psfstep, basis, nbasis, context):
    pixsize[0] = prefs.psf_pixsize[0]
    pixsize[1] = prefs.psf_pixsize[1]
    psf = psf_init(context, prefs.psf_size, psfstep, pixsize, set.nsample)

    psf.samples_total = set.nsample
    psf.samples_loaded = set.ngood
    psf.fwhm = set.fwhm
    psf_make(psf, set, 0.2)
    if (basis and nbasis):
        basis = psf.basis.copy()
        psf.nbasis = nbasis
    else :
        basistype = prefs.basis_type;
        if (basistype==BASIS_PIXEL_AUTO):
            if psf.fwhm < PSF_AUTO_FWHM:
                basistype = BASIS_NONE
            else:
                basistype = BASIS_PIXEL           
        
        psf_makebasis(psf, set, basistype, prefs.basis_number)
    
    psf_refine(psf, set)

    if (set.ngood>1):
        psf_clean(psf, set, 0.2)

        psf_make(psf, set, 0.1)
        psf_refine(psf, set)
    
    if (set.ngood>1) :
        psf_clean(psf, set, 0.1)

        psf_make(psf, set, 0.05)
        psf_refine(psf, set)
    

    if (set.ngood>1):
        psf_clean(psf, set, 0.05)
    
    psf.samples_accepted = set.ngood

    psf_make(psf, set, prefs.prof_accuracy)
    psf_refine(psf, set)

    psf_clip(psf)

    if set.ngood:
        psf.chi2 = psf_chi2(psf, set)
    else:
        psf.chi2 = 0.0
    
    return psf



def write_error(msg1, msg2):
    error = "%s%s" % (msg1,msg2)
    if (prefs.xml_flag):
        write_xmlerror(prefs.xml_name, error)
    
    end_xml()

    return
