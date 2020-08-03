import numpy as np
from define import *
from math import cos, sin, acos, asin, sqrt
from random import random
from enum import Enum
from astropy.wcs.docstrings import naxis

NAXIS = 2
DEG = (PI/180.0)
ARCMIN = (DEG/60.0)
ARCSEC = (DEG/3600.0)
MAS = (ARCSEC/1000.0)
MINUTE = 60.0
HOUR = (60.0*MINUTE)
DAY = (24.0*HOUR)
YEAR = (365.25*DAY)
MJD2000 = 51544.50000
JD2000 = (2400000.5+MJD2000)
MJD1950 = 33281.92346
JD1950 = (2400000.5+MJD1950)
JU2TROP = 1.0000214
WCS_NOCOORD = 1e31
WCS_NGRIDPOINTS = 12
WCS_NGRIDPOINTS2 = (WCS_NGRIDPOINTS*WCS_NGRIDPOINTS)
WCS_INVMAXDEG = 9
WCS_INVACCURACY = 0.001
WCS_NRANGEPOINTS = 32

class celsys_enum(Enum):
    CELSYS_NATIVE = 1
    CELSYS_PIXEL = 2
    CELSYS_EQUATORIAL = 3
    CELSYS_GALACTIC = 4
    CELSYS_ECLIPTIC = 5
    CELSYS_SUPERGALACTIC = 6
    
class radecsys_enum(Enum):
    RDSYS_ICRS = 1
    RDSYS_FK5 = 2
    RDSYS_FK4 = 3
    RDSYS_FK4_NO_E = 4
    RDSYS_GAPPT = 5
    
class wcs_struct():
    def __init__(self, naxis, naxisn, ctype, cunit, crval, cdelt, crpix,
                crder, csyer, cd, projp, nprojp, longpole, latpole,wcsmin, 
                wcsmax, wcsscale, wcsscalepos, wcsmaxradius, outmin, outmax,
                lat, lng, r0, lindet, chirality, pixscale, ap2000, dp2000,
                ap1950, dp1950, obsdate, equinox, epoch, radecsys, celsys,
                celsysmat, celsysconvflag, wcsprm, lin, cel, prj, tnx_latcor,
                tnx_lngcor, inv_x, inv_y):
        self.naxis = naxis
        self.naxisn = naxisn
        self.ctype = ctype
        self.cunit = cunit
        self.crval = crval
        self.cdelt = cdelt
        self.crpix = crpix
        self.crder = crder
        self.csyer = csyer
        self.cd = cd
        self.projp = projp
        self.nprojp = nprojp
        self.longpole = longpole
        self.latpole = latpole
        self.wcsmin = wcsmin
        self.wcsmax = wcsmax
        self.wcsscale = wcsscale
        self.wcsscalepos = wcsscalepos
        self.wcsmaxradius = wcsmaxradius
        self.outmin = outmin
        self.outmax = outmax
        self.lat = lat
        self.lng = lng
        self.r0 = r0
        self.lindet = lindet
        self.chirality = chirality
        self.pixscale = pixscale
        self.ap2000 = ap2000
        self.dp2000 = dp2000
        self.ap1950 = ap1950
        self.dp1950 = dp1950
        self.obsdate = obsdate
        self.equinox = equinox
        self.epoch = epoch
        self.radecsys = radecsys
        self.celsys = celsys
        self.celsysmat = celsysmat
        self.celsysconvflag = celsysconvflag
        self.wcsprm = wcsprm
        self.lin = lin
        self.cel = cel
        self.prj = prj
        self.tnx_latcor = tnx_latcor
        self.tnx_lngcor = tnx_lngcor
        self.inv_x = inv_x
        self.inv_y = inv_y


def copy_wcs(wcsin):
    wcsin = wcs.copy()
    wcsin.projp = wcs.projp.copy()

    wcs.wcsprm = None
    wcs.lin = None
    wcs.cel = None
    wcs.prj = None
    wcs.tnx_lngcor = copy_tnxaxis(wcsin.tnx_lngcor)
    wcs.tnx_latcor = copy_tnxaxis(wcsin.tnx_latcor)
    wcs.inv_x = wcs.inv_y = None

    wcsset(wcs.naxis,wcs.ctype, wcs.wcsprm)

    init_wcs(wcs)
    invert_wcs(wcs)
    range_wcs(wcs)

    return wcs


def create_wcs(ctype,crval, crpix, cdelt, naxisn, naxis):
    #    QCALLOC(wcs, wcsstruct, 1);

    wcs.naxis = naxis
    wcs.projp = np.zeros(naxis*100, dtype=np.float64)

    wcs.nprojp = 0

    wcs.longpole = wcs.latpole = 999.0
    for l in range(naxis):
        if naxisn:
            wcs.naxisn[l] = naxisn[l]
        else:
            wcs.naxisn[l] = 360.0
        if (ctype) :
            wcs.ctype[l] = ctype[l]
        elif (l==0):
            wcs.ctype[l] = "RA---AIT"
        elif (l==1):
            wcs.ctype[l] = "DEC--AIT"
        
        if crval:
            wcs.crval[l] = crval[l]
        else:
            wcs.crval[l] = 0.0
            
        if crpix:
            wcs.crpix[l] = crpix[l]
        else:
            wcs.crpix[l] = 0.0
        wcs.cdelt[l] = 1.0
        if cdelt:
            wcs.cd[l*(naxis+1)] = cdelt[l] 
        else:
            wcs.cd[l*(naxis+1)] = 1.0
    

    wcs.epoch = wcs.equinox = 2000.0
#    QCALLOC(wcs->wcsprm, struct wcsprm, 1);
    wcsset(wcs.naxis,wcs.ctype, wcs.wcsprm)
    init_wcs(wcs)
    invert_wcs(wcs)
    range_wcs(wcs)

    return wcs



def init_wcs(wcs):
    naxis = wcs.naxis
#    QCALLOC(wcs->lin, struct linprm, 1);
    wcs.lin.cdelt = np.zeros(naxis, dtype=np.float64)
    wcs.lin.crpix = np.zeros(naxis, dtype=np.float64)
    wcs.lin.pc = np.zeros(naxis**2, dtype=np.float64)
#    QCALLOC(wcs->cel, struct celprm, 1);
#    QCALLOC(wcs->prj, struct prjprm, 1);
    if (wcs.inv_x):
        poly_end(wcs.inv_x)
        wcs.inv_x = None
    
    if (wcs.inv_y):
        poly_end(wcs.inv_y)
        wcs.inv_y = None
    
    wcs.lin.flag = wcs.cel.flag = wcs.prj.flag = 0
    wcs.lin.naxis = naxis

    lng = wcs.lng = wcs.wcsprm.lng
    lat = wcs.lat = wcs.wcsprm.lat

    for l in range(naxis):
        wcs.lin.crpix[l] = wcs.crpix[l]
        wcs.lin.cdelt[l] = 1.0

    for l in range(naxis*naxis):
        wcs.lin.pc[l] = wcs.cd[l]

    if (lng>=0):
        wcs.cel.ref[0] = wcs.crval[lng]
        wcs.cel.ref[1] = wcs.crval[lat]
    else:
        wcs.cel.ref[0] = wcs.crval[0]
        wcs.cel.ref[1] = wcs.crval[1]
    
    wcs.cel.ref[2] = wcs.longpole
    wcs.cel.ref[3] = wcs.latpole

    wcs.prj.r0 = wcs.r0
    wcs.prj.tnx_lngcor = wcs.tnx_lngcor
    wcs.prj.tnx_latcor = wcs.tnx_latcor
    if (lng>=0):
        n = 0
        for l in range(100, -1,-1):
            wcs.prj.p[l] = wcs.projp[l+lat*100]
            wcs.prj.p[l+100] = wcs.projp[l+lng*100]
            if (not n and (wcs.prj.p[l] or wcs.prj.p[l+100])):
                n = l+1
            
        wcs.nprojp = n

    wcs.chirality = wcs_chirality(wcs)

    init_wcscelsys(wcs)

    return


def init_wcscelsys(wcs):
    lng = wcs.wcsprm.lng
    lat = wcs.wcsprm.lat
    if (lng==lat):
        wcs.celsysconvflag = 0
        return
    
    s = 0
    while(celsysname[s][0][0] and wcs.ctype[lng] == celsysname[s][0]):
        s+=1

    if (not s or not celsysname[s][0][0]):
        wcs.celsysconvflag = 0
        return
    
    wcs.celsys = s
    a0 = celsysorig[s][0]*DEG
    d0 = celsysorig[s][1]*DEG
    ap = celsyspole[s][0]*DEG
    dp = celsyspole[s][1]*DEG
    y = sin(ap - a0)

    ap2 = asin(cos(d0)*y)
    mat = wcs.celsysmat
    mat[0] = ap
    mat[1] = ap2
    mat[2] = cos(dp)
    mat[3] = sin(dp)

    wcs.celsysconvflag = 1
    return

def strpbrk(a, b):
    pos = -1
    for char in b:
        pos = max(pos, find(a, char))
    return pos

def read_wcs(tab):
    def  FITSREADF(buf, k, val, def_var):
        if (fitsread(buf,k, val, H_FLOAT,T_DOUBLE) != RETURN_OK):
            val = def_var
        
    def FITSREADI(buf, k, val, def_var):
        if (fitsread(buf,k, val, H_INT,T_LONG) != RETURN_OK):
            val = def_var
        
    def FITSREADS(buf, k, str, def_var):
        if (fitsread(buf,k,str, H_STRING,T_STRING) != RETURN_OK):
            str = def_var
    
    buf = tab.headbuf
    if tab.cat:
        filename = (tab.cat.filename)
    else:
         name = "internal header"
         filename = name

    FITSREADS(buf, "OBJECT  ", str, "Unnamed")
#    QCALLOC(wcs, wcsstruct, 1);
    if (tab.naxis > NAXIS):
        warning("Maximum number of dimensions supported by this version of the ", "software exceeded\n")
        tab.naxis = 2

    wcs.naxis = naxis = tab.naxis
    wcs.projp = np.zeros(naxis*100, np.float64)
    for l in range(naxis):
        wcs.naxisn[l] = tab.naxisn[l]
        str = "CTYPE%-3d"
        FITSREADS(buf, str, str, "")
        wcs.ctype[l] = str
        str = "CUNIT%-3d" % (l+1)
        FITSREADS(buf, str, str, "deg")
        wcs.cunit[l] = str
        str = "CRVAL%-3d" % (l+1)
        FITSREADF(buf, str, wcs.crval[l], 0.0)
        str = "CRPIX%-3d" % (l+1)
        FITSREADF(buf, str, wcs.crpix[l], 1.0)
        str = "CDELT%-3d" % (l+1)
        FITSREADF(buf, str, wcs.cdelt[l], 1.0)
        str = "CRDER%-3d" % (l+1)
        FITSREADF(buf, str, wcs.crder[l], 0.0)
        str = "CSYER%-3d" % (l+1)
        FITSREADF(buf, str, wcs.csyer[l], 0.0)
        if (fabs(wcs.cdelt[l]) < 1e-30):
            error(EXIT_FAILURE, "*Error*: CDELT parameters out of range in ",
                    filename)
        
    if (fitsfind(buf, "CD?_????")!=RETURN_ERROR):
        for l in range(naxis):
            for j in range(naxis):
                str = "CD%d_%d" % (l+1, j+1)
                if l==j:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 1.0)
                else:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 0.0)

    elif (fitsfind(buf, "PC?_????")!=RETURN_ERROR):
        for l in range(naxis):
            for j in range(naxis):
                str = "PC%d_%d" % (l+1, j+1)
                if l==j:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 1.0)
                else:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 0.0)
                
                wcs.cd[l*naxis+j] *= wcs.cdelt[l]
            
    elif (fitsfind(buf, "PC0??0??")!=RETURN_ERROR):
        for l in range(naxis):
            for j in range(naxis):
                str = "PC%03d%03d" % (l+1, j+1)
                if l==j:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 1.0)
                else:
                    FITSREADF(buf, str, wcs.cd[l*naxis+j], 0.0)

                wcs.cd[l*naxis+j] *= wcs.cdelt[l]
            
    else:
        FITSREADF(buf, "CROTA2  ", drota, 0.0)
        wcs.cd[2] = sin(drota*DEG)
        wcs.cd[1] = -wcs.cd[2]
        wcs.cd[0] *= wcs.cdelt[0]
        wcs.cd[2] *= wcs.cdelt[0]
        wcs.cd[1] *= wcs.cdelt[1]
        wcs.cd[3] *= wcs.cdelt[1]
#    QCALLOC(wcs->wcsprm, struct wcsprm, 1);
    
    if (not wcsset(wcs.naxis,wcs.ctype, wcs.wcsprm) and wcs.wcsprm.flag<999):

        FITSREADF(buf, "MJD-OBS ", date, -1.0)
        if (date<0.0):
            FITSREADF(buf, "MJDSTART", date, -1.0)
        
        if (date>0.0):
            wcs.obsdate = 2000.0 - (MJD2000 - date)/365.25
        else :
            FITSREADS(buf, "DATE-OBS", str, "")
            if (str):
                for l in range(3):
                    if l:
                        pstr = strtok_r(None,"/- ", ptr) 
                    else:
                        pstr = strtok_r(str,"/- ", ptr)
                    if not pstr:
                        break
                    dpar[l] = atoi(pstr)
                
                if (l<3 or not dpar[0] or not dpar[1] or not dpar[2]):
                    warning("Invalid DATE-OBS value in header: ", str)
                    dpar[0] = 2000
                    dpar[1] = 1
                    dpar[2] = 1
                elif (find(str, '/')>-1 and dpar[0]<32 and dpar[2]<100):
                    j = dpar[0]
                    dpar[0] = dpar[2]+1900
                    dpar[2] = j

                if (dpar[0]%4):
                    biss = 0
                else:
                    biss = 1
                if dpar[1]>2:
                    date = -678956 + (365*dpar[0]+dpar[0]/4) - biss + ((((dpar[1]+1)*30.6)-63+biss) + dpar[2])
                else:
                    date = -678956 + (365*dpar[0]+dpar[0]/4) - biss + ((((dpar[1]-1)*(63+biss))/2) + dpar[2])
                                                                          
                wcs.obsdate = 2000.0 - (MJD2000 - date)/365.25
            else:
                wcs.obsdate = 0.0

        FITSREADF(buf, "EPOCH", wcs.epoch, 2000.0)
        FITSREADF(buf, "EQUINOX", wcs.equinox, wcs.epoch)
        if (fitsread(buf, "RADESYS", str, H_STRING,T_STRING) != RETURN_OK):
            if wcs.equinox >= 2000.0:
                FITSREADS(buf, "RADECSYS", str, "ICRS")
            elif wcs.equinox<1984.0:
                FITSREADS(buf, "RADECSYS", str, "FK4")
            else:
                FITSREADS(buf, "RADECSYS", str, "FK5")
        if (str == "ICRS"):
            wcs.radecsys = RDSYS_ICRS
        elif (str == "FK5"):
            wcs.radecsys = RDSYS_FK5
        elif (str == "FK4"):
            if (wcs.equinox == 2000.0):
                FITSREADF(buf, "EPOCH  ", wcs.equinox, 1950.0)
                FITSREADF(buf, "EQUINOX", wcs.equinox, wcs.equinox)
            
            wcs.radecsys = RDSYS_FK4
            warning("FK4 precession formulae not yet implemented:\n","            Astrometry may be slightly inaccurate")
        elif (str == "FK4-NO-E"):
            if (wcs.equinox == 2000.0):
                FITSREADF(buf, "EPOCH", wcs.equinox, 1950.0)
                FITSREADF(buf, "EQUINOX", wcs.equinox, wcs.equinox)
            
            wcs.radecsys = RDSYS_FK4_NO_E
            warning("FK4 precession formulae not yet implemented:\n","            Astrometry may be slightly inaccurate")
        elif (str == "GAPPT"):
            wcs.radecsys = RDSYS_GAPPT
            warning("GAPPT reference frame not yet implemented:\n","            Astrometry may be slightly inaccurate")
        else:
            warning("Using ICRS instead of unknown astrometric reference frame: ",str)
            wcs.radecsys = RDSYS_ICRS
        
        if (wcs.wcsprm.pcode == "TNX"):
            if (fitsfind(buf, "WAT?????") != RETURN_ERROR):
                pstr = wstr1
                str = "WAT1_001"
                j=2
                while(fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK):
                    str = "WAT1_%03d" % (j)
                    pstr += len(pstr)
                    j+=1
                
                pstr = wstr2
                str = "WAT2_001"
                j=2
                while (fitsread(buf,str,pstr,H_STRINGS,T_STRING)==RETURN_OK):
                    str = "WAT2_%03d" % (j)
                    pstr += len(pstr)
                    j+=1
                
                pstr = max(find(wstr1, "longpole"), find(wstr2, "longpole"))
            
                if (pstr > -1):
                    pstr = strpbrk(pstr, "1234567890-+.");
                
                if pstr:
                    wcs.longpole = atof(pstr)
                else:
                    wcs.longpole = 999.0
                wcs.latpole = 999.0
                pstr = max(find(wstr1, "ro"), find(wstr2, "ro"))
                if pstr:
                    pstr = strpbrk(pstr, "1234567890-+.");
                if pstr:
                    wcs.r0 = atof(pstr)
                else:
                    wcs.r0 = 0.0
                pstr = max(find(wstr1, "lngcor"), find(wstr2, "lngcor"))
                if pstr:
                    wcs.tnx_lngcor = read_tnxaxis(pstr)
                
                if (not wcs.tnx_lngcor):
                    error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ", filename)
                
                pstr = max(find(wstr1, "latcor"), find(wstr2, "latcor"))
                if pstr:
                    wcs.tnx_latcor = read_tnxaxis(pstr)
                
                if (not wcs.tnx_latcor):
                    error(EXIT_FAILURE, "*Error*: incorrect TNX parameters in ", filename)
                
        else :
            if (fitsread(buf, "LONPOLE",wcs.longpole,H_FLOAT,T_DOUBLE) != RETURN_OK):
                FITSREADF(buf, "LONGPOLE", wcs.longpole, 999.0)
            
            FITSREADF(buf, "LATPOLE ", wcs.latpole, 999.0)
            if (fitsfind(buf, "PROJP???") != RETURN_ERROR):
                for j in range(10):
                    str = "PROJP%-3d" % (j);
                    FITSREADF(buf, str, wcs.projp[j], 0.0)
                
            if (fitsfind(buf, "PV?_????") != RETURN_ERROR):
                for l in range(naxis):
                    for j in range(100):
                        str = "PV%d_%d " % (l+1, j)
                        FITSREADF(buf, str, wcs.projp[j+l*100], 0.0)

    init_wcs(wcs)

    range_wcs(wcs)
    
    invert_wcs(wcs)

    return wcs


def write_wcs(tab, wcs):
    naxis = wcs.naxis
    addkeywordto_head(tab, "BITPIX  ", "Bits per pixel")
    fitswrite(tab.headbuf, "BITPIX  ", tab.bitpix, H_INT, T_LONG)
    addkeywordto_head(tab, "NAXIS   ", "Number of axes")
    fitswrite(tab.headbuf, "NAXIS   ", wcs.naxis, H_INT, T_LONG)
    for l in range(naxis):
        str = "NAXIS%-3d" % (l+1)
        addkeywordto_head(tab, str, "Number of pixels along this axis")
        fitswrite(tab.headbuf, str, wcs.naxisn[l], H_INT, T_LONG)
    
    addkeywordto_head(tab, "EQUINOX ", "Mean equinox")
    fitswrite(tab.headbuf, "EQUINOX ", wcs.equinox, H_FLOAT, T_DOUBLE)
    if (wcs.obsdate!=0.0):
        mjd = (wcs.obsdate-2000.0)*365.25 + MJD2000
        addkeywordto_head(tab, "MJD-OBS ", "Modified Julian date at start")
        fitswrite(tab.headbuf, "MJD-OBS ", mjd, H_EXPO,T_DOUBLE)
    
    addkeywordto_head(tab, "RADESYS ", "Astrometric system")
    if wcs.radecsys == RDSYS_ICRS:
        fitswrite(tab.headbuf, "RADESYS ", "ICRS", H_STRING, T_STRING)
        break
    elif wcs.radecsys == RDSYS_FK5:
        fitswrite(tab.headbuf, "RADESYS ", "FK5", H_STRING, T_STRING)
        break
    elif wcs.radecsys == RDSYS_FK4:
        fitswrite(tab.headbuf, "RADESYS ", "FK4", H_STRING, T_STRING)
        break
    elif wcs.radecsys == RDSYS_FK4_NO_E:
        fitswrite(tab.headbuf, "RADESYS ", "FK4-NO-E", H_STRING, T_STRING)
        break
    elif wcs.radecsys == RDSYS_GAPPT:
        fitswrite(tab.headbuf, "RADESYS ", "GAPPT", H_STRING, T_STRING)
        break
    else:
        error(EXIT_FAILURE, "*Error*: unknown RADESYS type in write_wcs()", "")
    
    for l in range(naxis):
        str = "CTYPE%-3d" % (l+1)
        addkeywordto_head(tab, str, "WCS projection type for this axis")
        fitswrite(tab.headbuf, str, wcs.ctype[l], H_STRING, T_STRING)
        str = "CUNIT%-3d" % (l+1)
        addkeywordto_head(tab, str, "Axis unit");
        fitswrite(tab.headbuf, str, wcs.cunit[l], H_STRING, T_STRING)
        str = "CRVAL%-3d" % (l+1)
        addkeywordto_head(tab, str, "World coordinate on this axis")
        fitswrite(tab.headbuf, str, wcs.crval[l], H_EXPO, T_DOUBLE)
        str = "CRPIX%-3d" % (l+1)
        addkeywordto_head(tab, str, "Reference pixel on this axis");
        fitswrite(tab.headbuf, str, wcs.crpix[l], H_EXPO, T_DOUBLE)
        for j in range(naxis):
            str = "CD%d_%d" % (l+1, j+1)
            addkeywordto_head(tab, str, "Linear projection matrix")
            fitswrite(tab.headbuf, str, wcs.cd[l*naxis+j], H_EXPO, T_DOUBLE)
        
        for j in range(100):
            if (abs(wcs.projp[j+100*l]) > TINY):
                str = "PV%d_%d" % (l+1, j)
                addkeywordto_head(tab, str, "Projection distortion parameter")
                fitswrite(tab.headbuf, str, wcs.projp[j+100*l], H_EXPO, T_DOUBLE)
            
    readbasic_head(tab)

    return


def wipe_wcs(tab):
    removekeywordfrom_head(tab, "CRVAL???")
    removekeywordfrom_head(tab, "CTYPE???")
    removekeywordfrom_head(tab, "CUNIT???")
    removekeywordfrom_head(tab, "CRPIX???")
    removekeywordfrom_head(tab, "CRDER???")
    removekeywordfrom_head(tab, "CSYER???")
    removekeywordfrom_head(tab, "CDELT???")
    removekeywordfrom_head(tab, "CROTA???")
    removekeywordfrom_head(tab, "CD?_????")
    removekeywordfrom_head(tab, "PROJP_??")
    removekeywordfrom_head(tab, "PV?_????")
    removekeywordfrom_head(tab, "PC?_????")
    removekeywordfrom_head(tab, "PC0??0??")
    removekeywordfrom_head(tab, "EQUINOX?")
    removekeywordfrom_head(tab, "RADESYS?")
    removekeywordfrom_head(tab, "RADECSYS")
    removekeywordfrom_head(tab, "LONPOLE?")
    removekeywordfrom_head(tab, "LONGPOLE")
    removekeywordfrom_head(tab, "LATPOLE?")
    removekeywordfrom_head(tab, "WAT?????")

    return



def end_wcs(wcs):
    if (wcs):
        free_tnxaxis(wcs.tnx_lngcor)
        free_tnxaxis(wcs.tnx_latcor)
        poly_end(wcs.inv_x)
        poly_end(wcs.inv_y)

    return


def wcs_supproj(name):
    projcode = ["TAN", "TPV", "AZP", "SIN", "STG", "ARC", "ZPN", "ZEA", "AIR", "CYP",
            "CAR", "MER", "CEA", "COP", "COD", "COE", "COO", "BON", "PCO", "GLS",
            "PAR", "AIT", "MOL", "CSC", "QSC", "TSC", "TNX", "NONE"]


    for code in projcode:
        if (name == code):
            return RETURN_OK
        
    return RETURN_ERROR



def invert_wcs(wcs):
    pixin = np.zeros(NAXIS, dtype=np.float64)
    raw = np.zeros(NAXIS, dtype=np.float64)
    rawmin = np.zeros(NAXIS, dtype=np.float64)

    group = [1, 1]
    lng = wcs.wcsprm.lng
    lat = wcs.wcsprm.lat
    if (wcs.wcsprm.pcode == "TNX"):
        tnxflag = 1;
    elif (wcs.wcsprm.pcode, "TAN" or wcs.wcsprm.pcode == "TPV") and (wcs.projp[1+lng*100] or wcs.projp[1+lat*100]):
        tnxflag = 0
    else:
        return

    lngstep = wcs.naxisn[lng]/(WCS_NGRIDPOINTS-1.0)
    latstep = wcs.naxisn[lat]/(WCS_NGRIDPOINTS-1.0)
    outpos = np.zeros(2*WCS_NGRIDPOINTS2, dtype=np.float64)
    lngpos = np.zeros(WCS_NGRIDPOINTS2, dtype=np.float64)
    latpos = np.zeros(WCS_NGRIDPOINTS2, dtype=np.float64)

    for i in range(wcs.naxis):
        raw[i] = rawmin[i] = 0.5
    
    outpost = outpos
    out_index = 0
    lngpost = lngpos
    lnglat_index = 0
    latpost = latpos
    for j in range(WCS_NGRIDPOINTS):
        raw[lng] = rawmin[lng]
        for i in range(WCS_NGRIDPOINTS):
            if (linrev(raw, wcs.lin, pixin)):
                error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ",wcs.wcsprm.pcode)
            
            lngpost[lnglat_index] = pixin[lng]
            latpost[lnglat_index] = pixin[lat]
            lnglat_index += 1
            if (tnxflag):
                outpost[out_index] = pixin[lng] +raw_to_tnxaxis(wcs.tnx_lngcor,pixin[lng],pixin[lat])
                out_index += 1
                outpost[out_index] = pixin[lat]+raw_to_tnxaxis(wcs.tnx_latcor,pixin[lng],pixin[lat])
                out_index +=1
            else:
                raw_to_pv(wcs.prj,pixin[lng],pixin[lat], outpost, outpost+1)
                outpost += 2
            raw[lng]+=lngstep
        raw[lat]+=latstep

    linrev(rawmin, wcs.lin, pixin)
    pixin[lng] += ARCSEC/DEG
    linfwd(pixin, wcs.lin, raw)
    rawsize = sqrt((raw[lng]-rawmin[lng])*(raw[lng]-rawmin[lng])+(raw[lat]-rawmin[lat])*(raw[lat]-rawmin[lat]))*DEG/ARCSEC;
    if (not rawsize):
        error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ",
                wcs.wcsprm.pcode)
    
    epsilon = WCS_INVACCURACY/rawsize
    poly = None
    maxflag = 1
    for deg in range(1,WCS_INVMAXDEG +1):
        if not maxflag:
            break
        if (deg>1):
            poly_end(poly)
        
        poly = poly_init(group, 2, deg, 1)
        poly_fit(poly, outpos, lngpos, None, WCS_NGRIDPOINTS2, None, 0.0)
        maxflag = 0
        outpost = outpos
        lngpost = lngpos
        for i in range(WCS_NGRIDPOINTS2):
            if (abs(poly_func(poly, outpost)-lngpost[i])>epsilon):
                maxflag = 1
                break
            
            outpost+=2
        
    wcs.prj.inv_x = wcs.inv_x = poly

    maxflagflag = maxflag

    linrev(rawmin, wcs.lin, pixin)
    pixin[lat] += ARCSEC/DEG
    linfwd(pixin, wcs.lin, raw)
    rawsize = sqrt((raw[lng]-rawmin[lng])*(raw[lng]-rawmin[lng])+(raw[lat]-rawmin[lat])*(raw[lat]-rawmin[lat]))*DEG/ARCSEC
    if (not rawsize):
        error(EXIT_FAILURE, "*Error*: incorrect linear conversion in ", wcs.wcsprm.pcode)
    
    epsilon = WCS_INVACCURACY/rawsize
    maxflag = 1
    for deg in range(1, WCS_INVMAXDEG+1):
        if not maxflag:
            break
        if (deg>1):
            poly_end(poly)
        
        poly = poly_init(group, 2, deg, 1)
        poly_fit(poly, outpos, latpos, None, WCS_NGRIDPOINTS2, None, 0.0)
        maxflag = 0
        outpost = outpos
        latpost = latpos
        for i in range(WCS_NGRIDPOINTS2):
            if (abs(poly_func(poly, outpost)-latpost[i])>epsilon):
                maxflag = 1
                break
            outpost = outpost[2:]

    wcs.prj.inv_y = wcs.inv_y = poly

    maxflagflag |= maxflag

    if (maxflagflag):
        warning("Significant inaccuracy likely to occur in projection","")
    return


def range_wcs(wcs):
    pixin = np.zeros(NAXIS, dtype=np.float64)
    raw = np.zeros(NAXIS, dtype=np.float64)
    rawmin = np.zeros(NAXIS, dtype=np.float64)
    world = np.zeros(NAXIS, dtype=np.float64)
    world2 = np.zeros(NAXIS, dtype=np.float64)
    linecount = np.zeros(NAXIS, dtype=np.int32)
    
    naxis = wcs.naxis
    npoints = 1
    worldmin = wcs.wcsmin
    worldmax = wcs.wcsmax
    lng = wcs.lng
    lat = wcs.lat
    for i in range(naxis):
        raw[i] = (wcs.naxisn[i]+1.0)/2.0
    
    if (raw_to_wcs(wcs, raw, world)):
        for j in range(100):
            for i in range(naxis):
                raw[i] += wcs.naxisn[i]/100.0*(0.5-random())
            
            if (not raw_to_wcs(wcs, raw, world)):
                break

    if (lng!=lat):
        lc = world[lng]
    else:
        lc = 0.0
        lng = -1
    
    scale = wcs.wcsscale
    for i in range(naxis):
        if ((i==lng or i==lat) and lng!=lat):
            wcs.pixscale = scale[i] = sqrt(wcs_scale(wcs, raw))
        else :
            raw[i] += 1.0
            raw_to_wcs(wcs, raw, world2)
            scale[i] = fabs(world2[i] - world[i])
            raw[i] -= 1.0
            if (lng==lat):
                wcs.pixscale = scale[i]
    
        wcs.wcsscalepos[i] = world[i]
    
    for i in range(naxis):
        raw[i] = rawmin[i] = 0.5
        step[i] = wcs.naxisn[i]/(WCS_NRANGEPOINTS-1.0)
        npoints *= WCS_NRANGEPOINTS
        worldmin[i] = 1e31
        worldmax[i] = -worldmin[i]
        linecount[i] = 0

    radmax = 0.0
    worldc = wcs.wcsscalepos

    for j in range(npoints):
        raw_to_wcs(wcs, raw, world)
        rad=wcs_dist(wcs, world, worldc)
        if (rad > radmax):
            radmax = rad
        
        for i in range(naxis):
            if (i==lng):
                world[i] -= lc
                if (world[i]>180.0):
                    world[i] -= 360.0
                elif (world[i] <= -180.0):
                    world[i] += 360.0
                
            if (world[i]<worldmin[i]) :
                worldmin[i] = world[i]
            
            if (world[i]>worldmax[i]) :
                worldmax[i] = world[i]

        for i in range (naxis):
            raw[i] += step[i]
            linecount[i]+=1
            if (linecount[i]<WCS_NRANGEPOINTS):
                break
            else:
                linecount[i] = 0
                raw[i] = rawmin[i]
            
    wcs.wcsmaxradius = radmax

    if (lng!=lat):
        worldmin[lng] = fmod_0_p360(worldmin[lng]+lc)
        worldmax[lng] = fmod_0_p360(worldmax[lng]+lc)
        if (worldmax[lat]<-90.0):
            worldmax[lat] = -90.0
        
        if (worldmax[lat]>90.0):
            worldmax[lat] = 90.0

    return


def frame_wcs(wcsin, wcsout):
    rawin = np.zeros(NAXIS, dtype=np.float64)
    rawout = np.zeros(NAXIS, dtype=np.float64)
    world = np.zeros(NAXIS, dtype=np.float64)
    linecount = np.zeros(NAXIS, dtype=np.int32)
    
    naxis = wcsin.naxis
    npoints = 1
    min = wcsin.outmin
    max = wcsin.outmax
    for i in range(naxis):
        rawin[i] = 0.5
        npoints *= WCS_NRANGEPOINTS
        min[i] = 1<<30
        max[i] = -min[i]
        linecount[i] = 0
    
    swapflag = (((wcsin.lng != wcsout.lng) or (wcsin.lat != wcsout.lat)) and (wcsin.lng != wcsin.lat) and (wcsout.lng != wcsout.lat))

    for j in range(npoints):
        if (not raw_to_wcs(wcsin, rawin, world)):
            if (swapflag):
                worldc = world[wcsout.lat]
                world[wcsout.lat] = world[wcsin.lat]
                world[wcsin.lat] = worldc
            
            if (not wcs_to_raw(wcsout, world, rawout)):
                for i in range(naxis):
                    out=int(rawout[i]+0.499)
                    if (out<min[i]):
                        min[i] = out
                    
                    if (out>max[i]):
                        max[i] = out

        for i in range(naxis):
            rawin[i] = 0.5 + 0.5*wcsin.naxisn[i]*(1-cos(PI*(linecount[i]+1.0)/(WCS_NRANGEPOINTS-1)))
            linecount[i]+=1
            if (linecount[i]<WCS_NRANGEPOINTS):
                break
            else:
                linecount[i] = 0
                rawin[i] =  0.5
          
    for i in range(naxis):
        if (min[i]>-2147483647):
            min[i] -= 2
        
        if (max[i]>2147483647):
            max[i] += 2
        
    overlapflag = 1
    for i in range(naxis):
        if (min[i]>wcsout.naxisn[i] or max[i]<0):
            overlapflag = 0
            break

    return overlapflag


def reaxe_wcs(wcs, lng, lat):
    olng = wcs.lng
    olat = wcs.lat
    if (lng<0 or lat<0 or olng<0 or olat<0) :
        return RETURN_ERROR

    ilng = wcs.naxisn[olng]
    ilat = wcs.naxisn[olat]
    wcs.naxisn[lng] = ilng
    wcs.naxisn[lat] = ilat
    strcpy(strlng, wcs.ctype[olng])
    strcpy(strlat, wcs.ctype[olat])
    strcpy(wcs.ctype[lng], strlng)
    strcpy(wcs.ctype[lat], strlat)
    dlng = wcs.crval[olng]
    dlat = wcs.crval[olat]
    wcs.crval[lng] = dlng
    wcs.crval[lat] = dlat
    naxis = wcs.naxis
    dlng =  wcs.cd[olng+olng*naxis]
    dlng2 = wcs.cd[olng+olat*naxis]
    dlat =  wcs.cd[olat+olat*naxis]
    dlat2 = wcs.cd[olat+olng*naxis]
    wcs.cd[lng+lng*naxis] = dlng2
    wcs.cd[lng+lat*naxis] = dlng
    wcs.cd[lat+lat*naxis] = dlat2
    wcs.cd[lat+lng*naxis] = dlat
    for l in range(100):
        dlng = wcs.projp[l+olng*100]
        dlat = wcs.projp[l+olat*100]
        wcs.projp[l+lng*100] = dlng
        wcs.projp[l+lat*100] = dlat

    wcsset(wcs.naxis,wcs.ctype, wcs.wcsprm)
    init_wcs(wcs)
    range_wcs(wcs)
    return RETURN_OK


def celsys_to_eq(wcs, wcspos):

    mat = wcs.celsysmat
    lng = wcs.wcsprm.lng
    a2 = wcspos[lng]*DEG - mat[1]
    lat = wcs.wcsprm.lat
    d2 = wcspos[lat]*DEG
    sd2 = sin(d2)
    cd2cp = cos(d2)*mat[2]
    sd = sd2*mat[3]-cd2cp*cos(a2)
    y = cd2cp*sin(a2)
    x = sd2 - sd*mat[3]
    wcspos[lng] = fmod((atan2(y,x) + mat[0])/DEG+360.0, 360.0)
    wcspos[lat] = asin(sd)/DEG

    return RETURN_OK


def eq_to_celsys(wcs, wcspos):
    mat = wcs.celsysmat
    lng = wcs.wcsprm.lng
    a = wcspos[lng]*DEG - mat[0]
    lat = wcs.wcsprm.lat
    d = wcspos[lat]*DEG
    sd = sin(d)
    cdcp = cos(d)*mat[2]
    sd2 = sd*mat[3]+cdcp*cos(a)
    y = cdcp*sin(a)
    x = sd2*mat[3]-sd
    wcspos[lng] = ((atan2(y,x) + mat[1])/DEG+360.0) % 360.0
    wcspos[lat] = asin(sd2)/DEG

    return RETURN_OK



def raw_to_wcs(wcs, pixpos, wcspos):
    imgcrd = np.zeros(NAXIS, dtype=np.float64)

    if (wcsrev(wcs.ctype, wcs.wcsprm, pixpos, wcs.lin,imgcrd, wcs.prj, phi, theta, wcs.crval, wcs.cel, wcspos)):
        for i in range(wcs.naxis):
            wcspos[i] = WCS_NOCOORD
        
        return RETURN_ERROR
    
    if (wcs.celsysconvflag):
        celsys_to_eq(wcs, wcspos)

    return RETURN_OK


def wcs_to_raw(wcs, wcspos, pixpos):
    imgcrd = np.zeros(NAXIS, dtype=np.float64)

    if (wcs.celsysconvflag):
        eq_to_celsys(wcs, wcspos)

    if (wcsfwd(wcs.ctype,wcs.wcsprm,wcspos, wcs.crval, wcs.cel,phi,theta,wcs.prj, imgcrd,wcs.lin,pixpos)):
        for i in range(wcs.naxis):
            pixpos[i] = WCS_NOCOORD
        
        return RETURN_ERROR
    
    return RETURN_OK


def red_to_raw(wcs, redpos, pixpos):
    wcsprm = wcs.wcsprm
    if (wcsprm and wcsprm.flag != WCSSET):
        if (wcsset(wcs.naxis, wcs.ctype, wcsprm)):
            return RETURN_ERROR
    if (wcsprm and wcsprm.flag != 999 and wcsprm.cubeface != -1):
        if wcs.prj.r0 == 0.0:
            offset = 90.0
        else:
            offset = wcs.prj.r0*PI/2.0

        if (redpos[wcs.lat] < -0.5*offset):
            redpos[wcs.lat] += offset
            redpos[wcsprm.cubeface] = 5.0
        elif (redpos[wcs.lat] > 0.5*offset):
            redpos[wcs.lat] -= offset
            redpos[wcsprm.cubeface] = 0.0
        elif (redpos[wcs.lng] > 2.5*offset):
            redpos[wcs.lng] -= 3.0*offset
            redpos[wcsprm.cubeface] = 4.0
        elif (redpos[wcs.lng] > 1.5*offset):
            redpos[wcs.lng] -= 2.0*offset
            redpos[wcsprm.cubeface] = 3.0
        elif (redpos[wcs.lng] > 0.5*offset):
            redpos[wcs.lng] -= offset
            redpos[wcsprm.cubeface] = 2.0
        else:
            redpos[wcsprm.cubeface] = 1.0
        
    if (linfwd(redpos, wcs.lin, pixpos)):
        return RETURN_ERROR

    return RETURN_OK


def raw_to_red(wcs, pixpos, redpos):
    wcsprm = wcs.wcsprm
    if (wcsprm and wcsprm.flag != WCSSET):
        if (wcsset(wcs.naxis, wcs.ctype, wcsprm)):
            return RETURN_ERROR
    
    if (linrev(pixpos, wcs.lin, redpos)) :
        return RETURN_ERROR

    if (wcsprm and wcsprm.flag != 999 and wcsprm.cubeface != -1):
        face = int(redpos[wcsprm.cubeface] + 0.5)
        if (abs(redpos[wcsprm.cubeface]-face) > 1e-10):
            return RETURN_ERROR
        if wcs.prj.r0 == 0.0:
            offset = 90.0
        else:
            offset = wcs.prj.r0*PI/2.0

        if face == 0:
            redpos[wcs.lat] += offset
            break
        elif face == 1:
            break
        elif face == 2:
            redpos[wcs.lng] += offset
            break
        elif face == 3:
            redpos[wcs.lng] += offset*2
            break
        elif face == 4:
            redpos[wcs.lng] += offset*3
            break
        elif face == 5:
            redpos[wcs.lat] -= offset
            break
        else:
            return RETURN_ERROR
    return RETURN_OK


def wcs_dist(wcs, wcspos1, wcspos2):
    lng = wcs.lng
    lat = wcs.lat
    if (lat!=lng):
        d = sin(wcspos1[lat]*DEG)*sin(wcspos2[lat]*DEG)+ cos(wcspos1[lat]*DEG)*cos(wcspos2[lat]*DEG)*cos((wcspos1[lng]-wcspos2[lng])*DEG)
        if d>-1.0:
            if d<1.0:
                return acos(d)/DEG
            else:
                return 0.0
        else:
            return 180.0
    
    d = 0.0
    for i in range(wcs.naxis):
        dp = wcspos1[i] - wcspos2[i]
        d += dp*dp
    
    return sqrt(d)


def wcs_scale(wcs, pixpos):
    wcspos = np.zeros(NAXIS, dtype=np.float64)
    wcspos1 = np.zeros(NAXIS, dtype=np.float64)
    wcspos2 = np.zeros(NAXIS, dtype=np.float64)
    pixpos2 = np.zeros(NAXIS, dtype=np.float64)

    if (raw_to_wcs(wcs, pixpos, wcspos)):
        return 0.0

    lng = wcs.lng
    lat = wcs.lat
    if (lng == lat):
        lng = 0
        lat = 1
    
    pixpos2[lng] = pixpos[lng] + 1.0
    pixpos2[lat] = pixpos[lat]
    if (raw_to_wcs(wcs, pixpos2, wcspos1)):
        return 0.0
    
    pixpos2[lng] -= 1.0
    pixpos2[lat] += 1.0
    if (raw_to_wcs(wcs, pixpos2, wcspos2)):
        return 0.0
    
    dpos1 = wcspos1[lng]-wcspos[lng]
    dpos2 = wcspos2[lng]-wcspos[lng]
    if (wcs.lng!=wcs.lat):
        if (dpos1>180.0):
            dpos1 -= 360.0
        elif (dpos1<-180.0):
            dpos1 += 360.0
        
        if (dpos2>180.0):
            dpos2 -= 360.0
        elif (dpos2<-180.0):
            dpos2 += 360.0
        
        return abs((dpos1*(wcspos2[lat]-wcspos[lat])-(wcspos1[lat]-wcspos[lat])*dpos2)*cos(wcspos[lat]*DEG))
    else :
        return abs((dpos1*(wcspos2[lat]-wcspos[lat])-(wcspos1[lat]-wcspos[lat])*dpos2))

def wcs_jacobian(wcs, pixpos, jacob):
    pixpos0 = np.zeros(NAXIS, dtype=np.float64)
    wcspos0 = np.zeros(NAXIS, dtype=np.float64)
    wcspos = np.zeros(NAXIS, dtype=np.float64)

    lng = wcs.lng
    lat = wcs.lat
    naxis = wcs.naxis
    for i in range(naxis):
        pixpos0[i] = pixpos[i]
    
    if (raw_to_wcs(wcs, pixpos0, wcspos0) == RETURN_ERROR):
        return -1.0
    
    for i in range(naxis):
        pixpos0[i] += 1.0
        if (raw_to_wcs(wcs, pixpos0, wcspos) == RETURN_ERROR):
            return -1.0
        
        pixpos0[i] -= 1.0
        for j in range(naxis):
            dpos = wcspos[j]-wcspos0[j]
            if (lng!=lat and j==lng):
                if (dpos>180.0):
                    dpos -= 360.0
                elif (dpos<-180.0):
                    dpos += 360.0
                dpos *= cos(wcspos0[lat]*DEG)
            
            jacob[j*naxis+i] = dpos

    if (lng==lat):
        lng = 0
        lat = 1

    return abs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat] - jacob[lat+naxis*lng]*jacob[lng+naxis*lat])


def wcs_rawtoraw(wcsin, wcsout, pixposin, pixposout, jacob):
    pixpos0 = np.zeros(NAXIS, dtype=np.float64)
    pixpos2 = np.zeros(NAXIS, dtype=np.float64)
    wcspos = np.zeros(NAXIS, dtype=np.float64)

    naxis = wcsin.naxis
    for i in range(naxis):
        pixpos0[i] = pixposin[i]

    if (raw_to_wcs(wcsin, pixposin, wcspos) == RETURN_ERROR):
        return -1.0
    
    if (wcs_to_raw(wcsout, wcspos, pixposout) == RETURN_ERROR):
        return -1.0
    
    if (not jacob):
        return 0.0

    lng = wcsin.lng
    lat = wcsin.lat
    for i in range(naxis):
        pixpos0[i] += 1.0
        if (raw_to_wcs(wcsin, pixpos0, wcspos) == RETURN_ERROR):
            return -1.0
        
        if (wcs_to_raw(wcsout, wcspos, pixpos2) == RETURN_ERROR):
            return -1.0
        
        pixpos0[i] -= 1.0
        for j in range(naxis):
            jacob[j*naxis+i] = pixpos2[j] - pixposout[j]

    if (lng==lat):
        lng = 0
        lat = 1
    
    return abs(jacob[lng+naxis*lng]*jacob[lat+naxis*lat] - jacob[lat+naxis*lng]*jacob[lng+naxis*lat])


def wcs_chirality(wcs):
    lng = wcs.lng
    lat = wcs.lat
    naxis = wcs.naxis
    if (lng==lat and naxis>=2):
        lng = 0
        lat = 1
    
    a = wcs.cd[lng*naxis+lng]*wcs.cd[lat*naxis+lat] - wcs.cd[lng*naxis+lat]*wcs.cd[lat*naxis+lng]
    if a>TINY:
        return 1
    if a<-TINY:
        return -1
    return 0


def precess_wcs(wcs, yearin, yearout):
    crval = np.zeros(NAXIS, dtype=np.float64)
    a = np.zeros(NAXIS**2, dtype=np.float64)
    b = np.zeros(NAXIS**2, dtype=np.float64)

    lng = wcs.lng
    lat = wcs.lat
    if (lat==lng or yearin==yearout):
        return
    naxis = wcs.naxis
    precess(yearin, wcs.crval[lng], wcs.crval[lat], yearout, crval[lng], crval[lat])

    dalpha = (crval[lng] - wcs.crval[lng])*DEG
    if (dalpha!=0.0 and (crval[lat] - wcs.crval[lat])*DEG != 0.0):
        angle = 180.0 - (atan2(sin(dalpha), cos(crval[lat]*DEG)*tan(wcs.crval[lat]*DEG) - sin(crval[lat]*DEG)*cos(dalpha))  + atan2(sin(dalpha), cos(wcs.crval[lat]*DEG)*tan(crval[lat]*DEG) - sin(wcs.crval[lat]*DEG)*cos(dalpha)))/DEG
    else:
        angle = 0.0

    c = wcs.cd

    cas = cos(angle*DEG)
    sas = sin(-angle*DEG)
    for i in range(naxis):
        b[i+i*naxis] = 1.0
    b[lng+lng*naxis] = cas
    b[lat+lng*naxis] = -sas
    b[lng+lat*naxis] = sas
    b[lat+lat*naxis] = cas
    at = a
    at_index = 0
    for j in range (naxis):
        for i in range(naxis):
            val = 0.0
            for k in range(naxis):
                val += c[k+j*naxis]*b[i+k*naxis]
            
            at[at_index] = val
            at_index+=1

    for i in range(naxis*naxis):
        c[i] = at[i]

    wcs.crval[lng] = crval[lng]
    wcs.crval[lat] = crval[lat]
    wcs.equinox = yearout

    init_wcs(wcs)
    range_wcs(wcs)
    invert_wcs(wcs)

    return


def precess(yearin, alphain, deltain, yearout, alphaout, deltaout):
    alphain *= DEG
    deltain *= DEG

    t1 = (yearin - 2000.0)/1000.0
    t2 = (yearout - yearin)/1000.0
    t1t1 = t1*t1
    t2t2 = t2*t2
    t2t2t2 = (t2t2)*t2
    theta = (97171.735e-06 - 413.691e-06*t1 - 1.052e-06 * t1t1) * t2 + (-206.846e-06 - 1.052e-06*t1) * t2t2 - 202.812e-06 * t2t2t2
    dzeta = (111808.609e-06 + 677.071e-06*t1 - 0.674e-06 * t1t1) * t2 + (146.356e-06 - 1.673e-06*t1) * t2t2 + 87.257e-06 * t2t2t2
    z = (111808.609e-06 +677.071e-06*t1 - 0.674e-06 * t1t1) * t2 + (530.716e-06 + 0.320e-06*t1) * t2t2 + 88.251e-06 * t2t2t2
    cdin=cos(deltain)
    cddsadz = cdin * sin(alphain+dzeta);
    sdin=sin(deltain)
    st=sin(theta)
    ct=cos(theta)
    caindz=cos(alphain+dzeta)
    cddcadz = -sdin*st +cdin*(ct)*(caindz)
    sdd = sdin*ct + cdin*st*caindz
    deltaout = asin(sdd)
    cdd = cos(deltaout)
    adz = asin(cddsadz/cdd)
    if (cddcadz<0.0):
        adz = PI - adz

    if (adz<0.0):
        adz += 2.0*PI

    adz += z
    alphaout = adz/DEG
    deltaout /= DEG

    return


def b2j(yearobs, alphain, deltain, alphaout, deltaout):
    a = [-1.62557e-6, -0.31919e-6, -0.13843e-6]
    ap = [1.245e-3, -1.580e-3, -0.659e-3]
    m = [[ 0.9999256782,     -0.0111820611,     -0.0048579477, 0.00000242395018, -0.00000002710663, -0.00000001177656],
        [0.0111820610,      0.9999374784,     -0.0000271765, 0.00000002710663,  0.00000242397878, -0.00000000006587],
        [0.0048579479,     -0.0000271474,      0.9999881997, 0.00000001177656, -0.00000000006582,  0.00000242410173],
        [-0.000551,        -0.238565,           0.435739, 0.99994704,         -0.01118251,        -0.00485767],
        [0.238514,        -0.002662,          -0.008541, 0.01118251,          0.99995883,        -0.00002718],
        [-0.435623,         0.012254,           0.002117, 0.00485767,      -0.00002714,         1.00000956]]
    a1 = np.zeros(3, dtype=np.float64)
    r = np.zeros(3, dtype=np.float64)
    ro = np.zeros(3, dtype=np.float64)
    r2 = np.zeros(3, dtype=np.float64)
    r2 = np.zeros(3, dtype=np.float64)
    v1 = np.zeros(3, dtype=np.float64)
    v = np.zeros(3, dtype=np.float64)
    t1 = (yearobs - 1950.0)/100.0

    alphain *= PI/180.0
    deltain *= PI/180.0
    cai = cos(alphain)
    sai = sin(alphain)
    cdi = cos(deltain)
    sdi = sin(deltain)
    ro[0] = cdi*cai
    ro[1] = cdi*sai
    ro[2] = sdi
    dotp = 0.0
    for i in range(3):
        a1[i] = a[i]+ap[i]*ARCSEC*t1
        dotp += a1[i]*ro[i]
    
    for i in range(3):
        r1[i] = ro[i] - a1[i] + dotp*ro[i]
        r[i] = v[i] = v1[i] = 0.0
    
    for j in range(6):
        for i in range(6):
            if (j<3):
                if i<3:
                    r[j] += m[j][i]*(r1[i])
                else:
                    r[j] += m[j][i]*(v1[i-3])
                    
            else:
                if i<3:
                    v[j-3] += m[j][i]*(r1[i])
                else:
                    v[j-3] += m[j][i]*(v1[i-3])
                
    rmod = 0.0
    for i in ragen(3):
        r2[i] = r[i]+v[i]*ARCSEC*(t1-0.5)
        rmod += r2[i]*r2[i]
    
    rmod = sqrt(rmod)
    delta = asin(r2[2]/rmod)
    alpha = acos(r2[0]/cos(delta)/rmod)
    if (r2[1]<0):
        alpha = 2*PI - alpha
    
    alphaout = alpha*180.0/PI
    deltaout = delta*180.0/PI

    return


def j2b(yearobs, alphain, deltain, alphaout, deltaout):
    a = [-1.62557e-6, -0.31919e-6, -0.13843e-6]
    ap = [1.245e-3, -1.580e-3, -0.659e-3]
    m[6][6] = [[0.9999256794678425,    0.01118148281196562,   0.004859003848996022,-2.423898417033081e-06,-2.710547600126671e-08,-1.177738063266745e-08],
               [-0.01118148272969232,   0.9999374849247641,   -2.717708936468247e-05,2.710547578707874e-08,-2.423927042585208e-06, 6.588254898401055e-11],
               [-0.00485900399622881,  -2.715579322970546e-05, 0.999988194643078, 1.177738102358923e-08, 6.582788892816657e-11,-2.424049920613325e-06],
               [-0.0005508458576414713, 0.2384844384742432,   -0.4356144527773499,0.9999043171308133,    0.01118145410120206,   0.004858518651645554],
               [-0.2385354433560954,   -0.002664266996872802,  0.01225282765749546,-0.01118145417187502,   0.9999161290795875,   -2.717034576263522e-05],
               [0.4357269351676567,   -0.008536768476441086,  0.002113420799663768,-0.004858518477064975, -2.715994547222661e-05, 0.9999668385070383]]
    a1 = np.zeros(3, dtype=np.float64)
    r = np.zeros(3, dtype=np.float64)
    ro = np.zeros(3, dtype=np.float64)
    r2 = np.zeros(3, dtype=np.float64)
    r2 = np.zeros(3, dtype=np.float64)
    v1 = np.zeros(3, dtype=np.float64)
    v = np.zeros(3, dtype=np.float64)
    
    t1 = ((yearobs - 2000.0) + (MJD2000 - MJD1950)/365.25)*JU2TROP/100.0
    alphain *= DEG
    deltain *= DEG
    cai = cos(alphain)
    sai = sin(alphain)
    cdi = cos(deltain)
    sdi = sin(deltain)
    r[0] = cdi*cai
    r[1] = cdi*sai
    r[2] = sdi
    for i in range(3):
        v[i] = r2[i] = v1[i] = 0.0
    
    for j in range(6):
        for i in range(6):
            if (j<3):
                if i<3:
                    r2[j] += m[j][i]*(r[i])
                else:
                    r2[j] += m[j][i]*(v[i-3])
                    
            else :
                if i<3:
                    v1[j-3] += m[j][i]*(r[i])
                else:
                    v1[j-3] += m[j][i]*(v[i-3])
                
    for i in range(3):
        r1[i] = r2[i]+v1[i]*ARCSEC*t1
    
    dotp = 0.0
    for i in range(3):
        a1[i] = a[i]+ap[i]*ARCSEC*t1
        dotp += a1[i]*(r1[i]+a1[i])
    
    dotp = 2.0/(sqrt(1+4.0*dotp)+1.0)
    rmod = 0.0
    for i in range(3):
        ro[i] = dotp*(r1[i]+a1[i])
        rmod += ro[i]*ro[i]
    
    rmod = sqrt(rmod)
    delta = asin(ro[2]/rmod)
    alpha = acos(ro[0]/cos(delta)/rmod)
    if (ro[1]<0):
        alpha = 2.0*PI - alpha
    
    alphaout = alpha/DEG
    deltaout = delta/DEG

    return


def degtosexal(alpha, str):

    if (alpha>=0.0 and alpha <360.0):
        hh = int(alpha/15.0)
        mm = int(60.0*(alpha/15.0 - hh))
        ss = 60.0*(60.0*(alpha/15.0 - hh) - mm)
    else :
        hh = mm = ss = 0.0
    
    str = "%02d:%02d:%05.2f" % (hh, mm, ss)

    return str


def degtosexde(delta, str):
    if delta<0.0:
        sign = '-'
    else:
        sign = '+'
    delta = abs(delta)
    if (delta>=-90.0 and delta <=90.0):
        dd = int(delta)
        dm = int(60.0*(delta - dd))
        ds = 60.0*abs(60.0*(delta - dd) - dm)
    else :
        dd = dm = ds = 0.0
    
    str = "%c%02d:%02d:%04.1f" % ( sign, dd, dm, ds)
    return str


def sextodegal(hms):
    val = float(strtok_r(hms, ": \t", ptr))*15.0 
    val += float(strtok_r(None, ": \t", ptr))/4.0
    val += float(strtok_r(None, ": \t", ptr))/240.0

    return val

def sextodegde(dms):

    str = strtok_r(dms, ": \t", ptr)
    if find(str, '-')>-1:
        sgn =-1.0
    else:
        sgn=1.0
    val = float(dms)
    val += float(strtok_r(None, ": \t", ptr))*sgn/60.0
    val += float(strtok_r(None, ": \t", ptr))*sgn/3600.0

    return val



def fmod_0_p360(angle):
    if angle>0.0:
        return angle%360.0
    return angle%360.0 + 360.0


def fmod_m90_p90(angle):
    if angle>0.0:
        return (angle+90.0)%180.0 -90.0 
    return (angle-90.0)%180.0 + 90.0


def fmod_0_pmod(angle, mod):
    if angle>0.0:
        return angle%mod
    return angle%mod+mod


def fmod_mmod_pmod(angle, mod):
    if angle>0.0:
        return (angle+mod)%(2.0*mod)-mod 
    return (angle-mod)%(2.0*mod)+mod


def fcmp_0_p360(anglep, anglem):
    dval = anglep - anglem

    return int((dval>0.0 and dval<180.0) or dval<-180.0)

