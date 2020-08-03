"""extern time_t        thetime,thetime2;  
extern pkeystruct    key[]
extern char        keylist[][32]
"""

from time import time
import numpy as np
from define import *

field_xml = None
nxml, nxmlmax = 0, 0

def init_xml(ncat):
    if (not ncat) :
        field_xml = None
#    else:
#        QMALLOC(field_xml, fieldstruct *, ncat);
    nxml = 0
    nxmlmax = ncat

    return EXIT_SUCCESS


def end_xml():
    return


def update_xml(field):
    field_xml[nxml] = field
    nxml+=1

    return EXIT_SUCCESS


def write_xml(filename):
    file = open(prefs.xml_name, "w")
    if  not file:
        return RETURN_ERROR

    write_xml_header(file)
    write_xml_meta(file, None)

    file.write("</RESOURCE>\n")
    file.wrifte("</VOTABLE>\n")

    file.close()

    return RETURN_OK


def write_xml_header(file):
    file.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    file.write("<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n" % (prefs.xsl_name))
    file.write("<VOTABLE version=\"1.1\"\n"
            " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
            " xsi:noNamespaceSchemaLocation=\"http://www.ivoa.net/xml/VOTable/v1.1\">\n")
    file.write("<DESCRIPTION>produced by %s</DESCRIPTION>\n" % (BANNER))
    file.write("<!-- VOTable description at "
            "http://www.ivoa.net/Documents/latest/VOT.html -.\n")
    file.write("<RESOURCE ID=\"%s\" name=\"%s\">\n" % (BANNER, BANNER))
    file.write(" <DESCRIPTION>Data related to %s"
            "</DESCRIPTION>\n" % (BANNER))
    file.write(" <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n")
    sysname = "ICRS"

    file.write(" <COOSYS ID=\"J2000\" equinox=\"J2000\""
            " epoch=\"2000.0\" system=\"ICRS\"/>\n")

    return RETURN_OK



def write_xml_meta(file, error):
    cp = np.zeros(CPLOT_NTYPES, dtype=np.int32)
    
    if (error):
        thetime = time()
        today = date.today()
        now = datetime.now()
        prefs.sdate_start = today.strftime("%Y-%m-%d")
        prefs.stime_start = now.strftime("%H:%M:%S")
        prefs.time_diff = thetime2 - thetime

    psuser = pspath = pshost = None
#ifdef HAVE_GETENV
    psuser=getenv("USERNAME")
    if not (psuser):
        psuser = getenv("LOGNAME")
    pspath = getenv("PWD")
    pshost = getenv("HOSTNAME")
#endif

    file.write( " <RESOURCE ID=\"MetaData\" name=\"MetaData\">\n")
    file.write( "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n" % (BANNER))
    file.write( "  <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n")
    file.write( "  <PARAM name=\"Software\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.title;meta.software\" value=\"%s\"/>\n" %
            (BANNER))
    file.write( "  <PARAM name=\"Version\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.version;meta.software\" value=\"%s\"/>\n" % 
            (MYVERSION))
    file.write( "  <PARAM name=\"Soft_URL\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n"
            % (WEBSITE))
    file.write( "  <PARAM name=\"Soft_Auth\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n" %
            ("Emmanuel Bertin"))
    file.write( "  <PARAM name=\"Soft_Ref\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n"
            % ("2006ASPC..351..112B"))
    file.write( "  <PARAM name=\"NThreads\" datatype=\"int\""
            " ucd=\"meta.number;meta.software\" value=\"%d\"/>\n" %
            (prefs.nthreads))
    file.write( "  <PARAM name=\"Date\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n" %
            (prefs.sdate_end))
    file.write( "  <PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n" %
            (prefs.stime_end))
    file.write( "  <PARAM name=\"Duration\" datatype=\"float\""
            " ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n" %
            (prefs.time_diff))

    file.write( "  <PARAM name=\"User\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.curation\" value=\"%s\"/>\n" %
            (psuser))
    file.write( "  <PARAM name=\"Host\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.curation\" value=\"%s\"/>\n" %
            (pshost))
    file.write( "  <PARAM name=\"Path\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta.dataset\" value=\"%s\"/>\n" %
            (pspath))

    if (error):
        file.write( "\n  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                "!!!!!!!!!!!!!!!!!!!! -.\n")
        file.write( "  <!-- !!!!!!!!!!!!!!!!!!!!!! an Error occured"
                " !!!!!!!!!!!!!!!!!!!!! -.\n")
        file.write( "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                "!!!!!!!!!!!!!!!!!!!! -.\n")
        file.write("  <PARAM name=\"Error_Msg\" datatype=\"char\" arraysize=\"*\""
                " ucd=\"meta\" value=\"%s\"/>\n" % (error))
        file.write( "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                "!!!!!!!!!!!!!!!!!!!! -.\n")
        file.write( "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                "!!!!!!!!!!!!!!!!!!!! -.\n\n")

    file.write( "  <TABLE ID=\"PSF_Fields\" name=\"PSF_Fields\">\n")
    file.write( "   <DESCRIPTION>PSF metadata and stats per field gathered by "
            "%s</DESCRIPTION>\n" % (BANNER))
    file.write( "   <!-- NFields may be 0"
            " if an error occurred early in the processing -.\n")
    file.write( "   <PARAM name=\"NFields\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n" % (nxmlmax))
    file.write( "   <!-- CurrField may differ from NFields"
            " if an error occurred -.\n")
    file.write( "   <PARAM name=\"CurrField\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n" % (nxml))
    if prefs.ncontext_group:
        file.write( "   <PARAM name=\"NSnapshots\" datatype=\"int\""
            " arraysize=\"%d\" ucd=\"meta.number;meta.dataset\""
            " value=\"%d" % 
            (prefs.ncontext_group, prefs.context_nsnap))
    else:
        file.write( "   <PARAM name=\"NSnapshots\" datatype=\"int\""
            " arraysize=\"%d\" ucd=\"meta.number;meta.dataset\""
            " value=\"%d" % 
            (1,1))
    for d in range(prefs.ncontext_group):
        file.write( " %d" % (prefs.context_nsnap))
    
    file.write( "\"/>\n")

    file.write( "   <FIELD name=\"Catalog_Name\" datatype=\"char\""
            " arraysize=\"*\" ucd=\"meta.id;meta.table;meta.file\"/>\n")
    file.write( "   <FIELD name=\"Image_Ident\" datatype=\"char\""
            " arraysize=\"*\" ucd=\"meta.id;obs.field\"/>\n")
    file.write( "   <FIELD name=\"NExtensions\" datatype=\"int\""
            " ucd=\"meta.number\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Total\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Min\" datatype=\"int\""
            " ucd=\"meta.number;stat.min;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Mean\" datatype=\"float\""
            " ucd=\"meta.number;stat.mean;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Max\" datatype=\"int\""
            " ucd=\"meta.number;stat.max;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Total\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Min\" datatype=\"int\""
            " ucd=\"meta.number;stat.min;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Mean\" datatype=\"int\""
            " ucd=\"meta.number;stat.mean;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Max\" datatype=\"int\""
            " ucd=\"meta.number;stat.max;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Min\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Max\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Min\" datatype=\"float\""
            " ucd=\"stat.param;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Mean\" datatype=\"float\""
            " ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Max\" datatype=\"float\""
            " ucd=\"stat.param;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Min\" datatype=\"float\""
            " ucd=\"stat.param;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Max\" datatype=\"float\""
            " ucd=\"stat.param;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Min\" datatype=\"float\""
            " ucd=\"phys.area;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Mean\" datatype=\"float\""
            " ucd=\"phys.area;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Max\" datatype=\"float\""
            " ucd=\"phys.area;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.min;instr.pixel\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.mean;instr.pixel\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.max;instr.pixel\"/>\n")


    file.write( "   <DATA><TABLEDATA>\n")
    for n in range(nxml):
        nloaded_min = naccepted_min = 2<<29
        nloaded_max = naccepted_max = nloaded_total = naccepted_total = 0
        minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min = ellipticity_min = ellipticity1_min = ellipticity2_min = beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min = pfellipticity_min = pfellipticity1_min = pfellipticity2_min = pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min = pixscale_wcs_min = BIG
        minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean = ellipticity_mean = ellipticity1_mean = ellipticity2_mean = beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean = pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean = pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean = pixscale_wcs_mean = nloaded_mean = naccepted_mean = 0.0
        minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max = ellipticity_max = ellipticity1_max = ellipticity2_max = beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max = pfellipticity_max = pfellipticity1_max = pfellipticity2_max = pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max = pixscale_wcs_max = -BIG
        neff = 0
        field = field_xml[n]
        for e in range(field.next):
            psf = field.psf[e]
            nloaded_total += psf.samples_loaded
            if (psf.samples_loaded < nloaded_min):
                nloaded_min = psf.samples_loaded
            nloaded_mean += float(psf.samples_loaded)
            if (psf.samples_loaded > nloaded_max):
                nloaded_max = psf.samples_loaded
            naccepted_total += psf.samples_accepted
            if (psf.samples_accepted < naccepted_min):
                naccepted_min = psf.samples_accepted
            naccepted_mean += float(psf.samples_accepted)
            if (psf.samples_accepted > naccepted_max):
                naccepted_max = psf.samples_accepted
            if (not psf.samples_accepted):
                continue
            neff+=1
            if (psf.fwhm < minrad_min):
                minrad_min = psf.fwhm
            minrad_mean += psf.fwhm
            if (psf.fwhm > minrad_max):
                minrad_max = psf.fwhm
            if (psf.pixstep < sampling_min):
                sampling_min = psf.pixstep
            sampling_mean += psf.pixstep
            if (psf.pixstep > sampling_max):
                sampling_max = psf.pixstep
            if (psf.chi2 < chi2_min):
                chi2_min = psf.chi2
            chi2_mean += psf.chi2
            if (psf.chi2 > chi2_max):
                chi2_max = psf.chi2
            if (psf.moffat_fwhm_min < fwhm_min):
                fwhm_min = psf.moffat_fwhm_min
            fwhm_mean += psf.moffat_fwhm
            if (psf.moffat_fwhm_max > fwhm_max):
                fwhm_max = psf.moffat_fwhm_max
            if (psf.moffat_fwhm_wcs_min < fwhm_wcs_min):
                fwhm_wcs_min = psf.moffat_fwhm_wcs_min
            fwhm_wcs_mean += psf.moffat_fwhm_wcs
            if (psf.moffat_fwhm_wcs_max > fwhm_wcs_max):
                fwhm_wcs_max = psf.moffat_fwhm_wcs_max
            if (psf.moffat_ellipticity_min < ellipticity_min):
                ellipticity_min = psf.moffat_ellipticity_min
            ellipticity_mean += psf.moffat_ellipticity;
            if (psf.moffat_ellipticity_max > ellipticity_max):
                ellipticity_max = psf.moffat_ellipticity_max
            if (psf.moffat_ellipticity1_min < ellipticity1_min):
                ellipticity1_min = psf.moffat_ellipticity1_min
            ellipticity1_mean += psf.moffat_ellipticity1;
            if (psf.moffat_ellipticity1_max > ellipticity1_max):
                ellipticity1_max = psf.moffat_ellipticity1_max
            if (psf.moffat_ellipticity2_min < ellipticity2_min):
                ellipticity2_min = psf.moffat_ellipticity2_min
            ellipticity2_mean += psf.moffat_ellipticity2
            if (psf.moffat_ellipticity2_max > ellipticity2_max):
                ellipticity2_max = psf.moffat_ellipticity2_max
            if (psf.moffat_beta_min < beta_min):
                beta_min = psf.moffat_beta_min
            beta_mean += psf.moffat_beta
            if (psf.moffat_beta_max > beta_max):
                beta_max = psf.moffat_beta_max
            if (psf.moffat_residuals_min < residuals_min):
                residuals_min = psf.moffat_residuals_min
            residuals_mean += psf.moffat_residuals
            if (psf.moffat_residuals_max > residuals_max):
                residuals_max = psf.moffat_residuals_max
            if (psf.pfmoffat_fwhm_min < pffwhm_min):
                pffwhm_min = psf.pfmoffat_fwhm_min
            pffwhm_mean += psf.pfmoffat_fwhm
            if (psf.pfmoffat_fwhm_max > pffwhm_max):
                pffwhm_max = psf.pfmoffat_fwhm_max
            if (psf.pfmoffat_fwhm_wcs_min < pffwhm_wcs_min):
                pffwhm_wcs_min = psf.pfmoffat_fwhm_wcs_min
            pffwhm_wcs_mean += psf.pfmoffat_fwhm_wcs
            if (psf.pfmoffat_fwhm_wcs_max > pffwhm_wcs_max):
                pffwhm_wcs_max = psf.pfmoffat_fwhm_wcs_max
            if (psf.pfmoffat_ellipticity_min < pfellipticity_min):
                pfellipticity_min = psf.pfmoffat_ellipticity_min
            pfellipticity_mean += psf.pfmoffat_ellipticity
            if (psf.pfmoffat_ellipticity_max > pfellipticity_max):
                pfellipticity_max = psf.pfmoffat_ellipticity_max
            if (psf.pfmoffat_ellipticity1_min < pfellipticity1_min):
                pfellipticity1_min = psf.pfmoffat_ellipticity1_min
            pfellipticity1_mean += psf.pfmoffat_ellipticity1
            if (psf.pfmoffat_ellipticity1_max > pfellipticity1_max):
                pfellipticity1_max = psf.pfmoffat_ellipticity1_max
            if (psf.pfmoffat_ellipticity2_min < pfellipticity2_min):
                pfellipticity2_min = psf.pfmoffat_ellipticity2_min
            pfellipticity2_mean += psf.pfmoffat_ellipticity2
            if (psf.pfmoffat_ellipticity2_max > pfellipticity2_max):
                pfellipticity2_max = psf.pfmoffat_ellipticity2_max
            if (psf.pfmoffat_beta_min < pfbeta_min):
                pfbeta_min = psf.pfmoffat_beta_min
            pfbeta_mean += psf.pfmoffat_beta
            if (psf.pfmoffat_beta_max > pfbeta_max):
                pfbeta_max = psf.pfmoffat_beta_max
            if (psf.pfmoffat_residuals_min < pfresiduals_min):
                pfresiduals_min = psf.pfmoffat_residuals_min
            pfresiduals_mean += psf.pfmoffat_residuals
            if (psf.pfmoffat_residuals_max > pfresiduals_max):
                pfresiduals_max = psf.pfmoffat_residuals_max
            if (psf.sym_residuals_min < symresiduals_min):
                symresiduals_min = psf.sym_residuals_min
            symresiduals_mean += psf.sym_residuals
            if (psf.sym_residuals_max > symresiduals_max):
                symresiduals_max = psf.sym_residuals_max
            if (psf.noiseqarea_min < noiseqarea_min):
                noiseqarea_min = psf.noiseqarea_min
            noiseqarea_mean += psf.noiseqarea
            if (psf.noiseqarea_max > noiseqarea_max):
                noiseqarea_max = psf.noiseqarea_max
            if (psf.pixscale_wcs_min < pixscale_wcs_min):
                pixscale_wcs_min = psf.pixscale_wcs_min
            pixscale_wcs_mean += psf.pixscale_wcs
            if (psf.pixscale_wcs_max > pixscale_wcs_max):
                pixscale_wcs_max = psf.pixscale_wcs_max

        if (field.next>1) :
            nloaded_mean /= field.next
            naccepted_mean /= field.next
        

        if (neff>1):
            minrad_mean /= neff;
            sampling_mean /= neff
            chi2_mean /= neff
            fwhm_mean /= neff
            fwhm_wcs_mean /= neff
            ellipticity_mean /= neff
            ellipticity1_mean /= neff
            ellipticity2_mean /= neff
            beta_mean /= neff
            residuals_mean /= neff
            pffwhm_mean /= neff
            pffwhm_wcs_mean /= neff
            pfellipticity_mean /= neff
            pfellipticity1_mean /= neff
            pfellipticity2_mean /= neff
            pfbeta_mean /= neff
            pfresiduals_mean /= neff
            symresiduals_mean /= neff
            noiseqarea_mean /= neff
            pixscale_wcs_mean /= neff
        elif (neff==0):
            minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min = ellipticity_min = ellipticity1_min = ellipticity2_min = beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min = pfellipticity_min = pfellipticity1_min = pfellipticity2_min= pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min= pixscale_wcs_min= minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean = ellipticity_mean = ellipticity1_mean = ellipticity2_mean = beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean= pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean= pixscale_wcs_mean= minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max= ellipticity_max  = ellipticity1_max = ellipticity2_max= beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max= pfellipticity_max = pfellipticity1_max = pfellipticity2_max = pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max = pixscale_wcs_max = 0.0

        file.write( "    <TR>\n"
                "     <TD>%s</TD><TD>%s</TD><TD>%d</TD>\n"
                "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
                "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n" %
                (field.rcatname, field.ident, field.next,
                nloaded_total, nloaded_min, nloaded_mean, nloaded_max,
                naccepted_total, naccepted_min, naccepted_mean, naccepted_max,
                minrad_min, minrad_mean, minrad_max,
                sampling_min, sampling_mean, sampling_max,
                chi2_min, chi2_mean, chi2_max,
                fwhm_min, fwhm_mean, fwhm_max,
                fwhm_wcs_min, fwhm_wcs_mean, fwhm_wcs_max,
                ellipticity_min, ellipticity_mean, ellipticity_max,
                ellipticity1_min, ellipticity1_mean, ellipticity1_max,
                ellipticity2_min, ellipticity2_mean, ellipticity2_max,
                beta_min, beta_mean, beta_max,
                residuals_min, residuals_mean, residuals_max,
                pffwhm_min, pffwhm_mean, pffwhm_max,
                pffwhm_wcs_min, pffwhm_wcs_mean, pffwhm_wcs_max,
                pfellipticity_min, pfellipticity_mean, pfellipticity_max,
                pfellipticity1_min, pfellipticity1_mean, pfellipticity1_max,
                pfellipticity2_min, pfellipticity2_mean, pfellipticity2_max,
                pfbeta_min, pfbeta_mean, pfbeta_max,
                pfresiduals_min, pfresiduals_mean, pfresiduals_max,
                symresiduals_min, symresiduals_mean, symresiduals_max,
                noiseqarea_min, noiseqarea_mean, noiseqarea_max,
                pixscale_wcs_min, pixscale_wcs_mean, pixscale_wcs_max))

        file.write( "    </TR>\n")
    file.write( "   </TABLEDATA></DATA>\n")
    file.write( "  </TABLE>\n")
    if nxml:
        next = field_xml[0].next
    else:
        next = 0
    file.write( "  <TABLE ID=\"PSF_Extensions\" name=\"PSF_Extensions\">\n")
    file.write( "   <DESCRIPTION>PSF metadata and stats per extension gathered"
            " by %s</DESCRIPTION>\n"% (BANNER))
    file.write( "   <!-- NExtensions may be 0"
            " if an error occurred early in the processing -.\n")
    file.write( "   <PARAM name=\"NExtensions\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n" % (next))
    file.write( "   <FIELD name=\"Extension\" datatype=\"int\""
            " ucd=\"meta.record\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Total\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Min\" datatype=\"int\""
            " ucd=\"meta.number;stat.min;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Mean\" datatype=\"float\""
            " ucd=\"meta.number;stat.mean;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Loaded_Max\" datatype=\"int\""
            " ucd=\"meta.number;stat.max;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Total\" datatype=\"int\""
            " ucd=\"meta.number;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Min\" datatype=\"int\""
            " ucd=\"meta.number;stat.min;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"NStars_Accepted_Mean\" datatype=\"float\""
            " ucd=\"meta.number;stat.mean;meta.dataset\"/>\n");
    file.write( "   <FIELD name=\"NStars_Accepted_Max\" datatype=\"int\""
            " ucd=\"meta.number;stat.max;meta.dataset\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_FromFluxRadius_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Sampling_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"arith.factor;instr.pixel;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Min\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Chi2_Max\" datatype=\"float\""
            " ucd=\"stat.fit.chi2;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Min\" datatype=\"float\""
            " ucd=\"stat.param;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Mean\" datatype=\"float\""
            " ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_Max\" datatype=\"float\""
            " ucd=\"stat.param;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Min\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Mean\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_Max\" unit=\"pix\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"FWHM_PixelFree_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"phys.size.diameter;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity1_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Min\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Ellipticity2_PixelFree_Max\" datatype=\"float\""
            " ucd=\"src.ellipticity;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Min\" datatype=\"float\""
            " ucd=\"stat.param;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"stat.param;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"MoffatBeta_PixelFree_Max\" datatype=\"float\""
            " ucd=\"stat.param;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Residuals_PixelFree_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Min\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Mean\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Asymmetry_Max\" datatype=\"float\""
            " ucd=\"stat.fit.residual;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Min\" datatype=\"float\""
            " ucd=\"phys.area;stat.min;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Mean\" datatype=\"float\""
            " ucd=\"phys.area;stat.mean;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"Area_Noise_Max\" datatype=\"float\""
            " ucd=\"phys.area;stat.max;instr.det.psf\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Min\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.min;instr.pixel\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Mean\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.mean;instr.pixel\"/>\n")
    file.write( "   <FIELD name=\"PixelScale_WCS_Max\" unit=\"arcsec\""
            " datatype=\"float\""
            " ucd=\"instr.scale;stat.max;instr.pixel\"/>\n")

    file.write( "   <DATA><TABLEDATA>\n")
    for e in range(next):
        nloaded_min = naccepted_min = 2<<29
        nloaded_max = naccepted_max = nloaded_total = naccepted_total = 0
        minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min = ellipticity_min = ellipticity1_min = ellipticity2_min = beta_min = residuals_min = pffwhm_min  = pfellipticity_min = pfellipticity1_min = pfellipticity2_min= pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min = pixscale_wcs_min = BIG;
        minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean = ellipticity_mean = ellipticity1_mean = ellipticity2_mean = beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean = pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean = pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean = pixscale_wcs_mean = nloaded_mean = naccepted_mean = 0.0
        minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max= ellipticity_max = ellipticity1_max = ellipticity2_max = beta_max = residuals_max = pffwhm_max = pfellipticity_max = pfellipticity1_max = pfellipticity2_max = pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max = pixscale_wcs_max = -BIG
        neff = 0
        for n in range(nxml):
            field = field_xml[n]
            psf = field.psf[e]
            nloaded_total += psf.samples_loaded
            if (psf.samples_loaded < nloaded_min):
                nloaded_min = psf.samples_loaded
            nloaded_mean += psf.samples_loaded
            if (psf.samples_loaded > nloaded_max):
                nloaded_max = psf.samples_loaded
            naccepted_total += psf.samples_accepted
            if (psf.samples_accepted < naccepted_min):
                naccepted_min = psf.samples_accepted
            naccepted_mean += psf.samples_accepted
            if (psf.samples_accepted > naccepted_max):
                naccepted_max = psf.samples_accepted

            if (not psf.samples_accepted):
                continue
            neff+=1
            if (psf.fwhm < minrad_min):
                minrad_min = psf.fwhm
            minrad_mean += psf.fwhm
            if (psf.fwhm > minrad_max):
                minrad_max = psf.fwhm
            if (psf.pixstep < sampling_min): 
                sampling_min = psf.pixstep
            sampling_mean += psf.pixstep
            if (psf.pixstep > sampling_max):
                sampling_max = psf.pixstep
            if (psf.chi2 < chi2_min):
                chi2_min = psf.chi2
            chi2_mean += psf.chi2
            if (psf.chi2 > chi2_max):
                chi2_max = psf.chi2
            if (psf.moffat_fwhm_min < fwhm_min):
                fwhm_min = psf.moffat_fwhm_min
            fwhm_mean += psf.moffat_fwhm
            if (psf.moffat_fwhm_max > fwhm_max):
                fwhm_max = psf.moffat_fwhm_max
            if (psf.moffat_fwhm_wcs_min < fwhm_wcs_min):
                fwhm_wcs_min = psf.moffat_fwhm_wcs_min
            fwhm_wcs_mean += psf.moffat_fwhm_wcs
            if (psf.moffat_fwhm_wcs_max > fwhm_wcs_max):
                fwhm_wcs_max = psf.moffat_fwhm_wcs_max
            if (psf.moffat_ellipticity_min < ellipticity_min):
                ellipticity_min = psf.moffat_ellipticity_min
            ellipticity_mean += psf.moffat_ellipticity
            if (psf.moffat_ellipticity_max > ellipticity_max):
                ellipticity_max = psf.moffat_ellipticity_max
            if (psf.moffat_ellipticity1_min < ellipticity1_min):
                ellipticity1_min = psf.moffat_ellipticity1_min
            ellipticity1_mean += psf.moffat_ellipticity1
            if (psf.moffat_ellipticity1_max > ellipticity1_max):
                ellipticity1_max = psf.moffat_ellipticity1_max
            if (psf.moffat_ellipticity2_min < ellipticity2_min):
                ellipticity2_min = psf.moffat_ellipticity2_min
            ellipticity2_mean += psf.moffat_ellipticity2
            if (psf.moffat_ellipticity2_max > ellipticity2_max):
                ellipticity2_max = psf.moffat_ellipticity2_max
            if (psf.moffat_beta_min < beta_min):
                beta_min = psf.moffat_beta_min
            beta_mean += psf.moffat_beta
            if (psf.moffat_beta_max > beta_max):
                beta_max = psf.moffat_beta_max
            if (psf.moffat_residuals_min < residuals_min):
                residuals_min = psf.moffat_residuals_min
            residuals_mean += psf.moffat_residuals
            if (psf.moffat_residuals_max > residuals_max):
                residuals_max = psf.moffat_residuals_max
            if (psf.pfmoffat_fwhm_min < pffwhm_min):
                pffwhm_min = psf.pfmoffat_fwhm_min
            pffwhm_mean += psf.pfmoffat_fwhm
            if (psf.pfmoffat_fwhm_max > pffwhm_max):
                pffwhm_max = psf.pfmoffat_fwhm_max
            if (psf.pfmoffat_fwhm_wcs_min < pffwhm_wcs_min):
                pffwhm_wcs_min = psf.pfmoffat_fwhm_wcs_min
            pffwhm_wcs_mean += psf.pfmoffat_fwhm_wcs
            if (psf.pfmoffat_fwhm_wcs_max > pffwhm_wcs_max):
                pffwhm_wcs_max = psf.pfmoffat_fwhm_wcs_max
            if (psf.pfmoffat_ellipticity_min < pfellipticity_min):
                pfellipticity_min = psf.pfmoffat_ellipticity_min
            pfellipticity_mean += psf.pfmoffat_ellipticity
            if (psf.pfmoffat_ellipticity_max > pfellipticity_max):
                pfellipticity_max = psf.pfmoffat_ellipticity_max
            if (psf.pfmoffat_ellipticity1_min < pfellipticity1_min):
                pfellipticity1_min = psf.pfmoffat_ellipticity1_min
            pfellipticity1_mean += psf.pfmoffat_ellipticity1
            if (psf.pfmoffat_ellipticity1_max > pfellipticity1_max):
                pfellipticity1_max = psf.pfmoffat_ellipticity1_max
            if (psf.pfmoffat_ellipticity2_min < pfellipticity2_min):
                pfellipticity2_min = psf.pfmoffat_ellipticity2_min
            pfellipticity2_mean += psf.pfmoffat_ellipticity2
            if (psf.pfmoffat_ellipticity2_max > pfellipticity2_max):
                pfellipticity2_max = psf.pfmoffat_ellipticity2_max
            if (psf.pfmoffat_beta_min < pfbeta_min):
                pfbeta_min = psf.pfmoffat_beta_min
            pfbeta_mean += psf.pfmoffat_beta
            if (psf.pfmoffat_beta_max > pfbeta_max):
                pfbeta_max = psf.pfmoffat_beta_max
            if (psf.pfmoffat_residuals_min < pfresiduals_min):
                pfresiduals_min = psf.pfmoffat_residuals_min
            pfresiduals_mean += psf.pfmoffat_residuals
            if (psf.pfmoffat_residuals_max > pfresiduals_max) :
                pfresiduals_max = psf.pfmoffat_residuals_max
            if (psf.sym_residuals_min < symresiduals_min):
                symresiduals_min = psf.sym_residuals_min
            symresiduals_mean += psf.sym_residuals
            if (psf.sym_residuals_max > symresiduals_max):
                symresiduals_max = psf.sym_residuals_max
            if (psf.noiseqarea_min < noiseqarea_min):
                noiseqarea_min = psf.noiseqarea_min
            noiseqarea_mean += psf.noiseqarea
            if (psf.noiseqarea_max > noiseqarea_max):
                noiseqarea_max = psf.noiseqarea_max
            if (psf.pixscale_wcs_min < pixscale_wcs_min):
                pixscale_wcs_min = psf.pixscale_wcs_min
            pixscale_wcs_mean += psf.pixscale_wcs
            if (psf.pixscale_wcs_max > pixscale_wcs_max):
                pixscale_wcs_max = psf.pixscale_wcs_max

        if (nxml>1):
            nloaded_mean /= nxml
            naccepted_mean /= nxml

        if (neff>1):
            minrad_mean /= neff
            sampling_mean /= neff
            chi2_mean /= neff
            fwhm_mean /= neff
            fwhm_wcs_mean /= neff
            ellipticity_mean /= neff
            ellipticity1_mean /= neff
            ellipticity2_mean /= neff
            beta_mean /= neff
            residuals_mean /= neff
            pffwhm_mean /= neff
            pffwhm_wcs_mean /= neff
            pfellipticity_mean /= neff
            pfellipticity1_mean /= neff
            pfellipticity2_mean /= neff
            pfbeta_mean /= neff
            pfresiduals_mean /= neff
            symresiduals_mean /= neff
            noiseqarea_mean /= neff
            pixscale_wcs_mean /= neff
        elif (neff==0) :
            minrad_min = sampling_min = chi2_min = fwhm_min = fwhm_wcs_min= ellipticity_min = ellipticity1_min = ellipticity2_min= beta_min = residuals_min = pffwhm_min = pffwhm_wcs_min = pfellipticity_min = pfellipticity1_min = pfellipticity2_min= pfbeta_min = pfresiduals_min = symresiduals_min = noiseqarea_min= pixscale_wcs_min = minrad_mean = sampling_mean = chi2_mean = fwhm_mean = fwhm_wcs_mean = ellipticity_mean = ellipticity1_mean = ellipticity2_mean= beta_mean = residuals_mean = pffwhm_mean = pffwhm_wcs_mean= pfellipticity_mean = pfellipticity1_mean = pfellipticity2_mean = pfbeta_mean = pfresiduals_mean = symresiduals_mean = noiseqarea_mean= pixscale_wcs_mean= minrad_max = sampling_max = chi2_max = fwhm_max = fwhm_wcs_max= ellipticity_max  = ellipticity1_max = ellipticity2_max = beta_max = residuals_max = pffwhm_max = pffwhm_wcs_max= pfellipticity_max = pfellipticity1_max = pfellipticity2_max= pfbeta_max = pfresiduals_max = symresiduals_max = noiseqarea_max= pixscale_wcs_max = 0.0
        
        file.write( "    <TR>\n"
                "     <TD>%d</TD>\n"
                "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
                "     <TD>%d</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n"
                "    </TR>\n" % 
                (e+1,
                nloaded_total, nloaded_min, nloaded_mean, nloaded_max,
                naccepted_total, naccepted_min, naccepted_mean, naccepted_max,
                minrad_min, minrad_mean, minrad_max,
                sampling_min, sampling_mean, sampling_max,
                chi2_min, chi2_mean, chi2_max,
                fwhm_min, fwhm_mean, fwhm_max,
                fwhm_wcs_min, fwhm_wcs_mean, fwhm_wcs_max,
                ellipticity_min, ellipticity_mean, ellipticity_max,
                ellipticity1_min, ellipticity1_mean, ellipticity1_max,
                ellipticity2_min, ellipticity2_mean, ellipticity2_max,
                beta_min, beta_mean, beta_max,
                residuals_min, residuals_mean, residuals_max,
                pffwhm_min, pffwhm_mean, pffwhm_max,
                pffwhm_wcs_min, pffwhm_wcs_mean, pffwhm_wcs_max,
                pfellipticity_min, pfellipticity_mean, pfellipticity_max,
                pfellipticity1_min, pfellipticity1_mean, pfellipticity1_max,
                pfellipticity2_min, pfellipticity2_mean, pfellipticity2_max,
                pfbeta_min, pfbeta_mean, pfbeta_max,
                pfresiduals_min, pfresiduals_mean, pfresiduals_max,
                symresiduals_min, symresiduals_mean, symresiduals_max,
                noiseqarea_min, noiseqarea_mean, noiseqarea_max,
                pixscale_wcs_min, pixscale_wcs_mean, pixscale_wcs_max))
    

    file.write( "   </TABLEDATA></DATA>\n")
    file.write( "  </TABLE>\n")

    file.write( "  <TABLE ID=\"Warnings\" name=\"Warnings\">\n")
    file.write(
            "   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n" % 
            (BANNER, WARNING_NMAX))
    file.write( "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta;time.event.end\"/>\n")
    file.write( "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta;time.event.end\"/>\n")
    file.write( "   <FIELD name=\"Msg\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"meta\"/>\n")
    file.write( "   <DATA><TABLEDATA>\n")
    str = warning_history()
    while(str[0]):
        file.write( "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n" %
                (str, str+11, str+22))
        str = warning_history()
    
    file.write( "   </TABLEDATA></DATA>\n")
    file.write( "  </TABLE>\n")

    file.write( "  <RESOURCE ID=\"Config\" name=\"Config\">\n")
    file.write( "   <DESCRIPTION>%s configuration</DESCRIPTION>\n" % (BANNER))
    file.write(
            "   <PARAM name=\"Command_Line\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"obs.param\" value=\"%s" % 
            (prefs.command_line[0]))
    for n in range(prefs.ncommand_line):
        file.write( " %s" % (prefs.command_line[n]))
    
    file.write( "\"/>\n")
    file.write(
            "   <PARAM name=\"Prefs_Name\" datatype=\"char\" arraysize=\"*\""
            " ucd=\"obs.param;meta.file\" value=\"%s\"/>\n" %
            (prefs.prefs_name))

    if (not error):
        write_xmlconfigparam(file, "Basis_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "Basis_Number", "", "meta.number","%d")
        write_xmlconfigparam(file, "Basis_Name", "", "meta.id;meta.file","%s")
        write_xmlconfigparam(file, "Basis_Scale", "", "arith.factor","%.6g")
        write_xmlconfigparam(file, "NewBasis_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "NewBasis_Number", "", "meta.number","%d")
        write_xmlconfigparam(file, "PSF_Sampling", "pix",
                "arith.factor;instr.pixel;inst.det.psf","%.6g")
        write_xmlconfigparam(file, "PSF_Accuracy", "",
                "obs.param;phot.flux.sb;arith.ratio;instr.det.psf","%.6g")
        write_xmlconfigparam(file, "PSF_Size", "pix",
                "meta.number;instr.pixel;instr.det.psf","%d")
        write_xmlconfigparam(file, "Center_Keys", "",
                "meta.id;src;instr.det.psf", "%s");
        write_xmlconfigparam(file, "PSF_Recenter", "", "meta.code","%c")
        write_xmlconfigparam(file, "PhotFlux_Key", "",
                "meta.id;src;instr.det.psf", "%s")
        write_xmlconfigparam(file, "PhotFluxErr_Key", "",
                "meta.id;src;instr.det.psf", "%s")
        write_xmlconfigparam(file, "MEF_Type", "", "meta.code","%s")

        write_xmlconfigparam(file, "PSFVar_Keys", "",
                "meta.id;src;instr.det.psf", "%s")
        write_xmlconfigparam(file, "PSFVar_Groups", "",
                "meta.id;stat.fit.param;instr.det.psf", "%d")
        write_xmlconfigparam(file, "PSFVar_Degrees", "",
                "stat.fit.param;instr.det.psf", "%d")
        write_xmlconfigparam(file, "PSFVar_NSnap", "",
                "stat.fit.param;instr.det.psf", "%d")
        write_xmlconfigparam(file, "HiddenMEF_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "Stability_Type", "", "meta.code","%s")

        write_xmlconfigparam(file, "Sample_AutoSelect", "", "meta.code","%c")
        write_xmlconfigparam(file, "SampleVar_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "Sample_FWHMRange", "pix",
                "phys.size.diameter;instr.det.psf","%.6g")
        write_xmlconfigparam(file, "Sample_Variability", "",
                "instr.det.psf;arith.ratio","%.6g")
        write_xmlconfigparam(file, "Sample_MinSN", "", "stat.snr;stat.min","%.6g")
        write_xmlconfigparam(file, "Sample_MaxEllip", "",
                "src.ellipticity;stat.max","%.6g")
        write_xmlconfigparam(file, "Sample_FlagMask", "", "meta.code","%d")
        write_xmlconfigparam(file, "BadPixel_Filter", "", "meta.code","%c")
        write_xmlconfigparam(file, "BadPixel_NMax", "",
                "meta.number;instr.pixel;stat.max","%d")

        write_xmlconfigparam(file, "HomoBasis_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "HomoBasis_Number", "", "meta.number","%d")
        write_xmlconfigparam(file, "HomoBasis_Scale", "", "arith.factor","%.6g")
        write_xmlconfigparam(file, "HomoPSF_Params", "", "stat.param","%.6g")
        write_xmlconfigparam(file, "HomoKernel_Dir", "", "meta.id;meta.file","%s")
        write_xmlconfigparam(file, "HomoKernel_Suffix", "",
                "meta.id;meta.file","%s")

        write_xmlconfigparam(file, "CheckPlot_Dev", "", "meta.code", "%s")
        write_xmlconfigparam(file, "CheckPlot_Res", "", "meta.number;meta", "%d")
        write_xmlconfigparam(file, "CheckPlot_AntiAlias", "", "meta.code", "%c")
        write_xmlconfigparam(file, "CheckPlot_Type", "", "meta.code", "%s")
        write_xmlconfigparam(file, "CheckPlot_Name", "", "meta.id;meta.file", "%s")

        write_xmlconfigparam(file, "CheckImage_Type", "", "meta.code", "%s")
        write_xmlconfigparam(file, "CheckImage_Name", "",
                "meta.id;meta.file;meta.fits", "%s")
        write_xmlconfigparam(file, "CheckImage_Cube", "", "meta.code","%c")

        write_xmlconfigparam(file, "PSF_Dir", "", "meta.id;meta.file","%s")
        write_xmlconfigparam(file, "PSF_Suffix", "", "meta.id;meta.file","%s")
        write_xmlconfigparam(file, "Verbose_Type", "", "meta.code","%s")
        write_xmlconfigparam(file, "Write_XML", "", "meta.code","%s")
        write_xmlconfigparam(file, "NThreads", "",
                "meta.number;meta.software", "%d")

    file.write( "  </RESOURCE>\n")
    file.write( " </RESOURCE>\n")

    return RETURN_OK


def write_xmlerror(filename, error):

    pipe_flag = 0;
    if (filename == "STDOUT"):
        file = stdout
        pipe_flag = 1
    else:
        file = open(filename, "w")
        if (not file): 
            return

    write_xml_header(file)
    write_xml_meta(file, error)

    file.write( "</RESOURCE>\n")
    file.write( "</VOTABLE>\n")

    if (not pipe_flag):
        file.close()
    
    return


def write_xmlconfigparam(file, name, unit, ucd, format):
    i=0
    
    while key[i].name[0] and name != key[i].name:
        i+=1
    
    if (not key[i].name[0]):
        return RETURN_ERROR
    
    if (unit):
        uunit = " unit=\"%s\"" % (unit)
    else:
        uunit = '\0'

    
    if key[i].type == P_FLOAT:
        value = format % (key[i].ptr[0])
        file.write( "   <PARAM name=\"%s\"%s datatype=\"double\""
                " ucd=\"%s\" value=\"%s\"/>\n" %
                (name, uunit, ucd, value))
        break
    elif key[i].type == P_FLOATLIST:
        n = key[i].nlistptr[0]
        if (n):
            value = format % (key[i].ptr[0])
            file.write( "   <PARAM name=\"%s\"%s datatype=\"double\""
                    " arraysize=\"%d\" ucd=\"%s\" value=\"%s" %
                    (name, uunit, n, ucd, value))
            for j in range(1, n):
                value = format % (key[i].ptr[j])
                file.write( " %s" % (value))
            file.write( "\"/>\n")
        else :
            file.write( "   <PARAM name=\"%s\"%s datatype=\"double\""
                    " ucd=\"%s\" value=\"\"/>\n" % 
                    (name, uunit, ucd))
        break
    elif key[i].type == P_INT:
        value = format % (int(key[i].ptr[0]))
        file.write( "   <PARAM name=\"%s\"%s datatype=\"int\""
                " ucd=\"%s\" value=\"%s\"/>\n" %
                (name, uunit, ucd, value))
        break
    elif key[i].type == P_INTLIST:
        n = key[i].nlistptr[0]
        if (n):
            value = format % (int(key[i].ptr)[0])
            file.write( "   <PARAM name=\"%s\"%s datatype=\"int\""
                    " arraysize=\"%d\" ucd=\"%s\" value=\"%s" %
                    (name, uunit, n, ucd, value))
            for j in range(n):
                value = format % (int(key[i].ptr[j]))
                file.write( " %s" % (value))
            file.write( "\"/>\n")
        else:
            file.write( "   <PARAM name=\"%s\"%s datatype=\"double\""
                    " ucd=\"%s\" value=\"\"/>\n" % 
                    (name, uunit, ucd))
        break
    elif key[i].type == P_BOOL:
        if key[i].ptr[0]:
            value = "%c" % ('T')
        else:
            value = "%c" % ('F')

        file.write( "   <PARAM name=\"%s\" datatype=\"boolean\""
                " ucd=\"%s\" value=\"%s\"/>\n" %
                (name, ucd, value))
        break
    elif key[i].type == P_BOOLLIST:
        n = key[i].nlistptr[0]
        if (n):
            if key[i].ptr[0]:
                value = "%c" % ('T')
            else:
                value = "%c" % ('F')
            file.write( "   <PARAM name=\"%s\" datatype=\"boolean\""
                    " arraysize=\"%d\" ucd=\"%s\" value=\"%s" % 
                    (name, n, ucd, value))
            for j in range(1, n):
                if key[i].ptr[j]:
                        value = "%c" % ('T')
                else:
                    value = "%c" % ('F')
                file.write( " %s" % (value))
            
            file.write( "\"/>\n")
        else:
            file.write( "   <PARAM name=\"%s\" datatype=\"boolean\""
                    " ucd=\"%s\" value=\"\"/>\n" %
                    (name, ucd))
        break
    elif key[i].type == P_STRING:
        value = "%s" % (key[i].ptr[0])
        if value[0]:
            file.write( "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
                " ucd=\"%s\" value=\"%s\"/>\n" % 
                (name, ucd, value))
        else:
            file.write( "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
                " ucd=\"%s\" value=\"%s\"/>\n" % 
                (name, ucd, " "))
        
        break
    elif key[i].type == P_STRINGLIST:
        n = key[i].nlistptr[0]
        if (n) :
            value = "%s" % (key[i].ptr[0])
            if value:
                file.write( "   <PARAM name=\"%s\" datatype=\"char\""
                        " arraysize=\"*\" ucd=\"%s\" value=\"%s" % 
                        (name, ucd,  value))
            else:
                file.write( "   <PARAM name=\"%s\" datatype=\"char\""
                    " arraysize=\"*\" ucd=\"%s\" value=\"%s" %
                    (name, ucd, " "))
            for j in range(1, n):
                value = "%s" % (key[i].ptr[j])
                if value[0]:
                    file.write( ",%s" % (value))
                else:
                    file.write( ",%s" % (" "))
                    

            file.write( "\"/>\n")
        else:
            file.write( "   <PARAM name=\"%s\" datatype=\"char\""
                    " arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n" % 
                    (name, ucd))
        break
    elif key[i].type == P_KEY:
        value = "%s" % (key[i].keylist[int(key[i].ptr)])
        file.write( "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
                " ucd=\"%s\" value=\"%s\"/>\n" % 
                (name, ucd, value))
        break
    elif key[i].type == P_KEYLIST:
        n = key[i].nlistptr[0]
        if (n) :
            value ="%s" % (key[i].keylist[(int(key[i].ptr)[0])])
            file.write( "   <PARAM name=\"%s\" datatype=\"char\""
                    " arraysize=\"*\" ucd=\"%s\" value=\"%s"%
                    (name, ucd, value))
            for j in range(1, n):
                value = "%s" % (key[i].keylist[(int(key[i].ptr)[j])]);
                file.write( ",%s" % (value))
            
            file.write( "\"/>\n")
        else:
            file.write( "   <PARAM name=\"%s\" datatype=\"char\""
                    " arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n" % 
                    (name, ucd))
        
        break
    else:
        error(EXIT_FAILURE, "*Internal Error*: Type Unknown",
                " in write_xmlconfigparam()")

    return RETURN_OK