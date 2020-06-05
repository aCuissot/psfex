/*
*				prefs.h
*
* Include file for prefs.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 1997-2015 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	PSFEx is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	PSFEx is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with PSFEx.  If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		16/12/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



#ifndef _PREFS_H_
#define _PREFS_H_

#include "structs.h"
#include "psf.h"

#include	<ctype.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include    <unistd.h>

#if defined(USE_THREADS) \
		&& (defined(__APPLE__) || defined(FREEBSD) || defined(NETBSD))	/* BSD, Apple */
#include	<sys/types.h>
#include	<sys/sysctl.h>
#elif defined(USE_THREADS) && defined(HAVE_MPCTL)		/* HP/UX */
#include	<sys/mpctl.h>
#endif

#ifdef HAVE_MKL
#include MKL_H
#endif


/* preflist */
#include "xml.h"

#ifdef  USE_THREADS
#define	THREADS_PREFMAX	THREADS_NMAX
#else
#define	THREADS_PREFMAX	65535
#endif

/*-------------------------------- initialization ---------------------------*/
int	idummy;
prefstruct	prefs;

pkeystruct key[] =
 {
  {"BADPIXEL_FILTER", P_BOOL, &prefs.badpix_flag},
  {"BADPIXEL_NMAX", P_INT, &prefs.badpix_nmax, 0,100000000},
  {"BASIS_NAME", P_STRING, prefs.basis_name},
  {"BASIS_NUMBER", P_INT, &prefs.basis_number, 0,10000},
  {"BASIS_SCALE", P_FLOAT, &prefs.basis_scale, 0,0, 0.0,1.0e3},
  {"BASIS_TYPE", P_KEY, &prefs.basis_type, 0,0, 0.0,0.0,
   {"NONE", "PIXEL", "GAUSS-LAGUERRE", "FILE", "PIXEL_AUTO", ""}},
  {"CENTER_KEYS", P_STRINGLIST, prefs.center_key, 0,0,0.0,0.0,
    {""}, 2, 2, &prefs.ncenter_key},
  {"CHECKIMAGE_CUBE", P_BOOL, &prefs.check_cubeflag},
  {"CHECKIMAGE_NAME", P_STRINGLIST, prefs.check_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncheck_name},
  {"CHECKIMAGE_TYPE", P_KEYLIST, prefs.check_type, 0,0, 0.0,0.0,
   {"NONE", "BASIS", "CHI", "PROTOTYPES", "RESIDUALS", "RESIDUALS_GRID",
	"SAMPLES", "SAMPLES_GRID",
	"SNAPSHOTS", "SNAPSHOTS_IMRES", "WEIGHTS",
	"MOFFAT", "-MOFFAT", "-SYMMETRICAL", "GREAT10", ""},
   0, MAXCHECK, &prefs.ncheck_type},
  {"CHECKPLOT_ANTIALIAS", P_BOOL, &prefs.cplot_antialiasflag},
  {"CHECKPLOT_DEV", P_KEYLIST, prefs.cplot_device, 0,0, 0.0,0.0,
    {"NULL", "XWIN", "TK", "XTERM", "PLMETA", "PS", "PSC", "XFIG", "LJIIP",
	"LJ_HPGL", "IMP", "PBM", "PNG", "JPEG", "PSTEX", "AQT", "PDF", "SVG",
	""},
    0, MAXCHECK, &prefs.ncplot_device},
  {"CHECKPLOT_NAME", P_STRINGLIST, prefs.cplot_name, 0,0,0.0,0.0,
    {""}, 0, MAXCHECK, &prefs.ncplot_name},
  {"CHECKPLOT_RES", P_INTLIST, prefs.cplot_res, 0,16384,0.0,0.0,
    {""}, 1, 2, &prefs.ncplot_res},
   {"CHECKPLOT_TYPE", P_KEYLIST, prefs.cplot_type, 0,0, 0.0,0.0,
    {"NONE", "FWHM", "ELLIPTICITY", "MOFFAT_RESIDUALS", "ASYMMETRY",
	"COUNTS", "COUNT_FRACTION", "CHI2", "RESIDUALS", "SELECTION_FWHM", ""},
    0, MAXCHECK, &prefs.ncplot_type},
  {"HIDDENMEF_TYPE", P_KEY, &prefs.hidden_mef_type, 0,0, 0.0,0.0,
	{"INDEPENDENT", "COMMON", ""}},
  {"HOMOBASIS_NUMBER", P_INT, &prefs.homobasis_number, 0,10000},
  {"HOMOBASIS_SCALE", P_FLOAT, &prefs.homobasis_scale, 0,0, 0.0,1.0e3},
  {"HOMOBASIS_TYPE", P_KEY, &prefs.homobasis_type, 0,0, 0.0,0.0,
   {"NONE", "GAUSS-LAGUERRE", ""}},
  {"HOMOKERNEL_DIR", P_STRING, prefs.homokernel_dir},
  {"HOMOKERNEL_SUFFIX", P_STRING, prefs.homokernel_suffix},
  {"HOMOPSF_PARAMS", P_FLOATLIST, prefs.homopsf_params, 0,0, 0.0,100.0, {""},
     2,2, &prefs.nhomopsf_params},
  {"MEF_TYPE", P_KEY, &prefs.psf_mef_type, 0,0, 0.0,0.0,
	{"INDEPENDENT", "COMMON", ""}},
  {"NEWBASIS_TYPE", P_KEY, &prefs.newbasis_type, 0,0, 0.0,0.0,
	{"NONE", "PCA_INDEPENDENT", "PCA_COMMON", ""}},
  {"NEWBASIS_NUMBER", P_INT, &prefs.newbasis_number, 0,1000},
  {"NTHREADS", P_INT, &prefs.nthreads, -THREADS_PREFMAX, THREADS_PREFMAX},
  {"OUTCAT_NAME", P_STRING, prefs.outcat_name},
  {"OUTCAT_TYPE", P_KEY, &prefs.outcat_type, 0,0, 0.0,0.0,
   {"NONE", "ASCII_HEAD", "ASCII", "ASCII_VOTABLE", "FITS_LDAC", ""}},
  {"PHOTFLUX_KEY", P_STRING, prefs.photflux_key},
  {"PHOTFLUXERR_KEY", P_STRING, prefs.photfluxerr_key},
  {"PSFVAR_DEGREES", P_INTLIST, prefs.group_deg, 0,32,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ngroup_deg},
  {"PSFVAR_KEYS", P_STRINGLIST, prefs.context_name, 0,0,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_name},
  {"PSFVAR_GROUPS", P_INTLIST, prefs.context_group, 1,MAXCONTEXT,0.0,0.0,
    {""}, 0, MAXCONTEXT, &prefs.ncontext_group},
  {"PSFVAR_NSNAP", P_INT, &prefs.context_nsnap, 1,256},
  {"PSF_ACCURACY", P_FLOAT, &prefs.prof_accuracy, 0,0, 0.0,1.0},
  {"PSF_DGEOCORRECT", P_BOOL, &prefs.dgeo_flag},
  {"PSF_DIR", P_STRING, prefs.psf_dir},
  {"PSF_PIXELSIZE", P_FLOATLIST, prefs.psf_pixsize, 0,0, 0.0,100.0, {""},
     1,2, &prefs.npsf_pixsize},
  {"PSF_RECENTER", P_BOOL, &prefs.recenter_flag},
  {"PSF_SAMPLING", P_FLOAT, &prefs.psf_step, 0,0, 0.0,1.0e3},
  {"PSF_SIZE", P_INTLIST, prefs.psf_size, 1,1024, 0.0,0.0, {""},
     1,2, &prefs.npsf_size},
  {"PSF_SUFFIX", P_STRING, prefs.psf_suffix},
  {"SAMPLE_AUTOSELECT", P_BOOL, &prefs.autoselect_flag},
  {"SAMPLE_FLAGMASK", P_INT, &prefs.flag_mask, 0,0xffff},
  {"SAMPLE_FWHMRANGE", P_FLOATLIST, prefs.fwhmrange, 0,0, 0.0,1e3, {""},
     2,2, &prefs.nfwhmrange},
  {"SAMPLE_IMAFLAGMASK", P_INT, &prefs.imaflag_mask, 0,0xff, 0.0,0.0},
  {"SAMPLE_MAXELLIP", P_FLOAT, &prefs.maxellip, 0,0, 0.0, 1.0},
  {"SAMPLE_MINSN", P_FLOAT, &prefs.minsn, 0,0, 1e-6,1e15},
//  {"SAMPLE_NMAX", P_INT, &prefs.nmax, 0,2147483648},
  {"SAMPLE_VARIABILITY", P_FLOAT, &prefs.maxvar, 0,0, 0.0, BIG},
  {"SAMPLE_WFLAGMASK", P_INT, &prefs.wflag_mask, 0,0xff, 0.0,0.0},
  {"SAMPLEVAR_TYPE", P_KEY, &prefs.var_type, 0,0, 0.0,0.0,
	{"NONE", "SEEING",""}},
  {"STABILITY_TYPE", P_KEY, &prefs.stability_type, 0,0, 0.0,0.0,
	{"EXPOSURE", "SEQUENCE", ""}},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","LOG","FULL",""}},
  {"XML_NAME", P_STRING, prefs.xml_name},
  {"XSL_URL", P_STRING, prefs.xsl_name},
  {"WRITE_XML", P_BOOL, &prefs.xml_flag},
  {""}
 };

char keylist[sizeof(key)/sizeof(pkeystruct)][32];
extern const char	notokstr[];

extern char *default_prefs[];



/*----------------------------- Internal constants --------------------------*/

#define	MAXCHARL	16384		/* max. nb of chars in a string list */
#define	MAXLIST		64		/* max. nb of list members */
#define	MAXLISTSIZE	2000000		/* max size of list */

/* NOTES:
One must have:	MAXLIST >= 1 (preferably >= 16!)
*/
/*------------------------------- preferences -------------------------------*/

//prefstruct
extern prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);
int	findkeys(char *str, char key[][32], int mode);

#endif
