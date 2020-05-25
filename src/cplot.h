/*
*				cplot.h
*
* Include file for cplot.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2008-2015 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		28/09/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



#ifndef _CPLOT_H_
#define _CPLOT_H_

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "fitswcs.h"
#include "field.h"

#include	<dlfcn.h>
#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	PLPLOT_H
#include	PLPLOTP_H

#include	"define.h"
#include	"fits/fitscat.h"
#include	"prefs.h"
#include	"psf.h"

/*------------------------------- constants ---------------------------------*/

#define	CPLOT_DEFRESX		800	/* Default X resol. for PNG and JPG */
#define	CPLOT_DEFRESY		600	/* Default X resol. for PNG and JPG */
#define	CPLOT_AAFAC		3	/* Anti-aliasing factor */
#define	CPLOT_NPOINTDEF		1024	/* Default number of points to plot */
#define	CPLOT_FGRIDLINES	9	/* Number of grid lines per axis */
#define	CPLOT_NDISTGRID		32	/* # of distort steps in each CCD dim*/
#define	CPLOT_ASTNSUBPLOTS	3	/* Number of subplot/dim/detector*/
#define	CPLOT_NTYPES		128	/* Number of CPLOT types (typedef below)*/
#define	CPLOT_NSHADES		32	/* Number of shading levels */
#define	CPLOT_SNRVSFWHM_NX	255	/* X resolution of SNR vs FWHM plots */
#define	CPLOT_SNRVSFWHM_NY	255	/* Y resolution of SNR vs FWHM plots */

/*---------------------------- return messages ------------------------------*/
/*-------------------------------- macros -----------------------------------*/
// Work-around to emulate the plwid() function replaced in later versions of
// the PLPlot library.
#ifdef HAVE_PLPLOT
 #ifndef __PLPLOT_H__
  #include	PLPLOT_H
 #endif

 #ifndef __PLPLOTP_H__
  #include	PLPLOTP_H
  #endif

 #ifdef plwidth
  #define	CPLOT_PLWID(wid)	plwidth((PLFLT)(wid))
 #else
  #define	CPLOT_PLWID		plwid
 #endif
#endif

/*--------------------------------- typedefs --------------------------------*/
typedef enum {CPLOT_NONE, CPLOT_FWHM, CPLOT_ELLIPTICITY, CPLOT_MOFFATRESI,
		CPLOT_ASYMRESI, CPLOT_COUNTS, CPLOT_COUNTFRAC, CPLOT_CHI2,
		CPLOT_MODRESI, CPLOT_SNRVSFWHM}	cplotenum;

typedef enum {CPLOT_NULL, CPLOT_XWIN, CPLOT_TK, CPLOT_XTERM, CPLOT_PLMETA,
	CPLOT_PS, CPLOT_PSC, CPLOT_XFIG, CPLOT_LJIIP, CPLOT_LJHPGL, CPLOT_IMP,
	CPLOT_PBM, CPLOT_PNG, CPLOT_JPEG, CPLOT_PSTEX, CPLOT_AQT, CPLOT_PDF,
	CPLOT_SVG} cplotdevenum;

typedef struct {cplotdevenum device; char *devname; char *extension;}
		devicestruct;

/*---------------------------------- svgp -----------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int		cplot_asymresi(fieldstruct *field),
			cplot_countfrac(fieldstruct *field),
			cplot_counts(fieldstruct *field),
			cplot_check(cplotenum cplottype),
			cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout),
			cplot_drawloccoordgrid(wcsstruct *wcs, double xmin,
					double xmax, double ymin, double ymax),
			cplot_ellipticity(fieldstruct *field),
			cplot_end(cplotenum cplottype),
			cplot_fwhm(fieldstruct *field),
			cplot_init(char *name, cplotenum cplottype),
			cplot_modchi2(fieldstruct *field),
			cplot_modresi(fieldstruct *field),
			cplot_moffatresi(fieldstruct *field),
			cplot_snrvsfwhm(fieldstruct *field, setstruct *set);
			
char			*cplot_degtosexal(char *str, double alpha,double step),
			*cplot_degtosexde(char *str, double delta,double step);
#endif

