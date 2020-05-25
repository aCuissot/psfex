/*
*				catout.h
*
* Include file for catout.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2014,2015 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	along with PSFEx. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		21/09/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _CATOUT_H_
#define _CATOUT_H_

#include "context.h"
#include "sample.h"

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "fits/fitscat.h"
#include "cathead.h"
#include "prefs.h"
#include "xml.h"

/*--------------------------------- constants -------------------------------*
/*--------------------------------- typedefs --------------------------------*/

typedef enum {CAT_NONE, CAT_ASCII_HEAD, CAT_ASCII,
		CAT_ASCII_VOTABLE, CAT_FITS_LDAC} cattypenum;

/*--------------------------- structure definitions -------------------------*/
typedef struct outsample
  {
  int		detindex;		/* Detection index */
  short		extindex;		/* Extension index */
  int		catindex;		/* Catalog index */
  int		badflag;		/* Rejection flags (0 = accepted) */
  double	context[MAXCONTEXT];	/* Context vector */
  int		ncontext;		/* Number of contexts */
  float		norm;			/* Normalisation flux */
  float		fwhm;			/* FWHM (derived from input catalog) */
  float		ellip;			/* Ellipticity */
  float		snr;			/* Signal-to-Noise ratio */
  double	x,y;			/* x,y position estimate in frame */
  float		dx,dy;			/* x,y shift / vignet center */
  float		chi2;			/* Chi2 of the fit */
  float		modresi;		/* Residual index */
  }	outsamplestruct;

typedef struct outcat
  {
  outsamplestruct	outsample;	/* Current output line */
  FILE			*ascfile;	/* Output ASCII file (if needed) */
  tabstruct		*objtab;	/* Output object table */
  keystruct		*objkeys;  	/* List of output catalog keys */
  char			*buf;		/* Line buffer */
  int			ncontext;	/* Number of contexts */
  }	outcatstruct;


  /* Output catalog fields */

  outsamplestruct	refoutsample;
  keystruct	refoutkey[] = {
    {"SOURCE_NUMBER", "Source index",
  	&refoutsample.detindex, H_INT, T_LONG,
  	"%10d", "", "meta.number", ""},
     {"EXTENSION", "Extension index",
  	&refoutsample.extindex, H_INT, T_SHORT,
  	"%4d", "", "meta.number", ""},
     {"CATALOG_NUMBER", "File index",
  	&refoutsample.catindex, H_INT, T_LONG,
  	"%7d", "", "meta.number", ""},
     {"FLAGS_PSF", "PSFEx rejection flags",
  	&refoutsample.badflag, H_INT, T_LONG,
  	"%7d", "", "meta.code", ""},
     {"VECTOR_CONTEXT", "Context vector",
  	&refoutsample.context, H_FLOAT, T_DOUBLE,
  	"%12.6g", "", "obs.param", "",
  	1, &refoutsample.ncontext},
     {"X_IMAGE", "Position along x image axis",
  	&refoutsample.x, H_FLOAT, T_DOUBLE,
  	"%11.4f", "pixel", "pos.cartesian.x", "pix"},
     {"Y_IMAGE", "Position along y image axis",
  	&refoutsample.y, H_FLOAT, T_DOUBLE,
  	"%11.4f", "pixel", "pos.cartesian.y", "pix"},
     {"DELTAX_IMAGE", "Position offset along x image axis",
  	&refoutsample.dx, H_FLOAT, T_FLOAT,
  	"%11.4f", "pixel", "pos.cartesian.x;arith.diff", "pix"},
     {"DELTAY_IMAGE", "Position offset along y image axis",
  	&refoutsample.dy, H_FLOAT, T_FLOAT,
  	"%11.4f", "pixel", "pos.cartesian.y;arith.diff", "pix"},
     {"NORM_PSF", "Source (inverse) normalization factor",
  	&refoutsample.norm, H_FLOAT, T_FLOAT,
  	"%12.6g", "count", "phot.flux;instr.det.psf", "ct"},
     {"FWHM_PSF", "Source (inverse) normalization factor",
  	&refoutsample.fwhm, H_FLOAT, T_FLOAT,
  	"%8.2f", "pixel", "phys.size.diameter;instr.det.psf", "pix",},
     {"ELLIP_PSF", "Source (inverse) normalization factor",
  	&refoutsample.ellip, H_FLOAT, T_FLOAT,
  	"%7.5f", "", "src.ellipticity;instr.det.psf", ""},
     {"SNR_PSF", "Signal-to-Noise ratio",
  	&refoutsample.snr, H_FLOAT, T_FLOAT,
  	"%10.4g", "", "stat.snr;instr.det.psf", ""},
     {"CHI2_PSF", "PSF fitting chi2/d.o.f.",
  	&refoutsample.chi2, H_FLOAT, T_FLOAT,
  	"%12.6g", "", "stat.fit.chi2;instr.det.psf", ""},
     {"RESI_PSF", "PSF fitting normalized residuals",
  	&refoutsample.modresi, H_FLOAT, T_FLOAT,
  	"%12.6g", "", "stat.fit.residual;instr.det.psf", ""},
    {""},
    };
/*-------------------------------- protos -----------------------------------*/

outcatstruct	*init_outcat(char *filename, int ncontext);

void		end_outcat(outcatstruct *outcat),
		write_outcat(outcatstruct *outcat, setstruct *set),
		write_vo_fields(FILE *file, tabstruct *objtab);

#endif
