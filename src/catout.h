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

#ifndef _CONTEXT_H_
#include "context.h"
#endif

#ifndef _SAMPLE_H_
#include "sample.h"
#endif

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "cathead.h"
#include "prefs.h"
#include "context.h"
#include "xml.h"

/*--------------------------------- constants -------------------------------*/
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

/*-------------------------------- protos -----------------------------------*/

outcatstruct	*init_outcat(char *filename, int ncontext);

void		end_outcat(outcatstruct *outcat),
		write_outcat(outcatstruct *outcat, setstruct *set),
		write_vo_fields(FILE *file, tabstruct *objtab);

#endif
