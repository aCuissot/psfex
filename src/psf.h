/*
*				psf.h
*
* Include file for psf.c
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
*	Last modified:		21/09/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/



#ifndef _PSF_H_
#define _PSF_H_
#include "structs.h"

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include "wcs/poly.h"

#include "sample.h"

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"prefs.h"
#include	"misc.h"
#include	"vignet.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
//#define MATSTORAGE_PACKED 1
#endif


/*----------------------------- Internal constants --------------------------*/

#define	PSF_NODIAG	0	/* Don't do diagnostics */
#define	PSF_DIAG	1	/* Do diagnostics */
#define	PSF_FREEDFACTOR	1.1	/* Margin against overfitting (10%) */
#define	PSF_NMASKDIM	3	/* Number of dimensions for PSF data */
#define	PSF_MAXSHIFT	3.0	/* Max shift from initial guess (pixels)*/
#define	PSF_MINSHIFT	1e-4	/* Min shift from previous guess (pixels)*/
#define PSF_NITER	40	/* Maximum number of iterations in fit */
#define	PSF_NSNAPMAX	16	/* Maximum number of PSF snapshots/dimension */
#define	GAUSS_LAG_OSAMP	3	/* Gauss-Laguerre oversampling factor */
#define	PSF_AUTO_FWHM	3.0	/* FWHM theshold for PIXEL-AUTO mode */
#define	PSF_NORTHOSTEP	16	/* Number of PSF orthonor. snapshots/dimension*/

/*----------------------------- Type definitions --------------------------*/

//basistypeenum

/*--------------------------- structure definitions -------------------------*/

// moffatstruct psfstruct

/*---------------------------------- protos --------------------------------*/
extern void	psf_build(psfstruct *psf, double *pos),
		psf_clip(psfstruct *psf),
		psf_end(psfstruct *psf),
		psf_make(psfstruct *psf, setstruct *set, double prof_accuracy),
		psf_makebasis(psfstruct *psf, setstruct *set,
			basistypenum basis_type,  int nvec),
		psf_makeresi(psfstruct *psf, setstruct *set, int centflag,
			double prof_accuracy),
		psf_makemask(psfstruct *psf, setstruct *set, double chithresh),
		psf_orthopoly(psfstruct *psf, setstruct *set),
		psf_save(psfstruct *psf,  char *filename, int ext, int next);

extern int	psf_pshapelet(float **shape, int w, int h, int nmax,
			double beta),
		psf_readbasis(psfstruct *psf, char *filename, int ext),
		psf_refine(psfstruct *psf, setstruct *set);

extern double	psf_chi2(psfstruct *psf, setstruct *set),
		psf_clean(psfstruct *psf, setstruct *set, double prof_accuracy);

extern psfstruct	*psf_copy(psfstruct *psf),
			*psf_inherit(contextstruct *context, psfstruct *psf),
			*psf_init(contextstruct *context, int *size,
				float psfstep, float *pixsize, int nsample),
			*psf_load(char *filename);

#endif

