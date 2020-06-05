/*
*				fitswcs.h
*
* Include file for fitswcs.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 1993-2016 IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		08/03/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#define _FITSWCS_H_
#include "structs.h"
#include "define.h"

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>



/*-------------------------------- macros -----------------------------------*/



/*----------------------------- Internal constants --------------------------*/

extern char	celsysname[][2][8];
extern double	celsysorig[][2], celsyspole[][2], celsyssign[];



#define		DEG	(PI/180.0)	/* 1 deg in radians */
#define		ARCMIN	(DEG/60.0)	/* 1 arcsec in radians */
#define		ARCSEC	(DEG/3600.0)	/* 1 arcsec in radians */
#define		MAS	(ARCSEC/1000.0)	/* 1 mas in radians */
#define		YEAR	(365.25*DAY)	/* 1 year in seconds */
#define		DAY	(24.0*HOUR)	/* 1 day in seconds */
#define		HOUR	(60.0*MINUTE)	/* 1 hour in seconds */
#define		MINUTE	60.0		/* 1 minute in seconds */
#define		MJD2000	51544.50000	/* Modified Julian date for J2000.0 */
#define		JD2000	(2400000.5+MJD2000)	/* Julian date for J2000.0 */
#define		MJD1950	33281.92346	/* Modified Julian date for B1950.0 */
#define		JD1950	(2400000.5+MJD1950)	/* Julian date for B1950.0 */
#define		JU2TROP	1.0000214	/* 1 Julian century in tropical units*/
#define		WCS_NOCOORD	1e31	/* Code for non-existing coordinates */

#define		WCS_NGRIDPOINTS	12	/* Number of WCS grid points / axis */
#define		WCS_NGRIDPOINTS2	(WCS_NGRIDPOINTS*WCS_NGRIDPOINTS)
#define		WCS_INVMAXDEG	9	/* Maximum inversion polynom degree */
#define		WCS_INVACCURACY	0.001	/* Maximum inversion error (pixels) */
#define		WCS_NRANGEPOINTS 32	/* Number of WCS range points / axis */

/*-------------------------------- typedefs ---------------------------------*/

//celsysenum

/*------------------------------- structures --------------------------------*/

//wcsstruct

/*------------------------------- functions ---------------------------------*/

extern wcsstruct	*create_wcs(char **ctype, double *crval, double *crpix,
				double *cdelt, int *naxisn, int naxis),
			*copy_wcs(wcsstruct *wcsin),
			*read_wcs(tabstruct *tab);

extern double		fmod_0_p360(double angle),
			fmod_0_pmod(double angle, double mod),
			fmod_m90_p90(double angle),
			fmod_mmod_pmod(double angle, double mod),
			sextodegal(char *hms),
			sextodegde(char *dms),
			wcs_dist(wcsstruct *wcs,
				double *wcspos1, double *wcspos2),
			wcs_jacobian(wcsstruct *wcs, double *pixpos,
				double *jacob),
			wcs_rawtoraw(wcsstruct *wcsin, wcsstruct *wcsout,
				double *pixposin, double *pixposout,
				double *jacob),
			wcs_scale(wcsstruct *wcs, double *pixpos);

extern int		celsys_to_eq(wcsstruct *wcs, double *wcspos),
			eq_to_celsys(wcsstruct *wcs, double *wcspos),
			fcmp_0_p360(double anglep, double anglem),
			frame_wcs(wcsstruct *wcsin, wcsstruct *wcsout),
			raw_to_red(wcsstruct *wcs,
				double *pixpos, double *redpos),
			raw_to_wcs(wcsstruct *wcs,
				double *pixpos, double *wcspos),
			reaxe_wcs(wcsstruct *wcs, int lng, int lat),
			red_to_raw(wcsstruct *wcs,
				double *redpos, double *pixpos),
			wcs_chirality(wcsstruct *wcs),
			wcs_supproj(char *name),
			wcs_to_raw(wcsstruct *wcs,
				double *wcspos, double *pixpos);

extern char		*degtosexal(double alpha, char *str),
			*degtosexde(double delta, char *str);

extern void		b2j(double yearobs, double alphain, double deltain,
				double *alphaout, double *deltaout),
			end_wcs(wcsstruct *wcs),
			init_wcs(wcsstruct *wcs),
			init_wcscelsys(wcsstruct *wcs),
			invert_wcs(wcsstruct *wcs),
			j2b(double yearobs, double alphain, double deltain,
				double *alphaout, double *deltaout),
			precess(double yearin, double alphain, double deltain,
				double yearout,
				double *alphaout, double *deltaout),
			precess_wcs(wcsstruct *wcs, double yearin,
				double yearout),
			range_wcs(wcsstruct *wcs),
			wipe_wcs(tabstruct *tab),
			write_wcs(tabstruct *tab, wcsstruct *wcs);

#endif
