/*
 *				define.h
 *
 * Global definitions.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	This file part of:	PSFEx
 *
 *	Copyright:		(C) 1997-2016 IAP/CNRS/UPMC
 *
 *	Author:			Emmanuel Bertin (IAP)
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
 *	Last modified:		19/10/2016
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _DEFINE_H_
#define _DEFINE_H_

/* Check if we are using a configure script here */
#include    "fits/fitscat_defs.h"
#include    "fits/fitscat.h"
#include    "wcs/wcs.h"
#include    "wcs/lin.h"
#include    "wcs/tnx.h"
#include    "wcs/poly.h"
#include "config.h"
/*------------------------ what, who, when and where ------------------------*/

#undef      BANNER
#define     BANNER          "PSFEx"

#undef      EXECUTABLE
#define     EXECUTABLE      "psfex"

#define		MYVERSION	VERSION

#undef      COPYRIGHT
#define     COPYRIGHT       "2011-2016 IAP/CNRS/UPMC"

#define		DISCLAIMER	BANNER " comes with ABSOLUTELY NO WARRANTY\n" \
		"You may redistribute copies of " BANNER "\n" \
		"under the terms of the GNU General Public License."

#undef      AUTHORS
#define		AUTHORS		"Emmanuel BERTIN <bertin@iap.fr>"

#undef      WEBSITE
#define		WEBSITE		"http://astromatic.net/software/psfex"

#undef      INSTITUTE
#define		INSTITUTE	"IAP http://www.iap.fr"

/*----------------------------- Internal constants --------------------------*/
#undef      OUTPUT
#define		OUTPUT		stderr		/* where all msgs are sent */

#define		BIG		1e+30		/* a huge number */
#define		MAXCHAR		512		/* max. number of characters */
#define		MAXFILE		32768		/* max number of input files */
#define		ALL_EXTENSIONS	(-1)		/* code for all MEF extensions*/
char gstr[MAXCHAR];

/*------------ Set defines according to machine's specificities -------------*/
#if _LARGEFILE_SOURCE
#define	FSEEKO	fseeko
#define	FTELLO	ftello
#else
#define	FSEEKO	fseek
#define	FTELLO	ftell
#endif

/*---------------------------- return messages ------------------------------*/

#define		RETURN_OK		0
#define		RETURN_ERROR		(-1)
#define		RETURN_FATAL_ERROR	(-2)

/*--------------------- in case of missing constants ------------------------*/

#ifndef PI
#define PI		3.1415926535898	/* never met before? */
#endif

#ifndef         SEEK_SET
#define         SEEK_SET        0
#endif
#ifndef         SEEK_CUR
#define         SEEK_CUR        1
#endif

#ifndef EXIT_SUCCESS
#define 	EXIT_SUCCESS	0
#endif
#ifndef EXIT_FAILURE
#define		EXIT_FAILURE	-1
#endif

/*------------------- a few definitions to read FITS parameters ------------*/

#define	FITSTOF(k, def)	((point = fitsnfind(buf, k, n))? \
					 atof(strncpy(st, &point[10], 70)) \
					:(def))
#define	FITSTOI(k, def)	((point = fitsnfind(buf, k, n))? \
					 atoi(strncpy(st, &point[10], 70)) \
					:(def))

#define	FITSTOS(k, str, def) \
                { if (fitsread(buf,k,str,H_STRING,T_STRING)!= RETURN_OK) \
                    strcpy(str, (def)); \
                }

/*------------------------------- Other Macros -----------------------------*/

#define	DEXP(x)	exp(2.30258509299*(x))	/* 10^x */

#undef QFREAD
#define QFREAD(ptr, size, afile, fname) \
		if (fread(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while reading ", fname)

#undef QFWRITE
#define QFWRITE(ptr, size, afile, fname) \
		if (fwrite(ptr, (size_t)(size), (size_t)1, afile)!=1) \
		  error(EXIT_FAILURE, "*Error* while writing ", fname)

#undef QFSEEK
#define	QFSEEK(afile, offset, pos, fname) \
		if (FSEEKO(afile, (offset), pos)) \
		  error(EXIT_FAILURE,"*Error*: file positioning failed in ", \
			fname)

#undef QFTELL
#define	QFTELL(afile, pos, fname) \
		if ((pos=FTELLO(afile))==-1) \
		  error(EXIT_FAILURE,"*Error*: file position unknown in ", \
			fname)

#undef QFREE
#define	QFREE(ptr) \
		{free(ptr); \
		ptr = NULL;}

#undef QCALLOC
#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%ld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#undef QMALLOC
#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%ld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#undef QREALLOC
#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))))\
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%ld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#undef QMEMCPY
#define QMEMCPY(ptrin, ptrout, typ, nel) \
		{if (ptrin) \
                  {if (!(ptrout = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		     { \
		     sprintf(gstr, #ptrout " (" #nel "=%ld elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		     error(EXIT_FAILURE,"Could not allocate memory for ",gstr);\
                     }; \
                   memcpy(ptrout, ptrin, (size_t)(nel)*sizeof(typ)); \
                   }; \
                 }

#define	RINT(x)	(int)(floor(x+0.5))

#define	NPRINTF		if (prefs.verbose_type == NORM) fprintf

#define	NFPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[1M> %s\n\33[1A",x); \
			else if (prefs.verbose_type == FULL) \
				fprintf(w, "%s.\n", x);}

#define	FPRINTF		if (prefs.verbose_type == FULL)	fprintf

#undef QPRINTF
#define QPRINTF	    if (prefs.verbose_type != QUIET)    fprintf

#define QIPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[7m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
				fprintf(w, "%s\n", x);}

#define QBPRINTF(w,x)	{if (prefs.verbose_type == NORM) \
				fprintf(w, "\33[01;31m%s\33[0m\n", x); \
			else if (prefs.verbose_type == LOG) \
			fprintf(w, "%s\n", x);}

/* Some defines from other files that I putted here */
#define		MAXCHECK	16
#define		NAXIS	2		/* Max number of FITS axes */
#define		MAXCONTEXT		8	/* max. # of context keys */

#endif
