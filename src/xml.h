/*
*				xml.h
*
* Include file for xml.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2005-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _XML_H_
#define _XML_H_
#include "structs.h"

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "field.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "cplot.h"
#include "prefs.h"
#include "psf.h"
/*----------------------------- Internal constants --------------------------*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif

#define         FIND_STRICT     0
#define         FIND_NOSTRICT   1

/*--------------------------------- typedefs --------------------------------*/
//struct pkeystruct
/*------------------------------- functions ---------------------------------*/

extern int	init_xml(int ncat),
		update_xml(fieldstruct *field),
		write_xml(char *filename),
		write_xml_header(FILE *file),
		write_xml_meta(FILE *file, char *error),
		write_xmlconfigparam(FILE *file, char *name, char *unit,
				char *ucd, char *format);
extern void	end_xml(void),
		write_xmlerror(char *filename, char *error);
#endif
