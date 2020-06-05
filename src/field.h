/*
*				field.h
*
* Include file for field.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		25/06/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _PSFMEF_H_
#define _PSFMEF_H_
#include "structs.h"

#include "fitswcs.h"
#include "psf.h"
#include "sample.h"

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>


#include	"misc.h"

/*----------------------------- Internal constants --------------------------*/
#define	COUNT_LOADED	1		/* Count detections that are loaded */
#define	COUNT_ACCEPTED	2		/* Count detections that are accepted */
//MAXCHECK IS DEFINED ALSO IN check.h BUT TO AVOID CIRCULAR INCLUDES I DEFINE IT HERE TOO FOR THE MOMENT => should set it in define.h

/*------------------------------ Type definitions ---------------------------*/
/*--------------------------- structure definitions -------------------------*/

//fieldstruct

/*---------------------------------- protos --------------------------------*/
extern fieldstruct	*field_init(char *catname);

extern void		field_count(fieldstruct **fields, setstruct *set,
				int counttype),
			field_end(fieldstruct *field),
			field_locate(fieldstruct *field),
			field_psfsave(fieldstruct *field, char *filename),
			field_stats(fieldstruct **fields, setstruct *set);

void			context_apply(contextstruct *context, psfstruct *psf,
				fieldstruct **fields, int ext, int catindex,
				int ncat),
			context_end(contextstruct *context);

#endif

