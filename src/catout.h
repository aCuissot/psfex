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

#include "structs.h"

#include "context.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "prefs.h"
#include "xml.h"

/*--------------------------------- constants -------------------------------*/

/*--------------------------------- typedefs --------------------------------*/

//cattypeenum
/*--------------------------- structure definitions -------------------------*/
//outsamplestruct outcatstruct
/* Output catalog fields */

outsamplestruct refoutsample;
extern keystruct refoutkey[];
/*-------------------------------- protos -----------------------------------*/

outcatstruct *
init_outcat (char *filename, int ncontext);

void
end_outcat (outcatstruct *outcat),
write_outcat (outcatstruct *outcat, setstruct *set),
write_vo_fields (FILE *file, tabstruct *objtab);

#endif
