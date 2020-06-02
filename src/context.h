/*
*				context.h
*
* Include file for context.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	PSFEx
*
*	Copyright:		(C) 2007-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		20/11/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include "structs.h"


#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "wcs/poly.h"

/*--------------------------------- constants -------------------------------*/

#define		CONTEXT_KEEPHIDDEN	0
#define		CONTEXT_REMOVEHIDDEN	1

/*--------------------------- structure definitions -------------------------*/

//contextstruct

/*-------------------------------- protos -----------------------------------*/

contextstruct	*context_init(char **names, int *group, int ndim, int *degree,
			int ngroup, int pcexflag);
void context_end(contextstruct *context);

#endif

