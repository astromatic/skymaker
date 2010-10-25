/*
*				image.h
*
* Include file for image.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SkyMaker is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SkyMaker is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SkyMaker. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		12/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif
#ifndef _LIST_H_
#include "list.h"
#endif

/*----------------------------- Internal constants --------------------------*/

#define	INTERP_TYPE		INTERP_LANCZOS3	/* Default interpolation */
#define	INTERP_MAXKERNELWIDTH	8	/* Max. range of kernel (pixels) */

/*--------------------------------- typedefs --------------------------------*/
typedef enum {INTERP_NEARESTNEIGHBOUR, INTERP_BILINEAR, INTERP_LANCZOS2,
		INTERP_LANCZOS3, INTERP_LANCZOS4}	interpenum;

/*----------------------------- Global variables ----------------------------*/

/*------------------------------- functions ---------------------------------*/
extern int	add_image(PIXTYPE *pix1, int w1, int h1, PIXTYPE *pix2,
			int w2, int h2, int ix, int iy, float amplitude),
		resample_image(PIXTYPE *pix1, int w1, int h1, objstruct *obj,
			double dx, double dy, double step2);

