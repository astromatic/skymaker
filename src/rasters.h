/*
*				rasters.h
*
* Include file for rasters.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		06/05/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SIMUL_H_
#include "simul.h"
#endif

#ifndef _LIST_H_
#include "list.h"
#endif

/*---------------------------- Internal constants ---------------------------*/

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern int	make_raster(simstruct *sim, objstruct *obj);

extern double	raster_raster(PIXTYPE *pix, int width, int height, double *jac,
		double size, double aspect, double posang);
