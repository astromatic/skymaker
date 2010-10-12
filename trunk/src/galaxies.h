/*
*				galaxies.h
*
* Include file for galaxies.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2010 IAP/CNRS/UPMC
*
*	Author:			Emmanuel Bertin (IAP)
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

#ifndef _SIMUL_H_
#include "simul.h"
#endif

#ifndef _LIST_H_
#include "list.h"
#endif

/*---------------------------- Internal constants ---------------------------*/
#define	VDKCUTRAD	5.0   /* van der Kruit disk truncation radius in r_h */

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	make_galaxy(simstruct *sim, objstruct *obj);
extern double	add_devauc(PIXTYPE *pix, double xcenter, double ycenter,
			int width, int height,
			double flux, double req, double aspect, double posang),
		add_expo(PIXTYPE *pix, double xcenter, double ycenter,
			int width, int height,
			double flux, double scale,
			double aspect, double posang);
	
extern PIXTYPE	trunc_prof(PIXTYPE *pix, double xcenter, double ycenter,
			int width, int height);
