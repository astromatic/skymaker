/*
*				corr.h
*
* Include file for corr.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2017 IAP/CNRS/UPMC
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
*	Last modified:		04/05/2017
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SIMUL_H_
#include "simul.h"
#endif

/*---------------------------- Internal constants --------------------------*/
#define	CORRFUNC_STEP		0.01	/* Corr. function integration step */

/*--------------------------- structure definitions ------------------------*/
typedef struct {
  corrfuncenum 	corrfunc_type;		/* Correlation interp. type */
  float		radius;			/* Correlation kernel size */
  PIXTYPE 	(*func)(float x);	/* Pointer to kernel function */
}	corrstruct;

/*---------------------------- External constants --------------------------*/
extern const corrstruct	coors[];

/*---------------------------------- protos --------------------------------*/
extern void	corr_conv(simstruct *sim, PIXTYPE **image),
		corr_convline(simstruct *sim, PIXTYPE *imagein,
			PIXTYPE *imageout, int y),
		corr_generate(simstruct *sim, corrfuncenum interp_type,
			float scale);

