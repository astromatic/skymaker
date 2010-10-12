/*
*				list.h
*
* Include file for list.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2005-2010 IAP/CNRS/UPMC
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
#define _LIST_H_
/*------------------------------ Internal constants -------------------------*/

#define	READLIST_DISPSTEP	1000

/*------------------------------ Type definitions ---------------------------*/

typedef	struct
	{
	int	type;			/* Type */
	double	pos[2];			/* Position of the object center */
	double	mag;			/* Catalog magnitude */
	double	flux;			/* Instrumental flux */
	double	bulge_ratio;		/* Bulge to total flux */
	double	bulge_req;		/* Bulge equiv. radius */
	double	bulge_posang;		/* Bulge position angle */
	double	bulge_ar;		/* Bulge aspect ratio */
	double	disk_scale;		/* Disk scale length */
	double	disk_ar;		/* Disk aspect ratio */
	double	disk_posang;		/* Disk position angle */
	double	z;			/* Galaxy redshift */
	PIXTYPE	*subimage;		/* Sub-frame with the final object */
	int	subsize[2];		/* Sub-frame dimensions */
	int	subpos[2];		/* Sub-frame position in image */
	double	subfactor;		/* Sub-frame amplitude factor */
	double	*maskbuf;		/* Interpolation kernel buffer */
	int	*nmaskbuf;		/* Interpolation kernel size */
	int	*startbuf;		/* Interpolation kernel start */
	int	buf1size;		/* Size of first buffers series */
	double	*buf2;			/* Second (intermediary) buffer */
	int	buf2size;		/* Size of second buffer */
	}	objstruct;

/*-------------------------------- protos -----------------------------------*/
extern void		closeoutlist(simstruct *sim),
			endobj(objstruct *obj),
			openoutlist(simstruct *sim),
			readlist(simstruct *sim),
			writeobj(simstruct *sim, objstruct *obj);
extern int		readobj(simstruct *sim, objstruct *obj, char *str,
				int proc);
#endif
