/*
*				list.h
*
* Include file for list.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2018 IAP/CNRS/UPMC
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
*	Last modified:		23/02/2018
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

typedef enum	{LISTCOORD_PIXEL, LISTCOORD_WORLD}	listcoordenum;

typedef	struct
	{
	int	type;			/* Source type */
	double	pos[2];			/* Position of object center in pixels*/
	double	wcspos[2];		/* Position of object center in WCS */
	float	mag;			/* Catalog magnitude */
	float	flux;			/* Instrumental flux */
	float	bulge_ratio;		/* Bulge to total flux */
	float	bulge_req;		/* Bulge equiv. radius */
	float	bulge_posang;		/* Bulge position angle */
	float	bulge_ar;		/* Bulge aspect ratio */
	float	bulge_sersicn;		/* Bulge sersic index */
	float	disk_scale;		/* Disk scale length */
	float	disk_ar;		/* Disk aspect ratio */
	float	disk_posang;		/* Disk position angle */
	double	z;			/* Galaxy redshift */
	float	hubble_type;		/* Hubble type */
        float	noiseqarea;		/* Galaxy image noise equivalent area */
	PIXTYPE	*subimage;		/* Sub-frame with the final object */
	int	subsize[2];		/* Sub-frame dimensions */
	int	subpos[2];		/* Sub-frame position in image */
	float	subfactor;		/* Sub-frame amplitude factor */
	float	*maskbuf;		/* Interpolation kernel buffer */
	int	*nmaskbuf;		/* Interpolation kernel size */
	int	*startbuf;		/* Interpolation kernel start */
	int	buf1size;		/* Size of first buffers series */
	float	*buf2;			/* Second (intermediary) buffer */
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
