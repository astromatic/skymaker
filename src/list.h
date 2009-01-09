 /*
 				list.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for list.c.
*
*	Last modify:	17/08/2006
**
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
	double	x,y;			/* Position of the object center */
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
	int	subx,suby;		/* Sub-frame position in image */
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
