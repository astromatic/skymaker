 /*
 				image.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for image.c.
*
*	Last modify:	23/09/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

