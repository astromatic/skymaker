 /*
 				galaxies.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for galaxies.c
*
*	Last modify:	15/08/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
