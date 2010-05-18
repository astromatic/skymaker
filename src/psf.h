 /*
 				psf.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for psf.c
*
*	Last modify:	18/05/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SIMUL_H_
#include "simul.h"
#endif

/*---------------------------- Internal constants ---------------------------*/
/*-------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	center_psf(simstruct *sim),
		freepsf(simstruct *sim),
		makeaureole(simstruct *sim),
		makepsf(simstruct *sim),
		readpsf(simstruct *sim);

extern PIXTYPE	*interp_psf(simstruct *sim, double *pos, double *dpos),
		*interp_dft(simstruct *sim, int order, double *pos,
			double *dpos);

extern int	pos_to_indices(simstruct *sim, double *pos,
			int *index, PIXTYPE *weight);

