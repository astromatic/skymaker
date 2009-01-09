 /*
 				alterimage.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP).
*
*	Contents:	Include for stars.c
*
*	Last modify:	20/09/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SIMUL_H_
#include "simul.h"
#endif

/*---------------------------- Internal constants ---------------------------*/

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	addaureole(simstruct *sim),
		addback(simstruct *sim),
		addnoise(simstruct *sim),
		cutborders(simstruct *sim),
		etoadu(simstruct *sim),
		microscan(simstruct *sim),
		saturate(simstruct *sim);
