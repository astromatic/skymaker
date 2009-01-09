/*
 				imaout.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	header structure and templates for catalog data.
*
*	Last modify:	28/09/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _PREFS_H_
#include "prefs.h"
#endif

#ifndef _SIMUL_H_
#include "simul.h"
#endif

extern char	prefsname[MAXCHAR];

/*--------------------------- structure definitions -------------------------*/

typedef struct structheadkey
  {
  char          name[80];               /* name */
  char          comment[80];            /* a comment */
  void          *ptr;                   /* pointer to the data */
  h_type        htype;                  /* standard ``h_type'' (display) */
  t_type        ttype;                  /* standard ``t_type'' (storage) */
  char		prefsname[80];		/* true name in the preferences */
  }             headkeystruct;

/*---------------------------------- protos --------------------------------*/
int	read_aschead(char *filename, int frameno, tabstruct *tab);
void	writeima(simstruct *sim);
