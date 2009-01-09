/*
                                  makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SkyMaker
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Main program.
*
*       Last modify:    27/09/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "alterimage.h"
#include "imaout.h"
#include "list.h"
#include "prefs.h"
#include "psf.h"
#include "simul.h"
#include "stars.h"

/******************************** makeit *************************************/
/*
Manage the whole stuff.
*/
void    makeit()

  {
   simstruct	*simul;
   time_t	thetime, thetime2;
   struct tm	*tm;

  thetime = time(NULL);	/* Record the time at beginning of sim */
  tm = localtime(&thetime);
  sprintf(prefs.sdate_start,"%04d-%02d-%02d",
	tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_start,"%02d:%02d:%02d",
	tm->tm_hour, tm->tm_min, tm->tm_sec);

  NFPRINTF(OUTPUT, "");
  QPRINTF(OUTPUT, "----- %s %s started on %s at %s with %d thread%s\n\n",
	BANNER, MYVERSION,
	prefs.sdate_start, prefs.stime_start,
	prefs.nthreads, prefs.nthreads>1? "s":"");

  NFPRINTF(OUTPUT, "Initializing simulation...");
  simul = sim_init();
  NFPRINTF(OUTPUT, "Creating the PSF...");
  if (simul->psftype == PSF_FILE)
    readpsf(simul);
  else
    makepsf(simul);
  QCALLOC(simul->image, PIXTYPE, simul->fimasize[0]*simul->fimasize[1]);
  openoutlist(simul);
  makestarfield(simul);
  readlist(simul);
  closeoutlist(simul);
  freepsf(simul);
  if (simul->aurange)
    {
    NFPRINTF(OUTPUT, "Adding diffuse component...");
    addaureole(simul);
    NFPRINTF(OUTPUT, "Removing margins...");
    cutborders(simul);
    }
  NFPRINTF(OUTPUT, "Adding background...");
  addback(simul);
  if (simul->imatype != SKY_NONOISE && simul->imatype != GRID_NONOISE)
    {
    NFPRINTF(OUTPUT, "Adding noise...");
    addnoise(simul);
    }
/*
  quantify(&simul);
*/
  if (simul->wellcap)
    {
    NFPRINTF(OUTPUT, "Simulating saturation...");
    saturate(simul);
    }
  NFPRINTF(OUTPUT, "Converting to ADUs...");
  etoadu(simul);
  NFPRINTF(OUTPUT, "Microscanning image...");
  microscan(simul);
  NFPRINTF(OUTPUT, "Writing image...");
  writeima(simul);
  NFPRINTF(OUTPUT, "Freeing memory...");
  sim_end(simul);

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = difftime(thetime2, thetime);

  return;
  }


