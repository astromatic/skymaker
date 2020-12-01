/*
*				makeit.c
*
* Main program
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		13/11/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include "define.h"
#include "globals.h"
#include "alterimage.h"
#include "corr.h"
#include "imaout.h"
#include "list.h"
#include "noise.h"
#include "prefs.h"
#include "psf.h"
#include "simul.h"
#include "stars.h"
#include "weight.h"

char		gstr[MAXCHAR];

/******************************** makeit *************************************/
/*
Manage the whole stuff.
*/
void    makeit() {

   simstruct	*simul;
   time_t	thetime, thetime2;
   struct tm	*tm;
   double	dtime;

/* Processing start date and time */
  thetime = time(NULL);
  dtime = counter_seconds();
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
  center_psf(simul);
  QCALLOC(simul->image, PIXTYPE, (size_t)simul->fimasize[0]*simul->fimasize[1]);
  if (prefs.weight_type != WEIGHT_NONE)
    weight_load(simul);
  list_openout(simul);
  list_read(simul);
  makestarfield(simul);
  list_closeout(simul);
  freepsf(simul);
  if (simul->aurange) {
    NFPRINTF(OUTPUT, "Adding diffuse component...");
    addaureole(simul);
  }

  if (prefs.weight_type == WEIGHT_NONE) {
    NFPRINTF(OUTPUT, "Adding background...");
    addback(simul);
  }

  if (simul->imatype != SKY_NONOISE && simul->imatype != GRID_NONOISE) {
/*-- Generate white noise */
    NFPRINTF(OUTPUT, "Generating noise...");
    noise_generate(simul);
    if (prefs.corr_type != CORR_NONE)
/*---- Correlate noise */
      corr_conv(simul, &simul->noise);
    noise_fix(simul);
  }

  if (prefs.corr_type == CORR_ALL)
/*-- Correlate the image itself */
    corr_conv(simul, &simul->image);

  if (simul->imatype != SKY_NONOISE && simul->imatype != GRID_NONOISE)
/*-- Apply noise to image */
    noise_apply(simul);


  NFPRINTF(OUTPUT, "Removing margins...");
  cutborders(simul, &simul->image);

//  quantify(&simul);

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
  imaout_write(simul);

  NFPRINTF(OUTPUT, "Freeing memory...");
  sim_end(simul);

/* Processing end date and time */
  thetime2 = time(NULL);
  tm = localtime(&thetime2);
  sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
  sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
  prefs.time_diff = counter_seconds() - dtime;

  return;
  }


/****** counter_seconds *******************************************************
PROTO	double counter_seconds(void)
PURPOSE	Count the number of seconds (with an arbitrary offset).
INPUT	-.
OUTPUT	Returns a number of seconds.
NOTES	Results are meaningful only for tasks that take one microsec or more.
AUTHOR	E. Bertin (IAP)
VERSION	24/09/2009
 ***/
double	counter_seconds(void)
  {
   struct timeval	tp;
   struct timezone	tzp;
   int			dummy;

  dummy = gettimeofday(&tp,&tzp);
  return (double) tp.tv_sec + (double) tp.tv_usec * 1.0e-6;
  }


