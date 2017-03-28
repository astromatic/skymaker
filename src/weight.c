/*
*				weight.c
*
* Manage input weight-maps.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2000-2017 IAP/CNRS/UPMC
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
*	Last modified:		28/03/2017
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "prefs.h"
#include "weight.h"

union {unsigned int i; float f;} anan = {0x7ff80000};

/******* load_weight *********************************************************
PROTO	int load_weight(simstruct *sim)
PURPOSE	Load a weight-map and convert it to a noise (RMS) map.
INPUT	Pointer to the simulation.
OUTPUT	RETURN_OK if no error, or RETURN_ERROR in case of non-fatal error(s).
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 28/03/2017
 ***/
int	load_weight(simstruct *sim)

  {
   catstruct	*cat;
   tabstruct	*tab;
   PIXTYPE	*pix,
		weight_fac, weight_thresh;
   char		lstr[MAXCHAR],
		*filename, *rfilename, *str,*str2;
   weightenum	wtype;
   int		p, npix;

  wtype = prefs.weight_type;

  filename = prefs.weight_name;
  if ((str = strrchr(filename, '[')))
    {
    *str = '\0';
    if ((str2 = strrchr(str, ']')))
      *str2 = '\0';
    }

/* A short, "relative" version of the filename */
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  sprintf(lstr,"Examining %s", rfilename);
  NFPRINTF(OUTPUT, lstr);
  if (!(cat=read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: no FITS data found in ", filename);

  if (str)
    {
    if (!(tab = name_to_tab(cat, str+1, 0)))
      tab = pos_to_tab(cat, atoi(str+1), 0);
    }
  else
    tab = cat->tab;

/*-- Force the data to be at least 2D */
  if (tab->naxis<2)
    error(EXIT_FAILURE, "*Error*: no 2D FITS data in ", filename);
  if (tab->naxisn[0] != sim->imasize[0] || tab->naxisn[1] != sim->imasize[1])
    error(EXIT_FAILURE, "*Error*: weight-map dimensions do not match in ", filename);
  npix = sim->imasize[0] * sim->imasize[1];
  sprintf(lstr,"Loading %s", rfilename);
  NFPRINTF(OUTPUT, lstr);
  QMALLOC16(sim->noise, PIXTYPE, npix);
  QFSEEK(tab->cat->file, tab->bodypos, SEEK_SET, tab->cat->filename);
  read_body(tab, sim->noise, npix);
  free_cat(&cat, 1);

  weight_fac = prefs.weight_factor;

  switch(prefs.weight_type)
    {
    case WEIGHT_NONE:
      return RETURN_ERROR;

    case WEIGHT_FROMRMSMAP:
      weight_thresh = prefs.nweight_thresh ? prefs.weight_thresh : BIG;
      pix = sim->noise;
      for (p=npix; p--; pix++)
//---- Negative variance means "bad pixel"
        *pix = *pix < weight_thresh ? weight_fac * *pix : anan.f;
      break;

    case WEIGHT_FROMVARMAP:
      weight_thresh = prefs.nweight_thresh ? prefs.weight_thresh : BIG;
      pix = sim->noise;
      for (p=npix; p--; pix++)
//---- Negative variance means "bad pixel"
        *pix = *pix < weight_thresh ? weight_fac * sqrtf(fabsf(*pix)) : anan.f;
      break;

    case WEIGHT_FROMWEIGHTMAP:
      weight_thresh = prefs.nweight_thresh ? prefs.weight_thresh : 0.0;
      pix = sim->noise;
      for (p=npix; p--; pix++)
//---- Negative variance means "bad pixel"
        *pix = *pix > weight_thresh ? weight_fac / sqrtf(fabs(*pix)) : anan.f;
      break;

    default:
      error(EXIT_FAILURE,
	"*Internal Error*: Unknown weight-map type in ", "load_weight()");
      break;
    }

  return RETURN_OK;
  }


