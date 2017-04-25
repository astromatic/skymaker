/*
*				corr.c
*
* Generate correlation square root function.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2017 IAP/CNRS/UPMC
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
*	Last modified:		25/04/2017
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

#include FFTW_H

#ifdef HAVE_MKL
 #include MKL_H
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "prefs.h"
#include "corr.h"
#include "image.h"
#include "simul.h"

static	PIXTYPE	corr_func_nearest(float x),
		corr_func_bilinear(float x),
		corr_func_lanczos2(float x),
		corr_func_lanczos3(float x),
		corr_func_lanczos4(float x);

/****** corr_generate **************************************************
PROTO	void corr_generate(simstruct *sim, corrinterpenum interp_type,
		float scale)
PURPOSE	Generate average correlation square root function for the provided
	interpolation function.
INPUT	Pointer to the simulation,
	interpolant type,
	relative pixel scale (>1 if the output grid oversamples the input grid)
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/04/2017
 ***/
void	corr_generate(simstruct *sim, corrinterpenum corrinterp_type,
		float scale) {

  static const corrstruct coors[6] = {
    {CORRINTERP_NONE, 0.0,  NULL},
    {CORRINTERP_NEAREST, 0.5, corr_func_nearest},
    {CORRINTERP_BILINEAR, 1.0, corr_func_bilinear},
    {CORRINTERP_LANCZOS2, 2.0, corr_func_lanczos2},
    {CORRINTERP_LANCZOS3, 3.0, corr_func_lanczos3},
    {CORRINTERP_LANCZOS4, 4.0, corr_func_lanczos4}
  };

   corrstruct	corr;
   PIXTYPE	(*corr_func)(float x),
		*kernel;
   double	dval, dsum;
   float	invscale, dx, dy;
   int		i, imax, ix, iy;

  if (corrinterp_type == CORRINTERP_NONE)
    return;

  corr = coors[corrinterp_type];
  imax = (int)(scale * (corr.radius + 0.499));
  sim->corrsize[0] = sim->corrsize[1] = imax * 2 + 1;
  invscale = 1.0 / scale;
  corr_func = corr.func;

/* Allocate memory for the correlation kernel */
  QMALLOC16(sim->corr, PIXTYPE, sim->corrsize[0] * sim->corrsize[1]);

/* Integrate over all kernel positions with respect to the pixel grid */
  kernel = sim->corr;
  dsum = 0.0;
  for (iy = -imax; iy <= imax; iy++)
    for (ix = -imax; ix <= imax; ix++) {
      dval = 0.0;
      for (dy = -0.5; dy <= 0.5; dy += CORRINTERP_STEP)
        for (dx = -0.5; dx <= 0.5; dx += CORRINTERP_STEP)
               dval += corr_func((ix+dx) * invscale) *
			corr_func((iy+dy) * invscale);
    *(kernel++) = (PIXTYPE)dval;
    dsum += dval;
    }

/* Normalize kernel sum to 1 */
  invscale = 1.0 / dsum;
  kernel = sim->corr;
  for (i = sim->corrsize[0] * sim->corrsize[1]; i--;)
    *(kernel++) *= invscale;

  return;
  }


/****** corr_func_nearest **************************************************
PROTO	PIXTYPE corr_func_nearest(float x)
PURPOSE	Generate value for the nearest neighbour interpolation function.
INPUT	x coordinate.
OUTPUT	Interpolant value at x.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2017
 ***/
PIXTYPE corr_func_nearest(float x) {

  return fabsf(x) < 0.5f ? 1.0f : 0.0f;
}


/****** corr_func_bilinear *************************************************
PROTO	PIXTYPE corr_func_bilinear(float x)
PURPOSE	Generate value for the linear interpolation function.
INPUT	x coordinate.
OUTPUT	Interpolant value at x.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2017
 ***/
PIXTYPE corr_func_bilinear(float x) {

  return fabsf(x) < 1.0f ? 1.0f - x : 0.0f;
}


/****** corr_func_lanczos2 *************************************************
PROTO	PIXTYPE corr_func_lanczos2(float x)
PURPOSE	Generate value for the Lanczos2 interpolation function.
INPUT	x coordinate.
OUTPUT	Interpolant value at x.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2017
 ***/
PIXTYPE corr_func_lanczos2(float x) {

  x = fabsf(x);
  if (x < 1e-5f)
    return 1.0f;
  else if (x < 2.0f) {
    x *= PI;
    return 2.0f * sinf(x) * sinf(0.5f * x) / (x * x);
  } else
    return 0.0f;
}


/****** corr_func_lanczos3 *************************************************
PROTO	PIXTYPE corr_func_lanczos3(float x)
PURPOSE	Generate value for the Lanczos2 interpolation function.
INPUT	x coordinate.
OUTPUT	Interpolant value at x.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2017
 ***/
PIXTYPE corr_func_lanczos3(float x) {

  x = fabsf(x);
  if (x < 1e-5f)
    return 1.0f;
  else if (x < 3.0f) {
    x *= PI;
    return 3.0f * sinf(x) * sinf(0.333333f * x) / (x * x);
  } else
    return 0.0f;
}


/****** corr_func_lanczos4 *************************************************
PROTO	PIXTYPE corr_func_lanczos4(float x)
PURPOSE	Generate value for the Lanczos4 interpolation function.
INPUT	x coordinate.
OUTPUT	Interpolant value at x.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/04/2017
 ***/
PIXTYPE corr_func_lanczos4(float x) {

  x = fabsf(x);
  if (x < 1e-5f)
    return 1.0f;
  else if (x < 4.0f) {
    x *= PI;
    return 4.0f * sinf(x) * sinf(0.25f * x) / (x * x);
  } else
    return 0.0f;
}


