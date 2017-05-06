/*
*				corr.c
*
* Manage noise correlation.
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
*	Last modified:		05/05/2017
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

#include FFTW_H

#ifdef HAVE_MKL
 #include MKL_H
#endif

#ifdef USE_THREADS
  #include <omp.h>
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

static const corrstruct corrs[] = {
  {CORRFUNC_NEAREST, 0.5, corr_func_nearest},
  {CORRFUNC_BILINEAR, 1.0, corr_func_bilinear},
  {CORRFUNC_LANCZOS2, 2.0, corr_func_lanczos2},
  {CORRFUNC_LANCZOS3, 3.0, corr_func_lanczos3},
  {CORRFUNC_LANCZOS4, 4.0, corr_func_lanczos4}
};


/****** corr_conv ******************************************************
PROTO	void corr_conv(simstruct *sim, PIXTYPE **image)
PURPOSE	Correlate (convolve) an image with some symmetric kernel.
INPUT	Pointer to the simulation,
	pointer to the input image pointer.
OUTPUT	-.
NOTES	The input image pointer is freed and replaced with that of the new
	correlated image. 
AUTHOR	E. Bertin (IAP)
VERSION	04/05/2017
 ***/
void	corr_conv(simstruct *sim, PIXTYPE **image) {

   PIXTYPE	*imageout;
   int		y;

  QMALLOC16(imageout, PIXTYPE, sim->fimasize[0]*sim->fimasize[1]);

#pragma omp parallel for num_threads(prefs.nthreads)
  for (y = 0; y < sim->fimasize[1]; y++)
    corr_convline(sim, *image, imageout, y);

  free(*image);
  *image = imageout;

  return;
  }


/****** corr_convline **************************************************
PROTO	void corr_convline(simstruct *sim, PIXTYPE *imagein, PIXTYPE *imageout,
		int y)
PURPOSE	Correlate an image line to replicate image interpolation effects
INPUT	Pointer to the simulation,
	pointer to the input image,
	pointer to the output correlated (convolved) image,
	y coordinate.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/05/2017
 ***/
void	corr_convline(simstruct *sim, PIXTYPE *imagein, PIXTYPE *imageout,
		int y) {

   int		mw,mw2,m0,me,m,mx,dmx, y0, dy, sw,sh;
   PIXTYPE	*corrline, *corrlinee, *corr,
		*s,*s0, *d,*de, mval;

  sw = sim->fimasize[0];
  mw = sim->corrsize[0];
  mw2 = mw/2;
  y0 = y - (sim->corrsize[1]/2);
  corrline = imageout + y * sw;
  if ((dy = -y0) > 0)
    {
    m0 = mw*dy;
    y0 = 0;
    }
  else
    m0 = 0;

  if ((dy = sim->fimasize[1] - y0) < sim->corrsize[1])
    me = mw * dy;
  else
    me = mw * sim->corrsize[1];

  corrlinee = corrline + sw;

  memset(corrline, 0, sw * sizeof(PIXTYPE));
  
  s0 = NULL;				/* To avoid gcc -Wall warnings */
  corr = sim->corr + m0;
  for (m = m0, mx = 0; m<me; m++, mx++) {
    if (mx==mw)
      mx = 0;
    if (!mx)
      s0 = imagein + sw * y0++;

    if ((dmx = mx - mw2) >= 0) {
      s = s0 + dmx;
      d = corrline;
      de = corrlinee - dmx;
    } else {
      s = s0;
      d = corrline - dmx;
      de = corrlinee;
    }

    mval = *(corr++);
    while (d < de)
      *(d++) += mval * *(s++);
  }

  return;
}


/****** corr_generate **************************************************
PROTO	void corr_generate(simstruct *sim, corrfuncenum corrfunc_type,
		float scale)
PURPOSE	Generate average correlation square root function for the provided
	interpolation function.
INPUT	Pointer to the simulation,
	correlation function type,
	relative pixel scale (>1 if the output grid oversamples the input grid)
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	05/05/2017
 ***/
void	corr_generate(simstruct *sim, corrfuncenum corrfunc_type,
		float scale) {

   corrstruct	corr;
   PIXTYPE	(*corr_func)(float x),
		*kernel;
   double	dval, dsum;
   float	invscale, dx, dy;
   int		i, imax, size, ix, iy;

  corr = corrs[corrfunc_type];
  imax = (int)(scale * (corr.radius + 0.499));
  sim->corrsize[0] = sim->corrsize[1] = size = imax * 2 + 1;
  invscale = 1.0 / scale;
  corr_func = corr.func;

/* Allocate memory for the correlation kernel */
  QMALLOC16(sim->corr, PIXTYPE, sim->corrsize[0] * sim->corrsize[1]);

/* Integrate over all kernel positions with respect to the pixel grid */
  kernel = sim->corr + imax * (size + 1);
  dsum = 0.0;

#pragma omp parallel for collapse(2) private(dval,dx,dy) reduction(+:dsum) num_threads(prefs.nthreads)
  for (iy = -imax; iy <= imax; iy++)
    for (ix = -imax; ix <= imax; ix++) {
      dval = 0.0;
      for (dy = -0.5; dy <= 0.5; dy += CORRFUNC_STEP)
        for (dx = -0.5; dx <= 0.5; dx += CORRFUNC_STEP)
               dval += corr_func((ix+dx) * invscale) *
			corr_func((iy+dy) * invscale);
      kernel[iy * size + ix] = (PIXTYPE)dval;
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

  x = fabsf(x);

  return x < 1.0f ? 1.0f - x : 0.0f;
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


