/*
*				noise.c
*
* Add noise to images.
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
*	Last modified:		24/04/2017
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
#include <time.h>

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
#include "corr.h"
#include "fft.h"
#include "noise.h"
#include "prefs.h"
#include "random.h"
#include "simul.h"

/****** noise_add **************************************************
PROTO	void noise_add(simstruct *sim)
PURPOSE	Add noise map to image.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2017
 ***/
void	noise_add(simstruct *sim)

  {
   PIXTYPE	*imapix, *rmspix;
   int		i;

  imapix = sim->image;
  rmspix = sim->noise;
  for (i = sim->imasize[0] * sim->imasize[1]; i--; imapix++, rmspix++)
    *imapix = isfinite(*rmspix) ? *imapix + *rmspix : 0.0;

  return;
  }


/****** noise_generate **************************************************
PROTO	void noise_generate(simstruct *sim)
PURPOSE	Generate a photon and read-out noise map.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2017
 ***/
void	noise_generate(simstruct *sim) {

#ifdef HAVE_MKL
   unsigned int	seed;
   int		p;
#endif

   int		y;

  QMALLOC16(sim->noise, PIXTYPE, sim->imasize[0]*sim->imasize[1]);

#ifdef HAVE_MKL
/* Allocate buffer memory */
  seed = (unsigned int)time(NULL);
  QCALLOC(sim->streams,  void *, prefs.nthreads);
  QCALLOC(sim->lambdabuf, double *, prefs.nthreads);
  QCALLOC(sim->poissonbuf, int *, prefs.nthreads);
  QCALLOC(sim->gaussbuf, float *, prefs.nthreads);
  for (p = 0; p<prefs.nthreads; p++) {
    vslNewStream((VSLStreamStatePtr *)&sim->streams[p], VSL_BRNG_MT2203,
		seed + p);
  QMALLOC16(sim->lambdabuf[p], double, sim->imasize[0]);
  QMALLOC16(sim->poissonbuf[p], int, sim->imasize[0]);
  QMALLOC16(sim->gaussbuf[p], float, sim->imasize[0]);
  }
#else
  init_random(0);
#endif

#pragma omp parallel for num_threads(prefs.nthreads)
  for (y = 0; y < sim->imasize[1]; y++)
    noise_generateline(sim, sim->noise, y);

#ifdef HAVE_MKL
/* Free buffer memory */
  for (p = 0; p<prefs.nthreads; p++) {
    vslDeleteStream((VSLStreamStatePtr *)&sim->streams[p]);
    free(sim->lambdabuf[p]);
    free(sim->poissonbuf[p]);
    free(sim->gaussbuf[p]);
  }
  free(sim->streams);
  free(sim->lambdabuf);
  free(sim->poissonbuf);
  free(sim->gaussbuf);
#endif

}


/****** noise_generateline ************************************************
PROTO	void noise_generateline(simstruct *sim, PIXTYPE *noise, int y)
PURPOSE	Generate lines of photon and read-out noise.
INPUT	Pointer to the simulation,
	pointer to the noise data,
	y coordinate.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2017
 ***/
void	noise_generateline(simstruct *sim, PIXTYPE *noise,  int y) {

   PIXTYPE		*pix, *imapix, *noisepix, *rmspix,
			ron;
   int			i, p, npix;

#ifdef HAVE_MKL
   PIXTYPE		*imapixt;
   double		*lambdabuf, *lambdabuft;
   int			*poissonbuf;
   float		*gaussbuf;
#endif

#ifdef USE_THREADS
  p = omp_get_thread_num();
#else
  p = 0;
#endif

  npix = sim->imasize[0];
  ron = (PIXTYPE)sim->ron;
  imapix = sim->image + y * npix;
  noisepix = sim->noise + y * npix;
  rmspix = sim->weight + y * npix;

#ifdef HAVE_MKL
  imapixt = imapix;
  lambdabuf = sim->lambdabuf[p];
  for (i=npix; i--;)
    *(lambdabuf++) = (double)*(imapixt++);
  viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM,
		(VSLStreamStatePtr)sim->streams[p], npix,
		sim->poissonbuf[p], sim->lambdabuf[p]);
  vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,
		(VSLStreamStatePtr)sim->streams[p], npix,
		sim->gaussbuf[p], 0.0, 1.0);
  imapixt = imapix;
  poissonbuf = sim->poissonbuf[p];
  gaussbuf = sim->gaussbuf[p];
  if (sim->weight) {
    for (i=npix; i--;
		rmspix++, noisepix++, imapix++, poissonbuf++, gaussbuf++) {
      if (isfinite(*rmspix) && *rmspix > SMALL)
        *noisepix = ((PIXTYPE)*poissonbuf - *imapix) / *rmspix
					+ (PIXTYPE)*gaussbuf;
    }
  } else
    for (i=npix; i--;)
      *(noisepix++) = (PIXTYPE)*(poissonbuf++) + ron * (PIXTYPE)*(gaussbuf++);
#else
  if (sim->weight)
    for (i=npix; i--; rmspix++, noisepix++, imapix++) {
      if (isfinite(*rmspix) && *rmspix > SMALL)
	*noisepix = (PIXTYPE)(random_poisson((double)*imapix, 0) - *imapix
		+ random_gauss(1.0, 0));
    }
  else
    for (i=npix; i--; imapix++)
      *(noisepix++) = (PIXTYPE)(random_poisson((double)*imapix, 0) - *imapix
		+ random_gauss(ron, 0));
#endif

  return;
}


/****** noise_corr **************************************************
PROTO	void noise_corr(simstruct *sim)
PURPOSE	Correlate (convolve) a noise map with some kernel.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	17/04/2017
 ***/
void	noise_corr(simstruct *sim) {

   PIXTYPE	*corrnoise, *corrnoisepix, *rmspix, *noisepix;
   int		i, y;

  QMALLOC16(corrnoise, PIXTYPE, sim->imasize[0]*sim->imasize[1]);

#pragma omp parallel for num_threads(prefs.nthreads)
  for (y = 0; y < sim->imasize[1]; y++) {
    noise_corrline(sim, corrnoise, y);
  }

  noisepix = sim->noise;
  corrnoisepix = corrnoise;
  if (sim->weight) {
    rmspix = sim->weight;
    for (i = sim->imasize[0] * sim->imasize[1]; i--;
		rmspix++, noisepix++, corrnoisepix++)
      if (isfinite(*rmspix) && *rmspix > SMALL)
        *noisepix = *corrnoisepix * *rmspix;
  } else {
    for (i = sim->imasize[0] * sim->imasize[1]; i--;)
      *(noisepix++) = *(corrnoisepix++);
  }

  free(corrnoise);

  return;
}


/****** noise_corrline **************************************************
PROTO	void noise_corrline(simstruct *sim, PIXTYPE *corrnoise, int y)
PURPOSE	Correlate a noise-map line to replicate image interpolation effects
INPUT	Pointer to the simulation,
	pointer to the correlated (convolved) data array,
	y coordinate.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2017
 ***/
void	noise_corrline(simstruct *sim, PIXTYPE *corrnoise, int y) {

   int		mw,mw2,m0,me,m,mx,dmx, y0, dy, sw,sh;
   PIXTYPE	*corrline, *corrlinee, *corr,
		*s,*s0, *d,*de, mval;

  sw = sim->imasize[0];
  mw = sim->corrsize[0];
  mw2 = mw/2;
  y0 = y - (sim->corrsize[1]/2);
  corrline = corrnoise + y * sw;
  if ((dy = -y0) > 0)
    {
    m0 = mw*dy;
    y0 = 0;
    }
  else
    m0 = 0;

  if ((dy = sim->imasize[1] - y0) < sim->corrsize[1])
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
      s0 = sim->noise + sw * y0++;

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

