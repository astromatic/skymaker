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
*	Last modified:		04/05/2017
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
VERSION	04/05/2017
 ***/
void	noise_add(simstruct *sim)

  {
   PIXTYPE	*imapix, *noisepix;
   int		i;

  imapix = sim->image;
  noisepix = sim->noise;
  for (i = sim->fimasize[0] * sim->fimasize[1]; i--; imapix++, noisepix++)
    if (isfinite(*noisepix))
      *imapix += *noisepix;

  return;
  }


/****** noise_generate **************************************************
PROTO	void noise_generate(simstruct *sim)
PURPOSE	Generate a photon and read-out noise map.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	04/05/2017
 ***/
void	noise_generate(simstruct *sim) {

#ifdef HAVE_MKL
   unsigned int	seed;
#endif

   PIXTYPE	*rmspix;
   int		i, p, y;

  QMALLOC16(sim->noise, PIXTYPE, sim->fimasize[0]*sim->fimasize[1]);
  QCALLOC(sim->weightbuf, PIXTYPE *, sim->fimasize[0] * prefs.nthreads);
  for (p = 0; p<prefs.nthreads; p++) {
    QMALLOC16(sim->weightbuf[p], PIXTYPE, sim->fimasize[0]);
    rmspix = sim->weightbuf[p];
    for (i = sim->fimasize[0]; i--;)
      *(rmspix++) = 1.0;
  }
    

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
    QMALLOC16(sim->lambdabuf[p], double, sim->fimasize[0]);
    QMALLOC16(sim->poissonbuf[p], int, sim->fimasize[0]);
    QMALLOC16(sim->gaussbuf[p], float, sim->fimasize[0]);
  }
#else
  init_random(0);
#endif

#pragma omp parallel for num_threads(prefs.nthreads)
  for (y = 0; y < sim->fimasize[1]; y++)
    noise_generateline(sim, sim->noise, y);

/* Free buffer memory */
  for (p = 0; p<prefs.nthreads; p++)
    free(sim->weightbuf[p]);
  free(sim->weightbuf);

#ifdef HAVE_MKL
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

  return;
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
VERSION	04/05/2017
 ***/
void	noise_generateline(simstruct *sim, PIXTYPE *noise,  int y) {

   PIXTYPE		*pix, *imapix, *noisepix, *rmspix,
			ron;
   int			i, p, npix;

#ifdef HAVE_MKL
   double		*lambdabuf, *lambdabuft;
   int			*poissonbuf;
   float		*gaussbuf;
#endif

#ifdef USE_THREADS
  p = omp_get_thread_num();
#else
  p = 0;
#endif

  npix = sim->fimasize[0];
  ron = (PIXTYPE)sim->ron;
  imapix = sim->image + y * npix;
  noisepix = sim->noise + y * npix;
  if (sim->weight) {
/* Weightmaps have no margin! */
    memcpy(sim->weightbuf[p] + sim->margin[0], sim->weight + y * npix, sim->imasize[0] * sizeof(PIXTYPE));
    rmspix = sim->weightbuf[p];
  }

#ifdef HAVE_MKL
  lambdabuf = sim->lambdabuf[p];
  for (i=npix; i--;)
    *(lambdabuf++) = (double)*(imapix++);
  viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM,
		(VSLStreamStatePtr)sim->streams[p], npix,
		sim->poissonbuf[p], sim->lambdabuf[p]);
  vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,
		(VSLStreamStatePtr)sim->streams[p], npix,
		sim->gaussbuf[p], 0.0, 1.0);
  imapix = sim->image + y * npix;
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
      *(noisepix++) = (PIXTYPE)*(poissonbuf++) - *(imapix++) +
			ron * (PIXTYPE)*(gaussbuf++);
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


/****** noise_corr ******************************************************
PROTO	void noise_corr(simstruct *sim)
PURPOSE	Correlate (convolve) a noise map with some kernel.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/05/2017
 ***/
void	noise_corr(simstruct *sim) {

   PIXTYPE	*rmspix, *noisepix;
   int		i;

  corr_conv(sim, &sim->noise);

  noisepix = sim->noise;
  if (sim->weight) {
    rmspix = sim->weight;
    for (i = sim->fimasize[0] * sim->fimasize[1]; i--; rmspix++, noisepix++)
      if (isfinite(*rmspix) && *rmspix > SMALL)
        *noisepix *= *rmspix;
  }

  return;
}


