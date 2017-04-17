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
*	Last modified:		17/04/2017
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

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fft.h"
#include "prefs.h"
#include "random.h"

#ifdef USE_THREADS
#include "threads.h"

static void	pthread_addnoise(simstruct *sim),
		*pthread_addnoiseline(void *arg);

static pthread_t	*thread;
static pthread_mutex_t	noisemutex;

static simstruct	*pthread_sim;
static int		pthread_line;

#ifdef HAVE_MKL
static VSLStreamStatePtr	*pthread_stream;
static double			**pthread_lambdabuf;
static int			**pthread_poissonbuf;
static float			**pthread_gaussbuf;
#endif
#endif

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


/****** noise_make **************************************************
PROTO	void noise_make(simstruct *sim)
PURPOSE	Generate a photon and read-out noise map.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2017
 ***/
void	noise_make(simstruct *sim)

  {
#ifndef USE_THREADS
   PIXTYPE		*pix, *imapix, *noisepix, *rmspix, *imapixt,
			ron;
   int			i, npix;
 #ifdef HAVE_MKL
   VSLStreamStatePtr	stream;
   PIXTYPE		*pixt;
   double		*lambdabuf, *lambdabuft;
   int			*poissonbuf, *poissonbuft,
			j;
   float		*gaussbuf, *gaussbuft;
 #endif
#endif

  init_random(0);

  QMALLOC16(sim->noise, PIXTYPE, sim->imasize[0]*sim->imasize[1]);

#ifdef USE_THREADS
  pthread_addnoise(sim);
#else
  ron = (PIXTYPE)sim->ron;
  imapix = sim->image;
  noisepix = sim->noise;
  rmspix = sim->weight;
 #ifdef HAVE_MKL
  vslNewStream(&stream, VSL_BRNG_MT19937, (unsigned int)time(NULL));
  npix = sim->imasize[0];
  QMALLOC(lambdabuf, double, npix);
  QMALLOC(poissonbuf, int, npix);
  QMALLOC(gaussbuf, float, npix);
  for (j=sim->imasize[1]; j--; imapix += npix) {
    imapixt = imapix;
    lambdabuft = lambdabuf;
    for (i=npix; i--;)
      *(lambdabuft++) = (double)*(imapixt++);
    viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM, stream, npix,
		poissonbuf, lambdabuf);
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npix,
		gaussbuf, 0.0, 1.0);
    imapixt = imapix;
    poissonbuft = poissonbuf;
    gaussbuft = gaussbuf;
    if (sim->weight)
      rmspixt = rmspix;
      for (i=npix; i--;
		rmspix++, noisepix++, imapixt++, poissonbuft++, gaussbuft++) {
        if (isfinite(*rmspix) && *rmspix > SMALL)
          *noisepix = ((PIXTYPE)*poissonbuft - *imapixt) / *rmspix
					+ (PIXTYPE)*gaussbuft;
      }
    else
      for (i=npix; i--;)
        *(noisepix++) = (PIXTYPE)*(poissonbuft++) + ron * (PIXTYPE)*(gaussbuft++);
  }
  vslDeleteStream(&stream);
  free(lambdabuf);
  free(poissonbuf);
  free(gaussbuf);
 #else
  npix =  sim->imasize[0] * sim->imasize[1];
  if (sim->weight)
    for (i=sim->imasize[0] * sim->imasize[1]; i--;
		rmspix++, noisepix++, imapix++) {
      if (isfinite(*rmspix) && *rmspix > SMALL)
	*noisepix = (PIXTYPE)(random_poisson((double)*imapix, 0) - *imapix
		+ random_gauss(1.0, 0));
    }
  else
    for (i=sim->imasize[0] * sim->imasize[1]; i--; imapix++)
      *(noisepix++) = (PIXTYPE)(random_poisson((double)*imapix, 0) - *imapix
		+ random_gauss(ron, 0));
 #endif
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
  for (y = 0; y < sim->imasize[1]; y++) {
    noise_corrline(sim, corrnoise + y * sim->imasize[0], y);
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
PROTO	void noise_corrline(simstruct *sim, PIXTYPE *corrline, int y)
PURPOSE	Correlate a noise-map line to replicate image interpolation effects
INPUT	Pointer to the simulation,
	Pointer to the correlated (convolved) line data,
	y coordinate.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	17/04/2017
 ***/
void	noise_corrline(simstruct *sim, PIXTYPE *corrline, int y) {

   int		mw,mw2,m0,me,m,mx,dmx, y0, dy, sw,sh;
   PIXTYPE	*corrlinee, *corr,
		*s,*s0, *d,*de, mval;

  sw = sim->imasize[0];
  mw = sim->corrsize[0];
  mw2 = mw/2;
  y0 = y - (sim->corrsize[1]/2);
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

#ifdef USE_THREADS

/****** pthread_addnoiseline **************************************************
PROTO	void *pthread_addnoiseline(void *arg)
PURPOSE	Noise generation thread.
INPUT	Pointer to the thread number.
OUTPUT	NULL void pointer.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2017
 ***/
static void	*pthread_addnoiseline(void *arg)
  {
   double	ron;
   PIXTYPE	*pix, *pixt, *imapix, *noisepix, *rmspix;
   int		i, line, npix, nlines, proc,
		wflag = (pthread_sim->noise != NULL);
#ifdef HAVE_MKL
   VSLStreamStatePtr	stream;
   double		*lambdabuf, *lambdabuft;
   float		*gaussbuf, *gaussbuft;
   int			*poissonbuf, *poissonbuft;
#endif

  proc = *((int *)arg);

#ifdef HAVE_MKL
  stream = pthread_stream[proc];
  lambdabuf = pthread_lambdabuf[proc];
  poissonbuf = pthread_poissonbuf[proc];
  gaussbuf = pthread_gaussbuf[proc];
#endif

  ron = pthread_sim->ron;
  npix = pthread_sim->imasize[0];
  nlines = pthread_sim->imasize[1];
/* Exit if the end of image has been reached */
  QPTHREAD_MUTEX_LOCK(&noisemutex);  
  while (pthread_line < nlines)
    {
    line = pthread_line++;
    QPTHREAD_MUTEX_UNLOCK(&noisemutex);  
    noisepix = pthread_sim->noise + line*npix;
    rmspix = pthread_sim->weight + line*npix;
    imapix = pthread_sim->image + line*npix;
#ifdef HAVE_MKL
    lambdabuft = lambdabuf;
    for (i=npix; i--;)
      *(lambdabuft++) = *(imapix++);
    viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM, stream, npix,
		poissonbuf, lambdabuf);
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npix,
		gaussbuf, 0.0, 1.0);
    imapix = pthread_sim->image + line*npix;
    poissonbuft = poissonbuf;
    gaussbuft = gaussbuf;
    if (pthread_sim->weight)
      for (i=npix; i--;
		noisepix++, rmspix++, imapix++, poissonbuft++, gaussbuft++) {
        if (isfinite(*rmspix) && *rmspix > SMALL)
          *noisepix = ((PIXTYPE)*poissonbuft - *imapix) / *rmspix
					+ (PIXTYPE)*gaussbuft;
      }
    else
      for (i=npix; i--;)
        *(noisepix++) = (PIXTYPE)*(poissonbuft++) - *(imapix++)
		+ ron * (PIXTYPE)*(gaussbuft++);
#else
    if (pthread_sim->weight)
      for (i=npix; i--; noisepix++, imapix++) {
        if (isfinite(*rmspix) && *rmspix > SMALL)
          *noisepix = ((PIXTYPE)random_poisson((double)*imapix, 0) - *imapix)
			/ *rmpix + random_gauss(1.0, 0));
      }
    else
      for (i=npix; i--; imapix++)
        *(noisepix++) = (PIXTYPE)(random_poisson((double)*imapix, 0) - *imapix
			+ random_gauss(ron, 0));
#endif
    QPTHREAD_MUTEX_LOCK(&noisemutex);  
    }

  QPTHREAD_MUTEX_UNLOCK(&noisemutex);  

  pthread_exit(NULL);
  return (void *)NULL;
  }


/****** pthread_addnoise ****************************************************
PROTO	void pthread_addnoise(simstruct *sim)
PURPOSE	Add photon photo and read-out noise to image using multithreads
INPUT	Pointer to the sim structure.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	10/03/2016
 ***/
static void	pthread_addnoise(simstruct *sim)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				p;

/* Number of active threads */
  nproc = prefs.nthreads;
/* Set up multi-threading stuff */
  QPTHREAD_MUTEX_INIT(&noisemutex, NULL);
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_sim = sim;
  pthread_line = 0;
/* Start the reading/generation threads */

#ifdef HAVE_MKL
  QMALLOC(pthread_stream, VSLStreamStatePtr, nproc);
  QMALLOC(pthread_lambdabuf, double *, nproc);
  QMALLOC(pthread_poissonbuf, int *, nproc);
  QMALLOC(pthread_gaussbuf, float *, nproc);
#endif

  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
#ifdef HAVE_MKL
    vslNewStream(&pthread_stream[p], VSL_BRNG_MT19937, p);
    QMALLOC(pthread_lambdabuf[p], double, sim->imasize[0]);
    QMALLOC(pthread_poissonbuf[p], int, sim->imasize[0]);
    QMALLOC(pthread_gaussbuf[p], float, sim->imasize[0]);
#endif
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_addnoiseline, &proc[p]);
    }
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
/* Clean up multi-threading stuff */
  QPTHREAD_MUTEX_DESTROY(&noisemutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(proc);
  free(thread);
#ifdef HAVE_MKL
  for (p=0; p<nproc; p++)
    {
    vslDeleteStream(&pthread_stream[p]);
    free(pthread_lambdabuf[p]);
    free(pthread_poissonbuf[p]);
    free(pthread_gaussbuf[p]);
    }
  free(pthread_stream);
  free(pthread_lambdabuf);
  free(pthread_poissonbuf);
  free(pthread_gaussbuf);
#endif

  return;
  }

#endif

