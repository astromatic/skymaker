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
*	Last modified:		17/03/2017
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

/****** add_noise **************************************************
PROTO	void make_noise(simstruct *sim)
PURPOSE	Add noise map to image.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2017
 ***/
void	add_noise(simstruct *sim)

  {
   PIXTYPE	*imapix, *rmspix;
   int		i;

  imapix = sim->image;
  rmspix = sim->noise;
  for (i = sim->imasize[0] * sim->imasize[1]; i--; imapix++, rmspix++)
    *imapix = isfinite(*rmspix) ? *imapix + *rmspix : 0.0;

  return;
  }


/****** make_noise **************************************************
PROTO	void make_noise(simstruct *sim)
PURPOSE	Generate a photon and read-out noise map.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	28/03/2017
 ***/
void	make_noise(simstruct *sim)

  {
#ifndef USE_THREADS
   PIXTYPE		*pix, *imapix, *rmspix, *imapixt, *pixt,
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
   int			wflag;

  init_random(0);

  if ((sim->noise))
    wflag = 1;
  else {
    wflag = 0;
    QMALLOC16(sim->noise, PIXTYPE, sim->imasize[0]*sim->imasize[1]);
  }
#ifdef USE_THREADS
  pthread_addnoise(sim);
#else
  ron = (PIXTYPE)sim->ron;
  imapix = sim->image;
  rmspix = sim->noise;
 #ifdef HAVE_MKL
  vslNewStream(&stream, VSL_BRNG_MT19937, (unsigned int)time(NULL));
  npix = sim->imasize[0];
  QMALLOC(lambdabuf, double, npix);
  QMALLOC(poissonbuf, int, npix);
  QMALLOC(gaussbuf, float, npix);
  for (j=sim->imasize[1]; j--; imapix += npix, rmspix += npix) {
    imapixt = imapix;
    lambdabuft = lambdabuf;
    for (i=npix; i--;)
      *(lambdabuft++) = (double)*(imapixt++);
    viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM, stream, npix,
		poissonbuf, lambdabuf);
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npix,
		gaussbuf, 0.0, 1.0);
    imapixt = imapix;
    pixt = rmspix;
    poissonbuft = poissonbuf;
    gaussbuft = gaussbuf;
    if (wflag)
      for (i=npix; i--; pixt++, imapixt++, poissonbuft++, gaussbuft++) {
        if (isfinite(*pixt))
          *pixt = (PIXTYPE)*poissonbuft - *imapixt
					+ *pixt * (PIXTYPE)*gaussbuft;
      }
    else
      for (i=npix; i--;)
        *(pixt++) = (PIXTYPE)*(poissonbuft++) + ron * (PIXTYPE)*(gaussbuft++);
  }
  vslDeleteStream(&stream);
  free(lambdabuf);
  free(poissonbuf);
  free(gaussbuf);
 #else
  npix =  sim->imasize[0] * sim->imasize[1];
  imapixt = imapix;
  pixt = rmspix;
  if (wflag)
    for (i=sim->imasize[0] * sim->imasize[1]; i--; pixt++, imapixt++) {
      if (isfinite(*pixt))
	*pixt = (PIXTYPE)(random_poisson((double)*imapixt, 0) - *imapixt
		+ random_gauss(*pixt, 0));
    }
  else
    for (i=sim->imasize[0] * sim->imasize[1]; i--; imapixt++)
      *(pixt++) = (PIXTYPE)(random_poisson((double)*imapixt, 0) - *imapixt
		+ random_gauss(ron, 0));
 #endif
#endif

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
   PIXTYPE	*pix, *pixt, *imapixt;
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
    imapixt = pthread_sim->image + line*npix;
    pixt = pthread_sim->noise + line*npix;
#ifdef HAVE_MKL
    lambdabuft = lambdabuf;
    for (i=npix; i--;)
      *(lambdabuft++) = *(imapixt++);
    viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM, stream, npix,
		poissonbuf, lambdabuf);
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npix,
		gaussbuf, 0.0, 1.0);
    imapixt = pthread_sim->image + line*npix;;
    poissonbuft = poissonbuf;
    gaussbuft = gaussbuf;
    if (wflag)
      for (i=npix; i--; pixt++, imapixt++, poissonbuft++, gaussbuft++) {
        if (isfinite(*pixt))
          *pixt = (PIXTYPE)*poissonbuft - *imapixt
					+ *pixt * (PIXTYPE)*gaussbuft;
      }
    else
      for (i=npix; i--;)
        *(pixt++) = (PIXTYPE)*(poissonbuft++) - *(imapixt++)
		+ ron * (PIXTYPE)*(gaussbuft++);
#else
    if (wflag)
      for (i=npix; i--; pixt++, imapixt++) {
        if (isfinite(*pixt))
          *pixt = (PIXTYPE)(random_poisson((double)*imapixt, 0) - *imapixt
			+ random_gauss(*pixt, 0));
      }
    else
      for (i=npix; i--; imapixt++)
        *(pixt++) = (PIXTYPE)(random_poisson((double)*imapixt, 0) - *imapixt
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

