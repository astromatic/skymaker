/*
*				alterimage.c
*
* Alter images as a whole.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2016 IAP/CNRS/UPMC
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
*	Last modified:		10/03/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include FFTW_H

#ifdef HAVE_MKL
 #include MKL_H
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fft.h"
#include "prefs.h"
#include "psf.h"
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

/*********************************** addback *********************************/
/*
Add a background to the image.
*/
void	addback(simstruct *sim)

  {
   PIXTYPE	*pix, back;
   int		i, nbpix;

  nbpix = sim->imasize[0]*sim->imasize[1];
  back = sim->pixscale[0]*sim->pixscale[1]
	*DEXP(-0.4*(sim->magback-sim->magzero2));
  pix = sim->image;
  for (i=nbpix; i--;)
    *(pix++) += back;

  return;
  }


/********************************** addnoise *********************************/
/*
Add photon and read-out noise to image. 
*/
void	addnoise(simstruct *sim)

  {
#ifndef USE_THREADS
   PIXTYPE	*pix;
   double	ron;
   int		i, nbpix;
#endif

  init_random(0);
#ifdef USE_THREADS
  pthread_addnoise(sim);
#else
  nbpix = sim->imasize[0]*sim->imasize[1];
  ron = sim->ron;
  pix = sim->image;
  for (i=0; i<nbpix; i++, pix++)
    *pix = (PIXTYPE)(random_poisson((double)*pix, 0)+random_gauss(ron, 0));
#endif

  return;
  }


/********************************* cutborders ********************************/
/*
Remove "optical" margins from image, for subsequent processing.
*/
void	cutborders(simstruct *sim)

  {
   PIXTYPE	*pix, *pix2;
   int		m,y, width,fwidth,height, marginstep, nmscan2;

/* No need to remove margins if they do not exist! */
  if (!sim->margin[0] && !sim->margin[1])
    return;
  width = sim->imasize[0]/sim->mscan[0];
  fwidth = sim->fimasize[0]/sim->mscan[0];
  height = sim->imasize[1]/sim->mscan[1];
  marginstep = fwidth*2*sim->margin[1]/sim->mscan[1];
  nmscan2 = sim->mscan[0]*sim->mscan[1];
/* Go on, on an SI-per-SI, line-per-line basis */
  pix = sim->image;
  pix2 = pix+(fwidth*sim->margin[1])/sim->mscan[1]+sim->margin[0]/sim->mscan[0];
  if (sim->margin[1])
    for (m=nmscan2; m--; pix2+=marginstep)
      for (y=height; y--; pix+=width, pix2+=fwidth)
        memcpy(pix, pix2, width*sizeof(PIXTYPE));
  else
/*-- If an overlap is possible, use a secure memmmove instead of memcpy */
    for (m=nmscan2; m--; pix2+=marginstep)
      for (y=height; y--; pix+=width, pix2+=fwidth)
        memmove(pix, pix2, width*sizeof(PIXTYPE));

/* Don't spoil memory */
  QREALLOC(sim->image, PIXTYPE, width*height*nmscan2);

  return;
  }


/*********************************** etoadu **********************************/
/*
Convert the image from photo-electron units to ADU
*/
void	etoadu(simstruct *sim)

  {
   PIXTYPE	*pix, gain, satlev;
   int		i, nbpix;

  
  nbpix = sim->imasize[0]*sim->imasize[1];
  gain = sim->gain;
  satlev = sim->satlev;
  pix = sim->image;
  for (i=nbpix; i--; pix++)
    {
    *pix /= gain;
    if (*pix>satlev)
      *pix = satlev;
    }

  return;
  }


/********************************** microscan ********************************/
/*
Built a microscanned image from a stack of sub-images.
*/
void	microscan(simstruct *sim)

  {
   PIXTYPE	*bmp, *bmpt, *pix;
   int		x,y, xm,ym, nmscanx,nmscany, width,height, hstep;

  if (sim->nmscan==1)
    return;

  QMALLOC(bmp, PIXTYPE, sim->imasize[0]*sim->imasize[1]);
  nmscanx = sim->mscan[0];
  nmscany = sim->mscan[1];
  width = sim->imasize[0]/nmscanx;
  height = sim->imasize[1]/nmscany;
  hstep = sim->imasize[0]*(nmscanx-1);
  pix = sim->image;
  for (ym=0; ym<nmscany; ym++)
    for (xm=0; xm<nmscanx; xm++)
      {
      bmpt = bmp+sim->imasize[0]*ym+xm;
      for (y=0; y<height; y++, bmpt += hstep)
        for (x=0; x<width; x++, bmpt += nmscanx)
          *bmpt = *(pix++);
      }


  free(sim->image);
  sim->image = bmp;    

  return;
  }

/*********************************** saturate ********************************/
/*
Add (CCD) saturation features to an image.
*/
void	saturate(simstruct *sim)

  {
   PIXTYPE	*pix, *pixy, *pixx;
   float	sum, overload, thresh;
   int		i,x,y, width,height, nbpix, nmscan2;

  
  width = sim->imasize[0]/sim->mscan[0];
  height = sim->imasize[1]/sim->mscan[1];
  nbpix = width*height;
  nmscan2 = sim->mscan[0]*sim->mscan[1];
  thresh = sim->wellcap;
  pix = sim->image;
  for (i=0; i<nmscan2; i++, pix+=nbpix)
    {
    pixx = pix;
    for (x=0; x<width; x++, pixx++)
      {
      pixy = pixx;
      sum = 0.0;
/*---- First pass: bottom to top */
      for (y=0; y<height; y++, pixy+=width)
        {
        if ((overload = *pixy - thresh)>0.0)
          {
          overload /= 2.0;
          *pixy -= (PIXTYPE)overload;
          sum += overload;
          }
        else if (sum>0.0)
          {
          if (-overload>sum)
            overload = -sum;
          *pixy -= (PIXTYPE)overload;
          sum += overload;
          }
        }
/*---- Second pass: top to bottom */
      pixy -= width;
      sum = 0.0;
      for (y=0; y<height; y++, pixy -= width)
        {
        if ((overload = *pixy - thresh)>0.0)
          {
          *pixy -= (PIXTYPE)overload;
          sum += overload;
          }
        else if (sum>0.0)
          {
          if (-overload>sum)
            overload = -sum;
          *pixy -= (PIXTYPE)overload;
          sum += overload;
          }
        }
      }
    }

  return;
  }


/********************************** addaureole *******************************/
/*
Add diffusion aureoles to images (heavy multi-convolution job).
*/
void	addaureole(simstruct *sim)

  {
   PIXTYPE	*pix,*mpix,*mima, *sipix,*sipix2,*simax,*simay,*sima,
		*chbuffer,*nhbuffer, *hbuffer1,*hbuffer2, *faureole;
   char		str[MAXCHAR];
   int		m,y, nmscan2, nsixm1,nsiym1, six,siy, nbpix, nbuffer,
		marginx,marginy, mwidth,mheight, fwidth,fheight,
		sistepx,sistepy,siwidth,siwidthu,siheight,siheightu,
		istep, nstep;

  if (!sim->margin[0] && !sim->margin[1])
    return;
  makeaureole(sim);
  mwidth = sim->aursize[0];
  mheight = sim->aursize[1];
  nmscan2 = sim->mscan[0]*sim->mscan[1];
  marginx = sim->margin[0]/sim->mscan[0];
  marginy = sim->margin[1]/sim->mscan[1];
  fwidth = sim->fimasize[0]/sim->mscan[0];
  fheight = sim->fimasize[1]/sim->mscan[1];
  nsixm1 = (sim->imasize[0]/sim->mscan[0])/(mwidth-2*marginx+1);
  nsiym1 = (sim->imasize[1]/sim->mscan[1])/(mheight-2*marginy+1);
  nbpix = fwidth*fheight;
  nbuffer = fwidth*marginy;
  sistepx = mwidth-2*marginx;
  sistepy = mheight-2*marginy;
  sima = sim->image;
  hbuffer1 = hbuffer2 = NULL;		/* to avoid gcc -Wall warnings */
  if (nsiym1>0)
    QMALLOC(hbuffer1, PIXTYPE, nbuffer);
  if (nsiym1>1)
    QMALLOC(hbuffer2, PIXTYPE, nbuffer);
  nstep = (nsixm1+1)*(nsiym1+1)*nmscan2;
  istep = 0;

  fft_init(prefs.nthreads);
  QFFTWF_CALLOC(mima, PIXTYPE, mwidth*mheight);
  faureole = fft_rtf(sim->aureole, mwidth, mheight);
  for (m=0; m<nmscan2; m++, sima+=nbpix)
    {
    siheight = mheight;
    siheightu = mheight-marginy;

/*-- Vertical loop over sub-images */
    simay = sima;
    for (siy=0; siy<=nsiym1; siy++, simay+=sistepy*fwidth)
      {
/*---- Save full-length lines of pixels of the 1/2 overlap before convolving */
      if (siy==nsiym1)
        siheight = siheightu = fheight - siy*sistepy;
      else
        {
        nhbuffer = (siy%2)?hbuffer2:hbuffer1;
        memcpy(nhbuffer, simay+sistepy*fwidth, nbuffer*sizeof(PIXTYPE));
        }
      if (siy)
        chbuffer = (siy%2)?hbuffer1:hbuffer2;
      else
        chbuffer = simay;
      siwidth = mwidth;
      siwidthu = mwidth-marginx;

/*---- Horizontal loop over sub-images */
      simax = simay;
      for (six=0; six<=nsixm1; six++, simax+=sistepx)
        {
        if (nstep>1)
          {
          sprintf(str, "Adding diffuse component... (%2.0f%% done)",
		(float)100.0*istep++/nstep);
          NFPRINTF(OUTPUT, str);
          }

        if (six==nsixm1)
          siwidth = siwidthu = fwidth - six*sistepx;

/*------ Copy the 1st sub-image of the row to the "Fourier-drome" */
        pix = chbuffer+six*sistepx;
        sipix2 = simax+nbuffer;
        mpix = mima;
        for (y=0; y<siheight; y++, pix+=fwidth, mpix+=mwidth)
          {
          if (y==marginy && siy)
            pix = sipix2;
          if (six)
            memcpy(mpix+marginx, pix+marginx,(siwidth-marginx)*sizeof(PIXTYPE));
          else
            memcpy(mpix, pix, siwidth*sizeof(PIXTYPE));
          }

/*------ Now the important call: fast convolution */
        fft_conv(mima, faureole, mwidth, mheight);

/*------ Copy back the relevant part from the Fourier-drome to the sub-image */
/*------ ...while preparing the next Fourier: "prise de tete" !! */
        pix = chbuffer+(six+1)*sistepx;
        sipix = simax;
        sipix2 = simax+nbuffer+sistepx;
        mpix = mima;
        for (y=0; y<siheight; y++, sipix+=fwidth, mpix+=mwidth)
          {
          if (!six && y<siheightu && (!siy || y>=marginy))
            memcpy(sipix, mpix, marginx*sizeof(PIXTYPE));
          if (six<nsixm1)
            {
            if (y==marginy && siy)
              pix = sipix2;
            memcpy(mpix, pix, marginx*sizeof(PIXTYPE));
            pix += fwidth;
            }
          if (y<siheightu && (!siy || y>=marginy))
            memcpy(sipix+marginx, mpix+marginx,
		(siwidthu-marginx)*sizeof(PIXTYPE));
          }
        }
      }
    }

  if (nsiym1>0)
    free(hbuffer1);
  if (nsiym1>1)
    free(hbuffer2);
  QFFTWF_FREE(mima);
  QFFTWF_FREE(faureole);
  QFFTWF_FREE(sim->aureole);
  fft_end(prefs.nthreads);

  return;
  }


#ifdef USE_THREADS

/****** pthread_addnoiseline ******************************************************
PROTO	void *pthread_addnoiseline(void *arg)
PURPOSE	Noise generation thread.
INPUT	Pointer to the thread number.
OUTPUT	NULL void pointer.
NOTES	Relies on global variables.
AUTHOR	E. Bertin (IAP)
VERSION	10/03/2016
 ***/
static void	*pthread_addnoiseline(void *arg)
  {
   double	ron;
   PIXTYPE	*pix, *pixt;
   int		i, line, npix, nlines, proc;
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
    pixt = pix = pthread_sim->image + line*npix;
#ifdef HAVE_MKL
    lambdabuft = lambdabuf;
    for (i=npix; i--;)
      *(lambdabuft++) = *(pixt++);
    viRngPoissonV(VSL_RNG_METHOD_POISSONV_POISNORM, stream, npix,
		poissonbuf, lambdabuf);
    vsRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, npix,
		gaussbuf, 0.0, ron);
    pixt = pix;
    poissonbuft = poissonbuf;
    gaussbuft = gaussbuf;
    for (i=npix; i--;)
      *(pixt++) = (PIXTYPE)*(poissonbuft++) + (PIXTYPE)*(gaussbuft++);
#else
    for (i=npix; i--; pixt++)
      *pixt = (PIXTYPE)(random_poisson((double)*pixt, proc)
		+random_gauss(ron, proc));
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

