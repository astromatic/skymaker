/*
				alterimage.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Routines for alterating the image as a whole.
*
*	Last modify:	15/01/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <fftw3.h>

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
Add photon photo and read-out noise to image. 
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
       /* memmove(pix, pix2, width*sizeof(PIXTYPE)); */
        memcpy(pix, pix2, width*sizeof(PIXTYPE));

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
  QFFTWCALLOC(mima, PIXTYPE, mwidth*mheight);
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
  QFFTWFREE(mima);
  QFFTWFREE(faureole);
  QFFTWFREE(sim->aureole);
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
VERSION	20/08/2006
 ***/
static void	*pthread_addnoiseline(void *arg)
  {
   double	ron;
   PIXTYPE	*pix;
   int		i, line, npix, nlines, proc;

  proc = *((int *)arg);
  ron = pthread_sim->ron;
  npix = pthread_sim->imasize[0];
  nlines = pthread_sim->imasize[1];
/* Exit if the end of image has been reached */
  QPTHREAD_MUTEX_LOCK(&noisemutex);  
  while (pthread_line < nlines)
    {
    line = pthread_line++;
    QPTHREAD_MUTEX_UNLOCK(&noisemutex);  
    pix = pthread_sim->image + line*npix;
    for (i=npix; i--; pix++)
      *pix = (PIXTYPE)(random_poisson((double)*pix, proc)
		+random_gauss(ron, proc));
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
VERSION	17/08/2006
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
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_addnoiseline, &proc[p]);
    }
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
/* Clean up multi-threading stuff */
  QPTHREAD_MUTEX_DESTROY(&noisemutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(proc);
  free(thread);

  return;
  }

#endif

