/*
				galaxies.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Routines for creating galaxy fields.
*
*	Last modify:	19/03/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <fftw3.h>

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "galaxies.h"
#include "image.h"
#include "list.h"
#include "prefs.h"
#include "psf.h"
#include "simul.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

#ifdef USE_THREADS
pthread_mutex_t		dftmutex;
extern pthread_mutex_t	fftmutex;
#endif

/****** make_galaxy *******************************************************
PROTO	void make_galaxy(simstruct *sim, objstruct *obj)
PURPOSE	Generate a galaxy at a definite position in the image.
INPUT	Pointer to the sim structure,
	pointer to the object structure.
OUTPUT	-.
NOTES	Writes to an allocated image buffer, not directly to the image to
	allow multithreading.
AUTHOR	E. Bertin (IAP)
VERSION	19/03/2008
 ***/
void	make_galaxy(simstruct *sim, objstruct *obj)

  {
   char		str[MAXCHAR];
   PIXTYPE	*bsub, *bsubt, *dsub, *sub, *subt;
   double	osamp, dratio, bsize,size,
		bflux,flux, dx,dy, beq, dscale,expo;
   int		i, subwidth,subheight, suborder,
		nsub,nsub2,nsubo,memnsub;

  osamp = sim->psfoversamp;
  beq = obj->bulge_req/sim->pixscale[0];
  dscale = obj->disk_scale/sim->pixscale[0];

/* Convert magnitude to linear units */
  if (!obj->flux)
    {
    expo = 0.4*(sim->magzero2-obj->mag);
    if (expo>-30.0)
      obj->flux = DEXP(expo);
    if (!obj->flux)
      return;
    }
/* No need to quantize below 1 sigma at the object scale: */
/* Set mask size limits from the bulge */
  if (obj->bulge_ratio)
    {
    if (beq < SMALL)
      {
      NFPRINTF(OUTPUT, "");
      sprintf(str, "Bulge radius too small for galaxy at (%.1f,%.1f): ",
	obj->x+1.0, obj->y+1.0);
      warning(str, "skipped");
      return;
      }
    bsize = log(obj->bulge_ratio*obj->flux*94.480
		/(obj->bulge_ar*beq*sim->minquant))/7.6693;
    size = bsize = bsize>0.0 ? beq*pow(bsize, 4.0) : 0.0;
    }
  else
    size = bsize = 0.0;

/* Set mask size limits from the disk */
  if ((dratio=(1.0-obj->bulge_ratio))>0.001) 
    {
    if (dscale < SMALL)
      {
      NFPRINTF(OUTPUT, "");
      sprintf(str, "Disk scale too small for galaxy at (%.1f,%.1f): ",
	obj->x+1.0,obj->y+1.0);
      warning(str, "skipped");
      return;
      }
    if (obj->disk_ar < SMALL)
      {
      NFPRINTF(OUTPUT, "");
      sprintf(str, "Disk A/R too small for galaxy at (%.1f,%.1f): ",
	obj->x+1.0, obj->y+1.0);
      warning(str, "skipped");
      return;
      }

/*-- Estimate the maximum extent of the disk */
    size = dscale*log(dratio*obj->flux
		/(2*PI*obj->disk_ar*dscale*sim->minquant));
    if (size<0.0)
      size = 0.0;
/*-- Don't go beyond van der Kruit cut-off */
/*
    if (size > (size2=dscale*(VDKCUTRAD+3.0)))
      size = size2;
*/
    }
  if (size<bsize)
    size = bsize;
  i = 2*(int)(osamp*sqrt(size*size
	+16.0*sim->psfarea/(sim->pixscale[0]*sim->pixscale[1])));
  for (suborder=1; i>>=1; suborder++);
  if (suborder>=PSF_NORDER)
    {
    osamp /= (double)(1<<(suborder-PSF_NORDER));
    suborder = PSF_NORDER;
    }
  subheight = subwidth = 1<<suborder;
/*
  sprintf(str,
	"Adding %dx%d galaxy at position (%.1f,%.1f) with magnitude %.2f",
	subwidth, subheight, obj->x+1.0,obj->y+1.0,obj->mag);

  NFPRINTF(OUTPUT, str);
*/
  nsub = subwidth*subheight;
  memnsub = (((subwidth>>1) + 1)<< 1) * subheight;/* Provide margin for FFTW */

/* Check that PSF DFTs already exist at that mask size */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&dftmutex);
#endif
  if (!sim->psfdft[suborder])
    makedft(sim, suborder);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&dftmutex);
#endif

  sub = NULL;			/* To avoid gcc -Wall warnings */
  flux = 0.0;			/* idem */
/* Bulge component */
  if (obj->bulge_ratio)
    {
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
    QFFTWMALLOC(bsub, PIXTYPE, memnsub);
#ifdef USE_THREADS
    QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
    memset(bsub, 0, memnsub*sizeof(PIXTYPE));
    flux = bflux = add_devauc(bsub, (double)(subwidth/2),(double)(subheight/2),
		subwidth, subheight,
		1e6, osamp*beq, obj->bulge_ar, obj->bulge_posang)
			/ obj->bulge_ratio;
    sub = bsub;
    }
  else
    {
    bflux = 0.0;		/* To avoid gcc -Wall warnings */
    bsub = NULL;
    }

/* Disk component */
  if (dratio>0.001) 
    {
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
    QFFTWMALLOC(dsub, PIXTYPE, memnsub);
#ifdef USE_THREADS
    QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif
    memset(dsub, 0, memnsub*sizeof(PIXTYPE));
    flux = add_expo(dsub, (double)(subwidth/2),(double)(subheight/2),
		subwidth, subheight,
		1e6, osamp*dscale, obj->disk_ar, obj->disk_posang) / dratio;
    sub = dsub;
    }
  else
    dsub = NULL;

/* Combine and normalize (flux will change later, though) */
  if (bsub && dsub)
    for (i=nsub,subt=sub,bsubt=bsub; i--; subt++)
      *subt = *subt/flux + *(bsubt++)/bflux;
  else
    for (i=nsub,subt=sub; i--;)
      *(subt++) /= flux;

/* Truncate to avoid introducing anistropy in the vignet corners */
  trunc_prof(sub, (double)(subwidth/2),(double)(subheight/2),
		subwidth, subheight);
/* Convolve with the oversampled PSF */
  fft_conv(sub, sim->psfdft[suborder], subwidth,subheight);
  dx = (obj->x + sim->margin[0])/sim->mscan[0];
  dy = (obj->y + sim->margin[1])/sim->mscan[1];
  dx -= (double)(obj->subx = (int)(dx+0.49999));
  dy -= (double)(obj->suby = (int)(dy+0.49999));
/* Truncate again */
  trunc_prof(sub, (double)(subwidth/2),(double)(subheight/2),
		subwidth, subheight);
/* Resample to lower resolution */
  nsubo = obj->subsize[0]*obj->subsize[1];
  obj->subsize[0]  = obj->subsize[1] = (int)(subwidth/osamp+1.0);
  nsub2 = obj->subsize[0]*obj->subsize[1];
  if (!nsubo)
    {
    QMALLOC(obj->subimage, PIXTYPE, nsub2);
    }
  else if (nsub2>nsubo)
    {
    QREALLOC(obj->subimage, PIXTYPE, nsub2);
    }
  resample_image(sub, subwidth, subheight, obj, -dx*osamp, -dy*osamp, osamp);
  flux = 0.0;
  for (i=nsub2,subt=obj->subimage; i--;)
    flux += (double)*(subt++);

  obj->subfactor = obj->flux/flux;

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWFREE(bsub);
  QFFTWFREE(dsub);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  return;
  }


/********************************* add_devauc ********************************/
/*
Add a 2D de Vaucouleurs (exp(r**(-1/4))) profile to a vignet.
*/
double	add_devauc(PIXTYPE *pix, double xcenter, double ycenter,
		int width, int height,
		double flux, double req, double aspect, double posang)

  {
   double	amp, cpos,spos, cxx,cyy,cxy,cxydy,cyydy2,
		a2,b2,dx,dy, sum;
   int		x,y;

  if (req == 0.0)
    return 0.0;

/* Preliminaries: compute the ellipse parameters */
  a2 = 1/(req*req);
  b2 = a2/(aspect*aspect);
  amp = flux*94.480/(req*req);
  cpos = cos(posang*DEG);
  spos = sin(posang*DEG);
  cxx = cpos*cpos*a2 + spos*spos*b2;
  cyy = spos*spos*a2 + cpos*cpos*b2;
  cxy = 2*cpos*spos*(a2-b2);

/* Compute and put the values */
  dy = -ycenter;
  sum = 0.0;
  for (y=height; y--; dy += 1.0)
    {
    cxydy = cxy*dy;
    cyydy2 = cyy*dy*dy;
    dx = -xcenter;
    for (x=width; x--; dx += 1.0)
      sum += (*(pix++) = (PIXTYPE)(amp*exp(-7.6693
	*pow(cxx*dx*dx + cyydy2 + cxydy*dx, 0.125))));
    }

  return sum;
  }


/********************************** add_expo *********************************/
/*
Add a 2D exponential profile to a vignet.
*/
double	add_expo(PIXTYPE *pix, double xcenter, double ycenter,
		int width, int height,
		double flux, double scale, double aspect, double posang)

  {
   double	amp, cpos,spos, cxx,cyy,cxy,cxydy,cyydy2,
		a2,b2,dx,dy, sum;
   int		x,y;

  if (scale == 0.0)
    return 0.0;

/* Preliminaries: compute the ellipse parameters */
  a2 = 1/(scale*scale);
  b2 = a2/(aspect*aspect);
  amp = flux/(2*PI*scale*scale);
  cpos = cos(posang*DEG);
  spos = sin(posang*DEG);
  cxx = cpos*cpos*a2 + spos*spos*b2;
  cyy = spos*spos*a2 + cpos*cpos*b2;
  cxy = 2*cpos*spos*(a2-b2);

/* Compute and put the values */
  dy = -ycenter;
  sum = 0.0;
  for (y=height; y--; dy += 1.0)
    {
    cxydy = cxy*dy;
    cyydy2 = cyy*dy*dy;
    dx = -xcenter;
    for (x=width; x--; dx += 1.0)
      sum += (*(pix++) = (PIXTYPE)(amp*exp(-sqrt(cxx*dx*dx+cyydy2+cxydy*dx))));
    }

  return sum;
  }


/********************************** trunc_prof *******************************/
/*
Truncate a 2D profile in a vignet.
*/
PIXTYPE	trunc_prof(PIXTYPE *pix, double xcenter, double ycenter,
		int width, int height)

  {
   double	dx,dy, r, rmin,rmin2, rmax,rmax2;
   PIXTYPE	*pixt, thresh, val;
   int		x,y, n;


 /* Find the shortest distance to a vignet border */
  rmax = width - xcenter;
  if (rmax > (r=xcenter+1.0))
    rmax = r;
  if (rmax > (r=height-ycenter))
    rmax = r;
  if (rmax > (r=ycenter+1.0))
    rmax = r;
  rmax -= 0.99;
  if (rmax<1.0)
    rmax = 1.0;
  rmax2 = rmax*rmax;
  rmin = rmax - 1.0;
  rmin2 = rmin*rmin;

/* Find best threshold (the max around the circle with radius rmax */
  dy = -ycenter;
  thresh = -BIG;
  pixt = pix;
  for (y=height; y--; dy += 1.0)
    {
    dx = -xcenter;
    for (x=width; x--; dx += 1.0)
      if ((val=*(pixt++))>thresh && (r=dx*dx+dy*dy)>rmin2 && r<rmax2)
        thresh = val;
    }

/* Threshold! */
  pixt = pix;
  for (n=width*height; n--; pixt++)
    if (*pixt<thresh)
      *pixt = 0.0;

  return thresh;
  }


