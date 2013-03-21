/*
*				galaxies.c
*
* Generate galaxy models.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		21/03/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include FFTW_H

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
extern pthread_mutex_t	fftmutex;
#endif

static double	gamm(double xx);

/****** make_galaxy *******************************************************
PROTO	void make_galaxy(simstruct *sim, objstruct *obj)
PURPOSE	Generate a galaxy at a definite position in the image.
INPUT	Pointer to the sim structure,
	pointer to the object structure.
OUTPUT	-.
NOTES	Writes to an allocated image buffer, not directly to the image to
	allow multithreading.
AUTHOR	E. Bertin (IAP)
VERSION	21/03/2013
 ***/
void	make_galaxy(simstruct *sim, objstruct *obj)

  {
   char		str[MAXCHAR];
   PIXTYPE	*bsub, *bsubt, *dsub, *sub, *subt, *psfdft;
   double	dpos[2],
		osamp, dratio, bsize,size,
		bflux,flux,flux2, dx,dy, beq, dscale,expo, n, bn, ampfac, dval;
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
	obj->pos[0]+1.0, obj->pos[1]+1.0);
      warning(str, "skipped");
      return;
      }
    n = obj->bulge_sersicn;
    bn = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
	+ 131.0/(1148175*n*n*n);        /* Ciotti & Bertin 1999 */
    ampfac = pow(bn, 2*n) / (PI * gamm(2*n+1));
    bsize = log(obj->bulge_ratio*obj->flux*ampfac
		/(obj->bulge_ar*beq*sim->minquant))/bn;
    size = bsize = bsize>0.0 ? beq*pow(bsize, n) : 0.0;
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
	obj->pos[0]+1.0,obj->pos[1]+1.0);
      warning(str, "skipped");
      return;
      }
    if (obj->disk_ar < SMALL)
      {
      NFPRINTF(OUTPUT, "");
      sprintf(str, "Disk A/R too small for galaxy at (%.1f,%.1f): ",
	obj->pos[0]+1.0, obj->pos[1]+1.0);
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
	subwidth, subheight, obj->pos[0]+1.0,obj->pos[1]+1.0,obj->mag);

  NFPRINTF(OUTPUT, str);
*/
  nsub = subwidth*subheight;
  memnsub = (((subwidth>>1) + 1)<< 1) * subheight;/* Provide margin for FFTW */

/* Compute (or retrieve) PSF DFT at this position */
  psfdft = interp_dft(sim, suborder, obj->pos, dpos);

  sub = NULL;			/* To avoid gcc -Wall warnings */
  flux = 0.0;			/* idem */
/* Bulge component */
  if (obj->bulge_ratio)
    {
    QFFTWF_MALLOC(bsub, PIXTYPE, memnsub);
    memset(bsub, 0, memnsub*sizeof(PIXTYPE));
    flux = bflux = make_sersic(bsub, subwidth, subheight, osamp*beq,
	obj->bulge_ar,obj->bulge_posang, n) / obj->bulge_ratio;
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
    QFFTWF_MALLOC(dsub, PIXTYPE, memnsub);
    memset(dsub, 0, memnsub*sizeof(PIXTYPE));
    flux = make_sersic(dsub, subwidth, subheight,  osamp*dscale*1.67835,
		obj->disk_ar, obj->disk_posang, 1.0) / dratio;
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
  fft_conv(sub, psfdft, subwidth,subheight);
  if (sim->npsf>1)
    free(psfdft);
  dx = (obj->pos[0] + sim->margin[0] - dpos[0])/sim->mscan[0];
  dy = (obj->pos[1] + sim->margin[1] - dpos[1])/sim->mscan[1];
  dx -= (double)(obj->subpos[0] = (int)(dx+0.49999));
  dy -= (double)(obj->subpos[1] = (int)(dy+0.49999));
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
  flux = flux2 = 0.0;
  for (i=nsub2,subt=obj->subimage; i--;)
    {
    dval = (double)*(subt++);
    flux += dval;
    flux2 += dval*dval;
    }

  obj->subfactor = obj->flux/flux;
  obj->noiseqarea = flux*flux / flux2;

  QFFTWF_FREE(bsub);
  QFFTWF_FREE(dsub);

  return;
  }


/****** make sersic **********************************************************
PROTO	double make_sersic(PIXTYPE *pix, int width, int height, double reff,
		double aspect, double posang, double n)
PURPOSE	Rasterize an unnormalized 2D Sersic (exp(r**(-1/n))) model.
INPUT	pointer to the raster,
	raster width,
	raster height,
	model effective radius,
	model aspect ratio,
	model position angle,
	Sersic index.
OUTPUT	Total flux.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/07/2011
 ***/
double	make_sersic(PIXTYPE *pix, int width, int height, double reff,
		double aspect, double posang, double n)

  {
   double	ddx1[36],ddx2[36], cpos,spos, a2,b2, bn, sum;
   float	k, invn,hinvn,krspinvn,ekrspinvn, krpinvn,dkrpinvn,
		r,r2,dr, rs,rs2,
		p1,p2,c0,c2,c3, x1c,x2c, x1,x10,x2, cxx,cyy,cxy,cxydy,cyydy2,
		dx1c,dx2c,a11,a12,a21,a22, ang,angstep, dca,dsa, selem, val;
   PIXTYPE	*pixt,*pixt2;
   int		a,i, ix1,ix2, nx2,npix2, nang;

/* No point source */
  if (reff == 0.0)
    return 0.0;

/* Preliminaries: compute the ellipse parameters */
  a2 = 1.0/(reff*reff);
  b2 = a2/(aspect*aspect);
  cpos = cos(posang*DEG);
  spos = sin(posang*DEG);
  cxx = cpos*cpos*a2 + spos*spos*b2;
  cyy = spos*spos*a2 + cpos*cpos*b2;
  cxy = 2*cpos*spos*(a2-b2);

/* Compute sharp/smooth transition radius */
  rs = SERSIC_SMOOTHR/(reff*aspect);
  if (rs<=0)
    rs = 1.0;
  rs2 = rs*rs;
  if (n==1.0)
    bn = 1.67835;
  else if (n==4.0)
    bn = 7.66924944;
  else
    bn = 2.0*n - 1.0/3.0 + 4.0/(405.0*n) + 46.0/(25515.0*n*n)
		+ 131.0/(1148175*n*n*n);	/* Ciotti & Bertin 1999 */
  invn = 1.0/n;
  hinvn = 0.5/n;
  k = -bn;
/* Compute central polynomial terms */
  krspinvn = (n==1.0)? k*rs : k*expf(logf(rs)*invn);
  ekrspinvn = expf(krspinvn);
  p2 = krspinvn*invn*invn;
  p1 = krspinvn*p2;
  c0 = (1+(1.0/6.0)*(p1+(1.0-5.0*n)*p2))*ekrspinvn;
  c2 = (-1.0/2.0)*(p1+(1.0-3.0*n)*p2)/rs2*ekrspinvn;
  c3 = (1.0/3.0)*(p1+(1.0-2.0*n)*p2)/(rs2*rs)*ekrspinvn;
  x1c = (float)(width/2);
  x2c = (float)(height/2);
  nx2 = height/2 + 1;

/* Compute the smooth part of the profile */
  x10 = -x1c;
  x2 = -x2c;
  pixt = pix;
  if (n==1.0)
    for (ix2=nx2; ix2--; x2+=1.0)
      {
      x1 = x10;
      cyydy2 = cyy*x2*x2;
      cxydy = cxy*x2;
      for (ix1=width; ix1--; x1+=1.0)
        *(pixt++) = (r2=cyydy2 + (cxx*x1+cxydy)*x1) < rs2?
		c0+r2*(c2+c3*sqrtf(r2)) : expf(-1.67835f*sqrtf(r2));
      }
  else
    for (ix2=nx2; ix2--; x2+=1.0)
      {
      x1 = x10;
      cyydy2 = cyy*x2*x2;
      cxydy = cxy*x2;
      for (ix1=width; ix1--; x1+=1.0)
        *(pixt++) = (r2=cyydy2 + (cxx*x1+cxydy)*x1)<rs2?
		c0+r2*(c2+c3*sqrtf(r2)) : expf(k*expf(logf(r2)*hinvn));
      }
/* Copy the symmetric part */
  if ((npix2=(height-nx2)*width) > 0)
    {
    pixt2 = pixt - width - 1;
    if (!(width&1))
      {
      *(pixt++) = 0.0;
      npix2--;
      }
    for (i=npix2; i--;)
      *(pixt++) = *(pixt2--);
    }

/* Compute the sharp part of the profile */
  dx1c = x1c + 0.4999999;
  dx2c = x2c + 0.4999999;
  a11 = reff*cpos;
  a12 = -reff*aspect*spos;
  a21 = reff*spos;
  a22 = reff*aspect*cpos;
  nang = 72 / 2;	/* 72 angles; only half of them are computed*/
  angstep = PI/nang;
  ang = 0.0;
  for (a=0; a<nang; a++)
    {
#ifdef HAVE_SINCOSF
    sincosf(ang, &dsa, &dca);
#else
    dsa = sinf(ang);
    dca = cosf(ang);
#endif
    ddx1[a] = a11*dsa+a12*dca;
    ddx2[a] = a21*dsa+a22*dca;
    ang += angstep;
    }
  r = DEXPF(-5.0);
  dr = DEXPF(0.05);
  selem = 0.5*angstep*(dr - 1.0/dr)*reff*reff*aspect;
  krpinvn = k*DEXPF(-5.0*invn);
  dkrpinvn = DEXPF(0.05*invn);
  for (; r<rs; r *= dr)
    {
    r2 = r*r;
    val = (expf(krpinvn) - (c0 + r2*(c2+c3*r)))*r2*selem;
    for (a=0; a<nang; a++)
      {
      ix1 = (int)(dx1c + r*ddx1[a]);
      ix2 = (int)(dx2c + r*ddx2[a]);
      if (ix1>=0 && ix1<width && ix2>=0 && ix2<height)
        pix[ix2*width+ix1] += val;
      ix1 = (int)(dx1c - r*ddx1[a]);
      ix2 = (int)(dx2c - r*ddx2[a]);
      if (ix1>=0 && ix1<width && ix2>=0 && ix2<height)
        pix[ix2*width+ix1] += val;
      }
    krpinvn *= dkrpinvn;
    }

/* Compute integral flux */
  pixt = pix;
  sum = 0.0;
  for (i=width*height; i--;)
    sum += *(pixt++);

  return sum;
  }


/****** trunc_prof **********************************************************
PROTO	PIXTYPE	trunc_prof(PIXTYPE *pix, double xcenter, double ycenter,
		int width, int height)
PURPOSE	Truncate a monotonously decreasing 2D model inside a raster using
	thresholding.
INPUT	pointer to the raster,
	model center x coordinate,
	model center y coordinate,
	raster width,
	raster height.
OUTPUT	Computed threshold.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	01/03/2011
 ***/
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


/****i* gamm ***************************************************************
PROTO   double gamm(double xx)
PURPOSE Returns the Gamma function (based on algorithm described in
	Numerical Recipes in C, chap 6.1).
INPUT   A double.
OUTPUT  Gamma function.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 21/12/2010
*/
static double	gamm(double xx)

  {
   double		x,tmp,ser;
   static double	cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
   int			j;

  tmp=(x=xx-1.0)+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<6;j++)
    ser += cof[j]/(x+=1.0);

  return 2.50662827465*ser*exp(-tmp);
  }

