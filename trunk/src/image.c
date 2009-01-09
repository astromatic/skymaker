 /*
 				image.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Function related to image manipulations.
*
*	Last modify:	07/04/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#define _GNU_SOURCE
#include <math.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "image.h"
#include "list.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

#ifdef USE_THREADS
pthread_mutex_t addimagemutex;
#endif

static void	make_kernel(double pos, double *kernel, interpenum interptype);

int		interp_kernwidth[5]={1,2,4,6,8};

/****** resample_image *******************************************************
PROTO	int resample_image(PIXTYPE *pix1, int w1, int h1, objstruct *obj,
			double dx, double dy, double step2)
PURPOSE Scale and shift a small image through interpolation. Image parts which
	lie outside boundaries are ignored.
INPUT   Pointer to input raster,
	input width,
	input height,
	pointer to output object (pointing to the output raster). 
        shift in x,
	shift in y,
	step in pixels.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 23/09/2006
 ***/
/*
Scale and shift a small image through sinc interpolation.
Image parts which lie outside boundaries are ignored.
*/
int	resample_image(PIXTYPE *pix1, int w1, int h1, objstruct *obj,
			double dx, double dy, double step2)
  {
   interpenum	interp_type;
   PIXTYPE	*pixin,*pixin0, *pixout,*pixout0;
   double	kernel[INTERP_MAXKERNELWIDTH], *kernelt,
		*dpixin,*dpixin0, *dpixout,*dpixout0, *maskt,
		xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, val;
   int		i,j,k,n,t, *startt, *nmaskt,
		ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
		ix,iy, ix1,iy1, w2,h2, interpw;

  interp_type = INTERP_TYPE;
  interpw = interp_kernwidth[interp_type];
  w2 = obj->subsize[0];
  h2 = obj->subsize[1];
  xc1 = (double)(w1/2);		/* Im1 center x-coord*/
  xc2 = (double)(w2/2);		/* Im2 center x-coord*/
  xs1 = xc1 + dx - xc2*step2;		/* Im1 start x-coord */
  if ((int)xs1 >= w1)
    return RETURN_ERROR;
  ixs2 = 0;			/* Int part of Im2 start x-coord */
  if (xs1<0.0)
    {
    dix2 = (int)(1-xs1/step2);
/*-- Simply leave here if the images do not overlap in x */
    if (dix2 >= w2)
      return RETURN_ERROR;
    ixs2 += dix2;
    xs1 += dix2*step2;
    }
  nx2 = (int)((w1-xs1)/step2);/* nb of interpolated Im2 pixels along x */
  if (nx2>(ix2=w2-ixs2))
    nx2 = ix2;
  if (!nx2)
    return RETURN_ERROR;
  yc1 = (double)(h1/2);		/* Im1 center y-coord */
  yc2 = (double)(h2/2);		/* Im2 center y-coord */
  ys1 = yc1 + dy - yc2*step2;		/* Im1 start y-coord */
  if ((int)ys1 >= h1)
    return RETURN_ERROR;
  iys2 = 0;			/* Int part of Im2 start y-coord */
  if (ys1<0.0)
    {
    diy2 = (int)(1-ys1/step2);
/*-- Simply leave here if the images do not overlap in y */
    if (diy2 >= h2)
      return RETURN_ERROR;
    iys2 += diy2;
    ys1 += diy2*step2;
    }
  ny2 = (int)((h1-ys1)/step2);/* nb of interpolated Im2 pixels along y */
  if (ny2>(iy2=h2-iys2))
    ny2 = iy2;
  if (!ny2)
    return RETURN_ERROR;

/* Set the yrange for the x-resampling with some margin for interpolation */
  iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
  hmh = hmw = interpw/2 - 1;	/* Interpolant start */
  if (iys1a<0 || ((iys1a -= hmh)< 0))
    iys1a = 0;
  ny1 = (int)(ys1+ny2*step2)+interpw-hmh;	/* Interpolated Im1 y size */
  if (ny1>h1)					/* with margin */
    ny1 = h1;
/* Express everything relative to the effective Im1 start (with margin) */
  ny1 -= iys1a;
  ys1 -= (double)iys1a;

/* Allocate interpolant stuff for the x direction */
  if (!obj->buf1size)
    {
    obj->buf1size = nx2;
    QMALLOC(obj->maskbuf, double, nx2*interpw);	/* Interpolation masks */
    QMALLOC(obj->nmaskbuf, int, nx2);		/* Interpolation mask sizes */
    QMALLOC(obj->startbuf, int, nx2);	/* Int part of Im1 conv starts */
    }
  else if (obj->buf1size < nx2)
    {
    obj->buf1size = nx2;
    QREALLOC(obj->maskbuf, double, nx2*interpw);
    QREALLOC(obj->nmaskbuf, int, nx2);
    QREALLOC(obj->startbuf, int, nx2);
    }

/* Compute the local interpolant and data starting points in x */
  x1 = xs1;
  maskt = obj->maskbuf;
  nmaskt = obj->nmaskbuf;
  startt = obj->startbuf;
  for (j=nx2; j--; x1+=step2)
    {
    ix = (ix1=(int)x1) - hmw;
    make_kernel(x1-ix1, kernel, interp_type);
    kernelt = kernel;
    if (ix < 0)
      {
      n = interpw+ix;
      kernelt -= ix;
      ix = 0;
      }
    else
      n = interpw;
    if (n>(t=w1-ix))
      n=t;
    *(startt++) = ix;
    *(nmaskt++) = n;
    for (i=n; i--;)
      *(maskt++) = *(kernelt++);
    }

  if (!obj->buf2size)
    {
    obj->buf2size = nx2*ny1;
    QCALLOC(obj->buf2, double, obj->buf2size);	/* Intermediary frame-buffer */
    }
  else if (obj->buf2size < nx2*ny1)
    {
    obj->buf2size = nx2*ny1;
    free(obj->buf2);
    QCALLOC(obj->buf2, double, obj->buf2size);	/* Intermediary frame-buffer */
    }

/* Make the interpolation in x (this includes transposition) */
  pixin0 = pix1+iys1a*w1;
  dpixout0 = obj->buf2;
  for (k=ny1; k--; pixin0+=w1, dpixout0++)
    {
    maskt = obj->maskbuf;
    nmaskt = obj->nmaskbuf;
    startt = obj->startbuf;
    dpixout = dpixout0;
    for (j=nx2; j--; dpixout+=ny1)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)*(double)*(pixin++);
      *dpixout = val;
      }
    }

/* Reallocate interpolant stuff for the y direction */
  if (obj->buf1size < ny2)
    {
    obj->buf1size = ny2;
    QREALLOC(obj->maskbuf, double, ny2*interpw);
    QREALLOC(obj->nmaskbuf, int, ny2);
    QREALLOC(obj->startbuf, int, ny2);
    }

/* Compute the local interpolant and data starting points in y */
  y1 = ys1;
  maskt = obj->maskbuf;
  nmaskt = obj->nmaskbuf;
  startt = obj->startbuf;
  for (j=ny2; j--; y1+=step2)
    {
    iy = (iy1=(int)y1) - hmh;
    make_kernel(y1-iy1, kernel, interp_type);
    kernelt = kernel;
    if (iy < 0)
      {
      n = interpw+iy;
      kernelt -= iy;
      iy = 0;
      }
    else
      n = interpw;
    if (n>(t=ny1-iy))
      n=t;
    *(startt++) = iy;
    *(nmaskt++) = n;
    for (i=n; i--;)
      *(maskt++) = *(kernelt++);
    }

/* Initialize destination buffer to zero */
  memset(obj->subimage, 0, (size_t)(w2*h2)*sizeof(PIXTYPE));

/* Make the interpolation in y  and transpose once again */
  dpixin0 = obj->buf2;
  pixout0 = obj->subimage+ixs2+iys2*w2;
  for (k=nx2; k--; dpixin0+=ny1, pixout0++)
    {
    maskt = obj->maskbuf;
    nmaskt = obj->nmaskbuf;
    startt = obj->startbuf;
    pixout = pixout0;
    for (j=ny2; j--; pixout+=w2)
      {
      dpixin = dpixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(dpixin++);
      *pixout += (PIXTYPE)val;
      }
    }

  return RETURN_OK;
  }


/****** make_kernel **********************************************************
PROTO	void make_kernel(double pos, double *kernel, interpenum interptype)
PURPOSE	Conpute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	07/04/2008
 ***/
static void	make_kernel(double pos, double *kernel, interpenum interptype)
  {
   double	x, val, sinx1,sinx2,sinx3,cosx1;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR)
    {
    *(kernel++) = 1.0-pos;
    *kernel = pos;
    }
  else if (interptype == INTERP_LANCZOS2)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/2.0*(pos+1.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/2.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0;
      val += (*kernel = cosx1/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else if (interptype == INTERP_LANCZOS3)
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/3.0*(pos+2.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx2=-0.5*sinx1-0.866025403785*cosx1)
				/ (x*x));
      x += PI/3.0;
      val += (*(kernel++) = (sinx3=-0.5*sinx1+0.866025403785*cosx1)
				/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx1/(x*x));
      x += PI/3.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/3.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }
  else
    {
    if (pos<1e-5 && pos>-1e5)
      {
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 1.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *(kernel++) = 0.0;
      *kernel = 0.0;
      }
    else
      {
      x = -PI/4.0*(pos+3.0);
#ifdef HAVE_SINCOS
      sincos(x, &sinx1, &cosx1);
#else
      sinx1 = sin(x);
      cosx1 = cos(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/4.0;
      val +=(*(kernel++) = -(sinx2=0.707106781186*(sinx1+cosx1))
				/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = cosx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -(sinx3=0.707106781186*(cosx1-sinx1))/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/4.0;
      val += (*kernel = sinx3/(x*x));
      val = 1.0/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
      }
    }

  return;
  }


/********************************* add_image *********************************/
/*
Copy a small part of the image.
*/
int     add_image(PIXTYPE *pix1, int w1, int h1, PIXTYPE *pix2, int w2, int h2,
		int ix, int iy, float amplitude)
  {
   PIXTYPE	*pix;
   int		x,y, xmin,xmax,ymin,ymax, w, dw1;

/* Set and check image boundaries */
  w = w1;
  ymin = iy-h1/2;
  ymax = ymin + h1;
  xmin = ix-w1/2;
  xmax = xmin + w1;
  if (ymax<=0 || xmax<=0 || xmin>=w2 || ymin>=h2)
    return RETURN_ERROR;

  if (ymin<0)
    {
    pix1 -= ymin*w1;
    ymin = 0;
    }
  if (ymax>h2)
    ymax = h2;

  if (xmax>w2)
    {
    w -= xmax-w2;
    xmax = w2;
    }
  if (xmin<0)
    {
    pix1 -= xmin;
    w += xmin;
    xmin = 0;
    }

  dw1 = w1-w;
/* Add the right pixels to the destination */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&addimagemutex);
#endif
  for (y=ymin; y<ymax; y++, pix1 += dw1)
    {
    pix = pix2+y*w2+xmin;
    for (x=w; x--;)
      *(pix++) += amplitude**(pix1++);
    }
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&addimagemutex);
#endif

   return RETURN_OK;
  }


