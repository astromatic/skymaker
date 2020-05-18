/*
*				image.c
*
* Manipulate image rasters.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		18/05/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _GNU_SOURCE
#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
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

static void	image_make_kernel(float pos, float *kernel,
			interpenum interptype);
static float	image_interpolate_pix(double *posin, PIXTYPE *pixin,
			int *naxisn, interpenum interptype);

int		interp_kernwidth[5]={1,2,4,6,8};

/****** image_add ***********************************************************
PROTO	void image_add(PIXTYPE *pix1, int w1, int h1,
		PIXTYPE *pix2, int w2, int h2,
		int ix, int iy, float amplitude)
PURPOSE	Copy and paste an image at a given position in another image.
INPUT   Pointer to input raster,
	input width,
	input height,
	pointer to output raster,
	output width,
	output height,
        shift in x (integer),
	shift in y (integer),
	flux scaling factor.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2020
 ***/
int     image_add(PIXTYPE *pix1, int w1, int h1,
		PIXTYPE *pix2, int w2, int h2,
		int ix, int iy, float amplitude) {

   PIXTYPE	*pix;
   int		x,y, xmin,xmax,ymin,ymax, w, dw1;

// Set and check image boundaries
  w = w1;
  ymin = iy - h1/2;
  ymax = ymin + h1;
  xmin = ix - w1/2;
  xmax = xmin + w1;
  if (ymax<=0 || xmax<=0 || xmin>=w2 || ymin>=h2)
    return RETURN_ERROR;

  if (ymin<0) {
    pix1 -= ymin*w1;
    ymin = 0;
  }
  if (ymax>h2)
    ymax = h2;

  if (xmax>w2) {
    w -= xmax-w2;
    xmax = w2;
  }
  if (xmin<0) {
    pix1 -= xmin;
    w += xmin;
    xmin = 0;
  }

  dw1 = w1-w;
// Add the right pixels to the destination
  for (y=ymin; y<ymax; y++, pix1 += dw1) {
    pix = pix2 + y*w2 + xmin;
    for (x=w; x--;)
      *(pix++) += amplitude * *(pix1++);
  }

  return RETURN_OK;
}


/****** image_resample_obj ***************************************************
PROTO	int image_resample_obj(PIXTYPE *pix1, int w1, int h1, objstruct *obj,
			double dx, double dy, double step2)
PURPOSE Scale and shift a small object image through interpolation.
	Image parts which lie outside boundaries are ignored.
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
VERSION 20/05/2020
 ***/
int	image_resample_obj(PIXTYPE *pix1, int w1, int h1, objstruct *obj,
			float dx, float dy, float step2) {

   interpenum	interp_type;
   PIXTYPE	*pixin,*pixin0, *pixout,*pixout0;
   float	kernel[INTERP_MAXKERNELWIDTH], *kernelt,
		*dpixin,*dpixin0, *dpixout,*dpixout0, *maskt,
		xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, val, ddix2,ddiy2;
   int		i,j,k,n,t, *startt, *nmaskt,
		ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
		ix,iy, ix1,iy1, w2,h2, interpw, imax;

  interp_type = INTERP_TYPE;
  interpw = interp_kernwidth[interp_type];
  w2 = obj->subsize[0];
  h2 = obj->subsize[1];
  xc1 = (double)(w1/2);		// Im1 center x-coord*/
  xc2 = (double)(w2/2);		// Im2 center x-coord*/
  xs1 = xc1 + dx - xc2*step2;	// Im1 start x-coord
  if ((int)xs1 >= w1)
    return RETURN_ERROR;
  ixs2 = 0;			// Int part of Im2 start x-coord
  if (xs1<0.0) {
    ddix2 = 1.0 - xs1/step2;
//-- Simply leave here if the images do not overlap in x
    if (ddix2 >= (float)w2)
      return RETURN_ERROR;
    dix2 = (int)ddix2;
    ixs2 += dix2;
    xs1 += dix2*step2;
  }
  nx2 = (int)((w1-xs1)/step2);	// nb of interpolated Im2 pixels along x
  if (nx2>(ix2=w2-ixs2))
    nx2 = ix2;
  if (!nx2)
    return RETURN_ERROR;
  yc1 = (double)(h1/2);		// Im1 center y-coord
  yc2 = (double)(h2/2);		// Im2 center y-coord
  ys1 = yc1 + dy - yc2*step2;		// Im1 start y-coord
  if ((int)ys1 >= h1)
    return RETURN_ERROR;
  iys2 = 0;			// Int part of Im2 start y-coord
  if (ys1<0.0) {
    ddiy2 = 1.0 - ys1/step2;
//-- Simply leave here if the images do not overlap in y
    if (ddiy2 >= (float)h2)
      return RETURN_ERROR;
    diy2 = (int)ddiy2;
    iys2 += diy2;
    ys1 += diy2*step2;
  }
  ny2 = (int)((h1-ys1)/step2);	// nb of interpolated Im2 pixels along y
  if (ny2>(iy2=h2-iys2))
    ny2 = iy2;
  if (!ny2)
    return RETURN_ERROR;

// Set the yrange for the x-resampling with some margin for interpolation
  iys1a = (int)ys1;		// Int part of Im1 start y-coord with margin
  hmh = hmw = interpw/2 - 1;	// Interpolant start
  if (iys1a<0 || ((iys1a -= hmh)< 0))
    iys1a = 0;
  ny1 = (int)(ys1+ny2*step2)+interpw-hmh;	// Interpolated Im1 y size
  if (ny1>h1)					// with margin
    ny1 = h1;
// Express everything relative to the effective Im1 start (with margin)
  ny1 -= iys1a;
  ys1 -= (double)iys1a;

// Allocate interpolant stuff for the x direction
  if (obj->buf1size < nx2) {
    if (!obj->buf1size) {
      free(obj->maskbuf);
      free(obj->nmaskbuf);
      free(obj->startbuf);
    }
    obj->buf1size = nx2;
    QMALLOC16(obj->maskbuf, float, nx2*interpw);	// Interpolation masks
    QMALLOC16(obj->nmaskbuf, int, nx2);		// Interpolation mask sizes
    QMALLOC16(obj->startbuf, int, nx2);		// Int part of Im1 conv starts
  }

// Compute the local interpolant and data starting points in x
  x1 = xs1;
  maskt = obj->maskbuf;
  nmaskt = obj->nmaskbuf;
  startt = obj->startbuf;
  for (j=nx2; j--; x1+=step2) {
    ix = (ix1=(int)x1) - hmw;
    image_make_kernel(x1-ix1, kernel, interp_type);
    kernelt = kernel;
    if (ix < 0) {
      n = interpw+ix;
      kernelt -= ix;
      ix = 0;
    } else
      n = interpw;
    if (n>(t=w1-ix))
      n=t;
    *(startt++) = ix;
    *(nmaskt++) = n;
    for (i=n; i--;)
      *(maskt++) = *(kernelt++);
  }

  if (obj->buf2size < nx2*ny1) {
    if (!obj->buf2size)
      free(obj->buf2);
    obj->buf2size = nx2*ny1;
    QCALLOC16(obj->buf2, float, obj->buf2size);	// Intermediary frame-buffer
  }

// Make the interpolation in x (this includes transposition)
  pixin0 = pix1+iys1a*w1;
  dpixout0 = obj->buf2;
  for (k=ny1; k--; pixin0+=w1, dpixout0++) {
    maskt = obj->maskbuf;
    nmaskt = obj->nmaskbuf;
    startt = obj->startbuf;
    dpixout = dpixout0;
    for (j=0; j<nx2; j++,dpixout+=ny1) {
      pixin = pixin0+*(startt++);
      val = 0.0;
      imax = *(nmaskt++);
#pragma omp simd reduction(+:val)
      for (i=0; i<imax; i++)
        val += *(maskt++)**(pixin++);
      *dpixout = val;
    }
  }

// Reallocate interpolant stuff for the y direction
  if (obj->buf1size < ny2) {
    obj->buf1size = ny2;
    free(obj->maskbuf);
    free(obj->nmaskbuf);
    free(obj->startbuf);
    QMALLOC16(obj->maskbuf, float, ny2*interpw);
    QMALLOC16(obj->nmaskbuf, int, ny2);
    QMALLOC16(obj->startbuf, int, ny2);
  }

// Compute the local interpolant and data starting points in y
  y1 = ys1;
  maskt = obj->maskbuf;
  nmaskt = obj->nmaskbuf;
  startt = obj->startbuf;
  for (j=ny2; j--; y1+=step2) {
    iy = (iy1=(int)y1) - hmh;
    image_make_kernel(y1-iy1, kernel, interp_type);
    kernelt = kernel;
    if (iy < 0) {
      n = interpw+iy;
      kernelt -= iy;
      iy = 0;
    } else
      n = interpw;
    if (n>(t=ny1-iy))
      n=t;
    *(startt++) = iy;
    *(nmaskt++) = n;
    for (i=n; i--;)
      *(maskt++) = *(kernelt++);
  }

// Initialize destination buffer to zero
  memset(obj->subimage, 0, (size_t)(w2*h2)*sizeof(PIXTYPE));

// Make the interpolation in y and transpose once again
  dpixin0 = obj->buf2;
  pixout0 = obj->subimage+ixs2+iys2*w2;
  for (k=nx2; k--; dpixin0+=ny1, pixout0++) {
    maskt = obj->maskbuf;
    nmaskt = obj->nmaskbuf;
    startt = obj->startbuf;
    pixout = pixout0;
    for (j=ny2; j--; pixout+=w2) {
      dpixin = dpixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(dpixin++);
      *pixout += (PIXTYPE)val;
    }
  }

  return RETURN_OK;
}


/****** image_resample_lin **************************************************
PROTO	int image_resample_lin(PIXTYPE *pix1, int w1, int h1,
			PIXTYPE *pix2, int w2, int h2, 
			double dx, double dy, double *lin)

PURPOSE Apply shifting and linear geometric transformation to a small image
	through interpolation.
INPUT   Pointer to input raster,
	input width,
	input height,
	pointer to output raster. 
	output width,
	output height,
	shift in x,
	shift in y,
	pointer to the linear transformation matrix.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/05/2020
 ***/

int	image_resample_lin(PIXTYPE *pix1, int w1, int h1,
			PIXTYPE *pix2, int w2, int h2, 
			double dx1, double dy1, double *lin) {

   double	invlin[4], pos1[2],
		det, invdet;
   int		dim1[2],
		interp_type, iy2, nx2, ny2;


// Invert 2x2 input transformation matrix
  det = lin[0] * lin[3] - lin[1] * lin[2];
  if (fabs(det) < SMALL)
    return RETURN_ERROR;

  invdet = 1.0 / det;
  invlin[0]  =  invdet * lin[3];
  invlin[1] *= -invdet;
  invlin[2] *= -invdet;
  invlin[3]  =  invdet * lin[0];

  interp_type = INTERP_TYPE;
  dim1[0] = w1;
  dim1[1] = h1;
// Incorporate initial x2 and y2 positions into the x1 and y2 shifts
  dx1 -= invlin[0] * (w2/2);
  dy1 -= invlin[2] * (w2/2);
// Loop on all output pixels (reverse mapping)
// Pixels outside of bounds are set to 0
  for (ny2=h2, iy2=-(h2/2); ny2--; iy2++) {
    pos1[0] = invlin[1] * iy2 + dx1;
    pos1[1] = invlin[3] * iy2 + dy1;
    for (nx2=w2; nx2--;) {
      pos1[0] += invlin[0];
      pos1[1] += invlin[2];
      *(pix2++) = image_interpolate_pix(pos1, pix1, dim1, interp_type);
    }
  }
      
  return RETURN_OK;
}


/****** image_interpolate_pix *************************************************
PROTO	float image_interpolate_pix(double *posin, PIXTYPE *pixin, int *naxisn,
		interpenum interptype)
PURPOSE	Interpolate pixel values at a given position in 2D raster.
INPUT	Pointer to input position vector,
	pointer to input image raster,
	pointer to output image raster,
	pointer to raster shape vector,
	interpolation type.
OUTPUT	Interpolated value.
NOTES	Pixels outside limits are set to 0.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2020
 ***/
static float	image_interpolate_pix(double *posin, PIXTYPE *pixin,
			int *naxisn, interpenum interptype) {

   float	buffer[INTERP_MAXKERNELWIDTH],
		kernel[INTERP_MAXKERNELWIDTH], fpos[2],
		*kvector, *pix, *pixout,
		val;
   int		fac, ival, kwidth, start, width, step,
		i,j, n;

  kwidth = interp_kernwidth[interptype];
  start = 0;
  fac = 1;
  for (n=0; n<2; n++) {
    val = *(posin++);
    width = naxisn[n];
//-- Get the integer part of the current coordinate or nearest neighbour
    ival = (interptype==INTERP_NEARESTNEIGHBOUR)? (int)(val-0.50001):(int)val;
//-- Store the fractional part of the current coordinate
    fpos[n] = (float)(val - ival);
//-- Check if interpolation start/end exceed image boundary...
    ival -= kwidth/2;
    if (ival < 0 || ival+kwidth <= 0 || ival+kwidth > width)
      return 0.0;
//-- Update starting pointer
    start += ival*fac;
//-- Update step between interpolated regions
    fac *= width;
  }

// First step: interpolate along NAXIS1 from the data themselves
  image_make_kernel(fpos[0], kernel, interptype);
  step = naxisn[0] - kwidth;
  pix = pixin + start;
  pixout = buffer;
  for (j=kwidth; j--;) {
    val = 0.0;
    kvector = kernel;
    for (i=kwidth; i--;)
      val += *(kvector++)**(pix++);
    *(pixout++) = val;
    pix += step;
  }

// Second step: interpolate along NAXIS2 from the interpolation buffer
  image_make_kernel(fpos[1], kernel, interptype);
  pix = buffer;
  val = 0.0;
  kvector = kernel;
  for (i=kwidth; i--;)
    val += *(kvector++)**(pix++);

  return val;
}


/****** image_make_kernel ****************************************************
PROTO	void image_make_kernel(float pos, float *kernel, interpenum interptype)
PURPOSE	Compute interpolation-kernel data
INPUT	Position,
	Pointer to the output kernel data,
	Interpolation method.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2020
 ***/
static void	image_make_kernel(float pos, float *kernel,
				interpenum interptype) {

   float	x, val, sinx1,sinx2,sinx3,cosx1;

  if (interptype == INTERP_NEARESTNEIGHBOUR)
    *kernel = 1;
  else if (interptype == INTERP_BILINEAR) {
    *(kernel++) = 1.0f-pos;
    *kernel = pos;
  } else if (interptype == INTERP_LANCZOS2) {
    if (pos<1e-5f && pos>-1e-5f) {
      *(kernel++) = 0.0f;
      *(kernel++) = 1.0f;
      *(kernel++) = 0.0f;
      *kernel = 0.0f;
    } else {
      x = -PI/2.0f*(pos+1.0f);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/2.0f;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/2.0f;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/2.0f;
      val += (*kernel = cosx1/(x*x));
      val = 1.0f/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
    }
  } else if (interptype == INTERP_LANCZOS3) {
    if (pos<1e-5f && pos>-1e-5f) {
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 1.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *kernel = 0.0f;
    } else {
      x = -PI/3.0f*(pos+2.0f);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/3.0f;
      val += (*(kernel++) = (sinx2=-0.5f*sinx1-0.866025403785f*cosx1)
				/ (x*x));
      x += PI/3.0f;
      val += (*(kernel++) = (sinx3=-0.5f*sinx1+0.866025403785f*cosx1)
				/(x*x));
      x += PI/3.0f;
      val += (*(kernel++) = sinx1/(x*x));
      x += PI/3.0f;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/3.0f;
      val += (*kernel = sinx3/(x*x));
      val = 1.0f/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
    }
  } else if (interptype == INTERP_LANCZOS4) {
    if (pos<1e-5f && pos>-1e-5f) {
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 1.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *(kernel++) = 0.0f;
      *kernel = 0.0f;
    } else {
      x = -PI/4.0f*(pos+3.0f);
#ifdef HAVE_SINCOSF
      sincosf(x, &sinx1, &cosx1);
#else
      sinx1 = sinf(x);
      cosx1 = cosf(x);
#endif
      val = (*(kernel++) = sinx1/(x*x));
      x += PI/4.0f;
      val +=(*(kernel++) = -(sinx2=0.707106781186f*(sinx1+cosx1))
				/(x*x));
      x += PI/4.0f;
      val += (*(kernel++) = cosx1/(x*x));
      x += PI/4.0f;
      val += (*(kernel++) = -(sinx3=0.707106781186f*(cosx1-sinx1))/(x*x));
      x += PI/4.0f;
      val += (*(kernel++) = -sinx1/(x*x));
      x += PI/4.0f;
      val += (*(kernel++) = sinx2/(x*x));
      x += PI/4.0f;
      val += (*(kernel++) = -cosx1/(x*x));
      x += PI/4.0f;
      val += (*kernel = sinx3/(x*x));
      val = 1.0f/val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *(kernel--) *= val;
      *kernel *= val;
    }
  } else
    error(EXIT_FAILURE, "*Internal Error*: Unknown interpolation type in ",
		"make_kernel()");

  return;
}


