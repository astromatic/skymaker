/*
*				rasters.c
*
* Generate images from external rasters.
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

#ifdef HAVE_SINCOSF
#define _GNU_SOURCE 
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include FFTW_H

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fft.h"
#include "galaxies.h"
#include "image.h"
#include "list.h"
#include "prefs.h"
#include "psf.h"
#include "rasters.h"
#include "simul.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

#ifdef USE_THREADS
extern pthread_mutex_t	fftmutex;
#endif

static double	gamm(double xx);

/****** make_raster *******************************************************
PROTO	int make_raster(simstruct *sim, objstruct *obj)
PURPOSE	Render a raster at a definite position in the image.
INPUT	Pointer to the sim structure,
	pointer to the object structure.
OUTPUT	-.
NOTES	Writes to an allocated image buffer, not directly to the image to
	allow multithreading.
AUTHOR	E. Bertin (IAP)
VERSION	18/05/2020
 ***/
int	make_raster(simstruct *sim, objstruct *obj)

  {
   catstruct	*cat;
   tabstruct	*tab;
   char		str[MAXCHAR];
   PIXTYPE	*sub, *subt, *psfdft, *raster,
		invflux;
   double	dpos[2], jac[4],
		osamp, size,
		flux,flux2, dx,dy, beq, dscale,expo, n, bn, ampfac, dval;
   int		i, subwidth,subheight, suborder,
		rasterwidth, rasterheight, rastersize,
		nsub,nsub2,nsubo,memnsub;

  osamp = sim->psfoversamp;
  size = obj->raster_size / sim->pixscale[0];
  n = 0.0;	// avoid gcc -Wall warnings

// Convert magnitude to linear units
  if (!obj->flux)
    {
    expo = 0.4*(sim->magzero2-obj->mag);
    if (expo>-30.0)
      obj->flux = DEXP(expo);
    if (!obj->flux)
      return RETURN_ERROR;
    }

// Set mask size limits
  if (size < SMALL) {
    NFPRINTF(OUTPUT, "");
    sprintf(str, "Final image size too small for raster at (%.1f,%.1f): ",
	obj->pos[0], obj->pos[1]);
    warning(str, "skipped");
    return RETURN_ERROR;
  }

// Read file
  if (obj->type==300) {
    if (obj->raster_index <= 0) {
      error(EXIT_FAILURE, "Raster with negative index in ", sim->inlistname);
    }
    sprintf(str, "%s_%09d.fits", prefs.raster_prefix, obj->raster_index);
    if (!(cat = read_cat(str))) {
      sprintf(gstr, "*Error*: %s not found", str);
      error(EXIT_FAILURE, gstr,"");
    }
  }
  tab = cat->tab;
  if (tab->naxis<2)
    error(EXIT_FAILURE, "*Error*: Not an image in ", cat->filename);

  rasterwidth = tab->naxisn[0];
  rasterheight = tab->naxisn[1];
  rastersize =  rasterwidth * rasterheight;
  QMALLOC(raster, PIXTYPE, rastersize);
  read_body(tab, raster, rastersize);

  i = 2*(int)(osamp*sqrt(size*size
	+16.0*sim->psfarea/(sim->pixscale[0]*sim->pixscale[1])));
  for (suborder=1; i>>=1; suborder++);
  if (suborder>=PSF_NORDER) {
    osamp /= (double)(1<<(suborder-PSF_NORDER));
    suborder = PSF_NORDER;
  }
  subheight = subwidth = 1<<suborder;
  nsub = subwidth*subheight;
  memnsub = (((subwidth>>1) + 1)<< 1) * subheight; // Provide margin for FFTW

// Compute (or retrieve) PSF DFT at this position
  psfdft = interp_dft(sim, suborder, obj->pos, dpos);

// Compute Jacobian
  if (sim->wcs) {
     double	invpixscale = 1.0 / (sim->pixscale[0] * ARCSEC / DEG);
    wcs_jacobian(sim->wcs, obj->pos, jac);
    for (i=0; i<4; i++)
      jac[i] *= invpixscale;
  }

// Render
  sub = NULL;
  QFFTWF_CALLOC(sub, PIXTYPE, memnsub);
  invflux = obj->raster_aspect / raster_raster(sub, subwidth, subheight,
	sim->wcs? jac : NULL,
	osamp*obj->raster_size / sim->pixscale[0], obj->raster_aspect,
	obj->raster_posang);

// Combine and normalize (flux will change later, though)
    for (i=nsub,subt=sub; i--;)
      *(subt++) *= invflux;

/* Truncate to avoid introducing anistropy in the vignet corners */
//  trunc_prof(sub, (double)(subwidth/2),(double)(subheight/2),
//		subwidth, subheight);
/* Convolve with the oversampled PSF */
  fft_conv(sub, psfdft, subwidth,subheight);
  if (sim->npsf>1)
    QFFTWF_FREE(psfdft);
  dx = (obj->pos[0] - 1.0 + sim->margin[0] - dpos[0])/sim->mscan[0];
  dy = (obj->pos[1] - 1.0 + sim->margin[1] - dpos[1])/sim->mscan[1];
  dx -= (double)(obj->subpos[0] = (int)(dx+0.49999));
  dy -= (double)(obj->subpos[1] = (int)(dy+0.49999));
// Truncate again
//  trunc_prof(sub, (double)(subwidth/2),(double)(subheight/2),
//		subwidth, subheight);
// Resample to lower resolution
  nsubo = obj->subsize[0]*obj->subsize[1];
  obj->subsize[0]  = obj->subsize[1] = (int)(subwidth/osamp+1.0);
  nsub2 = obj->subsize[0]*obj->subsize[1];
  if (!nsubo) {
    QMALLOC(obj->subimage, PIXTYPE, nsub2);
  }
  else if (nsub2>nsubo) {
    QREALLOC(obj->subimage, PIXTYPE, nsub2);
  }
  image_resample_obj(sub, subwidth, subheight, obj, -dx*osamp, -dy*osamp, osamp);
  flux = flux2 = 0.0;
  for (i=nsub2,subt=obj->subimage; i--;) {
    dval = (double)*(subt++);
    flux += dval;
    flux2 += dval*dval;
  }

  obj->subfactor = obj->flux/flux;
  obj->noiseqarea = flux*flux / flux2;

  QFFTWF_FREE(sub);
  free(raster);

  return RETURN_OK;
  }


/****** raster_raster ********************************************************
PROTO	double raster_raster(PIXTYPE *pix, int width, int height, double *jac,
		double size, double aspect, double posang)
PURPOSE	Rasterize an unnormalized 2D image.
INPUT	pointer to the raster,
	raster width,
	raster height,
	Jacobian array of the local astrometric deprojection (NULL if none),
	raster size (scale),
	raster aspect ratio,
	raster position angle.
OUTPUT	Total flux.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/05/2020
 ***/
double	raster_raster(PIXTYPE *pix, int width, int height, double *jac,
		double size, double aspect, double posang) {

  
  return	1.0;

};
