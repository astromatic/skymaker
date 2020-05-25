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
*	Last modified:		25/05/2020
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
VERSION	25/05/2020
 ***/
int	make_raster(simstruct *sim, objstruct *obj) {

   catstruct	*cat;
   tabstruct	*tab;
   char		str[MAXCHAR];
   PIXTYPE	*sub, *subt, *psfdft, *raster,
		invflux;
   double	dpos[2], geo[4], jac[4], lin[4], invlin[4],
		osamp, size, cang, sang, det, invdet,
		flux, flux2, dx,dy, scale, expo, n, dval;
   int		i, subwidth,subheight, suborder,
		rasterwidth, rasterheight, rastersize,
		nsub, nsub2, nsubo, memnsub, oversamp;

  osamp = sim->psfoversamp;
  size = obj->raster_size / sim->pixscale[0];
  n = 0.0;	// avoid gcc -Wall warnings

// Convert magnitude to linear units
  if (!obj->flux) {
    expo = 0.4*(sim->magzero2-obj->mag);
    if (expo>-30.0)
      obj->flux = DEXP(expo);
    if (!obj->flux)
      return RETURN_ERROR;
  }

// Set mask size limits
  if (size < SMALL) {
    NFPRINTF(OUTPUT, "");
    sprintf(str, "Final image size too small for raster #%d at (%.1f,%.1f): ",
	obj->raster_index, obj->pos[0], obj->pos[1]);
    warning(str, "skipped");
    return RETURN_ERROR;
  }

// Read file
  if (obj->type==300) {
    if (obj->raster_index <= 0) {
      error(EXIT_FAILURE, "Raster with negative index in ", sim->inlistname);
    }
    sprintf(str, prefs.raster_pattern, obj->raster_index);
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
  QFSEEK(cat->file, tab->bodypos, SEEK_SET, cat->filename);
  read_body(tab, raster, rastersize);
  free_cat(&cat, 1);

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
#ifdef HAVE_SINCOS
  sincos(obj->raster_posang * DEG, &sang, &cang);
#else
  cang = cos(obj->raster_posang * DEG);
  sang = sin(obj->raster_posang * DEG)
#endif
  scale = osamp * (double)size / rasterwidth;
  lin[0] = geo[0] = scale * cang;
  lin[1] = geo[1] = - scale * obj->raster_aspect * sang;
  lin[2] = geo[2] = scale * sang;
  lin[3] = geo[3] = scale * obj->raster_aspect * cang;

  if (sim->wcs) {
//-- Include Jacobian from the WCS transformation if in world coordinates
     double	invpixscale = 1.0 / (sim->pixscale[0] * ARCSEC / DEG);
    wcs_jacobian(sim->wcs, obj->pos, jac);
    for (i=0; i<4; i++)
      jac[i] *= invpixscale;
    lin[0] = jac[0] * geo[0] + jac[1] * geo[2];
    lin[1] = jac[0] * geo[1] + jac[1] * geo[3];
    lin[2] = jac[2] * geo[0] + jac[3] * geo[2];
    lin[3] = jac[2] * geo[1] + jac[3] * geo[3];
  }

// Compute determinant and invert 2x2 input transformation matrix
  det = fabs(lin[0] * lin[3] - lin[1] * lin[2]);
  if (det < SMALL) {
    NFPRINTF(OUTPUT, "");
    sprintf(str, "Null matrix determinant for raster #%d at (%.1f,%.1f): ",
	obj->raster_index, obj->pos[0], obj->pos[1]);
    warning(str, "skipped");
    return RETURN_ERROR;
  }

  invdet = 1.0 / det;
  invlin[0] =  invdet * lin[3];
  invlin[1] = -invdet * lin[1];
  invlin[2] = -invdet * lin[2];
  invlin[3] =  invdet * lin[0];

// Set oversampling factor
  oversamp = (int)sqrt(invdet);
  if (oversamp < 1)
    oversamp = 1;
// Render
  QFFTWF_CALLOC(sub, PIXTYPE, memnsub);
  flux = fabs(image_resample_lin(raster, rasterwidth, rasterheight,
			sub, subwidth, subheight,
			obj->raster_aspect, obj->raster_posang, invlin,
			oversamp));

  if (fabs(flux) <= SMALL) {
    NFPRINTF(OUTPUT, "");
    sprintf(str, "Null average for raster #%d at (%.1f,%.1f): ",
	obj->raster_index, obj->pos[0], obj->pos[1]);
    warning(str, "normalization dropped");
  } else {
//-- Combine and normalize (flux will change later, though)
    invflux = 1.0 / flux;
    for (i=n,subt=sub; i--;)
      *(subt++) *= invflux;
  }

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

  if (fabs(flux) <= SMALL) {
    obj->subfactor = 1.0;
    obj->noiseqarea = 0.0;
  } else {
    obj->subfactor = fabs(obj->flux/flux);
    obj->noiseqarea = flux*flux / flux2;
  }

  QFFTWF_FREE(sub);
  free(raster);

  return RETURN_OK;
  }


