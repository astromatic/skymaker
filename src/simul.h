/*
*				simul.h
*
* Include file for simul.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 2003-2015 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		06/07/2015
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SIMUL_H_
#define _SIMUL_H_

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

/*----------------------------- Internal constants --------------------------*/
#define	PSF_INTERPW	2		/* Footprint of PSF interpolant */
#define	PSF_NINTERP	(PSF_INTERPW*PSF_INTERPW)	/* Interpolant area */
#define	PSF_NVARSNAP	11		/* Number of PSF snapshots per axis */
#define	PSF_NVARORDER	3		/* Number of PSF variation orders (-1)*/
#define	PSF_VARTHRESH	1e-9		/* Aberration significance threshold */

#define	QUANT_ACCURACY	1.0   /* flux quantization in units of noise RMS */

/*------------------------------ Type definitions ---------------------------*/

typedef enum	{PUPIL_REAL,PUPIL_IMAGINARY,PUPIL_MODULUS,PUPIL_PHASE,PUPIL_MTF,
		PSF_MTF,PSF_FULLRES,PSF_FINALRES, SKY_NONOISE, SKY,
		GRID_NONOISE, GRID} imatypenum;

typedef enum	{PSF_INTERNAL, PSF_FILE} psftypenum;

typedef enum	{NO_SEEING, LONG_EXPOSURE, SHORT_EXPOSURE} seeingtypenum;

typedef enum	{CENTER_UPPERHALF, CENTER_LOWERHALF, CENTER_HALF,
		CENTER_CENTROID, CENTER_CENTROID_COMMON,
		CENTER_PEAK}	centertypenum;

typedef enum	{NO_TRACKERR, DRIFT, JITTER}	tracktypenum;

typedef enum	{CORREL_NONE, CORREL_NOISE, CORREL_ALL}	correltypenum;

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
/* Image itself */
  char		filename[MAXCHAR];	/* Image filename */
  char		headname[MAXCHAR];	/* Header filename */
  catstruct	*cat;			/* cat (FITS file) structure */ 
  char		inlistname[MAXCHAR];	/* Input list filename */
  FILE		*inlistfile;		/* Input list file */
  char		outlistname[MAXCHAR];	/* Output list filename */
  FILE		*outlistfile;		/* Output list file */
  PIXTYPE	*image;			/* Pointer to the image pixel map */
  imatypenum	imatype;		/* Image type */
/*------ Astrometry */
  int		wcsflag;		/* Use WCS coordinates? */
  wcsstruct	*wcs;			/* World Coordinate System */
  int		imasize[4];		/* Dimension of the image */
  int		nimasize;		/* Number of arguments */
  int		fimasize[2];		/* Dimension of the FULL pixmap */
  int		margin[2];		/* Margin around the actual image */
  int		mscan[2];		/* Number of microscanning steps */
  int		nmscan;			/* Number of arguments */
  double	pixscale[2];		/* Pixel scale (in arcsec) */
  int		npixscale;		/* Number of arguments */
/*------ Photometry */
  double	lambdaeq;		/* Central wavelength (in microns) */
  double	wellcap;		/* Full well capacity (e-) */
  double	satlev;			/* Saturation level (ADU) */
  double	gain;			/* Gain (e-/ADU) */
  double	ron;			/* Read-out noise (e-) */
  double	magzero;		/* Mag. zero-point "ADU per second" */
  double	magzero2;		/* Mag. zero-point "e- per exposure" */
  double	expotime;		/* Exposure time (s) */
  double	magback;		/* Background surface brightness */
  double	minquant;		/* Minimum quantization step */

/*------ PSF */
  psftypenum	psftype;		/* PSF type */
  char		psfname[MAXCHAR];	/* PSF file name */
  double	seeing;			/* FWHM of the PSF (in arcsec) */
  double	ro;			/* Fried parameter (in m) */
  seeingtypenum	psfseeingtype;		/* Seeing type */
  double	psfmotion;		/* PSF motion amplitude in arcsec */
  double	*psfmot[2];		/* PSF motion vector */
  double	psfhalosb;		/* SB at 1' for a 0-mag star */
  double	psfoversamp;		/* Oversampling of the PSF */
  int		npsf;			/* Number of PSF components */
  PIXTYPE	*psf;			/* Pointer to the PSF pixel map(s) */
  PIXTYPE	**psfdft;		/* Pointers to DFTs of the PSF */
  int		psfsize[5];		/* Dimensions of the PSF */
  int		npsfsize;		/* Number of arguments */
  double	*psfdpos[2];		/* Relative PSF center coordinates*/
  int		psfnarms;		/* Number of spider arms */
  double	psfdm1;			/* Diameter of the primary mirror (m)*/
  double	psfdm2;			/* Diam. of the 2nd mir. support (m) */
  double	psfarmw;		/* Thickness of the spider arms (m)*/
  double	psfarmang;		/* Pos. angle of spider arms (deg) */
  double	psfd80defoc[PSF_NVARORDER];/* FWHM induced by defocus (arcsec)*/
  double	psfdefocc[2];		/* Defocalisation center coords */
  double	psfd80spher[PSF_NVARORDER];/* Spherical d80% diameter (arcsec)*/
  double	psfspherc[2];		/* Spherical aberr. center coords */
  double	psfd80comax[PSF_NVARORDER];/* X-coma d80% diameter (arcsec) */
  double	psfd80comay[PSF_NVARORDER];/* Y-coma d80% diameter (arcsec) */
  double	psfcomac[2];		/* Coma center coords */
  double	psfd80ast00[PSF_NVARORDER];/* 0 deg. astigmatism d80% (arcsec)*/
  double	psfd80ast45[PSF_NVARORDER];/* 45 deg astigmatism d80% (arcsec)*/
  double	psfastc[2];		/* Astigmatism center coords */
  double	psfd80tri00[PSF_NVARORDER];/* 0 deg. triangular d80% (arcsec) */
  double	psfd80tri30[PSF_NVARORDER];/* 30 deg. triangular d80% (arcsec)*/
  double	psftric[2];		/* Triangular aberr center coords */
  double	psfd80qua00[PSF_NVARORDER];/* 0 deg. quadratic d80% (arcsec) */
  double	psfd80qua22[PSF_NVARORDER];/* 22.5 deg quadratic d80% (arcsec)*/
  double	psfquac[2];		/* Quadratic aberr center coords */
  tracktypenum	psftracktype;		/* Tracking type */
  double	psftrackmaj;		/* Maximum RMS tracking error (") */
  double	psftrackmin;		/* Minimum RMS tracking error (") */
  double	psftrackang;		/* CC angle of maj. axis (deg) */
  double	psfarea;		/* Typical area of the PSF */
/*------ Aureole */
  PIXTYPE	*aureole;		/* Pointer to the aureole pixel map */
  int		aursize[2];		/* Dim (/nmscan) of the aureole map */
  int		aurange;		/* Maximum radius of the aurole */
/*------ Stellar field */
  double	scountdens;		/* Star nb per sq. deg. at mag.lim. */
  double	scountslope;		/* Diff. star count slope (dexp/mag) */
  double	maglim[2];		/* Brightest and faintest mag. allowed*/
  int		nmaglim;		/* Number of arguments */

/*------ Random generator */
  int		psfmotionseed;		/* Seed for PSF motion */
  int		starposseed;		/* Seed for star positions */
/*------ Grid parameters */
  int		gridindex;		/* Position in grid */
  int		gridstep;		/* Step between objects */
  int		ngrid[2];		/* Number of objects per axis */
  int		gridoffset[2];		/* Grid offset on each axis */
/* Misc */
        }	simstruct;

/*-------------------------------- protos -----------------------------------*/
extern simstruct	*sim_init(void);

extern void		sim_end(simstruct *sim);

#endif
