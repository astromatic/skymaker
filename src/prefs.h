/*
*				prefs.h
*
* Include file for prefs.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2020 IAP/CNRS/SorbonneU
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

#ifndef _PSF_H_
#include "psf.h"
#endif

#ifndef _LIST_H_
#include "list.h"
#endif

#ifndef _SIMUL_H_
#include "simul.h"
#endif

#ifndef _WEIGHT_H_
#include "weight.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

//----------------------------- Internal constants --------------------------

#define		MAXLIST		32	// max. nb of list members 

//------------------------------- preferences -------------------------------
typedef struct
  {
//------ Image itself 
  char		filename[MAXCHAR];	// Image filename 
  char		headname[MAXCHAR];	// Header filename 
  char		inlistname[MAXCHAR];	// Input list filename 
  char		outlistname[MAXCHAR];	// Output list filename 
  imatypenum	imatype;		// Image type 
//----- Weighting 
  char		weight_name[MAXCHAR];	// Input weight filename 
  weightenum	weight_type;		// Weighting scheme 
  double	weight_thresh;		// Bad pixel weight threshold 
  int		nweight_thresh;		// Number of arguments 
  double	weight_factor;		// Noise scaling factor 
  int		weightgain_flag;	// Weight proportional to gain? 
//------ Astrometry 
  listcoordenum	listcoord_type;		// Use World coordinates? 
  int		imasize[2];		// Dimension of the image 
  int		nimasize;		// Number of arguments 
  int		mscan[2];		// Number of microscanning steps 
  int		nmscan;			// Number of arguments 
  double	pixscale[2];		// Pixel scale (in arcsec) 
  int		npixscale;		// Number of arguments 
//------ Source list 
  double	listmag_limits[2];	// Brightest and faintest mag. allowed
  int		nlistmag_limits;	// Number of arguments 
//------ Detector 
  double	lambdaeq;		// Central wavelength (in microns) 
  double	wellcap;		// Full well capacity (e-) 
  double	satlev;			// Saturation level (ADU) 
  double	gain;			// Gain (e-/ADU) 
  double	ron;			// Read-out noise (e-) 
  double	magzero;		// Mag. zero-point "ADU per second" 
  double	expotime;		// Exposure time (s) 
  double	magback;		// Background surface brightness 
//------ Noise correlation 
  corrtypenum	corr_type;		// Correlation type 
  corrfuncenum  corr_func_type;		// Correlation function type 
  double	corr_scale;		// Scaling of the correl. function 
//------ PSF 
  psftypenum	psftype;		// PSF type 
  char		psfname[MAXCHAR];	// PSF file name 
  double	seeing;			// FWHM of the PSF (in arcsec) 
  seeingtypenum	psfseeingtype;		// Seeing type 
  centertypenum	psfcentertype;		// PSF centering type 
  double	psfhalosb;		// SB at 1' for a 0-mag star 
  double	psfoversamp;		// Oversampling of the PSF 
  int		psfsize[2];		// Dimensions of the PSF 
  int		npsfsize;		// Number of arguments 
  int		psfnarms;		// Number of spider arms 
  double	psfdm1;			// Diameter of the primary mirror (m)
  double	psfdm2;			// Diam. of the 2nd mir. support (m) 
  double	psfarmw;		// Thickness of the spider arms (m)
  double	psfarmang;		// Pos. angle of spider arms (deg) 
  double	psfd80defoc[PSF_NVARORDER];// FWHM induced by defocus (arcsec)
  int		npsfd80defoc;		// Number of arguments 
  double	psfdefocc[2];		// PSF defocus center 
  int		npsfdefocc;		// Number of arguments 
  double	psfd80spher[PSF_NVARORDER];// Spherical d80% diameter (arcsec)
  int		npsfd80spher;		// Number of arguments 
  double	psfspherc[2];		// PSF spherical aber. center 
  int		npsfspherc;		// Number of arguments 
  double	psfd80comax[PSF_NVARORDER];// X-coma d80% diameter (arcsec) 
  int		npsfd80comax;		// Number of arguments 
  double	psfd80comay[PSF_NVARORDER];// Y-coma d80% diameter (arcsec) 
  int		npsfd80comay;		// Number of arguments 
  double	psfcomac[2];		// PSF coma center 
  int		npsfcomac;		// Number of arguments 
  double	psfd80ast00[PSF_NVARORDER];// 0 deg. astigmatism d80% (arcsec)
  int		npsfd80ast00;		// Number of arguments 
  double	psfd80ast45[PSF_NVARORDER];// 45 deg astigmatism d80% (arcsec)
  int		npsfd80ast45;		// Number of arguments 
  double	psfastc[2];		// PSF astigmatism center 
  int		npsfastc;		// Number of arguments 
  double	psfd80tri00[PSF_NVARORDER];// 0 deg. triangular d80% (arcsec) 
  int		npsfd80tri00;		// Number of arguments 
  double	psfd80tri30[PSF_NVARORDER];// 30 deg. triangular d80% (arcsec)
  int		npsfd80tri30;		// Number of arguments 
  double	psftric[2];		// PSF triangular aber center 
  int		npsftric;		// Number of arguments 
  double	psfd80qua00[PSF_NVARORDER];// 0 deg. quadratic d80% (arcsec) 
  int		npsfd80qua00;		// Number of arguments 
  double	psfd80qua22[PSF_NVARORDER];// 22.5 deg quadratic d80% (arcsec)
  int		npsfd80qua22;		// Number of arguments 
  double	psfquac[2];		// PSF quadratic aberration center 
  int		npsfquac;		// Number of arguments 
  tracktypenum	psftracktype;		// Tracking type 
  double	psftrackmaj;		// Maximum RMS tracking error (") 
  double	psftrackmin;		// Minimum RMS tracking error (") 
  double	psftrackang;		// CC angle of maj. axis (deg) 
//------ Aureole 
  int		aurange;		// Maximum radius of the aurole 
//------ Stellar field 
  double	scountdens;		// Star nb per sq. deg. at mag.lim. 
  double	scountslope;		// Diff. star count slope (dexp/mag) 
  double	maglim[2];		// Brightest and faintest mag. allowed
  int		nmaglim;		// Number of arguments 
//------ Rasters 
  char		raster_pattern[MAXCHAR];	// File pattern for rasters
//------ Random generator 
  int		psfmotionseed;		// Seed for PSF motion 
  int		starposseed;		// Seed for star positions 
//------ Multithreading 
  int		nthreads;		// Number of active threads 
//------ Grid parameters 
  int		grid_size;		// Grid size (pixels) in GRID mode 
// Misc 
  enum {QUIET, NORM, LOG, FULL} verbose_type;	// display type 
  enum {NONE, ASTROMATIC}  xml_type;		// XML output type 
  char		xml_name[MAXCHAR];		// XML file name 
  char 		sdate_start[12];		// SkyMaker start date 
  char		stime_start[12];		// SkyMaker start time 
  char		sdate_end[12];			// SkyMaker end date 
  char		stime_end[12];			// SkyMaker end time 
  double	time_diff; 			// Execution time 
  int		nobj;				// Number of sources added 
        }	prefstruct;

prefstruct	prefs;

//-------------------------------- protos -----------------------------------
extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);

#endif
