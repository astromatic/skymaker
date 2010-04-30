 /*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SkyMaker
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for prefs.c.
*
*	Last modify:	30/04/2010
**
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _SIMUL_H_
#include "simul.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXLIST		32	/* max. nb of list members */

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
/* Image itself */
  char		filename[MAXCHAR];	/* Image filename */
  char		headname[MAXCHAR];	/* Header filename */
  char		inlistname[MAXCHAR];	/* Input list filename */
  char		outlistname[MAXCHAR];	/* Output list filename */
  imatypenum	imatype;		/* Image type */
  int		imasize[2];		/* Dimension of the image */
  int		nimasize;		/* Number of arguments */
  int		mscan[2];		/* Number of microscanning steps */
  int		nmscan;			/* Number of arguments */
  double	pixscale[2];		/* Pixel scale (in arcsec) */
  int		npixscale;		/* Number of arguments */
  double	lambdaeq;		/* Central wavelength (in microns) */
  double	wellcap;		/* Full well capacity (e-) */
  double	satlev;			/* Saturation level (ADU) */
  double	gain;			/* Gain (e-/ADU) */
  double	ron;			/* Read-out noise (e-) */
  double	magzero;		/* Mag. zero-point "ADU per second" */
  double	expotime;		/* Exposure time (s) */
  double	magback;		/* Background surface brightness */

/*------ PSF */
  psftypenum	psftype;		/* PSF type */
  char		psfname[MAXCHAR];	/* PSF file name */
  double	seeing;			/* FWHM of the PSF (in arcsec) */
  seeingtypenum	psfseeingtype;		/* Seeing type */
  centertypenum	psfcentertype;		/* PSF centering type */
  double	psfhalosb;		/* SB at 1' for a 0-mag star */
  double	psfoversamp;		/* Oversampling of the PSF */
  int		psfsize[2];		/* Dimensions of the PSF */
  int		npsfsize;		/* Number of arguments */
  int		psfnarms;		/* Number of spider arms */
  double	psfdm1;			/* Diameter of the primary mirror (m)*/
  double	psfdm2;			/* Diam. of the 2nd mir. support (m) */
  double	psfarmw;		/* Thickness of the spider arms (m)*/
  double	psfarmang;		/* Pos. angle of spider arms (deg) */
  double	psfd80defoc;		/* FWHM induced by defocus (arcsec) */
  double	psfd80spher;		/* Spherical d80% diameter (arcsec) */
  double	psfd80comax;		/* X-coma d80% diameter (arcsec) */
  double	psfd80comay;		/* Y-coma d80% diameter (arcsec) */
  double	psfd80ast00;		/* 0 deg. astigmatism d80% (arcsec) */
  double	psfd80ast45;		/* 45 deg. astigmatism d80% (arcsec) */
  double	psfd80tri00;		/* 0 deg. triangular d80% (arcsec) */
  double	psfd80tri30;		/* 30 deg. triangular d80% (arcsec) */
  double	psfd80qua00;		/* 0 deg. quadratic d80% (arcsec) */
  double	psfd80qua22;		/* 22.5 deg. quadratic d80% (arcsec) */
  tracktypenum	psftracktype;		/* Tracking type */
  double	psftrackmaj;		/* Maximum RMS tracking error (") */
  double	psftrackmin;		/* Minimum RMS tracking error (") */
  double	psftrackang;		/* CC angle of maj. axis (deg) */
/*------ Aureole */
  int		aurange;		/* Maximum radius of the aurole */
/*------ Stellar field */
  double	scountdens;		/* Star nb per sq. deg. at mag.lim. */
  double	scountslope;		/* Diff. star count slope (dexp/mag) */
  double	maglim[2];		/* Brightest and faintest mag. allowed*/
  int		nmaglim;		/* Number of arguments */

/*------ Random generator */
  int		psfmotionseed;		/* Seed for PSF motion */
  int		starposseed;		/* Seed for star positions */
/* Multithreading */
  int		nthreads;		/* Number of active threads */
/*------ Grid parameters */
  int		grid_size;		/* Grid size (pixels) in GRID mode */
/* Misc */
  enum {QUIET, NORM, LOG, FULL} verbose_type;	/* display type */
  enum {NONE, ASTROMATIC}  xml_type;		/* XML output type */
  char		xml_name[MAXCHAR];		/* XML file name */
  char 		sdate_start[12];		/* SkyMaker start date */
  char		stime_start[12];		/* SkyMaker start time */
  char		sdate_end[12];			/* SkyMaker end date */
  char		stime_end[12];			/* SkyMaker end time */
  double	time_diff; 			/* Execution time */
  int		nobj;				/* Number of sources added */
        }	prefstruct;

prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);

#endif
