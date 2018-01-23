/*
*				imaout.c
*
* Manage data input/output.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2018 IAP/CNRS/UPMC
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
*	Last modified:		23/01/2018
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <pwd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "imaout.h"
#include "key.h"
#include "prefs.h"
#include "simul.h"

/*--------------------------- structure definitions -------------------------*/
headkeystruct	headparam[] = {
  {"OBJECT", "Input object list ",
	prefs.inlistname, H_STRING, T_STRING},
  {"CONFFILE", "Configuration filename",
	prefsname, H_STRING, T_STRING},
  {"IMATYPE", "Image type ",
	&prefs.imatype, H_STRING, T_STRING, "IMAGE_TYPE"},
  {"GAIN", "Detector gain or conversion factor (e-/ADU)",
	&prefs.gain, H_FLOAT, T_DOUBLE},
  {"WELLCAPA", "Full well capacity (e-)",
	&prefs.wellcap, H_EXPO, T_DOUBLE},
  {"SATURATE", "Saturation level (ADU)",
	&prefs.satlev, H_FLOAT, T_DOUBLE},
  {"RON", "Read-out noise (ADU)",
	&prefs.ron, H_FLOAT, T_DOUBLE},
  {"EXPTIME", "Exposure time (s)",
	&prefs.expotime, H_FLOAT, T_DOUBLE},
  {"MAGZERO", "Magnitude Zero-point (mag/s)",
	&prefs.magzero, H_FLOAT, T_DOUBLE},
  {"PIXSIZE", "Pixel size (arcsec)",
	&prefs.pixscale[0], H_FLOAT, T_DOUBLE},
  {"NMSCAN", "Number of microscanning steps",
	&prefs.mscan[0], H_INT, T_LONG},
  {"PSFTYPE", "PSF type",
	&prefs.psftype, H_STRING, T_STRING, "PSF_TYPE"},
  {"PSFNAME", "PSF filename",
	prefs.psfname, H_STRING, T_STRING},
  {"SEEING", "FWHM of atmospheric seeing (arcsec)",
	&prefs.seeing, H_FLOAT, T_DOUBLE},
  {"SNGTYPE", "Seeing type",
	&prefs.psfseeingtype, H_STRING, T_STRING, "SEEING_TYPE"},
  {"AURESB", "Aureole SB at 1' (mag/arcsec2)",
	&prefs.psfhalosb, H_FLOAT, T_DOUBLE},
  {"AURERAD", "Aureole radius",
	&prefs.aurange, H_FLOAT, T_DOUBLE},
  {"PSFSAMP", "Oversampling factor for PSF",
	&prefs.psfoversamp, H_FLOAT, T_DOUBLE},
  {"PSFSIZE", "PSF mask size",
	&prefs.psfsize[0], H_INT, T_LONG},
  {"TRCKTYPE", "Tracking error type",
	&prefs.psftracktype, H_STRING, T_STRING, "TRACKERROR_TYPE"},
  {"TRACKMAJ", "Tracking RMS error: major axis (arcsec)",
	&prefs.psftrackmaj, H_FLOAT, T_DOUBLE},
  {"TRACKMIN", "Tracking RMS error: minor axis (arcsec)",
	&prefs.psftrackmin, H_FLOAT, T_DOUBLE},
  {"TRACKANG", "Tracking error: pos. angle (deg, CCW/horiz.)",
	&prefs.psftrackang, H_FLOAT, T_DOUBLE},
  {"DIAMM1", "M1 diameter (m)",
	&prefs.psfdm1, H_FLOAT, T_DOUBLE},
  {"DIAMM2", "M2 obstruction diameter (m)",
	&prefs.psfdm2, H_FLOAT, T_DOUBLE},
  {"ARMCOUNT", "Number of spider arms",
	&prefs.psfnarms, H_INT, T_LONG},
  {"ARMTHICK", "Thickness of spider arms (mm)",
	&prefs.psfarmw, H_FLOAT, T_DOUBLE},
  {"ARMANG", "Position angle of arms (deg, CCW/horiz.)",
	&prefs.psfarmang, H_FLOAT, T_DOUBLE},
  {"PSFDEFOC", "D80 of defocus component (arcsec)",
	&prefs.psfd80defoc, H_FLOAT, T_DOUBLE},
  {"PSFSPHER", "D80 of spherical aberration (arcsec)",
	&prefs.psfd80spher, H_FLOAT, T_DOUBLE},
  {"PSFCOMAX", "D80 of coma along x (arcsec)",
	&prefs.psfd80comax, H_FLOAT, T_DOUBLE},
  {"PSFCOMAY", "D80 of coma along y (arcsec)",
	&prefs.psfd80comay, H_FLOAT, T_DOUBLE},
  {"PSFAST00", "D80 of astigmatism at 0 deg. (arcsec)",
	&prefs.psfd80ast00, H_FLOAT, T_DOUBLE},
  {"PSFAST45", "D80 of astigmatism at 45 deg. (arcsec)",
	&prefs.psfd80ast45, H_FLOAT, T_DOUBLE},
  {"PSFTRI00", "D80 of triang. aber. at 0 deg. (arcsec)",
	&prefs.psfd80tri00, H_FLOAT, T_DOUBLE},
  {"PSFTRI30", "D80 of triang. aber. at30 deg. (arcsec)",
	&prefs.psfd80tri30, H_FLOAT, T_DOUBLE},
  {"PSFQUA00", "D80 of quad. aber. at 0 deg. (arcsec)",
	&prefs.psfd80qua00, H_FLOAT, T_DOUBLE},
  {"PSFQUA22", "D80 of quad. aber. at 22.5 deg. (arcsec)",
	&prefs.psfd80qua22, H_FLOAT, T_DOUBLE},
  {"WAVELEN", "Equivalent wavelength (microns)",
	&prefs.lambdaeq, H_EXPO, T_DOUBLE},
  {"BACKMAG", "Sky background SB (mag/arcsec2)",
	&prefs.magback, H_FLOAT, T_DOUBLE},
  {"SCOUNTZP", "Total number of added stars",
	&prefs.scountdens, H_EXPO, T_DOUBLE},
  {"SCOUNTSL", "Slope of star counts (dlogN/dm)",
	&prefs.scountslope, H_FLOAT, T_DOUBLE},
  {"SMAGMIN", "Minimum mag. for simulated stars (mag)",
	&prefs.maglim[0], H_FLOAT, T_DOUBLE},
  {"SMAGMAX", "Maximum mag. for simulated stars (mag)",
	&prefs.maglim[1], H_FLOAT, T_DOUBLE},
  {"SEEDPMOT", "Seed for randomized PSF motion (0=rand)",
	&prefs.psfmotionseed, H_INT, T_LONG},
  {"SEEDSPOS", "Seed for randomized star pos (0=rand)",
	&prefs.starposseed, H_INT, T_LONG},
  {""}
  };


extern pkeystruct	key[];
extern char		keylist[][32];
extern time_t		thetime;


/****** imaout_readaschead ****************************************************
PROTO	int imaout_readaschead(char *filename, int frameno, tabstruct *tab)
PURPOSE	Read an ASCII header file and update the current field's tab
INPUT	Name of the ASCII file,
	Frame number (if extensions),
	Tab structure.
OUTPUT	RETURN_OK if the file was found, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/03/2011
 ***/
int	imaout_readaschead(char *filename, int frameno, tabstruct *tab)
  {
   char		keyword[88],data[88],comment[88], str[88];
   FILE		*file;
   h_type	htype;
   t_type	ttype;
   int		i;

  if ((file=fopen(filename, "r")))
    {
/*- Skip previous ENDs in multi-FITS extension headers */
    for (i=frameno; i--;)
      while (fgets(str, MAXCHAR, file)
		&& strncmp(str,"END ",4)
		&& strncmp(str,"END\n",4));
    memset(str, ' ', 80);
    while (fgets(str, 81, file) && strncmp(str,"END ",4)
		&& strncmp(str,"END\n",4))
      {
      fitspick(str, keyword, data, &htype, &ttype, comment);
/*---- Block critical keywords */
      if (!wstrncmp(keyword, "SIMPLE  ", 8)
	||!wstrncmp(keyword, "BITPIX  ", 8)
	||!wstrncmp(keyword, "NAXIS   ", 8)
	||!wstrncmp(keyword, "NAXIS1  ", 8)
	||!wstrncmp(keyword, "NAXIS2  ", 8)
	||!wstrncmp(keyword, "NAXIS3  ", 8)
	||!wstrncmp(keyword, "NAXIS4  ", 8)
	||!wstrncmp(keyword, "NAXIS5  ", 8)
	||!wstrncmp(keyword, "BSCALE  ", 8)
	||!wstrncmp(keyword, "BZERO   ", 8))
        continue;
      addkeywordto_head(tab, keyword, comment);
      fitswrite(tab->headbuf, keyword, data, htype, ttype);
      memset(str, ' ', 80);
      }
    fclose(file);
/*-- Update the tab data */
    readbasic_head(tab);
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }


/******* imaout_inithead ******************************************************
PROTO	catstruct *imaout_inithead(simstruct *sim)
PURPOSE	Initialize the header of the simulated image and the FITScat structure.
INPUT	Pointer to the simulation.
OUTPUT	Pointer to the initialized FITScat structure.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/05/2012
 ***/
catstruct	*imaout_inithead(simstruct *sim)

  {
   catstruct		*cat;
   tabstruct		*tab;
   headkeystruct	*headkey;
   void			*ptr;
   char			str[MAXCHAR],
			*pstr;
   double		dval;
   int			hflag;

  cat = new_cat(1);
  init_cat(cat);
  tab = cat->tab;
  tab->bitpix = BP_FLOAT;
  tab->bytepix = 4;
  tab->naxis = sim->imasize[2]>1? (sim->imasize[3]>1? 4 : 3): 2;
  QCALLOC(tab->naxisn, int, 4);
  tab->naxisn[0] = sim->imasize[0];
  tab->naxisn[1] = sim->imasize[1];
  tab->naxisn[2] = sim->imasize[2];
  tab->naxisn[3] = sim->imasize[3];
  update_head(tab);

  addkeywordto_head(tab, "COMMENT ",
	"Created by " BANNER " version " MYVERSION " (" DATE ")");
/* Update optional FITS keywords */
  dval = 1.0;
  addkeywordto_head(tab, "BSCALE  ", "True value = BSCALE*PIXEL+BZERO");
  fitswrite(tab->headbuf, "BSCALE  ",&dval, H_FLOAT, T_DOUBLE);
  dval = 0.0;
  addkeywordto_head(tab, "BZERO   ", "True value = BSCALE*PIXEL+BZERO");
  fitswrite(tab->headbuf, "BZERO   ",&dval, H_FLOAT, T_DOUBLE);
  hflag = 1;

/* See if the user supplies a header or not */
  if (cistrcmp(sim->headname, "INTERNAL", 0))
    {
    if (imaout_readaschead(sim->headname, 0, tab)==RETURN_OK)
      hflag = 0;
    else
      warning(sim->headname,
	" could not be opened\n> I shall therefore write a minimum header.");
    }

/* If no external header is used, add SkyMaker custom parameters */
  if (hflag)
    {
/*--- User name */
#ifdef HAVE_GETENV
  if (!(pstr=getenv("USERNAME")))       /* Cygwin,... */
    if ((pstr=getenv("LOGNAME")))       /* Linux,... */
      {
      addkeywordto_head(tab, "AUTHOR  ", "Whom it comes from");
      fitswrite(tab->headbuf, "AUTHOR  ",pstr, H_STRING, T_STRING);
      }
#endif
/*--- Host name */
    if (!gethostname(str, 80))
      {
      addkeywordto_head(tab, "ORIGIN  ", "Where it comes from");
      fitswrite(tab->headbuf, "ORIGIN  ", str, H_STRING, T_STRING);
      }
/*--- Obs dates */
    addkeywordto_head(tab, "DATE-OBS", "When it was started");
    fitswrite(tab->headbuf, "DATE-OBS", prefs.stime_start, H_STRING, T_STRING);
/*--- Reference */
    addkeywordto_head(tab, "REFERENC", "Software web page");
    fitswrite(tab->headbuf, "REFERENC", "http://astromatic.net/software/skymaker",
		H_STRING, T_STRING);

    for (headkey=headparam;*(headkey->name); headkey++)
      {
      addkeywordto_head(tab, headkey->name, headkey->comment);
/*---- The case for keywords */
      ptr = (*headkey->prefsname)?
          key[findkeys(headkey->prefsname, keylist, FIND_STRICT)].keylist[
						*((int *)headkey->ptr)]
        : headkey->ptr;
      fitswrite(tab->headbuf, headkey->name, ptr,headkey->htype,headkey->ttype);
      }
    }

  if (sim->wcsflag)
/*-- Initialize WCS structure */
    sim->wcs = read_wcs(tab);

  return cat;
  }


/******* imaout_write *********************************************************
PROTO	void imaout_write(simstruct *sim)
PURPOSE	Write the simulated image to a FITS file.
INPUT	Pointer to the simulation.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/01/2018
 ***/
void	imaout_write(simstruct *sim) {

   tabstruct	*tab;

  tab = sim->cat->tab;
  tab->bodybuf = (void *)sim->image;
  tab->tabsize = (KINGSIZE_T)tab->naxisn[0]*tab->naxisn[1] * tab->bytepix;
  if (tab->naxis > 2) {
    tab->tabsize *= tab->naxisn[2];
    if (tab->naxis > 3)
      tab->tabsize *= tab->naxisn[4];
  }

  save_cat(sim->cat, sim->filename);
  sim->cat->tab->bodybuf = NULL;

  return;
}

