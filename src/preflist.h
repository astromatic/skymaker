/*
*				preflist.h
*
* Configuration keyword definitions.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2016 IAP/CNRS/UPMC
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
*	Last modified:		16/03/2016
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "key.h"

#ifndef _PREFS_H_
#include "prefs.h"
#endif

#ifdef	USE_THREADS
#define	THREADS_PREFMAX THREADS_NMAX
#else
#define	THREADS_PREFMAX 65535
#endif

int idummy;

pkeystruct key[] =
 {
  {"ARM_COUNT", P_INT, &prefs.psfnarms, 0, 16},
  {"ARM_THICKNESS", P_FLOAT, &prefs.psfarmw, 0,0, 0.0, 10000.0},
  {"ARM_POSANGLE", P_FLOAT, &prefs.psfarmang, 0,0, -360.0, 360.0},
  {"AST00_D80", P_FLOATLIST, prefs.psfd80ast00, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80ast00},
  {"AST45_D80", P_FLOATLIST, prefs.psfd80ast45, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80ast45},
  {"AST_CENTER", P_FLOATLIST, prefs.psfastc, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsfastc},
  {"AUREOLE_SB", P_FLOAT, &prefs.psfhalosb, 0,0, 12.0,50.0},
  {"AUREOLE_RADIUS", P_INT, &prefs.aurange, 0, 16384},
  {"BACK_MAG", P_FLOAT, &prefs.magback, 0,0, -100.0, 100.0},
  {"COMA_CENTER", P_FLOATLIST, prefs.psfcomac, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsfcomac},
  {"COMAX_D80", P_FLOATLIST, prefs.psfd80comax, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80comax},
  {"COMAY_D80", P_FLOATLIST, prefs.psfd80comay, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80comay},
  {"CORRELATION_LENGTH", P_FLOAT, &prefs.correlation_length, 0,0, 0.0, 1024.0},
  {"CORRELATION_TYPE", P_KEY, &prefs.correlation_type, 0,0, 0.0,0.0,
   {"NONE","NOISE", "ALL",""}},
  {"DEFOC_D80", P_FLOATLIST, prefs.psfd80defoc, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80defoc},
  {"DEFOC_CENTER", P_FLOATLIST, prefs.psfdefocc, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsfdefocc},
  {"EXPOSURE_TIME", P_FLOAT, &prefs.expotime, 0,0, 0.0, 1e18},
  {"GAIN", P_FLOAT, &prefs.gain, 0,0, 0.0,1e18},
  {"GRID_SIZE", P_INT, &prefs.grid_size, 3, 16384},
  {"IMAGE_HEADER", P_STRING, prefs.headname},
  {"IMAGE_NAME", P_STRING, prefs.filename},
  {"IMAGE_TYPE", P_KEY, &prefs.imatype, 0,0, 0.0,0.0,
   {"PUPIL_REAL","PUPIL_IMAGINARY","PUPIL_MODULUS","PUPIL_PHASE","PUPIL_MTF",
	"PSF_MTF","PSF_FULLRES","PSF_FINALRES","SKY_NONOISE","SKY",
	"GRID_NONOISE","GRID",""}},
  {"IMAGE_SIZE", P_INTLIST, prefs.imasize, 1, 100000, 0.0,0.0,
   {""}, 1, 2, &prefs.nimasize},
  {"LISTCOORD_TYPE", P_KEY, &prefs.listcoord_type, 0,0, 0.0,0.0,
   {"PIXEL","WORLD",""}},
  {"LISTMAG_LIMITS", P_FLOATLIST, prefs.listmag_limits, 0,0, -100.0, 100.0,
   {""}, 2,2, &prefs.nlistmag_limits},
  {"M1_DIAMETER", P_FLOAT, &prefs.psfdm1, 0,0, 0.0, 100.0},
  {"M2_DIAMETER", P_FLOAT, &prefs.psfdm2, 0,0, 0.0, 100.0},
  {"MAG_LIMITS", P_FLOATLIST, prefs.maglim, 0,0, -30.0, 50.0,
   {""}, 2,2, &prefs.nmaglim},
  {"MAG_ZEROPOINT", P_FLOAT, &prefs.magzero, 0,0, -100.0,100.0},
  {"MICROSCAN_NSTEP", P_INTLIST, prefs.mscan, 1,16, 0.0,0.0,
   {""}, 1,2, &prefs.nmscan},
  {"NTHREADS", P_INT, &prefs.nthreads, -THREADS_PREFMAX, THREADS_PREFMAX},
  {"PIXEL_SIZE", P_FLOATLIST, &prefs.pixscale, 0,0, 1e-12,1e12,
   {""}, 1, 2, &prefs.npixscale},
  {"PSFCENTER_TYPE", P_KEY, &prefs.psfcentertype, 0,0, 0.0,0.0,
	{"UPPERHALF", "LOWERHALF", "HALF", "CENTROID", "CENTROID_COMMON",
	"PEAK",""}},
  {"PSF_OVERSAMP", P_FLOAT, &prefs.psfoversamp, 0,0, 0.01, 100.0},
  {"PSF_MAPSIZE", P_INTLIST, prefs.psfsize, 1, 8192, 0.0,0.0,
   {""}, 1, 2, &prefs.npsfsize},
  {"PSF_NAME", P_STRING, prefs.psfname},
  {"PSF_TYPE", P_KEY, &prefs.psftype, 0,0, 0.0,0.0,
	{"INTERNAL","FILE",""}},
  {"QUA00_D80", P_FLOATLIST, prefs.psfd80qua00, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80qua00},
  {"QUA22_D80", P_FLOATLIST, prefs.psfd80qua22, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80qua22},
  {"QUA_CENTER", P_FLOATLIST, prefs.psfquac, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsfquac},
  {"READOUT_NOISE", P_FLOAT, &prefs.ron, 0,0, 0.0,1e9},
  {"SATUR_LEVEL", P_FLOAT, &prefs.satlev, 0,0, -1e18,1e18},
  {"SEED_MOTION", P_INT, &prefs.psfmotionseed, 0, 0x7fffffffL},
  {"SEED_STARPOS", P_INT, &prefs.starposseed, 0, 0x7fffffffL},
  {"SEEING_FWHM", P_FLOAT, &prefs.seeing, 0,0, 0.0, 1e12},
  {"SEEING_TYPE", P_KEY, &prefs.psfseeingtype, 0,0, 0.0,0.0,
	{"NONE","LONG_EXPOSURE","SHORT_EXPOSURE",""}},
  {"SPHER_D80", P_FLOATLIST, prefs.psfd80spher, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80spher},
  {"SPHER_CENTER", P_FLOATLIST, prefs.psfspherc, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsfspherc},
  {"STARCOUNT_SLOPE", P_FLOAT, &prefs.scountslope, 0,0, 1e-3, 10.0},
  {"STARCOUNT_ZP", P_FLOAT, &prefs.scountdens, 0,0, 0.0, 1e12},
  {"TRACKERROR_ANG", P_FLOAT, &prefs.psftrackang, 0,0, -360.0, 360.0},
  {"TRACKERROR_MAJ", P_FLOAT, &prefs.psftrackmaj, 0,0, 0.0, 100.0},
  {"TRACKERROR_MIN", P_FLOAT, &prefs.psftrackmin, 0,0, 0.0, 100.0},
  {"TRACKERROR_TYPE", P_KEY, &prefs.psftracktype, 0,0, 0.0, 0.0,
	{"NONE","DRIFT","JITTER",""}},
  {"TRI00_D80", P_FLOATLIST, prefs.psfd80tri00, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80tri00},
  {"TRI30_D80", P_FLOATLIST, prefs.psfd80tri30, 0,0, -1e3, 1e3,
   {""}, 1, 3, &prefs.npsfd80tri30},
  {"TRI_CENTER", P_FLOATLIST, prefs.psftric, 0,0, -1e18, 1e18,
   {""}, 2, 2, &prefs.npsftric},
  {"VERBOSE_TYPE", P_KEY, &prefs.verbose_type, 0,0, 0.0,0.0,
   {"QUIET","NORMAL","FULL",""}},
  {"WAVELENGTH", P_FLOAT, &prefs.lambdaeq, 0,0, 1e-12, 1e12},
  {"WELL_CAPACITY", P_FLOAT, &prefs.wellcap, 0,0, 0.0,1e18},
  {""}
 };

char			keylist[sizeof(key)/sizeof(pkeystruct)][32];
static const char	notokstr[] = {" \t=,;\n\r\""};

char *default_prefs[] =
 {
"# Default configuration file for " BANNER " " MYVERSION,
"# EB " DATE,
"#",
" ",
"#--------------------------------- Image -------------------------------------",
" ",
"IMAGE_NAME         sky.fits     # Name of the output frame",
"IMAGE_SIZE         1024         # Width,[height] of the output frame",
"IMAGE_TYPE         SKY          # PUPIL_REAL,PUPIL_IMAGINARY,PUPIL_MODULUS,",
"                                # PUPIL_PHASE,PUPIL_MTF,PSF_MTF,PSF_FULLRES,",
"                                # PSF_FINALRES,SKY_NONOISE,SKY,GRID",
"                                # or GRID_NONOISE",
"*GRID_SIZE          64           # Distance between objects in GRID mode",
"IMAGE_HEADER       INTERNAL     # File name or INTERNAL",
"LISTCOORD_TYPE     PIXEL        # Coordinates in input lists: PIXEL or WORLD",
"LISTMAG_LIMITS     -99.0,99.0   # Magnitude range restriction for input lists",
" ",
"#-------------------------------- Detector -----------------------------------",
" ",
"GAIN               1.0          # gain (e-/ADU)",
"WELL_CAPACITY      0            # full well capacity in e- (0 = infinite)",
"SATUR_LEVEL        65535        # saturation level (ADU)",
"READOUT_NOISE      1.0          # read-out noise (e-)",
"EXPOSURE_TIME      300.0        # total exposure time (s)",
"MAG_ZEROPOINT      26.0         # magnitude zero-point (\"ADU per second\")",
" ",
"#-------------------------------- Sampling -----------------------------------",
" ",
"PIXEL_SIZE         0.200        # pixel size in arcsec.",
"MICROSCAN_NSTEP    1            # number of microscanning steps (1=no mscan)",
"CORRELATION_TYPE   NONE         # pixel correlation (NONE, NOISE or ALL)",
"CORRELATION_LENGTH 2            # autocorrelation length in pixels",
" ",
"#---------------------------------- PSF --------------------------------------",
" ",
"PSF_TYPE           INTERNAL     # INTERNAL or FILE",
"PSF_NAME           psf.fits     # Name of the FITS image containing the PSF",
"*PSFCENTER_TYPE     UPPERHALF    # UPPERHALF, LOWERHALF, HALF, CENTROID,",
"*                                # CENTROID_COMMON or PEAK",
"SEEING_TYPE       LONG_EXPOSURE # (NONE, LONG_EXPOSURE or SHORT_EXPOSURE)",
"SEEING_FWHM        0.7          # FWHM of seeing in arcsec (incl. motion)",
"AUREOLE_RADIUS     200          # Range covered by aureole (pix) 0=no aureole",
"AUREOLE_SB         16.0         # SB (mag/arcsec2) at 1' from a 0-mag star",
"PSF_OVERSAMP       5            # Oversampling factor / final resolution",
"PSF_MAPSIZE        1024         # PSF mask size (pixels): must be a power of 2",
"TRACKERROR_TYPE    NONE         # Tracking error model: NONE, DRIFT or JITTER",
"TRACKERROR_MAJ     0.0          # Tracking RMS error (major axis) (in arcsec)",
"TRACKERROR_MIN     0.0          # Tracking RMS error (minor axis) (in arcsec)",
"TRACKERROR_ANG     0.0          # Tracking angle (in deg, CC/horizontal)",
" ",
"#----------------------------- Pupil features --------------------------------",
" ",
"M1_DIAMETER        3.6          # Diameter of the primary mirror (in meters)",
"M2_DIAMETER        1.0          # Obstruction diam. from the 2nd mirror in m.",
"ARM_COUNT          4            # Number of spider arms (0 = none)",
"ARM_THICKNESS      20.0         # Thickness of the spider arms (in mm)",
"ARM_POSANGLE       0.0          # Position angle of the spider pattern / AXIS1",
"DEFOC_D80          0.0          # Defocusing d80% diameter (arcsec)",
"*DEFOC_CENTER       0.5,0.5      # Relative center of PSF focus variations",
"SPHER_D80          0.0          # Spherical d80% diameter (arcsec)",
"*SPHER_CENTER       0.5,0.5      # Center of PSF spherical aber. variations",
"COMAX_D80          0.0          # Coma along X d80% diameter (arcsec)",
"COMAY_D80          0.0          # Coma along Y d80% diameter (arcsec)",
"*COMA_CENTER        0.5,0.5      # Center of PSF coma variations",
"AST00_D80          0.0          # 0 deg. astigmatism d80% diameter (arcsec)",
"AST45_D80          0.0          # 45 deg. astigmatism d80% diameter (arcsec)",
"*AST_CENTER         0.5,0.5      # Center of PSF astigmatism variations",
"TRI00_D80          0.0          # 0 deg. triangular d80% diameter (arcsec)",
"TRI30_D80          0.0          # 30 deg. triangular d80% diameter (arcsec)",
"*TRI_CENTER         0.5,0.5      # Center of PSF triangular aber. variations",
"QUA00_D80          0.0          # 0 deg. quadratic d80% diameter (arcsec)",
"QUA22_D80          0.0          # 22.5 deg. quadratic d80% diameter (arcsec)",
"*QUA_CENTER         0.5,0.5      # Center of PSF quad. aber. variations",
" ",
"#--------------------------------- Signal ------------------------------------",
" ",
"WAVELENGTH         0.8          # average wavelength analysed (microns)",
"BACK_MAG           20.0         # background surface brightness (mag/arcsec2)",
" ",
"#------------------------------ Stellar field --------------------------------",
" ",
"STARCOUNT_ZP       3e4          # nb of stars /deg2 brighter than MAG_LIMITS",
"STARCOUNT_SLOPE    0.2          # slope of differential star counts (dexp/mag)",
"MAG_LIMITS         17.0,26.0    # stellar magnitude range",
" ",
"#------------------------------ Random Seeds ---------------------------------",
" ",
"SEED_MOTION        0            # rand. seed for PSF turbulent motion (0=time)",
"SEED_STARPOS       0            # random seed for star positions (0=time)",
" ",
"#----------------------------- Miscellaneous ---------------------------------",
" ",
"VERBOSE_TYPE       NORMAL       # QUIET, NORMAL or FULL",
#ifdef USE_THREADS
"NTHREADS           0            # Number of simultaneous threads for",
"                                # the SMP version of " BANNER,
#else
"NTHREADS           1            # 1 single thread",
#endif
""};
