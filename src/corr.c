/*
*				corr.c
*
* Generate correlation square root function.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SkyMaker
*
*	Copyright:		(C) 1998-2017 IAP/CNRS/UPMC
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
*	Last modified:		18/04/2017
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#include FFTW_H

#ifdef HAVE_MKL
 #include MKL_H
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "prefs.h"
#include "corr.h"
#include "simul.h"

/****** corr_generate **************************************************
PROTO	void corr_generate(simstruct *sim, corrinterpenum interp_type,
		float scale)
PURPOSE	Generate average correlation square root function for the provided
	interpolation function.
INPUT	Pointer to the simulation,
	interpolant type,
	relative pixel scale (>1 if the output grid oversamples the input grid)
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	18/04/2017
 ***/
void	corr_generate(simstruct *sim, corrinterpenum corrinterp_type,
		float scale)

  {
   int		i;

  sim->corrsize[0] = sim->corrsize[1] = 7;
  QMALLOC16(sim->corr, PIXTYPE, sim->corrsize[0] * sim->corrsize[1]);
  for (i = sim->corrsize[0] * sim->corrsize[1]; i--;)
    sim->corr[i] = 1.0 / (sim->corrsize[0] * sim->corrsize[1]);

  return;
  }


