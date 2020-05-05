#
#				Makefile.am
#
# src Makefile. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SkyMaker
#
#	Copyright:		(C) 2003-2020 IAP/CNRS/SorbonneU
#
#	License:		GNU General Public License
#
#	SkyMaker is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SkyMaker is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SkyMaker. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		05/05/2020
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBDIRS		= fits wcs

bin_PROGRAMS	= sky

sky_SOURCES	= alterimage.c corr.c fft.c fitswcs.c galaxies.c image.c \
		  imaout.c list.c main.c makeit.c noise.c poly.c prefs.c psf.c \
		  random.c rasters.c simul.c stars.c weight.c \
		  alterimage.h corr.h define.h fft.h fitswcs.h galaxies.h \
		  globals.h image.h imaout.h key.h list.h noise.h poly.h \
		  preflist.h prefs.h psf.h random.h rasters.h simul.h stars.h \
		  threads.h \
		  types.h wcscelsys.h weight.h

sky_LDADD	= $(srcdir)/fits/libfits.a \
		  $(srcdir)/wcs/libwcs_c.a

DATE=`date +"%Y-%m-%d"`

