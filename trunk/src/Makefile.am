# Program Makefile for SkyMaker
# Copyright (C) 1995-2010 Emmanuel Bertin.
SUBDIRS		= fits

bin_PROGRAMS	= sky

sky_SOURCES	= alterimage.c fft.c galaxies.c image.c imaout.c list.c main.c \
		  makeit.c poly.c prefs.c psf.c random.c simul.c stars.c \
		  alterimage.h define.h fft.h galaxies.h globals.h image.h \
		  imaout.h key.h list.h poly.h preflist.h prefs.h psf.h \
		  random.h simul.h stars.h threads.h types.h

sky_LDADD	= $(top_builddir)/src/fits/libfits.a

DATE=`date +"%Y-%m-%d"`
