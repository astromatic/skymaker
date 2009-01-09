# Main Makefile for SkyMaker
# Copyright (C) 1995-2008 Emmanuel Bertin.
AUTOMAKE_OPTIONS	= foreign no-dependencies
SUBDIRS			= man src
EXTRA_DIST		= config doc AUTHORS BUGS ChangeLog COPYING \
			  HISTORY INSTALL README THANKS \
			  acx_fftw.m4 acx_pthread.m4 acx_prog_cc_optim.m4 
RPM_ROOTDIR		= `rpmbuild --nobuild -E %_topdir`
RPM_SRCDIR		= $(RPM_ROOTDIR)/SOURCES
dist-hook:
	rm -rf `find $(distdir) -name .svn`

rpm:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	USE_ICC="1" rpmbuild -ba --clean --nodeps $(PACKAGE_NAME).spec

rpm-opteron:	dist
	echo @PLPLOT_CFLAGS@
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -march=opteron" TPXLIBS="-lpthread -lfftw3f -lm" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXFLAGS="-O3 -g -pthread -funroll-loops -fomit-frame-pointer -Wall -march=opteron" TPXLIBS="-lpthread -lfftw3f_threads -lfftw3f -lm" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME)-mp.spec

rpm-athlon:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lpthread -lfftw3f -lm" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="gcc" TPXFLAGS="-O3 -g -pthread -funroll-loops -fomit-frame-pointer -Wall -m32 -march=i686 -msse -mfpmath=sse -mtune=athlon" TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lpthread -lfftw3f_threads -lfftw3f -lm" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME)-mp.spec

rpm-opteron-icc:	dist
	echo @PLPLOT_CFLAGS@
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -axWP -ip -unroll" TPXLIBS="-lpthread -lfftw3f -lm" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME).spec
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc" TPXFLAGS="-O3 -axWP -ip -unroll -pthread" TPXLIBS="-lpthread -lfftw3f_threads -lfftw3f -lm" rpmbuild -ba --target=x86_64 --clean $(PACKAGE_NAME)-mp.spec

rpm-athlon-icc:	dist
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -ip -unroll" TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lpthread -lfftw3f -lm" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME).spec
	cp -f $(PACKAGE_NAME)-$(PACKAGE_VERSION).tar.gz $(RPM_SRCDIR)
	TPXCC="icc32" TPXFLAGS="-O3 -ip -unroll -pthread" TPXLIBS="-L/usr/local/lib/ -L/usr/lib/ -lpthread -lfftw3f_threads -lfftw3f -lm" rpmbuild -ba --target=i686 --clean $(PACKAGE_NAME)-mp.spec
