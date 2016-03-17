#
#				configure.ac
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SkyMaker
#
#	Copyright:		(C) 2003-2016 IAP/CNRS/UPMC
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
#	along with SkyMaker.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		17/03/2016
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Process this file with autoconf to produce a configure script.

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT([skymaker],[3.12.1],[bertin@iap.fr])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "Emmanuel Bertin")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_fftw.m4)
sinclude(acx_mkl.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)

# Display pakage and version number
AC_MSG_RESULT([*********** Configuring:  $PACKAGE $VERSION ($date) **********])

AC_MSG_CHECKING([if compilation flags are set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags \
(off by default)])],
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether we should use the INTEL compiler])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (off by default)])],
	CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Provide special options for INTEL MKL
# We force the use of icc
AC_MSG_CHECKING([whether we should use INTEL's MKL])
AC_ARG_ENABLE(mkl,
	[AS_HELP_STRING([--enable-mkl],
	[Use INTEL's MKL for solvers and FFTs (off by default)])],
	enable_icc="yes"
	CC="icc"
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Checks for programs.
# GCC is chosen last because it is likely to yield less optimized code
AC_LANG(C)

AC_SEARCH_LIBS([strerror],[cposix])
if test "$enable_auto_flags" = "yes"; then
  CFLAGS=""
  LDFLAGS =""
  ACX_PROG_CC_OPTIM
fi
LT_INIT
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h sys/mman.h unistd.h])
if test "$enable_icc" = "yes"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv memcpy memmove memset munmap sincos sincosf \
		strstr sysconf])
 
# Provide special options for FFTW
AC_ARG_WITH(fftw-libdir,
        [AS_HELP_STRING([--with-fftw-libdir=<FFTW library path>],
        [Alternative directory for the FFTW library])])
AC_ARG_WITH(fftw-incdir,
        [AS_HELP_STRING([--with-fftw-incdir=<FFTW include dir>],
        [Alternative directory for FFTW include files])])

# Set flags for multithreading
n_pthreads=1024
AC_ARG_ENABLE(threads,
        [AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 16 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (off by default)])]
	AC_MSG_RESULT([yes]),
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (off by default)])])

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$PTHREAD_LIBS $LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

#################### handle the INTEL MKL library (FFTW + LAPACKe) ###########
if test "$enable_mkl" = "yes"; then
  convlibs="${srcdir}/wcs/libwcs_c.a"
  ACX_MKL(,,$enable_best_link,$convlibs)
  if test "$MKL_WARN" != ""; then
    AC_MSG_WARN([$MKL_WARN])
  fi
  AM_CFLAGS="$AM_CFLAGS $MKL_CFLAGS "
  AM_LDFLAGS="$AM_LDFLAGS $MKL_LDFLAGS "
  LIBS="$LIBS $MKL_LIBS"
else
################ handle the FFTW library (Fourier transforms) ################
  ACX_FFTW($with_fftw_libdir,$with_fftw_incdir,$use_pthreads,yes,
	[LIBS="$FFTW_LIBS $LIBS"]
	if test "$FFTW_WARN" != ""; then
	  AC_MSG_WARN([$FFTW_WARN])
	fi,
	AC_MSG_ERROR([$FFTW_ERROR Exiting.]))
fi

# Compile with profiling option
if test "$enable_profiling" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  enable_best_link="no"
fi

# "Best" linking option
if test "$enable_best_link" = "yes"; then
  if test "$enable_icc" = "yes"; then
    AM_LDFLAGS="-static-intel -qopenmp-link static -shared-libgcc \
	-static-libtool-libs -no-intel-extensions  -avoid-version $AM_LDFLAGS"
  else
    AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
  fi
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
man/Makefile
src/Makefile
src/fits/Makefile
src/wcs/Makefile
skymaker.spec
man/sky.1])
AC_OUTPUT

