# configure.in for SkyMaker
# (C) E.Bertin 1995-2009
# Process this file with autoconf to produce a configure script.
# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT([skymaker],[3.3.2],[bertin@iap.fr])
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
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)

# Display pakage and version number
AC_MSG_RESULT([*********** Configuring:  $PACKAGE $VERSION ($date) **********])

# Initialize the list of compilers to consider
cclist="cc gcc"
# Backup and reset the input CFLAGS
mycflags="$CFLAGS"
CFLAGS=""

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([for Linux Intel C compiler mode])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING(--enable-icc,Enable special mode for compilation with the Intel compiler \
(off by default))],
        use_icc="yes"
        cclist="icc $cclist"
	AC_MSG_RESULT([yes]),
        use_icc="no"
	AC_MSG_RESULT([no]))

# Checks for programs.
# GCC is chosen last because it is likely to yield less optimized code
AC_LANG(C)
AC_PROG_CC([$cclist])
# C Compiler: Check that it is ANSI C and POSIX-compliant
AM_PROG_CC_STDC
AC_ISC_POSIX
ACX_PROG_CC_OPTIM
AC_DISABLE_STATIC
#AC_DISABLE_SHARED
AC_PROG_LIBTOOL
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h sys/mman.h])
echo

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv memcpy memmove memset munmap strstr sysconf sincos])

# Provide special options for FFTW
AC_ARG_WITH(fftw,
        [AC_HELP_STRING([--with-fftw=<FFTW library path>],
        [Provide an alternative path to the FFTW library])],
    fftw_libdir=$withval,
    fftw_libdir=""
    )
AC_ARG_WITH(fftw-incdir,
        [AC_HELP_STRING([--with-fftw-incdir=<FFTW include dir>],
        [Provide an alternative path to the FFTW include directory])],
    fftw_incdir=$withval,
    fftw_incdir=""
    )

# Set flags for multithreading
n_pthreads=16
AC_ARG_ENABLE(threads,
        [AC_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
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

# Provide special option for debugging
AC_MSG_CHECKING([for debugger mode])
AC_ARG_ENABLE(debug,
	[AS_HELP_STRING(--enable-debug,Enable special mode for debugging \
(off by default))],
	use_debug="yes"
	AC_MSG_RESULT([yes]),
	use_debug="no"
	AC_MSG_RESULT([no]))

# Provide special option for gprof profiling
AC_MSG_CHECKING([for gprof profiler mode])
AC_ARG_ENABLE(gprof,
	[AC_HELP_STRING([--enable-gprof],
	[Enable special mode for compilation with the gprof profiler \
(off by default)])],
	use_gprof="yes"
	AC_MSG_RESULT([yes]),
	use_gprof="no"
	AC_MSG_RESULT([no]))

# Enable static linking
AC_MSG_CHECKING([for static linking])
AC_ARG_ENABLE(static,
	[AC_HELP_STRING([--enable-static],
	[Enable static linking (off by default)])],
	use_static="$enableval"
	AC_MSG_RESULT([yes]),
	use_static="no"
	AC_MSG_RESULT([no]))

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [CFLAGS="$CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$PTHREAD_LIBS $LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

################ handle the FFTW library (Fourier transforms) ################
ACX_FFTW($fftw_libdir,$fftw_incdir,$use_pthreads,yes,
	[use_fftw=yes],[use_fftw=no])
if test "$use_fftw" = "yes"; then
  LIBS="$FFTW_LIBS $LIBS"
else
  AC_MSG_ERROR([$FFTW_ERROR Exiting.])
fi

# Debugging
if test "$use_debug" = "yes"; then
  CFLAGS="-g $CFLAGS"
fi

# Link with gprof option
if test "$use_gprof" = "yes"; then
  if test "$use_icc" = "yes"; then
    CFLAGS="$CFLAGS -pq"
  else
    CFLAGS="$CFLAGS -pg"
  fi
  use_static="no"
fi

# Static linking option
if test "$use_static" = "yes"; then
  LDFLAGS="-static -shared-libgcc $LDFLAGS"
fi

# Override automatic CFLAGS and LDFLAGS with those of user
#if test -n "$mycflags"; then
#CFLAGS="$mycflags"
#fi

# Display compiler and linker flags
AC_MSG_RESULT([I will compile using:  $CC $CFLAGS])
AC_MSG_RESULT([I will link using:     $CC $LDFLAGS $LIBS])
AC_CONFIG_FILES([Makefile man/Makefile src/Makefile src/fits/Makefile \
skymaker.spec man/sky.1])
AC_OUTPUT

