#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pkdgrav2, 2.2.10, bugs@pkdgrav2.org)
AM_INIT_AUTOMAKE(pkdgrav2,2.2.10)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET

AC_ARG_WITH(mpicc,
  [  --with-mpicc=mpicc      specify the name of the MPI C compiler],
  [MPICC="${withval}"],[AC_CHECK_PROG([MPICC],mpicc,mpicc,no)])
AC_SUBST(MPICC)

# Cooling
AC_ARG_ENABLE(cooling,
  [  --enable-cooling             enable cooling],
  [COOLING=yes
   AC_DEFINE([OOOLING], 1,
 	[Define if cooling module should be used.])])
AM_CONDITIONAL(COOLING, test "x$COOLING" = "xyes")

if test "x$COOLING" = "xyes" ; then
#AC_PROG_FC
#AC_FC_MAIN
#AC_FC_DUMMY_MAIN
##AC_FC_WRAPPERS
#FC_MODULE_WRAPPERS
AC_ARG_WITH(mpif90,
  [  --with-mpif90=mpif90    specify the name of the MPI f90 compiler],
  [MPIFC="${withval}"],[AC_CHECK_PROG([MPIFC],mpif90,mpif90,no)])
AC_SUBST(MPIFC)
else
AC_DEFINE([FC_MAIN],main,
	[C main function (varies if FORTRAN is being used).])
fi

# The floating point data type - This should be detected and --with-float added
AC_DEFINE([FLOAT],double,
	[This defines the base floating point type for particles.])
AC_DEFINE([FLOAT_MAXVAL],HUGE,
	[This defines the maximum value that a FLOAT can represent.])

# This adds processor specific timing to the code - default disabled
AC_ARG_ENABLE(instrument,
  [  --disable-instrument           Disable instrumentation code],
  [INSTRUMENT="$enableval"],[INSTRUMENT="yes"])
if test "x$INSTRUMENT" = "xyes" ; then
   AC_DEFINE([INSTRUMENT], 1,
 	[Instrument the code.])
fi

# This adds backtracing for crashed processes -- default enabled
AC_ARG_ENABLE(bt,
  [  --disable-bt           Disable backtracing code],
  [BT="$enableval"],[BT="yes"])
if test "x$BT" = "xyes" ; then
   AC_DEFINE([USE_BT], 1,
 	[Print a backtrace on crash.])
fi
AM_CONDITIONAL(USE_BT, test "x$USE_BT" = "xyes")

# This adds FFTW analysis - default disabled
AC_ARG_ENABLE(fftw,
  [  --enable-fftw                  Enable FFT support modules],
  [FFTW="$enableval"],[FFTW="no"])
if test "x$FFTW" = "xyes" ; then
   AC_DEFINE([MDL_FFTW], 1,
 	[Include support for FFTW.])
fi
AM_CONDITIONAL(FFTW, test "x$FFTW" = "xyes")

# Local expansions - default is enabled
AC_ARG_ENABLE(local-expansion,
  [  --disable-local-expansion      disable local expansion code],
  [LOCAL_EXPANSION="$enableval"],[LOCAL_EXPANSION="yes"])
# Local expansions... this is the preferred flag
AC_ARG_ENABLE(fmm,
  [  --disable-fmm                  disable local expansion code],
  [LOCAL_EXPANSION="$enableval"],[LOCAL_EXPANSION="yes"])
if test "x$LOCAL_EXPANSION" != "xno" ; then
   AC_DEFINE([LOCAL_EXPANSION], 1,
 	[Define if local expansions should be used.])
fi
AM_CONDITIONAL(LOCAL_EXPANSION, test "x$LOCAL_EXPANSION" != "xno")

# Electric fence
AC_ARG_WITH(efence,
  [  --with-efence=PATH      enable electric fence debugging],
  [EFENCE_DIR="${withval}"],[EFENCE_DIR=no])
AM_CONDITIONAL(USE_EFENCE, test "x$EFENCE_DIR" != "xno")
AM_CONDITIONAL(HAVE_EFENCE_PATH, test "x$EFENCE_DIR" != "xyes")
AC_SUBST([EFENCE_PATH],[$EFENCE_DIR])

# Python
AC_ARG_WITH(python,
  [  --with-python=VERSION   enable python scripting],
  [PYTHON_EXE="${withval}"],[PYTHON_EXE=no])
AM_CONDITIONAL(USE_PYTHON, test "x$PYTHON_EXE" != "xno")
if test "x$PYTHON_EXE" != "xno" ; then
AC_DEFINE([USE_PYTHON],1,
  [Defined if python support should be compiled into pkdgrav2])
if test "x$PYTHON_EXE" = "xyes" ; then
  PYTHON_EXE=python
fi
PYTHON_VER=`$PYTHON_EXE -c "import sys; import distutils.sysconfig as ds; sys.stdout.write(ds.get_python_version())"`
PYTHON_INC=`$PYTHON_EXE -c "import sys; import distutils.sysconfig as ds; sys.stdout.write(ds.get_python_inc())"`
PYTHON_LIB=`$PYTHON_EXE -c "import os; import sys; import distutils.sysconfig as ds; sys.stdout.write(os.path.join(ds.get_python_lib(plat_specific=1,standard_lib=1), 'config'))"`
AC_SUBST([PYTHON_VER],[$PYTHON_VER])
AC_SUBST([PYTHON_INC],[$PYTHON_INC])
AC_SUBST([PYTHON_LIB],[$PYTHON_LIB])
fi

# HDF5 I/O
AC_ARG_ENABLE(hdf5,
  [  --enable-hdf5                  enable support for I/O using HDF5 libraries],
  [USE_HDF5="$enableval"])
AC_ARG_WITH(hdf5,
  [  --with-hdf5=PATH        path to HDF5 installation],
  [CFLAGS="$CFLAGS -I${withval}/include"
   LDFLAGS="$LDFLAGS -L${withval}/lib64"
   if test -n "${need_dash_r}"; then
	LDFLAGS="$LDFLAGS -R${withval}/lib64"
   fi
   USE_HDF5=yes
   ])
if test "x$USE_HDF5" = "xyes" ; then
  AC_CHECK_LIB([hdf5],[H5Fopen],,[AC_MSG_ERROR([The HDF5 library cannot be located])])
elif test "x$USE_HDF5" != "xno" ; then
  USE_HDF5="yes"
  AC_CHECK_LIB([hdf5],[H5Fopen],,[USE_HDF5="no"])
fi
if test "x$USE_HDF5" = "xyes" ; then
   AC_CHECK_LIB([z], [inflate],,[AC_MSG_ERROR([The HDF5 library needs libz])])
   AC_CHECK_LIB([gpfs], [gpfs_fcntl])
   AC_DEFINE([USE_HDF5], 1,
 	[Define if HDF5 I/O should be available.])
   AC_DEFINE([H5_USE_16_API],1,
        [Use the HDF5 1.6 API instead of later versions])
fi

AC_ARG_WITH(mdl,
  [  --with-mdl=PATH         base path of MDL source],
  [if test -d "${withval}" ; then
     MDL_PATH="${withval}"
   else
     AC_MSG_ERROR([MDL does not exist at ${withval}])
   fi
  ],
  [for D in mdl2 mdl ../mdl2 ../mdl ; do
     if test -d $D ; then
       MDL_PATH="$D"
       break
     fi
   done
   if test -z "${MDL_PATH}" ; then
     AC_MSG_ERROR([Cannot locate MDL, try --with-mdl=PATH])
   fi])
AC_SUBST([MDL_PATH],[$MDL_PATH])

AC_ARG_WITH(grafic,
  [  --with-grafic=PATH      base path of GRAFIC IC library],
  [AC_SUBST([GRAFIC_PATH],[${withval}])
   USE_GRAFIC="$enableval"
   AC_DEFINE([USE_GRAFIC], 1,
 	[Define GRAFICIC IC should be available.])],
  [USE_GRAFIC=no])
dnl if test "$USE_GRAFIC" != "yes" -a -d grafic ; then
dnl   USE_GRAFIC=yes
dnl   GRAFIC_PATH=grafic
dnl   AC_CONFIG_SUBDIRS(grafic)
dnl fi
AM_CONDITIONAL(USE_GRAFIC,test "x$USE_GRAFIC" == "xyes")
			       
# MDL dedicated I/O
AC_ARG_ENABLE(mdl-io,
  [  --enable-mdl-io                enable dedicated I/O processors (MPI/pthread)],
  [USE_MDL_IO=yes
   AC_DEFINE([USE_MDL_IO], 1,
 	[Define if MDL style I/O should be used.])
   AC_DEFINE([MDL_MAX_IO_PROCS],128,[The maximum number of IO processors])])
AM_CONDITIONAL(USE_MDL_IO, test "x$USE_MDL_IO" = "xyes")
AM_CONDITIONAL(USE_HDF5, test "x$USE_HDF5" = "xyes")

# SISCI MDL
AC_ARG_WITH(sisci,
  [  --with-sisci[=/opt/DIS]   path to SISCI installation],
  [if test "x$enableval" = "xyes" ; then
    AC_SUBST([SISCI_PATH],/opt/DIS)
    USE_SISCI=yes
   elif test "x$enableval" = "xno" ; then
    USE_SISCI=no
   else
    AC_SUBST([SISCI_PATH],$withval)
    USE_SISCI=yes
   fi
   CFLAGS="$CFLAGS -I${withval}/include"
   LDFLAGS="$LDFLAGS -L${withval}/lib64"
   if test -n "${need_dash_r}"; then
	LDFLAGS="$LDFLAGS -R${withval}/lib64"
   fi
   ])
AM_CONDITIONAL(HAVE_SISCI, test "x$USE_SISCI" = "xyes")

# PNG Output
AC_ARG_ENABLE(png,
  [  --enable-png                   enable PNG outputs for MDL IO)],
  [USE_PNG=yes
   AC_DEFINE([USE_PNG], 1,
 	[Define if PNG output is supported.])])
AM_CONDITIONAL(USE_PNG, test "x$USE_PNG" = "xyes")
AC_ARG_WITH(gd,
  [  --with-gd=PATH          path to GD installation],
  [CFLAGS="$CFLAGS -I${withval}/include"
   LDFLAGS="$LDFLAGS -L${withval}/lib"
   if test -n "${need_dash_r}"; then
	LDFLAGS="$LDFLAGS -R${withval}/lib"
   fi
])
if test "x$USE_PNG" = "xyes" ; then
  AC_CHECK_LIB([png],[png_create_write_struct],,[AC_MSG_ERROR([The png library cannot be located])])
  AC_CHECK_LIB([gd],[gdImageCreateTrueColor],,[AC_MSG_ERROR([The gd library cannot be located])])
fi

# Lustre support
AC_ARG_ENABLE(lustre,
  [  --enable-lustre                enable support for lustre],
  [LUSTRE="${enableval}"])
if test "x$LUSTRE" != "xno" ; then
  AC_MSG_CHECKING([for lustre])
  lustre_ok=yes
  saved_LIBS="$LIBS"
  LIBS="$LIBS -llustreapi"
  AC_TRY_LINK([], [llapi_file_create()], [AC_DEFINE(USE_LUSTRE,1,[Define if lustre support should be included.])], [lustre_ok=no])
  LIBS="$saved_LIBS"
  AC_MSG_RESULT($lustre_ok)
  if test "x$LUSTRE" = "xyes" -a "x$lustre_ok" != "xyes" ; then
      AC_MSG_ERROR([Lustre support was requested, but it doesn't seem to be available.])
  fi
fi

AM_CONDITIONAL(LUSTRE, test "x$lustre_ok" = "xyes")

# Enable BSC Tracing
AC_ARG_ENABLE(bsc-trace,
  [  --enable-bsc-trace             enable tracing at BSC],
  [USE_BSC="$enableval"])
if test "x$USE_BSC" = "xyes" ; then
  AC_DEFINE([USE_BSC], 1,
    [Define if BSC tracing should be used.])
fi
AM_CONDITIONAL(USE_BSC, test "x$USE_BSC" = "xyes")

# Enable SIMD version (SSE,Altivec)
AC_ARG_ENABLE(simd,
  [  --disable-simd                 disable SIMD optimizations],
  [USE_SIMD="$enableval"])

USE_SIMD_MOMR="no"
AC_ARG_ENABLE(simd-momr,
  [  --enable-simd-momr             enable MOMR SIMD optimizations],
  [USE_SIMD_MOMR="$enableval"])
if test "x$USE_SIMD_MOMR" = "xyes" ; then
  USE_SIMD=yes
fi

USE_SIMD_PP="$USE_SIMD"
AC_ARG_ENABLE(simd-pp,
  [  --disable-simd-pp              disable PP SIMD optimizations],
  [USE_SIMD_PP="$enableval"])
if test "x$USE_SIMD_PP" != "xno" ; then
  USE_SIMD="$USE_SIMD_PP"
fi

# The floating point data type for moments
if test "x$USE_SIMD_MOMR" = "xyes" -o "x$USE_SIMD_LOCR" = "xyes" ; then
AC_DEFINE([MOMFLOAT],float,
	[This defines the base floating point type for moments.])
else
AC_DEFINE([MOMFLOAT],double,
	[This defines the base floating point type for moments.])
fi

# Change softening
CHANGESOFT=yes
AC_ARG_ENABLE(change-softening,
  [  --disable-change-softening     disable softening change],
  [if test "x$enableval" = "xno" ; then CHANGESOFT=no ; fi])
if test "x$CHANGESOFT" = "xyes" ; then
AC_DEFINE([CHANGESOFT], 1,
 	[Define if softening can be changed.])
fi
AM_CONDITIONAL(CHANGESOFT, test "x$CHANGESOFT" = "xyes")

# HSHRINK M4 Kernel
AC_ARG_ENABLE(hshrink,
  [  --enable-hshrink               enable HSHRINK M4 Kernel],
  [HSHRINK=yes
   AC_DEFINE([HSHRINK], 1,
 	[Define if SHRINK H4 kernel should be used.])])
AM_CONDITIONAL(HSHRINK, test "x$HSHRINK" = "xyes")

# Relaxation
AC_ARG_ENABLE(relaxation,
  [  --enable-relaxation            enable relaxation],
  [RELAXATION=yes
   AC_DEFINE([RELAXATION], 1,
 	[Define if RELAXATION should be used.])])
AM_CONDITIONAL(RELAXATION, test "x$RELAXATION" = "xyes")

# Mass-weighted softening
AC_ARG_ENABLE(soft-mass,
  [  --disable-soft-mass            Do not weight softening by mass],
  [if test "x$enableval" = "xyes" ; then
     SOFTMASS="yes"
   elif test "x$enableval" = "xno" ; then
     SOFTMASS="no"
   else
     SOFTMASS="yes"
   fi])
if test "x$SOFTMASS" = "xno" ; then
  AC_DEFINE([SOFTENING_NOT_MASS_WEIGHTED],1, 
 	[Define if softening should NOT be weighted by mass.])
fi

# Soft Linear
AC_ARG_ENABLE(soft-linear,
  [  --enable-soft-linear           enable soft linear],
  [SOFTLINEAR=yes
   AC_DEFINE([SOFTLINEAR], 1,
 	[Define if SOFTLINEAR should be used.])])
AM_CONDITIONAL(SOFTLINEAR, test "x$SOFTLINEAR" = "xyes")

# Relaxation
AC_ARG_ENABLE(soft-square,
  [  --enable-soft-square           enable soft square],
  [SOFTSQUARE=yes
   AC_DEFINE([SOFTSQUARE], 1,
 	[Define if SOFTSQUARE should be used.])])
AM_CONDITIONAL(SOFTSQUARE, test "x$SOFTSQUARE" = "xyes")

# PLANETS 
AC_ARG_ENABLE(planets,
  [  --enable-planets               enable planets],
  [PLANETS=yes
   AC_DEFINE([PLANETS], 1,
 	[Define if planets should be used.])])
AM_CONDITIONAL(PLANETS, test "x$PLANETS" = "xyes")

# HERMITE
AC_ARG_ENABLE(hermite,
  [  --enable-hermite               enable hermite],
  [HERMITE=yes
   AC_DEFINE([HERMITE], 1,
 	[Define if hermite should be used.])])
AM_CONDITIONAL(HERMITE, test "x$HERMITE" = "xyes")

# SYMBA
AC_ARG_ENABLE(symba,
  [  --enable-symba               enable symba],
  [SYMBA=yes
   AC_DEFINE([SYMBA], 1,
 	[Define if symba should be used.])])
AM_CONDITIONAL(SYMBA, test "x$SYMBA" = "xyes")

# PSD
AC_ARG_ENABLE(psd,
  [  --enable-psd                 enable phase-space options],
  [PSD=yes
   AC_DEFINE([USE_PSD], 1,
 	[Define if phase-space coordinates/functions should be available.])])
AM_CONDITIONAL(PSD, test "x$PSD" = "xyes")

# Check if pthread is available
AC_CHECK_LIB([pthread], [pthread_create],[HAVE_PTHREAD=yes],[HAVE_PTHREAD=no])
if test "x$HAVE_PTHREAD" = "xyes" ; then
  AC_CHECK_HEADER([pthread.h],[HAVE_PTHREAD=yes],[HAVE_PTHREAD=no])
fi
AM_CONDITIONAL(HAVE_PTHREAD,test "x$HAVE_PTHREAD" = "xyes")

# Check if mpi is available
if test "x$MPICC" != "xno" ; then
  HAVE_MPI=yes
#  AC_CHECK_LIB([mpi], [MPI_Init],[HAVE_MPI=yes],[HAVE_MPI=no])
#  if test "x$HAVE_MPI" = "xyes" ; then
#    AC_CHECK_HEADER([mpi.h],[HAVE_MPI=yes],[HAVE_MPI=no])
#  fi
fi
AM_CONDITIONAL(HAVE_MPI,test "x$HAVE_MPI" = "xyes")

AC_CHECK_LIB([m], [fabs])
# FIXME: Replace `main' with a function in `-lrt':
#AC_CHECK_LIB([rt], [main])

if test "x$USE_SIMD" != "xno" ; then
  AC_CHECK_FUNCS([posix_memalign],HAVE_POSIX_MEMALIGN=yes)

  AC_MSG_CHECKING([if the compiler accepts -march=native])
  oldCFLAGS="$CFLAGS"; CFLAGS="-march=native"
  AC_TRY_COMPILE([],[],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-march=native",
	AC_MSG_RESULT(no);)
  CFLAGS="$oldCFLAGS"
  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_CHECKING([if the compiler accepts -msse3])
    oldCFLAGS="$CFLAGS"; CFLAGS="-msse3"
    AC_TRY_COMPILE([#include <pmmintrin.h>],[__m128 r;],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-msse3",
	AC_MSG_RESULT(no);)
    CFLAGS="$oldCFLAGS"
  fi
  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_CHECKING([if the compiler accepts -msse2])
    oldCFLAGS="$CFLAGS"; CFLAGS="-msse2"
    AC_TRY_COMPILE([#include <emmintrin.h>],[__m128 r;],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-msse2",
	AC_MSG_RESULT(no);)
    CFLAGS="$oldCFLAGS"
  fi
  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_CHECKING([if the compiler accepts -msse])
    oldCFLAGS="$CFLAGS"; CFLAGS="-msse"
    AC_TRY_COMPILE([#include <emmintrin.h>],[__m128 r;],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-msse",
	AC_MSG_RESULT(no);)
    CFLAGS="$oldCFLAGS"
  fi

  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_CHECKING([if the compiler accepts -qaltivec])
    oldCFLAGS="$CFLAGS"; CFLAGS="-qarch=ppc970 -qaltivec"
    AC_TRY_COMPILE([#include <altivec.h>],[void foo(void) {vector float r;}],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-qarch=ppc970 -qaltivec",
	AC_MSG_RESULT(no);)
    CFLAGS="$oldCFLAGS"
  fi
  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_CHECKING([if the compiler accepts -maltivec])
    oldCFLAGS="$CFLAGS"; CFLAGS="-maltivec"
    AC_TRY_COMPILE([#include <altivec.h>],[vector float r;],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-maltivec",
	AC_MSG_RESULT(no);)
    CFLAGS="$oldCFLAGS"
  fi
  AC_MSG_CHECKING([for SSE2/3 support])
  oldCFLAGS="$CFLAGS"; CFLAGS="$SIMD_CFLAGS"
  AC_TRY_COMPILE([
#if defined(__SSE__) && defined(__SSE2__)
#include <xmmintrin.h>
#include <emmintrin.h>
#ifdef __SSE3__
#include <pmmintrin.h>
#endif
#else
#error 'NO SSE'
#endif
],[__m128 x;],
	AC_MSG_RESULT(yes);SIMD_OK=yes,
	AC_MSG_RESULT(no);)
  CFLAGS="$oldCFLAGS"
  if test "x$SIMD_OK" = "x" ; then
    oldCFLAGS="$CFLAGS"; CFLAGS="$SIMD_CFLAGS"
    AC_MSG_CHECKING([for Altivec support])
    AC_TRY_COMPILE([
#if defined(__ALTIVEC__)
#include <altivec.h>
#else
#error 'No Altivec'
#endif
],[vector float x;],
	AC_MSG_RESULT(yes);SIMD_OK=yes,
	AC_MSG_RESULT(no);)
  CFLAGS="$oldCFLAGS"
fi
  if test "x$SIMD_OK" != "xyes" ; then
    if test "x$USE_SIMD" = "xyes" ; then
      AC_MSG_ERROR([SIMD was requested, but the compiler doesn't seem to support it.])
    else
      USE_SIMD=no
      USE_SIMD_PP=no
      USE_SIMD_MOMR=no
    fi
  fi
  AC_SUBST(SIMD_CFLAGS)
fi

if test "x$USE_SIMD" != "xno" ; then
  AC_DEFINE([USE_SIMD], 1, [Define if SIMD optimizations should be used.])
fi
if test "x$USE_SIMD_PP" != "xno" ; then
  AC_DEFINE([USE_SIMD_PP],1, 
 	[Define if SIMD optimizations should be used for PP interactions.])
fi
if test "x$USE_SIMD_MOMR" = "xyes" ; then
  AC_DEFINE([USE_SIMD_MOMR],1, 
 	[Define if SIMD optimizations should be used for MOMRs.])
fi

AM_CONDITIONAL(USE_SIMD, test "x$USE_SIMD" != "xno")
AM_CONDITIONAL(USE_SIMD_PP,test "x$USE_SIMD_PP" != "xno")
AM_CONDITIONAL(USE_SIMD_MOMR,test "x$USE_SIMD_MOMR" = "xyes")

if test "x$HAVE_POSIX_MEMALIGN" == "xyes" ; then
  AC_DEFINE([HAVE_POSIX_MEMALIGN],1, 
 	[Define if posix_memalign should be used instead of malloc for SIMD.])
  AC_DEFINE([_GNU_SOURCE],1,
	[This is needed for posix_memalign to be defined in stdlib.h])
fi

AC_CHECK_LIB([bz2],[BZ2_bzWriteOpen])

# Checks for header files.
AC_FUNC_ALLOCA
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/param.h sys/time.h unistd.h ansidecl.h fenv.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SIGNAL
AC_HEADER_STDBOOL
AC_HEADER_TIME
AC_SYS_LARGEFILE

if test "x${INSTRUMENT}" == "xyes" ; then
AC_CHECK_HEADERS([c_asm.h intrinsics.h mach/mach_time.h])
AC_C_VOLATILE
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])

AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
dnl Cray UNICOS _rtc() (real-time clock) intrinsic
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
AC_MSG_RESULT($rtc_ok)
fi

# Checks for library functions.
AC_CHECK_FUNCS([wordexp wordfree glob globfree gethostname getpagesize mkdir])
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strchr strdup strrchr strstr],,[AC_MSG_ERROR([not found])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
