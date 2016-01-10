#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pkdgrav3, 0.9.18, bugs@pkdgrav.org)
AM_INIT_AUTOMAKE
AM_MAINTAINER_MODE([disable])
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET
AM_PROG_CC_C_O

AC_ARG_WITH(mpicc,
  [  --with-mpicc=mpicc      specify the name of the MPI C compiler],
  [MPICC="${withval}"],[MPICC=""])
if test "$MPICC" == "no" ; then
  MPICC=""
else
  AC_LANG_PUSH(C)
  if test "$MPICC" = "yes" ; then
    MPICC=""
  fi
  oldCC="$CC"
  for CC in cc mpicc $oldCC ; do
    test -n "$MPICC" && break
    AC_TRY_COMPILE([#include "mpi.h"],[MPI_Init(0,0)],[MPICC="$CC"],[])
  done
  CC="$oldCC"
  AC_LANG_POP(C)
fi
AC_SUBST(MPICC)

AC_DEFINE([HAVE_ROMBERG], 1, [Use the pkdgrav3 Romberg integrator.])
AC_DEFINE([NAIVE_DOMAIN_DECOMP], 1, [Use a simplified domain decomposition.])

# Cooling
AC_ARG_ENABLE(cooling,
  [  --enable-cooling             enable cooling],
  [COOLING=yes
   AC_DEFINE([COOLING], 1,
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
AC_ARG_WITH(fftw,
  [  --with-fftw                    Enable FFT support modules],
  [FFTW="$enableval"],[FFTW="no"])
if test "x$FFTW" != "xno" ; then
  AC_DEFINE([USE_SINGLE], 1,
        [Use single precision FFT.])
  AC_CHECK_LIB([fftw3f], [fftwf_execute], [],[AC_MSG_ERROR([Cannot find FFTW3])])
  AC_CHECK_LIB([fftw3f_mpi], [fftwf_mpi_init],[],[AC_MSG_ERROR([Cannot find FFTW3])])
  AC_CHECK_LIB([fftw3f_threads], [fftwf_init_threads],[],[AC_MSG_ERROR([Cannot find FFTW3])])
  AC_DEFINE([MDL_FFTW], 1, [Include support for FFTW.])
  AM_PATH_GSL(1.8,[],[AC_MSG_ERROR([Unable to locate GNU Scientific Library (gsl)])])
fi
AM_CONDITIONAL(FFTW, test "x$FFTW" != "xno")

# Electric fence
AC_ARG_WITH(efence,
  [  --with-efence=PATH      enable electric fence debugging],
  [EFENCE_DIR="${withval}"],[EFENCE_DIR=no])
AM_CONDITIONAL(USE_EFENCE, test "x$EFENCE_DIR" != "xno")
AM_CONDITIONAL(HAVE_EFENCE_PATH, test "x$EFENCE_DIR" != "xyes")
AC_SUBST([EFENCE_PATH],[$EFENCE_DIR])

# Intel VTUNE
AC_ARG_WITH(itt,
  [  --with-itt=PATH      enable Vtune debugging],
  [ITT_DIR="${withval}"],[ITT_DIR=no])
AM_CONDITIONAL(USE_ITT, test "x$ITT_DIR" != "xno")
AM_CONDITIONAL(HAVE_ITT_PATH, test "x$ITT_DIR" != "xyes")
AC_SUBST([ITT_PATH],[$ITT_DIR])
if test "x$ITT_DIR" != "xno" ; then
AC_DEFINE([USE_ITT],1,
  [Defined if Intel vtune support should be compiled into pkdgrav3])
fi

AC_CHECK_HEADERS([libaio.h],[USE_LIBAIO=yes],[USE_LIBAIO=no])
AC_CHECK_HEADERS([aio.h],[USE_AIO=yes],[USE_AIO=no])
AM_CONDITIONAL(USE_LIBAIO, test "$USE_LIBAIO" != "no")
AM_CONDITIONAL(USE_AIO, test "$USE_AIO" != "no")

# Python
AC_ARG_WITH(python,
  [  --with-python=VERSION   enable python scripting],
  [PYTHON_EXE="${withval}"],[PYTHON_EXE=no])
AM_CONDITIONAL(USE_PYTHON, test "x$PYTHON_EXE" != "xno")
if test "x$PYTHON_EXE" != "xno" ; then
AC_DEFINE([USE_PYTHON],1,
  [Defined if python support should be compiled into pkdgrav3])
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
#elif test "x$USE_HDF5" != "xno" ; then
#  USE_HDF5="yes"
#  AC_CHECK_LIB([hdf5],[H5Fopen],,[USE_HDF5="no"])
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

AM_CONDITIONAL(USE_HDF5, test "x$USE_HDF5" = "xyes")

AC_ARG_ENABLE(integer-positions,
  [  --enable-integer-positions     enable Integerized positions to save memory)],
  [AC_DEFINE([INTEGER_POSITION], 1,
        [Define if Integerized positions are to be used.])])

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

# Python
AC_ARG_WITH(cuda,
  [  --with-cuda=PATH               enable cude acceleration],
  [CUDA_DIR="${withval}"],[CUDA_DIR=no])
AM_CONDITIONAL(USE_CUDA, test "x$CUDA_DIR" != "xno")
if test "x$CUDA_DIR" != "xno" ; then
AC_DEFINE([USE_CUDA],1,
  [Use CUDA to accelerate calculations])
if test "x$CUDA_DIR" = "xyes" ; then
  NVCC=nvcc
  CUDA_LIBS=""
else
  NVCC=$CUDA_DIR/bin/nvcc
  CUDA_LIBS="-L$CUDA_DIR/lib64 -L$CUDA_DIR/lib"
fi
NVCCFLAGS="-arch compute_35"
AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)
AC_SUBST(CUDA_LIBS)
fi

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

USE_SIMD_EWALD="$USE_SIMD"
AC_ARG_ENABLE(simd-ewald,
  [  --disable-simd-ewald           disable Ewald SIMD optimizations],
  [USE_SIMD_EWALD="$enableval"])
if test "x$USE_SIMD_EWALD" != "xno" ; then
  USE_SIMD="$USE_SIMD_EWALD"
fi

USE_SIMD_OPEN="$USE_SIMD"
AC_ARG_ENABLE(simd-open,
  [  --disable-simd-open            disable Opening Criteria SIMD optimizations],
  [USE_SIMD_OPEN="$enableval"])
if test "x$USE_SIMD_OPEN" != "xno" ; then
  USE_SIMD="$USE_SIMD_OPEN"
fi

USE_SIMD_PP="$USE_SIMD"
AC_ARG_ENABLE(simd-pp,
  [  --disable-simd-pp              disable PP SIMD optimizations],
  [USE_SIMD_PP="$enableval"])
if test "x$USE_SIMD_PP" != "xno" ; then
  USE_SIMD="$USE_SIMD_PP"
fi

USE_SIMD_PC="$USE_SIMD"
AC_ARG_ENABLE(simd-pc,
  [  --disable-simd-pc              disable PC SIMD optimizations],
  [USE_SIMD_PC="$enableval"])
if test "x$USE_SIMD_PC" != "xno" ; then
  USE_SIMD="$USE_SIMD_PC"
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

# Check if mpi is available
if test -n "$MPICC" ; then
  HAVE_MPI=yes
fi
AM_CONDITIONAL(HAVE_MPI,test "x$HAVE_MPI" = "xyes")
if test "x$HAVE_MPI" = "xyes" ; then
  AC_DEFINE([HAVE_MPI], 1, [Set to 1 if MPI is available.])
fi
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
      USE_SIMD_OPEN=no
      USE_SIMD_EWALD=no
      USE_SIMD_PP=no
      USE_SIMD_PC=no
      USE_SIMD_MOMR=no
    fi
  fi
  AC_SUBST(SIMD_CFLAGS)
fi

if test "x$USE_SIMD" != "xno" ; then
  AC_DEFINE([USE_SIMD], 1, [Define if SIMD optimizations should be used.])
fi
if test "x$USE_SIMD_OPEN" != "xno" ; then
  AC_DEFINE([USE_SIMD_OPEN],1, 
 	[Define if SIMD optimizations should be used for the opening criteria.])
fi
if test "x$USE_SIMD_PP" != "xno" ; then
  AC_DEFINE([USE_SIMD_PP],1, 
 	[Define if SIMD optimizations should be used for PP interactions.])
fi
if test "x$USE_SIMD_PC" != "xno" ; then
  AC_DEFINE([USE_SIMD_PC],1, 
 	[Define if SIMD optimizations should be used for PC interactions.])
fi
if test "x$USE_SIMD_EWALD" != "xno" ; then
  AC_DEFINE([USE_SIMD_EWALD],1,
        [Define if SIMD optimizations should be used for Ewald.])
fi
if test "x$USE_SIMD_MOMR" = "xyes" ; then
  AC_DEFINE([USE_SIMD_MOMR],1, 
 	[Define if SIMD optimizations should be used for MOMRs.])
fi

AM_CONDITIONAL(USE_SIMD, test "x$USE_SIMD" != "xno")
AM_CONDITIONAL(USE_SIMD_OPEN,test "x$USE_SIMD_OPEN" != "xno")
AM_CONDITIONAL(USE_SIMD_PP,test "x$USE_SIMD_PP" != "xno")
AM_CONDITIONAL(USE_SIMD_PC,test "x$USE_SIMD_PC" != "xno")
AM_CONDITIONAL(USE_SIMD_MOMR,test "x$USE_SIMD_MOMR" = "xyes")

if test "x$HAVE_POSIX_MEMALIGN" == "xyes" ; then
  AC_DEFINE([HAVE_POSIX_MEMALIGN],1, 
 	[Define if posix_memalign should be used instead of malloc for SIMD.])
  AC_DEFINE([_GNU_SOURCE],1,
	[This is needed for posix_memalign to be defined in stdlib.h])
fi

AC_MSG_CHECKING([if the linker accepts -rdynamic])
oldLDFLAGS="$LDFLAGS"; LDFLAGS="-Xlinker -rdynamic"
AC_TRY_LINK([],[],
AC_MSG_RESULT(yes); HAVE_RDYNAMIC="yes",
AC_MSG_RESULT(no);)
LDFLAGS="$oldLDFLAGS"
AM_CONDITIONAL(HAVE_RDYNAMIC,test "x$HAVE_RDYNAMIC" = "xyes")

AC_CHECK_LIB([bz2],[BZ2_bzWriteOpen])

# Checks for header files.
AC_FUNC_ALLOCA
AC_FUNC_FSEEKO
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/param.h sys/time.h unistd.h ansidecl.h fenv.h signal.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_HEADER_STDBOOL
AC_HEADER_TIME
AC_SYS_LARGEFILE

if test "x${INSTRUMENT}" == "xyes" ; then
AC_CHECK_HEADERS([c_asm.h intrinsics.h mach/mach_time.h])
AC_C_VOLATILE
AC_CHECK_TYPE([hrtime_t],[AC_DEFINE(HAVE_HRTIME_T, 1, [Define to 1 if hrtime_t is defined in <sys/time.h>])],,[#if HAVE_SYS_TIME_H
#include <sys/time.h>
#endif])

AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time atexit memmove memset malloc_stats])
dnl Cray UNICOS _rtc() (real-time clock) intrinsic
AC_MSG_CHECKING([for _rtc intrinsic])
rtc_ok=yes
AC_TRY_LINK([#ifdef HAVE_INTRINSICS_H
#include <intrinsics.h>
#endif], [_rtc()], [AC_DEFINE(HAVE__RTC,1,[Define if you have the UNICOS _rtc() intrinsic.])], [rtc_ok=no])
AC_MSG_RESULT($rtc_ok)
fi

# Checks for library functions.
AC_CHECK_FUNCS([wordexp wordfree glob globfree gethostname getpagesize mkdir strverscmp])
AC_CHECK_FUNCS([floor gettimeofday pow sqrt strchr strdup strrchr strstr],,[AC_MSG_ERROR([not found])])
AX_PTHREAD
AM_CONDITIONAL(HAVE_PTHREAD,test "x$ax_pthread_ok" = "xyes")
AM_CONDITIONAL(HAVE_MPIPT,false)

AC_CONFIG_FILES([Makefile])
AC_CONFIG_SUBDIRS([openpa])
#if test "x$USE_GRAFIC" != "xno" ; then
#  AC_CONFIG_SUBDIRS([grafic])
#fi

AC_OUTPUT
