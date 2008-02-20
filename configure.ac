#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pkdgrav2m, 2.2, bugs@pkdgrav2.org)
AM_INIT_AUTOMAKE(pkdgrav2m,2.2)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET
AC_CHECK_PROG([MPICC],mpicc,mpicc,no)

# The floating point data type - This should be detected and --with-float added
AC_DEFINE([FLOAT],double,
	[This defines the base floating point type for particles.])
AC_DEFINE([FLOAT_MAXVAL],HUGE,
	[This defines the maximum value that a FLOAT can represent.])

# Local expansions
AC_ARG_ENABLE(local-expansion,
  [  --enable-local-expansion       enable local expansion code],
  [LOCAL_EXPANSION=yes])
# Local expansions
AC_ARG_ENABLE(fmm,
  [  --enable-fmm                   enable local expansion code],
  [LOCAL_EXPANSION=yes])
if test "x$LOCAL_EXPANSION" = "xyes" ; then
   AC_DEFINE([LOCAL_EXPANSION], 1,
 	[Define if local expansions should be used.])
fi
AM_CONDITIONAL(LOCAL_EXPANSION, test "x$LOCAL_EXPANSION" = "xyes")

# With particle-cell interactions
AC_ARG_ENABLE(particle-cell,
  [  --enable-particle-cell         enable particle-cell interactions],
  [WITH_PARTICLE_CELL=yes
   AC_DEFINE([WITH_PARTICLE_CELL], 1,
 	[Define if particle-cell interactions should be used.])])
AM_CONDITIONAL(WITH_PARTICLE_CELL, test "x$WITH_PARTICLE_CELL" = "xyes")

# HDF5 I/O
AC_ARG_ENABLE(hdf5,
  [  --enable-hdf5                  enable I/O using HDF5 libraries],
  [USE_HDF5=yes
   AC_DEFINE([USE_HDF5], 1,
 	[Define if HDF5 I/O should be available.])])
AC_ARG_WITH(hdf5,
  [  --with-hdf5=PATH               path to HDF5 installation],
  [CFLAGS="$CFLAGS -I${withval}/include"
   LDFLAGS="$LDFLAGS -L${withval}/lib"
   if test -n "${need_dash_r}"; then
	LDFLAGS="$LDFLAGS -R${withval}/lib"
   fi
   USE_HDF5=yes
   AC_DEFINE([USE_HDF5], 1,
 	[Define if HDF5 I/O should be available.])])


if test "${cross_compiling}" = "no";
then
AC_ARG_WITH(mdl,
  [  --with-mdl=PATH                base path of MDL source],
  [AC_CHECK_FILE([${withval}],[MDL_PATH="${withval}"],[AC_MSG_ERROR([MDL does not exist at ${withval}])])],
  [AC_CHECK_FILE([../mdl2],
	[AC_SUBST([MDL_PATH],../mdl2)],
	[AC_CHECK_FILE([../mdl],[AC_SUBST([MDL_PATH],../mdl)],[AC_MSG_ERROR([Cannot locate MDL, try --with-mdl=PATH])])]
   )])
else
AC_ARG_WITH(mdl,
  [  --with-mdl=PATH                base path of MDL source],
  [MDL_PATH="${withval}"]
  [MDL_PATH="../mdl2"])
fi

AC_ARG_WITH(grafic,
  [  --with-grafic=PATH             base path of GRAFIC IC library],
  [AC_SUBST([GRAFIC_PATH],[${withval}])
   USE_GRAFIC=yes
   AC_DEFINE([USE_GRAFIC], 1,
 	[Define GRAFICIC IC should be available.])],
  [USE_GRAFIC=no])
AM_CONDITIONAL(USE_GRAFIC,test "x$USE_GRAFIC" == "xyes")

# MDL dedicated I/O
AC_ARG_ENABLE(mdl-io,
  [  --enable-mdl-io                enable dedicated I/O processors (MPI/pthread)],
  [USE_MDL_IO=yes
   USE_HDF5=yes
   AC_DEFINE([USE_HDF5], 1,
 	[Define if HDF5 I/O should be available.])
   AC_DEFINE([USE_MDL_IO], 1,
 	[Define if MDL style I/O should be used.])
   AC_DEFINE([MDL_MAX_IO_PROCS],128,[The maximum number of IO processors])])
AM_CONDITIONAL(USE_MDL_IO, test "x$USE_MDL_IO" = "xyes")
AM_CONDITIONAL(USE_HDF5, test "x$USE_HDF5" = "xyes")


# SISCI MDL
AC_ARG_WITH(sisci,
  [  --with-sisci[=/opt/DIS]        path to SISCI installation],
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
  [  --with-gd=PATH                path to GD installation],
  [CFLAGS="$CFLAGS -I${withval}/include"
   LDFLAGS="$LDFLAGS -L${withval}/lib"
   if test -n "${need_dash_r}"; then
	LDFLAGS="$LDFLAGS -R${withval}/lib"
   fi
])

# Enable BSC Tracing
AC_ARG_ENABLE(bsc-trace,
  [  --enable-bsc-trace             enable tracing at BSC],
  [USE_BSC=yes
   AC_DEFINE([USE_BSC], 1,
 	[Define if BSC tracing should be used.])])
AM_CONDITIONAL(USE_BSC, test "x$USE_BSC" = "xyes")

# Enable SIMD version (SSE,Altivec)
AC_ARG_ENABLE(simd,
  [  --enable-simd                  enable SIMD optimizations],
  [USE_SIMD=yes
   AC_DEFINE([USE_SIMD], 1,
 	[Define if SIMD optimizations should be used.])])

USE_SIMD_LOCR="no"
AC_ARG_ENABLE(simd-locr,
  [  --enable-simd-locr            enable LOCR SIMD optimizations (experimental)],
  [if test "x$enableval" = "xyes" ; then
     USE_SIMD_LOCR="yes"
   elif test "x$enableval" = "xno" ; then
     USE_SIMD_LOCR="no"
   fi])
if test "x$USE_SIMD_LOCR" == "xyes" -a "x$LOCAL_EXPANSION" = "xyes" ; then
  USE_SIMD=yes
  AC_DEFINE([USE_SIMD_LOCR],1, 
 	[Define if SIMD optimizations should be used for LOCRs.])
fi
AM_CONDITIONAL(USE_SIMD_LOCR,test "x$USE_SIMD_LOCR" != "xno")

USE_SIMD_MOMR="no"
AC_ARG_ENABLE(simd-momr,
  [  --enable-simd-momr             enable MOMR SIMD optimizations],
  [if test "x$enableval" = "xyes" ; then
     USE_SIMD_MOMR="yes"
   elif test "x$enableval" = "xno" ; then
     USE_SIMD_MOMR="no"
   fi])
if test "x$USE_SIMD_MOMR" = "xyes" ; then
  USE_SIMD=yes
  AC_DEFINE([USE_SIMD_MOMR],1, 
 	[Define if SIMD optimizations should be used for MOMRs.])
fi
AM_CONDITIONAL(USE_SIMD_MOMR,test "x$USE_SIMD_MOMR" != "xno")

USE_SIMD_PP="$USE_SIMD"
AC_ARG_ENABLE(simd-pp,
  [  --disable-simd-pp              disable PP SIMD optimizations],
  [if test "x$enableval" = "xyes" ; then
     USE_SIMD_PP="yes"
   elif test "x$enableval" = "xno" ; then
     USE_SIMD_PP="no"
   fi])
if test "x$USE_SIMD_PP" == "xyes" ; then
  USE_SIMD=yes
  AC_DEFINE([USE_SIMD_PP],1, 
 	[Define if SIMD optimizations should be used for PP interactions.])
fi
AM_CONDITIONAL(USE_SIMD_PP,test "x$USE_SIMD_PP" == "xyes")

AM_CONDITIONAL(USE_SIMD, test "x$USE_SIMD" = "xyes")

AC_ARG_ENABLE(simd-locr,
  [  --disable-simd-locr            disable LOCR SIMD optimizations],
  [if test "x$enableval" = "xyes" ; then
     USE_SIMD_LOCR="yes"
   elif test "x$enableval" = "xno" ; then
     USE_SIMD_LOCR="no"
   else
     USE_SIMD_LOCR="$USE_SIMD"
   fi])


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

# Checks for libraries.
if test "x$USE_HDF5" = "xyes" ; then
	save_CC="$CC"
	 CC=mpicc
	 AC_CHECK_LIB([z], [inflate])
	 AC_CHECK_LIB([gpfs], [gpfs_fcntl])
	 AC_CHECK_LIB([hdf5], [H5open])
	CC="$save_CC"
fi
AC_CHECK_LIB([m], [fabs])
# FIXME: Replace `main' with a function in `-lrt':
#AC_CHECK_LIB([rt], [main])
# FIXME: Replace `main' with a function in `-lz':
#AC_CHECK_LIB([z], [main])

if test "x$USE_SIMD" == "xyes" ; then
  AC_CHECK_FUNCS([posix_memalign],HAVE_POSIX_MEMALIGN=yes)

  AC_MSG_CHECKING([if the compiler accepts -msse3])
  oldCFLAGS="$CFLAGS"; CFLAGS="-msse3"
  AC_TRY_COMPILE([#include <pmmintrin.h>],[__m128 r;],
	AC_MSG_RESULT(yes); SIMD_CFLAGS="-msse3",
	AC_MSG_RESULT(no);)
  CFLAGS="$oldCFLAGS"
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
  if test "x$SIMD_CFLAGS" = "x" ; then
    AC_MSG_ERROR([SIMD was requested, but the compiler doesn't seem to support it.])
  fi
  AC_SUBST(SIMD_CFLAGS)
fi
if test "x$HAVE_POSIX_MEMALIGN" == "xyes" ; then
  AC_DEFINE([HAVE_POSIX_MEMALIGN],1, 
 	[Define if posix_memalign should be used instead of malloc for SIMD.])
  AC_DEFINE([_GNU_SOURCE],1,
	[This is needed for posix_memalign to be defined in stdlib.h])
fi

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_SIGNAL
AC_HEADER_STDBOOL
AC_HEADER_TIME
AC_SYS_LARGEFILE

# Checks for library functions.
AC_CHECK_FUNCS([floor gethostname gettimeofday memset pow sqrt strchr strdup strrchr strstr ])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
