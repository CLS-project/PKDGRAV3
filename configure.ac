#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pkdgrav2m, 2.1.1, bugs@pkdgrav2.org)
AM_INIT_AUTOMAKE(pkdgrav2m,2.1)
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
  [  --enable-local-expansion       enable local expansion code (experimental)],
  [LOCAL_EXPANSION=yes
   AC_DEFINE([LOCAL_EXPANSION], 1,
 	[Define if local expansions should be used.])])
AM_CONDITIONAL(LOCAL_EXPANSION, test "x$LOCAL_EXPANSION" = "xyes")

# HDF5 I/O
AC_ARG_ENABLE(hdf5,
  [  --enable-hdf5                  enable I/O using HDF5 libraries],
  [USE_HDF5=yes
   AC_DEFINE([USE_HDF5], 1,
 	[Define if HDF5 I/O should be available.])])

AC_CHECK_FILE([../mdl2],
	AC_SUBST([MDL_PATH],../mdl2),
	AC_SUBST([MDL_PATH],../mdl))

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

# Enable SIMD version (SSE,Altivec)
AC_ARG_ENABLE(simd,
  [  --enable-simd                  enable SIMD optimizations],
  [USE_SIMD=yes
   AC_DEFINE([USE_SIMD], 1,
 	[Define if SIMD optimizations should be used.])])
AM_CONDITIONAL(USE_SIMD, test "x$USE_SIMD" = "xyes")

AC_ARG_ENABLE(simd-momr,
  [  --disable-simd-momr            disable MOMR SIMD optimizations],
  [if test "x$enableval" = "xyes" ; then
     USE_SIMD_MOMR="yes"
   elif test "x$enableval" = "xno" ; then
     USE_SIMD_MOMR="no"
   else
     USE_SIMD_MOMR="$USE_SIMD"
   fi])
if test "x$USE_SIMD_MOMR" != "xno" -a "x$USE_SIMD" = "xyes" ; then
  AC_DEFINE([USE_SIMD_MOMR],1, 
 	[Define if SIMD optimizations should be used for MOMRs.])
fi
AM_CONDITIONAL(USE_SIMD_MOMR,test "x$USE_SIMD_MOMR" != "xno")

# The floating point data type for moments
if test "x$USE_SIMD" = "xyes" ; then
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
  AC_CHECK_LIB([hdf5], [H5open])
fi
AC_CHECK_LIB([m], [fabs])
# FIXME: Replace `main' with a function in `-lrt':
#AC_CHECK_LIB([rt], [main])
# FIXME: Replace `main' with a function in `-lz':
#AC_CHECK_LIB([z], [main])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_STDC
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h stddef.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_SYS_LARGEFILE

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gethostname gettimeofday pow sqrt strchr strdup])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
