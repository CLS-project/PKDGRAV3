#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pkdgrav2m, 2.1, bugs@pkdgrav.com)
AM_INIT_AUTOMAKE(pkdgrav2m,2.1)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC([pgcc icc gcc cc])
AC_PROG_MAKE_SET

# Local expansions
AC_ARG_ENABLE(local-expansion,
  [  --enable-local-expansion       enable local expansion code (experimental)],
  [LOCAL_EXPANSION=true
   AC_DEFINE([LOCAL_EXPANSION], 1,
 	[Define if local expansions should be used.])])
AM_CONDITIONAL(LOCAL_EXPANSION, test "x$LOCAL_EXPANSION" = "xtrue")

# MDL dedicated I/O
AC_ARG_ENABLE(mdl-io,
  [  --enable-mdl-io                enable dedicated I/O processors (MPI/pthread)],
  [USE_MDL_IO=true
   AC_DEFINE([USE_MDL_IO], 1,
 	[Define if MDL style I/O should be used.])])
AM_CONDITIONAL(USE_MDL_IO, test "x$USE_MDL_IO" = "xtrue")

# Enable SIMD version (SSE,Altivec)
AC_ARG_ENABLE(simd,
  [  --enable-simd                  enable SIMD optimizations],
  [USE_SIMD=true
   AC_DEFINE([USE_SIMD], 1,
 	[Define if SIMD optimizations should be used.])])
AM_CONDITIONAL(USE_SIMD, test "x$USE_SIMD" = "xtrue")

# Relaxation
AC_ARG_ENABLE(relaxation,
  [  --enable-relaxation            enable relaxation],
  [RELAXATION=true
   AC_DEFINE([RELAXATION], 1,
 	[Define if RELAXATION should be used.])])
AM_CONDITIONAL(RELAXATION, test "x$RELAXATION" = "xtrue")

# Soft Linear
AC_ARG_ENABLE(soft-linear,
  [  --enable-soft-linear           enable soft linear],
  [SOFTLINEAR=true
   AC_DEFINE([SOFTLINEAR], 1,
 	[Define if SOFTLINEAR should be used.])])
AM_CONDITIONAL(SOFTLINEAR, test "x$SOFTLINEAR" = "xtrue")

# Relaxation
AC_ARG_ENABLE(soft-square,
  [  --enable-soft-square           enable soft square],
  [SOFTSQUARE=true
   AC_DEFINE([SOFTSQUARE], 1,
 	[Define if SOFTSQUARE should be used.])])
AM_CONDITIONAL(SOFTSQUARE, test "x$SOFTSQUARE" = "xtrue")

# Checks for libraries.
AC_CHECK_LIB([hdf5], [H5open])
AC_CHECK_LIB([m], [fabs])
AC_CHECK_LIB([mpi], [MPI_Init])
AC_CHECK_LIB([pthread], [pthread_create])
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

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([floor gethostname gettimeofday pow sqrt strchr strdup])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
