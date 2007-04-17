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

# Optional features
AC_ARG_ENABLE(mdl-io, 
  [  --enable-mdl-io                enable dedicated I/O processors (MPI/pthread)],
  [AC_DEFINE([USE_MDL_IO], 1,
 	[Define if MDL style I/O should be used.])])

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
