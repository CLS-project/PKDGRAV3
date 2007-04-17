bin_PROGRAMS = pkdgrav2m_mpi pkdgrav2m_pthread

# Pseudo targets for 'make pthread' and 'make mpi'
pthread		: pkdgrav2m_pthread
mpi		: pkdgrav2m_mpi

SOURCES = cosmo.c master.c param.c pst.c smooth.c tree.c ewald.c\
 main.c moments.c outtype.c pkd.c romberg.c smoothfcn.c\
 ../bsc/grav2.c ../bsc/walk2.c
if USE_MDL_IO
SOURCES += io.c
endif
if USE_SIMD
SOURCES += ../bsc/momGenLocrAddSIMDMomr.c
endif

AM_CFLAGS = -msse3 -DSINGLE

pkdgrav2m_mpi_SOURCES = $(SOURCES) ../bsc/mpi/mdl.c
pkdgrav2m_mpi_CC = mpicc
pkdgrav2m_mpi_CFLAGS = -DMDL_MPI -I../bsc/mpi/ $(AM_CFLAGS)

pkdgrav2m_pthread_SOURCES = $(SOURCES) ../mdl/pthread/mdl.c
pkdgrav2m_pthread_CFLAGS = -DMDL_PTHREAD -I../mdl/pthread/ $(AM_CFLAGS)

