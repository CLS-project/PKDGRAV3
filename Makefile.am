bin_PROGRAMS = pkdgrav2m_mpi pkdgrav2m_pthread

# Pseudo targets for 'make pthread' and 'make mpi'
pthread		: pkdgrav2m_pthread
mpi		: pkdgrav2m_mpi

SOURCES = cosmo.c grav2.c io.c master.c momGenLocrAddSIMDMomr.c param.c pst.c\
 smooth.c tree.c walk.c ewald.c grav.c main.c moments.c outtype.c pkd.c\
 romberg.c smoothfcn.c walk2.c

AM_CFLAGS = -msse3

pkdgrav2m_mpi_SOURCES = $(SOURCES) ../bsc/mpi/mdl.c
pkdgrav2m_mpi_CFLAGS = -DMDL_MPI -I../bsc/mpi/

pkdgrav2m_pthread_SOURCES = $(SOURCES) ../mdl/pthread/mdl.c
pkdgrav2m_pthread_CFLAGS = -DMDL_PTHREAD -I../mdl/pthread/

