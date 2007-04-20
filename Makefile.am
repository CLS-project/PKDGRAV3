bin_PROGRAMS = pkdgrav2m_mpi pkdgrav2m_pthread

# Pseudo targets for 'make pthread' and 'make mpi'
pthread		: pkdgrav2m_pthread
mpi		: pkdgrav2m_mpi

SOURCES  = cosmo.c master.c param.c pst.c smooth.c tree.c ewald.c
SOURCES += main.c moments.c outtype.c pkd.c romberg.c smoothfcn.c
if LOCAL_EXPANSION
SOURCES += ../bsc/grav2.c ../bsc/walk2.c
if USE_SIMD
SOURCES += ../bsc/momGenLocrAddSIMDMomr.c
endif
else
SOURCES += grav.c walk.c
endif
if USE_MDL_IO
SOURCES += io.c
MDL=../bsc
else
MDL=../mdl
endif
AM_CFLAGS = -msse3 -DSINGLE

pkdgrav2m_mpi_SOURCES = $(SOURCES) $(MDL)/mpi/mdl.c
pkdgrav2m_mpi_CFLAGS = -DMDL_MPI -I$(MDL)/mpi/ $(AM_CFLAGS)

pkdgrav2m_pthread_SOURCES = $(SOURCES) $(MDL)/pthread/mdl.c
pkdgrav2m_pthread_CFLAGS = -DMDL_PTHREAD -I../mdl/pthread/ $(AM_CFLAGS)

