EXTRA_PROGRAMS = pkdgrav2m_null
TARGETS = null

if HAVE_PTHREAD
EXTRA_PROGRAMS	+= pkdgrav2m_pthread
TARGETS	+= pthread
pthread		: pkdgrav2m_pthread
endif

if HAVE_MPI
EXTRA_PROGRAMS	+= pkdgrav2m_mpi
TARGETS	+= mpi
mpi:
	$(MAKE) CC=mpicc $(AM_MAKEFLAGS) pkdgrav2m_mpi
endif

null	: pkdgrav2m_null

all-local: $(TARGETS)

clean-local:
	rm -f $(EXTRA_PROGRAMS)

SOURCES  = cosmo.c master.c param.c pst.c smooth.c tree.c ewald.c
SOURCES += main.c moments.c outtype.c pkd.c romberg.c smoothfcn.c
if LOCAL_EXPANSION
SOURCES += ../bsc/grav2.c ../bsc/walk2.c
if USE_SIMD
AM_CFLAGS = -msse3 -I../bsc/
SOURCES += ../bsc/momGenLocrAddSIMDMomr.c
SOURCES += ../bsc/momGenLocrAddVMomr.c
if USE_SIMD_MOMR
SOURCES += ../bsc/momEvalSIMDMomr.c
endif
endif
else
SOURCES += grav.c walk.c
endif
if USE_HDF5
SOURCES += iohdf5.c
endif
if USE_MDL_IO
SOURCES += io.c
MDL=../bsc
else
MDL=@MDL_PATH@
endif

pkdgrav2m_null_SOURCES = $(SOURCES) $(MDL)/null/mdl.c
pkdgrav2m_null_CFLAGS = -I$(MDL)/null/ $(AM_CFLAGS)

if HAVE_MPI
pkdgrav2m_mpi_SOURCES = $(SOURCES) $(MDL)/mpi/mdl.c
pkdgrav2m_mpi_CFLAGS = -I$(MDL)/mpi/ $(AM_CFLAGS)
endif

if HAVE_PTHREAD
pkdgrav2m_pthread_SOURCES = $(SOURCES) $(MDL)/pthread/mdl.c
pkdgrav2m_pthread_CFLAGS = -I$(MDL)/pthread/ $(AM_CFLAGS)
pkdgrav2m_pthread_LDADD = -lpthread
endif

