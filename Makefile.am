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

if USE_HDF5
EXTRA_PROGRAMS += tip2hdf5 hdf52tip hdf5cat
endif

null	: pkdgrav2m_null

all-local: $(TARGETS)

clean-local:
	rm -f $(EXTRA_PROGRAMS)

SOURCES  = cosmo.c master.c param.c pst.c smooth.c tree.c ewald.c
SOURCES += main.c moments.c outtype.c pkd.c romberg.c smoothfcn.c
if PLANETS
SOURCES += ssio.c collision.c kepler.c
endif
if LOCAL_EXPANSION
SOURCES += grav2.c walk2.c
if USE_SIMD
AM_CFLAGS = -I../bsc/ @SIMD_CFLAGS@
if USE_SIMD_LOCR
SOURCES += ../bsc/momGenLocrAddSIMDMomr.c
endif
if USE_SIMD_MOMR
SOURCES += ../bsc/momEvalSIMDMomr.c
endif
if USE_SIMD_PP
SOURCES += ../bsc/PPInteractSIMD.c
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
if USE_PNG
SOURCES += png.c
endif

pkdgrav2m_null_SOURCES = $(SOURCES) $(MDL)/null/mdl.c
pkdgrav2m_null_CFLAGS = -I$(MDL)/null/ $(AM_CFLAGS)
pkdgrav2m_null_LDADD =
if USE_PNG
pkdgrav2m_null_LDADD += -lgd -lpng
endif
if USE_HDF5
pkdgrav2m_null_LDADD += -lhdf5
endif


if HAVE_MPI
pkdgrav2m_mpi_LDADD =
pkdgrav2m_mpi_SOURCES = $(SOURCES) $(MDL)/mpi/mdl.c
pkdgrav2m_mpi_CFLAGS = -I$(MDL)/mpi/ $(AM_CFLAGS)
if USE_BSC
pkdgrav2m_mpi_LDADD += -L/gpfs/apps/CEPBATOOLS/mpitrace-devel/64/lib -ldl -lmpitrace_burst
pkdgrav2m_mpi_CFLAGS += -I/gpfs/apps/CEPBATOOLS/mpitrace-devel/64/include
endif
if USE_PNG
pkdgrav2m_mpi_LDADD += -lgd -lpng
endif
if USE_HDF5
pkdgrav2m_mpi_LDADD += -lhdf5
endif
endif

if HAVE_PTHREAD
pkdgrav2m_pthread_SOURCES = $(SOURCES) $(MDL)/pthread/mdl.c
pkdgrav2m_pthread_CFLAGS = -I$(MDL)/pthread/ $(AM_CFLAGS)
pkdgrav2m_pthread_LDADD = -lpthread
if USE_PNG
pkdgrav2m_pthread_LDADD += -lgd -lpng
endif
if USE_HDF5
pkdgrav2m_pthread_LDADD += -lhdf5
endif
endif

if USE_HDF5
tip2hdf5_SOURCES = tip2hdf5.c iohdf5.c
hdf52tip_SOURCES = hdf52tip.c iohdf5.c
hdf5cat_SOURCES = hdf5cat.c iohdf5.c
endif

