SUBDIRS=openpa
ACLOCAL_AMFLAGS  = -I m4

SUFFIXES = .cu

.cu.o:
	$(NVCC) $(CFLAGS) -I. -DHAVE_CONFIG_H -o $@ -c $< $(NVCCFLAGS) -Imdl2/mpi -Imdl2 -Iopenpa/src


#if USE_GRAFIC
#SUBDIRS=grafic
#endif

EXTRA_DIST = server.py
EXTRA_PROGRAMS = pkdgrav2_null
TARGETS = null

if HAVE_PTHREAD
EXTRA_PROGRAMS	+= pkdgrav2_pthread
TARGETS	+= pthread
pthread		: pkdgrav2_pthread
endif

if HAVE_MPI
EXTRA_PROGRAMS	+= pkdgrav2_mpi
TARGETS	+= mpi

mpi:
	$(MAKE) CC=@MPICC@ FC=@MPIFC@ $(AM_MAKEFLAGS) pkdgrav2_mpi
endif

EXTRA_PROGRAMS += tostd tdiff

null	: pkdgrav2_null

all-local: $(TARGETS)

clean-local:
	rm -f $(EXTRA_PROGRAMS)

SOURCES = main.c tree.c analysis.c romberg.c bt.c
noinst_HEADERS = bt.h parameters.h grav.h walk.h cycle.h tipsydefs.h
SOURCES        += cosmo.c master.c param.c pst.c smooth.c ewald.c
noinst_HEADERS += cosmo.h master.h param.h pst.h smooth.h ewald.h meval.h qeval.h
SOURCES        += moments.c intype.c outtype.c pkd.c smoothfcn.c rbtree.c fio.c listcomp.c illinois.c
noinst_HEADERS += moments.h intype.h outtype.h pkd.h smoothfcn.h rbtree.h fio.h listcomp.h illinois.h
if COOLING
#SOURCES        += cooling.c cooling_module.F90
#noinst_HEADERS += cooling.h
endif

SOURCES += psdtree.c psd.c unbind.c knn6d.c
noinst_HEADERS += psdtree.h psd.h unbind.h knn6d.h qsort.h

SOURCES += hop.c
noinst_HEADERS += hop.h

SOURCES += $(MDL)/mdlbase.c
noinst_HEADERS += $(MDL)/mdlbase.h

if PLANETS
SOURCES += ssio.c collision.c kepler.c
noinst_HEADERS += ssio.h collision.h kepler.h
endif
if LOCAL_EXPANSION
SOURCES += grav2.c walk2.c
SOURCES += lst.c ilp.c ilc.c cl.c
noinst_HEADERS += lst.h ilp.h ilc.h cl.h
if USE_SIMD
AM_CFLAGS = @SIMD_CFLAGS@
noinst_HEADERS += simd.h sse_mathfun.h avx_mathfun.h
endif
else
SOURCES += grav.c walk.c
endif
MDL=@MDL_PATH@
if USE_PNG
SOURCES += png.c
noinst_HEADERS += png.h
endif
if USE_PYTHON
SOURCES += pkdpython.c
noinst_HEADERS += pkdpython.h
endif
if USE_CUDA
SOURCES += cudautil.cu cudapp.cu cudapc.cu
noinst_HEADERS += cudautil.h
endif

pkdgrav2_null_SOURCES = $(SOURCES) $(MDL)/null/mdl.c
noinst_HEADERS += $(MDL)/null/mdl.h
pkdgrav2_null_CFLAGS =  -I$(MDL)/null/ -I$(MDL) $(AM_CFLAGS)
pkdgrav2_null_LDADD = $(FLIBS)
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav2_null_LDADD += -L @EFENCE_PATH@
endif
pkdgrav2_null_LDADD += -lefence
endif
if USE_CUDA
pkdgrav2_null_LDADD += @CUDA_LIBS@
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav2_null_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav2_null_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lrt -lutil
pkdgrav2_null_CFLAGS += -I$(PYTHON_INC)
endif
if LUSTRE
pkdgrav2_null_LDADD += -llustreapi
endif

if HAVE_MPI
pkdgrav2_mpi_LDADD = -lpthread
pkdgrav2_mpi_SOURCES = $(SOURCES) $(MDL)/mpi/mdl.c openpa/src/opa_primitives.c openpa/src/opa_queue.c
noinst_HEADERS += $(MDL)/mpi/mdl.h
pkdgrav2_mpi_CFLAGS = -I$(MDL)/mpi/ -I$(MDL) -Iopenpa/src $(AM_CFLAGS)
if USE_BSC
pkdgrav2_mpi_LDADD += -L /gpfs/apps/CEPBATOOLS/mpitrace-devel/64/lib -ldl -lmpitrace_burst
pkdgrav2_mpi_CFLAGS += -I/gpfs/apps/CEPBATOOLS/mpitrace-devel/64/include
endif
if FFTW
pkdgrav2_mpi_LDADD += -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
#if USE_GRAFIC
#pkdgrav2_mpi_CFLAGS += -Igrafic/
#pkdgrav2_mpi_LDADD += grafic/libgrafic.a
##pkdgrav2_mpi_LDADD += -L @GRAFIC_PATH@ -lgrafic -lsrfftw_mpi -lsfftw_mpi -lsrf#ftw -lsfftw -lgfortran
##CCLD=mpicxx
#endif
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav2_mpi_LDADD += -L @EFENCE_PATH@
endif
pkdgrav2_mpi_LDADD += -lefence
endif
if USE_CUDA
pkdgrav2_mpi_LDADD += @CUDA_LIBS@
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav2_mpi_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav2_mpi_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lutil
pkdgrav2_mpi_CFLAGS += -I$(PYTHON_INC)
endif
if LUSTRE
pkdgrav2_mpi_LDADD += -llustreapi
endif
endif

if HAVE_PTHREAD
pkdgrav2_pthread_SOURCES = $(SOURCES) $(MDL)/pthread/mdl.c
noinst_HEADERS += $(MDL)/pthread/mdl.h
pkdgrav2_pthread_CFLAGS = -DUSE_PTHREAD -I$(MDL)/pthread/ -I$(MDL) $(AM_CFLAGS)
pkdgrav2_pthread_LDADD = -lpthread
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav2_pthread_LDADD += -L @EFENCE_PATH@
endif
pkdgrav2_pthread_LDADD += -lefence
endif
if USE_CUDA
pkdgrav2_pthread_LDADD += @CUDA_LIBS@
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav2_pthread_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav2_pthread_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lutil
pkdgrav2_pthread_CFLAGS += -I$(PYTHON_INC)
endif
if LUSTRE
pkdgrav2_pthread_LDADD += -llustreapi
endif
endif

tostd_SOURCES = tostd.c fio.c romberg.c
tdiff_SOURCES = tdiff.c fio.c romberg.c
