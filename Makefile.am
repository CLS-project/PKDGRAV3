SUBDIRS=openpa
ACLOCAL_AMFLAGS  = -I m4

SUFFIXES = .cu

.cu.o:
	$(NVCC) $(CFLAGS) -I. -DHAVE_CONFIG_H -lineinfo -o $@ -c $< $(NVCCFLAGS) -Imdl2/mpi -Imdl2 -Iopenpa/src


#if USE_GRAFIC
#SUBDIRS=grafic
#endif

EXTRA_DIST = server.py
EXTRA_PROGRAMS = pkdgrav3_null
TARGETS = null

if HAVE_PTHREAD
EXTRA_PROGRAMS	+= pkdgrav3_pthread
TARGETS	+= pthread
pthread		: pkdgrav3_pthread
endif

if HAVE_MPI
EXTRA_PROGRAMS	+= pkdgrav3_mpi
TARGETS	+= mpi

mpi:
	$(MAKE) CC=@MPICC@ FC=@MPIFC@ $(AM_MAKEFLAGS) pkdgrav3_mpi
endif

EXTRA_PROGRAMS += tostd tdiff hpb2fits lcp2tip

null	: pkdgrav3_null

all-local: $(TARGETS)

clean-local:
	rm -f $(EXTRA_PROGRAMS)

SOURCES = main.c tree.c analysis.c romberg.c bt.c
noinst_HEADERS = bt.h parameters.h grav.h walk.h cycle.h tipsydefs.h basetype.h vqsort.h
SOURCES        += cosmo.c master.c param.c pst.c smooth.c ewald.c
noinst_HEADERS += cosmo.h master.h param.h pst.h smooth.h ewald.h meval.h qeval.h
SOURCES        += moments.c intype.c outtype.c pkd.c smoothfcn.c rbtree.c fio.c listcomp.c illinois.c iomodule.c output.c
noinst_HEADERS += moments.h intype.h outtype.h pkd.h smoothfcn.h rbtree.h fio.h listcomp.h illinois.h iomodule.h output.h
if COOLING
#SOURCES        += cooling.c cooling_module.F90
#noinst_HEADERS += cooling.h
endif

SOURCES += healpix.c
noinst_HEADERS += healpix.h

SOURCES += hop.c fof.c group.c groupstats.c
noinst_HEADERS += hop.h fof.h group.h groupstats.h

SOURCES += $(MDL)/mdlbase.c
noinst_HEADERS += $(MDL)/mdlbase.h

SOURCES += grav2.c walk2.c
SOURCES += lst.c ilp.c ilc.c cl.c
noinst_HEADERS += lst.h ilp.h ilc.h cl.h
if USE_SIMD
SOURCES        += vmoments.cxx
AM_CFLAGS = @SIMD_CFLAGS@
AM_CXXFLAGS = @SIMD_CFLAGS@
noinst_HEADERS += simd.h sse_mathfun.h avx_mathfun.h
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
SOURCES += cudautil.cu cudapppc.cu cudaewald.cu
noinst_HEADERS += cudautil.h
endif

pkdgrav3_null_SOURCES = $(SOURCES) $(MDL)/null/mdl.c
noinst_HEADERS += $(MDL)/null/mdl.h
pkdgrav3_null_CFLAGS =  -I$(MDL)/null/ -I$(MDL) $(AM_CFLAGS)
pkdgrav3_null_LDADD = $(FLIBS)
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav3_null_LDADD += -L @EFENCE_PATH@
endif
pkdgrav3_null_LDADD += -lefence
endif
if USE_CUDA
pkdgrav3_null_LDADD += @CUDA_LIBS@
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav3_null_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav3_null_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lrt -lutil
pkdgrav3_null_CFLAGS += -I$(PYTHON_INC)
endif

if HAVE_MPI
pkdgrav3_mpi_LDADD = -lpthread
pkdgrav3_mpi_SOURCES = $(SOURCES) $(MDL)/mpi/mdl.c openpa/src/opa_primitives.c openpa/src/opa_queue.c
noinst_HEADERS += $(MDL)/mpi/mdl.h
pkdgrav3_mpi_CFLAGS = -I$(MDL)/mpi/ -I$(MDL) -Iopenpa/src $(AM_CFLAGS)
if FFTW
pkdgrav3_mpi_CFLAGS += $(GSL_CFLAGS)
pkdgrav3_mpi_LDADD += $(GSL_LIBS)
pkdgrav3_mpi_SOURCES += ic.c RngStream.c
noinst_HEADERS += ic.h RngStream.h
endif
if USE_LIBAIO
pkdgrav3_mpi_LDADD += -laio
endif
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav3_mpi_LDADD += -L @EFENCE_PATH@
endif
pkdgrav3_mpi_LDADD += -lefence
endif
if USE_ITT
if HAVE_ITT_PATH
pkdgrav3_mpi_LDADD += -L @ITT_PATH@/lib64
pkdgrav3_mpi_CFLAGS += -I @ITT_PATH@/include
endif
pkdgrav3_mpi_LDADD += -littnotify
endif
if USE_CUDA
pkdgrav3_mpi_LDADD += @CUDA_LIBS@ -lcudart -lnvToolsExt
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav3_mpi_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav3_mpi_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lutil
pkdgrav3_mpi_CFLAGS += -I$(PYTHON_INC)
endif
endif

if HAVE_PTHREAD
pkdgrav3_pthread_SOURCES = $(SOURCES) $(MDL)/pthread/mdl.c
noinst_HEADERS += $(MDL)/pthread/mdl.h
pkdgrav3_pthread_CFLAGS = -DUSE_PTHREAD -I$(MDL)/pthread/ -I$(MDL) $(AM_CFLAGS)
pkdgrav3_pthread_LDADD = -lpthread
if USE_EFENCE
if HAVE_EFENCE_PATH
pkdgrav3_pthread_LDADD += -L @EFENCE_PATH@
endif
pkdgrav3_pthread_LDADD += -lefence
endif
if USE_CUDA
pkdgrav3_pthread_LDADD += @CUDA_LIBS@ -lcudart -lnvToolsExt
endif
if USE_PYTHON
if HAVE_RDYNAMIC
pkdgrav3_pthread_LDFLAGS = -Xlinker -rdynamic 
endif
pkdgrav3_pthread_LDADD += -L$(PYTHON_LIB) -lpython$(PYTHON_VER) -ldl -lutil
pkdgrav3_pthread_CFLAGS += -I$(PYTHON_INC)
endif
endif

tostd_SOURCES = tostd.c fio.c romberg.c
tdiff_SOURCES = tdiff.c fio.c romberg.c
hpb2fits_SOURCES = hpb2fits.c
lcp2tip_SOURCES = lcp2tip.c

if HAVE_CFITSIO_PATH
hpb2fits_LDADD = -L@CFITSIO_PATH@/lib -lcfitsio -lm
hpb2fits_CFLAGS = -I @CFITSIO_PATH@/include
endif
