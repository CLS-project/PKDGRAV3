#
# Makefile for pkdgrav2
#

#MY_CC = gcc
#MY_CFLAGS = -O3 -Wall #-mcmodel=medium
MY_CC = pgcc
MY_CFLAGS = -fastsse -tp athlonxp
#MY_CC = icc
#MY_CFLAGS = -D__GNUC__ -D_REENTRANT -O3

EXE = pkdgrav2.$(MY_CC).null32
#EXE = pkdgrav2.$(MY_CC).mpi32
#EXE = pkdgrav2.$(MY_CC).mpi64

CODEDEF = -DCHANGESOFT
#CODEDEF = -DRELAXATION -DGROUPFIND 
#CODEDEF = -DCHANGESOFT -DSOFTLINEAR
#CODEDEF = -DCOLLISIONS -DSLIDING_PATCH
#CODEDEF = -DRELAXATION
#CODEDEF = -DDEBUG -HSHRINK

#ev6
#CC = ccc -arch ev67

#
#       NULL defines
#
NULL_MDL		= ../mdl/null
NULL_CFLAGS		= $(MY_CFLAGS) -I$(NULL_MDL) $(CODEDEF)

#NULL_LD_FLAGS	= -Wl,-s
NULL_XOBJ		= 
NULL_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       SGI defines
#
SGI_MDL			= ../mdl/mpi
SGI_CFLAGS		= -O2 -I$(SGI_MDL) $(CODEDEF) -mips4 -64 -r10000
SGI_LD_FLAGS	= -mips4 -64 -r10000
SGI_XOBJ		=
SGI_LIBMDL		= $(SGI_MDL)/mdl.o -lmpi -lm
SGI_MDL_CFLAGS	= -g2 -O2 -mips4 -64 -r10000

#
#       LINUX ScaMPI defines [SCALI]
#
SCA_DIR			= /opt/scali
SCA_MDL			= ../mdl/mpi
SCA_CFLAGS		=  -D_REENTRANT -O3 -I$(SCA_DIR)/include -I$(SCA_MDL) $(CODEDEF) 
SCA_LD_FLAGS	= -L$(SCA_DIR)/lib -Wl,-R/opt/scali/lib -lmpi
SCA_XOBJ        =
SCA_LIBMDL      = $(SCA_MDL)/mdl.o 
# SCA_LIBMDL      = -L$(SCA_DIR)/lib $(SCA_MDL)/mdl.o $CRT_BEGIN -lmpi $CRT_END -lm -I$(SCA_DIR)/include
SCA_MDL_CFLAGS  = $(SCA_CFLAGS)	

#
#       LINUX LAM MPI defines
#
LAM_DIR			= ../lam/lam-6.3.2
#LAM_DIR			= /net/lam-gcc
LAM_MDL			= ../mdl/mpi
#LAM_CFLAGS		= -O3 -malign-double -mstack-align-double -mpentiumpro -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX -I$(LAM_DIR)/include
LAM_CFLAGS		= -fast -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX -I$(LAM_DIR)/include
LAM_LD_FLAGS		=  
LAM_XOBJ                =
LAM_LIBMDL              = -L$(LAM_DIR)/lib $(LAM_MDL)/mdl.o -lmpi -ltstdio -lt -largs -ltrillium -ltstdio -lmpi++ -lm -I$(LAM_DIR)/include
#LAM_MDL_CFLAGS = -O3 -malign-double -mstack-align-double -mpentiumpro -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX  -I$(LAM_DIR)/include 
LAM_MDL_CFLAGS  = -fast -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX  -I$(LAM_DIR)/include 

#
#       MPI defines
#
SPX_MDL			= ../mdl/mpi
SPX_CFLAGS		= $(MY_CFLAGS) -I$(SPX_MDL) $(CODEDEF)
SPX_LD_FLAGS	=
SPX_XOBJ		= 
SPX_LIBMDL		= $(SPX_MDL)/mdl.o -lm
SPX_MDL_CFLAGS	= $(MY_CFLAGS)

#
#		PVM defines
#
PVMDIR	= $(PVM_ROOT)
PVMLIB	= $(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	= $(HOME)/pvm3/bin
XDIR	= $(BDIR)/$(PVM_ARCH)

PVM_MDL		= ../mdl/pvm
PVM_CFLAGS	= -O3 -I$(PVMDIR)/include -I$(PVM_MDL) $(CODEDEF)
#PVM_CFLAGS	= -mips4 -g -I$(PVMDIR)/include -I$(PVM_MDL) $(CODEDEF)
PVM_XOBJ	= 
PVM_LIBMDL	= $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm

#
#       PTHREAD defines
#
PTHREAD_MDL			= ../mdl/pthread
PTHREAD_CFLAGS		= -fast -I$(PTHREAD_MDL) $(CODEDEF)
PTHREAD_LD_FLAGS 	= 
PTHREAD_XOBJ		=
PTHREAD_LIBMDL 		= $(PTHREAD_MDL)/mdl.o -lm -lpthread

#
#       PTHREAD_SGI defines
#
PTHREAD_SGI_MDL			= ../mdl/pthread
PTHREAD_SGI_CFLAGS		= -O2 -D_REENTRANT -I$(PTHREAD_SGI_MDL) $(CODEDEF) -mips4 -64 -r10000
PTHREAD_SGI_LD_FLAGS 	= -mips4 -64 -r10000
PTHREAD_SGI_XOBJ		= 
PTHREAD_SGI_LIBMDL 		= $(PTHREAD_SGI_MDL)/mdl.o -lm -lpthread
PTHREAD_SGI_MDL_CFLAGS	= -O2 -mips4 -64 -r10000

#
#       T3D MPP defines
#
T3D_MDL	= ../mdl/mpp
RPC		= ../rpc

T3D_CFLAGS		= -O3 -g -DCRAY_T3D -I$(T3D_MDL) -I$(RPC) $(CODEDEF)
T3D_XOBJ		= hyperlib.o
T3D_LIBMDL		= -O3 -g -L $(RPC) $(T3D_MDL)/mdl.o -lmpi -lrpc -lm
T3D_LD_FLAGS	=

#
#       T3DMPI MPP defines
#
T3DMPI_MDL	= ../mdl/t3dmpi
RPC			= ../rpc

T3DMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3DMPI_MDL) -I$(RPC) $(CODEDEF)
T3DMPI_XOBJ		= hyperlib.o
T3DMPI_LIBMDL	= -O3 -L $(RPC) $(T3DMPI_MDL)/mdl.o -lmpi -lrpc -lm
T3DMPI_LD_FLAGS	=

#
#       T3EMPI MPP defines
#
T3EMPI_MDL		= ../mdl/t3empi
T3EMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3EMPI_MDL) -I$(RPC) $(CODEDEF)
T3EMPI_XOBJ		= hyperlib.o
T3EMPI_LIBMDL	= $(T3EMPI_MDL)/mdl.o -lmpi -lm
T3DMPI_LD_FLAGS	=

#
#       KSR1 defines
#
KSR_MDL			= ../mdl/ksr
KSR_CFLAGS		= -O2 -w2 -I$(KSR_MDL) $(CODEDEF)
KSR_LD_FLAGS	= -para
KSR_XOBJ		=
KSR_LIBMDL		= $(KSR_MDL)/mdl.o -lm -lrpc

OBJ	= main.o master.o param.o outtype.o pkd.o pst.o grav.o tree.o \
	  ewald.o walk.o smooth.o smoothfcn.o moments.o cosmo.o romberg.o

EXTRA_OBJ =

default:	
	@echo "Please tell me what architecture to make."
	@echo "Choices are null, sgi, pvm, pthread, spx, sca_mpi, t3d, t3dmpi, t3empi, ksr."

install:
	@echo "No installation rules."

clean:
	-rm -f $(OBJ) $(EXTRA_OBJ)

spotless: clean
	-rm -f $(EXE)

depend:
	makedepend -Y -- $(CODEDEF) -- *.c

$(XDIR):
	-mkdir $(BDIR)
	-mkdir $(XDIR)

null:
	cd $(NULL_MDL); rm *.o; make "CC=$(MY_CC)" "CFLAGS=$(NULL_CFLAGS)"
	make $(EXE) "CC=$(MY_CC)" "CFLAGS=$(NULL_CFLAGS)" "LD_FLAGS=$(NULL_LD_FLAGS)"\
		"MDL=$(NULL_MDL)" "XOBJ=$(NULL_XOBJ)" "LIBMDL=$(NULL_LIBMDL)"

sgi:
	cd $(SGI_MDL); make CC=cc "CC_FLAGS=$(SGI_MDL_CFLAGS)"
	make $(EXE) CC=cc "CFLAGS=$(SGI_CFLAGS)" "LD_FLAGS=$(SGI_LD_FLAGS)"\
		"MDL=$(SGI_MDL)" "XOBJ=$(SGI_XOBJ)" "LIBMDL=$(SGI_LIBMDL)"

pvm:
	cd $(PVM_MDL); aimk	
	make $(EXE) "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)"\
		"MDL=$(PVM_MDL)" "XOBJ=$(PVM_XOBJ)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f $(EXE) $(XDIR)

pthread:
	cd $(PTHREAD_MDL); rm *.o; make "CC=$(MY_CC)" "CFLAGS=$(CFLAGS)"
	make $(EXE) "CC=$(MY_CC)" "CFLAGS=$(PTHREAD_CFLAGS)" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

pthread_dec:
	cd $(PTHREAD_MDL); make CC=ccc
	make $(EXE) CC=ccc "CFLAGS=$(PTHREAD_CFLAGS) -fast -arch ev6" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

pthread_sgi:
	cd $(PTHREAD_MDL); make CC=cc "CC_FLAGS=$(PTHREAD_SGI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(PTHREAD_SGI_CFLAGS)" "LD_FLAGS=$(PTHREAD_SGI_LD_FLAGS)"\
		"MDL=$(PTHREAD_SGI_MDL)" "XOBJ=$(PTHREAD_SGI_XOBJ)" "LIBMDL=$(PTHREAD_SGI_LIBMDL)"

sca_mpi:
	cd $(SCA_MDL); make "CC=$(CC)" "CC_FLAGS=$(SCA_MDL_CFLAGS)"
	make $(EXE) "CC=$(CC)" "CFLAGS=$(SCA_CFLAGS)" "LD_FLAGS=$(SCA_LD_FLAGS)"\
		"MDL=$(SCA_MDL)" "XOBJ=$(SCA_XOBJ)" "LIBMDL=$(SCA_LIBMDL)"

lam_mpi:
	cd $(LAM_MDL); make CC=pgcc "CC_FLAGS=$(LAM_MDL_CFLAGS)"
	make $(EXE) CC=pgcc "CFLAGS=$(LAM_CFLAGS)" "LD_FLAGS=$(LAM_LD_FLAGS)"\
		"MDL=$(LAM_MDL)" "XOBJ=$(LAM_XOBJ)" "LIBMDL=$(LAM_LIBMDL)"


mpi: spx

spx:
	cd $(SPX_MDL); rm *o; make CC="mpicc -cc=$(MY_CC)" "CC_FLAGS=$(SPX_MDL_CFLAGS)"
	make $(EXE) CC="mpicc -cc=$(MY_CC)" "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)"\
		"MDL=$(SPX_MDL)" "XOBJ=$(SPX_XOBJ)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make $(EXE) "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)"\
		"MDL=$(T3D_MDL)" "XOBJ=$(T3D_XOBJ)" "LIBMDL=$(T3D_LIBMDL)"

t3dmpi:
	cd $(T3DMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3DMPI_CFLAGS)" "LD_FLAGS=$(T3DMPI_LD_FLAGS)"\
		"MDL=$(T3DMPI_MDL)" "XOBJ=$(T3DMPI_XOBJ)" "LIBMDL=$(T3DMPI_LIBMDL)"

t3empi:
	cd $(T3EMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3EMPI_CFLAGS)" "LD_FLAGS=$(T3EMPI_LD_FLAGS)"\
		"MDL=$(T3EMPI_MDL)" "XOBJ=$(T3EMPI_XOBJ)" "LIBMDL=$(T3EMPI_LIBMDL)"

ksr:
	cd $(KSR_MDL); make
	make $(EXE) "CFLAGS=$(KSR_CFLAGS)" "LD_FLAGS=$(KSR_LD_FLAGS)"\
		"MDL=$(KSR_MDL)" "XOBJ=$(KSR_XOBJ)" "LIBMDL=$(KSR_LIBMDL)"

$(EXE): $(OBJ) $(XOBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o $@ $(OBJ) $(XOBJ) $(LIBMDL)

$(OBJ) $(EXTRA_OBJ): Makefile

# DO NOT DELETE

collision.o: pkd.h floattype.h
ewald.o: ewald.h pkd.h floattype.h meval.h qeval.h
grav.o: pkd.h floattype.h grav.h meval.h qeval.h
main.o: pst.h pkd.h floattype.h smoothfcn.h master.h param.h
main.o: parameters.h outtype.h
master.o: master.h param.h pst.h pkd.h floattype.h smoothfcn.h
master.o: parameters.h tipsydefs.h outtype.h
outtype.o: pkd.h floattype.h outtype.h
param.o: param.h
pkd.o: pkd.h floattype.h ewald.h grav.h walk.h tipsydefs.h
pst.o: pst.h pkd.h floattype.h smoothfcn.h outtype.h smooth.h
qqsmooth.o: smooth.h pkd.h floattype.h smoothfcn.h
smooth.o: smooth.h pkd.h floattype.h smoothfcn.h
smoothfcn.o: smoothfcn.h pkd.h floattype.h
walk.o: walk.h pkd.h floattype.h




