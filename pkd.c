#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <errno.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#ifdef __linux__
#include <sys/resource.h>
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif
#include "cudautil.h"
#include "pkd.h"
#include "ewald.h"
#include "walk.h"
#include "grav.h"
#include "mdl.h"
#include "tipsydefs.h"
#include "outtype.h"
#include "parameters.h"
#include "cosmo.h"

double pkdGetTimer(PKD pkd,int iTimer) {
    return(pkd->ti[iTimer].sec);
    }

double pkdGetSystemTimer(PKD pkd,int iTimer) {
    return(pkd->ti[iTimer].system_sec);
    }

double pkdGetWallClockTimer(PKD pkd,int iTimer) {
    return(pkd->ti[iTimer].wallclock_sec);
    }


void pkdClearTimer(PKD pkd,int iTimer) {
   int i;

    if (iTimer >= 0) {
	pkd->ti[iTimer].sec = 0.0;
	pkd->ti[iTimer].system_sec = 0.0;
	pkd->ti[iTimer].wallclock_sec = 0.0;
	pkd->ti[iTimer].iActive = 0;
	}
    else {
	for (i=0;i<MAX_TIMERS;++i) {
	    pkd->ti[i].sec = 0.0;
	    pkd->ti[i].system_sec = 0.0;
	    pkd->ti[i].wallclock_sec = 0.0;
	    pkd->ti[i].iActive = 0;
	    }
	}
    }


void pkdStartTimer(PKD pkd,int iTimer) {
    struct timeval tv;

    pkd->ti[iTimer].iActive++;

    if (pkd->ti[iTimer].iActive == 1) {
#ifdef _MSC_VER
        FILETIME ft;
        uint64_t clock;
	GetSystemTimeAsFileTime(&ft);
	clock = ft.dwHighDateTime;
	clock <<= 32;
	clock |= ft.dwLowDateTime;
	/* clock is in 100 nano-second units */
	pkd->ti[iTimer].wallclock_stamp = clock / 10000000.0;
#else
	gettimeofday(&tv,NULL);
	pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
#endif
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
#ifdef __linux__
	{
	    struct rusage ru;

	    getrusage(RUSAGE_SELF,&ru);
	    pkd->ti[iTimer].system_stamp = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
	    }
#endif
	}
    }


void pkdStopTimer(PKD pkd,int iTimer) {
    double sec;
#ifdef _MSC_VER
    FILETIME ft;
    uint64_t clock;
#else
    struct timeval tv;
#endif
    sec = -pkd->ti[iTimer].stamp;
    pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
    sec += pkd->ti[iTimer].stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].sec += sec;

    sec = -pkd->ti[iTimer].wallclock_stamp;

#ifdef _MSC_VER
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    pkd->ti[iTimer].wallclock_stamp = clock / 10000000.0;
#else
    gettimeofday(&tv,NULL);
    pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
#endif
    sec += pkd->ti[iTimer].wallclock_stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].wallclock_sec += sec;

#ifdef __linux__
	{
	struct rusage ru;

	sec = -pkd->ti[iTimer].system_stamp;
	getrusage(RUSAGE_SELF,&ru);
	pkd->ti[iTimer].system_stamp = ((double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec);
	sec += pkd->ti[iTimer].system_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].system_sec += sec;
	}
#endif
    pkd->ti[iTimer].iActive--;
    }

/* Add a NODE structure: assume double alignment */
static int pkdNodeAddStruct(PKD pkd,int n) {
    int iOffset = pkd->iTreeNodeSize;
    mdlassert( pkd->mdl, pkd->kdNodeListPRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(double)-1)) == 0 );
    pkd->iTreeNodeSize += n;
    return iOffset;
    }
/* Add n doubles to the node structure */
static int pkdNodeAddDouble(PKD pkd,int n) {
    int iOffset = pkd->iTreeNodeSize;
    mdlassert( pkd->mdl, pkd->kdNodeListPRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(double)-1)) == 0 );
    pkd->iTreeNodeSize += sizeof(double) * n;
    return iOffset;
    }
/* Add n floats to the node structure */
static int pkdNodeAddFloat(PKD pkd,int n) {
    int iOffset = pkd->iTreeNodeSize;
    mdlassert( pkd->mdl, pkd->kdNodeListPRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(float)-1)) == 0 );
    pkd->iTreeNodeSize += sizeof(float) * n;
    return iOffset;
    }
/* Add n 64-bit integers to the node structure */
static int pkdNodeAddInt64(PKD pkd,int n) {
    int iOffset = pkd->iTreeNodeSize;
    mdlassert( pkd->mdl, pkd->kdNodeListPRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(int64_t)-1)) == 0 );
    pkd->iTreeNodeSize += sizeof(int64_t) * n;
    return iOffset;
    }
/* Add n 32-bit integers to the node structure */
static int pkdNodeAddInt32(PKD pkd,int n) {
    int iOffset = pkd->iTreeNodeSize;
    mdlassert( pkd->mdl, pkd->kdNodeListPRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(int32_t)-1)) == 0 );
    pkd->iTreeNodeSize += sizeof(int32_t) * n;
    return iOffset;
    }

/* Add a structure: assume double alignment */
static int pkdParticleAddStruct(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(double)-1)) == 0 );
    pkd->iParticleSize += n;
    return iOffset;
    }

/* Add n doubles to the particle structure */
static int pkdParticleAddDouble(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(double)-1)) == 0 );
    pkd->iParticleSize += sizeof(double) * n;
    return iOffset;
    }

/* Add n floats to the particle structure */
static int pkdParticleAddFloat(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(float)-1)) == 0 );
    pkd->iParticleSize += sizeof(float) * n;
    return iOffset;
    }

/* Add n 64-bit integers to the particle structure */
static int pkdParticleAddInt64(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(int64_t)-1)) == 0 );
    pkd->iParticleSize += sizeof(int64_t) * n;
    return iOffset;
    }

/* Add n 32-bit integers to the particle structure */
static int pkdParticleAddInt32(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(int32_t)-1)) == 0 );
    pkd->iParticleSize += sizeof(int32_t) * n;
    return iOffset;
    }

/* Add n 16-bit integers to the particle structure */
static int pkdParticleAddInt16(PKD pkd,int n) {
    int iOffset = pkd->iParticleSize;
    mdlassert( pkd->mdl, pkd->pStorePRIVATE == NULL );
    mdlassert( pkd->mdl, (iOffset & (sizeof(int16_t)-1)) == 0 );
    pkd->iParticleSize += sizeof(int16_t) * n;
    return iOffset;
    }

/* Extend the tree by adding more nodes */
void pkdExtendTree(PKD pkd) {
    if ( pkd->nTreeTiles >= (1<<pkd->nTreeBitsHi) ) {
	fprintf(stderr, "ERROR: insufficent nodes available in tree build"
	    "-- Increase nTreeBitsLo and/or nTreeBitsHi\n"
	    "nTreeBitsLo=%d nTreeBitsHi=%d\n",
	    pkd->nTreeBitsLo, pkd->nTreeBitsHi);
	assert( pkd->nTreeTiles < (1<<pkd->nTreeBitsHi) );
	}
    pkd->kdNodeListPRIVATE[pkd->nTreeTiles] = mdlMalloc(pkd->mdl,(1<<pkd->nTreeBitsLo)*pkd->iTreeNodeSize);
    mdlassert(pkd->mdl,pkd->kdNodeListPRIVATE[pkd->nTreeTiles] != NULL);
    ++pkd->nTreeTiles;
    pkd->nMaxNodes = (1<<pkd->nTreeBitsLo) * pkd->nTreeTiles;
    }

void pkdInitialize(
    PKD *ppkd,MDL mdl,int nStore,int nBucket,int nGroup,int nTreeBitsLo, int nTreeBitsHi,
    int iCacheSize,int iWorkQueueSize,int iCUDAQueueSize,FLOAT *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,
    uint64_t mMemoryModel, int nMaxDomainRungs) {
    PKD pkd;
    PARTICLE *p;
    uint32_t pi;
    int j,ism;

#define RANDOM_SEED 1
    srand(RANDOM_SEED);

    pkd = (PKD)SIMD_malloc(sizeof(struct pkdContext));
    mdlassert(mdl,pkd != NULL);
    pkd->mdl = mdl;
    pkd->idSelf = mdlSelf(mdl);
    pkd->nThreads = mdlThreads(mdl);
    pkd->kdNodeListPRIVATE = NULL;
    pkd->pStorePRIVATE = NULL;
    pkd->pStorePRIVATE2 = NULL;
    pkd->nStore = nStore;
    pkd->nLocal = 0;
    pkd->nDark = nDark;
    pkd->nGas = nGas;
    pkd->nStar = nStar;
    pkd->nRejects = 0;
    for (j=0;j<3;++j) {
	pkd->fPeriod[j] = fPeriod[j];
	}

    pkd->uMinRungActive  = 0;
    pkd->uMaxRungActive  = 255;
    pkd->uRungVeryActive = 255;

    pkd->psGroupTable.nGroups = 0;
    pkd->psGroupTable.pGroup = NULL;


    /*
    ** Calculate the amount of memory (size) of each particle.  This is the
    ** size of a base particle (PARTICLE), plus any extra fields as defined
    ** by the current memory model.  Fields need to be added in order of
    ** descending size (i.e., doubles & int64 and then float & int32)
    */
    pkd->iParticleSize = sizeof(PARTICLE);
    pkd->iTreeNodeSize = sizeof(KDN);

#ifndef INTEGER_POSITION
    pkd->oPosition = pkdParticleAddDouble(pkd,3);
#endif
    if ( mMemoryModel & PKD_MODEL_PARTICLE_ID )
	pkd->oParticleID = pkdParticleAddInt64(pkd,1);
    else
	pkd->oParticleID = 0;

    pkd->oVelocity = 0;
    if ( mMemoryModel & PKD_MODEL_VELOCITY ) {
	if (sizeof(vel_t) == sizeof(double)) {
	    pkd->oVelocity = pkdParticleAddDouble(pkd,3);
	    }
	}
#ifdef INTEGER_POSITION
    pkd->oPosition = pkdParticleAddInt32(pkd,3);
#endif
    if ( mMemoryModel & PKD_MODEL_RELAXATION )
	pkd->oRelaxation = pkdParticleAddDouble(pkd,1);
    else
	pkd->oRelaxation = 0;

    if ( mMemoryModel & PKD_MODEL_SPH )
	pkd->oSph = pkdParticleAddStruct(pkd,sizeof(SPHFIELDS));
    else
	pkd->oSph = 0;

    if ( mMemoryModel & PKD_MODEL_STAR )
	pkd->oStar = pkdParticleAddStruct(pkd,sizeof(STARFIELDS));
    else
	pkd->oStar = 0;

    if ( mMemoryModel & PKD_MODEL_VELSMOOTH )
	pkd->oVelSmooth = pkdParticleAddStruct(pkd,sizeof(VELSMOOTH));
    else
	pkd->oVelSmooth = 0;
    if ( mMemoryModel & PKD_MODEL_VELOCITY ) {
	if (sizeof(vel_t) == sizeof(float)) {
	    pkd->oVelocity = pkdParticleAddFloat(pkd,3);
	    }
	}
    if ( mMemoryModel & PKD_MODEL_ACCELERATION )
	pkd->oAcceleration = pkdParticleAddFloat(pkd,3);
    else
	pkd->oAcceleration = 0;

    if ( mMemoryModel & PKD_MODEL_POTENTIAL )
	pkd->oPotential = pkdParticleAddFloat(pkd,1);
    else
	pkd->oPotential = 0;

    if ( mMemoryModel & PKD_MODEL_MASS )
	pkd->oMass = pkdParticleAddFloat(pkd,1);
    else
	pkd->oMass = 0;

    if ( mMemoryModel & PKD_MODEL_SOFTENING )
	pkd->oSoft = pkdParticleAddFloat(pkd,1);
    else
	pkd->oSoft = 0;

    if ( mMemoryModel & PKD_MODEL_SPH )
	pkd->oBall = pkdParticleAddFloat(pkd,1);
    else pkd->oBall = 0;
    if ( mMemoryModel & PKD_MODEL_SPH )
	pkd->oDensity = pkdParticleAddFloat(pkd,1);
    else pkd->oDensity = 0;

    if ( mMemoryModel & PKD_MODEL_GROUPS ) {
	pkd->oGroup = pkdParticleAddInt32(pkd,1);
	}
    else {
	pkd->oGroup = 0;
	}

    if ( mMemoryModel & PKD_MODEL_RUNGDEST ) {
	pkd->nMaxDomainRungs = nMaxDomainRungs;
	pkd->oRungDest = pkdParticleAddInt16(pkd,pkd->nMaxDomainRungs);
	}
    else {
	pkd->nMaxDomainRungs = 0;
	pkd->oRungDest = 0;
	}

    /*
    ** Tree node memory models
    */
    if ( mMemoryModel & PKD_MODEL_NODE_BND ) {
        pkd->oNodeBnd  = pkdNodeAddStruct(pkd,sizeof(BND));
    }
    else {
	pkd->oNodeBnd  = 0;
    }

    if ( mMemoryModel & PKD_MODEL_NODE_VBND ) {
	pkd->oNodeVBnd  = pkdNodeAddStruct(pkd,sizeof(BND));
    }
    else {
	pkd->oNodeVBnd = 0;
    }

    pkd->oNodeVelocity = 0;
    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(double))
	    pkd->oNodeVelocity = pkdNodeAddDouble(pkd,3);
    /*
    ** Three extra bounds are required by the fast gas SPH code.
    */
    if ( mMemoryModel & PKD_MODEL_NODE_SPHBNDS ) {
	pkd->oNodeSphBounds = pkdNodeAddStruct(pkd,sizeof(SPHBNDS));
    }
    else
	pkd->oNodeSphBounds = 0;

    if ( mMemoryModel & PKD_MODEL_NODE_MOMENT )
	pkd->oNodeMom = pkdNodeAddStruct(pkd,sizeof(FMOMR));
    else
	pkd->oNodeMom = 0;

    /* The acceleration is required for the new time step criteria */
    if ( mMemoryModel & PKD_MODEL_NODE_ACCEL )
	pkd->oNodeAcceleration = pkdNodeAddFloat(pkd,3);
    else
	pkd->oNodeAcceleration = 0;

    if ( (mMemoryModel & PKD_MODEL_NODE_VEL) && sizeof(vel_t) == sizeof(float))
	    pkd->oNodeVelocity = pkdNodeAddFloat(pkd,3);
    /*
    ** N.B.: Update pkdMaxNodeSize in pkd.h if you add fields.  We need to
    **       know the size of a node when setting up the pst.
    */
    assert(pkdNodeSize(pkd) > 0);
    if (pkdNodeSize(pkd) > pkdMaxNodeSize()) {
	fprintf(stderr, "Node size is too large. Node size=%llu, max node size=%llu\n", pkdNodeSize(pkd), pkdMaxNodeSize());
	}
    assert(pkdNodeSize(pkd)<=pkdMaxNodeSize());

    /*
    ** Allocate the main particle store.
    ** Need to use mdlMalloc() since the particles will need to be
    ** visible to all other processors thru mdlAcquire() later on.
    **
    ** We need one EXTRA storage location at the very end to use for
    ** calculating acceleration on arbitrary positions in space, for example
    ** determining the force on the sun. The easiest way to do this is to
    ** allocate one hidden particle, which won't interfere with the rest of
    ** the code (hopefully). pkd->pStore[pkd->nStore] is this particle.
    **
    ** We also allocate a temporary particle used for swapping.  We need to do
    ** this now because the outside world can no longer know the size of a
    ** particle.
    */
    pkd->iParticleSize = (pkd->iParticleSize + sizeof(double) - 1 ) & ~(sizeof(double)-1);
    pkd->pStorePRIVATE = mdlMallocArray(pkd->mdl,nStore+1,pkdParticleSize(pkd));
    mdlassert(mdl,pkd->pStorePRIVATE != NULL);
    if ( mMemoryModel & PKD_MODEL_RUNGDEST ) {
	pkd->pStorePRIVATE2 = mdlMalloc(pkd->mdl,(nStore+1)*pkdParticleSize(pkd));
	mdlassert(mdl,pkd->pStorePRIVATE2 != NULL);
	}
    pkd->pTempPRIVATE = malloc(pkdParticleSize(pkd));
    mdlassert(mdl,pkd->pTempPRIVATE != NULL);

    /* Create a type for our particle -- we use an opaque type */
#ifdef xMPI_VERSION
    mdlTypeContiguous(pkd->mdl, pkd->iParticleSize, MDL_BYTE, &pkd->typeParticle );
    mdlTypeCommit(pkd->mdl,&pkd->typeParticle);
#endif

#ifdef MDL_CACHE_SIZE
    if ( iCacheSize > 0 ) mdlSetCacheSize(pkd->mdl,iCacheSize);
#endif
    // This is cheeserific - chooses the largest specified

#ifdef USE_CUDA
    mdlSetCudaBufferSize(pkd->mdl,MAX_EWALD_PARTICLES*sizeof(double)*4,MAX_EWALD_PARTICLES*sizeof(double)*5);
    mdlSetCudaBufferSize(pkd->mdl,PP_CUDA_MEMORY_LIMIT,PP_CUDA_MEMORY_LIMIT);
#endif
    mdlSetWorkQueueSize(pkd->mdl,iWorkQueueSize,iCUDAQueueSize);
    /*
    ** Initialize neighbor list pointer to NULL if present.
    */
    if (pkd->oSph) {
	for (pi=0;pi<(pkd->nStore+1);++pi) {
	    p = pkdParticle(pkd,pi);
	    *pkd_pNeighborList(pkd,p) = NULL;
	}
    }

    /*
    ** We support up to 256 classes
    */
    pkd->pClass = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlassert(mdl,pkd->pClass != NULL);
    for (j=0;j<PKD_MAX_CLASSES;j++) {
	pkd->pClass[j].fMass = pkd->pClass[j].fSoft = -1.0;
	pkd->pClass[j].eSpecies = FIO_SPECIES_LAST;
	}
    pkd->nClasses = 0;

    pkd->fSoftFix = -1.0;
    pkd->fSoftFac = 1.0;
    pkd->fSoftMax = HUGE;
    /*
    ** Now we setup the node storage for the tree.  This storage is no longer
    ** continguous as the MDL now supports non-contiguous arrays.  We allocate
    ** a single "tile" for the tree.  If this is not sufficient, then additional
    ** tiles are allocated dynamically.  The default parameters allow for 2^32
    ** nodes total which is the integer limit anyway.
    */
    pkd->iTreeNodeSize = (pkd->iTreeNodeSize + sizeof(double) - 1 ) & ~(sizeof(double)-1);
    pkd->nTreeBitsLo = nTreeBitsLo;
    pkd->nTreeBitsHi = nTreeBitsHi;
    pkd->iTreeMask = (1<<pkd->nTreeBitsLo) - 1;
    pkd->kdNodeListPRIVATE = mdlMalloc(pkd->mdl,(1<<pkd->nTreeBitsHi)*sizeof(KDN *));
    mdlassert(mdl,pkd->kdNodeListPRIVATE != NULL);
    pkd->kdNodeListPRIVATE[0] = mdlMalloc(pkd->mdl,(1<<pkd->nTreeBitsLo)*pkd->iTreeNodeSize);
    mdlassert(mdl,pkd->kdNodeListPRIVATE[0] != NULL);
    pkd->nTreeTiles = 1;
    pkd->nMaxNodes = (1<<pkd->nTreeBitsLo) * pkd->nTreeTiles;
    /*
    ** pLite particles are also allocated and are quicker when sorting particle
    ** type operations such as tree building and domain decomposition are being
    ** performed.
    */
    pkd->pLite = mdlMallocArray(pkd->mdl,nStore+1,EPHEMERAL_BYTES);
    mdlassert(mdl,pkd->pLite != NULL);
    pkd->nNodes = 0;
    /*
    ** Ewald stuff!
    */
    pkd->ew.nMaxEwhLoop = 0;
    *ppkd = pkd;
    /*
    ** Tree walk stuff.
    */
    ilpInitialize(&pkd->ilp);
    ilcInitialize(&pkd->ilc);
    /*
    ** Allocate Checklist.
    */
    pkd->clFreeList.list = NULL;
    pkd->clFreeList.nRefs = 0;
    pkd->clFreeList.nTiles = 0;
    clInitialize(&pkd->cl,&pkd->clFreeList);
    clInitialize(&pkd->clNew,&pkd->clFreeList);
    /*
    ** Allocate the stack.
    */
    pkd->nMaxStack = 30;
    pkd->S = malloc(pkd->nMaxStack*sizeof(CSTACK));
    assert(pkd->S != NULL);
    for (ism=0;ism<pkd->nMaxStack;++ism) {
	clInitialize(&pkd->S[ism].cl,&pkd->clFreeList);
	}
    /*
    ** Allocate initial particle pointer arrays for active/inactive particles.
    */
    pkd->nMaxBucketActive = 1000;
    pkd->piActive = malloc(pkd->nMaxBucketActive*sizeof(PARTICLE *));
    mdlassert(mdl,pkd->piActive != NULL);
    pkd->piInactive = malloc(pkd->nMaxBucketActive*sizeof(PARTICLE *));
    mdlassert(mdl,pkd->piInactive != NULL);

    pkd->profileBins = NULL;
    pkd->groupBin = NULL;

    pkd->grid = NULL;
    pkd->gridData = NULL;

    pkd->tmpHopGroups = NULL;
    pkd->hopGroups = NULL;
    pkd->hopNumRoots = NULL;
    pkd->hopRootIndex = NULL;
    pkd->hopRoots = NULL;

#ifdef COOLING
    pkd->Cool = CoolInit();
#endif
    assert(pkdNodeSize(pkd) > 0);

#ifdef xNOxUSE_CUDA
	{
	int sizeILP = sizeof(ILP_BLK)*pkd->ilp->lst.nBlocksPerTile;
	int sizeILC = sizeof(ILC_BLK)*pkd->ilc->lst.nBlocksPerTile;
	pkd->cudaCtx = CUDA_initialize(mdlCore(pkd->mdl),
	    iCUDAQueueSize,
	    sizeILP>sizeILC ? sizeILP : sizeILC,
	    nGroup*sizeof(PINFOIN),
	    nGroup*sizeof(PINFOOUT)
	    * (pkd->ilp->lst.nBlocksPerTile>pkd->ilc->lst.nBlocksPerTile
		? pkd->ilp->lst.nBlocksPerTile : pkd->ilc->lst.nBlocksPerTile) );
	}
#endif
    }


void pkdFinish(PKD pkd) {
    PARTICLE *p;
    char **ppCList;
    uint32_t pi;
    int ism;
    int i;

#ifdef xNOxUSE_CUDA
    CUDA_finish(pkd->cudaCtx);
#endif

    if (pkd->kdNodeListPRIVATE) {
	/*
	** Close caching space and free up nodes.
	*/
	if (pkd->nNodes > 0)
	    mdlFinishCache(pkd->mdl,CID_CELL);
	for( i=0; i<pkd->nTreeTiles; i++)
	    mdlFree(pkd->mdl,pkd->kdNodeListPRIVATE[i]);
	mdlFree(pkd->mdl,pkd->kdNodeListPRIVATE);
	}
    /*
    ** Free Interaction lists.
    */
    ilpFinish(pkd->ilp);
    ilcFinish(pkd->ilc);
    /*
    ** Free checklist.
    */
    clFinish(pkd->cl);
    clFinish(pkd->clNew);
    /*
    ** Free Stack.
    */
    for (ism=0;ism<pkd->nMaxStack;++ism) {
	clFinish(pkd->S[ism].cl);
	}
    free(pkd->S);
    if (pkd->ew.nMaxEwhLoop) {
	SIMD_free(pkd->ewt.hx.f);
	SIMD_free(pkd->ewt.hy.f);
	SIMD_free(pkd->ewt.hz.f);
	SIMD_free(pkd->ewt.hCfac.f);
	SIMD_free(pkd->ewt.hSfac.f);
	}

    free(pkd->pClass);
    /*
    ** Free any neighbor lists that were left hanging around.
    */
    if (pkd->oSph) {
	for (pi=0;pi<(pkd->nStore+1);++pi) {
	    p = pkdParticle(pkd,pi);
	    ppCList = pkd_pNeighborList(pkd,p);
	    if (*ppCList) {
		free(*ppCList);
		*ppCList = NULL;
	    }
	}
    }
#ifdef xMPI_VERSION
    mdlTypeFree(pkd->mdl,&pkd->typeParticle);
#endif
    mdlFreeArray(pkd->mdl,pkd->pStorePRIVATE);
    if (pkd->pStorePRIVATE2)
	mdlFree(pkd->mdl,pkd->pStorePRIVATE2);
    free(pkd->pTempPRIVATE);
    mdlFreeArray(pkd->mdl,pkd->pLite);
    free(pkd->piActive);
    free(pkd->piInactive);
    csmFinish(pkd->param.csm);
    SIMD_free(pkd);
    }

size_t pkdClCount(PKD pkd) {
    size_t nCount = clCount(pkd->cl);
    int i;
    for(i=0; i<pkd->nMaxStack; ++i)
	nCount += clCount(pkd->S[i].cl);
    return nCount;
    }

size_t pkdClMemory(PKD pkd) {
    return clMemory(pkd->cl);
    }

size_t pkdIlpMemory(PKD pkd) {
    return ilpMemory(pkd->ilp);
    }

size_t pkdIlcMemory(PKD pkd) {
    return ilcMemory(pkd->ilc);
    }

size_t pkdTreeMemory(PKD pkd) {
    return pkd->nTreeTiles * (1<<pkd->nTreeBitsLo) * pkd->iTreeNodeSize;
    }

static void getClass( PKD pkd, float fMass, float fSoft, FIO_SPECIES eSpecies, PARTICLE *p ) {
    int i;

    if ( pkd->oMass ) {
	float *pMass = pkdField(p,pkd->oMass);
	*pMass = fMass;
	fMass = 0.0;
	}
    if ( pkd->oSoft ) {
	float *pSoft = pkdField(p,pkd->oSoft);
	*pSoft = fSoft;
	fSoft = 0.0;
	}
    /* NOTE: The above can both be true, in which case a "zero" class is recorded */
    /* NOTE: Species is always part of the class table, so there will be at least one class per species */

    /* TODO: This is a linear search which is fine for a small number of classes */
    for ( i=0; i<pkd->nClasses; i++ )
	if ( pkd->pClass[i].fMass == fMass && pkd->pClass[i].fSoft == fSoft && pkd->pClass[i].eSpecies==eSpecies )
	    break;

    if ( i == pkd->nClasses ) {
	assert( pkd->nClasses < PKD_MAX_CLASSES );
	i = pkd->nClasses++;
	pkd->pClass[i].fSoft    = fSoft;
	pkd->pClass[i].fMass    = fMass;
	pkd->pClass[i].eSpecies = eSpecies;
	}
    p->iClass = i;
    }

int pkdGetClasses( PKD pkd, int nMax, PARTCLASS *pClass ) {
    int i;
    for ( i=0; i<pkd->nClasses; i++ )
	pClass[i] = pkd->pClass[i];
    return pkd->nClasses;
    }

void pkdSetClasses( PKD pkd, int n, PARTCLASS *pClass, int bUpdate ) {
    uint8_t map[PKD_MAX_CLASSES];
    PARTICLE *p;
    int i,j;

    if ( bUpdate ) {
	/* Build a map from the old class to the new class */
	assert( n >= pkd->nClasses );
	for ( i=0; i<pkd->nClasses; i++ ) {
	    for ( j=0; j<n; j++ )
		if ( pClass[j].fMass==pkd->pClass[i].fMass && pClass[j].fSoft==pkd->pClass[i].fSoft && pClass[j].eSpecies==pkd->pClass[i].eSpecies )
		    break;
	    assert(j<n);
	    map[i] = j;
	    }

	/* Now update the class with the new value */
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    assert( p->iClass <= pkd->nClasses );
	    p->iClass = map[p->iClass];
	    }
	}

    /* Finally, set the new class table */
    for ( i=0; i<n; i++ )
	pkd->pClass[i] = pClass[i];
    pkd->nClasses = n;
    }

void pkdSeek(PKD pkd,FILE *fp,uint64_t nStart,int bStandard,int bDoublePos) {
#ifndef HAVE_FSEEKO
    off_t MAX_OFFSET = 2147483640;
    int iErr;
#endif
    off_t lStart;

    /*
    ** Seek according to true XDR size structures when bStandard is true.
    ** This may be a bit dicey, but it should work as long
    ** as no one changes the tipsy binary format!
    */
    if (bStandard) lStart = 32;
    else lStart = sizeof(struct dump);
    if (nStart > pkd->nGas) {
	if (bStandard) lStart += pkd->nGas*(bDoublePos?60:48);
	else lStart += pkd->nGas*sizeof(struct gas_particle);
	nStart -= pkd->nGas;
	if (nStart > pkd->nDark) {
	    if (bStandard) lStart += pkd->nDark*(bDoublePos?48:36);
	    else lStart += pkd->nDark*sizeof(struct dark_particle);
	    nStart -= pkd->nDark;
	    if (bStandard) lStart += nStart*(bDoublePos?56:44);
	    else lStart += nStart*sizeof(struct star_particle);
	    }
	else {
	    if (bStandard) lStart += nStart*(bDoublePos?48:36);
	    else lStart += nStart*sizeof(struct dark_particle);
	    }
	}
    else {
	if (bStandard) lStart += nStart*(bDoublePos?60:48);
	else lStart += nStart*sizeof(struct gas_particle);
	}

#ifdef HAVE_FSEEKO
    fseeko(fp,lStart,SEEK_SET);
#else
    /*fseek fails for offsets >= 2**31; this is an ugly workaround;*/
    if (lStart > MAX_OFFSET) {
	iErr = fseek(fp,0,SEEK_SET);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	while (lStart > MAX_OFFSET) {
	    fseek(fp,MAX_OFFSET,SEEK_CUR);
	    lStart -= MAX_OFFSET;
	    }
	iErr = fseek(fp,lStart,SEEK_CUR);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	}
    else {
	iErr = fseek(fp,lStart,SEEK_SET);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	}
#endif
    }


#ifdef USE_GRAFIC
void pkdGenerateIC(PKD pkd, GRAFICCTX gctx,  int iDim,
		   double fSoft, double fMass, int bComove) {
    PARTICLE *p;
    int i, j, k, d, pi, n1, n2, n3;
    double dx;
    double dvFac, a;
    vel_t *v;

    graficGenerate(gctx, iDim, 1 );

    pkd->nLocal = pkd->nActive = graficGetLocal(gctx);

    pi = 0;
    n1 = graficGetLocalDim(gctx,1);
    n2 = graficGetLocalDim(gctx,2);
    n3 = graficGetLocalDim(gctx,3);

    dx = 1.0 / n1;
    d = iDim - 1;
    a = graficGetExpansionFactor(gctx);
    dvFac = bComove ? a*a : 1.0;

    for ( i=0; i<n1; i++ ) {
	for ( j=0; j<n2; j++ ) {
	    for ( k=0; k<n3; k++ ) {
		p = pkdParticle(pkd,pi);
		v = pkdVel(pkd,p);
		p->uRung = p->uNewRung = 0;
		p->bSrcActive = p->bDstActive = 1;
		p->fDensity = 0.0;
		p->fBall = 0.0;
		/*
		** Clear the accelerations so that the timestepping calculations
		** do not get funny uninitialized values!
		*/
		p->a[0] = p->a[1] = p->a[2] = 0.0;
		p->r[d] = graficGetPosition(gctx,i,j,k,d) - 0.5;
		if ( p->r[d] < -0.5 ) p->r[d] += 1.0;
		if ( p->r[d] >= 0.5 ) p->r[d] -= 1.0;
		assert( p->r[d] >= -0.5 && p->r[d] < 0.5 );
		v[d] = graficGetVelocity(gctx,i,j,k) * dvFac;
		getClass(pkd,fMass,fSoft,p);
		p->iOrder = 0; /* FIXME */
		p->iClass = 0;
		pi++;
		}
	    }
	}

    assert( pi == pkd->nLocal );
    }

#endif

void pkdReadFIO(PKD pkd,FIO fio,uint64_t iFirst,int nLocal,double dvFac, double dTuFac) {
    int i,j;
    PARTICLE *p;
    STARFIELDS *pStar;
    SPHFIELDS *pSph;
    float *pPot, dummypot;
    double r[3];
    double vel[3];
    float fMass, fSoft,fDensity;
    FIO_SPECIES eSpecies;
    uint64_t iParticleID;

    mdlassert(pkd->mdl,fio != NULL);

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Read");
     __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    if (pkd->oStar) {
	/* Make sure star class established -- how do all procs know of these classes? How do we ensure they agree on the class identifiers? */
	p = pkdParticle(pkd,pkd->nLocal);
	getClass(pkd,0,0,FIO_SPECIES_STAR,p);
	}

    fioSeek(fio,iFirst,FIO_SPECIES_ALL);
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,pkd->nLocal+i);
	/*
	** General initialization.
	*/
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	pkdSetDensity(pkd,p,0.0);
	if (pkd->oBall) pkdSetBall(pkd,p,0.0);
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	if ( pkd->oAcceleration ) {
	    float *a = pkdAccel(pkd,p);
	    for (j=0;j<3;++j) a[j] = 0.0;
	    }
	if ( pkd->oPotential) pPot = pkdPot(pkd,p);
	else pPot = &dummypot;

	/* Initialize SPH fields if present */
	if (pkd->oSph) {
	    pSph = pkdField(p,pkd->oSph);
	    pSph->u = pSph->uPred = pSph->uDot = pSph->c = pSph->divv = pSph->BalsaraSwitch
		= pSph->fMetals = pSph->diff = pSph->fMetalsPred = pSph->fMetalsDot = 0.0;
	    }
	else pSph = NULL;

	/* Initialize Star fields if present */
	if (pkd->oStar) {
	    pStar = pkdField(p,pkd->oStar);
	    pStar->fTimer = 0;
/*	    pStar->iGasOrder = IORDERMAX;*/
	    }
	else pStar = NULL;

	eSpecies = fioSpecies(fio);
	switch(eSpecies) {
	case FIO_SPECIES_SPH:
	    assert(pSph); /* JW: Could convert to dark ... */
	    assert(dTuFac>0.0);
	    fioReadSph(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,
			     &fDensity/*?*/,&pSph->u,&pSph->fMetals);
	    pkdSetDensity(pkd,p,fDensity);
	    pSph->u *= dTuFac; /* Can't do precise conversion until density known */
	    pSph->uPred = pSph->u;
	    pSph->fMetalsPred = pSph->fMetals;
	    pSph->vPred[0] = vel[0]*dvFac;
	    pSph->vPred[1] = vel[1]*dvFac;
	    pSph->vPred[2] = vel[2]*dvFac; /* density, divv, BalsaraSwitch, c set in smooth */
	    break;
	case FIO_SPECIES_DARK:
	    fioReadDark(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,&fDensity);
	    pkdSetDensity(pkd,p,fDensity);
	    break;
	case FIO_SPECIES_STAR:
	    assert(pStar && pSph);
	    fioReadStar(fio,&iParticleID,r,vel,&fMass,&fSoft,pPot,&fDensity,
			      &pSph->fMetals,&pStar->fTimer);
	    pkdSetDensity(pkd,p,fDensity);
	    pSph->vPred[0] = vel[0]*dvFac;
	    pSph->vPred[1] = vel[1]*dvFac;
	    pSph->vPred[2] = vel[2]*dvFac;
	    break;
	default:
	    fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
	    assert(0);
	    }

	for (j=0;j<3;++j) pkdSetPos(pkd,p,j,r[j]);
	if (pkd->oVelocity) {
	    for (j=0;j<3;++j) pkdVel(pkd,p)[j] = vel[j]*dvFac;
	    }

	p->iOrder = iFirst++;
	if (pkd->oParticleID) *pkdParticleID(pkd,p) = iParticleID;

	getClass(pkd,fMass,fSoft,eSpecies,p);
	}
    
    pkd->nLocal += nLocal;
    pkd->nActive += nLocal;

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
    }

void pkdCalcBound(PKD pkd,BND *pbnd) {
    double r[3],dMin[3],dMax[3];
    PARTICLE *p;
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->nLocal > 0);
    p = pkdParticle(pkd,i);
    for (j=0;j<3;++j) {
	dMin[j] = dMax[j] = pkdPos(pkd,p,j);
	}
    for (++i;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) r[j] = pkdPos(pkd,p,j);
	pkdMinMax(r,dMin,dMax);
	}
    for (j=0;j<3;++j) {
	pbnd->fCenter[j] = pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pbnd->fMax[j] = pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
    }

void pkdCalcVBound(PKD pkd,BND *pbnd) {
    double dMin[3],dMax[3];
    PARTICLE *p;
    vel_t *v;
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->nLocal > 0);
    p = pkdParticle(pkd,i);
    v = pkdVel(pkd,p);
    for (j=0;j<3;++j) {
	dMin[j] = v[j];
	dMax[j] = v[j];
	}
    for (++i;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	v = pkdVel(pkd,p);
	pkdMinMax(v,dMin,dMax);
	}
    for (j=0;j<3;++j) {
	pbnd->fCenter[j] = pkd->vbnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pbnd->fMax[j] = pkd->vbnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
    }


void pkdEnforcePeriodic(PKD pkd,BND *pbnd) {
    PARTICLE *p;
    double r;
    int i,j;

    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) {
	    r = pkdPos(pkd,p,j);
	    if (r < pbnd->fCenter[j] - pbnd->fMax[j]) r += 2*pbnd->fMax[j];
	    else if (r >= pbnd->fCenter[j] + pbnd->fMax[j]) r -= 2*pbnd->fMax[j];
	    pkdSetPos(pkd,p,j,r);
	    /*
	    ** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the 
	    ** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
	    ** to the specified periodic box in a gross way.
	    */
	    mdlassert(pkd->mdl,((r >= pbnd->fCenter[j] - pbnd->fMax[j])&&
	    (r < pbnd->fCenter[j] + pbnd->fMax[j])));
	    }
	}
    }


/*
** x and y must have range [1,2) !
** returns key in the range [0,2^63-1]
*/
uint64_t hilbert2d(float x,float y) {
    uint64_t s = 0;
    uint32_t m,ux,uy,ut;

    union {
	float    f;
	uint32_t u;
	} punner;

    punner.f = x; ux = punner.u >> 2;
    punner.f = y; uy = punner.u >> 2;
   
    m = 0x00100000;

    while (m) {
	s = s << 2;
	if (ux&m) {
	    if (uy&m) {
		s |= 2;
		}
	    else {
		ut = ux;
		ux = ~uy;
		uy = ~ut;
		s |= 3;
		}
	    }
	else {
	    if (uy&m) {
		s |= 1;
		}
	    else {
		ut = ux;
		ux = uy;
		uy = ut;
		}
	    }
	m = m >> 1;
	}
    return s;
    }

/*
** x, y and z must have range [1,2) !
*/
uint64_t hilbert3d(float x,float y,float z) {
    uint64_t s = 0;
    uint32_t m,ux,uy,uz,ut;

    union {
	float    f;
	uint32_t u;
	} punner;

    punner.f = x; ux = punner.u >> 2;
    punner.f = y; uy = punner.u >> 2;
    punner.f = z; uz = punner.u >> 2;
    /* Or: ux = (uint32_t)((x-1.0f) * 0x00200000)*/

    m = 0x00100000;
    while (m) {
	s = s << 3;

	if (ux&m) {
	    if (uy&m) {
		if (uz&m) {
		    ut = ux;
		    ux = uy;
		    uy = ~uz;
		    uz = ~ut;
		    s |= 5;
		    }
		else {
		    ut = uz;
		    uz = ux;
		    ux = uy;
		    uy = ut;
		    s |= 2;
		    }
		}
	    else {
		ux = ~ux;
		uy = ~uy;
		if (uz&m) {
		    s |= 4;
		    }
		else {
		    s |= 3;
		    }
		}
	    }
	else {
	    if (uy&m) {
		if (uz&m) {
		    ut = ux;
		    ux = uy;
		    uy = ~uz;
		    uz = ~ut;
		    s |= 6;
		    }
		else {
		    ut = uz;
		    uz = ux;
		    ux = uy;
		    uy = ut;
		    s |= 1;
		    }
		}
	    else {
		if (uz&m) {
		    ut = uy;
		    uy = ux;
		    ux = ~uz;
		    uz = ~ut;
		    s |= 7;
		    }
		else {
		    ut = uy;
		    uy = ux;
		    ux = uz;
		    uz = ut;
		    s |= 0;
		    }
		}
	    }
	m = m >> 1;
	}
    return s;
    }

#ifdef MPI_VERSION
typedef struct {
    int64_t lKey;
    int32_t  i;
    } PLITEDD;

typedef struct {
    int64_t r,s,t;
    double fr, fs;
    } ILLINOISDD;

static int cmpPeanoHilbert(const void *pva,const void *pvb) {
    PLITEDD *pa = (PLITEDD *)pva;
    PLITEDD *pb = (PLITEDD *)pvb;
    if (pa->lKey < pb->lKey ) return -1;
    else if (pa->lKey > pb->lKey ) return 1;
    else return 0;
    }

/* Handy macros for: rung index, active index, inactive index */
#define RI(domain,rung) ((domain)*nRungs+(rung))
#define AI(domain,rung) (RI(domain,rung)*2)
#define II(domain,rung) (AI(domain,rung)+1)

/*
** Methods:
**  0 - Split particles by active, but use the same domains for inactive
**  1 - Split particles by active and by inactive
**  2 - Split particles by active, inactive on every rung
**  3 - Split particles by active on every rung
**  4 - Split on all rungs
*/
void pkdPeanoHilbertDecomp(PKD pkd, int nRungs, int iMethod) {
    PLITEDD *pl = (PLITEDD *)pkd->pLite;
    PARTICLE *p;
    float x,y,z;
    int nDomains = mdlThreads(pkd->mdl);
    int nThreads = mdlThreads(pkd->mdl);
    int nLocal   = pkd->nLocal;
    int nID      = mdlSelf(pkd->mdl);
    int i, iRung, iRungMin, iRungMax, iActiveRung, iKey;
    int64_t lKey;
    int L, U, M;
    int32_t *nActiveBelow;
    int64_t *nRungTotal;
    int64_t *splitKeys, *nTarget;
    int64_t *nBelowLocal, *nBelowGlobal;
    int64_t *myBound;
    int *counts, *piIndex;
    int iIter, nIter;

    /* Root finder state */
    int64_t *r, *s, *t;
    double *fr, *fs, ft;

    /* Root finder state - active/inactive for each rung */
    r  = malloc(nRungs * 2 * sizeof(*r));
    s  = malloc(nRungs * 2 * sizeof(*s));
    t  = malloc(nRungs * 2 * sizeof(*t));
    fr = malloc(nRungs * 2 * sizeof(*fr));
    fs = malloc(nRungs * 2 * sizeof(*fs));


    /* Number of particles active on the given rung for each particle */
    nActiveBelow = malloc((nLocal+1) * nRungs*2 * sizeof(*nActiveBelow));
    assert(nActiveBelow != NULL);

    /* Number of split keys to send to each processor */
    counts = malloc(sizeof(*counts) * nThreads);
    assert(counts != NULL);
    for(i=0; i<nDomains; ++i) counts[i] = nRungs*2;

    /* Domain indexes for each rung */
    piIndex = malloc(sizeof(*piIndex) * nRungs*2);
    assert(piIndex != NULL);

    /* Number of active particles on each rung */
    nRungTotal = malloc(sizeof(*nRungTotal) * nRungs*2);
    assert(nRungTotal != NULL);

    /* The targets for each rung  */
    nTarget = malloc(sizeof(*nTarget) * nRungs*2);
    assert(nTarget!=NULL);

    /* What we actually found */
    nBelowGlobal = malloc(sizeof(*nBelowGlobal) * nRungs*2);
    assert(nBelowGlobal!=NULL);

    myBound = malloc(sizeof(*myBound) * nRungs * 2);
    assert(myBound!=NULL);

    /* The split keys - active and inactive for each rung */
    splitKeys = malloc(sizeof(*splitKeys) * nDomains * nRungs*2);
    assert(splitKeys!=NULL);

    /* The results - active and inactive below the split */
    nBelowLocal = malloc(sizeof(*nBelowLocal) * nDomains * nRungs*2);
    assert(nBelowLocal!=NULL);

    /* Create a Peano-Hilbert key for each particle, then sort them into key order */
    /* Let's take this opportunity to calculate the first active rung */
    iRungMin = MAX_RUNG;
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	if (p->uRung < iRungMin) iRungMin = p->uRung;
	if (p->uRung > iRungMax) iRungMax = p->uRung;
	x = pkdPos(pkd,p,0) + 1.5;
	if (x < 1.0) x = 1.0;
	else if (x >= 2.0) x = 2.0;
	y = pkdPos(pkd,p,1) + 1.5;
	if (y < 1.0) y = 1.0;
	else if (y >= 2.0) y = 2.0;
	z = pkdPos(pkd,p,2) + 1.5;
	if (z < 1.0) z = 1.0;
	else if (z >= 2.0) z = 2.0;

#if PEANO_HILBERT_KEY_MAX > 0x3ffffffffffll
	pl[i].lKey = hilbert3d(x,y,z);
#else
	pl[i].lKey = hilbert2d(x,y);
#endif
	assert(pl[i].lKey >= 0 && pl[i].lKey <= PEANO_HILBERT_KEY_MAX);
	pl[i].i = i;
	}
    mdlAllreduce( pkd->mdl, &iRungMin, &pkd->iFirstDomainRung, 1, MDL_INT, MDL_MIN);
    mdlAllreduce( pkd->mdl, &iRungMax, &pkd->iLastDomainRung, 1, MDL_INT, MDL_MAX);
    qsort(pl,nLocal,sizeof(PLITEDD),cmpPeanoHilbert);
    pl[nLocal].lKey = PEANO_HILBERT_KEY_MAX;
    pl[nLocal].i = -1;

    /*
    ** Calculate the number of active particles to the left of the given particle, for every particle
    ** and for every rung. The number of inactive is this number subtracted from the particle index.
    ** The numbers to the right can also easily be calculated.
    */
    p = pkdParticle(pkd,pl[0].i);
    for( iRung=0; iRung<nRungs; iRung++) nActiveBelow[AI(0,iRung)] = nActiveBelow[II(0,iRung)] = 0;
    for (i=1;i<=nLocal;++i) {
	p = pkdParticle(pkd,pl[i-1].i);
	iActiveRung = p->uRung - pkd->iFirstDomainRung + 1;
	if (iActiveRung>nRungs) iActiveRung = nRungs;
	for( iRung=0; iRung<nRungs; iRung++) {
	    nActiveBelow[AI(i,iRung)] = nActiveBelow[AI(i-1,iRung)];
	    nActiveBelow[II(i,iRung)] = nActiveBelow[II(i-1,iRung)];
	    }
	if ( iMethod==3 || iMethod==4) {
	    ++nActiveBelow[AI(i,iActiveRung-1)];
	    }
	else {
	    for( iRung=0; iRung<iActiveRung; iRung++) {
		++nActiveBelow[AI(i,iRung)];
		}
	    }
	if (iMethod==2 || iMethod==4 ) {
	    if (iActiveRung<nRungs)
		++nActiveBelow[II(i,iActiveRung)];
	    }
	else {
	    for( iRung=iActiveRung; iRung<nRungs; iRung++ ) {
		++nActiveBelow[II(i,iRung)];
		}
	    }
	}

    /*
    ** We need to calculate a target for our processor which is simply our share of the active
    ** and inactive particles. For that, we need to know the total number of active on each rung.
    ** This is the sum of the active counts for the last particle.
    */
    for( iRung=0; iRung<nRungs; iRung++) {
	nTarget[AI(0,iRung)] = nActiveBelow[AI(nLocal,iRung)];
	nTarget[II(0,iRung)] = nActiveBelow[II(nLocal,iRung)];
	}
    mdlAllreduce( pkd->mdl, nTarget, nRungTotal, nRungs*2, MDL_LONG_LONG, MDL_SUM);

    /* This is the number of particles we want to the left of our split
    ** We need P-1 roots, so the root processor doesn't need to root find,
    ** so we set the target to zero which indicates we don't want to root find.
    ** Also, if there are no particles in our domain, we skip root finding as well.
    */
    for( iRung=0; iRung<nRungs*2; iRung++) {
	nTarget[iRung] = (int64_t)((nID+0.0) / nDomains  * nRungTotal[iRung]);
	}

    /*
    ** It could be at any key, so we set the range to [0,MAX] and evaluation our
    ** cost function which is known at the boundaries. Set for each rung, and for
    ** each rung both the active and inactive counts.
    */
    for( iRung=0; iRung<nRungs*2; iRung++) {
	r[iRung]  = 0ll;
	s[iRung]  = PEANO_HILBERT_KEY_MAX;
	fr[iRung] = -nTarget[iRung];
	fs[iRung] = nRungTotal[iRung] - nTarget[iRung];
	/*
	** Setting t to zero indicates we have already converged and we should
	** no longer root find.
	*/
	myBound[iRung] = t[iRung] = nTarget[iRung] > 0;
	}

    /* Keep root finding until converged */
    iIter = nIter = 0;
    do {
	iIter++;

	/* If we haven't converged, then we calculate the next guess. */
	for( iRung=0; iRung<nRungs*2; iRung++) {
	    if ( t[iRung]!=0 ) {
#ifdef USE_ILLINOIS_ROOT_FINDER
		int64_t lAdjust = (s[iRung] - r[iRung]) * (fr[iRung]/(fs[iRung] - fr[iRung]));
		myBound[iRung] = t[iRung] = r[iRung] - lAdjust;
#else
		/* Simple bisection */
		myBound[iRung] = t[iRung] = r[iRung]/2 + s[iRung]/2;
		if ( (r[iRung]&1) && (s[iRung]&1) ) myBound[iRung] = ++t[iRung];
#endif
		assert(t[iRung] >= 0 && t[iRung] <= PEANO_HILBERT_KEY_MAX);
		}
	    }

	/* Now we must distribute the guess "t" to all processors */
	mdlAllGather(pkd->mdl,t,nRungs*2,MDL_LONG_LONG, splitKeys, nRungs*2, MDL_LONG_LONG);

	/* Check if we are done */
	for(i=0; i<nDomains*nRungs*2; i++) if (splitKeys[i]!=0) break;
	if (i==nDomains*nRungs*2) break;

	nBelowLocal[0] = 0;
	for(iKey=1; iKey<nDomains; iKey++) { /* O(P-1) */
	    for(iRung=0; iRung<nRungs*2; ++iRung) {
		lKey = splitKeys[iKey*nRungs*2 + iRung];
		if (lKey == 0)
		    nBelowLocal[iKey*nRungs*2 + iRung] = 0;
		else {
		    L = 0;
		    U = nLocal - 1;
		    while(L<U) { /* O(log(N/P/P)) */
			M = (U+L) / 2;
			if ( pl[M].lKey > lKey ) U = M;
			else L = M+1;
			}
		    if (pl[L].lKey <= lKey) L++;
		    nBelowLocal[iKey*nRungs*2 + iRung] = nActiveBelow[L*nRungs*2 + iRung];
		    assert(pl[L-1].lKey <= lKey);
		    assert(pl[L].lKey > lKey);
		    }
		}
	    }

	/* Sum the results and scatter back */
	mdlReduceScatter(pkd->mdl, nBelowLocal, nBelowGlobal, counts, MDL_LONG_LONG, MDL_SUM);

	/* Now do the root finding for each rung */
	for( iRung=0; iRung<nRungs*2; iRung++) {
	    if ( t[iRung]>0 ) {
		ft = nBelowGlobal[iRung] - nTarget[iRung];
		if (fabs(ft) <= 0.0 || labs(r[iRung]-s[iRung]) <= 1 || iIter > 80) { /* Criterian for exit */
		    t[iRung] = 0;
		    nIter = iIter;
		    }
#ifdef USE_ILLINOIS_ROOT_FINDER
		else if (ft*fs[iRung] < 0) {
		    /*
		    ** Unmodified step.
		    */
		    r[iRung] = s[iRung];
		    s[iRung] = t[iRung];
		    fr[iRung] = fs[iRung];
		    fs[iRung] = ft;
		    }
		else {
		    /*
		    ** Modified step to make sure we do not retain the 
		    ** endpoint r indefinitely.
		    */
		    double phis = ft/fs[iRung];
		    double phir = ft/fr[iRung];
		    double gamma = 1.0 - (phis/(1.0-phir));  /* method 3 */
		    if (gamma < 0) gamma = 0.5; /* illinois */
		    fr[iRung] *= gamma;
		    s[iRung] = t[iRung];
		    fs[iRung] = ft;
		    }
#else
		else if (ft<0.0) {
		    r[iRung] = t[iRung];
		    fr[iRung] = ft;
		    }
		else {
		    s[iRung] = t[iRung];
		    fs[iRung] = ft;
		    }
#endif
		}
	    }
	} while(1);

    /* Send the bounds to each processor so we can assign particles to domains */
    mdlAllGather(pkd->mdl,myBound,nRungs*2,MDL_LONG_LONG, splitKeys, nRungs*2, MDL_LONG_LONG);

    for(iKey=1; iKey<nDomains; iKey++) {
	for(iRung=0; iRung<nRungs*2; ++iRung) {
	    if( splitKeys[(iKey-1)*nRungs*2 + iRung] > splitKeys[iKey*nRungs*2 + iRung]) {
		printf( "%d: rung %d - domain %d - %llx <= %llx <= %llx !!<=!! %llx <= %llx\n", nID, iRung, iKey,
		    splitKeys[(iKey-3)*nRungs*2 + iRung], 
		    splitKeys[(iKey-2)*nRungs*2 + iRung], 
		    splitKeys[(iKey-1)*nRungs*2 + iRung],
		    splitKeys[iKey*nRungs*2 + iRung],
		    splitKeys[(iKey+1)*nRungs*2 + iRung] );
		}
	    assert( splitKeys[(iKey-1)*nRungs*2 + iRung] <= splitKeys[iKey*nRungs*2 + iRung]);
	    splitKeys[(iKey-1)*nRungs*2 + iRung] = splitKeys[iKey*nRungs*2 + iRung];
	    }
	}
    for(iRung=0; iRung<nRungs*2; ++iRung) {
	splitKeys[(nDomains-1)*nRungs*2 + iRung] = PEANO_HILBERT_KEY_MAX;
	}


    /* Check that things are valid */
    for(iKey=1; iKey<nDomains; iKey++) { /* O(P-1) */
	for(iRung=0; iRung<nRungs*2; ++iRung) {
	    assert(splitKeys[(iKey-1)*nRungs*2 + iRung] <= splitKeys[iKey*nRungs*2 + iRung]);
	    }
	}

    /*
    ** The splits are now in order for each rung so we can optimize the update.
    */
    for(iRung=0; iRung<nRungs*2; ++iRung) {
	piIndex[iRung] = 0;
	}
    for(i=0; i<nLocal; i++) {
	uint16_t *pRungDest;
	p = pkdParticle(pkd,pl[i].i);
	pRungDest = pkdRungDest(pkd,p);
	if (iMethod==0) {
	    for(iRung=0; iRung<nRungs; ++iRung) {
		while(pl[i].lKey > splitKeys[piIndex[iRung*2]*nRungs*2 + iRung*2])
		    piIndex[iRung*2]++;
		pRungDest[iRung] = piIndex[iRung*2];
		assert(pRungDest[iRung]<nDomains);
		}
	    }
	else {
	    iActiveRung = p->uRung-pkd->iFirstDomainRung >= nRungs ? nRungs : p->uRung-pkd->iFirstDomainRung+1;
	    if ( iMethod==3 || iMethod==4) {
		while(pl[i].lKey > splitKeys[piIndex[(iActiveRung-1)*2]*nRungs*2 + (iActiveRung-1)*2])
		    piIndex[(iActiveRung-1)*2]++;
		for(iRung=0; iRung<iActiveRung; ++iRung) {
		    pRungDest[iRung] = piIndex[(iActiveRung-1)*2];
		    }
		}
	    else {
		for(iRung=0; iRung<iActiveRung; ++iRung) {
		    while(pl[i].lKey > splitKeys[piIndex[iRung*2]*nRungs*2 + iRung*2])
			piIndex[iRung*2]++;
		    pRungDest[iRung] = piIndex[iRung*2];
		    assert(pRungDest[iRung]<nDomains);
		    }
		}
	    if ( iMethod==2 || iMethod==4) {
		if (iActiveRung<nRungs) {
		    while(pl[i].lKey > splitKeys[piIndex[iActiveRung*2+1]*nRungs*2 + iActiveRung*2+1])
			piIndex[iActiveRung*2+1]++;
		    for(iRung=iActiveRung; iRung<nRungs; ++iRung) {
			pRungDest[iRung] = piIndex[iActiveRung*2+1];
			}
		    }
		}
	    else {
		for(iRung=iActiveRung; iRung<nRungs; ++iRung) {
		    while(pl[i].lKey > splitKeys[piIndex[iRung*2+1]*nRungs*2 + iRung*2+1])
			piIndex[iRung*2+1]++;
		    pRungDest[iRung] = piIndex[iRung*2+1];
		    assert(pRungDest[iRung]<nDomains);
		    }
		}
	    }
	}

    free(r); free(s); free(t);
    free(fr); free(fs);

    free(myBound);
    free(nBelowGlobal);
    free(nTarget);
    free(nActiveBelow);
    free(counts);
    free(piIndex);
    free(nBelowLocal);
    free(splitKeys);
    free(nRungTotal);
    }

/*
** This does an ORB decomposition without moving any particles
*/
typedef struct {
    double r[3];
    int32_t i;
    int     uRung;
    } PLITEORB;


static inline void swapPLITEORB(PLITEORB *pl,int lb, int ub) {
    PLITEORB tmp;
    tmp = pl[lb];
    pl[lb] = pl[ub];
    pl[ub] = tmp;
    }

void pkdOrbSplit(PKD pkd, int iDomain) {
    int nNodes = mdlThreads(pkd->mdl);
    int i,j;

    /* Gather which domains to split */
    mdlAllGather(pkd->mdl,&iDomain,1,MDL_INT, pkd->counts,1,MDL_INT);

    for(i=nNodes-1;i>=0;--i) {
	if (pkd->counts[i]>=0) {
	    j = pkd->counts[i];
	    pkd->iFirstActive[j+1] = pkd->iFirstActive[i+1];
	    pkd->iFirstActive[j]   = pkd->iSplitActive[i];
	    pkd->iFirstActive[i+1] = pkd->iSplitActive[i];
	    pkd->iFirstInActive[j+1] = pkd->iFirstInActive[i+1];
	    pkd->iFirstInActive[j]   = pkd->iSplitInActive[i];
	    pkd->iFirstInActive[i+1] = pkd->iSplitInActive[i];
	    }
	}
    }

void pkdOrbBegin(PKD pkd, int nRungs) {
    int nNodes = mdlThreads(pkd->mdl);
    PLITEORB *pl = (PLITEORB *)pkd->pLite;
    int i, j;
    PARTICLE *p;
    int iRungMin, iRungMax;

    pkd->iFirstActive  = malloc(sizeof(*pkd->iFirstActive)   * (nNodes+1));
    pkd->iFirstInActive= malloc(sizeof(*pkd->iFirstInActive) * (nNodes+1));
    pkd->iSplitActive  = malloc(sizeof(*pkd->iSplitActive)   * (nNodes+1));
    pkd->iSplitInActive= malloc(sizeof(*pkd->iSplitInActive) * (nNodes+1));
    pkd->counts        = malloc(sizeof(*pkd->counts)         * nNodes);
    pkd->rdisps        = malloc(sizeof(*pkd->rdisps)         * nNodes);
    pkd->cSplitDims    = malloc(sizeof(*pkd->cSplitDims)     * nNodes);
    pkd->dSplits       = malloc(sizeof(*pkd->dSplits)        * nNodes);
    pkd->pDomainCountsLocal = malloc(sizeof(*pkd->pDomainCountsLocal) * nNodes);

    /* These are the "particles" we partition */
    iRungMin = MAX_RUNG;
    iRungMax = 0;
    for (i=0;i<pkd->nLocal;++i) {
        p = pkdParticle(pkd,i);
	if (p->uRung < iRungMin) iRungMin = p->uRung;
	if (p->uRung > iRungMax) iRungMax = p->uRung;
        for(j=0;j<3;++j) pl[i].r[j] = pkdPos(pkd,p,j);
        pl[i].i = i;
        pl[i].uRung = p->uRung;
        }
    mdlAllreduce( pkd->mdl, &iRungMin, &pkd->iFirstDomainRung, 1, MDL_INT, MDL_MIN);
    mdlAllreduce( pkd->mdl, &iRungMax, &pkd->iLastDomainRung,  1, MDL_INT, MDL_MAX);
    }

int pkdOrbSelectRung(PKD pkd,int uRung) {
    PLITEORB *pl = (PLITEORB *)pkd->pLite;
    int j, lb, ub;

    pkd->iDomainRung = uRung;
    uRung += pkd->iFirstDomainRung;

    /* We start with all particles */
    pkd->iFirstActive[0] = 0;
    pkd->iFirstActive[1] = pkd->nLocal;
    pkd->iFirstInActive[0] = pkd->iFirstInActive[1] = pkd->nLocal;
    pkd->iSplitActive[1] = pkd->iSplitInActive[1] = pkd->nLocal;

    /* Separate them into active and inactive */
    for(j=0;j<3;++j) pl[pkd->nLocal].r[j] = HUGE_VAL; /* sentinal node */
    if ( uRung > 0 ) {
        lb = pkd->iFirstActive[0];
        ub = pkd->iFirstActive[1] - 1;
        PARTITION(lb<ub,lb<=ub,++lb,--ub,swapPLITEORB(pl,lb,ub),
            pl[lb].uRung >= uRung,pl[ub].uRung < uRung);
        pkd->iFirstActive[1] = pkd->iFirstInActive[0] = lb;
        }
    else lb = pkd->nLocal;
    return lb;
    }

void pkdOrbUpdateRung(PKD pkd) {
    int nNodes = mdlThreads(pkd->mdl);
    int iRung = pkd->iDomainRung;
    PLITEORB *pl = (PLITEORB *)pkd->pLite;
    PARTICLE *p;
    int i, iDomain;

    for(iDomain=0; iDomain<nNodes; ++iDomain) {
	for( i=pkd->iFirstActive[iDomain]; i<pkd->iFirstActive[iDomain+1]; ++i) {
	    p = pkdParticle(pkd,pl[i].i);
	    pkdRungDest(pkd,p)[iRung] = iDomain;
	    }
	}
    for(iDomain=0; iDomain<nNodes; ++iDomain) {
	for( i=pkd->iFirstInActive[iDomain]; i<pkd->iFirstInActive[iDomain+1]; ++i) {
	    p = pkdParticle(pkd,pl[i].i);
	    pkdRungDest(pkd,p)[iRung] = iDomain;
	    }
	}
    }

void pkdOrbFinish(PKD pkd) {
    free(pkd->iFirstActive);
    free(pkd->iFirstInActive);
    free(pkd->iSplitActive);
    free(pkd->iSplitInActive);
    free(pkd->counts);
    free(pkd->rdisps);
    free(pkd->cSplitDims);
    free(pkd->dSplits);
    }


int pkdOrbRootFind(
    PKD pkd, double dFraction, uint64_t nLowerMax, uint64_t nUpperMax,
    double dReserveFraction, BND *bnd, double *dSplitOut, int *iDim) {
    /**/
    static int nMaxIter = 60;
    /**/

    int nNodes = mdlThreads(pkd->mdl);
    PLITEORB *pl = (PLITEORB *)pkd->pLite;

    uint8_t cSplitDim;
    double dSplit, dSplitMin, dSplitMax;
    double dFracActive;

    int i, j, lb, ub, d;
    int iIter, nUnbalancedTarget;
    int nDomainsActive, iDomain, iActive, iWork;

    ORBCOUNT DomainCountGlobal;

    if (dSplitOut) *dSplitOut = HUGE_VAL;

    /* Now figure out which dimension to split along -- the longest */
    if (dFraction>0.0) {
        cSplitDim = 0;
	for(j=1;j<3;++j) if (bnd->fMax[j]>bnd->fMax[cSplitDim]) cSplitDim=j;
        dSplit    = bnd->fCenter[cSplitDim];
        dSplitMin = dSplit - bnd->fMax[cSplitDim];
        dSplitMax = dSplit + bnd->fMax[cSplitDim];
        }
    else cSplitDim = 255;
    if (iDim) *iDim = cSplitDim;

    /* Split dimension or 255 if we are not currently directing a split */
    mdlAllGather(pkd->mdl,&cSplitDim,1,MDL_BYTE, pkd->cSplitDims,1,MDL_BYTE);

    pkd->rdisps[0] = 0;
    pkd->counts[0] = 1;
    nDomainsActive = 0;
    for(i=0;i<nNodes;++i) {
	if (pkd->cSplitDims[i]<255) ++nDomainsActive;
	}
    if (nDomainsActive==0) return 0;

    nUnbalancedTarget = 0;
    iIter = 0;
    do {
	iWork = 0;

	pkd->rdisps[0] = 0;
	pkd->counts[0] = 1;
	for(i=0;i<nNodes;++i) {
	    if (pkd->cSplitDims[i]<255) pkd->counts[i] = 1;
	    else pkd->counts[i] = 0;
	    pkd->rdisps[i] = pkd->rdisps[i-1] + pkd->counts[i-1];
	    }

        mdlAllGatherv(pkd->mdl,&dSplit,dFraction>0.0,MDL_DOUBLE,
	    pkd->dSplits, pkd->counts, pkd->rdisps, MDL_DOUBLE);

        /* Now partition about this split and send back the counts */
	iActive = -1;
        for(iDomain=0; iDomain<nNodes; ++iDomain) {
            pkd->counts[iDomain] = 4;
	    if ( (d=pkd->cSplitDims[iDomain]) == 255 ) {
		pkd->pDomainCountsLocal[iDomain].nActiveBelow = 0;
		pkd->pDomainCountsLocal[iDomain].nActiveAbove = 0;
		pkd->pDomainCountsLocal[iDomain].nTotalBelow = 0;
		pkd->pDomainCountsLocal[iDomain].nTotalAbove = 0;
		continue;
		}
	    else ++iActive;

            if ( pkd->dSplits[iActive] == HUGE_VAL ) {
                continue;
                }

            iWork++;

            /* Only active particles here */
            lb = pkd->iFirstActive[iDomain];
            ub = pkd->iFirstActive[iDomain+1] - 1;
	    if (lb <= ub) {
		PARTITION(lb<ub,lb<=ub,++lb,--ub,swapPLITEORB(pl,lb,ub),
		    pl[lb].r[d] < pkd->dSplits[iActive],pl[ub].r[d] >= pkd->dSplits[iActive]);
		}
	    else lb = pkd->iFirstActive[iDomain+1];
            pkd->iSplitActive[iDomain] = lb;
            pkd->pDomainCountsLocal[iDomain].nActiveBelow = lb - pkd->iFirstActive[iDomain];
            pkd->pDomainCountsLocal[iDomain].nActiveAbove = pkd->iFirstActive[iDomain+1] - lb;

	    for(i=pkd->iFirstActive[iDomain]; i<pkd->iSplitActive[iDomain]; ++i)
		assert(pl[i].r[d]<pkd->dSplits[iActive]);
	    for(i=pkd->iSplitActive[iDomain]; i<pkd->iFirstActive[iDomain+1]; ++i)
		assert(pl[i].r[d]>=pkd->dSplits[iActive]);

            /* Now partition inactive and accumulate totals */
            lb = pkd->iFirstInActive[iDomain];
            ub = pkd->iFirstInActive[iDomain+1] - 1;
	    if (lb <= ub) {
		PARTITION(lb<ub,lb<=ub,++lb,--ub,swapPLITEORB(pl,lb,ub),
		    pl[lb].r[d] < pkd->dSplits[iActive],pl[ub].r[d] >= pkd->dSplits[iActive]);
		}
	    else lb = pkd->iFirstInActive[iDomain+1];
	    pkd->iSplitInActive[iDomain] = lb;
	    pkd->pDomainCountsLocal[iDomain].nTotalBelow = lb - pkd->iFirstInActive[iDomain]
		+ pkd->pDomainCountsLocal[iDomain].nActiveBelow;
	    pkd->pDomainCountsLocal[iDomain].nTotalAbove = pkd->iFirstInActive[iDomain+1] - lb
		+ pkd->pDomainCountsLocal[iDomain].nActiveAbove;
            }

        mdlReduceScatter(pkd->mdl, pkd->pDomainCountsLocal, &DomainCountGlobal,
	    pkd->counts, MDL_LONG_LONG, MDL_SUM);

        if (dFraction>0.0 && dSplit < HUGE_VAL) {
	    uint64_t nTotal = DomainCountGlobal.nTotalBelow + DomainCountGlobal.nTotalAbove;
	    uint64_t nActive= DomainCountGlobal.nActiveBelow + DomainCountGlobal.nActiveAbove;
	    uint64_t nUnbalanced = llabs(DomainCountGlobal.nActiveBelow - dFraction*nActive)
		+ llabs(DomainCountGlobal.nActiveAbove - (1.0-dFraction)*nActive);
	    uint64_t nLowerLimit, nUpperLimit;

	    /* Here we calculate the real limits, based on how much we should reserve */
	    if ( nLowerMax==0 || nUpperMax==0 ) {
		nLowerLimit = nLowerMax;
		nUpperLimit = nUpperMax;
		}
	    else {
		uint64_t nMax = nLowerMax + nUpperMax;
		uint64_t nReserve = (nMax - nTotal) * dReserveFraction;
		assert(nMax>=nTotal);
		nLowerLimit = nLowerMax - nReserve*(1.0*nLowerMax/nMax);
		nUpperLimit = nUpperMax - nReserve*(1.0*nUpperMax/nMax);
		}

	    if (nActive<2) {
		if (nTotal<2) dFracActive = dFraction;
		else dFracActive=1.0*DomainCountGlobal.nTotalBelow/nTotal;
		}
            else dFracActive=1.0*DomainCountGlobal.nActiveBelow/nActive;
	    assert(DomainCountGlobal.nTotalBelow+DomainCountGlobal.nTotalAbove <= nLowerLimit+nUpperLimit);
	    if (DomainCountGlobal.nTotalBelow > nLowerLimit) {
		dSplitMax = dSplit;
		dSplit = (dSplitMin+dSplitMax) * 0.5;
		}
	    else if (DomainCountGlobal.nTotalAbove > nUpperLimit) {
		dSplitMin = dSplit;
		dSplit = (dSplitMin+dSplitMax) * 0.5;
		}
	    else if (DomainCountGlobal.nTotalBelow == nLowerLimit && dFracActive <= dFraction) {
		if (dSplitOut) *dSplitOut = dSplit;
		dSplit = HUGE_VAL;
		}
	    else if (DomainCountGlobal.nTotalAbove == nUpperLimit && dFracActive >= dFraction ) {
		if (dSplitOut) *dSplitOut = dSplit;
		dSplit = HUGE_VAL;
		}
	    else if ( nUnbalanced <= nUnbalancedTarget ) {
		if (dSplitOut) *dSplitOut = dSplit;
		dSplit = HUGE_VAL;
		}
            else {
		++nUnbalancedTarget;
                if ( dFracActive > dFraction ) dSplitMax = dSplit;
                else dSplitMin = dSplit;
		if (dSplitOut) *dSplitOut = dSplit;
                dSplit = (dSplitMin+dSplitMax) * 0.5;
                }
            }
        ++iIter;
        } while(iIter < nMaxIter && iWork);
    assert(iIter<nMaxIter);
    return nDomainsActive;
    }



void pkdRungOrder(PKD pkd, int iRung, total_t *nMoved) {
    int nDomains = mdlThreads(pkd->mdl);
    int nLocal   = pkd->nLocal;
    int nID      = mdlSelf(pkd->mdl);
    int *scounts, *rcounts, *sdispls, *rdispls, *ioffset;
    int nSelf, iSelf, iTarget;
    int i;
    PARTICLE *p, *p2;
    uint16_t *pRungDest;

    if (iRung>pkd->iLastDomainRung) iRung = pkd->iLastDomainRung;
    iRung -= pkd->iFirstDomainRung;
    if (iRung<0) iRung = 0;
    else if (iRung >= pkd->nMaxDomainRungs) iRung = pkd->nMaxDomainRungs - 1;

    scounts = malloc(sizeof(*scounts) * nDomains);
    assert(scounts != NULL);
    rcounts = malloc(sizeof(*rcounts) * nDomains);
    assert(rcounts != NULL);
    sdispls = malloc(sizeof(*sdispls) * nDomains);
    assert(sdispls != NULL);
    rdispls = malloc(sizeof(*rdispls) * nDomains);
    assert(rdispls != NULL);
    ioffset = malloc(sizeof(*ioffset) * nDomains);
    assert(ioffset != NULL);

    /*
    ** Count the number of particles destined for each processor and tell
    ** the other processor this information.
    */
    for(i=0; i<nDomains; ++i) scounts[i] = 0;
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	pRungDest = pkdRungDest(pkd,p);
	scounts[pRungDest[iRung]]++;
	}
    nSelf = scounts[nID];
    scounts[nID] = 0;
    mdlAlltoall( pkd->mdl, scounts, 1, MDL_INT, rcounts, 1, MDL_INT );

    ioffset[0] = sdispls[0] = rdispls[0] = 0;
    for(i=1; i<nDomains; i++) {
	ioffset[i] = sdispls[i] = sdispls[i-1] + scounts[i-1];
	rdispls[i] = rdispls[i-1] + rcounts[i-1];
	}
    assert(sdispls[nDomains-1] + scounts[nDomains-1] == nLocal-nSelf);
    if(rdispls[nDomains-1] + rcounts[nDomains-1] + nSelf > pkd->nStore) {
	printf("%d <= %d\n",
	    rdispls[nDomains-1] + rcounts[nDomains-1] + nSelf, pkd->nStore);
	}
    assert(rdispls[nDomains-1] + rcounts[nDomains-1] + nSelf <= pkd->nStore);


    /* Put the particles into domain order */
    iSelf = 0;
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	pRungDest = pkdRungDest(pkd,p);
	iTarget = pRungDest[iRung];
	if (iTarget!=nID) {
	    p2 = pkdParticle2(pkd,ioffset[iTarget]++);
	    pkdCopyParticle(pkd,p2,p);
	    }
	else if (iSelf==i) iSelf++;
	else  {
	    p2 = pkdParticle(pkd,iSelf++);
	    pkdCopyParticle(pkd,p2,p);
	    }
	}
    assert(iSelf == nSelf);
    if (*nMoved) *nMoved = nLocal - nSelf;

    /* Send the particles to their correct processors and update our local count */
    mdlAlltoallv(pkd->mdl,pkdParticle2(pkd,0), scounts, sdispls, pkd->typeParticle,
	pkdParticle(pkd,nSelf), rcounts, rdispls, pkd->typeParticle);
    nLocal = pkd->nLocal = rdispls[nDomains-1] + rcounts[nDomains-1] + nSelf;

    /* Check that this worked */
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	pRungDest = pkdRungDest(pkd,p);
	iTarget = pRungDest[iRung];
	assert(iTarget==nID);
	}

    /* The bounds have likely changed */
    pkdCalcBound(pkd,&pkd->bnd);

    free(scounts);
    free(rcounts);
    free(sdispls);
    free(rdispls);
    free(ioffset);
    }
#endif

/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeight(PKD pkd,int d,FLOAT fSplit,int iSplitSide,int iFrom,int iTo,
	      int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh) {
    int i,iPart;
    FLOAT fLower,fUpper;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
	iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
	*pnLow = pkdLocal(pkd)-iPart;
	*pnHigh = iPart;
	}
    else {
	iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    /*
    ** Calculate the lower weight and upper weight BETWEEN the particles
    ** iFrom to iTo!
    */
    fLower = 0.0;
    for (i=iFrom;i<iPart;++i) {
	fLower += 1.0;
	}
    fUpper = 0.0;
    for (i=iPart;i<=iTo;++i) {
	fUpper += 1.0;
	}
    if (iSplitSide) {
	*pfLow = fUpper;
	*pfHigh = fLower;
	}
    else {
	*pfLow = fLower;
	*pfHigh = fUpper;
	}
    return(iPart);
    }


void pkdCountVA(PKD pkd,int d,FLOAT fSplit,int *pnLow,int *pnHigh) {
    PARTICLE *p;
    int i;

    *pnLow = 0;
    *pnHigh = 0;
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsVeryActive(pkd,p)) {
	    if (pkdPos(pkd,p,d) < fSplit) *pnLow += 1;
	    else *pnHigh += 1;
	    }
	}
    }

/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeightWrap(PKD pkd,int d,FLOAT fSplit,FLOAT fSplit2,int iSplitSide,int iVASplitSide,
		  int iFrom,int iTo,int *pnLow,int *pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (!iSplitSide) {
	iPart = pkdLowerPartWrap(pkd,d,fSplit,fSplit2,iVASplitSide,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    else {
	iPart = pkdUpperPartWrap(pkd,d,fSplit,fSplit2,iVASplitSide,iFrom,iTo);
	*pnHigh = iPart;
	*pnLow = pkdLocal(pkd)-iPart;
	}
    return(iPart);
    }


int pkdOrdWeight(PKD pkd,uint64_t iOrdSplit,int iSplitSide,int iFrom,int iTo,
		 int *pnLow,int *pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
	iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
	*pnLow = pkdLocal(pkd)-iPart;
	*pnHigh = iPart;
	}
    else {
	iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    return(iPart);
    }


int pkdLowerPart(PKD pkd,int d,FLOAT fSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	pkdSwapParticle(pkd,pi,pj),
	pkdPos(pkd,pi,d) >= fSplit,pkdPos(pkd,pj,d) < fSplit);
    return(i);
    }


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	pkdSwapParticle(pkd,pi,pj),
	pkdPos(pkd,pi,d) < fSplit,pkdPos(pkd,pj,d) >= fSplit);
    return(i);
    }


int pkdLowerPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE *pi = pkdParticle(pkd,i);
    PARTICLE *pj = pkdParticle(pkd,j);

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 || pkdPos(pkd,pi,d) >= fSplit1) &&
		       !pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) >= fSplit2 && pkdPos(pkd,pj,d) < fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 || pkdPos(pkd,pi,d) >= fSplit1) ||
		pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) >= fSplit2 && pkdPos(pkd,pj,d) < fSplit1) &&
		!pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 || pkdPos(pkd,pi,d) >= fSplit1),
	    (pkdPos(pkd,pj,d) >= fSplit2 && pkdPos(pkd,pj,d) < fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 && pkdPos(pkd,pi,d) >= fSplit1) &&
		!pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) >= fSplit2 || pkdPos(pkd,pj,d) < fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 && pkdPos(pkd,pi,d) >= fSplit1) ||
		pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) >= fSplit2 || pkdPos(pkd,pj,d) < fSplit1) &&
		!pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) < fSplit2 && pkdPos(pkd,pi,d) >= fSplit1),
	    (pkdPos(pkd,pj,d) >= fSplit2 || pkdPos(pkd,pj,d) < fSplit1));
	    }
	}
    return(i);
    }


int pkdUpperPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE *pi = pkdParticle(pkd,i);
    PARTICLE *pj = pkdParticle(pkd,j);

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 && pkdPos(pkd,pi,d) < fSplit1) ||
		pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) < fSplit2 || pkdPos(pkd,pj,d) >= fSplit1) &&
		!pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 && pkdPos(pkd,pi,d) < fSplit1) &&
		!pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) < fSplit2 || pkdPos(pkd,pj,d) >= fSplit1) ||
		pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 && pkdPos(pkd,pi,d) < fSplit1),
	    (pkdPos(pkd,pj,d) < fSplit2 || pkdPos(pkd,pj,d) >= fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 || pkdPos(pkd,pi,d) < fSplit1) ||
		pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) < fSplit2 && pkdPos(pkd,pj,d) >= fSplit1) &&
		!pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 || pkdPos(pkd,pi,d) < fSplit1) &&
		!pkdIsVeryActive(pkd,pi),
	    (pkdPos(pkd,pj,d) < fSplit2 && pkdPos(pkd,pj,d) >= fSplit1) ||
		pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		pkdSwapParticle(pkd,pi,pj),
	    (pkdPos(pkd,pi,d) >= fSplit2 || pkdPos(pkd,pi,d) < fSplit1),
	    (pkdPos(pkd,pj,d) < fSplit2 && pkdPos(pkd,pj,d) >= fSplit1));
	    }
	}
    return(i);
    }


int pkdLowerOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	       pkdSwapParticle(pkd,pi,pj),
	       pi->iOrder >= nOrdSplit,pj->iOrder < nOrdSplit);
    return(i);
    }


int pkdUpperOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	       pkdSwapParticle(pkd,pi,pj),
	       pi->iOrder < nOrdSplit,pj->iOrder >= nOrdSplit);
    return(i);
    }


int pkdActiveOrder(PKD pkd) {
    int i=0;
    int j=pkdLocal(pkd)-1;
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	       pkdSwapParticle(pkd,pi,pj),
	       pkdIsActive(pkd,pi),!pkdIsActive(pkd,pj));
    return (pkd->nActive = i);
    }


int pkdColRejects(PKD pkd,int nSplit) {
    int iRejects,i;

    mdlassert(pkd->mdl,pkd->nRejects == 0);

    pkd->nRejects = pkdLocal(pkd) - nSplit;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    /*
    ** Move rejects to High memory.
    */
    if (pkdLocal(pkd) != pkdFreeStore(pkd)) {
	for (i=pkd->nRejects-1;i>=0;--i)
	    pkdCopyParticle(pkd,pkdParticle(pkd,iRejects+i),pkdParticle(pkd,nSplit+i));
	}
    pkd->nLocal = nSplit;
    return(pkd->nRejects);
    }


int pkdSwapRejects(PKD pkd,int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;

    if (idSwap != -1) {
	nBuf = (pkdSwapSpace(pkd))*pkdParticleSize(pkd);
	nOutBytes = pkd->nRejects*pkdParticleSize(pkd);
	mdlassert(pkd->mdl,pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
	mdlSwap(pkd->mdl,idSwap,nBuf,pkdParticle(pkd,pkdLocal(pkd)),
		nOutBytes,&nSndBytes,&nRcvBytes);
	pkd->nLocal += nRcvBytes/pkdParticleSize(pkd);
	pkd->nRejects -= nSndBytes/pkdParticleSize(pkd);
	}
    return(pkd->nRejects);
    }

void pkdSwapAll(PKD pkd, int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;

    /*
    ** Move particles to High memory.
    */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkdCopyParticle(pkd,pkdParticle(pkd,iBuf+i),pkdParticle(pkd,i));
    nBuf = pkdFreeStore(pkd)*pkdParticleSize(pkd);
    nOutBytes = pkdLocal(pkd)*pkdParticleSize(pkd);
    mdlSwap(pkd->mdl,idSwap,nBuf,pkdParticleBase(pkd), nOutBytes,
	    &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,nSndBytes/pkdParticleSize(pkd) == pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/pkdParticleSize(pkd);
    }

int pkdSwapSpace(PKD pkd) {
    return(pkdFreeStore(pkd) - pkdLocal(pkd));
    }


int pkdFreeStore(PKD pkd) {
    return(pkd->nStore);
    }

int pkdActive(PKD pkd) {
    return(pkd->nActive);
    }

int pkdInactive(PKD pkd) {
    return(pkd->nLocal - pkd->nActive);
    }

int pkdLocal(PKD pkd) {
    return(pkd->nLocal);
    }

int pkdNodes(PKD pkd) {
    return(pkd->nNodes);
    }

/*
** Returns a pointer to the i'th KDN in the tree.  Used for fetching
** cache element.  Normal code should call pkdTreeNode().
*/
void *pkdTreeNodeGetElement(void *vData,int i,int iDataSize) {
    PKD pkd = vData;
    return pkdTreeNode(pkd,i);
    }

int pkdNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    int i, n;
    for (n=0,i=0;i<pkdLocal(pkd);++i)
	if ( pkdIsSrcActive(pkdParticle(pkd,i),uRungLo,uRungHi) ) n++;
    return n;
    }

int pkdNumDstActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    int i, n;
    for (n=0,i=0;i<pkdLocal(pkd);++i)
	if ( pkdIsDstActive(pkdParticle(pkd,i),uRungLo,uRungHi) ) n++;
    return n;
    }

int pkdColOrdRejects(PKD pkd,uint64_t nOrdSplit,int iSplitSide) {
    int nSplit;
    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    return pkdColRejects(pkd,nSplit);
    }

int cmpParticles(const void *pva,const void *pvb) {
    PARTICLE *pa = (PARTICLE *)pva;
    PARTICLE *pb = (PARTICLE *)pvb;

    return(pa->iOrder - pb->iOrder);
    }


void pkdLocalOrder(PKD pkd,uint64_t iMinOrder, uint64_t iMaxOrder) {
    int i;
    assert(pkd->nLocal == iMaxOrder - iMinOrder + 1);
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p1 = pkdParticle(pkd,i);
	assert(p1->iOrder >= iMinOrder && p1->iOrder <= iMaxOrder);
	while(p1->iOrder - iMinOrder !=  i) {
	    PARTICLE *p2 = pkdParticle(pkd,p1->iOrder-iMinOrder);
	    pkdSwapParticle(pkd,p1,p2);
	    }
	}
    /* Above replaces: qsort(pkdParticleBase(pkd),pkdLocal(pkd),pkdParticleSize(pkd),cmpParticles); */
    }

static void writeParticle(PKD pkd,FIO fio,double dvFac,BND *bnd,PARTICLE *p) {
    STARFIELDS *pStar;
    SPHFIELDS *pSph;
    float *pPot, dummypot;
    double v[3],r[3];
    float fMass, fSoft, fDensity;
    uint64_t iParticleID;
    int j;

    dummypot = 0.0;

    if ( pkd->oPotential) pPot = pkdPot(pkd,p);
    else pPot = &dummypot;
    if (pkd->oVelocity) {
	vel_t *pV = pkdVel(pkd,p);
	v[0] = pV[0] * dvFac;
	v[1] = pV[1] * dvFac;
	v[2] = pV[2] * dvFac;
	}
    else v[0] = v[1] = v[2] = 0.0;
 
    /* Initialize SPH fields if present */
    if (pkd->oSph) pSph = pkdField(p,pkd->oSph);
    else pSph = NULL;
    if (pkd->oStar) pStar = pkdField(p,pkd->oStar);
    else pStar = NULL;
    fMass = pkdMass(pkd,p);
    fSoft = pkdSoft0(pkd,p);
    if (pkd->oParticleID) iParticleID = *pkdParticleID(pkd,p);
    else iParticleID = p->iOrder;
    if (pkd->oDensity) fDensity = pkdDensity(pkd,p);
    else fDensity = 0.0;

    r[0] = pkdPos(pkd,p,0);
    r[1] = pkdPos(pkd,p,1);
    r[2] = pkdPos(pkd,p,2);
    /* Enforce periodic boundaries */
    for (j=0;j<3;++j) {
	if (r[j] < bnd->fCenter[j] - bnd->fMax[j]) r[j] += 2*bnd->fMax[j];
	else if (r[j] >= bnd->fCenter[j] + bnd->fMax[j]) r[j] -= 2*bnd->fMax[j];
	/*
	** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the 
	** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
	** to the specified periodic box in a gross way.
	*/
	mdlassert(pkd->mdl,((r[j] >= bnd->fCenter[j] - bnd->fMax[j])&&
		(r[j] < bnd->fCenter[j] + bnd->fMax[j])));
	
	}

    switch(pkdSpecies(pkd,p)) {
    case FIO_SPECIES_SPH:
	assert(pSph);
	assert(pkd->param.dTuFac>0.0);
	    {
	    double T;
#ifdef COOLING
	    COOLPARTICLE cp;
	    if (pkd->param.bGasCooling) {
		double E = pSph->u;
		CoolTempFromEnergyCode( pkd->Cool, 
		    &cp, &E, &T, p->fDensity, pSph->fMetals );
		}
	    else T = pSph->u/pkd->param.dTuFac;
#else
	    T = pSph->u/pkd->param.dTuFac;
#endif
	    fioWriteSph(fio,iParticleID,r,v,fMass,fSoft,*pPot,
		fDensity,T,pSph->fMetals);
	    }
	break;
    case FIO_SPECIES_DARK:
	fioWriteDark(fio,iParticleID,r,v,fMass,fSoft,*pPot,fDensity);
	break;
    case FIO_SPECIES_STAR:
	assert(pStar && pSph);
	fioWriteStar(fio,iParticleID,r,v,fMass,fSoft,*pPot,fDensity,
	    pSph->fMetals,pStar->fTimer);
	break;
    default:
	fprintf(stderr,"Unsupported particle type: %d\n",pkdSpecies(pkd,p));
	assert(0);
	}

    }

struct packWriteCtx {
    PKD pkd;
    FIO fio;
    BND *bnd;
    double dvFac;
    int iIndex;
    };

static int unpackWrite(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packWriteCtx *ctx = (struct packWriteCtx *)vctx;
    PKD pkd = ctx->pkd;
    PARTICLE *p = (PARTICLE *)vBuff;
    int n = nSize / pkdParticleSize(pkd);
    int i;
    assert( n*pkdParticleSize(pkd) == nSize);
    for(i=0; i<n; ++i) {
	writeParticle(pkd,ctx->fio,ctx->dvFac,ctx->bnd,pkdParticleGet(pkd,p,i));
	}
    return 1;
    }

void pkdWriteFromNode(PKD pkd,int iNode, FIO fio,double dvFac,BND *bnd) {
    struct packWriteCtx ctx;
    ctx.pkd = pkd;
    ctx.fio = fio;
    ctx.bnd = bnd;
    ctx.dvFac = dvFac;
    ctx.iIndex = 0;
#ifdef MPI_VERSION
    mdlRecv(pkd->mdl,iNode,unpackWrite,&ctx);
#endif
    }

static int packWrite(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packWriteCtx *ctx = (struct packWriteCtx *)vctx;
    PKD pkd = ctx->pkd;
    int nLeft = pkd->nLocal - ctx->iIndex;
    int n = nSize / pkdParticleSize(pkd);
    if ( n > nLeft ) n = nLeft;
    nSize = n*pkdParticleSize(pkd);
    memcpy(vBuff,pkdParticle(pkd,ctx->iIndex), nSize );
    ctx->iIndex += n;
    return nSize;
    }

/* Send all particled data to the specified node for writing */
void pkdWriteViaNode(PKD pkd, int iNode) {
    struct packWriteCtx ctx;
    ctx.pkd = pkd;
    ctx.fio = NULL;
    ctx.dvFac = 1.0;
    ctx.iIndex = 0;
#ifdef MPI_VERSION
    mdlSend(pkd->mdl,iNode,packWrite, &ctx);
#endif
    }

uint32_t pkdWriteFIO(PKD pkd,FIO fio,double dvFac,BND *bnd) {
    PARTICLE *p;
    int i;
    uint32_t nCount;
    nCount = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;  /* JW: Ack! */
	writeParticle(pkd,fio,dvFac,bnd,p);
	nCount++;
	}
    return nCount;
    }

void pkdSetSoft(PKD pkd,double dSoft) {
    pkd->fSoftFix = dSoft;
    }

void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul) {
    pkd->fSoftFac = dFac;
    pkd->fSoftMax = bSoftMaxMul ? HUGE : dSoftMax;
    }

void
pkdGravAll(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
    int bKickClose,int bKickOpen,double *dtClose,double *dtOpen,
    double dAccFac,double dTime,int nReps,int bPeriodic,
    int iOrder,int bEwald,int nGroup,int iRoot1, int iRoot2,
    double fEwCut,double fEwhCut,double dThetaMin,
    uint64_t *pnActive,
    double *pdPart,double *pdPartNumAccess,double *pdPartMissRatio,
    double *pdCell,double *pdCellNumAccess,double *pdCellMissRatio,
    double *pdFlop,uint64_t *pnRung) {

    double dActive;
    double dPartSum;
    double dCellSum;
    int i;

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("Gravity");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif

    /*
    ** Clear all the rung counters to be safe.
    */
    for (i=0;i<=IRUNGMAX;++i) pkd->nRung[i] = 0;
     
    pkdClearTimer(pkd,1);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    mdlTimeReset(pkd->mdl);
#endif

    /*
    ** Set up Ewald tables and stuff.
    */
    if (bPeriodic && bEwald) {
	pkdEwaldInit(pkd,nReps,fEwCut,fEwhCut);	/* ignored in Flop count! */
	}
    /*
    ** Start particle caching space (cell cache already active).
    */
    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),
	       pkdLocal(pkd));
    /*
    ** Calculate newtonian gravity, including replicas if any.
    */
    *pdFlop = 0.0;
    dPartSum = 0.0;
    dCellSum = 0.0;
    pkdStartTimer(pkd,1);
    *pnActive = pkdGravWalk(pkd,uRungLo,uRungHi,bKickClose,bKickOpen,dtClose,dtOpen,dAccFac,dTime,nReps,bPeriodic && bEwald,nGroup,iRoot1,iRoot2,0,dThetaMin,pdFlop,&dPartSum,&dCellSum);
    pkdStopTimer(pkd,1);

    dActive = (double)(*pnActive);
    if (*pnActive) {
	*pdPart = dPartSum/dActive;
	*pdCell = dCellSum/dActive;
	}
    else {
	assert(dPartSum == 0 && dCellSum == 0);
	*pdPart = 0;  /* for the statistics we don't count this processor, see pstGravity(). */
	*pdCell = 0;
	}
    /*
    ** Get caching statistics.
    */
    if (*pnActive) {
	*pdCellNumAccess = mdlNumAccess(pkd->mdl,CID_CELL)/dActive;
	*pdPartNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE)/dActive;
	}
    else {
	*pdCellNumAccess = 0;
	*pdPartNumAccess = 0;
	}
    *pdCellMissRatio = 100.0*mdlMissRatio(pkd->mdl,CID_CELL);      /* as a percentage */
    *pdPartMissRatio = 100.0*mdlMissRatio(pkd->mdl,CID_PARTICLE);  /* as a percentage */
    /*
    ** Output flops count in GFlops!
    */
    *pdFlop *= 1e-9;
    /*
    ** Stop particle caching space.
    */
    mdlFinishCache(pkd->mdl,CID_PARTICLE);

    for (i=0;i<=IRUNGMAX;++i) pnRung[i] = pkd->nRung[i];

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
    }

void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double *L,double *F,double *W) {
    /* L is calculated with respect to the origin (0,0,0) */

    PARTICLE *p;
    vel_t *v;
    float *a;
    FLOAT rx,ry,rz,vx,vy,vz;
    float fMass;
    int i,n;

    n = pkdLocal(pkd);
    *T = 0.0;
    *U = 0.0;
    *Eth = 0.0;
    L[0] = L[1] = L[2] = 0.0;
    F[0] = F[1] = F[2] = 0.0;
    *W = 0.0;
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	fMass = pkdMass(pkd,p);
	if (pkd->oPotential) *U += 0.5*fMass*(*(pkdPot(pkd,p)));
	if (pkd->oAcceleration) {
	    a = pkdAccel(pkd,p);
	    *W += fMass*(pkdPos(pkd,p,0)*a[0] + pkdPos(pkd,p,1)*a[1] + pkdPos(pkd,p,2)*a[2]);
	    F[0] += fMass*a[0];
	    F[1] += fMass*a[1];
	    F[2] += fMass*a[2];
	}
	if (pkd->oSph && pkdIsGas(pkd,p)) *Eth += fMass*pkdSph(pkd,p)->u;
	if (pkd->oVelocity) {
	    v = pkdVel(pkd,p);
	    rx = pkdPos(pkd,p,0); ry = pkdPos(pkd,p,1); rz = pkdPos(pkd,p,2);
	    vx = v[0]; vy = v[1]; vz = v[2];
	    L[0] += fMass*(ry*vz - rz*vy);
	    L[1] += fMass*(rz*vx - rx*vz);
	    L[2] += fMass*(rx*vy - ry*vx);
	    }
	}
    *T += pkd->dEnergyT;
    if (!pkd->oPotential) *U = pkd->dEnergyU;
    }


void pkdScaleVel(PKD pkd,double dvFac) {
    PARTICLE *p;
    vel_t *v;
    int i,j,n;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	v = pkdVel(pkd,p);
	for (j=0;j<3;++j) v[j] *= dvFac;
	}
    }

/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd,double dDelta,double dDeltaVPred,double dDeltaUPred,uint8_t uRungLo,uint8_t uRungHi) {
    PARTICLE *p;
    vel_t *v;
    float *a;
    SPHFIELDS *sph;
    int i,j,n;
    double r[3],dMin[3],dMax[3];

    mdlDiag(pkd->mdl, "Into pkdDrift\n");
    assert(pkd->oVelocity);

    for (j=0;j<3;++j) {
	dMin[j] = pkd->bnd.fCenter[j] - pkd->bnd.fMax[j];
	dMax[j] = pkd->bnd.fCenter[j] + pkd->bnd.fMax[j];
	}
    n = pkdLocal(pkd);
    /*
    ** Update particle positions
    */
    if (pkd->param.bDoGas) {
	assert(pkd->oSph);
	assert(pkd->oAcceleration);
	for (i=0;i<n;++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
		v = pkdVel(pkd,p);
		if (pkdIsGas(pkd,p)) {
		    a = pkdAccel(pkd,p);
		    sph = pkdSph(pkd,p);
		    for (j=0;j<3;++j) { /* NB: Pred quantities must be done before std. */
			sph->vPred[j] += a[j]*dDeltaVPred;
			}
		    sph->uPred += sph->uDot*dDeltaUPred;
		    sph->fMetalsPred += sph->fMetalsDot*dDeltaUPred;
		    }
		for (j=0;j<3;++j) {
		    pkdSetPos(pkd,p,j,r[j] = pkdPos(pkd,p,j) + dDelta*v[j]);
		    }
		pkdMinMax(r,dMin,dMax);
		}
	    }
	}
    else {
	for (i=0;i<n;++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
		v = pkdVel(pkd,p);
		for (j=0;j<3;++j) {
		    pkdSetPos(pkd,p,j,r[j] = pkdPos(pkd,p,j) + dDelta*v[j]);
		    }
		pkdMinMax(r,dMin,dMax);
		}
	    }
	}
    for (j=0;j<3;++j) {
	pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
    }


void pkdGravityVeryActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int bEwald,int nGroup,int nReps,
			  double dStep,double dTheta) {
    int nActive;
    double dFlop,dPartSum,dCellSum;

    /*
    ** Calculate newtonian gravity for the very active particles ONLY, including replicas if any.
    */
    dFlop = 0.0;
    dPartSum = 0.0;
    dCellSum = 0.0;
    nActive = pkdGravWalk(pkd,uRungLo,uRungHi,0,0,NULL,NULL,1.0,dTime,nReps,bEwald,nGroup,ROOT,0,VAROOT,dTheta,&dFlop,&dPartSum,&dCellSum);
    }


void pkdStepVeryActiveKDK(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dStep, double dTime, double dDelta,
			  int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double dThetaMin,
			  int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass) {
    uint64_t nRungCount[256];
    double dDriftFac;

    if (iAdjust && (iRung < pkd->param.iMaxRung)) {

	/*
	** The following should be replaced with a single call which sets the rungs of all particles.
	*/
	pkdActiveRung(pkd, iRung, 1);
	if (pkd->param.bAccelStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dVelFac = 1.0/(a*a);
	    double dAccFac = 1.0/(a*a*a);
	    double dhMinOverSoft = 0;
	    pkdAccelStep(pkd,uRungLo,uRungHi,pkd->param.dEta, dVelFac,dAccFac,pkd->param.bDoGravity,
			 pkd->param.bEpsAccStep,dhMinOverSoft);
	    }
	*pnMaxRung = pkdUpdateRung(pkd,iRung,pkd->param.iMaxRung,
				   iRung,pkd->param.iMaxRung, nRungCount);


	if (pkd->param.bVDetails) {
	    printf("%*cAdjust at iRung: %d, nMaxRung:%d nRungCount[%d]=%lld\n",
		   2*iRung+2,' ',iRung,*pnMaxRung,*pnMaxRung,nRungCount[*pnMaxRung]);
	    }

	}
    if (iRung > iRungVeryActive) {	/* skip this if we are
					   entering for the first
					   time: Kick is taken care of
					   in master().
					*/
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickOpen  at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	pkdKickKDKOpen(pkd,dTime,0.5*dDelta,iRung,iRung);
	}
    if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
	pkdStepVeryActiveKDK(pkd,uRungLo,uRungHi,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
			     dThetaMin,pnMaxRung,aSunInact,adSunInact,dSunMass);
	dStep += 1.0/(2 << iRung);
	dTime += 0.5*dDelta;

	pkdActiveRung(pkd,iRung,0);   /* is this needed? */

	pkdStepVeryActiveKDK(pkd,uRungLo,uRungHi,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
			     dThetaMin,pnMaxRung,aSunInact,adSunInact,dSunMass);
	}
    else {
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive Drift at iRung: %d, drifting %d and higher with dDelta: %g\n",
		   2*iRung+2,' ',iRung,iRungVeryActive+1,dDelta);
	    }
	/*
	** We need to account for cosmological drift factor here!
	** Normally this is done at the MASTER level in msrDrift.
	** Note that for kicks we have written new "master-like" functions
	** KickOpen and KickClose which do this same job at PKD level.
	*/
	if (pkd->param.csm->bComove) {
	    dDriftFac = csmComoveDriftFac(pkd->param.csm,dTime,dDelta);
	    }
	else {
	    dDriftFac = dDelta;
	    }
	/*
	** This should drift *all* very actives!
	*/
	pkdDrift(pkd,dDriftFac,0,0,iRungVeryActive+1,MAX_RUNG);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);

	if (iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master().
						*/

	    if (pkd->param.bVDetails) {
		printf("%*cGravityVA: iRung %d Gravity for rungs %d to %d ... ",
		       2*iRung+2,' ',iRung,iKickRung,*pnMaxRung);
		}

	    pkdActiveRung(pkd,iKickRung,1);
//	    pkdVATreeBuild(pkd,pkd->param.nBucket);
	    pkdGravityVeryActive(pkd,uRungLo,uRungHi,dTime,pkd->param.bEwald && pkd->param.bPeriodic,pkd->param.nGroup,
				 pkd->param.nReplicas,dStep,dThetaMin);

	    }
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	}
    if (iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master().
						*/
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickClose at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	pkdKickKDKClose(pkd,dTime,0.5*dDelta,iRung,iRung);
	}
    }


/*
 * Stripped down versions of routines from master.c
 */
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    if (pkd->param.csm->bComove) {
	dDelta = csmComoveKickFac(pkd->param.csm,dTime,dDelta);
    }
    pkdKick(pkd,dTime,dDelta,0,0,0,uRungLo,uRungHi);
    }

void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    if (pkd->param.csm->bComove) {
	dDelta = csmComoveKickFac(pkd->param.csm,dTime,dDelta);
    }
    pkdKick(pkd,dTime,dDelta,0,0,0,uRungLo,uRungHi);
    }


void pkdKick(PKD pkd,double dTime,double dDelta,double dDeltaVPred,double dDeltaU,double dDeltaUPred,uint8_t uRungLo,uint8_t uRungHi) {
    PARTICLE *p;
    vel_t *v;
    float *a;
    SPHFIELDS *sph;
    int i,j,n;

    assert(pkd->oVelocity);
    assert(pkd->oAcceleration);

    pkdClearTimer(pkd,1);
    pkdStartTimer(pkd,1);

    if (pkd->param.bDoGas) {
	assert(pkd->oSph);
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
		a = pkdAccel(pkd,p);
		v = pkdVel(pkd,p);
		if (pkdIsGas(pkd,p)) {
		    sph = pkdSph(pkd,p);
		    for (j=0;j<3;++j) { /* NB: Pred quantities must be done before std. */
			sph->vPred[j] = v[j] + a[j]*dDeltaVPred;
			}
		    sph->uPred = sph->u + sph->uDot*dDeltaUPred;
		    sph->u += sph->uDot*dDeltaU;
		    sph->fMetalsPred = sph->fMetals + sph->fMetalsDot*dDeltaUPred;
		    sph->fMetals += sph->fMetalsDot*dDeltaU;
		    }
		for (j=0;j<3;++j) {
		    v[j] += a[j]*dDelta;
		    }
		}
	    }
	}
    else {
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
		a = pkdAccel(pkd,p);
		v = pkdVel(pkd,p);
		for (j=0;j<3;++j) {
		    v[j] += a[j]*dDelta;
		    }
		}
	    }
	}


    pkdStopTimer(pkd,1);
    mdlDiag(pkd->mdl, "Done pkdkick\n");
    }

/* Kick the tree at iRoot. */
void pkdKickTree(PKD pkd,double dTime,double dDelta,double dDeltaVPred,double dDeltaU,double dDeltaUPred,int iRoot) {
    KDN *c;
    PARTICLE *p;
    vel_t *v;
    float *a;
    int i,j;

    /* Skip to local tree */
    c = pkdTreeNode(pkd,iRoot);
    while(c->bRemote) c = pkdTreeNode(pkd,iRoot = c->iLower);

    /* Now just kick all of the particles in the tree */
    for(i=c->pLower; i<=c->pUpper; ++i) {
	p = pkdParticle(pkd,i);
	a = pkdAccel(pkd,p);
	v = pkdVel(pkd,p);
	for (j=0;j<3;++j) {
	    v[j] += a[j]*dDelta;
	    a[j] = 0.0;
	    }
	}
    }

void pkdInitStep(PKD pkd, struct parameters *p, CSM csm) {
    pkd->param = *p;
    /*
    ** Need to be careful to correctly copy the cosmo
    ** parameters. This is very ugly!
    */
    csmInitialize(&pkd->param.csm);
    *pkd->param.csm = *csm;
    }


void pkdSetRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung) {
    PARTICLE *p;
    int i;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsDstActive(p,uRungLo,uRungHi) ) continue;
	p->uRung = p->uNewRung = uRung;
	}
    }

void pkdZeroNewRung(PKD pkd,uint8_t uRungLo, uint8_t uRungHi, uint8_t uRung) {  /* JW: Ugly -- need to clean up */
    PARTICLE *p;
    int i;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsActive(pkd,p) ) continue;
	p->uNewRung = 0;
	}
    }

void pkdActiveRung(PKD pkd, int iRung, int bGreater) {
    pkd->uMinRungActive = iRung;
    pkd->uMaxRungActive = bGreater ? 255 : iRung;
    }

int pkdCurrRung(PKD pkd,uint8_t uRung) {
    PARTICLE *p;
    int i;
    int iCurrent;

    iCurrent = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (p->uRung == uRung) {
	    iCurrent = 1;
	    break;
	    }
	}
    return iCurrent;
    }

void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
		  double dEta,double dVelFac,double dAccFac,
		  int bDoGravity,int bEpsAcc,double dhMinOverSoft) {
    PARTICLE *p;
    float *a, *pPot;
    vel_t *v;
    int i,uNewRung;
    double vel;
    double acc;
    int j;
    double dT;
    FLOAT fSoft;

    assert(pkd->oVelocity);
    assert(pkd->oAcceleration);

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    v = pkdVel(pkd,p);
	    a = pkdAccel(pkd,p);
	    fSoft = pkdSoft(pkd,p);
	    vel = 0;
	    acc = 0;
	    for (j=0;j<3;j++) {
		vel += v[j]*v[j];
		acc += a[j]*a[j];
		}
	    mdlassert(pkd->mdl,vel >= 0);
	    vel = sqrt(vel)*dVelFac;
	    mdlassert(pkd->mdl,acc >= 0);
	    acc = sqrt(acc)*dAccFac;
	    dT = FLOAT_MAXVAL;
	    if (acc>0) {
		if (bEpsAcc) {
		    dT = dEta*sqrt(fSoft/acc);
		    }
		}
	    uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung);
	    if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
	    }
	}
    }


void pkdSphStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,double dAccFac) {
    PARTICLE *p;
    float *a, uDot;
    int i,j,uNewRung;
    double acc;
    double dtNew;
    int u1,u2,u3;

    assert(pkd->oAcceleration);
    assert(pkd->oSph);

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    if (pkdIsGas(pkd,p)) {
		u1 = p->uNewRung;
		a = pkdAccel(pkd,p);
		acc = 0;
		for (j=0;j<3;j++) {
		    acc += a[j]*a[j];
		    }
		acc = sqrt(acc)*dAccFac;
		dtNew = FLOAT_MAXVAL;
		if (acc>0) dtNew = pkd->param.dEta*sqrt(pkdBall(pkd,p)/acc);
		u2 = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung);
		uDot = *pkd_uDot(pkd,p);
		u3=0;
		if (uDot < 0) {
		    double dtemp = pkd->param.dEtaUDot*(*pkd_u(pkd,p))/fabs(uDot);
		    if (dtemp < dtNew) dtNew = dtemp;
		    u3 = pkdDtToRung(dtemp,pkd->param.dDelta,pkd->param.iMaxRung);
		    }
		uNewRung = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung);
		if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
		if (!(p->iOrder%10000) || (p->uNewRung > 5 && !(p->iOrder%1000))) {
		    SPHFIELDS *sph = pkdSph(pkd,p);
#ifdef COOLING
		    double T, E = sph->u;
		    if (pkd->param.bGasIsothermal) T = E/pkd->param.dTuFac;
		    else {
			COOLPARTICLE cp;
			CoolTempFromEnergyCode( pkd->Cool, &cp, &E, &T, p->fDensity, sph->fMetals );
		    }

#else
		    /*T = E/pkd->param.dTuFac;*/
#endif
		    }
		}
	    }
	}
    }

void pkdStarForm(PKD pkd, double dRateCoeff, double dTMax, double dDenMin,
		 double dDelta, double dTime,
		 double dInitStarMass, double dESNPerStarMass, double dtCoolingShutoff,
		 double dtFeedbackDelay,  double dMassLossPerStarMass,    
		 double dZMassPerStarMass, double dMinGasMass,
		 int bdivv,
		 int *nFormed, /* number of stars formed */
		 double *dMassFormed,	/* mass of stars formed */
		 int *nDeleted) /* gas particles deleted */ {

    PARTICLE *p;
#ifdef COOLING
    COOLPARTICLE cp;
#endif
    SPHFIELDS *sph;
    double T, E, dmstar, dt, prob;
    PARTICLE *starp;
    int i;
    
    assert(pkd->oStar);
    assert(pkd->oSph);
    assert(pkd->oMass);

    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;
    starp = (PARTICLE *) malloc(pkdParticleSize(pkd));
    assert(starp != NULL);

    printf("pkdSF calc dTime %g\n",dTime);
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	
	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    sph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung); /* Actual Rung */
	    pkdStar(pkd,p)->totaltime += dt;
	    if (pkdDensity(pkd,p) < dDenMin || (bdivv && sph->divv >= 0.0)) continue;
	    E = sph->uPred;
#ifdef COOLING
	    if (pkd->param.bGasCooling) 
		CoolTempFromEnergyCode( pkd->Cool, &cp, &E, &T, p->fDensity, sph->fMetals );
	    else T=E/pkd->param.dTuFac;
#else
	    T=E/pkd->param.dTuFac;
#endif
	    if (T > dTMax) continue;
	    
            /* Note: Ramses allows for multiple stars per step -- but we have many particles
	      and he has one cell that may contain many times m_particle */
	    if (pkd->param.bGasCooling) {
		if (fabs(pkdStar(pkd,p)->totaltime-dTime) > 1e-3*dt) {
		    printf("total time error: %lu,  %g %g %g\n",p->iOrder,pkdStar(pkd,p)->totaltime,dTime,dt);
		    assert(0);
		    }
		}

	    dmstar = dRateCoeff*sqrt(pkdDensity(pkd,p))*pkdMass(pkd,p)*dt;
	    prob = 1.0 - exp(-dmstar/dInitStarMass); 
	    
	    /* Star formation event? */
	    if (rand()<RAND_MAX*prob) {
		float *starpMass = pkdField(starp,pkd->oMass);
		float *pMass = pkdField(p,pkd->oMass);
		pkdCopyParticle(pkd, starp, p);	/* grab copy */
		*pMass -= dInitStarMass;
		*starpMass = dInitStarMass;
/*		pkdStar(pkd,starp)->iGasOrder = p->iOrder;*/
		if (*pMass < 0) {
		    *starpMass += *pMass;
		    *pMass = 0;
		    }
	        if (*pMass < dMinGasMass) {
		    pkdDeleteParticle(pkd, p);
		    (*nDeleted)++;
		    }

                /* Time formed  
		   -- in principle it could have formed any time between dTime-dt and dTime 
		   so dTime-0.5*dt is good -- just check it's less that dtFB */
		if (dt < dtFeedbackDelay) pkdStar(pkd,starp)->fTimer = dTime-dt*.5;
		else pkdStar(pkd,starp)->fTimer = dTime-0.5*dtFeedbackDelay;
		pkdSph(pkd,starp)->u = 1; /* no FB yet */
		
		getClass(pkd,pkdMass(pkd,starp),pkdSoft(pkd,starp),FIO_SPECIES_STAR,starp); /* How do I make a new particle? -- this is bad it rewrites mass and soft for particle */
		/* JW: If class doesn't exist this is very bad -- what is the soft? 
		   For now force softening to exist to get around this */
		(*nFormed)++;
		*dMassFormed += *starpMass;
		pkdNewParticle(pkd, starp);    
		}
	    }
	}

    free(starp);
}


void pkdCooling(PKD pkd, double dTime, double z, int bUpdateState, int bUpdateTable, int bIterateDt, int bIsothermal )
    {
    PARTICLE *p;
    int i;
    SPHFIELDS *sph;
#ifdef COOLING
    COOLPARTICLE cp;  /* Dummy: Not yet fully implemented */
#endif
    double E,dt,ExternalHeating;
  
    pkdClearTimer(pkd,1);
    pkdStartTimer(pkd,1);
  
    assert(pkd->oSph);

    if (bIsothermal)  {
	for (i=0;i<pkdLocal(pkd);++i) {
	    p = pkdParticle(pkd,i);
	    pkdSph(pkd,p)->uDot = 0;
	    }
	}
    else {

#ifdef COOLING
	CoolSetTime( pkd->Cool, dTime, z, bUpdateTable );
#endif	
	if (bIterateDt) { /* Iterate Cooling & dt for each particle */
	    for (i=0;i<pkdLocal(pkd);++i) {
		p = pkdParticle(pkd,i);
		if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
		    if (pkdStar(pkd,p)->fTimer > dTime) continue;
		    sph = pkdSph(pkd,p);
		    ExternalHeating = sph->uDot;
		    for (;;) {
			double uDot;
			
			E = sph->u;
			dt = pkd->param.dDelta/(1<<p->uNewRung); /* Rung Guess */
#ifdef COOLING
			CoolIntegrateEnergyCode(pkd->Cool, &cp, &E, ExternalHeating, p->fDensity, sph->fMetals, p->r, dt);
#endif
			uDot = (E-sph->u)/dt; 
			if (uDot < 0) {
			    double dtNew;
			    int uNewRung;
			    dtNew = pkd->param.dEtaUDot*sph->u/fabs(uDot);
			    uNewRung = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung);
			    if (uNewRung > p->uNewRung) {
				p->uNewRung = uNewRung;
				continue;
				}
			    }
			sph->uDot = uDot;
			break;
			}
		    }
		}
	    }
	else {
	    for (i=0;i<pkdLocal(pkd);++i) {
		p = pkdParticle(pkd,i);
		if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
		    if (pkdStar(pkd,p)->fTimer > dTime) {
			continue;
			}
		    sph = pkdSph(pkd,p);
		    ExternalHeating = sph->uDot;
		    E = sph->u;
		    dt = pkd->param.dDelta/(1<<p->uRung); /* Actual Rung */
#ifdef COOLING
		    CoolIntegrateEnergyCode(pkd->Cool, &cp, &E, ExternalHeating, p->fDensity, sph->fMetals, p->r, dt);
#endif
		    sph->uDot = (E-sph->u)/dt; /* To let us interpolate/extrapolate uPred */
		    }
		}
	    }
	}
    pkdStopTimer(pkd,1);
    }

void pkdCorrectEnergy(PKD pkd, double dTuFac, double z, double dTime, int iDirection )
    {
    PARTICLE *p;
    SPHFIELDS *sph;
    int i;
    double T,E;
#ifdef COOLING
    COOL *cl;
    COOLPARTICLE cp; /* Dummy for now */
#endif

#ifdef COOLING
    cl = pkd->Cool;
    CoolSetTime( cl, dTime, z, 1 );
#endif
    switch(iDirection)  {
    case CORRECTENERGY_IN:
#ifdef COOLING
	for(i=0;i<pkdLocal(pkd);++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsGas(pkd,p)) {
		sph = pkdSph(pkd,p);
		T = sph->u/dTuFac;
		CoolEnergyCodeFromTemp( cl, &cp, &E, &T, p->fDensity, sph->fMetals );
		sph->u = E;
		sph->uPred = E;
		pkdStar(pkd,p)->totaltime = dTime;
		}
	    }
#endif
	break;
	/* Careful using this -- it permanenty converts the thermal energy */
    case CORRECTENERGY_OUT: 
#ifdef COOLING
	for(i=0;i<pkdLocal(pkd);++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsGas(pkd,p)) {
		sph = pkdSph(pkd,p);
		E = sph->u;
		CoolTempFromEnergyCode( cl, &cp, &E, &T, p->fDensity, sph->fMetals );
		sph->u = T*dTuFac;
		sph->uPred = T*dTuFac;
		}
	    }
#endif
	break;
    case CORRECTENERGY_SPECIAL:
#ifdef COOLING
	for(i=0;i<pkdLocal(pkd);++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsGas(pkd,p)) {
		sph = pkdSph(pkd,p);
		T = sph->u/dTuFac; 
		CoolInitEnergyCode( cl, &cp, &E, &T, p->fDensity, sph->fMetals );
		sph->u = E;
		sph->uPred = E;
		}
	    }
#endif
	break;
    default:
	assert(0);
	break;
	}
    }
    
void pkdDensityStep(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, double dEta, double dRhoFac) {
    PARTICLE *p;
    int i;
    double dT;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    dT = dEta/sqrt(pkdDensity(pkd,p)*dRhoFac);
	    p->uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung);
	    }
	}
    }

uint8_t pkdDtToRung(double dT, double dDelta, uint8_t uMaxRung) {
    double dRung;

    assert(dT>0.0);
    dRung = log(dDelta/dT) / M_LN2;
    dRung = (dRung > 0)?dRung:0;
    if (dRung > (double)uMaxRung) return(uMaxRung);
    else return((uint8_t)floor(dRung));
    }


void pkdUpdateRungByTree(PKD pkd,int iRoot,uint8_t uMinRung,int iMaxRung,
uint64_t *nRungCount) {
    KDN *c = pkdTreeNode(pkd,iRoot);
    int i;
    for (i=0;i<=iMaxRung;++i) nRungCount[i] = 0;
    for (i=c->pLower; i<=c->pUpper; ++i) {
	PARTICLE *p = pkdParticle(pkd,i);
	if ( p->uNewRung > iMaxRung ) p->uNewRung = iMaxRung;
	else if (p->uNewRung < uMinRung) p->uNewRung = uMinRung;
	if ( p->uNewRung > p->uRung ) ++p->uRung;
	else if ( p->uNewRung < p->uRung ) --p->uRung;
	nRungCount[p->uRung] += 1;
	}
    }


int pkdUpdateRung(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
		  uint8_t uRung,int iMaxRung,uint64_t *nRungCount) {
    PARTICLE *p;
    int i;
    int iTempRung;
    for (i=0;i<iMaxRung;++i) nRungCount[i] = 0;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( pkdIsActive(pkd,p) ) {
	    if ( p->uNewRung > iMaxRung ) p->uNewRung = iMaxRung;
	    if ( p->uNewRung >= uRung ) p->uRung = p->uNewRung;
	    else if ( p->uRung > uRung) p->uRung = uRung;
	    }
	/*
	** Now produce a count of particles in rungs.
	*/
	nRungCount[p->uRung] += 1;
	}
    iTempRung = iMaxRung;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
    }

void pkdDeleteParticle(PKD pkd, PARTICLE *p) {
    /* p->iOrder = -2 - p->iOrder; JW: Not needed -- just preserve iOrder */
    getClass(pkd,pkdMass(pkd,p),pkdSoft(pkd,p),FIO_SPECIES_LAST,p); /* Special "DELETED" class == FIO_SPECIES_LAST */
    }

void pkdNewParticle(PKD pkd, PARTICLE *p) {
    PARTICLE *newp;

    mdlassert(pkd->mdl,pkd->nLocal < pkd->nStore);
    newp = pkdParticle(pkd,pkd->nLocal);
    pkdCopyParticle(pkd,newp,p);
    newp->iOrder = IORDERMAX;
    pkd->nLocal++;
    }

void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
		  int *nDeltaStar) {
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;

    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for (pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	p = pkdParticle(pkd,pi);
	if (pj < pi)
	    pkdCopyParticle(pkd,pkdParticle(pkd,pj),p);
	if (pkdIsNew(pkd,p)) {
	    ++pj;
	    ++nNew;
	    if (pkdIsGas(pkd, p))
		++ndGas;
	    else if (pkdIsDark(pkd, p))
		++ndDark;
	    else if (pkdIsStar(pkd, p))
		++ndStar;
	    else
		mdlassert(pkd->mdl,0);
	    if (pkdIsActive(pkd,p))
		++pkd->nActive;
	    continue;
	    }
	else if (pkdIsDeleted(pkd,p)) {
	    --newnLocal; /* no idea about type now -- type info lost */
	    --ndGas; /* JW: Hack: assume only gas deleted fix this! */
/*	    if (pkdIsGas(pkd, p))
		--ndGas;
	    else if (pkdIsDark(pkd, p))
		--ndDark;
	    else if (pkdIsStar(pkd, p))
		--ndStar;
	    else
	    mdlassert(pkd->mdl,0);*/
	    if (pkdIsActive(pkd,p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->nLocal = newnLocal;
    }

void pkdNewOrder(PKD pkd,int nStart) {
    PARTICLE *p;
    int pi;

    for (pi=0;pi<pkdLocal(pkd);pi++) {
	p = pkdParticle(pkd,pi);
	if (p->iOrder == IORDERMAX) {
	    p->iOrder = nStart++;
	    }
	}
    }

void
pkdGetNParts(PKD pkd, struct outGetNParts *out )
{
    int pi;
    int n;
    int nGas;
    int nDark;
    int nStar;
    total_t iMaxOrder;
    total_t iOrder;
    PARTICLE *p;
    
    n = 0;
    nGas = 0;
    nDark = 0;
    nStar = 0;
    iMaxOrder = 0;
    for(pi = 0; pi < pkdLocal(pkd); pi++) {
	p = pkdParticle(pkd,pi);
	iOrder = p->iOrder;
	if (iOrder>iMaxOrder) iMaxOrder = iOrder;
	n++;
	if(pkdIsGas(pkd, p)) {
	    ++nGas;
	    }
	else if(pkdIsDark(pkd, p)) {
	    ++nDark;
	    }
	else if(pkdIsStar(pkd, p)) {
	    ++nStar;
	    }
	}
    
    out->n  = n;
    out->nGas = nGas;
    out->nDark = nDark;
    out->nStar = nStar;
    out->nMaxOrder = iMaxOrder;
}


void pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar) {
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    }


void pkdSetRungVeryActive(PKD pkd, int iRung) {
    /* Remember, the first very active particle is at iRungVeryActive + 1 */
    pkd->uRungVeryActive = iRung;
    }

int pkdIsGas(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_SPH;
    }

int pkdIsDark(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_DARK;
    }

int pkdIsStar(PKD pkd,PARTICLE *p) {
    return pkdSpecies(pkd,p) == FIO_SPECIES_STAR;
    }

void pkdInitRelaxation(PKD pkd) {
    PARTICLE *p;
    double *pRelax;
    int i;

    assert(pkd->oRelaxation);
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	pRelax = pkdField(p,pkd->oRelaxation);
	*pRelax = 0.0;
	}
    }

double pkdTotalMass(PKD pkd) {
    PARTICLE *p;
    double m;
    int i,n;

    m = 0.0;
    n = pkdLocal(pkd);
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	m += pkdMass(pkd,p);
	}
    return m;
    }

/*
** This function checks the predicate and returns a new value based on the flags.
** setIfTrue:    >0 -> return true if the predicate is true
**               <0 -> "clear if true"
** clearIfFalse: >0 -> return false if the predicate is false
**               <0 -> "set if false"
** A value of zero for either results in no action for the "IfTrue" or "IfFalse" flags.
** Conflicting options (e.g., setIfTrue and setIfFalse) result in a toggle.
*/
static inline int isSelected( int predicate, int setIfTrue, int clearIfFalse, int value ) {
    int s = (predicate&(setIfTrue>0)) | (~predicate&(clearIfFalse<0));
    int c = (predicate&(setIfTrue<0)) | (~predicate&(clearIfFalse>0));
    return (~s&~c&value) | (s&~(c&value));
    }

int pkdSelSrcAll(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    for( i=0; i<n; i++ ) pkdParticle(pkd,i)->bSrcActive = 1;
    return n;
    }
int pkdSelDstAll(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    for( i=0; i<n; i++ ) pkdParticle(pkd,i)->bDstActive = 1;
    return n;
    }

int pkdSelSrcGas(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	if (pkdIsGas(pkd,p)) p->bSrcActive = 1; else p->bSrcActive = 0;
	}
    return n;
    }

int pkdSelDstGas(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	if (pkdIsGas(pkd,p)) p->bDstActive = 1; else p->bDstActive = 0;
	}
    return n;
    }

int pkdSelSrcStar(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	if (pkdIsStar(pkd,p)) p->bSrcActive = 1; else p->bSrcActive = 0;
	}
    return n;
    }

int pkdSelDstStar(PKD pkd, int bFB, double dTimeFB) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;

    if (bFB) {
	for( i=0; i<n; i++ ) {
	    p=pkdParticle(pkd,i);
	    if (pkdIsStar(pkd,p) && pkdIsActive(pkd,p)) {
		double dtp = pkd->param.dDelta/(1<<p->uRung);
		if (dTimeFB-dtp < *pkd_Timer(pkd,p)) p->bDstActive = 1; 
		else p->bDstActive = 0;
		}
	    else p->bDstActive = 0;
	    }
	}
    else {
	for( i=0; i<n; i++ ) {
	    p=pkdParticle(pkd,i);
	    if (pkdIsStar(pkd,p)) p->bDstActive = 1; else p->bDstActive = 0;
	    }
	}
    return n;
    }

int pkdSelSrcDeleted(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	if (pkdIsDeleted(pkd,p)) p->bSrcActive = 1; else p->bSrcActive = 0;
	}
    return n;
    }

int pkdSelDstDeleted(PKD pkd) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	if (pkdIsDeleted(pkd,p)) p->bDstActive = 1; else p->bDstActive = 0;
	}
    return n;
    }

int pkdSelSrcMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double m;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	m = pkdMass(pkd,p);
	p->bSrcActive = isSelected((m >= dMinMass && m <=dMaxMass),setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelDstMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double m;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	m = pkdMass(pkd,p);
	p->bDstActive = isSelected((m >= dMinMass && m <=dMaxMass),setIfTrue,clearIfFalse,p->bDstActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelSrcById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	p->bSrcActive = isSelected((p->iOrder >= idStart && p->iOrder <= idEnd),setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelDstById(PKD pkd,uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	p->bDstActive = isSelected((p->iOrder >= idStart && p->iOrder <= idEnd),setIfTrue,clearIfFalse,p->bDstActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }


int pkdSelSrcPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    VELSMOOTH *pvel;
    float density;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	pvel = pkdField(p,pkd->oVelSmooth);
	density = pkdDensity(pkd,p) * pow(pvel->veldisp2,-1.5);
	p->bSrcActive = isSelected((density >= dMinDensity && density <=dMaxDensity),setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelDstPhaseDensity(PKD pkd,double dMinDensity, double dMaxDensity, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    VELSMOOTH *pvel;
    float density;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	pvel = pkdField(p,pkd->oVelSmooth);
	density = pkdDensity(pkd,p) * pow(pvel->veldisp2,-1.5);
	p->bDstActive = isSelected((density >= dMinDensity && density <=dMaxDensity),setIfTrue,clearIfFalse,p->bDstActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelSrcBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,j,n,nSelected;
    int predicate;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	predicate = 1;
	for(j=0; j<3; j++ ) {
	    double dx = dCenter[j] - pkdPos(pkd,p,j);
	    predicate = predicate && dx < dSize[j] && dx >= -dSize[j];
	    }
	p->bSrcActive = isSelected(predicate,setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelDstBox(PKD pkd,double *dCenter, double *dSize, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    int i,j,n,nSelected;
    int predicate;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	predicate = 1;
	for(j=0; j<3; j++ )
	    predicate = predicate && fabs(dCenter[j] - pkdPos(pkd,p,j)) <= dSize[j];
	p->bDstActive = isSelected(predicate,setIfTrue,clearIfFalse,p->bDstActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelSrcSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double d2,dx,dy,dz,dRadius2;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    dRadius2 = dRadius*dRadius;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	dx = r[0] - pkdPos(pkd,p,0);
	dy = r[1] - pkdPos(pkd,p,1);
	dz = r[2] - pkdPos(pkd,p,2);

	d2 = dx*dx + dy*dy + dz*dz;
	p->bSrcActive = isSelected((d2<=dRadius2),setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelDstSphere(PKD pkd,double *r, double dRadius, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double d2,dx,dy,dz,dRadius2;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    dRadius2 = dRadius*dRadius;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	dx = r[0] - pkdPos(pkd,p,0);
	dy = r[1] - pkdPos(pkd,p,1);
	dz = r[2] - pkdPos(pkd,p,2);

	d2 = dx*dx + dy*dy + dz*dz;
	p->bDstActive = isSelected((d2<=dRadius2),setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }

/*
** Select all particles that fall inside a cylinder between two points P1 and P2
** with radius dRadius.
*/
int pkdSelSrcCylinder(PKD pkd,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double dCyl[3], dPart[3];
    double dLength2, dRadius2, dL2;
    double pdotr;
    int i,j,n,nSelected;
    int predicate;

    dRadius2 = dRadius*dRadius;
    dLength2 = 0.0;
    for( j=0;j<3;j++ ) {
	dCyl[j] = dP2[j] - dP1[j];
	dLength2 += dCyl[j] * dCyl[j];
	}
    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	pdotr = 0.0;
	for( j=0;j<3;j++ ) {
	    dPart[j] = pkdPos(pkd,p,j) - dP1[j];
	    pdotr += dPart[j] * dCyl[j];
	    }

	if ( pdotr < 0.0 || pdotr > dLength2 ) predicate = 0;
	else {
	    dL2 = dPart[0]*dPart[0] + dPart[1]*dPart[1] + dPart[2]*dPart[2] - pdotr*pdotr/dLength2;
	    predicate = (dL2 <= dRadius2);
	    }
	p->bSrcActive = isSelected(predicate,setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bSrcActive ) nSelected++;
	}
    return nSelected;
    }

/*
** Select all particles that fall inside a cylinder between two points P1 and P2
** with radius dRadius.
*/
int pkdSelDstCylinder(PKD pkd,double *dP1, double *dP2, double dRadius,
		      int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double dCyl[3], dPart[3];
    double dLength2, dRadius2, dL2;
    double pdotr;
    int i,j,n,nSelected;
    int predicate;

    dRadius2 = dRadius*dRadius;
    dLength2 = 0.0;
    for( j=0;j<3;j++ ) {
	dCyl[j] = dP2[j] - dP1[j];
	dLength2 += dCyl[j] * dCyl[j];
	}

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = pkdParticle(pkd,i);
	pdotr = 0.0;
	for( j=0;j<3;j++ ) {
	    dPart[j] = pkdPos(pkd,p,j) - dP1[j];
	    pdotr += dPart[j] * dCyl[j];
	    }

	if ( pdotr < 0.0 || pdotr > dLength2 ) predicate = 0;
	else {
	    dL2 = dPart[0]*dPart[0] + dPart[1]*dPart[1] + dPart[2]*dPart[2] - pdotr*pdotr/dLength2;
	    predicate = (dL2 <= dRadius2);
	    }
	p->bDstActive = isSelected(predicate,setIfTrue,clearIfFalse,p->bSrcActive);
	if ( p->bDstActive ) nSelected++;
	}
    return nSelected;
    }

int pkdSelSrcGroup(PKD pkd, int iGroup) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	p->bSrcActive = *pkdGroup(pkd,p)==iGroup;
	}
    return n;
    }

int pkdSelDstGroup(PKD pkd, int iGroup) {
    int i;
    int n=pkdLocal(pkd);
    PARTICLE *p;
    for( i=0; i<n; i++ ) {
	p=pkdParticle(pkd,i);
	p->bDstActive = *pkdGroup(pkd,p)==iGroup;
	}
    return n;
    }

/*
**  Find the source particle with the deepest potential
*/
int pkdDeepestPot(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
    double *r, float *fPot) {
    int i,n,nChecked;
    PARTICLE *p, *pLocal;
    float *pPot, *pPotLocal;

    assert(pkd->oPotential);

    n = pkdLocal(pkd);
    pLocal = pkdParticle(pkd,0);
    pPotLocal = pkdPot(pkd,pLocal);
    nChecked = 0;
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsSrcActive(p,uRungLo,uRungHi)) {
	    nChecked++;
	    pPot = pkdPot(pkd,p);
	    if ( *pPot < *pPotLocal ) {
		pLocal = p;
		pPotLocal = pkdPot(pkd,pLocal);
		}
	    }
	}
    r[0] = pkdPos(pkd,pLocal,0);
    r[1] = pkdPos(pkd,pLocal,1);
    r[2] = pkdPos(pkd,pLocal,2);
    *fPot= *pPotLocal;
    return nChecked;
    }

void pkdOutPsGroup(PKD pkd,char *pszFileName,int iType)
{
    FILE *fp;
    int i,j,nout,lStart;

    if (iType == OUT_PSGROUP_STATS) {
	fp = fopen(pszFileName,"a+");
	assert(fp != NULL);
	struct psGroup *gd = pkd->psGroupTable.pGroup;

	for (i=1;i<pkd->psGroupTable.nGroups;++i)
	{
	    if (gd[i].iPid != pkd->idSelf) continue;
	    fprintf(fp,"%d",gd[i].iGlobalId);
	    fprintf(fp," %10llu",gd[i].nTotal);
	    fprintf(fp," %12.8e",gd[i].fMass);
	    fprintf(fp," %12.8e",gd[i].fRMSRadius);
	    fprintf(fp," %12.8e",gd[i].r[0]);
	    fprintf(fp," %12.8e",gd[i].r[1]);
	    fprintf(fp," %12.8e",gd[i].r[2]);
	    fprintf(fp," %12.8e",gd[i].v[0]);
	    fprintf(fp," %12.8e",gd[i].v[1]);
	    fprintf(fp," %12.8e",gd[i].v[2]);
#if 0
	    fprintf(fp,"%.11g ",pkd->groupData[i].rcom[0]);
	    fprintf(fp,"%.11g ",pkd->groupData[i].rcom[1]);
	    fprintf(fp,"%.11g ",pkd->groupData[i].rcom[2]);
	    fprintf(fp,"%.11g ",pkd->groupData[i].r[0]);
	    fprintf(fp,"%.11g ",pkd->groupData[i].r[1]);
	    fprintf(fp,"%.11g ",pkd->groupData[i].r[2]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[0]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[1]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[2]);
#endif
	    fprintf(fp,"\n");
	}
	if (fclose(fp) == EOF)
	{
	    perror("pkdOutGroup: could not close file");
	    exit(1);
	}
    }
    else
	assert(0);
}

