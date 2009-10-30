#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <inttypes.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <sys/time.h>
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "pkd.h"
#include "ewald.h"
#include "walk.h"
#include "grav.h"
#include "mdl.h"
#include "tipsydefs.h"
#include "ssio.h"

#include "parameters.h"
#include "cosmo.h"

#ifdef USE_HDF5
#include "iohdf5.h"
#endif

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif

const char *pkd_module_id = "$Id$";
const char *pkd_h_module_id = PKD_H_MODULE_ID;

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
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
	gettimeofday(&tv,NULL);
	pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
	    {
	    struct rusage ru;

	    getrusage(0,&ru);
	    pkd->ti[iTimer].system_stamp = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
	    }
	}
    }


void pkdStopTimer(PKD pkd,int iTimer) {
    double sec;
    struct timeval tv;

    sec = -pkd->ti[iTimer].stamp;
    pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
    sec += pkd->ti[iTimer].stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].sec += sec;

    sec = -pkd->ti[iTimer].wallclock_stamp;
    gettimeofday( &tv, NULL );
    pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double)tv.tv_usec;
    sec += pkd->ti[iTimer].wallclock_stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].wallclock_sec += sec;

#ifndef _CRAYMPP
	{
	struct rusage ru;

	sec = -pkd->ti[iTimer].system_stamp;
	getrusage(0,&ru);
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

void pkdAllocateTopTree(PKD pkd,int nCell) {
    if (pkd->kdTopPRIVATE != NULL) free(pkd->kdTopPRIVATE);
    pkd->kdTopPRIVATE = malloc(nCell*pkd->iTreeNodeSize);
    assert(pkd->kdTopPRIVATE != NULL);
    }

void pkdInitialize(
    PKD *ppkd,MDL mdl,int nStore,int nBucket,int nTreeBitsLo, int nTreeBitsHi,
    int iCacheSize,FLOAT *fPeriod,uint64_t nDark,uint64_t nGas,uint64_t nStar,
    uint64_t mMemoryModel) {
    PKD pkd;
    int j,ism;

#define RANDOM_SEED 1
    srand(RANDOM_SEED);

    pkd = (PKD)malloc(sizeof(struct pkdContext));
    mdlassert(mdl,pkd != NULL);
    pkd->mdl = mdl;
    pkd->idSelf = mdlSelf(mdl);
    pkd->nThreads = mdlThreads(mdl);
    pkd->kdNodeListPRIVATE = NULL;
    pkd->pStorePRIVATE = NULL;
    pkd->nStore = nStore;
    pkd->nLocal = 0;
    pkd->nDark = nDark;
    pkd->nGas = nGas;
    pkd->nStar = nStar;
/*    pkd->nMaxOrderGas = nGas;
    pkd->nMaxOrderDark = nGas + nDark;
    pkd->nMaxOrder = nGas+nDark+nStar;  JW: Deprecate this: Probably wrong if Order in input file */
    pkd->nRejects = 0;
    for (j=0;j<3;++j) {
	pkd->fPeriod[j] = fPeriod[j];
	}

    pkd->uMinRungActive  = 0;
    pkd->uMaxRungActive  = 255;
    pkd->uRungVeryActive = 255;


    /*
    ** Calculate the amount of memory (size) of each particle.  This is the
    ** size of a base particle (PARTICLE), plus any extra fields as defined
    ** by the current memory model.  Fields need to be added in order of
    ** descending size (i.e., doubles & int64 and then float & int32)
    */
    pkd->iParticleSize = sizeof(PARTICLE);
    pkd->iTreeNodeSize = sizeof(KDN);

    if ( mMemoryModel & PKD_MODEL_VELOCITY )
	pkd->oVelocity = pkdParticleAddDouble(pkd,3);
    else
	pkd->oVelocity = 0;

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

    if ( mMemoryModel & PKD_MODEL_HERMITE )
	pkd->oHermite = pkdParticleAddStruct(pkd,sizeof(HERMITEFIELDS));
    else
	pkd->oHermite = 0;

    if ( mMemoryModel & PKD_MODEL_VELSMOOTH )
	pkd->oVelSmooth = pkdParticleAddStruct(pkd,sizeof(VELSMOOTH));
    else
	pkd->oVelSmooth = 0;

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

    if ( mMemoryModel & PKD_MODEL_GROUPS ) {
	pkd->oGroup = pkdParticleAddInt32(pkd,1);
	}
    else {
	pkd->oGroup = 0;
	}

    /*
    ** Tree node memory models
    */
    if ( mMemoryModel & PKD_MODEL_TREE_MOMENT )
	pkd->oNodeMom = pkdNodeAddStruct(pkd,sizeof(MOMR));
    else
	pkd->oNodeMom = 0;

    if ( mMemoryModel & PKD_MODEL_VELOCITY )
	pkd->oNodeVelocity = pkdNodeAddDouble(pkd,3);
    else
	pkd->oNodeVelocity = 0;

    /* The acceleration is required for the new time step criteria */
#ifdef LOCAL_EXPANSION
    pkd->oNodeAcceleration = pkdNodeAddDouble(pkd,3);
#else
    pkd->oNodeAcceleration = 0;
#endif

    /*
    ** N.B.: Update pkdMaxNodeSize in pkd.h if you add fields.  We need to
    **       know the size of a node when setting up the pst.
    */
    assert(pkdNodeSize(pkd)<=pkdMaxNodeSize());

    /*
    ** Allocate the main particle store.
    ** Need to use mdlMalloc() since the particles will need to be
    ** visible to all other processors thru mdlAquire() later on.
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
    pkd->pStorePRIVATE = mdlMalloc(pkd->mdl,(nStore+1)*pkdParticleSize(pkd));
    mdlassert(mdl,pkd->pStorePRIVATE != NULL);
    pkd->pTempPRIVATE = malloc(pkdParticleSize(pkd));
    mdlassert(mdl,pkd->pTempPRIVATE != NULL);

#ifdef MDL_CACHE_SIZE
    if ( iCacheSize > 0 ) mdlSetCacheSize(pkd->mdl,iCacheSize);
#endif

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
    pkd->pLite = malloc((nStore+1)*sizeof(PLITE));
    mdlassert(mdl,pkd->pLite != NULL);
    pkd->nNodes = 0;
    pkd->kdTopPRIVATE = NULL;
    /*
    ** Ewald stuff!
    */
    pkd->ew.nMaxEwhLoop = 100;
    pkd->ew.ewt = malloc(pkd->ew.nMaxEwhLoop*sizeof(EWT));
    mdlassert(mdl,pkd->ew.ewt != NULL);
    *ppkd = pkd;
    /*
    ** Tree walk stuff.
    */
#ifdef LOCAL_EXPANSION
    ilpInitialize(&pkd->ilp);
    ilcInitialize(&pkd->ilc);
#else
    pkd->nMaxPart = 10000;
    pkd->ilp = malloc(pkd->nMaxPart*sizeof(ILP));
    assert(pkd->ilp != NULL);
    pkd->nMaxCell = 1000;
    pkd->ilc = malloc(pkd->nMaxCell*sizeof(ILC));
    assert(pkd->ilc != NULL);
#endif
    /*
    ** Allocate Checklist.
    */
    pkd->nMaxCheck = 10000;
    pkd->Check = malloc(pkd->nMaxCheck*sizeof(CELT));
    assert(pkd->Check != NULL);
    /*
    ** Allocate the stack.
    */
    pkd->nMaxStack = 30;
    pkd->S = malloc(pkd->nMaxStack*sizeof(CSTACK));
    assert(pkd->S != NULL);
    for (ism=0;ism<pkd->nMaxStack;++ism) {
	pkd->S[ism].Check = malloc(pkd->nMaxCheck*sizeof(CELT));
	assert(pkd->S[ism].Check != NULL);
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

    pkd->Cool = CoolInit();
    }


void pkdFinish(PKD pkd) {
    int ism;
    int i;

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
#ifdef LOCAL_EXPANSION
    ilpFinish(pkd->ilp);
    ilcFinish(pkd->ilc);
#else
    free(pkd->ilp);
    free(pkd->ilc);
#endif
    /*
    ** Free checklist.
    */
    free(pkd->Check);
    /*
    ** Free Stack.
    */
    for (ism=0;ism<pkd->nMaxStack;++ism) {
	free(pkd->S[ism].Check);
	}
    free(pkd->S);
    if (pkd->kdTopPRIVATE) free(pkd->kdTopPRIVATE);
    free(pkd->ew.ewt);
    free(pkd->pClass);
    mdlFree(pkd->mdl,pkd->pStorePRIVATE);
    free(pkd->pLite);
    free(pkd->piActive);
    free(pkd->piInactive);
    csmFinish(pkd->param.csm);
    free(pkd);
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
	fprintf(stderr,"New class %d: %g %g %d\n",i,fMass,fSoft,eSpecies);
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
#endif
    off_t lStart;
    int iErr;
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


void IOCheck(int nout) {
    if (nout != 1) {
	perror("IOCheck failed");
	exit(errno);
	}
    }


#ifdef USE_GRAFIC
void pkdGenerateIC(PKD pkd, GRAFICCTX gctx,  int iDim,
		   double fSoft, double fMass, int bComove) {
    PARTICLE *p;
    int i, j, k, d, pi, n1, n2, n3;
    double dx;
    double dvFac, a;

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
		p->v[d] = graficGetVelocity(gctx,i,j,k) * dvFac;
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
    double *v, dummyv[3];
    float fMass, fSoft;
    FIO_SPECIES eSpecies;
    uint64_t iOrder;

    mdlassert(pkd->mdl,fio != NULL);

    if (pkd->oStar) {
	/* Make sure star class established -- how do all procs know of these classes? How do we ensure they agree on the class identifiers? */
	p = pkdParticle(pkd,pkd->nLocal);
	getClass(pkd,0,0,FIO_SPECIES_STAR,p);
	fprintf(stderr,"INFO: Dummy star to establish star class\n");
	}

    fioSeek(fio,iFirst,FIO_SPECIES_ALL);
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,pkd->nLocal+i);
	/*
	** General initialization.
	*/
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	p->fDensity = 0.0;
	p->fBall = 0.0;
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
	if (pkd->oVelocity) v = pkdVel(pkd,p);
	else v = dummyv;

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
//	    pStar->iGasOrder = IORDERMAX;
	    }
	else pStar = NULL;

	eSpecies = fioSpecies(fio);
	switch(eSpecies) {
	case FIO_SPECIES_SPH:
	    assert(pSph); /* JW: Could convert to dark ... */
	    assert(dTuFac>0.0);
	    fioReadSph(fio,&iOrder,p->r,v,&fMass,&fSoft,pPot,
		       &p->fDensity/*?*/,&pSph->u,&pSph->fMetals);
/*	    if ((iOrder%1000)==0) printf("%d: %g %g %g\n",iOrder,pSph->u,dTuFac,pSph->u*dTuFac);*/
	    pSph->u *= dTuFac; /* Can't do precise conversion until density known */
	    pSph->uPred = pSph->u;
	    pSph->fMetalsPred = pSph->fMetals;
	    pSph->vPred[0] = v[0]*dvFac;
	    pSph->vPred[1] = v[1]*dvFac;
	    pSph->vPred[2] = v[2]*dvFac; /* density, divv, BalsaraSwitch, c set in smooth */
	    break;
	case FIO_SPECIES_DARK:
	    fioReadDark(fio,&iOrder,p->r,v,&fMass,&fSoft,pPot);
	    break;
	case FIO_SPECIES_STAR:
	    assert(pStar && pSph);
	    fioReadStar(fio,&iOrder,p->r,v,&fMass,&fSoft,pPot,
			&pSph->fMetals,&pStar->fTimer);
	    pSph->vPred[0] = v[0]*dvFac;
	    pSph->vPred[1] = v[1]*dvFac;
	    pSph->vPred[2] = v[2]*dvFac;
	    break;
	default:
	    fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
	    assert(0);
	    }
	p->iOrder = iOrder;
	for(j=0;j<3;++j) v[j] *= dvFac;
	getClass(pkd,fMass,fSoft,eSpecies,p);
	}
    
    pkd->nLocal += nLocal;
    pkd->nActive += nLocal;
    }

#ifdef USE_MDL_IO
void pkdIOInitialize( PKD pkd, int nLocal) {
    int i, j;
    PARTICLE *p;

    pkd->nLocal = pkd->nActive = nLocal;

    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	p->iClass = 0;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	if ( pkd->oAcceleration ) {
	    float *a = pkdAccel(pkd,p);
	    for (j=0;j<3;++j) {
		a[j] = 0.0;
		}
	    }
	}

    }
#endif

void pkdCalcBound(PKD pkd,BND *pbnd) {
    double dMin[3],dMax[3];
    PARTICLE *p;
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->nLocal > 0);
    p = pkdParticle(pkd,i);
    for (j=0;j<3;++j) {
	dMin[j] = p->r[j];
	dMax[j] = p->r[j];
	}
    for (++i;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	pkdMinMax(p->r,dMin,dMax);
	}
    for (j=0;j<3;++j) {
	pbnd->fCenter[j] = pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pbnd->fMax[j] = pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
    }


void pkdEnforcePeriodic(PKD pkd,BND *pbnd) {
    PARTICLE *p;
    int i,j;

    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) {
	    if (p->r[j] < pbnd->fCenter[j] - pbnd->fMax[j]) p->r[j] += 2*pbnd->fMax[j];
	    else if (p->r[j] >= pbnd->fCenter[j] + pbnd->fMax[j]) p->r[j] -= 2*pbnd->fMax[j];
	    /*
	    ** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the 
	    ** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
	    ** to the specified periodic box in a gross way.
	    */
	    mdlassert(pkd->mdl,((p->r[j] >= pbnd->fCenter[j] - pbnd->fMax[j])&&
				(p->r[j] < pbnd->fCenter[j] + pbnd->fMax[j])));

	    }
	}
    }


/*
** x, y and z must have range [1,2) !
*/
uint64_t hilbert3d(float x,float y,float z) {
    uint64_t s = 0;
    uint32_t m,ux,uy,uz,ut;

    ux = (UNION_CAST(x,float,uint32_t))>>2;
    uy = (UNION_CAST(y,float,uint32_t))>>2;
    uz = (UNION_CAST(z,float,uint32_t))>>2;

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

#if 0
void pkdPeanoHilbertCount(PKD pkd) {
    PARTICLE *p;
    PLITEDD *pl = (PLITEDD *)pkd->pLite;
    uint64_t uMask;
    float x,y,z;
    int i,j,bits,iShift;

    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	/*
	** For now we just assume the particles are coming from a standard
	** cosmological box. We scale the volume by a factor of 0.99 to be
	** certain that fast particles are still captured by the domain
	** decomposition.
	*/
	x = 0.99*p->r[0] + 1.5;
	if (x < 1.0) x = 1.0;
	else if (x >= 2.0) x = 2.0;
	y = 0.99*p->r[1] + 1.5;
	if (y < 1.0) y = 1.0;
	else if (y >= 2.0) y = 2.0;
	z = 0.99*p->r[2] + 1.5;
	if (z < 1.0) z = 1.0;
	else if (z >= 2.0) z = 2.0;
	pl[i].uKey = hilbert3d(x,y,z);
	pl[i].i = i;
	}
    /*
    ** Now create a local count of particles in the PeanoHilbert cells.
    */
    iShift = 63 - bits;
    uMask = (1<<(bits+1)) - 1;
    for (j=0;j<pkd->nPHCount;++j) pkd->auPHCount[j] = 0;
    for (i=0;i<pkd->nLocal;++i) {

	}

    /*
    ** Example of mdlReduce from Doug's code.
    mdlReduce(pkd->mdl,limg,img,N,MPI_FLOAT,MPI_SUM,0);
    */
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
	    if (p->r[d] < fSplit) *pnLow += 1;
	    else *pnHigh += 1;
	    }
	}
    }

double pkdTester(PKD pkd,int i) {
    PARTICLE *p = pkdParticle(pkd,i);
    return p->r[2];
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
	       pi->r[d] >= fSplit,pj->r[d] < fSplit);
    return(i);
    }


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j) {
    PARTICLE *pi, *pj;
    pi = pkdParticle(pkd,i);
    pj = pkdParticle(pkd,j);
    PARTITION(pi<pj,pi<=pj,
	       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
	       pkdSwapParticle(pkd,pi,pj),
	       pi->r[d] < fSplit,pj->r[d] >= fSplit);
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
		       (pi->r[d] < fSplit2 || pi->r[d] >= fSplit1) &&
		       !pkdIsVeryActive(pkd,pi),
		       (pj->r[d] >= fSplit2 && pj->r[d] < fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] < fSplit2 || pi->r[d] >= fSplit1) ||
		       pkdIsVeryActive(pkd,pi),
		       (pj->r[d] >= fSplit2 && pj->r[d] < fSplit1) &&
		       !pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] < fSplit2 || pi->r[d] >= fSplit1),
		       (pj->r[d] >= fSplit2 && pj->r[d] < fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] < fSplit2 && pi->r[d] >= fSplit1) &&
		       !pkdIsVeryActive(pkd,pi),
		       (pj->r[d] >= fSplit2 || pj->r[d] < fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] < fSplit2 && pi->r[d] >= fSplit1) ||
		       pkdIsVeryActive(pkd,pi),
		       (pj->r[d] >= fSplit2 || pj->r[d] < fSplit1) &&
		       !pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] < fSplit2 && pi->r[d] >= fSplit1),
		       (pj->r[d] >= fSplit2 || pj->r[d] < fSplit1));
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
		       (pi->r[d] >= fSplit2 && pi->r[d] < fSplit1) ||
		       pkdIsVeryActive(pkd,pi),
		       (pj->r[d] < fSplit2 || pj->r[d] >= fSplit1) &&
		       !pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] >= fSplit2 && pi->r[d] < fSplit1) &&
		       !pkdIsVeryActive(pkd,pi),
		       (pj->r[d] < fSplit2 || pj->r[d] >= fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] >= fSplit2 && pi->r[d] < fSplit1),
		       (pj->r[d] < fSplit2 || pj->r[d] >= fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] >= fSplit2 || pi->r[d] < fSplit1) ||
		       pkdIsVeryActive(pkd,pi),
		       (pj->r[d] < fSplit2 && pj->r[d] >= fSplit1) &&
		       !pkdIsVeryActive(pkd,pj));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] >= fSplit2 || pi->r[d] < fSplit1) &&
		       !pkdIsVeryActive(pkd,pi),
		       (pj->r[d] < fSplit2 && pj->r[d] >= fSplit1) ||
		       pkdIsVeryActive(pkd,pj));
	    }
	else {
	    PARTITION(pi<pj,pi<=pj,
		       pi=pkdParticle(pkd,++i),pj=pkdParticle(pkd,--j),
		       pkdSwapParticle(pkd,pi,pj),
		       (pi->r[d] >= fSplit2 || pi->r[d] < fSplit1),
		       (pj->r[d] < fSplit2 && pj->r[d] >= fSplit1));
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


int pkdColRejects_Active_Inactive(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide) {
    int nSplit,nSplitInactive,iRejects,i,j;

    mdlassert(pkd->mdl,pkd->nRejects == 0);
    if (iSplitSide) {
	nSplit = pkdLowerPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
	}
    else {
	nSplit = pkdUpperPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
	}
    if (iSplitSide) {
	nSplitInactive = pkdLowerPart(pkd,d,fSplitInactive,
				      pkdActive(pkd),pkdLocal(pkd)-1);
	}
    else {
	nSplitInactive = pkdUpperPart(pkd,d,fSplitInactive,
				      pkdActive(pkd),pkdLocal(pkd)-1);
	}
    nSplitInactive -= pkdActive(pkd);
    /*
    ** Now do some fancy rearrangement.
    */
    i = nSplit;
    j = nSplit+nSplitInactive;
    while (j < pkdActive(pkd) + nSplitInactive) {
	pkdSwapParticle(pkd,pkdParticle(pkd,i),pkdParticle(pkd,j));
	++i; ++j;
	}
    pkd->nRejects = pkdLocal(pkd) - nSplit - nSplitInactive;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nActive = nSplit;
    pkd->nLocal = nSplit + nSplitInactive;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkdCopyParticle(pkd,pkdParticle(pkd,iRejects+i),pkdParticle(pkd,pkd->nLocal+i));
    return(pkd->nRejects);
    }


int pkdColRejects(PKD pkd,int nSplit) {
    int iRejects,i;

    mdlassert(pkd->mdl,pkd->nRejects == 0);

    pkd->nRejects = pkdLocal(pkd) - nSplit;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nLocal = nSplit;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkdCopyParticle(pkd,pkdParticle(pkd,iRejects+i),pkdParticle(pkd,pkd->nLocal+i));
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
    int nSplit,iRejects,i;

    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    pkd->nRejects = pkdLocal(pkd) - nSplit;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nLocal = nSplit;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkdCopyParticle(pkd,pkdParticle(pkd,iRejects+i),pkdParticle(pkd,pkd->nLocal+i));
    return(pkd->nRejects);
    }


int cmpParticles(const void *pva,const void *pvb) {
    PARTICLE *pa = (PARTICLE *)pva;
    PARTICLE *pb = (PARTICLE *)pvb;

    return(pa->iOrder - pb->iOrder);
    }


void pkdLocalOrder(PKD pkd) {
    qsort(pkdParticleBase(pkd),pkdLocal(pkd),pkdParticleSize(pkd),cmpParticles);
    }

int pkdUnpackIO(PKD pkd,
		PIO *io,
		int nMax,
		local_t *iIndex,
		total_t iMinOrder,
		total_t iMaxOrder,
		double dvFac ) {
    int i,d;

    assert( pkd != NULL );

    for ( i=0; i<nMax; i++ ) {
	local_t I = *iIndex + i;
	PARTICLE *p = pkdParticle(pkd,I);
	float *a, dummya[3];
	float *pPot, dummypot;
	double *v, dummyv[3];

	if ( pkd->oPotential) pPot = pkdPot(pkd,p);
	else pPot = &dummypot;
	if (pkd->oVelocity) v = pkdVel(pkd,p);
	else v = dummyv;
	if ( pkd->oAcceleration ) a = pkdAccel(pkd,p);
	else a = dummya;

	for ( d=0; d<3; d++ ) {
	    p->r[d]  = io[i].r[d];
	    v[d]  = io[i].v[d] * dvFac; /*FIXME: times??*/
	    a[d]  = 0.0;
	    }
	p->iOrder = io[i].iOrder;
	getClass(pkd,io[i].fMass,io[i].fSoft,FIO_SPECIES_DARK,p);
	p->fDensity = io[i].fDensity;
	*pPot = io[i].fPot;
	}
    *iIndex += nMax;

    return pkd->nLocal - *iIndex;
    }

/*
** This routine will pack at most nMax particles into an array of packed
** particles (io).  It scans the particle list from *iIndex to the end
** packing only particles with iOrder in the range [iMinOrder,iMaxOrder).
** When this routine is done, it will return 0 for the number of particles
** packed (and continues to return 0 if called again).
*/
int pkdPackIO(PKD pkd,
	      PIO *io,
	      int nMax,
	      local_t *iIndex,
	      total_t iMinOrder,
	      total_t iMaxOrder,
	      double dvFac ) {
    PARTICLE *p;
    float *pPot, dummypot;
    double *v, dummyv[3];
    local_t i;
    int nCopied, d;

    dummyv[0] = dummyv[1] = dummyv[2] = 0.0;
    dummypot = 0.0;

    mdlassert(pkd->mdl,*iIndex<=pkd->nLocal);

    for ( i=*iIndex,nCopied=0; nCopied < nMax && i < pkd->nLocal; i++ ) {
	p = pkdParticle(pkd,i);
	if ( pkd->oPotential) pPot = pkdPot(pkd,p);
	else pPot = &dummypot;
	if (pkd->oVelocity) v = pkdVel(pkd,p);
	else v = dummyv;

	/* Not a particle of interest? */
	if ( p->iOrder<iMinOrder || p->iOrder>=iMaxOrder)
	    continue;

	/* We should have certain special cases here */
	mdlassert( pkd->mdl, pkdIsDark(pkd,p) );

	for ( d=0; d<3; d++ ) {
	    io[nCopied].r[d] = p->r[d];
	    io[nCopied].v[d] = v[d] * dvFac;
	    }
	io[nCopied].iOrder= p->iOrder;
	io[nCopied].fMass = pkdMass(pkd,p);
	io[nCopied].fSoft = pkdSoft(pkd,p);
	io[nCopied].fDensity = p->fDensity;
	io[nCopied].fPot = *pPot;
	nCopied++;
	}
    *iIndex = i;
    return nCopied;
    }

#ifdef USE_HDF5
void pkdWriteHDF5(PKD pkd, IOHDF5 io, IOHDF5V ioDen, IOHDF5V ioPot, double dvFac) {
    PARTICLE *p;
    float *pPot, dummypot;
    FLOAT v[3], fSoft, fMass;
    double *pv, dummyv[3];
    int i;

    dummyv[0] = dummyv[1] = dummyv[2] = 0.0;
    dummypot = 0.0;
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	if ( pkd->oPotential) pPot = pkdPot(pkd,p);
	else pPot = &dummypot;
	if (pkd->oVelocity) pv = pkdVel(pkd,p);
	else pv = dummyv;

	v[0] = pv[0] * dvFac;
	v[1] = pv[1] * dvFac;
	v[2] = pv[2] * dvFac;
	fSoft = pkdSoft0(pkd,p);
	fMass = pkdMass(pkd,p);
	if (pkdIsDark(pkd,p)) {
	    ioHDF5AddDark(io,p->iOrder,p->r,v,
			  fMass,fSoft,*pPot );
	    }
	else if (pkdIsGas(pkd,p)) {
	    assert(0);
	    /* Why are temp and metals always set to zero? */
	    ioHDF5AddGas( io,p->iOrder,p->r,v,
			  fMass,fSoft,*pPot,0.0,0.0);
	    }
	else if (pkdIsStar(pkd,p)) {
	    assert(0);
	    /* Why are metals and tform always set to zero? */
	    ioHDF5AddStar(io, p->iOrder, p->r, v,
			  fMass,fSoft,*pPot,0.0,0.0);
	    }
	else mdlassert(pkd->mdl,0);

	ioHDF5AddVector( ioDen, p->iOrder, p->fDensity );
	ioHDF5AddVector( ioPot, p->iOrder, *pPot );
	}
    }


#endif

uint32_t pkdWriteTipsy(PKD pkd,char *pszFileName,uint64_t iFirst,
		       int bStandard,double dvFac,int bDoublePos) {
    FIO fio;
    PARTICLE *p;
    STARFIELDS *pStar;
    SPHFIELDS *pSph;
    float *pPot, dummypot;
    double v[3];
    float fMass, fSoft;
    uint32_t nCount;
    uint64_t iOrder;
    int i;

    fio = fioTipsyAppend(pszFileName,bDoublePos,bStandard);
    if (fio==NULL) {
	fprintf(stderr,"ERROR: unable to reopen tipsy file for output\n");
	perror(pszFileName);
	mdlassert(pkd->mdl,fio!=NULL);
	}
    fioSeek(fio,iFirst,FIO_SPECIES_ALL);

    v[0] = v[1] = v[2] = 0.0;
    dummypot = 0.0;
    nCount = 0;
    /*printf("Writing %lu\n", pkdLocal(pkd));*/
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;  /* JW: Ack! */
	if ( pkd->oPotential) pPot = pkdPot(pkd,p);
	else pPot = &dummypot;
	if (pkd->oVelocity) {
	    double *pV = pkdVel(pkd,p);
	    v[0] = pV[0] * dvFac;
	    v[1] = pV[1] * dvFac;
	    v[2] = pV[2] * dvFac;
	    }
	/* Initialize SPH fields if present */
	if (pkd->oSph) pSph = pkdField(p,pkd->oSph);
	else pSph = NULL;
	if (pkd->oStar) pStar = pkdField(p,pkd->oStar);
	else pStar = NULL;
	fMass = pkdMass(pkd,p);
	fSoft = pkdSoft0(pkd,p);
	iOrder = p->iOrder;
	nCount++;
	switch(pkdSpecies(pkd,p)) {
	case FIO_SPECIES_SPH:
	    assert(pSph);
	    assert(pkd->param.dTuFac>0.0);
		{
		double T, E;
		COOLPARTICLE cp;
		if (pkd->param.bGasCooling) {
		    E = pSph->u;
		    CoolTempFromEnergyCode( pkd->Cool, 
					    &cp, &E, &T, p->fDensity, pSph->fMetals );
		    }
		else T = pSph->u/pkd->param.dTuFac;
		fioWriteSph(fio,iOrder,p->r,v,fMass,fSoft,*pPot,
			    p->fDensity,T,pSph->fMetals);
		}
	    break;
	case FIO_SPECIES_DARK:
	    fioWriteDark(fio,iOrder,p->r,v,fMass,fSoft,*pPot);
	    break;
	case FIO_SPECIES_STAR:
	    assert(pStar && pSph);
	    fioWriteStar(fio,iOrder,p->r,v,fMass,fSoft,*pPot,
		       pSph->fMetals,pStar->fTimer);
	    break;
	default:
	    fprintf(stderr,"Unsupported particle type: %d\n",pkdSpecies(pkd,p));
	    assert(0);
	    }

	}
    fioClose(fio);
    return nCount;
    }


void pkdSetSoft(PKD pkd,double dSoft) {
    pkd->fSoftFix = dSoft;
    }

void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul) {
    pkd->fSoftFac = dFac;
    pkd->fSoftMax = bSoftMaxMul ? dSoftMax : HUGE;
    }

#ifdef USE_BSC_trace
static int foo = 0;
#endif

void
pkdGravAll(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bPeriodic,
	   int iOrder,int bEwald,double fEwCut,double fEwhCut,int *nActive,
	   double *pdPartSum, double *pdCellSum,CASTAT *pcs, double *pdFlop) {
    int bVeryActive = 0;

    pkdClearTimer(pkd,1);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    mdlTimeReset(pkd->mdl);
#endif

#ifdef USE_BSC_trace
    foo++;
    if ( foo == 1 ) {
	MPItrace_restart();
	}
#endif
#ifdef USE_BSC
    MPItrace_event(10000,3);
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
    *pdPartSum = 0.0;
    *pdCellSum = 0.0;
    pkdStartTimer(pkd,1);
    *nActive = pkdGravWalk(pkd,uRungLo,uRungHi,dTime,nReps,bPeriodic && bEwald,bVeryActive,pdFlop,pdPartSum,pdCellSum);
    pkdStopTimer(pkd,1);

#ifdef USE_BSC
    MPItrace_event(10001, *nActive);
    MPItrace_event(10000, 0 );
#endif
#ifdef USE_BSC_trace
    if ( foo == 1 ) {
	MPItrace_shutdown();
	}
#endif

    /*
    ** Get caching statistics.
    */
    pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
    pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
    pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
    pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
    pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
    pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
    pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
    pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
    /*
    ** Stop particle caching space.
    */
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    }


void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double L[]) {
    /* L is calculated with respect to the origin (0,0,0) */

    PARTICLE *p;
    float *pPot;
    double *v;
    FLOAT rx,ry,rz,vx,vy,vz,fMass;
    int i,n;

    n = pkdLocal(pkd);
    *T = 0.0;
    *U = 0.0;
    *Eth = 0.0;
    L[0] = L[1] = L[2] = 0;
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	fMass = pkdMass(pkd,p);
	if (pkd->oPotential) *U += 0.5*fMass*(*(pkdPot(pkd,p)));
	if (pkd->oSph && pkdIsGas(pkd,p)) *Eth += fMass*pkdSph(pkd,p)->u;
	if (pkd->oVelocity) {
	    v = pkdVel(pkd,p);
	    rx = p->r[0]; ry = p->r[1]; rz = p->r[2];
	    vx = v[0]; vy = v[1]; vz = v[2];
	    *T += 0.5*fMass*(vx*vx + vy*vy + vz*vz);
	    L[0] += fMass*(ry*vz - rz*vy);
	    L[1] += fMass*(rz*vx - rx*vz);
	    L[2] += fMass*(rx*vy - ry*vx);
	    }
	}
    }


/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd,double dTime,double dDelta,double dDeltaVPred,double dDeltaUPred,uint8_t uRungLo,uint8_t uRungHi) {
    PARTICLE *p;
    double *v;
    float *a;
    SPHFIELDS *sph;
    int i,j,n;
    double dMin[3],dMax[3];

    mdlDiag(pkd->mdl, "Into pkdDrift\n");
    assert(pkd->oVelocity);

#ifdef USE_BSC
    MPItrace_event(10000,4);
#endif
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
		    p->r[j] += dDelta*v[j];
		    }
		pkdMinMax(p->r,dMin,dMax);
		if (p->iOrder == 32769) printf("SFTEST %i: drift %20.14g   %g\n",p->iOrder,dTime,dDeltaUPred);
		}
	    }
	}
    else {
	for (i=0;i<n;++i) {
	    p = pkdParticle(pkd,i);
	    if (pkdIsRungRange(p,uRungLo,uRungHi)) {
		v = pkdVel(pkd,p);
		for (j=0;j<3;++j) {
		    p->r[j] += dDelta*v[j];
		    }
		pkdMinMax(p->r,dMin,dMax);
		}
	    }
	}

    for (j=0;j<3;++j) {
	pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
#ifdef USE_BSC
    MPItrace_event(10000,0);
#endif
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
    }


void pkdGravityVeryActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int bEwald,int nReps,double dStep) {
    int nActive;
    int bVeryActive = 1;
    double dFlop,dPartSum,dCellSum;

    /*
    ** Calculate newtonian gravity for the very active particles ONLY, including replicas if any.
    */
    dFlop = 0.0;
    dPartSum = 0.0;
    dCellSum = 0.0;
    nActive = pkdGravWalk(pkd,uRungLo,uRungHi,dTime,nReps,bEwald,bVeryActive,&dFlop,&dPartSum,&dCellSum);
    }


void pkdStepVeryActiveKDK(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dStep, double dTime, double dDelta,
			  int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
			  int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass) {
    int nRungCount[256];
    double dDriftFac;
#ifdef PLANETS
    double aSun[3], adSun[3];
    int j;
#endif

    if (iAdjust && (iRung < pkd->param.iMaxRung-1)) {

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
			 pkd->param.bEpsAccStep,pkd->param.bSqrtPhiStep,dhMinOverSoft);
	    }
	*pnMaxRung = pkdUpdateRung(pkd,iRung,pkd->param.iMaxRung-1,
				   iRung,pkd->param.iMaxRung-1, nRungCount);


	if (pkd->param.bVDetails) {
	    printf("%*cAdjust at iRung: %d, nMaxRung:%d nRungCount[%d]=%d\n",
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
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	dStep += 1.0/(2 << iRung);
	dTime += 0.5*dDelta;

	pkdActiveRung(pkd,iRung,0);   /* is this needed? */

	pkdStepVeryActiveKDK(pkd,uRungLo,uRungHi,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
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
	pkdDrift(pkd,dTime,dDriftFac,0,0,iRungVeryActive+1,MAX_RUNG);
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
	    pkdVATreeBuild(pkd,pkd->param.nBucket,diCrit2,dTime);
	    pkdGravityVeryActive(pkd,uRungLo,uRungHi,dTime,pkd->param.bEwald && pkd->param.bPeriodic,pkd->param.nReplicas,dStep);

#ifdef PLANETS
	    /* Sun's gravity */
	    if (pkd->param.bHeliocentric) {
		/* Sun's indirect gravity due to very active particles*/
		pkdActiveRung(pkd,iRungVeryActive+1,1);
		pkdSunIndirect(pkd,aSun,adSun,1);
		for (j=0;j<3;++j) {
		    aSun[j] += aSunInact[j];
		    adSun[j] += adSunInact[j];
		    }
		pkdActiveRung(pkd,iKickRung,1);
		pkdGravSun(pkd,aSun,adSun,dSunMass);
		}
#endif
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

#ifdef HERMITE
void
pkdStepVeryActiveHermite(PKD pkd, double dStep, double dTime, double dDelta,
			 int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
			 int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass) {
    int nRungCount[256];
    double dDriftFac;
    int i;

    double aSun[3], adSun[3];
    int j;

    if (iAdjust && (iRung < pkd->param.iMaxRung-1)) {
	pkdActiveRung(pkd, iRung, 1);
	if (pkd->param.bAarsethStep) {
	    pkdAarsethStep(pkd,pkd->param.dEta);
	    }
	if (pkd->param.bAccelStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dVelFac = 1.0/(a*a);
	    double dAccFac = 1.0/(a*a*a);
	    double dhMinOverSoft = 0;
	    pkdAccelStep(pkd,uRungLo,uRungHi,pkd->param.dEta, dVelFac,dAccFac,pkd->param.bDoGravity,
			 pkd->param.bEpsAccStep,pkd->param.bSqrtPhiStep,dhMinOverSoft);
	    }
	*pnMaxRung = pkdUpdateRung(pkd,iRung,pkd->param.iMaxRung-1,
				   iRung,pkd->param.iMaxRung-1, 1, nRungCount);

	if (pkd->param.bVDetails) {
	    printf("%*cAdjust at iRung: %d, nMaxRung:%d nRungCount[%d]=%d\n",
		   2*iRung+2,' ',iRung,*pnMaxRung,*pnMaxRung,nRungCount[*pnMaxRung]);
	    }

	}

    if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
	pkdStepVeryActiveHermite(pkd,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
				 diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	dStep += 1.0/(2 << iRung);
	dTime += 0.5*dDelta;
	pkdActiveRung(pkd,iRung,0);
	pkdStepVeryActiveHermite(pkd,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
				 diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	}
    else {
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive Predictor at iRung: %d, drifting %d and higher with dDelta: %g\n",
		   2*iRung+2,' ',iRung,iRungVeryActive+1,dDelta);
	    }
	/*
	** This should predict *all* very actives!
	*/
	pkdActiveRung(pkd,iRungVeryActive+1,1);
	/*
	** We need to account for cosmological drift factor here!
	** Normally this is done at the MASTER level in msrDrift.
	** Note that for kicks we have written new "master-like" functions
	** KickOpen and KickClose which do this same job at PKD level.
	*/
	/*if (pkd->param.csm->bComove) {
	    dDriftFac = csmComoveDriftFac(pkd->param.csm,dTime,dDelta);
	    }
	else {
	    dDriftFac = dDelta;
	    }*/

	dTime += dDelta;
	dStep += 1.0/(1 << iRung);

	pkdPredictor(pkd,dTime);

	if (iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master().
						*/

	    if (pkd->param.bVDetails) {
		printf("%*cGravityVA: iRung %d Gravity for rungs %d to %d ...\n ",
		       2*iRung+2,' ',iRung,iKickRung,*pnMaxRung);
		}

	    pkdActiveRung(pkd,iKickRung,1);
	    pkdVATreeBuild(pkd,pkd->param.nBucket,diCrit2,dTime);
	    pkdGravityVeryActive(pkd,dTime,pkd->param.bEwald && pkd->param.bPeriodic,pkd->param.nReplicas,dStep);


#ifdef PLANETS
	    /* Sun's gravity */
	    if (pkd->param.bHeliocentric) {
		/* Sun's indirect gravity due to very active particles*/
		pkdActiveRung(pkd,iRungVeryActive+1,1);
		pkdSunIndirect(pkd,aSun,adSun,1);
		for (j=0;j<3;++j) {
		    aSun[j] += aSunInact[j];
		    adSun[j] += adSunInact[j];
		    }
		pkdActiveRung(pkd,iKickRung,1);
		pkdGravSun(pkd,aSun,adSun,dSunMass);
		}
#endif

	    if (pkd->param.bVDetails) {
		printf("%*cVeryActive pkdCorrector at iRung: %d, 0.5*dDelta: %g\n",
		       2*iRung+2,' ',iRung,0.5*dDelta);
		}
	    pkdCorrector(pkd,dTime);

#ifdef PLANETS
	    if (pkd->param.bHeliocentric) {
		int nite = 0;
		int nitemax = 3;
		do {
		    nite += 1;
		    pkdSunCorrector(pkd,dTime,dSunMass);
		    }
		while (nite < nitemax);
		}
	    if (pkd->param.bCollision && pkd->iCollisionflag) pkdDoCollisionVeryActive(pkd,dTime);
#endif
	    pkdCopy0(pkd,dTime);
	    }
	}
    }

void
pkdCopy0(PKD pkd,double dTime) {
    PARTICLE *p;
    int i,j,n;

    mdlDiag(pkd->mdl, "Into pkdCopy0\n");

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    for (j=0;j<3;++j) {
		p->r0[j] = p->r[j];
		p->v0[j] = p->v[j];
		p->a0[j] = p->a[j];
		p->ad0[j] = p->ad[j];
		}
	    p->dTime0 = dTime;
#ifdef PLANETS
	    if (pkd->param.bCollision) {
		p->iColflag = 0;  /* just in case */
		}
#endif
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdCopy0\n");
    }

void
pkdPredictor(PKD pkd,double dTime) {
    PARTICLE *p;
    int i,j,n;
    double dt; /* time since last evaluation of force */

    mdlDiag(pkd->mdl, "Into pkdPredictor\n");

    if (pkd->param.bVDetails) printf("Into pkdPredictor\n");

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    dt =  dTime - p->dTime0;
	    for (j=0;j<3;++j) {
		p->r[j] = p->r0[j] + dt*(p->v0[j] + 0.5*dt*(p->a0[j]+dt*p->ad0[j]/3.0));
		p->v[j] = p->v0[j] + dt*(p->a0[j]+0.5*dt*p->ad0[j]);
		p->rp[j] = p->r[j];
		p->vp[j] = p->v[j];
		}
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdPredictor\n");
    }

void pkdCorrector(PKD pkd,double dTime) {
    PARTICLE *p;
    int i,j,n;
    double dt, add, addd, am;
    double a0d2, a1d2, a2d2, a3d2, dti, dt2i;
    /*double alpha = 7.0/6.0;*/

    mdlDiag(pkd->mdl, "Into pkdCorrector\n");

    if (pkd->param.bVDetails) printf("Into pkdCorrector\n");

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    dt =  dTime - p->dTime0;
	    if (pkd->param.bAarsethStep) {
		a0d2  = 0.0;
		a1d2  = 0.0;
		a2d2  = 0.0;
		a3d2  = 0.0;
		dti = 1.0/dt;
		dt2i = dti*dti;
		}

	    for (j=0;j<3;++j) {
		am = p->a0[j]-p->a[j];

		if (pkd->param.bAarsethStep) {
		    add = 2.0*(3.0*am*dt2i+(2.0*p->ad0[j]+p->ad[j])*dti);
		    addd = 6.0*(2.0*am*dti+(p->ad0[j]+p->ad[j]))*dt2i;
		    add += dt*addd;

		    a0d2 += p->a[j]*p->a[j];
		    a1d2 += p->ad[j]*p->ad[j];
		    a2d2 += add*add;
		    a3d2 += addd*addd;
		    }

		add = 16.0*am+(13.0*p->ad0[j]+3.0*p->ad[j])*dt;
		addd = 6.0*am+(5.0*p->ad0[j]+p->ad[j])*dt;
		p->r[j] = p->rp[j] - add/120.0*dt*dt;
		p->v[j] = p->vp[j] - addd/12.0*dt;

		}
	    if (pkd->param.bAarsethStep) {
		p->dtGrav = (sqrt(a0d2*a2d2)+a1d2)/(sqrt(a1d2*a3d2)+a2d2);
		}
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdCorrector\n");
    }

void
pkdSunCorrector(PKD pkd,double dTime,double dSunMass) {
    /* iteratively correct only sun's direct gravity */
    PARTICLE *p;
    int i,j,n;
    double dt;
    double add, addd, r2,r3i,r5i,rv, am;
    /*double alpha = 7.0/6.0;*/

    mdlDiag(pkd->mdl, "Into pkdSunCorrector\n");

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    dt =  dTime - p->dTime0;
	    r2  = 0.0;
	    rv  = 0.0;
	    for (j=0;j<3;++j) {
		r2 += p->r[j]*p->r[j];
		rv += p->r[j]*p->v[j];
		}
	    r2 = 1.0/r2;
	    r3i = r2*sqrt(r2)*dSunMass;
	    r5i = 3.0*r3i*r2*rv;

	    for (j=0;j<3;++j) {
		p->a[j] = -p->r[j]*r3i + p->app[j];
		p->ad[j] = -(p->v[j]*r3i - p->r[j]*r5i) + p->adpp[j];
		am = p->a0[j]-p->a[j];
		add = 16.0*am+(13.0*p->ad0[j]+3.0*p->ad[j])*dt;
		addd = 6.0*am+(5.0*p->ad0[j]+p->ad[j])*dt;
		p->r[j] = p->rp[j] - add/120.0*dt*dt;
		p->v[j] = p->vp[j] - addd/12.0*dt;

		}
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdSunCorrector\n");
    }

void
pkdPredictorInactive(PKD pkd,double dTime) {
    PARTICLE *p;
    int i,j,n;
    double dt; /* time since last evaluation of force */

    mdlDiag(pkd->mdl, "Into pkdPredictorInactive\n");

    if (pkd->param.bVDetails) printf("Into pkdPredictorInactive\n");

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) continue;
	dt =  dTime - p->dTime0;
	for (j=0;j<3;++j) {
	    p->r[j] = p->r0[j] + dt*(p->v0[j] + 0.5*dt*(p->a0[j]+dt*p->ad0[j]/3.0));
	    p->v[j] = p->v0[j] + dt*(p->a0[j]+0.5*dt*p->ad0[j]);
	    p->rp[j] = p->r[j];
	    p->vp[j] = p->v[j];
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdPredictorInactive\n");
    }

void pkdAarsethStep(PKD pkd,double dEta) {
    PARTICLE *p;
    double dT;
    int i;

    for (i=0;i<pkdLocal(pkd);i++) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    mdlassert(pkd->mdl, p->dtGrav > 0);
	    dT = dEta*sqrt(p->dtGrav);
	    if (dT < p->dt)
		p->dt = dT;
	    }
	}
    }

void pkdFirstDt(PKD pkd) {
    int i,j;
    PARTICLE *p;
    double a0d2,a1d2;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	a0d2 = 0.0;
	a1d2 = 0.0;
	if (pkdIsActive(pkd,p))
	    for (j=0;j<3;++j) {
		a0d2 += p->a0[j]*p->a0[j];
		a1d2 += p->ad0[j]*p->ad0[j];
		}
	p->dtGrav = (a0d2/a1d2);
#ifdef PLANETS
	if (pkd->param.bCollision) {
	    p->iColflag = 0; /* initial reset of collision flag */
	    }
#endif
	}
    }
#endif /* Hermite */

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
    double *v;
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
		if (p->iOrder == 32769) printf("SFTEST %i: kick %g   %g %g\n",p->iOrder,dTime,dDeltaU,dDeltaU*2);
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
		  int bDoGravity,int bEpsAcc,int bSqrtPhi,double dhMinOverSoft) {
    PARTICLE *p;
    float *a, *pPot;
    double *v;
    int i,uNewRung;
    double vel;
    double acc;
    int j;
    double dT;
    FLOAT fSoft;

    assert(pkd->oVelocity);
    assert(pkd->oAcceleration);
    assert(pkd->oPotential);

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
	    mdlassert(pkd->mdl,acc > 0);
	    dT = FLOAT_MAXVAL;
	    if (bEpsAcc && acc>0) {
		dT = dEta*sqrt(fSoft/acc);
		}
	    if (bSqrtPhi) {
		/*
		** NOTE: The factor of 3.5 keeps this criterion in sync
		** with DensityStep. The nominal value of dEta for both
		** cases is then 0.02-0.03.
		*/
		pPot = pkdPot(pkd,p);
		double dtemp =
		    dEta*3.5*sqrt(dAccFac*fabs(*pPot))/acc;
		if (dtemp < dT)
		    dT = dtemp;
		}
	    uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
	    if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
	    }
	}
    }


void pkdSphStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,double dAccFac) {
    PARTICLE *p;
    float *a, uDot;
    int i,j,uNewRung;
    double acc;
    double dtNew,dtC;
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
		if (acc>0) dtNew = pkd->param.dEta*sqrt(p->fBall/acc);
		u2 = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung-1);
		uDot = *pkd_uDot(pkd,p);
		u3=0;
		if (uDot < 0) {
		    double dtemp = pkd->param.dEtaUDot*(*pkd_u(pkd,p))/fabs(uDot);
		    if (dtemp < dtNew) dtNew = dtemp;
		    u3 = pkdDtToRung(dtemp,pkd->param.dDelta,pkd->param.iMaxRung-1);
		    }
/*	    if (p->uRung > 5) printf("%d %d: %g %g %g *%g* %g %g,\n",i,p->iOrder,*pkd_u(pkd,p),*pkd_uDot(pkd,p),*pkd_divv(pkd,p),1/(*pkd_divv(pkd,p)),pkd->param.dDelta,dtNew); */

		uNewRung = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung-1);
		if (uNewRung > p->uNewRung) p->uNewRung = uNewRung;
		if (!(p->iOrder%10000) || (p->uNewRung > 5 && !(p->iOrder%1000))) {
		    SPHFIELDS *sph = pkdSph(pkd,p);
		    double T, E = sph->u;
		    COOLPARTICLE cp;
		    if (pkd->param.bGasIsothermal) T = E/pkd->param.dTuFac;
		    else CoolTempFromEnergyCode( pkd->Cool, &cp, &E, &T, p->fDensity, sph->fMetals );
//		    printf("RUNG %d: grav+sph %d acc %d udot %d final %d (prev %d)\nRUNG %d: dens %16.10g temp %g r %g\n",p->iOrder,u1,u2,u3,(int) p->uNewRung,p->uRung, p->iOrder, p->fDensity, T, sqrt(p->r[0]*p->r[0]+p->r[1]*p->r[1]));
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
    COOLPARTICLE cp;
    SPHFIELDS *sph;
    int n = pkdLocal(pkd);
    double T, E, dmstar, dt, prob;
    PARTICLE *starp;
    int i,j;
    
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
	
	if (p->iOrder==6667) {
	    dt = pkd->param.dDelta/(1<<p->uRung); /* Actual Rung */
	    printf("SF %d: %d %d %d %g\n",p->iOrder,pkdIsActive(pkd,p),pkdIsGas(pkd,p),p->uRung,dt,pkdStar(pkd,p)->totaltime);
	    }
	if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
	    sph = pkdSph(pkd,p);
	    dt = pkd->param.dDelta/(1<<p->uRung); /* Actual Rung */
	    pkdStar(pkd,p)->totaltime += dt;
	    if (p->iOrder==6667) {
		printf("SF %d: %d %d %d %g AFTER PLUS\n",p->iOrder,pkdIsActive(pkd,p),pkdIsGas(pkd,p),p->uRung,dt,pkdStar(pkd,p)->totaltime);
		}


	    if (p->fDensity < dDenMin || (bdivv && sph->divv >= 0.0)) continue;
	    E = sph->uPred;
	    if (pkd->param.bGasCooling) 
		CoolTempFromEnergyCode( pkd->Cool, &cp, &E, &T, p->fDensity, sph->fMetals );
	    else T=E/pkd->param.dTuFac;
	    if (T > dTMax) continue;
	    
            /* Note: Ramses allows for multiple stars per step -- but we have many particles
	      and he has one cell that may contain many times m_particle */
	    if (pkd->param.bGasCooling) {
		if (fabs(pkdStar(pkd,p)->totaltime-dTime) > 1e-3*dt) {
		    printf("total time error: %i,  %g %g %g\n",p->iOrder,pkdStar(pkd,p)->totaltime,dTime,dt);
		    assert(0);
		    }
		}

	    dmstar = dRateCoeff*sqrt(p->fDensity)*pkdMass(pkd,p)*dt;
	    prob = 1.0 - exp(-dmstar/dInitStarMass); 
//	    if (!(p->iOrder%1000)) printf("SF %d: %g %g %g  %g\n",p->iOrder,dt,dmstar,dInitStarMass,prob);
	    
	    /* Star formation event? */
	    if (rand()<RAND_MAX*prob) {
		float *starpMass = pkdField(starp,pkd->oMass);
		float *pMass = pkdField(p,pkd->oMass);
		pkdCopyParticle(pkd, starp, p);	/* grab copy */
		*pMass -= dInitStarMass;
		*starpMass = dInitStarMass;
//		pkdStar(pkd,starp)->iGasOrder = p->iOrder;
		if (*pMass < 0) {
		    *starpMass += *pMass;
		    *pMass = 0;
		    }
	        if (*pMass < dMinGasMass) {
//		    printf("Delete* particle %d %d %d %d, %g\n",(int) p->iOrder,(int) pkdIsGas(pkd,p),pkdSpecies( pkd, p ), FIO_SPECIES_SPH,*pMass);
		    pkdDeleteParticle(pkd, p);
//		    printf("Deleted particle %d %d %d %d, %g\n",(int) p->iOrder,(int) pkdIsGas(pkd,p),pkdSpecies( pkd, p ),FIO_SPECIES_LAST,*pMass);
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
    int i,n;
    SPHFIELDS *sph;
    COOLPARTICLE cp;  /* Dummy: Not yet fully implemented */
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

	CoolSetTime( pkd->Cool, dTime, z, bUpdateTable );
	
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
			CoolIntegrateEnergyCode(pkd->Cool, &cp, &E, ExternalHeating, p->fDensity, sph->fMetals, p->r, dt);
			uDot = (E-sph->u)/dt; 
			if (uDot < 0) {
			    double dtNew;
			    int uNewRung;
			    dtNew = pkd->param.dEtaUDot*sph->u/fabs(uDot);
			    uNewRung = pkdDtToRung(dtNew,pkd->param.dDelta,pkd->param.iMaxRung-1);
			    if (uNewRung > p->uNewRung) {
				p->uNewRung = uNewRung;
				continue;
				}
			    }
			sph->uDot = uDot;
			break;
			}
/*	    printf("%d %d: %g %g %g *%g* %g %g,\n",i,p->iOrder,*pkd_u(pkd,p),*pkd_uDot(pkd,p),*pkd_divv(pkd,p),*pkd_divv(pkd,p)*pkd->param.dDelta,pkd->param.dDelta,dT); */
		    }
		}
	    }
	else {
	    for (i=0;i<pkdLocal(pkd);++i) {
		p = pkdParticle(pkd,i);
		if (pkdIsActive(pkd,p) && pkdIsGas(pkd,p)) {
		    if (pkdStar(pkd,p)->fTimer > dTime) {
//		    printf("COOLING shut off %d: %g %g %g  %g %g\n",p->iOrder,pkdStar(pkd,p)->fTimer,dTime,(dTime-pkdStar(pkd,p)->fTimer)*1.7861e+18/(365.*60*60*24.)/1e6,pkdSph(pkd,p)->u,pkdSph(pkd,p)->uPred);
			continue;
			}
		    sph = pkdSph(pkd,p);
		    ExternalHeating = sph->uDot;
		    E = sph->u;
		    dt = pkd->param.dDelta/(1<<p->uRung); /* Actual Rung */
		    CoolIntegrateEnergyCode(pkd->Cool, &cp, &E, ExternalHeating, p->fDensity, sph->fMetals, p->r, dt);
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
    COOL *cl;
    double T,E;
    COOLPARTICLE cp; /* Dummy for now */

    cl = pkd->Cool;
    CoolSetTime( cl, dTime, z, 1 );

    switch(iDirection)  {
    case CORRECTENERGY_IN:
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
	break;
	/* Careful using this -- it permanenty converts the thermal energy */
    case CORRECTENERGY_OUT: 
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
	break;
    case CORRECTENERGY_SPECIAL:
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
	    dT = dEta/sqrt(p->fDensity*dRhoFac);
	    p->uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
	    }
	}
    }

uint8_t pkdDtToRung(double dT, double dDelta, uint8_t uMaxRung) {
    int iSteps;
    uint8_t uRung;

    assert(dT>0.0);
    iSteps = dDelta/dT;
    uRung = 0;
    while (iSteps) {
	++uRung;
	iSteps >>= 1;
	}
    if ( uRung > uMaxRung ) uRung = uMaxRung;
    return uRung;
    }

int pkdUpdateRung(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,
		  uint8_t uRung,int iMaxRung,int *nRungCount) {
    PARTICLE *p;
    int i;
    int iTempRung;
    for (i=0;i<iMaxRung;++i) nRungCount[i] = 0;
	
//    printf("RUNG UPDATES %d %d %d\n",uRungLo,uRungHi,uRung);
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if ( pkdIsActive(pkd,p) ) {
	    if ( p->uNewRung >= iMaxRung )
		p->uRung = iMaxRung-1;
	    else if ( p->uNewRung >= uRung )
		p->uRung = p->uNewRung;
	    else if ( p->uRung > uRung)
		p->uRung = uRung;

//	    if (!(p->iOrder%10000) || (p->uRung > 5 && !(p->iOrder%1000))) printf("RUNG %d: UPDATED %d %d\n",p->iOrder,p->uNewRung,p->uRung);
	    }
	/*
	** Now produce a count of particles in rungs.
	*/
	nRungCount[p->uRung] += 1;
	}
    iTempRung = iMaxRung-1;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
    }

void pkdDeleteParticle(PKD pkd, PARTICLE *p) {

#ifdef PLANETS
    /* the following operation is an emergent treatment for
       to-be-deleted VeryActive particles. */

    int j;
    assert(0);
    if (pkd->param.bGravStep) {
	p->dtGrav = 0.00001;
	}
#ifdef HERMITE
    if (pkd->param.bAarsethStep) {
	p->dtGrav = 10000;
	}
#endif
#ifdef SYMBA
    p->drmin = 1000.0;
    p->n_VA=0;
#endif
    for (j=0;j<3;j++) {
	p->r[j] = 100.0*p->iOrder;
	p->v[j] = 0.0;
#ifdef HERMITE
	p->r0[j] = 100.0*p->iOrder;
	p->v0[j] = 0.0;
	p->a0[j] =  0.000001;
	p->ad0[j] = 0.000001;
	p->a[j] =  0.000001;
	p->ad[j] = 0.000001;
#endif
	}
#endif

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
//	    printf("DELETE %d: %d %d  %g %g\n",p->iOrder,pkdIsGas(pkd,p),pkdIsStar(pkd,p),pkdMass(pkd,p),pkdStar(pkd,p)->fTimer);
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

//    printf("NEWORDER: New particles from iOrder=%d\n",nStart);
    for (pi=0;pi<pkdLocal(pkd);pi++) {
	p = pkdParticle(pkd,pi);
	if (p->iOrder == IORDERMAX) {
//	    printf("NEWORDER: Assigning new particle %d\n",nStart);
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
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    int iOrder;
    PARTICLE *p;
    
    n = 0;
    nGas = 0;
    nDark = 0;
    nStar = 0;
    iMaxOrderGas = -1;
    iMaxOrderDark = -1;
    iMaxOrderStar = -1;
    for(pi = 0; pi < pkdLocal(pkd); pi++) {
	p = pkdParticle(pkd,pi);
	iOrder = p->iOrder;
	n++;
	if(pkdIsGas(pkd, p)) {
	    ++nGas;
	    if (iOrder > iMaxOrderGas) iMaxOrderGas = iOrder;
	    }
	else if(pkdIsDark(pkd, p)) {
	    ++nDark;
	    if (iOrder > iMaxOrderDark) iMaxOrderDark = iOrder;
	    }
	else if(pkdIsStar(pkd, p)) {
	    ++nStar;
	    if (iOrder > iMaxOrderStar) iMaxOrderStar = iOrder;
	    }
	}
    
    out->n  = n;
    out->nGas = nGas;
    out->nDark = nDark;
    out->nStar = nStar;
    out->iMaxOrderGas = iMaxOrderGas;
    out->iMaxOrderDark = iMaxOrderDark;
    out->iMaxOrderStar = iMaxOrderStar;
}


void pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar,int nMaxOrderGas,
		  int nMaxOrderDark, int nMaxOrder) {
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
/* JW: Depr.
    pkd->nMaxOrder = nMaxOrder;
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;*/
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

#ifdef PLANETS
void
pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal) {
    SSIO ssio;
    SSDATA data;
    PARTICLE *p;
    int i,j;

    pkd->nLocal = nLocal;
    pkd->nActive = nLocal;
    /*
     ** General initialization (modeled after pkdReadTipsy()).
     */
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	p->iRung = 0;
	p->fDensity = 0.0;
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	for (j=0;j<3,++j) {
	    p->a[j] = 0.0;
	    }
	}
    /*
     ** Seek past the header and up to nStart.
     */
    if (ssioOpen(pszFileName,&ssio,SSIO_READ))
	mdlassert(pkd->mdl,0); /* unable to open ss file */
    if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
	mdlassert(pkd->mdl,0); /* unable to seek in ss file */
    /*
     ** Read Stuff!
     */
    for (i=0;i<nLocal;++i) {
	p = pkdParticle(pkd,i);
	p->iOrder = nStart + i;
	if (ssioData(&ssio,&data))
	    mdlassert(pkd->mdl,0); /* error during read in ss file */
	p->iOrgIdx = data.org_idx;
	getClass(data.mass,data.radius,p);
	for (j=0;j<3;++j) p->r[j] = data.pos[j];
	for (j=0;j<3;++j) p->v[j] = data.vel[j];
	for (j=0;j<3;++j) p->w[j] = data.spin[j];
	p->iColor = data.color;
#ifdef NEED_VPRED
	for (j=0;j<3;++j) p->vPred[j] = p->v[j];
#endif
	}
    if (ssioClose(&ssio))
	mdlassert(pkd->mdl,0); /* unable to close ss file */
    }

void
pkdWriteSS(PKD pkd,char *pszFileName,int nStart) {
    SSIO ssio;
    SSDATA data;
    PARTICLE *p;
    int i,j;

    /*
     ** Seek past the header and up to nStart.
     */
    if (ssioOpen(pszFileName,&ssio,SSIO_UPDATE))
	mdlassert(pkd->mdl,0); /* unable to open ss file */
    if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
	mdlassert(pkd->mdl,0); /* unable to seek in ss file */
    /*
     ** Write Stuff!
     */
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (!pkdIsDark(pkd,p))
	    mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
	data.org_idx = p->iOrgIdx;
	data.mass = pkdMass(pkd,p);
	data.radius = pkdSoft(pkd,p);
	for (j=0;j<3;++j) data.pos[j]  = p->r[j];
	for (j=0;j<3;++j) data.vel[j]  = p->v[j];
	for (j=0;j<3;++j) data.spin[j] = p->w[j];
	data.color = p->iColor;
	if (ssioData(&ssio,&data))
	    mdlassert(pkd->mdl,0); /* unable to write in ss file */
	}
    if (ssioClose(&ssio))
	mdlassert(pkd->mdl,0); /* unable to close ss file */
    }

void pkdSunIndirect(PKD pkd,double aSun[],double adSun[],int iFlag) {
    PARTICLE *p;
    FLOAT r2,r1i,r3i,r5i,rv,fMass;
    int i,j,n;

    for (j=0;j<3;++j) {
	aSun[j] = 0;
	adSun[j] = 0;
	}
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	fMass = pkdMass(pkd,p);
	if (iFlag == 2) {
	    if (pkdIsActive(pkd,p)) continue;  /* inactive */
	    }
	else if (iFlag == 1) {
	    if (!pkdIsActive(pkd,p)) continue; /* active */
	    }

	r2 = 0;
	rv = 0;
	for (j=0;j<3;++j) {
	    r2 += p->r[j]*p->r[j];
	    rv += p->v[j]*p->r[j];
	    }
	r1i = (r2 == 0 ? 0 : 1/sqrt(r2));
	r3i = fMass*r1i*r1i*r1i;
	r5i = 3.0*rv*r3i*r1i*r1i;
	for (j=0;j<3;++j) {
	    aSun[j] += p->r[j]*r3i;
	    adSun[j] += p->v[j]*r3i-p->r[j]*r5i;
	    }
	}
    }

void pkdGravSun(PKD pkd,double aSun[],double adSun[],double dSunMass) {
    PARTICLE *p;
    double r2,v2,r1i,r3i;
    double aai,aa3i,idt2;
    double rv, r5i;
    int i,j,n;
    double hx,hy,hz,h2,E,e,sum;

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    r2 = 0;
	    v2 = 0;
	    rv = 0;
	    for (j=0;j<3;++j) {
		r2 += p->r[j]*p->r[j];
		v2 += p->v[j]*p->v[j];
		rv += p->v[j]*p->r[j];
		}

	    r1i = (r2 == 0 ? 0 : 1/sqrt(r2)); /*gravity at origin = zero */
	    p->fPot -= dSunMass*r1i;
	    r3i = dSunMass*r1i*r1i*r1i;
	    r5i = 3.0*rv*r3i*r1i*r1i;
	    /* time step is determined by semimajor axis, not the heliocentric distance*/
	    if (pkd->param.bGravStep && p->fMass > 0) {
		/* E and h are normalized by the reduced mass */
		sum = dSunMass + p->fMass;
		hx =  p->r[1]*p->v[2] - p->r[2]*p->v[1];
		hy =  p->r[2]*p->v[0] - p->r[0]*p->v[2];
		hz =  p->r[0]*p->v[1] - p->r[1]*p->v[0];
		h2 = hx*hx + hy*hy + hz*hz;
		E = 0.5*v2 - sum*r1i;
		e = sqrt(1.0+2.0*E*h2/sum/sum);
		aai = -2.0*E/sum/(1.0-e);
		aai = aai*aai*aai;
		idt2 = sum*aai;
		/*if (p->dtSun > p->dtGrav) p->dtGrav = p->dtSun;*/
		if (idt2 > p->dtGrav) p->dtGrav = idt2;
		}
	    for (j=0;j<3;++j) {
#ifdef HERMITE
		if (pkd->param.bHermite) {
		    p->app[j] = p->a[j] - aSun[j]; /* perturbation force*/
		    p->adpp[j] = p->ad[j] - adSun[j];
		    p->ad[j] -= (adSun[j] + p->v[j]*r3i-p->r[j]*r5i);
		    }
#endif
		p->a[j] -= (aSun[j] + p->r[j]*r3i);
		}
	    }
	}
    }

void pkdHandSunMass(PKD pkd, double dSunMass) {
    pkd->dSunMass = dSunMass;
    }

#ifdef SYMBA
void
pkdStepVeryActiveSymba(PKD pkd, double dStep, double dTime, double dDelta,
		       int iRung, int iKickRung, int iRungVeryActive,
		       int iAdjust, double diCrit2,
		       int *pnMaxRung, double dSunMass, int multiflag) {
    int nRungCount[256];
    double dDriftFac;
    double ddDelta, ddStep;


    if (iRung ==0) {
	/* get pointa of interacting opponents from thier iOrder's
	   multiflag > 0 if close encounters with more than 2 particles exist*/

	pkd->nVeryActive = pkdSortVA(pkd, iRung);
	multiflag = pkdGetPointa(pkd);
	/*printf("Time %e, nVA %d,  multiflag %d \n",dTime, pkd->nVeryActive, multiflag);*/
	}
    if (iRung > 0) {
	/* at iRung = 0 we just call pkdStepVeryActiveSymba three times*/

	pkdActiveRung(pkd,iRung,1);
	if (iAdjust && (iRung < pkd->param.iMaxRung-1)) {
	    /* determine KickRung from current position */
	    *pnMaxRung = pkdDrminToRungVA(pkd,iRung,pkd->param.iMaxRung,multiflag);
	    if (pkd->param.bVDetails) {
		printf("%*cAdjust, iRung: %d\n",2*iRung+2,' ',iRung);
		}
	    }
	if (iAdjust == 0) {
	    /*if(iRung >= 2) pkdSortVA(pkd, iRung); maybe added later */
	    pkdGravVA(pkd,iRung);
	    if (pkd->param.bVDetails) {
		printf("%*cpkdGravVA, iRung: %d\n",2*iRung+2,' ',iRung);
		}

	    }
	/* following kick is for particles at iRung and iRung + 1 */
	pkdKickVA(pkd,0.5*dDelta);

	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickOpen  at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }

	/* if particles exist at the current iRung */
	pkdActiveRung(pkd,iRung,0);
	pkdKeplerDrift(pkd,dDelta,dSunMass,1); /* 1 means VA */
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKeplerDrift  at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	/* if drmin druing drift is less than R_(iRung+1),
	   this interacting pair is sent to iRung + 1*/
	if (iRung < pkd->param.iMaxRung-1) {
	    *pnMaxRung = pkdCheckDrminVA(pkd,iRung,multiflag,*pnMaxRung);
	    }
	}

    if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
	ddDelta = dDelta/3.0;
	ddStep = 1.0/pow(3.0, iRung); /* this is currently unnecessary */

	pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iRung+1,iRungVeryActive,0,
			       diCrit2,pnMaxRung,dSunMass,multiflag);
	dStep += ddStep;
	dTime += ddDelta;
	pkdActiveRung(pkd,iRung,0);
	pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iRung+1,iRungVeryActive,1,
			       diCrit2,pnMaxRung,dSunMass,multiflag);
	dStep += ddStep;
	dTime += ddDelta;
	pkdActiveRung(pkd,iRung,0);
	pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iKickRung,iRungVeryActive,1,
			       diCrit2,pnMaxRung,dSunMass,multiflag);

	/* move time back to 1 step */
	dTime -= dDelta;
	}
    if (iRung > 0) {
	dTime += 0.5*dDelta;

	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickClose at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }

	/* Kick due to F_iRung for particles at the current or higher iRungs */
	pkdActiveRung(pkd,iRung,1);
	pkdGravVA(pkd,iRung);
	pkdKickVA(pkd,0.5*dDelta);

	if (pkd->param.bCollision && pkd->iCollisionflag) pkdDoCollisionVeryActive(pkd,dTime);

	}
    }

int pkdSortVA(PKD pkd, int iRung) {
    int i,j,n;
    PARTICLE *pi, *pj;
    PARTICLE t;
    n = pkd->nLocal;
    int nSort = 0;

    if (iRung == 0) {
	i = 0;
	}
    else {
	i = n-(pkd->nVeryActive);
	}
    j = n - 1;

    while (i <= j) {
	pi = pkdParticle(pkd,i);
	if ( pi->iRung <= iRung ) ++i;
	else break;
	}
    while (i <= j) {
	pj = pkdParticle(pkd,j);
	if ( pj->iRung > iRung ) --j;
	else break;
	}

    if (i < j) {
	pi = pkdParticle(pkd,i);
	pj = pkdParticle(pkd,j);
	pkdSwapParticle(pkd,pi,pj);
	while (1) {
	    while ((p[++i].iRung <= iRung));
	    while (p[--j].iRung > iRung);
	    if (i < j) {
		t = p[i];
		p[i] = p[j];
		p[j] = t;
		}
	    else break;
	    }
	}
    nSort = n - i;
    return(nSort);
    }


void pkdGravVA(PKD pkd, int iRung) {
    int i, j, k, n, iStart;
    double x, y, z;
    double fac, d2, d1;
    PARTICLE *pi, *pj;
    double r_k = 3.0/pow(2.08,iRung-1);
    double r_kk = r_k/2.08;
    double r_kkk = r_kk/2.08;
    double fourh, dir, dir2;

    assert(iRung >= 1);
    iStart = pkd->nLocal - pkd->nVeryActive;

    n = pkdLocal(pkd);

    /* reset */
    for (i=iStart;i<n;++i) {
	pi = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    p->drmin = 1000.0;
	    for (k=0;k<3;k++) {
		p->a_VA[k] = 0.0;
		}
	    }
	}

    for (i=iStart;i<n;++i) {
	pi = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    for (k=0;k<pi->n_VA;k++) {
		j = pi->i_VA[k];
		pj = pkdParticle(pkd,j);
		if (pi->iOrder<pj->iOrder) { /* after three body encounter and a collision
						we sometimes obtain i<j, pi->iOrder = pj->iOrder*/
		    x = pi->r[0] - pj->r[0];
		    y = pi->r[1] - pj->r[1];
		    z = pi->r[2] - pj->r[2];
		    d2 = x*x + y*y +z*z;
		    d1 = sqrt(d2);

		    fourh = pi->fSoft + pj->fSoft;

		    if (d1 > fourh) {
			dir2 = 1.0/(d1*d2);
			}
		    else {

			if (pkd->param.bCollision) {

			    pkd->iCollisionflag = 1;
			    pi->iColflag = 1;
			    pi->iOrderCol = pj->iOrder;
			    pi->dtCol = 1.0*pi->iOrgIdx;
			    printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i piRung %d iRung %d\n",
				   d1,fourh,pi->iOrgIdx,pj->iOrgIdx, pi->iRung, iRung);
			    printf("i %d j %d inVA %d, jnVA %d \n",i, j, pi->n_VA, pj->n_VA);
			    }

			dir = 1.0/fourh;
			dir2 = dir*dir;
			d2 *= dir2;
			dir2 *= dir;
			d2 = 1.0 - d2;
			dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
			dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
			}


		    d1 /= pi->hill_VA[k]; /* distance normalized by hill */

		    if (d1 < r_kkk) {
			fac = 0.0;
			}
		    else if (d1 < r_kk && d1 >r_kkk) {
			fac = symfac*(1.0-d1/r_kk);
			fac *= fac*(2.0*fac -3.0);
			fac += 1.0;
			}
		    else if (d1 < r_k && d1 >r_kk) {
			fac = symfac*(1.0-d1/r_k);
			fac *= -fac*(2.0*fac -3.0);
			}
		    else if (d1 > r_k) {
			fac = 0.0;
			}

		    pi->drmin = (d1 < pi->drmin)?d1:pi->drmin;
		    pj->drmin = (d1 < pj->drmin)?d1:pj->drmin;

		    dir2 *= fac;
		    /*printf("iOrder %d %d, d2 %e iRung %d\n",
		      pi->iOrder,pj->iOrder,d2,iRung);*/
		    pi->a_VA[0] -= x*dir2*pj->fMass;
		    pi->a_VA[1] -= y*dir2*pj->fMass;
		    pi->a_VA[2] -= z*dir2*pj->fMass;
		    pj->a_VA[0] += x*dir2*pi->fMass;
		    pj->a_VA[1] += y*dir2*pi->fMass;
		    pj->a_VA[2] += z*dir2*pi->fMass;
		    }
		}
	    }
	}
    }

int pkdCheckDrminVA(PKD pkd, int iRung, int multiflag, int nMaxRung) {
    int i, j, k, n, iStart, nloop;
    double x, y, z;
    double d2;
    int iTempRung;
    PARTICLE *pi, *pj;
    double r_k = 3.0/pow(2.08,iRung);
    int iupdate = 0;

    assert(iRung >= 1);
    iStart = pkd->nLocal - pkd->nVeryActive;
    n = pkdLocal(pkd);

    for (i=iStart;i<n;++i) {
	pi = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,pi)) {
	    for (k=0;k<pi->n_VA;k++) {
		j = pi->i_VA[k];
		pj = pkdParticle(pkd,j);
		if (i<j) {
		    x = pi->r[0] - pj->r[0];
		    y = pi->r[1] - pj->r[1];
		    z = pi->r[2] - pj->r[2];
		    d2 = x*x + y*y +z*z;
		    d2 = sqrt(d2);
		    d2 /= pi->hill_VA[k]; /* distance normalized by hill */

		    /* d2 should be replaced by min.dist. during drift */

		    if (d2 < r_k) {
			iupdate = 1;
			pi->iRung = iRung +1;
			pj->iRung = iRung +1;
			nMaxRung = ((iRung+1) >= nMaxRung)?(iRung+1):nMaxRung;
			/* retrive */
			for (k=0;k<3;k++) {
			    pi->r[k] = pi->rb[k];
			    pi->v[k] = pi->vb[k];
			    pj->r[k] = pj->rb[k];
			    pj->v[k] = pj->vb[k];
			    }
			}
		    }
		}
	    }
	}
    /* If a close encounter with more than 2 particles exists, we synchronize
       iRungs of interacting particles to the highest iRung in them */
    if (multiflag && iupdate) {
	/*nloop = multiflag;*/
	nloop = 2;
	do {
	    for (i=iStart;i<n;++i) {
		if (pkdIsActive(pkd,&p[i])) {
		    if (pi->n_VA >= 2) {
			for (k=0;k<pi->n_VA;k++) {
			    iTempRung = p[pi->i_VA[k]].iRung;
			    if (iTempRung > pi->iRung) {
				if (iTempRung != (iRung +1)) {
				    printf("p_iOrder %d pi_n_VA %d pj_nVA %d iTempRung %d, iRung+1 %d \n",
					   pi->iOrder,pi->n_VA,p[pi->i_VA[k]].n_VA,iTempRung,iRung +1);
				    printf("too many particles in 3 hill radius !! \n");
				    assert(0);
				    }
				/* assert(iTempRung == (iRung +1));*/
				pi->iRung = iRung + 1;
				for (k=0;k<3;k++) {
				    pi->r[k] = pi->rb[k];
				    pi->v[k] = pi->vb[k];
				    }
				}
			    }
			}
		    }
		}
	    nloop --;
	    }
	while (nloop > 0);
	}
    return(nMaxRung);
    }

void pkdKickVA(PKD pkd, double dt) {
    int i, k, iStart;
    PARTICLE *p;
    iStart = pkd->nLocal - pkd->nVeryActive;

    for (i=iStart;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    for (k=0;k<3;k++) {
		p->v[k] += p->a_VA[k]*dt;
		}
	    }
	}
    }

int pkdDrminToRung(PKD pkd, int iRung, int iMaxRung, int *nRungCount) {
    int i, j, iTempRung;;
    PARTICLE *p ;

    /* R_k = 3.0/(2.08)^(k-1) with (k= 1,2, ...)*/
    /* note iMaxRung = pkd.param->iMaxRung */
    for (i=0;i<iMaxRung;++i) nRungCount[i] = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    iTempRung = iRung;
	    if (p->drmin > 3.0) {
		iTempRung = 0;
		}
	    else {
		iTempRung = floor(log(3.0/p->drmin)/log(2.08)) + 1;
		}
	    if (iTempRung >= iMaxRung) {
		iTempRung = iMaxRung-1;
		}

	    /* if min. dist. during drift is less than 3 Hill radius,
	    set iRung = 1 */
	    if (iTempRung == 0 && p->drmin2 < 3.0) {
		iTempRung = 1;
		}

	    /* retrive position and velocity of active particle to
	    those before drift*/
	    if (iTempRung > 0) {
		for (j=0;j<3;j++) {
		    p->r[j] = p->rb[j];
		    p->v[j] = p->vb[j];
		    }
		}
	    /*
	    ** Now produce a count of particles in rungs.
	    */
	    nRungCount[iTempRung] += 1;
	    p->iRung = iTempRung;

	    /*printf("iorder %d, drmin1 %e, drmin2 %e, \n",
	    p->iOrder, p->drmin, p->drmin2);*/
	    }

	}
    iTempRung = iMaxRung-1;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
    }

int pkdGetPointa(PKD pkd) {
    int i, j, k, n, iStart, iOrderTemp, n_VA;
    int multiflag = 0;
    int iTempRung, nloop;
    int bRung = 0;
    PARTICLE *p;
    char c;
    int iMaxRung = pkd->param.iMaxRung;
    int nRungCount[iMaxRung-1];/* 1 to iMaxRung-1 */

    iStart = pkd->nLocal - pkd->nVeryActive;
    n = pkdLocal(pkd);

    if (bRung) {
	for (i=1;i<iMaxRung;++i) nRungCount[i] = 0;
	}

    for (i=iStart;i<n;++i) {
	p = pkdParticle(pkd,i);
	assert(pkdIsActive(pkd,p));
	if (bRung) nRungCount[p->iRung] += 1;
	n_VA = p->n_VA;
	assert(n_VA >= 1);

	if (n_VA >= 2) multiflag = (multiflag < n_VA)?n_VA:multiflag; /* multiflag = largest n_VA */

	for (k=0;k<n_VA;k++) {
	    iOrderTemp = p->iOrder_VA[k];

	    for (j=iStart;j<n;++j) {
		if (i==j)continue;
		pj = pkdParticle(pkd,j);
		assert(pj->iOrder != p->iOrder);
		if (pj->iOrder == iOrderTemp) {
		    p->i_VA[k] = j;
		    break;
		    }
		}
	    }
	}

    if (bRung) {
	for (i=1;i<iMaxRung;++i) {
	    if (nRungCount[i] == 0) continue;
	    else c = ' ';
	    printf(" %c rung:%d %d\n",c,i,nRungCount[i]);
	    }
	printf("\n");
	}

    /* If a close encounter with more than 2 particles exists, we synchronize
       iRungs of interacting particles to the highest iRung in them */
    if (multiflag) {
	nloop = 2;
	do {
	    for (i=iStart;i<n;++i) {
		if (pkdIsActive(pkd,&p[i])) {
		    if (p[i].n_VA >= 2) {
			for (k=0;k<p[i].n_VA;k++) {
			    j = p[i].i_VA[k];
			    iTempRung = p[j].iRung;
			    p[i].iRung = (iTempRung >= p[i].iRung)?iTempRung:p[i].iRung;
			    p[j].iRung = p[i].iRung;
			    }
			}
		    }
		}
	    nloop --;
	    }
	while (nloop > 0);
	}
    return(multiflag);
    }

int pkdDrminToRungVA(PKD pkd, int iRung, int iMaxRung, int multiflag) {
    int i, j, k, iTempRung, iStart;
    int nMaxRung = 0;
    PARTICLE *pi, *pj;
    int nloop;
    /* note iMaxRung = pkd.param->iMaxRung -1 */

    int bRung = 0;
    char c;
    int nRungCount[iMaxRung-1];/* 1 to iMaxRung-1*/
    if (bRung) {
	for (i=1;i<iMaxRung;++i) nRungCount[i] = 0;
	}

    iStart = pkd->nLocal - pkd->nVeryActive;

    for (i=iStart;i<pkdLocal(pkd);++i) {
	pi = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)) {
	    /*assert(pi->n_VA >= 1); n_VA might be 0 after collision */
	    iTempRung = floor(log(3.0/pi->drmin)/log(2.08)) + 1;
	    /*if(iTempRung >= iMaxRung){
	    printf("iRung %d for particle %d larger than Max iRung %d\n",
	       iTempRung, pi->iOrder,iMaxRung-1);
	    iTempRung = iMaxRung;
		}*/
	    iTempRung = (iTempRung >= iMaxRung)?(iMaxRung-1):iTempRung;
	    /* iTempRung = (iTempRung >= 0)?iTempRung:0;
	       pi->iKickRung = iTempRung; */

	    iTempRung = (iTempRung >= iRung)?iTempRung:iRung;
	    pi->iRung = iTempRung;
	    nMaxRung = (iTempRung >= nMaxRung)?iTempRung:nMaxRung;
	    }
	if (bRung) nRungCount[pi->iRung] += 1;
	}

    if (bRung) {
	printf("nVA %d iRung %d\n",pkd->nVeryActive, iRung);
	for (i=1;i<iMaxRung;++i) {
	    if (nRungCount[i] == 0) continue;
	    else c = ' ';
	    printf(" %c rung:%d %d\n",c,i,nRungCount[i]);
	    }
	printf("\n");
	}
    /* If a close encounter with more than 2 particles exists, we synchronize
       iRungs of interacting particles to the highest iRung in them */
    if (multiflag) {
	nloop = 2;
	do {
	    for (i=iStart;i<pkdLocal(pkd);++i) {
		if (pkdIsActive(pkd,&p[i])) {
		    if (pi->n_VA >= 2) {
			for (k=0;k<pi->n_VA;k++) {
			    j = pi->i_VA[k];
			    pj= pkdParticle(pkd,j);
			    iTempRung = pj->iRung;
			    pi->iRung = (iTempRung >= pi->iRung)?iTempRung:pi->iRung;
			    pj->iRung = pi->iRung;
			    }
			}
		    }
		}
	    nloop --;
	    }
	while (nloop > 0);
	}
    return(nMaxRung);
    }


void pkdMomSun(PKD pkd,double momSun[]) {
    PARTICLE *p;
    int i,j,n;

    for (j=0;j<3;++j) {
	momSun[j] = 0.0;
	}
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) {
	    momSun[j] -= p->fMass*p->v[j];
	    }
	}
    }

void pkdDriftSun(PKD pkd,double vSun[],double dt) {
    PARTICLE *p;
    int i,j,n;

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) {
	    p->r[j] += vSun[j]*dt;
	    }
	}
    }

void pkdKeplerDrift(PKD pkd,double dt,double mu, int tag_VA) {
    /*
     * Use the f and g functions to advance an unperturbed orbit.
     */
    PARTICLE *p;
    int  i,j,n, iStart;
    int iflg = 0;
    double dttmp;
    mdlDiag(pkd->mdl, "Into pkdKeplerDrift \n");
    n = pkdLocal(pkd);

    if (tag_VA) {
	iStart = pkd->nLocal - pkd->nVeryActive;
	}
    else {
	iStart = 0;
	}

    for (i=iStart;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsActive(pkd,p)&&(p->iOrder>=0)) {

	    /* copy r and v before drift */
	    for (j=0;j<3;++j) {
		p->rb[j] = p->r[j];
		p->vb[j] = p->v[j];
		}

	    /*printf("before: iorder = %d, (x,y,z)= (%e,%e,%e), (vx,vy,vz)= (%e,%e,%e),ms =%e \n",
	    p->iOrder, p->r[0],p->r[1],p->r[2],
	    p->v[0],p->v[1],p->v[2], mu);*/

	    iflg = drift_dan(mu,p->r,p->v,dt); /* see kepler.c */

	    /* printf("exit drift_dan, iflg = %d, (x,y,z)= (%e,%e,%e),(vx,vy,vz)= (%e,%e,%e), ms =%e \n",
	    iflg, p->r[0],p->r[1],p->r[2],p->v[0],p->v[1],p->v[2],mu);*/

	    if (iflg != 0) {
		/*printf("call drift_dan*10, iflg = %d\n",iflg);*/
		dttmp = 0.1*dt;
		for (j=0;j<10;++j) {
		    iflg = drift_dan(mu,p->r,p->v,dttmp);
		    if (iflg != 0) {
			printf("exit drift_dan, iflg = %d, (x,y,z)= (%e,%e,%e),(vx,vy,vz)= (%e,%e,%e), m = %e, ms =%e \n", iflg, p->r[0],p->r[1],p->r[2],p->v[0],p->v[1],p->v[2],p->fMass,mu);
			}
		    assert(iflg == 0);
		    }
		}
	    }
	}
    }

#endif /* SYMBA */
#endif /* PLANETS */

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
	density = p->fDensity * pow(pvel->veldisp2,-1.5);
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
	density = p->fDensity * pow(pvel->veldisp2,-1.5);
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
	    double dx = dCenter[j] - p->r[j];
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
	    predicate = predicate && fabs(dCenter[j] - p->r[j]) <= dSize[j];
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
	dx = r[0] - p->r[0];
	dy = r[1] - p->r[1];
	dz = r[2] - p->r[2];

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
	dx = r[0] - p->r[0];
	dy = r[1] - p->r[1];
	dz = r[2] - p->r[2];

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
	    dPart[j] = p->r[j] - dP1[j];
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
	    dPart[j] = p->r[j] - dP1[j];
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
    r[0] = pLocal->r[0];
    r[1] = pLocal->r[1];
    r[2] = pLocal->r[2];
    *fPot= *pPotLocal;
    return nChecked;
    }
