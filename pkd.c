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
#include <time.h> /* added MZ */

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

double Zeit() { /* added MZ */
    struct timeval tv;

    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }

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

#define returnWarning(x) if (!(x)) { printf("assert failed: %s\n",#x); return 0; }

int pkdVerify(PKD pkd) {
    returnWarning(pkd->pStore[-1].iOrder==0x23456789abc);
    returnWarning(pkd->pStore[pkd->nStore].iOrder==0x23456789abc);
    returnWarning(pkd->pStore[-1].uRung==0x2d);
    returnWarning(pkd->pStore[pkd->nStore].uRung==0x2d);
    returnWarning(pkd->pStore[-1].bSrcActive==0);
    returnWarning(pkd->pStore[pkd->nStore].bSrcActive==0);
    returnWarning(pkd->pStore[-1].bDstActive==0);
    returnWarning(pkd->pStore[pkd->nStore].bDstActive==0);
    returnWarning(pkd->pStore[-1].iClass==0xdb);
    returnWarning(pkd->pStore[pkd->nStore].iClass==0xdb);
    returnWarning(pkd->pStore[-1].fBall==2.0);
    returnWarning(pkd->pStore[pkd->nStore].fBall==2.0);
    returnWarning(pkd->pStore[-1].fPot==4.0);
    returnWarning(pkd->pStore[pkd->nStore].fPot==4.0);
    returnWarning(pkd->pStore[-1].a[0]==8.0);
    returnWarning(pkd->pStore[pkd->nStore].a[0]==8.0);
    returnWarning(pkd->pStore[-1].a[1]==16.0);
    returnWarning(pkd->pStore[pkd->nStore].a[1]==16.0);
    returnWarning(pkd->pStore[-1].a[2]==32.0);
    returnWarning(pkd->pStore[pkd->nStore].a[2]==32.0);
    returnWarning(pkd->pStore[-1].fDensity==64.0);
    returnWarning(pkd->pStore[pkd->nStore].fDensity==64.0);
    returnWarning(pkd->pStore[-1].v[0]==128.0);
    returnWarning(pkd->pStore[pkd->nStore].v[0]==128.0);
    returnWarning(pkd->pStore[-1].v[1]==256.0);
    returnWarning(pkd->pStore[pkd->nStore].v[1]==256.0);
    returnWarning(pkd->pStore[-1].v[2]==512.0);
    returnWarning(pkd->pStore[pkd->nStore].v[2]==512.0);
    return 1;
    }

void pkdInitialize(PKD *ppkd,MDL mdl,int nStore,int nBucket,FLOAT *fPeriod,
		   uint64_t nDark,uint64_t nGas,uint64_t nStar) {
    PKD pkd;
    int j,ism;

    pkd = (PKD)malloc(sizeof(struct pkdContext));
    mdlassert(mdl,pkd != NULL);
    pkd->mdl = mdl;
    pkd->idSelf = mdlSelf(mdl);
    pkd->nThreads = mdlThreads(mdl);
    pkd->nStore = nStore;
    pkd->nLocal = 0;
    pkd->nDark = nDark;
    pkd->nGas = nGas;
    pkd->nStar = nStar;
    pkd->nMaxOrderGas = nGas;
    pkd->nMaxOrderDark = nGas + nDark;
    pkd->nRejects = 0;
    for (j=0;j<3;++j) {
	pkd->fPeriod[j] = fPeriod[j];
	}

    pkd->uMinRungActive  = 0;
    pkd->uMaxRungActive  = 255;
    pkd->uRungVeryActive = 255;

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
    */
    pkd->pStore = mdlMalloc(pkd->mdl,(nStore+2)*sizeof(PARTICLE));
    mdlassert(mdl,pkd->pStore != NULL);
    pkd->pStore++;
    /*
    ** NOTE: the particle store now had an extra field below the first
    ** position that we use as a sentinal.  That, combined with the one
    ** at the end is used to verify that there are no memory bounds issues.
    ** The free call MUST take into account the pointer magic.
    */
    pkd->pStore[-1].iOrder     = pkd->pStore[nStore].iOrder     = 0x23456789abc;
    pkd->pStore[-1].uRung      = pkd->pStore[nStore].uRung      = 0x2d;
    pkd->pStore[-1].bSrcActive = pkd->pStore[nStore].bSrcActive = 0;
    pkd->pStore[-1].bDstActive = pkd->pStore[nStore].bDstActive = 0;
    pkd->pStore[-1].iClass     = pkd->pStore[nStore].iClass     = 0xdb;
    pkd->pStore[-1].fBall      = pkd->pStore[nStore].fBall      = 2.0;
    pkd->pStore[-1].fPot       = pkd->pStore[nStore].fPot       = 4.0;
    pkd->pStore[-1].a[0]       = pkd->pStore[nStore].a[0]       = 8.0;
    pkd->pStore[-1].a[1]       = pkd->pStore[nStore].a[1]       = 16.0;
    pkd->pStore[-1].a[2]       = pkd->pStore[nStore].a[2]       = 32.0;
    pkd->pStore[-1].fDensity   = pkd->pStore[nStore].fDensity   = 64.0;
    pkd->pStore[-1].v[0]       = pkd->pStore[nStore].v[0]       = 128.0;
    pkd->pStore[-1].v[1]       = pkd->pStore[nStore].v[1]       = 256.0;
    pkd->pStore[-1].v[2]       = pkd->pStore[nStore].v[2]       = 512.0;
    assert(pkdVerify(pkd));
    /*
    ** We support up to 256 classes
    */
    pkd->pClass = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
    mdlassert(mdl,pkd->pClass != NULL);
    for (j=0;j<PKD_MAX_CLASSES;j++)
	pkd->pClass[j].fMass = pkd->pClass[j].fSoft = pkd->pClass[j].fSoft0 = -1.0;
    pkd->nClasses = 0;
    /*
    ** Now also allocate all node storage here.
    ** We guess that the upper bound is based on the number of particles in
    ** a bucket. The mean number of particles per bucket is always somewhat
    ** less than nBucket, and roughly given by nBucket-sqrt(nBucket).  For
    ** small numbers of particles, we must correct for the minimum cell size.
    */
    /*pkd->nMaxNodes = (int)ceil(3.0*nStore/floor(nBucket - sqrt(nBucket)));*/
    /* j is an estimate of the lower limit of the total number of cells */
    j = 2.0 / (PKD_MAX_CELL_SIZE*PKD_MAX_CELL_SIZE*PKD_MAX_CELL_SIZE*mdlThreads(mdl));
    pkd->nMaxNodes = (int)ceil(15.0*nStore/floor(nBucket-sqrt(nBucket))) + j;
    if ( pkd->nMaxNodes < j ) pkd->nMaxNodes = j;
    pkd->kdNodes = mdlMalloc(pkd->mdl,pkd->nMaxNodes*sizeof(KDN));
    mdlassert(mdl,pkd->kdNodes != NULL);
    /*
    ** pLite particles are also allocated and are quicker when sorting particle
    ** type operations such as tree building and domain decomposition are being
    ** performed.
    */
    pkd->pLite = malloc((nStore+1)*sizeof(PLITE));
    mdlassert(mdl,pkd->pLite != NULL);
    pkd->nNodes = 0;
    pkd->kdTop = NULL;
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
    }


void pkdFinish(PKD pkd) {
    int ism;

    if (pkd->kdNodes) {
	/*
	** Close caching space and free up nodes.
	*/
	if (pkd->nNodes > 0)
	    mdlFinishCache(pkd->mdl,CID_CELL);
	mdlFree(pkd->mdl,pkd->kdNodes);
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
    if (pkd->kdTop) free(pkd->kdTop);
    free(pkd->ew.ewt);
    free(pkd->pClass);
    mdlFree(pkd->mdl,pkd->pStore-1);
    free(pkd->pLite);
    free(pkd->piActive);
    free(pkd->piInactive);
    csmFinish(pkd->param.csm);
    free(pkd);
    }

static uint8_t getClass( PKD pkd, FLOAT fMass, FLOAT fSoft ) {
    int i;

    /* TODO: This is a linear search which is fine for a small number of classes */
    for ( i=0; i<pkd->nClasses; i++ )
	if ( pkd->pClass[i].fMass == fMass && pkd->pClass[i].fSoft0 == fSoft )
	    return i;
    assert( pkd->nClasses < PKD_MAX_CLASSES );
    i = pkd->nClasses++;
    pkd->pClass[i].fSoft = pkd->pClass[i].fSoft0 = fSoft;
    pkd->pClass[i].fMass = fMass;
    return i;
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
		if ( pClass[j].fMass==pkd->pClass[i].fMass && pClass[j].fSoft==pkd->pClass[i].fSoft )
		    break;
	    assert(j<n);
	    map[i] = j;
	    }

	/* Now update the class with the new value */
	for (i=0;i<pkd->nLocal;++i) {
	    p = &pkd->pStore[i];
	    assert( p->iClass <= pkd->nClasses );
	    p->iClass = map[p->iClass];
	    }
	}

    /* Finally, set the new class table */
    for ( i=0; i<n; i++ )
	pkd->pClass[i] = pClass[i];
    pkd->nClasses = n;
    }

void pkdSeek(PKD pkd,FILE *fp,uint64_t nStart,int bStandard,int bDoublePos,int bNoHeader) {
    off_t MAX_OFFSET = 2147483640;
    off_t lStart;
    int iErr;
    /*
    ** Seek according to true XDR size structures when bStandard is true.
    ** This may be a bit dicey, but it should work as long
    ** as no one changes the tipsy binary format!
    */
    if (bNoHeader) lStart = 0;
    else if (bStandard) lStart = 32;
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

    pkd->nClasses = 1;
    pkd->pClass[0].fSoft = pkd->pClass[0].fSoft0 = fSoft;
    pkd->pClass[0].fMass = fMass;

    for ( i=0; i<n1; i++ ) {
	for ( j=0; j<n2; j++ ) {
	    for ( k=0; k<n3; k++ ) {
		p = &pkd->pStore[pi];
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
		p->iOrder = 0; /* FIXME */
		p->iClass = 0;
		pi++;
		}
	    }
	}

    assert( pi == pkd->nLocal );
    }

#endif


#ifdef USE_HDF5
void pkdReadHDF5(PKD pkd, IOHDF5 io, double dvFac,
		 uint64_t nStart, int nLocal ) {
    PARTICLE *p;
    FLOAT dT1, dT2;
    FLOAT fSoft, fMass;
    uint64_t iOrder;
    int i, j;
    IOHDF5V ioPot;

//    ioPot  = ioHDFF5NewVector( io, "potential",IOHDF5_SINGLE );


    ioHDF5SeekDark( io, nStart );
//    ioHDF5SeekVector( ioPot, nStart );

    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[pkd->nLocal+i];
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	for (j=0;j<3;++j) {
	    p->a[j] = 0.0;
	    }
	}

    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[pkd->nLocal+i];
	p->iOrder = nStart + i;

	if (pkdIsDark(pkd,p)) {
	    ioHDF5GetDark( io, &iOrder, p->r, p->v,
			   &fMass, &fSoft, &p->fPot );
	    }
	else if (pkdIsGas(pkd,p)) {
	    ioHDF5GetGas( io, &iOrder, p->r, p->v,
			  &fMass, &fSoft, &p->fPot,
			  &dT1, &dT2 );
	    }
	else if (pkdIsStar(pkd,p)) {
	    ioHDF5GetStar( io, &iOrder, p->r, p->v,
			   &fMass, &fSoft, &p->fPot,
			   &dT1, &dT2 );
	    }
	else mdlassert(pkd->mdl,0);
	if (fSoft < sqrt(2.0e-38)) { /* set minimum softening */
	    fSoft = sqrt(2.0e-38);
	    }
	for (j=0;j<3;++j) p->v[j] *= dvFac;
	p->iClass = getClass(pkd,fMass,fSoft);
	p->iOrder = iOrder;
//	p->fPot = ioHDF5GetVector(ioPot);
	p->fPot = 0.0;
	}

    pkd->nLocal += nLocal;
    pkd->nActive += nLocal;



    }
#endif


#ifdef USE_MDL_IO
void pkdIOInitialize( PKD pkd, int nLocal) {
    int i, j;
    PARTICLE *p;

    pkd->nLocal = pkd->nActive = nLocal;

    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[i];
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	p->iClass = 0;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	for (j=0;j<3;++j) {
	    p->a[j] = 0.0;
	    }
	}

    }
#endif

void pkdReadTipsy(PKD pkd,char *pszFileName, uint64_t nStart,int nLocal,
		  int bStandard,double dvFac,int bDoublePos,int bNoHeader) {
    FILE *fp;
    int i,j;
    PARTICLE *p;
    struct dark_particle dp;
    struct gas_particle gp;
    struct star_particle sp;
    float fTmp;
    FLOAT fMass, fSoft;
    double dTmp;
    float mass=0.0;

    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[pkd->nLocal+i];
	p->uRung = p->uNewRung = 0;
	p->bSrcActive = p->bDstActive = 1;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	/*
	** Clear the accelerations so that the timestepping calculations do not
	** get funny uninitialized values!
	*/
	for (j=0;j<3;++j) {
	    p->a[j] = 0.0;
	    }
	}
    /*
    ** Seek past the header and up to nStart.
    */
    fp = fopen(pszFileName,"r");
    mdlassert(pkd->mdl,fp != NULL);
    /*
    ** Seek to right place in file.
    */
    pkdSeek(pkd,fp,nStart,bStandard,bDoublePos,bNoHeader);
    /*
    ** Read Stuff!
    */
    if (bStandard) {
	FLOAT vTemp;
	XDR xdrs;
	xdrstdio_create(&xdrs,fp,XDR_DECODE);
	for (i=0;i<nLocal;++i) {
	    p = &pkd->pStore[pkd->nLocal+i];
	    p->iOrder = nStart + i;
	    if (pkdIsDark(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		fMass = fTmp;
		mass += fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;
		    }
		xdr_float(&xdrs,&fTmp);
		fSoft = fTmp;
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else if (pkdIsGas(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		fMass = fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;
		    }

		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		fSoft = fTmp;
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else if (pkdIsStar(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		fMass = fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;
		    }
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		fSoft = fTmp;
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else mdlassert(pkd->mdl,0);
	    p->iClass = getClass(pkd,fMass,fSoft);
	    }
	xdr_destroy(&xdrs);
	}
    else {
	for (i=0;i<nLocal;++i) {
	    p = &pkd->pStore[pkd->nLocal+i];
	    p->iOrder = nStart + i;
	    if (pkdIsDark(pkd,p)) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = dp.pos[j];
		    p->v[j] = dvFac*dp.vel[j];
		    }
		fMass = dp.mass;
		mass += dp.mass;
		fSoft = dp.eps;
		p->fPot = dp.phi;
		}
	    else if (pkdIsGas(pkd,p)) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = gp.pos[j];
		    p->v[j] = dvFac*gp.vel[j];
		    }
		fMass = gp.mass;
		fSoft = gp.hsmooth;
		p->fPot = gp.phi;
		}
	    else if (pkdIsStar(pkd,p)) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = sp.pos[j];
		    p->v[j] = dvFac*sp.vel[j];
		    }
		fMass = sp.mass;
		fSoft = sp.eps;
		p->fPot = sp.phi;
		}
	    else mdlassert(pkd->mdl,0);
	    p->iClass = getClass(pkd,fMass,fSoft);
	    }
	}

    pkd->nLocal += nLocal;
    pkd->nActive += nLocal;

    fclose(fp);
    }


void pkdCalcBound(PKD pkd,BND *pbnd) {
    double dMin[3],dMax[3];
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->nLocal > 0);
    for (j=0;j<3;++j) {
	dMin[j] = pkd->pStore[i].r[j];
	dMax[j] = pkd->pStore[i].r[j];
	}
    for (++i;i<pkd->nLocal;++i) {
	pkdMinMax(pkd->pStore[i].r,dMin,dMax);
	}
    for (j=0;j<3;++j) {
	pbnd->fCenter[j] = pkd->bnd.fCenter[j] = 0.5*(dMin[j] + dMax[j]);
	pbnd->fMax[j] = pkd->bnd.fMax[j] = 0.5*(dMax[j] - dMin[j]);
	}
    }


void pkdEnforcePeriodic(PKD pkd,BND *pbnd) {
    int i,j;

    for (i=0;i<pkd->nLocal;++i) {
	for (j=0;j<3;++j) {
	    if (pkd->pStore[i].r[j] < pbnd->fCenter[j] - pbnd->fMax[j]) pkd->pStore[i].r[j] += 2*pbnd->fMax[j];
	    else if (pkd->pStore[i].r[j] >= pbnd->fCenter[j] + pbnd->fMax[j]) pkd->pStore[i].r[j] -= 2*pbnd->fMax[j];
	    /*
	    ** If it still doesn't lie in the "unit" cell then something has gone quite wrong with the 
	    ** simulation. Either we have a super fast particle or the initial condition is somehow not conforming
	    ** to the specified periodic box in a gross way.
	    */
	    mdlassert(pkd->mdl,((pkd->pStore[i].r[j] >= pbnd->fCenter[j] - pbnd->fMax[j])&&
				(pkd->pStore[i].r[j] < pbnd->fCenter[j] + pbnd->fMax[j])));

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
    PARTICLE *p = pkd->pStore;
    PLITEDD *pl = (PLITEDD *)pkd->pLite;
    uint64_t uMask;
    float x,y,z;
    int i,j,bits,iShift;

    for (i=0;i<pkd->nLocal;++i) {
	/*
	** For now we just assume the particles are coming from a standard
	** cosmological box. We scale the volume by a factor of 0.99 to be
	** certain that fast particles are still captured by the domain
	** decomposition.
	*/
	x = 0.99*p[i].r[0] + 1.5;
	if (x < 1.0) x = 1.0;
	else if (x >= 2.0) x = 2.0;
	y = 0.99*p[i].r[1] + 1.5;
	if (y < 1.0) y = 1.0;
	else if (y >= 2.0) y = 2.0;
	z = 0.99*p[i].r[2] + 1.5;
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
    int i;

    *pnLow = 0;
    *pnHigh = 0;
    for (i=0;i<pkd->nLocal;++i) {
	if (pkdIsVeryActive(pkd,pkd->pStore+i)) {
	    if (pkd->pStore[i].r[d] < fSplit) *pnLow += 1;
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
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].r[d] >= fSplit,
	      pkd->pStore[j].r[d] < fSplit);
    return(i);
    }


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j) {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].r[d] < fSplit,
	      pkd->pStore[j].r[d] >= fSplit);
    return(i);
    }


int pkdLowerPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE pTemp;

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1));
	    }
	}
    return(i);
    }


int pkdUpperPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE pTemp;

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1));
	    }
	}
    return(i);
    }


int pkdLowerOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].iOrder >= nOrdSplit,
	      pkd->pStore[j].iOrder < nOrdSplit);
    return(i);
    }


int pkdUpperOrdPart(PKD pkd,uint64_t nOrdSplit,int i,int j) {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].iOrder < nOrdSplit,
	      pkd->pStore[j].iOrder >= nOrdSplit);
    return(i);
    }


int pkdActiveOrder(PKD pkd) {
    PARTICLE pTemp;
    int i=0;
    int j=pkdLocal(pkd)-1;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkdIsActive(pkd,&(pkd->pStore[i])),
	      !pkdIsActive(pkd,&(pkd->pStore[j])));
    return (pkd->nActive = i);
    }


int pkdColRejects_Active_Inactive(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide) {
    PARTICLE pTemp;
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
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;
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
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
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
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
    return(pkd->nRejects);
    }


int pkdSwapRejects(PKD pkd,int idSwap) {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;

    if (idSwap != -1) {
	nBuf = (pkdSwapSpace(pkd))*sizeof(PARTICLE);
	nOutBytes = pkd->nRejects*sizeof(PARTICLE);
	mdlassert(pkd->mdl,pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
	mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkdLocal(pkd)],
		nOutBytes,&nSndBytes,&nRcvBytes);
	pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
	pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
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
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
	    &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,nSndBytes/sizeof(PARTICLE) == pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/sizeof(PARTICLE);
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

int pkdNumSrcActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    int i, n;
    for (n=0,i=0;i<pkdLocal(pkd);++i)
	if ( pkdIsSrcActive(&pkd->pStore[i],uRungLo,uRungHi) ) n++;
    return n;
    }

int pkdNumDstActive(PKD pkd,uint8_t uRungLo,uint8_t uRungHi) {
    int i, n;
    for (n=0,i=0;i<pkdLocal(pkd);++i)
	if ( pkdIsDstActive(&pkd->pStore[i],uRungLo,uRungHi) ) n++;
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
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
    return(pkd->nRejects);
    }


int cmpParticles(const void *pva,const void *pvb) {
    PARTICLE *pa = (PARTICLE *)pva;
    PARTICLE *pb = (PARTICLE *)pvb;

    return(pa->iOrder - pb->iOrder);
    }


void pkdLocalOrder(PKD pkd) {
    qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
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
    assert( pkd->pStore != NULL );

    for ( i=0; i<nMax; i++ ) {
	local_t I = *iIndex + i;
	PARTICLE *p = pkd->pStore + I;

	for ( d=0; d<3; d++ ) {
	    p->r[d]  = io[i].r[d];
	    p->v[d]  = io[i].v[d] * dvFac; //FIXME: times??
	    p->a[d]  = 0.0;
	    }
	p->iOrder = io[i].iOrder;
	p->iClass = getClass(pkd,io[i].fMass,io[i].fSoft);
	p->fDensity = io[i].fDensity;
	p->fPot = io[i].fPot;
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
    local_t i;
    int nCopied, d;

    mdlassert(pkd->mdl,*iIndex<=pkd->nLocal);

    for ( i=*iIndex,nCopied=0; nCopied < nMax && i < pkd->nLocal; i++ ) {
	/* Not a particle of interest? */
	if ( pkd->pStore[i].iOrder<iMinOrder || pkd->pStore[i].iOrder>=iMaxOrder)
	    continue;

	/* We should have certain special cases here */
	mdlassert( pkd->mdl, pkdIsDark(pkd,pkd->pStore+i) );

	for ( d=0; d<3; d++ ) {
	    io[nCopied].r[d] = pkd->pStore[i].r[d];
	    io[nCopied].v[d] = pkd->pStore[i].v[d] * dvFac;
	    }
	io[nCopied].iOrder= pkd->pStore[i].iOrder;
	io[nCopied].fMass = pkdMass(pkd,&pkd->pStore[i]);
	io[nCopied].fSoft = pkdSoft(pkd,&pkd->pStore[i]);
	io[nCopied].fDensity = pkd->pStore[i].fDensity;
	io[nCopied].fPot = pkd->pStore[i].fPot;
	nCopied++;
	}
    *iIndex = i;
    return nCopied;
    }

#ifdef USE_HDF5
void pkdWriteHDF5(PKD pkd, IOHDF5 io, IOHDF5V ioDen, IOHDF5V ioPot, double dvFac) {
    PARTICLE *p;
    FLOAT v[3], fSoft, fMass;
    int i;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = &pkd->pStore[i];
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;

	v[0] = p->v[0] * dvFac;
	v[1] = p->v[1] * dvFac;
	v[2] = p->v[2] * dvFac;
	fSoft = pkdSoft0(pkd,p);
	fMass = pkdMass(pkd,p);
	if (pkdIsDark(pkd,p)) {
	    ioHDF5AddDark(io,p->iOrder,p->r,v,
			  fMass,fSoft,p->fPot );
	    }
	else if (pkdIsGas(pkd,p)) {
	    assert(0);
	    /* Why are temp and metals always set to zero? */
	    ioHDF5AddGas( io,p->iOrder,p->r,v,
			  fMass,fSoft,p->fPot,0.0,0.0);
	    }
	else if (pkdIsStar(pkd,p)) {
	    assert(0);
	    /* Why are metals and tform always set to zero? */
	    ioHDF5AddStar(io, p->iOrder, p->r, v,
			  fMass,fSoft,p->fPot,0.0,0.0);
	    }
	else mdlassert(pkd->mdl,0);

	ioHDF5AddVector( ioDen, p->iOrder, p->fDensity );
	ioHDF5AddVector( ioPot, p->iOrder, p->fPot );
	}
    }


#endif

uint32_t pkdWriteTipsy(PKD pkd,char *pszFileName,uint64_t nStart,
		       int bStandard,double dvFac,int bDoublePos) {
    PARTICLE *p;
    FILE *fp;
    int i,j;
    struct dark_particle dp;
    struct gas_particle gp;
    struct star_particle sp;
    int nout;
    float fTmp;
    double dTmp;
    uint32_t nCount;

    /*
    ** Seek past the header and up to nStart.
    */
    fp = fopen(pszFileName,"r+");
    mdlassert(pkd->mdl,fp != NULL);
    pkdSeek(pkd,fp,nStart,bStandard,bDoublePos,0);

    nCount = 0;
    if (bStandard) {
	FLOAT vTemp;
	XDR xdrs;
	/*
	** Write Stuff!
	*/
	xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	for (i=0;i<pkdLocal(pkd);++i) {
	    p = &pkd->pStore[i];
	    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	    nCount++;
	    if (pkdIsDark(pkd,p)) {
		fTmp = pkdMass(pkd,p);
		assert(fTmp > 0.0);
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
		fTmp = pkdSoft0(pkd,p);
		assert(fTmp > 0.0);
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else if (pkdIsGas(pkd,p)) {
		fTmp = pkdMass(pkd,p);
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
		fTmp = p->fDensity;
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = pkdSoft0(pkd,p);
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else if (pkdIsStar(pkd,p)) {
		fTmp = pkdMass(pkd,p);
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = pkdSoft0(pkd,p);
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else mdlassert(pkd->mdl,0);
	    }
	xdr_destroy(&xdrs);
	}
    else {
	/*
	** Write Stuff!
	*/
	for (i=0;i<pkdLocal(pkd);++i) {
	    p = &pkd->pStore[i];
	    if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	    nCount++;
	    if (pkdIsDark(pkd,p)) {
		for (j=0;j<3;++j) {
		    dp.pos[j] = p->r[j];
		    dp.vel[j] = dvFac*p->v[j];
		    }
		dp.mass = pkdMass(pkd,p);
		dp.eps = pkdSoft0(pkd,p);
		dp.phi = p->fPot;
		IOCheck(fwrite(&dp,sizeof(struct dark_particle),1,fp));
		}
	    else if (pkdIsGas(pkd,p)) {
		for (j=0;j<3;++j) {
		    gp.pos[j] = p->r[j];
		    gp.vel[j] = dvFac*p->v[j];
		    }
		gp.mass = pkdMass(pkd,p);
		gp.hsmooth = pkdSoft0(pkd,p);
		gp.phi = p->fPot;
		gp.rho = p->fDensity;
		gp.temp = 0.0;
		gp.metals = 0.0;
		IOCheck(fwrite(&gp,sizeof(struct gas_particle),1,fp));
		}
	    else if (pkdIsStar(pkd,p)) {
		for (j=0;j<3;++j) {
		    sp.pos[j] = p->r[j];
		    sp.vel[j] = dvFac*p->v[j];
		    }
		sp.mass = pkdMass(pkd,p);
		sp.eps = pkdSoft0(pkd,p);
		sp.phi = p->fPot;
		sp.metals = 0.0;
		sp.tform = 0.0;
		IOCheck(fwrite(&sp,sizeof(struct star_particle),1,fp));
		}
	    else mdlassert(pkd->mdl,0);
	    }
	}
    nout = fclose(fp);
    mdlassert(pkd->mdl,nout == 0);
    return nCount;
    }


void pkdSetSoft(PKD pkd,double dSoft) {
    PARTICLE *p;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
    if (dSoft < sqrt(2.0e-38)) { /* set minimum softening */
	dSoft = sqrt(2.0e-38);
	}
    for (i=0;i<pkd->nClasses;i++)
	pkd->pClass[i].fSoft = dSoft;
    }

#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul) {
    PARTICLE *p;
    int i,n;

    p = pkd->pStore;

    n = pkd->nClasses;
    mdlassert(pkd->mdl,dFac > 0);
    if (bSoftMaxMul) {
	for (i=0;i<n;++i) {
	    mdlassert(pkd->mdl,pkd->pClass[i].fSoft0 > 0);
	    pkd->pClass[i].fSoft = pkd->pClass[i].fSoft0*dFac;
	    mdlassert(pkd->mdl,pkd->pClass[i].fSoft > 0);
	    }
	}
    else {
	mdlassert(pkd->mdl,dSoftMax > 0);
	for (i=0;i<n;++i) {
	    mdlassert(pkd->mdl,pkd->pClass[i].fSoft0 > 0);
	    pkd->pClass[i].fSoft = pkd->pClass[i].fSoft0*dFac;
	    if (pkd->pClass[i].fSoft > dSoftMax) pkd->pClass[i].fSoft = dSoftMax;
	    mdlassert(pkd->mdl,pkd->pClass[i].fSoft > 0);
	    }
	}
    }
#endif

#ifdef USE_BSC_trace
static int foo = 0;
#endif

void
pkdGravAll(PKD pkd,uint8_t uRungLo,uint8_t uRungHi,double dTime,int nReps,int bPeriodic,
	   int iOrder,int bEwald,double fEwCut,double fEwhCut,int *nActive,
	   double *pdPartSum, double *pdCellSum,CASTAT *pcs, double *pdFlop) {
    int bVeryActive = 0;

    pkdClearTimer(pkd,1);
#ifdef INSTRUMENT
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
    mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
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
    FLOAT rx,ry,rz,vx,vy,vz,fMass;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
    *T = 0.0;
    *U = 0.0;
    *Eth = 0.0;
    L[0] = L[1] = L[2] = 0;
    for (i=0;i<n;++i) {
	fMass = pkdMass(pkd,&p[i]);
	rx = p[i].r[0]; ry = p[i].r[1]; rz = p[i].r[2];
	vx = p[i].v[0]; vy = p[i].v[1]; vz = p[i].v[2];
	*T += 0.5*fMass*(vx*vx + vy*vy + vz*vz);
	*U += 0.5*fMass*p[i].fPot;
	L[0] += fMass*(ry*vz - rz*vy);
	L[1] += fMass*(rz*vx - rx*vz);
	L[2] += fMass*(rx*vy - ry*vx);
	}
    }


/*
** Drift particles whose Rung falls between uRungLo (large step) and uRungHi (small step) inclusive,
** and those whose destination activity flag is set.
**
** Note that the drift funtion no longer wraps the particles around the periodic "unit" cell. This is
** now done by Domain Decomposition only.
*/
void pkdDrift(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    int i,j,n;
    double dMin[3],dMax[3];

    mdlDiag(pkd->mdl, "Into pkdDrift\n");
#ifdef USE_BSC
    MPItrace_event(10000,4);
#endif
    for (j=0;j<3;++j) {
	dMin[j] = pkd->bnd.fCenter[j] - pkd->bnd.fMax[j];
	dMax[j] = pkd->bnd.fCenter[j] + pkd->bnd.fMax[j];
	}
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsDstActive(&pkd->pStore[i],uRungLo,uRungHi)) {
	    /*
	    ** Update particle positions
	    */
	    for (j=0;j<3;++j) {
		pkd->pStore[i].r[j] += dDelta*pkd->pStore[i].v[j];
		}
	    pkdMinMax(pkd->pStore[i].r,dMin,dMax);
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

    double time1,time2; /* added MZ 1.6.2006 */

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
	pkdDrift(pkd,dTime,dDriftFac,iRungVeryActive+1,MAX_RUNG);
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

	    time1 = Zeit(); /* added MZ 1.6.2006 */

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

	    time2 = Zeit();
	    if (pkd->param.bVDetails)
		printf("Time: %g\n",time2-time1);
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

    double time1,time2; /* added MZ 1.6.2006 */

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

	    time1 = Zeit(); /* added MZ 1.6.2006 */

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

	    time2 = Zeit();
	    if (pkd->param.bVDetails)
		printf("Time: %g\n",time2-time1);

	    }
	}
    }

void
pkdCopy0(PKD pkd,double dTime) {
    PARTICLE *p;
    int i,j,n;

    mdlDiag(pkd->mdl, "Into pkdCopy0\n");

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    for (j=0;j<3;++j) {
		p[i].r0[j] = p[i].r[j];
		p[i].v0[j] = p[i].v[j];
		p[i].a0[j] = p[i].a[j];
		p[i].ad0[j] = p[i].ad[j];
		}
	    p[i].dTime0 = dTime;
#ifdef PLANETS
	    if (pkd->param.bCollision) {
		p[i].iColflag = 0;  /* just in case */
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

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    dt =  dTime - p[i].dTime0;
	    for (j=0;j<3;++j) {
		p[i].r[j] = p[i].r0[j] + dt*(p[i].v0[j] + 0.5*dt*(p[i].a0[j]+dt*p[i].ad0[j]/3.0));
		p[i].v[j] = p[i].v0[j] + dt*(p[i].a0[j]+0.5*dt*p[i].ad0[j]);
		p[i].rp[j] = p[i].r[j];
		p[i].vp[j] = p[i].v[j];
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

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    dt =  dTime - p[i].dTime0;
	    if (pkd->param.bAarsethStep) {
		a0d2  = 0.0;
		a1d2  = 0.0;
		a2d2  = 0.0;
		a3d2  = 0.0;
		dti = 1.0/dt;
		dt2i = dti*dti;
		}

	    for (j=0;j<3;++j) {
		am = p[i].a0[j]-p[i].a[j];

		if (pkd->param.bAarsethStep) {
		    add = 2.0*(3.0*am*dt2i+(2.0*p[i].ad0[j]+p[i].ad[j])*dti);
		    addd = 6.0*(2.0*am*dti+(p[i].ad0[j]+p[i].ad[j]))*dt2i;
		    add += dt*addd;

		    a0d2 += p[i].a[j]*p[i].a[j];
		    a1d2 += p[i].ad[j]*p[i].ad[j];
		    a2d2 += add*add;
		    a3d2 += addd*addd;
		    }

		add = 16.0*am+(13.0*p[i].ad0[j]+3.0*p[i].ad[j])*dt;
		addd = 6.0*am+(5.0*p[i].ad0[j]+p[i].ad[j])*dt;
		p[i].r[j] = p[i].rp[j] - add/120.0*dt*dt;
		p[i].v[j] = p[i].vp[j] - addd/12.0*dt;

		}
	    if (pkd->param.bAarsethStep) {
		p[i].dtGrav = (sqrt(a0d2*a2d2)+a1d2)/(sqrt(a1d2*a3d2)+a2d2);
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

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    dt =  dTime - p[i].dTime0;
	    r2  = 0.0;
	    rv  = 0.0;
	    for (j=0;j<3;++j) {
		r2 += p[i].r[j]*p[i].r[j];
		rv += p[i].r[j]*p[i].v[j];
		}
	    r2 = 1.0/r2;
	    r3i = r2*sqrt(r2)*dSunMass;
	    r5i = 3.0*r3i*r2*rv;

	    for (j=0;j<3;++j) {
		p[i].a[j] = -p[i].r[j]*r3i + p[i].app[j];
		p[i].ad[j] = -(p[i].v[j]*r3i - p[i].r[j]*r5i) + p[i].adpp[j];
		am = p[i].a0[j]-p[i].a[j];
		add = 16.0*am+(13.0*p[i].ad0[j]+3.0*p[i].ad[j])*dt;
		addd = 6.0*am+(5.0*p[i].ad0[j]+p[i].ad[j])*dt;
		p[i].r[j] = p[i].rp[j] - add/120.0*dt*dt;
		p[i].v[j] = p[i].vp[j] - addd/12.0*dt;

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

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) continue;
	dt =  dTime - p[i].dTime0;
	for (j=0;j<3;++j) {
	    p[i].r[j] = p[i].r0[j] + dt*(p[i].v0[j] + 0.5*dt*(p[i].a0[j]+dt*p[i].ad0[j]/3.0));
	    p[i].v[j] = p[i].v0[j] + dt*(p[i].a0[j]+0.5*dt*p[i].ad0[j]);
	    p[i].rp[j] = p[i].r[j];
	    p[i].vp[j] = p[i].v[j];
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdPredictorInactive\n");
    }

void pkdAarsethStep(PKD pkd,double dEta) {
    double dT;
    int i;

    for (i=0;i<pkdLocal(pkd);i++) {
	if (pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    mdlassert(pkd->mdl, pkd->pStore[i].dtGrav > 0);
	    dT = dEta*sqrt(pkd->pStore[i].dtGrav);
	    if (dT < pkd->pStore[i].dt)
		pkd->pStore[i].dt = dT;

	    }
	}
    }

void pkdFirstDt(PKD pkd) {
    int i,j;
    PARTICLE *p ;
    double a0d2,a1d2;

    p = pkd->pStore;
    for (i=0;i<pkdLocal(pkd);++i) {
	a0d2 = 0.0;
	a1d2 = 0.0;
	if (pkdIsActive(pkd,&p[i]))
	    for (j=0;j<3;++j) {
		a0d2 += p[i].a0[j]*p[i].a0[j];
		a1d2 += p[i].ad0[j]*p[i].ad0[j];
		}
	p[i].dtGrav = (a0d2/a1d2);
#ifdef PLANETS
	if (pkd->param.bCollision) {
	    p[i].iColflag = 0; /* initial reset of collision flag */
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
    pkdKick(pkd,dTime,dDelta,uRungLo,uRungHi);
    }

void pkdKickKDKClose(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    if (pkd->param.csm->bComove) {
	dDelta = csmComoveKickFac(pkd->param.csm,dTime,dDelta);
    }
    pkdKick(pkd,dTime,dDelta,uRungLo,uRungHi);
    }


void pkdKick(PKD pkd,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    int i,j,n;

    pkdClearTimer(pkd,1);
    pkdStartTimer(pkd,1);

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsDstActive(&pkd->pStore[i],uRungLo,uRungHi)) {
	    for (j=0;j<3;++j) {
		pkd->pStore[i].v[j] += pkd->pStore[i].a[j]*dDelta;
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
    int i;

    for (i=0;i<pkdLocal(pkd);++i) {
	if ( !pkdIsDstActive(&pkd->pStore[i],uRungLo,uRungHi) ) continue;
	pkd->pStore[i].uRung = pkd->pStore[i].uNewRung = uRung;
	}
    }

void pkdActiveRung(PKD pkd, int iRung, int bGreater) {
    pkd->uMinRungActive = iRung;
    pkd->uMaxRungActive = bGreater ? 255 : iRung;
    }

int pkdCurrRung(PKD pkd,uint8_t uRung) {
    int i;
    int iCurrent;

    iCurrent = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkd->pStore[i].uRung == uRung) {
	    iCurrent = 1;
	    break;
	    }
	}
    return iCurrent;
    }

void pkdAccelStep(PKD pkd, uint8_t uRungLo,uint8_t uRungHi,
		  double dEta,double dVelFac,double dAccFac,
		  int bDoGravity,int bEpsAcc,int bSqrtPhi,double dhMinOverSoft) {
    int i;
    double vel;
    double acc;
    int j;
    double dT;
    FLOAT fSoft;

    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkdIsDstActive(&(pkd->pStore[i]),uRungLo,uRungHi)) {
	    fSoft = pkdSoft(pkd,&pkd->pStore[i]);
	    vel = 0;
	    acc = 0;
	    for (j=0;j<3;j++) {
		vel += pkd->pStore[i].v[j]*pkd->pStore[i].v[j];
		acc += pkd->pStore[i].a[j]*pkd->pStore[i].a[j];
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
		double dtemp =
		    dEta*3.5*sqrt(dAccFac*fabs(pkd->pStore[i].fPot))/acc;
		if (dtemp < dT)
		    dT = dtemp;
		}
	    pkd->pStore[i].uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
	    }
	}
    }


void pkdDensityStep(PKD pkd, uint8_t uRungLo, uint8_t uRungHi, double dEta, double dRhoFac) {
    int i;
    double dT;

    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkdIsDstActive(&(pkd->pStore[i]),uRungLo,uRungHi)) {
	    dT = dEta/sqrt(pkd->pStore[i].fDensity*dRhoFac);
	    pkd->pStore[i].uNewRung = pkdDtToRung(dT,pkd->param.dDelta,pkd->param.iMaxRung-1);
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
    int i;
    int iTempRung;
    for (i=0;i<iMaxRung;++i) nRungCount[i] = 0;
    for (i=0;i<pkdLocal(pkd);++i) {
	if ( pkdIsDstActive(&pkd->pStore[i],uRungLo,uRungHi) ) {
	    if ( pkd->pStore[i].uNewRung >= iMaxRung )
		pkd->pStore[i].uRung = iMaxRung-1;
	    else if ( pkd->pStore[i].uNewRung >= uRung )
		pkd->pStore[i].uRung = pkd->pStore[i].uNewRung;
	    else if ( pkd->pStore[i].uRung > uRung)
		pkd->pStore[i].uRung = uRung;
	    }
	/*
	** Now produce a count of particles in rungs.
	*/
	nRungCount[pkd->pStore[i].uRung] += 1;
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

    p->iOrder = -2 - p->iOrder;

    p->iClass = getClass(pkd,0.0,0.0); /* Special "DELETED" class */
    }

void pkdNewParticle(PKD pkd, PARTICLE *p) {
    assert(0);
    mdlassert(pkd->mdl,pkd->nLocal < pkd->nStore);
    pkd->pStore[pkd->nLocal] = *p;
    pkd->pStore[pkd->nLocal].iOrder = -1;
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

    assert(0);
    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for (pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if (pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if (p->iOrder == -1) {
	    ++pj;
	    ++nNew;
	    ++ndDark;
	    if (pkdIsActive(pkd,p))
		++pkd->nActive;
	    continue;
	    }
	else if (p->iOrder < -1) {
	    --newnLocal;
	    p->iOrder = -2 - p->iOrder;
	    if (pkdIsGas(pkd, p))
		--ndGas;
	    else if (pkdIsDark(pkd, p))
		--ndDark;
	    else if (pkdIsStar(pkd, p))
		--ndStar;
	    else
		mdlassert(pkd->mdl,0);
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
    int pi;

    assert(0);
    for (pi=0;pi<pkdLocal(pkd);pi++) {
	if (pkd->pStore[pi].iOrder == -1) {
	    pkd->pStore[pi].iOrder = nStart++;
	    }
	}
    }

void pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar,int nMaxOrderGas,
	     int nMaxOrderDark) {
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;
    }


void pkdSetRungVeryActive(PKD pkd, int iRung) {
    /* Remember, the first very active particle is at iRungVeryActive + 1 */
    pkd->uRungVeryActive = iRung;
    }

int pkdIsGas(PKD pkd,PARTICLE *p) {
    if (p->iOrder < pkd->nMaxOrderGas) return 1;
    else return 0;
    }

int pkdIsDark(PKD pkd,PARTICLE *p) {
    if (p->iOrder >= pkd->nMaxOrderGas && p->iOrder < pkd->nMaxOrderDark)
	return 1;
    else return 0;
    }

int pkdIsStar(PKD pkd,PARTICLE *p) {
    if (p->iOrder >= pkd->nMaxOrderDark) return 1;
    else return 0;
    }

#ifdef RELAXATION
void pkdInitRelaxation(PKD pkd) {
    int i,j;

    for (i=0;i<pkdLocal(pkd);++i) {
	pkd->pStore[i].fRelax = 0.0;
	}
    }

#endif /* RELAXATION */

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
	p = &pkd->pStore[i];
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
	p = &pkd->pStore[i];
	p->iOrder = nStart + i;
	if (ssioData(&ssio,&data))
	    mdlassert(pkd->mdl,0); /* error during read in ss file */
	p->iOrgIdx = data.org_idx;
	p->iClass = getClass(data.mass,data.radius);
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
	p = &pkd->pStore[i];
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
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	fMass = pkdMass(pkd,&p[i]);
	if (iFlag == 2) {
	    if (pkdIsActive(pkd,&p[i])) continue;  /* inactive */
	    }
	else if (iFlag == 1) {
	    if (!pkdIsActive(pkd,&p[i])) continue; /* active */
	    }

	r2 = 0;
	rv = 0;
	for (j=0;j<3;++j) {
	    r2 += p[i].r[j]*p[i].r[j];
	    rv += p[i].v[j]*p[i].r[j];
	    }
	r1i = (r2 == 0 ? 0 : 1/sqrt(r2));
	r3i = fMass*r1i*r1i*r1i;
	r5i = 3.0*rv*r3i*r1i*r1i;
	for (j=0;j<3;++j) {
	    aSun[j] += p[i].r[j]*r3i;
	    adSun[j] += p[i].v[j]*r3i-p[i].r[j]*r5i;
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

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&(p[i]))) {
	    r2 = 0;
	    v2 = 0;
	    rv = 0;
	    for (j=0;j<3;++j) {
		r2 += p[i].r[j]*p[i].r[j];
		v2 += p[i].v[j]*p[i].v[j];
		rv += p[i].v[j]*p[i].r[j];
		}

	    r1i = (r2 == 0 ? 0 : 1/sqrt(r2)); /*gravity at origin = zero */
	    p[i].fPot -= dSunMass*r1i;
	    r3i = dSunMass*r1i*r1i*r1i;
	    r5i = 3.0*rv*r3i*r1i*r1i;
	    /* time step is determined by semimajor axis, not the heliocentric distance*/
	    if (pkd->param.bGravStep && p[i].fMass > 0) {
		/* E and h are normalized by the reduced mass */
		sum = dSunMass + p[i].fMass;
		hx =  p[i].r[1]*p[i].v[2] - p[i].r[2]*p[i].v[1];
		hy =  p[i].r[2]*p[i].v[0] - p[i].r[0]*p[i].v[2];
		hz =  p[i].r[0]*p[i].v[1] - p[i].r[1]*p[i].v[0];
		h2 = hx*hx + hy*hy + hz*hz;
		E = 0.5*v2 - sum*r1i;
		e = sqrt(1.0+2.0*E*h2/sum/sum);
		aai = -2.0*E/sum/(1.0-e);
		aai = aai*aai*aai;
		idt2 = sum*aai;
		/*if (p[i].dtSun > p[i].dtGrav) p[i].dtGrav = p[i].dtSun;*/
		if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
		}
	    for (j=0;j<3;++j) {
#ifdef HERMITE
		if (pkd->param.bHermite) {
		    p[i].app[j] = p[i].a[j] - aSun[j]; /* perturbation force*/
		    p[i].adpp[j] = p[i].ad[j] - adSun[j];
		    p[i].ad[j] -= (adSun[j] + p[i].v[j]*r3i-p[i].r[j]*r5i);
		    }
#endif
		p[i].a[j] -= (aSun[j] + p[i].r[j]*r3i);
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
    PARTICLE *p = pkd->pStore;
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
	if ( p[i].iRung <= iRung ) ++i;
	else break;
	}
    while (i <= j) {
	if ( p[j].iRung > iRung ) --j;
	else break;
	}

    if (i < j) {
	t = p[i];
	p[i] = p[j];
	p[j] = t;
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
    PARTICLE *p;
    double r_k = 3.0/pow(2.08,iRung-1);
    double r_kk = r_k/2.08;
    double r_kkk = r_kk/2.08;
    double fourh, dir, dir2;

    assert(iRung >= 1);
    iStart = pkd->nLocal - pkd->nVeryActive;

    p = pkd->pStore;
    n = pkdLocal(pkd);

    /* reset */
    for (i=iStart;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    p[i].drmin = 1000.0;
	    for (k=0;k<3;k++) {
		p[i].a_VA[k] = 0.0;
		}
	    }
	}

    for (i=iStart;i<n;++i) {

	if (pkdIsActive(pkd,&p[i])) {
	    for (k=0;k<p[i].n_VA;k++) {
		j = p[i].i_VA[k];
		if (p[i].iOrder<p[j].iOrder) { /* after three body encounter and a collision
						we sometimes obtain i<j, p[i].iOrder = p[j].iOrder*/
		    x = p[i].r[0] - p[j].r[0];
		    y = p[i].r[1] - p[j].r[1];
		    z = p[i].r[2] - p[j].r[2];
		    d2 = x*x + y*y +z*z;
		    d1 = sqrt(d2);

		    fourh = p[i].fSoft + p[j].fSoft;

		    if (d1 > fourh) {
			dir2 = 1.0/(d1*d2);
			}
		    else {

			if (pkd->param.bCollision) {

			    pkd->iCollisionflag = 1;
			    p[i].iColflag = 1;
			    p[i].iOrderCol = p[j].iOrder;
			    p[i].dtCol = 1.0*p[i].iOrgIdx;
			    printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i piRung %d iRung %d\n",
				   d1,fourh,p[i].iOrgIdx,p[j].iOrgIdx, p[i].iRung, iRung);
			    printf("i %d j %d inVA %d, jnVA %d \n",i, j, p[i].n_VA, p[j].n_VA);
			    }

			dir = 1.0/fourh;
			dir2 = dir*dir;
			d2 *= dir2;
			dir2 *= dir;
			d2 = 1.0 - d2;
			dir *= 1.0 + d2*(0.5 + d2*(3.0/8.0 + d2*(45.0/32.0)));
			dir2 *= 1.0 + d2*(1.5 + d2*(135.0/16.0));
			}


		    d1 /= p[i].hill_VA[k]; /* distance normalized by hill */

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

		    p[i].drmin = (d1 < p[i].drmin)?d1:p[i].drmin;
		    p[j].drmin = (d1 < p[j].drmin)?d1:p[j].drmin;

		    dir2 *= fac;
		    /*printf("iOrder %d %d, d2 %e iRung %d\n",
		      p[i].iOrder,p[j].iOrder,d2,iRung);*/
		    p[i].a_VA[0] -= x*dir2*p[j].fMass;
		    p[i].a_VA[1] -= y*dir2*p[j].fMass;
		    p[i].a_VA[2] -= z*dir2*p[j].fMass;
		    p[j].a_VA[0] += x*dir2*p[i].fMass;
		    p[j].a_VA[1] += y*dir2*p[i].fMass;
		    p[j].a_VA[2] += z*dir2*p[i].fMass;
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
    PARTICLE *p;
    double r_k = 3.0/pow(2.08,iRung);
    int iupdate = 0;

    assert(iRung >= 1);
    p = pkd->pStore;
    iStart = pkd->nLocal - pkd->nVeryActive;
    n = pkdLocal(pkd);

    for (i=iStart;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    for (k=0;k<p[i].n_VA;k++) {
		j = p[i].i_VA[k];
		if (i<j) {
		    x = p[i].r[0] - p[j].r[0];
		    y = p[i].r[1] - p[j].r[1];
		    z = p[i].r[2] - p[j].r[2];
		    d2 = x*x + y*y +z*z;
		    d2 = sqrt(d2);
		    d2 /= p[i].hill_VA[k]; /* distance normalized by hill */

		    /* d2 should be replaced by min.dist. during drift */

		    if (d2 < r_k) {
			iupdate = 1;
			p[i].iRung = iRung +1;
			p[j].iRung = iRung +1;
			nMaxRung = ((iRung+1) >= nMaxRung)?(iRung+1):nMaxRung;
			/* retrive */
			for (k=0;k<3;k++) {
			    p[i].r[k] = p[i].rb[k];
			    p[i].v[k] = p[i].vb[k];
			    p[j].r[k] = p[j].rb[k];
			    p[j].v[k] = p[j].vb[k];
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
		    if (p[i].n_VA >= 2) {
			for (k=0;k<p[i].n_VA;k++) {
			    iTempRung = p[p[i].i_VA[k]].iRung;
			    if (iTempRung > p[i].iRung) {
				if (iTempRung != (iRung +1)) {
				    printf("p_iOrder %d pi_n_VA %d pj_nVA %d iTempRung %d, iRung+1 %d \n",
					   p[i].iOrder,p[i].n_VA,p[p[i].i_VA[k]].n_VA,iTempRung,iRung +1);
				    printf("too many particles in 3 hill radius !! \n");
				    assert(0);
				    }
				/* assert(iTempRung == (iRung +1));*/
				p[i].iRung = iRung + 1;
				for (k=0;k<3;k++) {
				    p[i].r[k] = p[i].rb[k];
				    p[i].v[k] = p[i].vb[k];
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

    p = pkd->pStore;
    for (i=iStart;i<pkdLocal(pkd);++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    for (k=0;k<3;k++) {
		p[i].v[k] += p[i].a_VA[k]*dt;
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
    p = pkd->pStore;
    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    iTempRung = iRung;
	    if (p[i].drmin > 3.0) {
		iTempRung = 0;
		}
	    else {
		iTempRung = floor(log(3.0/p[i].drmin)/log(2.08)) + 1;
		}
	    if (iTempRung >= iMaxRung) {
		iTempRung = iMaxRung-1;
		}

	    /* if min. dist. during drift is less than 3 Hill radius,
	    set iRung = 1 */
	    if (iTempRung == 0 && p[i].drmin2 < 3.0) {
		iTempRung = 1;
		}

	    /* retrive position and velocity of active particle to
	    those before drift*/
	    if (iTempRung > 0) {
		for (j=0;j<3;j++) {
		    p[i].r[j] = p[i].rb[j];
		    p[i].v[j] = p[i].vb[j];
		    }
		}
	    /*
	    ** Now produce a count of particles in rungs.
	    */
	    nRungCount[iTempRung] += 1;
	    p[i].iRung = iTempRung;

	    /*printf("iorder %d, drmin1 %e, drmin2 %e, \n",
	    p[i].iOrder, p[i].drmin, p[i].drmin2);*/
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

    p = pkd->pStore;
    iStart = pkd->nLocal - pkd->nVeryActive;
    n = pkdLocal(pkd);

    if (bRung) {
	for (i=1;i<iMaxRung;++i) nRungCount[i] = 0;
	}

    for (i=iStart;i<n;++i) {
	assert(pkdIsActive(pkd,&p[i]));
	if (bRung) nRungCount[p[i].iRung] += 1;
	n_VA = p[i].n_VA;
	assert(n_VA >= 1);

	if (n_VA >= 2) multiflag = (multiflag < n_VA)?n_VA:multiflag; /* multiflag = largest n_VA */

	for (k=0;k<n_VA;k++) {
	    iOrderTemp = p[i].iOrder_VA[k];

	    for (j=iStart;j<n;++j) {
		if (i==j)continue;
		assert(p[j].iOrder != p[i].iOrder);
		if (p[j].iOrder == iOrderTemp) {
		    p[i].i_VA[k] = j;
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
    PARTICLE *p;
    int nloop;
    /* note iMaxRung = pkd.param->iMaxRung -1 */

    int bRung = 0;
    char c;
    int nRungCount[iMaxRung-1];/* 1 to iMaxRung-1*/
    if (bRung) {
	for (i=1;i<iMaxRung;++i) nRungCount[i] = 0;
	}

    iStart = pkd->nLocal - pkd->nVeryActive;

    p = pkd->pStore;
    for (i=iStart;i<pkdLocal(pkd);++i) {
	if (pkdIsActive(pkd,&p[i])) {
	    /*assert(p[i].n_VA >= 1); n_VA might be 0 after collision */
	    iTempRung = floor(log(3.0/p[i].drmin)/log(2.08)) + 1;
	    /*if(iTempRung >= iMaxRung){
	    printf("iRung %d for particle %d larger than Max iRung %d\n",
	       iTempRung, p[i].iOrder,iMaxRung-1);
	    iTempRung = iMaxRung;
		}*/
	    iTempRung = (iTempRung >= iMaxRung)?(iMaxRung-1):iTempRung;
	    /* iTempRung = (iTempRung >= 0)?iTempRung:0;
	       p[i].iKickRung = iTempRung; */

	    iTempRung = (iTempRung >= iRung)?iTempRung:iRung;
	    p[i].iRung = iTempRung;
	    nMaxRung = (iTempRung >= nMaxRung)?iTempRung:nMaxRung;
	    }
	if (bRung) nRungCount[p[i].iRung] += 1;
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
    return(nMaxRung);
    }


void pkdMomSun(PKD pkd,double momSun[]) {
    PARTICLE *p;
    int i,j,n;

    for (j=0;j<3;++j) {
	momSun[j] = 0.0;
	}
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	for (j=0;j<3;++j) {
	    momSun[j] -= p[i].fMass*p[i].v[j];
	    }
	}
    }

void pkdDriftSun(PKD pkd,double vSun[],double dt) {
    PARTICLE *p;
    int i,j,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	for (j=0;j<3;++j) {
	    p[i].r[j] += vSun[j]*dt;
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
    p = pkd->pStore;
    n = pkdLocal(pkd);

    if (tag_VA) {
	iStart = pkd->nLocal - pkd->nVeryActive;
	}
    else {
	iStart = 0;
	}

    for (i=iStart;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])&&(p[i].iOrder>=0)) {

	    /* copy r and v before drift */
	    for (j=0;j<3;++j) {
		p[i].rb[j] = p[i].r[j];
		p[i].vb[j] = p[i].v[j];
		}

	    /*printf("before: iorder = %d, (x,y,z)= (%e,%e,%e), (vx,vy,vz)= (%e,%e,%e),ms =%e \n",
	    p[i].iOrder, p[i].r[0],p[i].r[1],p[i].r[2],
	    p[i].v[0],p[i].v[1],p[i].v[2], mu);*/

	    iflg = drift_dan(mu,p[i].r,p[i].v,dt); /* see kepler.c */

	    /* printf("exit drift_dan, iflg = %d, (x,y,z)= (%e,%e,%e),(vx,vy,vz)= (%e,%e,%e), ms =%e \n",
	    iflg, p[i].r[0],p[i].r[1],p[i].r[2],p[i].v[0],p[i].v[1],p[i].v[2],mu);*/

	    if (iflg != 0) {
		/*printf("call drift_dan*10, iflg = %d\n",iflg);*/
		dttmp = 0.1*dt;
		for (j=0;j<10;++j) {
		    iflg = drift_dan(mu,p[i].r,p[i].v,dttmp);
		    if (iflg != 0) {
			printf("exit drift_dan, iflg = %d, (x,y,z)= (%e,%e,%e),(vx,vy,vz)= (%e,%e,%e), m = %e, ms =%e \n", iflg, p[i].r[0],p[i].r[1],p[i].r[2],p[i].v[0],p[i].v[1],p[i].v[2],p[i].fMass,mu);
			}
		    assert(iflg == 0);
		    }
		}
	    }
	}
    }

#endif /* SYMBA */
#endif /* PLANETS */

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
    PARTICLE *p;
    int i;
    int n=pkdLocal(pkd);
    for( i=0; i<n; i++ ) pkd->pStore[i].bSrcActive = 1;
    return n;
    }
int pkdSelDstAll(PKD pkd) {
    PARTICLE *p;
    int i;
    int n=pkdLocal(pkd);
    for( i=0; i<n; i++ ) pkd->pStore[i].bDstActive = 1;
    return n;
    }

int pkdSelSrcMass(PKD pkd,double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse ) {
    PARTICLE *p;
    double m;
    int i,n,nSelected;

    n = pkdLocal(pkd);
    nSelected = 0;
    for( i=0; i<n; i++ ) {
	p = &pkd->pStore[i];
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
	p = &pkd->pStore[i];
	m = pkdMass(pkd,p);
	p->bDstActive = isSelected((m >= dMinMass && m <=dMaxMass),setIfTrue,clearIfFalse,p->bDstActive);
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
	p = &pkd->pStore[i];
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
	p = &pkd->pStore[i];
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
	p = &pkd->pStore[i];

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
	p = &pkd->pStore[i];

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
    int iLocal;

    n = pkdLocal(pkd);
    iLocal = 0;
    nChecked = 0;
    for (i=0;i<n;++i) {
	if (pkdIsSrcActive(&pkd->pStore[i],uRungLo,uRungHi)) {
	    nChecked++;
	    if ( pkd->pStore[i].fPot < pkd->pStore[iLocal].fPot )
		iLocal = i;
	    }
	}
    r[0] = pkd->pStore[iLocal].r[0];
    r[1] = pkd->pStore[iLocal].r[1];
    r[2] = pkd->pStore[iLocal].r[2];
    *fPot= pkd->pStore[iLocal].fPot;
    return nChecked;
    }
