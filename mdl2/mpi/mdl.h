#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include "mpi.h"
#ifdef MDL_FFTW
#include <srfftw_mpi.h>
#endif
#ifdef INSTRUMENT
#include "cycle.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define SRV_STOP		0

#define MDL_CACHE_SIZE		150000000
#define MDL_CACHELINE_BITS	4
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)
#define MDL_CHECK_MASK  	0x7f

#ifndef MDL_MAX_IO_PROCS
#define MDL_MAX_IO_PROCS        128
#endif

/* Maximum number of communicators */
#define MDL_MAX_COMM 10

/*
** A MDL Key must be large enough to hold the largest unique particle key.
** It must be (several times) larger than the total number of particles.
** An "unsigned long" is normally 32 bits on a 32 bit machine and
** 64 bits on a 64 bit machine.  We use uint64_t to be sure.
*/
#ifndef MDLKEY
#define MDLKEY uint64_t
#endif
typedef MDLKEY mdlkey_t;
static const mdlkey_t MDL_INVALID_KEY = (mdlkey_t)(-1);

typedef struct cacheTag {
    mdlkey_t iKey;
    int nLock;
    int iLink;
    } CTAG;


/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct cacheHeader {
    int cid;
    int mid;
    int id;
    int iLine;
    } CAHEAD;


typedef struct cacheSpace {
    int iType;
    void *pData;
    int iDataSize;
    int nData;
    int iLineSize;
    int nLines;
    int iLine;
    int nTrans;
    int iKeyShift;
    int iInvKeyShift;
    int iLastVictim;
    mdlkey_t iTransMask;
    mdlkey_t iIdMask;
    int *pTrans;
    CTAG *pTag;
    char *pLine;
    int nCheckIn;
    int nCheckOut;
    CAHEAD caReq;
    void *ctx;
    void (*init)(void *,void *);
    void (*combine)(void *,void *,void *);
    void * (*getElt)(void *pData,int i,int iDataSize);
    /*
     ** Statistics stuff.
     */
    int nAccess;
    int nAccHigh;
    long nMiss;
    long nColl;
    long nMin;
    int nKeyMax;
    } CACHE;

typedef struct serviceRec {
    int nInBytes;
    int nOutBytes;
    void *p1;
    void (*fcnService)(void *,void *,int,void *,int *);
    } SERVICE;


typedef struct mdlContext {
    int nThreads;
    int nIO;
    int commCount;
    int cacheSize;
    MPI_Comm commMDL;  /* Current active communicator */
    MPI_Comm commList[MDL_MAX_COMM];
    /*MPI_Comm commWork;*/
    /*MPI_Comm commPeer;*/
    int idSelf;
    int bDiag;
    FILE *fpDiag;
    int dontcare;
    int allgrp;
    /*
     ** Services stuff!
     */
    int nMaxServices;
    int nMaxSrvBytes;
    SERVICE *psrv;
    char *pszIn;
    char *pszOut;
    char *pszBuf;
    /*
     ** Swapping buffer.
     */
    char *pszTrans;
    /*
     ** Caching stuff!
     */
    int iMaxDataSize;
    int iCaBufSize;
    char *pszRcv;
    int *pmidRpl;
    MPI_Request *pReqRpl;
    MPI_Request ReqRcv;
    char **ppszRpl;
    char *pszFlsh;
    int nMaxCacheIds;
    CACHE *cache;
    char nodeName[MPI_MAX_PROCESSOR_NAME];
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks;
    double dWaiting;
    double dComputing;
    double dSynchronizing;
#endif
    } * MDL;


/*
 * MDL debug and Timer macros and prototypes
 */
/*
 * Compile time mdl debugging options
 *
 * mdl asserts: define MDLASSERT
 * Probably should always be on unless you want no mdlDiag output at all
 *
 * NB: defining NDEBUG turns off all asserts so MDLASSERT will not assert
 * however it will output using mdlDiag and the code continues.
 */
#define MDLASSERT
/*
 * Debug functions active: define MDLDEBUG
 * Adds debugging mdldebug prints and mdldebugassert asserts
 */
#define MDLDEBUG
/*
 * Timer functions active: define MDLTIMER
 * Makes mdl timer functions active
 */
#define MDLTIMER


void mdlprintf( MDL mdl, const char *format, ... );

#ifdef MDLASSERT
#ifdef __ANSI_CPP__
#define mdlassert(mdl,expr) \
    { \
      if (!(expr)) { \
	     mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, # expr ); \
	     assert( expr ); \
	     } \
    }
#else
#define mdlassert(mdl,expr) \
    { \
      if (!(expr)) { \
	     mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, "expr" ); \
	     assert( expr ); \
	     } \
    }
#endif
#else
#define mdlassert(mdl,expr)  assert(expr)
#endif

#ifdef MDLDEBUG
#define mdldebugassert(mdl,expr)   mdlassert(mdl,expr)
void mdldebug( MDL mdl, const char *format, ... );
#else
#define mdldebug
#define mdldebugassert
#endif

typedef struct {
    double wallclock;
    double cpu;
    double system;
    } mdlTimer;

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl,mdlTimer *);
void mdlGetTimer(MDL mdl,mdlTimer *,mdlTimer *);
void mdlPrintTimer(MDL mdl,char *message,mdlTimer *);
#else
#define mdlZeroTimer
#define mdlGetTimer
#define mdlPrintTimer
#endif

/*
 ** General Functions
 */
double mdlCpuTimer(MDL);
int mdlInitialize(MDL *,char **,void (*)(MDL),void (*)(MDL));
void mdlFinish(MDL);
int  mdlSplitComm(MDL mdl, int nProcs);
void mdlSetComm(MDL mdl, int iComm);
void mdlStop(MDL);
int mdlThreads(MDL);
int mdlIO(MDL);
int mdlSelf(MDL);
int mdlOldSelf(MDL);
const char *mdlName(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
typedef int (*mdlPack)(void *,int *,size_t,void*);
void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx);
void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx);
void mdlDiag(MDL,char *);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
		   int,int);
void mdlReqService(MDL,int,int,void *,int);
void mdlGetReply(MDL,int,void *,int *);
void mdlHandler(MDL);

/*
** Collective operations
*/
#define MDL_BAND  MPI_BAND
#define MDL_BOR MPI_BOR
#define MDL_BXOR MPI_BXOR
#define MDL_LAND MPI_LAND
#define MDL_LOR MPI_LOR
#define MDL_LXOR MPI_LXOR
#define MDL_MAX MPI_MAX
#define MDL_MAXLOC MPI_MAXLOC
#define MDL_MIN MPI_MIN
#define MDL_MINLOC MPI_MINLOC
#define MDL_PROD MPI_PROD
#define MDL_REPLACE MPI_REPLACE
#define MDL_SUM MPI_SUM

#define MDL_BYTE MPI_BYTE
#define MDL_INT MPI_INT
#define MDL_LONG MPI_LONG
#define MDL_LONG_LONG_INT MPI_LONG_LONG_INT
#define MDL_SHORT MPI_SHORT
#define MDL_UNSIGNED MPI_UNSIGNED
#define MDL_UNSIGNED_LONG MPI_UNSIGNED_LONG
#define MDL_UNSIGNED_LONG_LONG MPI_UNSIGNED_LONG_LONG
#define MDL_UNSIGNED_SHORT MPI_UNSIGNED_SHORT
#define MDL_DOUBLE_INT MPI_DOUBLE_INT
#define MDL_FLOAT_INT MPI_FLOAT_INT
#define MDL_LONG_INT MPI_LONG_INT
#define MDL_LONG_DOUBLE_INT MPI_LONG_DOUBLE_INT
#define MDL_SHORT_INT MPI_SHORT_INT
#define MDL_2INT MPI_2INT
#define MDL_COMPLEX  MPI_COMPLEX
#define MDL_DOUBLE MPI_DOUBLE
#define MDL_DOUBLE_PRECISION MPI_DOUBLE_PRECISION
#define MDL_FLOAT MPI_FLOAT
#define MDL_LONG_DOUBLE MPI_LONG_DOUBLE
#define MDL_REAL MPI_REAL
#define MDL_LOGICAL  MPI_LOGICAL

typedef MPI_Op MDL_Op;
typedef MPI_Datatype MDL_Datatype;
#define MDL_Op_create(f,c,o) MPI_Op_create(f,c,o) 
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );
int mdlAllreduce( MDL mdl, void *sendbuf, void *recvbuf, int count,
		  MDL_Datatype datatype, MDL_Op op );
int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
		 void *recvbuf, int rcount, MDL_Datatype rtype);

/*
** FFT Operations
*/
#ifdef MDL_FFTW
typedef struct mdlFFTContext {
    MDL mdl;
    rfftwnd_mpi_plan fplan;
    rfftwnd_mpi_plan iplan;

    int n1, n2, n3; /* Real dimensions */
    int a1r, a1k;   /* Actual size in r and k space */
    int sz, nz;     /* Start z and number of z */
    int sy, ny;     /* Transposed start and number */
    int nlocal;     /* Number of local elements */

    uint32_t *rsz;  /* Starting z for each processor */
    uint32_t *rsy;
    uint32_t *rnz;  /* Number of z slabs on each processor */
    uint32_t *rny;

    uint32_t *zid;  /* Which processor has this z */
    uint32_t *yid;

    } * MDLFFT;

size_t mdlFFTInitialize(MDL mdl,MDLFFT *fft,
			int nx,int ny,int nz,int bMeasure);
void mdlFFTFinish( MDLFFT fft );
fftw_real *mdlFFTMAlloc( MDLFFT fft );
void mdlFFTFree( MDLFFT fft, void *p );
void mdlFFT( MDLFFT fft, fftw_real *data, int bInverse );

/*
** This gives the processor on which the given "z" slab can
** be found when in r-space (z is the major dimension).
*/
static inline int mdlFFTrId(MDLFFT fft, uint32_t x, uint32_t y, uint32_t z) {
    return fft->zid[z];
    }

/*
** This gives the processor on which the given "y" slab can
** be found when in k-space (y is the major dimension).
*/
static inline int mdlFFTkId(MDLFFT fft, uint32_t x, uint32_t y, uint32_t z) {
    return fft->yid[y];
    }

/*
** This returns the index into the array on the appropriate processor
** for the specified cell when in r-space.
*/
static inline int mdlFFTrIdx(MDLFFT fft, uint32_t x, uint32_t y, uint32_t z) {
    z -= fft->rsz[fft->zid[z]]; /* Make "z" zero based for its processor */
    return x + fft->a1r*(y + fft->n2*z); /* Local index */
    }

/*
** This returns the index into the array on the appropriate processor
** for the specified cell when in k-space.
*/
static inline int mdlFFTkIdx(MDLFFT fft, uint32_t x, uint32_t y, uint32_t z) {
    y -= fft->rsy[fft->yid[y]]; /* Make "y" zero based for its processor */
    return x + fft->a1k*(z + fft->n3*y); /* Local index */
    }
#endif
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
void mdlSetCacheSize(MDL,int);
void mdlROcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData);
void mdlCOcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData,
		void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *));
void mdlFinishCache(MDL,int);
void mdlCacheCheck(MDL);
void mdlCacheBarrier(MDL,int);

void *mdlDoMiss(MDL mdl, int cid, int iIndex, int id, mdlkey_t iKey, int lock);

static inline void *mdlAquire(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 1;  /* we always lock in aquire */
    CACHE *c = &mdl->cache[cid];
    char *pLine;
    mdlkey_t iKey;
    int iElt;
    int i;

    ++c->nAccess;
    if (!(c->nAccess & MDL_CHECK_MASK))
	mdlCacheCheck(mdl);
    /*
     ** Determine memory block key value and cache line.
     */
    iKey = ((iIndex&MDL_INDEX_MASK) << c->iKeyShift)| id;
    i = c->pTrans[iKey & c->iTransMask];
    /*
     ** Check for a match!
     */
    if (c->pTag[i].iKey == iKey) {
	++c->pTag[i].nLock;
	pLine = &c->pLine[i*c->iLineSize];
	iElt = iIndex & MDL_CACHE_MASK;
	return(&pLine[iElt*c->iDataSize]);
	}
    i = c->pTag[i].iLink;
    /*
     ** Collision chain search.
     */
    while (i) {
	++c->nColl;
	if (c->pTag[i].iKey == iKey) {
	    ++c->pTag[i].nLock;
	    pLine = &c->pLine[i*c->iLineSize];
	    iElt = iIndex & MDL_CACHE_MASK;
	    return(&pLine[iElt*c->iDataSize]);
	    }
	i = c->pTag[i].iLink;
	}
    /*
     ** Is it a local request? This should not happen in pkdgrav2.
     */
    if (id == mdl->idSelf) {
	return (*c->getElt)(c->pData,iIndex,c->iDataSize);
	}
    return(mdlDoMiss(mdl, cid, iIndex, id, iKey, lock));
    }


void mdlRelease(MDL,int,void *);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);
double mdlMinRatio(MDL,int);

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
void mdlTimeReset(MDL);
double mdlTimeComputing(MDL);
double mdlTimeSynchronizing(MDL);
double mdlTimeWaiting(MDL);
#endif

#ifdef __cplusplus
}
#endif
#endif
