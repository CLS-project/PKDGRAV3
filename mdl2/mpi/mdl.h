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
#ifdef INSTRUMENT
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
typedef MPI_Op MDL_Op;
typedef MPI_Datatype MDL_Datatype;
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );

/*
** FFT Operations
*/
#ifdef MDL_FFTW
typedef struct mdlFFTContext {
    MDL mdl;
    rfftwnd_mpi_plan fplan;
    rfftwnd_mpi_plan iplan;

    int rx, ry, rz; /* Real dimensions */
    int sz, nz;     /* Start z and number of z */
    int sy, ny;     /* Transposed start and number */
    int nlocal;     /* Number of local elements */
    } * MDLFFT;

size_t mdlFFTInitialize(MDL mdl,MDLFFT *fft,
			int nx,int ny,int nz,int bMeasure);
void mdlFFTFinish( MDLFFT fft );
void mdlFFT( MDLFFT fft, fftw_real *data, int bInverse );
#endif
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
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

#ifdef INSTRUMENT
void mdlTimeReset(MDL);
double mdlTimeComputing(MDL);
double mdlTimeSynchronizing(MDL);
double mdlTimeWaiting(MDL);
#endif

#endif
