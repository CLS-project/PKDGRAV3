#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "mdlbase.h"
#include <stdio.h>
#include <stdint.h>
#include <pthread.h>
#include <assert.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef MDL_FFTW
#include <srfftw_threads.h>
#endif
#ifdef INSTRUMENT
#include "cycle.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __osf__
#define vsnprintf(a,b,c,d) vsprintf((a),(c),(d))
#endif

#define MDL_CACHE_SIZE		2000000
#define MDL_CACHELINE_BITS	3
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)

#define MDL_MBX_RING_SZ	8

#define SRV_STOP		0

typedef int (*mdlWorkFunction)(void *ctx);

typedef struct cacheTag {
    mdlkey_t iKey;
    int nLock;
    int nLast;
    int iLink;
    } CTAG;

typedef struct cacheHeader {
    int cid;
    int mid;
    int id;
    int iLine;
    int iSeq;
    } CAHEAD;

typedef struct mbxStruct {
    pthread_mutex_t mux;
    pthread_cond_t sigReq;
    pthread_cond_t sigRpl;
    pthread_cond_t sigRel;
    int bReq;
    int bRpl;
    int bRel;
    int sid;
    int nBytes;
    char *pszIn;
    char *pszOut;
    } MBX;

typedef struct swxStruct {
    pthread_mutex_t mux;
    pthread_cond_t sigRel;
    pthread_cond_t sigRec;
    pthread_cond_t sigSnd;
    int bRel;
    int bRec;
    size_t nInBytes;
    size_t nOutBufBytes;
    char *pszBuf;
    } SWX;

typedef struct cacheSpace {
    void *pData;
    int iType;
    int iDataSize;
    int nData;
    int iLineSize;
    int nLines;
    int nTrans;
    int iKeyShift;
    int iInvKeyShift;

    mdlkey_t iTransMask;
    mdlkey_t iIdMask;
    int *pTrans;
    CTAG *pTag;
    char *pLine;

    int nCheckOut;

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
    } CACHE;



typedef struct {
    int next;
    int prev;
    void *ctx;
    mdlWorkFunction checkFcn;
    mdlWorkFunction doFcn;
    mdlWorkFunction doneFcn;
    } wqNode;

typedef struct mdlContext {
    mdlBASE base;
    int cacheSize;
    pthread_t *pt;
    struct mdlContext **pmdl;
    /*
     ** Services stuff!
     */
    MBX mbxOwn;
    /*
     ** Swapping Box.
     */
    SWX swxOwn;
    /*
     ** Caching stuff!
     */
    MBX mbxCache[MDL_MBX_RING_SZ];
    pthread_mutex_t muxRing;
    int iRingHd;
    int iRingTl;
    int iRecSeq[128];
    int iSndSeq[128];
    unsigned long uRand;
    int iCaBufSize;
    int nMaxCacheIds;
    int iMaxDataSize;
    CACHE *cache;
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks;
    double dWaiting;
    double dComputing;
    double dSynchronizing;
#endif

    /* The work queue */
    pthread_mutex_t  wqMux;
    int              wqSize;
    int              cudaSize;
    wqNode           *wq;
    wqNode           *wqCUDA;
    int              wqFree; /* Unused queue entries */
    int              wqWait; /* Waiting to be processed */
    int              wqDone; /* Finished entries */
    int              freeCUDA;
    int              busyCUDA;
    } * MDL;

/*
 ** General Functions
 */
int mdlLaunch(int,char **,int (*)(MDL,int,char **),void (*)(MDL));
void mdlFinish(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
		   int,int);
void mdlCommitServices(MDL mdl);
void mdlReqService(MDL,int,int,void *,int);
void mdlGetReply(MDL,int,void *,int *);
void mdlHandler(MDL);

/*
** Collective operations
*/
typedef int MDL_Op;

#define MDL_BAND ((MDL_Op)0x40000001)
#define MDL_BOR ((MDL_Op)0x40000002)
#define MDL_BXOR ((MDL_Op)0x40000003)
#define MDL_LAND ((MDL_Op)0x40000004)
#define MDL_LOR ((MDL_Op)0x40000005)
#define MDL_LXOR ((MDL_Op)0x40000006)
#define MDL_MAX ((MDL_Op)0x40000007)
#define MDL_MAXLOC ((MDL_Op)0x40000008)
#define MDL_MIN ((MDL_Op)0x40000009)
#define MDL_MINLOC ((MDL_Op)0x4000000a)
#define MDL_PROD ((MDL_Op)0x4000000b)
#define MDL_REPLACE ((MDL_Op)0x4000000c)
#define MDL_SUM ((MDL_Op)0x4000000d)

typedef int MDL_Datatype;
#define MDL_FLOAT ((MDL_Datatype)0x50000001)
#define MDL_DOUBLE ((MDL_Datatype)0x50000002)
#define MDL_BYTE ((MDL_Datatype)0x50000003)
#define MDL_INT ((MDL_Datatype)0x50000004)
#define MDL_LONG_LONG ((MDL_Datatype)0x50000005)

int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );
int mdlAllreduce( MDL mdl, void *sendbuf, void *recvbuf, int count,
		  MDL_Datatype datatype, MDL_Op op );
int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype rtype);
int mdlAlltoallv( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype sendtype,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype recvtype);
int mdlAllGather( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype recvtype);
int mdlReduceScatter( MDL mdl, void* sendbuf, void* recvbuf, int *recvcounts,
    MDL_Datatype datatype, MDL_Op op);
int mdlTypeContiguous(MDL mdl,int count, MDL_Datatype old_type, MDL_Datatype *newtype);
int mdlTypeCommit(MDL mdl, MDL_Datatype *datatype );
int mdlTypeFree (MDL mdl, MDL_Datatype *datatype );


/*
** Grid Operations
*/

typedef struct mdlGridContext {
    int n1, n2, n3; /* Real dimensions */
    int a1;         /* Actual size of dimension 1 */
    int s, n;       /* Start and number of slabs */
    int nlocal;     /* Number of local elements */

    uint32_t *rs;  /* Starting slab for each processor */
    uint32_t *rn;  /* Number of slabs on each processor */
    uint32_t *id;  /* Which processor has this slab */
    } * MDLGRID;


/*
** Allocate a MDLGRID context.  This has no actual data, but only describes
** the grid geometry.  The global geometry is set.
*/
void mdlGridInitialize(MDL mdl,MDLGRID *pgrid,int n1,int n2,int n3,int a1);
/*
** Free all memory associated with a MDLGRID context.
*/
void mdlGridFinish(MDL mdl, MDLGRID grid);
/*
** Sets the local geometry (i.e., what is on this processor) of this grid.
*/
void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, int nlocal);
/*
** Share the local geometry between processors.
*/
void mdlGridShare(MDL mdl,MDLGRID grid);
/*
** Allocate the local elements.  The size of a single element is
** given and the local GRID information is consulted to determine
** how many to allocate.
*/
void *mdlGridMalloc(MDL mdl,MDLGRID grid,int nEntrySize);
void mdlGridFree( MDL mdl, MDLGRID grid, void *p );
/*
** This gives the processor on which the given slab can be found.
*/
static inline int mdlGridId(MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    assert(z>=0&&z<grid->n3);
    return grid->id[z];
    }
/*
** This returns the index into the array on the appropriate processor.
*/
static inline int mdlGridIdx(MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    assert(x>=0&&x<grid->a1&&y>=0&&y<grid->n2&&z>=0&&z<grid->n3);
    z -= grid->rs[grid->id[z]]; /* Make "z" zero based for its processor */
    return x + grid->a1*(y + grid->n2*z); /* Local index */
    }

/*
** FFT Operations
*/
#ifdef MDL_FFTW
typedef struct mdlFFTContext {
    MDLGRID rgrid;
    MDLGRID kgrid;
    rfftwnd_plan fplan;
    rfftwnd_plan iplan;
    } * MDLFFT;

size_t mdlFFTInitialize(MDL mdl,MDLFFT *fft,
			int nx,int ny,int nz,int bMeasure);
void mdlFFTFinish( MDL mdl, MDLFFT fft );
fftw_real *mdlFFTMAlloc( MDL mdl, MDLFFT fft );
void mdlFFTFree( MDL mdl, MDLFFT fft, void *p );
void mdlFFT( MDL mdl, MDLFFT fft, fftw_real *data, int bInverse );

/* Grid accessors: r-space */
#define mdlFFTrId(fft,x,y,z) mdlGridId((fft)->rgrid,x,y,z)
#define mdlFFTrIdx(fft,x,y,z) mdlGridIdx((fft)->rgrid,x,y,z)

/* Grid accessors: k-space (note permuted indices) */
#define mdlFFTkId(fft,x,y,z) mdlGridId((fft)->kgrid,x,z,y)
#define mdlFFTkIdx(fft,x,y,z) mdlGridIdx((fft)->kgrid,x,z,y)




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
void *mdlAquire(MDL,int,int,int);
void mdlPrefetch(MDL,int,int,int);
void mdlRelease(MDL,int,void *);
void mdlCacheBarrier(MDL,int);
void mdlCacheCheck(MDL);
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

void mdlSetWorkQueueSize(MDL,int,int);
void mdlAddWork(MDL mdl, void *ctx, mdlWorkFunction initWork, mdlWorkFunction checkWork, mdlWorkFunction doWork, mdlWorkFunction doneWork);

#ifdef __cplusplus
    }
#endif

#endif
