#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef INSTRUMENT
#include "cycle.h"
#endif
#include <pthread.h>
#include "mdlbase.h"
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#include <stdint.h>
#include "mpi.h"
#ifdef MDL_FFTW
#include <fftw3-mpi.h>
#ifdef USE_SINGLE
#define FFTW3(name) fftwf_ ## name
typedef float fftwf_real;
#else
#define FFTW3(name) fftw_ ## name
typedef double fftw_real;
#endif
#endif
#include "opa_queue.h"

typedef struct {
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    } MDLserviceElement;

typedef struct  {
    struct CacheDataBucket *next;      /* for doubly linked list */
    struct CacheDataBucket *prev;      /* for doubly linked list */
    } CacheDataLinks;

typedef struct CacheDataBucket {
    /*
    ** When we are flushing, the element has been removed from all lists,
    ** and from the hash table. We can "send" the CDB to the MPI node and
    ** get it returned to us later.
    */
    union {
	CacheDataLinks  links; /* When in the ARC cache */
	MDLserviceElement svc; /* When sending to/from the MPI thread */
	} hdr;
    union {
	struct CacheDataBucket *coll; /* collision chain for hash table */
	struct cacheSpace *cache; /* when flushing need to know which cache */
	} extra;
    uint64_t *data;      /* page's location in cache */
    uint32_t uId;       /* upper 4 bits encode ARC_where and dirty bit */
    uint32_t uIndex;    /* page's ID number */
} CDB;


typedef struct ArcContext {
    CDB **Hash;
    CDB *cdbBase;
    uint64_t *dataBase;
    uint64_t *dataLast;
    CDB *T1;
    CDB *B1;
    CDB *T2;
    CDB *B2;
    CDB *Free;
    uint32_t nHash;
    uint32_t uHashMask;
    uint32_t nCache;
    uint32_t uDataSize;
    uint32_t T1Length;
    uint32_t B1Length;
    uint32_t T2Length;
    uint32_t B2Length;
    uint32_t target_T1;
    struct cacheSpace *cache;
} * ARC;

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MPI_VERSION
#define MPI_VERSION 1
#endif

#define SRV_STOP		0

#define MDL_CACHE_SIZE		15000000
#define MDL_CACHELINE_ELTS	(32)
#define MDL_CHECK_MASK  	0x7f

typedef struct {
    MDLserviceElement svc;
    void *buf;
    int count;
    int target;
    int tag;
    MPI_Datatype datatype;
    } MDLserviceSend;

typedef int (*mdlWorkFunction)(void *ctx);

typedef struct mdl_wq_node {
    /* We can put this on different types of queues */
    union {
	OPA_Queue_element_hdr_t hdr;
	} q;
    int iCoreOwner;
    void *ctx;
    mdlWorkFunction doFcn;
    mdlWorkFunction doneFcn;
#ifdef USE_CUDA
    double *pHostBuf;
    double *pCudaBuf;
#endif
    } MDLwqNode;

typedef struct cacheTag {
    mdlkey_t iKey;
    int nLock;
    int iLink;
    } CTAG;


/*
** This structure should be "maximally" aligned.
*/
typedef struct cacheHeader {
    uint8_t cid;
    uint8_t mid;
    uint16_t nItems;
    int32_t idFrom;
    int32_t idTo;
    int32_t iIndex;
    } CAHEAD;

typedef struct cacheHeaderNew {
    uint8_t mid;     /*  2: Message ID: request, flush, reply */
    uint8_t cid;     /*  6: Cache for this entry */
    union {
	uint8_t tidFrom; /* 12: thread id that made the request */
	uint8_t xxx;
	};
    uint8_t tidTo;   /* 12: thread id that gets the reply */
    uint32_t iIndex; /* Index of the element */
    } CAHEADnew;

typedef struct {
    struct mdl_flush_buffer *next, *prev;
    MPI_Request request;
    } mdl_flush_links;

typedef struct mdl_flush_buffer {
    union {
	mdl_flush_links mpi;
	OPA_Queue_element_hdr_t hdr;
	MDLserviceElement svc; /* When sending to/from the MPI thread */
	} hdr;
    MPI_Request request;
    uint32_t nBufferSize;
    uint32_t nBytes;
    uint32_t iRankTo;
    } MDLflushBuffer;
/* followed by CAHEAD, element, CAHEAD, element, etc. */

#define MDL_CACHE_DATA_SIZE (8000)
typedef struct cache_reply_data {
    union {
	struct cache_reply_data *next;
	};
    MPI_Request mpiRequest;
    int nBytes;
    } MDLcacheReplyData;
    /* followed by CAHEAD and data */

typedef struct {
    MDLserviceElement svc;
    void *pLine;
    CAHEAD caReq;
    MPI_Request request;
    } MDLserviceCacheReq;

typedef struct cacheSpace {
    void *pData;
    uint16_t iType;
    uint16_t iCID;
    int iDataSize;
    int nData;
    int nLineElements;
    int iLineSize;
    ARC arc;
    char *pOneLine;

    MDLserviceCacheReq cacheRequest;
    void *ctx;
    void (*init)(void *,void *);
    void (*combine)(void *,void *,void *);
    void * (*getElt)(void *pData,int i,int iDataSize);
    /*
     ** Statistics stuff.
     */
    uint64_t nAccess;
    uint64_t nMiss;
    uint64_t nColl;
    } CACHE;

typedef struct {
    MPI_Comm commMDL;             /* Current active communicator */
    OPA_Queue_info_t queueMPI;
    OPA_Queue_info_t localFlushBuffers;
    MDLflushBuffer **flushBuffersByRank;
    MDLflushBuffer **flushBuffersByCore;
    MDLflushBuffer flushHeadBusy;
    MDLflushBuffer flushHeadFree;
    MDLflushBuffer flushHeadSent;
    int flushBusyCount, flushBuffCount;
    int *pRequestTargets;
    int nRequestTargets;
    int iRequestTarget;
    int nSendRecvReq;
    int nActiveCores;
    int nOpenCaches;
    int iCacheBufSize;  /* Cache input buffer size */
    int iReplyBufSize;  /* Cache reply buffer size */
    MPI_Request *pSendRecvReq;
    MDLcacheReplyData *pReqRcv;
    MDLserviceSend **pSendRecvBuf;
    MDLserviceCacheReq **pThreadCacheReq;
    MDLcacheReplyData *freeCacheReplies;
    MDLcacheReplyData *busyCacheReplies, **busyCacheRepliesTail;
    } mdlContextMPI;

typedef struct mdlContext {
    mdlBASE base;
    struct mdlContext **pmdl;
    pthread_t *threadid;
    pthread_barrier_t barrier;
    void * (*fcnWorker)(struct mdlContext *mdl);

    OPA_Queue_info_t *inQueue;
    MDLserviceElement inMessage;
    void *pvMessageData; /* These two are for the collective malloc */
    size_t nMessageData;
    int iCoreMPI;             /* Core that handles MPI requests */
    int cacheSize;

    /* Work Queues */
    OPA_Queue_info_t wq;     /* Work for us to do */
    OPA_Queue_info_t wqDone; /* Completed work from other threads */
    OPA_Queue_info_t wqFree; /* Free work queue nodes */
    int wqMaxSize;
    uint16_t wqAccepting;
    uint16_t wqLastHelper;
    OPA_int_t wqCurSize;

    /*
     ** Services stuff!
     */
    int nMaxSrvBytes;
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
    int nMaxCacheIds;
    CACHE *cache;
    MDLflushBuffer *coreFlushBuffer;
    OPA_Queue_info_t coreFlushBuffers; /* Buffers we can destage to */
    OPA_Queue_info_t wqCacheFlush;

    int nFlushOutBytes;

    mdlContextMPI *mpi;
    MDLserviceSend sendRequest;
    MDLserviceSend recvRequest;

#ifdef USE_CUDA
    void *cudaCtx;
    int inCudaBufSize, outCudaBufSize;
#endif
    } * MDL;


/*
 ** General Functions
 */
void mdlLaunch(int,char **,void * (*)(MDL),void * (*)(MDL));
void mdlFinish(MDL);
int  mdlSplitComm(MDL mdl, int nProcs);
void mdlSetComm(MDL mdl, int iComm);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
typedef int (*mdlPack)(void *,int *,size_t,void*);
void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx);
void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
		   int,int);
void mdlCommitServices(MDL mdl);
int  mdlReqService(MDL, int, int, void *, int);
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
#define MDL_LONG_LONG MPI_LONG_LONG
#define MDL_FLOAT MPI_FLOAT
#define MDL_DOUBLE MPI_DOUBLE
#define MDL_DATATYPE_NULL MPI_DATATYPE_NULL

typedef MPI_Op MDL_Op;
typedef MPI_Datatype MDL_Datatype;
#define MDL_Op_create(f,c,o) MPI_Op_create(f,c,o) 
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );
int mdlExscan   ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op );
int mdlAllreduce( MDL mdl, void *sendbuf, void *recvbuf, int count,
		  MDL_Datatype datatype, MDL_Op op );
int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype rtype);
int mdlAlltoallv( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype sendtype,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype recvtype);
int mdlAlltoallw( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype *stypes,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype *rtypes);
int mdlAllGather( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype recvtype);
int mdlAllGatherv( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype recvtype);
int mdlReduceScatter( MDL mdl, void* sendbuf, void* recvbuf, int *recvcounts,
    MDL_Datatype datatype, MDL_Op op);
int mdlTypeContiguous(MDL mdl,int count, MDL_Datatype old_type, MDL_Datatype *newtype);
int mdlTypeIndexed(MDL mdl, int count,
    int array_of_blocklengths[], int array_of_displacements[],
    MDL_Datatype oldtype, MDL_Datatype *newtype);
int mdlTypeCommit(MDL mdl, MDL_Datatype *datatype );
int mdlTypeFree (MDL mdl, MDL_Datatype *datatype );

/*
** Grid Operations
*/

typedef struct mdlGridContext {
    uint32_t n1,n2,n3;     /* Real dimensions */
    uint32_t a1;           /* Actual size of dimension 1 */
    uint32_t sSlab, nSlab; /* Start and number of slabs */
    uint64_t nLocal;       /* Number of local elements */
    uint32_t *rs;  /* Starting slab for each processor */
    uint32_t *rn;  /* Number of slabs on each processor */
    uint32_t *id;  /* Which processor has this slab */
    } * MDLGRID;

typedef struct {
    MDLGRID grid;
    uint64_t I;
    int x, y, z; /* Real coordinate */
    int i;       /* Index into local array */
    } mdlGridCoord;

static inline int mdlGridCoordCompare(const mdlGridCoord *a,const mdlGridCoord *b) {
    return a->x==b->x && a->y==b->y && a->z==b->z; 
    }

static inline mdlGridCoord *mdlGridCoordIncrement(mdlGridCoord *a) {
    ++a->i;
    ++a->I;
    if ( ++a->x == a->grid->n1 ) {
	a->i += a->grid->a1 - a->grid->n1;
	a->x = 0;
	if ( ++a->y == a->grid->n2 ) {
	    a->y = 0;
	    ++a->z;
	    }
	}
    return a;
    }

void mdlGridCoordFirstLast(MDL mdl,MDLGRID grid,mdlGridCoord *f,mdlGridCoord *l);

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
void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, uint64_t nLocal);
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
static inline int mdlGridId(MDL mdl,MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    assert(z<grid->n3);
    return mdlProcToThread(mdl,grid->id[z]);
    }
/*
** This returns the index into the array on the appropriate processor.
*/
static inline int mdlGridIdx(MDL mdl,MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    assert(x<=grid->a1 && y<grid->n2 && z<grid->n3);
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
    FFTW3(plan) fplan, iplan;
    } * MDLFFT;

size_t mdlFFTlocalCount(MDL mdl,int n1,int n2,int n3,int *nz,int *sz,int *ny,int*sy);
size_t mdlFFTInitialize(MDL mdl,MDLFFT *pfft,int nx,int ny,int nz,int bMeasure,FFTW3(real) *data);
void mdlFFTFinish( MDL mdl, MDLFFT fft );
FFTW3(real) *mdlFFTMalloc( MDL mdl, MDLFFT fft );
void mdlFFTFree( MDL mdl, MDLFFT fft, void *p );
void mdlFFT( MDL mdl, MDLFFT fft, FFTW3(real) *data);
void mdlIFFT( MDL mdl, MDLFFT fft, FFTW3(complex) *data);

/* Grid accessors: r-space */
#define mdlFFTrId(mdl,fft,x,y,z) mdlGridId(mdl,(fft)->rgrid,x,y,z)
#define mdlFFTrIdx(mdl,fft,x,y,z) mdlGridIdx(mdl,(fft)->rgrid,x,y,z)

/* Grid accessors: k-space (note permuted indices) */
#define mdlFFTkId(mdl,fft,x,y,z) mdlGridId(mdl,(fft)->kgrid,x,z,y)
#define mdlFFTkIdx(mdl,fft,x,y,z) mdlGridIdx(mdl,(fft)->kgrid,x,z,y)

#endif
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
void *mdlMallocArray(MDL mdl,size_t nmemb,size_t size,size_t minSize);
void *mdlSetArray(MDL mdl,size_t nmemb,size_t size,void *vdata);
void mdlFreeArray(MDL,void *);
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
void mdlPrefetch(MDL mdl,int cid,int iIndex, int id);
void *mdlAcquire(MDL mdl,int cid,int iIndex,int id);
void *mdlFetch(MDL mdl,int cid,int iIndex,int id);
void *mdlVirtualAcquire(MDL mdl,int cid,int iIndex,int id,int bLock);
void mdlRelease(MDL,int,void *);
void mdlFlushCache(MDL,int);
void mdlThreadBarrier(MDL);
void mdlCompleteAllWork(MDL);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);

void mdlSetWorkQueueSize(MDL,int,int);
void mdlSetCudaBufferSize(MDL,int,int);
void mdlAddWork(MDL mdl, void *ctx,
    int (*initWork)(void *ctx,void *vwork),
    int (*checkWork)(void *ctx,void *vwork),
    mdlWorkFunction doWork,
    mdlWorkFunction doneWork);

#ifdef __cplusplus
}
#endif
#endif
