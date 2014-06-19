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
#include <srfftw_mpi.h>
#endif
#include "opa_queue.h"
/*#define USE_ARC*/
#ifdef USE_ARC
typedef struct CacheDataBucket {
    uint32_t uId;       /* upper 4 bits encode ARC_where and dirty bit */
    uint32_t uIndex;    /* page's ID number */
    struct CacheDataBucket *coll; /* collision chain for hash table */
    struct CacheDataBucket *next;      /* for doubly linked list */
    struct CacheDataBucket *prev;      /* for doubly linked list */
    uint64_t *data;      /* page's location in cache */
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
} * ARC;
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef MPI_VERSION
#define MPI_VERSION 1
#endif

#define SRV_STOP		0

#define MDL_CACHE_SIZE		150000000
#define MDL_CACHELINE_BITS	4
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)
#define MDL_CHECK_MASK  	0x7f

/* Maximum number of communicators */
#define MDL_MAX_COMM 10


typedef struct {
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    } MDLserviceElement;

typedef struct {
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    void *buf;
    int count;
    int target;
    int tag;
    MPI_Datatype datatype;
    } MDLserviceSend;

typedef int (*mdlWorkFunction)(void *ctx);

typedef struct {
    OPA_Queue_element_hdr_t hdr;
    int iCoreOwner;
    void *ctx;
    mdlWorkFunction checkFcn;
    mdlWorkFunction doFcn;
    mdlWorkFunction doneFcn;
    } MDLwqNode;

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
    int16_t cid;
    int16_t mid;
    int32_t idFrom;
    int32_t idTo;
    int32_t iLine;
    } CAHEAD;

typedef struct {
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    void *pLine;
    CAHEAD caReq;
    } MDLserviceCacheReq;

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
    int iLastVictim;
    mdlkey_t iTransMask;
    mdlkey_t iIdMask;
#ifdef USE_ARC
    ARC arc;
    void *pOneLine;
#endif
    int *pTrans;
    CTAG *pTag;
    char *pLine;
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
    int nKeyMax;
    } CACHE;

typedef struct {
    MPI_Comm commMDL;             /* Current active communicator */
    OPA_Queue_info_t queueMPI;
    int *pRequestTargets;
    int nRequestTargets;
    int iRequestTarget;
    int nSendRecvReq;

    int nActiveCores;
    int nOpenCaches;
    int iMaxDataSize; /* Largest cache element */
    int iCaBufSize;  /* Cache buffer size */
    char *pszRcv; /* Cache receive buffer */

    int *pmidRpl;
    MPI_Request *pReqRpl;

    MPI_Request *pSendRecvReq;
    MPI_Request ReqRcv;
    MDLserviceSend **pSendRecvBuf;
    MDLserviceCacheReq **pThreadCacheReq;
    char **ppszRpl;
    } mdlContextMPI;

typedef struct mdlContext {
    mdlBASE base;
    struct mdlContext **pmdl;
    pthread_t *threadid;
    void * (*fcnWorker)(struct mdlContext *mdl);

    OPA_Queue_info_t *inQueue;
    MDLserviceElement inMessage;
    int iCoreMPI;             /* Core that handles MPI requests */
    int cacheSize;

    /* Work Queues */
    OPA_Queue_info_t wq;     /* Work for us to do */
    OPA_Queue_info_t wqDone; /* Completed work from other threads */
    OPA_Queue_info_t wqFree; /* Free work queue nodes */
    int wqMaxSize;
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
    char *pszFlsh;
    int nMaxCacheIds;
    CACHE *cache;

    mdlContextMPI mpi;
    MDLserviceSend sendRequest;
    MDLserviceSend recvRequest;
    } * MDL;


/*
 ** General Functions
 */
int mdlLaunch(int,char **,void * (*)(MDL),void * (*)(MDL));
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
int mdlBcast ( MDL mdl, void *buf, int count, MDL_Datatype datatype, int root );
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );
int mdlScan   ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op );
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
    rfftwnd_mpi_plan fplan;
    rfftwnd_mpi_plan iplan;
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
void mdlCacheCheck(MDL);
void mdlCacheBarrier(MDL,int);
void mdlPrefetch(MDL mdl,int cid,int iIndex, int id);
void *mdlAquire(MDL mdl,int cid,int iIndex,int id);
void mdlRelease(MDL,int,void *);
void mdlFlushCache(MDL,int);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);

void mdlSetWorkQueueSize(MDL,int,int);
void mdlAddWork(MDL mdl, void *ctx, mdlWorkFunction initWork, mdlWorkFunction checkWork, mdlWorkFunction doWork, mdlWorkFunction doneWork);

#ifdef __cplusplus
}
#endif
#endif
