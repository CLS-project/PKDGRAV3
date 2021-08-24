/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#include "mdlbase.h"
#include <pthread.h>
#ifdef __APPLE__
#include "pthread_barrier.h"
#endif
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#include <stdint.h>
#include "mpi.h"
#include "mdlfft.h"
#ifdef __cplusplus
#include "arc.h"
#include "mdlmessages.h"
#ifdef USE_CUDA
#include "mdlcuda.h"
#endif
#include "rwlock.h"
#include <memory>
#include <tuple>
#include <vector>
#include <list>
#include <map>
#endif

#ifndef MPI_VERSION
#define MPI_VERSION 1
#endif

#define SRV_STOP		0

/*                                   MQ   MPI */
#define MDL_TAG_BARRIER        	1 /* Yes  Yes */
#define MDL_TAG_SWAPINIT 	2 /* NO   Yes */
#define MDL_TAG_CACHE_FLUSH 	2 /* Yes  NO  */
#define MDL_TAG_SWAP		3 /* NO   Yes */
#define MDL_TAG_REQ	   	4 /* NO   Yes */
#define MDL_TAG_RPL		5 /* This is treated specially */
#define MDL_TAG_SEND            6 /* NO   Yes */
#define MDL_TAG_CACHECOM	7 /* Yes  Yes */

#define MDL_TAG_MAX             8

typedef int (*mdlWorkFunction)(void *ctx);
typedef int (*mdlPack)(void *,int *,size_t,void*);

typedef void * MDL;

#ifdef __cplusplus
namespace mdl {
    static const auto mdl_cache_size = 15'000'000;
    static const auto mdl_check_mask = 0x7f;

typedef struct mdl_wq_node {
    /* We can put this on different types of queues */
    union {
	OPA_Queue_element_hdr_t hdr;
	};
    int iCoreOwner;
    void *ctx;
    mdlWorkFunction doFcn;
    mdlWorkFunction doneFcn;
    } MDLwqNode;

// typedef struct cacheHeaderNew {
//     uint32_t mid     :  2; /*  2: Message ID: request, flush, reply */
//     uint32_t cid     :  6; /*  6: Cache for this entry */
//     uint32_t tidFrom : 12; /* 12: thread id that made the request */
//     uint32_t tidTo   : 12; /* 12: thread id that gets the reply */
//     uint32_t iLine;        /* Index of the cache line */
//     } CacheHeadernew;

#define MDL_CACHE_DATA_SIZE (512)

//*****************************************************************************

// This tell the cach how to marshall data
//     pack - pack data from a table to an outgoing buffer
//   unpack - expand data from a buffer to a cache element
//     init - For a combiner cache perform the initialization of the element
//            This is also called for a virtual fetch (pack/unpack skipped)
//    flush - pack data from the cache into a buffer for flushing
//  combine - expand data from the buffer into the original table
class CACHEhelper {
    friend class mdlClass;
    friend class mpiClass;
    friend class CACHE;
protected:
    uint32_t nData;   // Size of a data element in bytes
    bool     bModify; // If this cache will be modified (and thus flushed)
    virtual void    pack(void *dst, const void *src) { if (src) memcpy(dst,src,nData); else memset(dst,0,nData); }
    virtual void  unpack(void *dst, const void *src, const void *key) { memcpy(dst,src,nData); }
    virtual void    init(void *dst) {}
    virtual void   flush(void *dst, const void *src)                  { memcpy(dst,src,nData); }
    virtual void combine(void *dst, const void *src, const void *key) { memcpy(dst,src,nData); }
    virtual void *create(uint32_t size, const void *pKey) { return nullptr; }
    virtual uint32_t getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) {return uId;}
    virtual bool  modify() {return bModify; }
    virtual uint32_t pack_size() {return nData;}
    virtual uint32_t flush_size() {return nData;}
public:
    explicit CACHEhelper(uint32_t nData, bool bModify=false) : nData(nData), bModify(bModify) {}
    virtual ~CACHEhelper() = default;
    };

class CACHE : private ARChelper {
    friend class mdlClass;
    friend class mpiClass;
    std::shared_ptr<CACHEhelper> cache_helper;
    std::unique_ptr<GARC> arc_cache;
private:
    virtual uint32_t  getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey)                            override;
    virtual void * getLocalData(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey)                            override;
    virtual void   flushElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) override;
    virtual void *invokeRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual)             override;
    virtual void *finishRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual, void *dst, const void *src) override;
    virtual void combineElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) override;
    hash::GHASH *hash_table = nullptr;

protected:
    auto lookup(uint32_t uHash, const void *pKey) {
    	return hash_table->lookup(uHash,pKey);
	}

public:
    auto fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual)
	{return arc_cache->fetch(uIndex,uId,bLock,bModify,bVirtual);}
    void *fetch(uint32_t uHash, void *pKey, int bLock,int bModify,bool bVirtual);
    void release(void *p) {return arc_cache->release(p);}
    void clear() {arc_cache->clear();}
protected:
    static void *getArrayElement(void *vData,int i,int iDataSize);
protected:
    class mdlClass * const mdl; // MDL is needed for cache operations
    mdlMessageCacheRequest CacheRequest;
public:
    void initialize(uint32_t cacheSize,
	void * (*getElt)(void *pData,int i,int iDataSize),
	void *pData,int iDataSize,int nData,
	std::shared_ptr<CACHEhelper> helper);

    void initialize_advanced(uint32_t cacheSize,hash::GHASH *hash,int iDataSize,std::shared_ptr<CACHEhelper> helper);

protected:
    void * (*getElt)(void *pData,int i,int iDataSize);
    void *pData;
    uint16_t iCID;
    // Fast lock to allow local updates while cache is active
    RWLock rwlock;
public:
    int iDataSize;
    int nData;
    uint32_t nLineBits;
    uint32_t getLineElementCount() const {return 1 << nLineBits;}
    uint32_t getLineMask()         const {return getLineElementCount()-1; }
    int iLineSize;
    std::vector<char> OneLine;
    /*
     ** Statistics stuff.
     */
    uint64_t nAccess;
    uint64_t nMiss;
public:
    explicit CACHE(mdlClass * mdl,uint16_t iCID);
    void close();
    bool isActive() {return cache_helper.get() != nullptr; }
    bool modify() {assert(cache_helper); return cache_helper->modify();}
    uint32_t key_size() {return arc_cache->key_size();}
    void *getElement(int i) {return (*getElt)(pData,i,iDataSize);}
    void *WriteLock(int iIndex) {rwlock.lock_write(); return getElement(iIndex);}
    void WriteUnlock(const void *p) {rwlock.unlock_write();}
    void *ReadLock(int iIndex) {rwlock.lock_read(); return getElement(iIndex);}
    void ReadUnlock(const void *p) {rwlock.unlock_read();}
    };

//*****************************************************************************

class mdlClass : public mdlBASE {
    int exit_code; // Only used by the worker thread
    friend class CACHE;
protected:
    friend class mdlMessageFlushToCore;
    void MessageFlushToCore(mdlMessageFlushToCore *message);
public:
    class mdlClass **pmdl;
    class mpiClass *mpi;
    void * (*fcnWorkerInit)(MDL mdl);
    void   (*fcnWorkerDone)(MDL mdl, void *ctx);
    int    (*fcnMaster)(    MDL mdl, void *ctx);

    void *worker_ctx;
    mdlMessageQueue threadBarrierQueue;

    std::vector<mdlMessageQueue> queueReceive; // Receive "Send/Ssend"
    void *pvMessageData; /* These two are for the collective malloc */
    size_t nMessageData;
    int iCoreMPI;             /* Core that handles MPI requests */
    int cacheSize;

    /* Work Queues */
    mdlMessageQueue wq;     /* Work for us to do */
    mdlMessageQueue wqDone; /* Completed work from other threads */
    mdlMessageQueue wqFree; /* Free work queue nodes */
    int wqMaxSize;
    uint16_t wqAccepting;
    uint16_t wqLastHelper;
    OPA_int_t wqCurSize;
    /*
     ** Services stuff!
     */
    int nMaxSrvBytes;
    std::vector<char> input_buffer;
    /*
     ** Swapping buffer.
     */
    std::vector<char> pszTrans;
    /*
     ** Caching stuff!
     */

    // Flush messages are added here. Eventually we send them to the MPI thread.
    mdlMessageQueue coreFlushBuffers; // Available buffers
    mdlMessageFlushFromCore *coreFlushBuffer; // Active buffer

    mdlMessageQueue queueCacheReply; // Replies to cache requests
    std::vector< std::unique_ptr<CACHE> >cache;
    mdlMessageQueue wqCacheFlush;

    int nFlushOutBytes;

#if defined(USE_CUDA) || defined(USE_CL)
    int inCudaBufSize, outCudaBufSize;
#endif
#ifdef USE_CL
    void *clCtx;
#endif
protected:
    void CommitServices();
    void Handler();
    int run_master();
    void mdl_MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag);
    int mdl_MPI_Sendrecv(
	void *sendbuf, int sendcount, MPI_Datatype sendtype,
	int dest, int sendtag, void *recvbuf, int recvcount,
	MPI_Datatype recvtype, int source, int recvtag, int *nReceived);
    void mdl_start_MPI_Ssend(mdlMessageSend &M, mdlMessageQueue &replyTo);
    int mdl_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, int *nBytes);
    virtual int checkMPI();
    void drainMPI();
    static void *mdlWorkerThread(void *vmdl); // Called by pthread_create with an mdlClass *
    void *WorkerThread();
    void combine_all_incoming();

    int DoSomeWork();
    void bookkeeping();
    void *finishCacheRequest(int cid,uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual, void *dst, const void *src);
  
protected:
    void flush_core_buffer();
    basicMessage & waitQueue(basicQueue &wait);
    void enqueue(mdlMessage &M);
    void enqueue(const mdlMessage &M, basicQueue &replyTo, bool bWait=false);
    void enqueueAndWait(const mdlMessage &M);
public:
#ifdef USE_CUDA
    int nCUDA;
    cudaMessageQueue cudaDone;
    int flushCompletedCUDA(); // Returns how many are still outstanding
    void enqueue(const cudaMessage &M);
    void enqueue(const cudaMessage &M, basicQueue &replyTo);
    void enqueueAndWait(const cudaMessage &M);
#endif
private:
    void init(bool bDiag = false);

public:
    explicit mdlClass(class mpiClass *mpi, int iMDL);
    explicit mdlClass(class mpiClass *mpi,
		int (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),
		int argc=0, char **argv=0);
    virtual ~mdlClass();

    CACHE *CacheInitialize(int cid,
	void * (*getElt)(void *pData,int i,int iDataSize),
	void *pData,int iDataSize,int nData,
	std::shared_ptr<CACHEhelper> helper);
    CACHE *AdvancedCacheInitialize(int cid,hash::GHASH *hash,int iDataSize,std::shared_ptr<CACHEhelper> helper);

    void *AcquireWrite(int cid,int iIndex) {
	return cache[cid]->WriteLock(iIndex);
	}

    void ReleaseWrite(int cid,void *p) {
	cache[cid]->WriteUnlock(p);
	}

    void FlushCache(int cid);
    void FinishCache(int cid);
    int ReqService(int id,int sid,void *vin=nullptr,int nInBytes=0);
    int GetReply(int rID,void *vout=nullptr);
    void Send(int id,mdlPack pack, void *ctx);
    void Recv(int id,mdlPack unpack, void *ctx);
    int Swap(int id,size_t nBufBytes,void *vBuf,size_t nOutBytes, size_t *pnSndBytes,size_t *pnRcvBytes);

    void CacheCheck();
    void CacheBarrier(int cid);
    void ThreadBarrier(bool bGlobal=false);
    void CompleteAllWork();
    bool isCudaActive();
    void Backtrace() {show_backtrace();}

    void *Access(int cid, uint32_t uIndex,  int uId, bool bLock,bool bModify,bool bVirtual);
    void *Access(int cid, uint32_t uHash,void *pKey, bool bLock,bool bModify,bool bVirtual);

    size_t FFTlocalCount(int n1,int n2,int n3,int *nz,int *sz,int *ny,int*sy);
    MDLFFT FFTNodeInitialize(int n1,int n2,int n3,int bMeasure,FFTW3(real) *data);
    void GridShare(MDLGRID grid);

    void FFT( MDLFFT fft, FFTW3(real) *data );
    void IFFT( MDLFFT fft, FFTW3(complex) *kdata );
    void Alltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps);

    };

class mpiClass : public mdlClass {
public:
    mdlMessageQueue queueWORK;
    mdlMessageQueue queueREGISTER;
protected:
    MPI_Comm commMDL;             /* Current active communicator */
    pthread_barrier_t barrier;
    std::vector<pthread_t> threadid;
#ifdef USE_CUDA
    CUDA cuda;
#endif
    mdlMessageQueue queueMPInew;     // Queue of work sent to MPI task
    std::vector<mdlMessageCacheRequest*> CacheRequestMessages;

    // Used to buffer incoming flush requests before sending them to each core
    mdlMessageQueue localFlushBuffers;
    std::vector<mdlMessageFlushToCore *> flushBuffersByCore;

    typedef std::list<mdlMessageFlushToRank *> FlushToRankList;
    FlushToRankList flushHeadFree; // Unused buffers
    FlushToRankList flushHeadBusy; // Buffers currently being filled (in flushBuffersByRank)
    std::vector<FlushToRankList::iterator> flushBuffersByRank;

    int flushBuffCount;
    std::vector<int> pRequestTargets;
    int iRequestTarget;
    int nActiveCores;
    int nOpenCaches;
    int iCacheBufSize;  /* Cache input buffer size */
    int iReplyBufSize;  /* Cache reply buffer size */
    std::vector<MPI_Request>    SendReceiveRequests;
    std::vector<MPI_Status>     SendReceiveStatuses;
    std::vector<int>            SendReceiveIndices;
    std::vector<mdlMessageMPI*> SendReceiveMessages;
#ifndef NDEBUG
    uint64_t nRequestsCreated, nRequestsReaped;
#endif
    MPI_Request *newRequest(mdlMessageMPI*message);

    mdlMessageCacheReceive *pReqRcv;
    std::list<mdlMessageCacheReply *> freeCacheReplies;

    // Cached FFTW plans
    struct fft_plan_information {
	ptrdiff_t nz, sz, ny, sy, nLocal;
	FFTW3(plan) fplan, iplan;
	};
    typedef std::tuple<ptrdiff_t,ptrdiff_t,ptrdiff_t> fft_plan_key; 
    std::map<fft_plan_key,fft_plan_information> fft_plans;

protected:
    // These are functions that are called as a result of a worker thread sending
    // us a message to perform a certain task.
    friend class mdlMessageSTOP;
    void MessageSTOP(mdlMessageSTOP *message);
    friend class mdlMessageBarrierMPI;
    void MessageBarrierMPI(mdlMessageBarrierMPI *message);
    friend class mdlMessageFlushFromCore;
    void MessageFlushFromCore(mdlMessageFlushFromCore *message);
    friend class mdlMessageFlushToRank;
    void MessageFlushToRank(mdlMessageFlushToRank *message);
    void FinishFlushToRank(mdlMessageFlushToRank *message);
    friend class mdlMessageCacheReply;
    void MessageCacheReply(mdlMessageCacheReply *message);
    void FinishCacheReply(mdlMessageCacheReply *message);
    friend class mdlMessageCacheReceive;
    void MessageCacheReceive(mdlMessageCacheReceive *message);
    void FinishCacheReceive(mdlMessageCacheReceive *message, const MPI_Status &status);
    void CacheReceiveRequest(int count, const CacheHeader *ph);
    void CacheReceiveReply(int count, const CacheHeader *ph);
    void CacheReceiveFlush(int count, CacheHeader *ph);
    friend class mdlMessageCacheOpen;
    void MessageCacheOpen(mdlMessageCacheOpen *message);
    friend class mdlMessageCacheClose;
    void MessageCacheClose(mdlMessageCacheClose *message);
    friend class mdlMessageCacheFlushOut;
    void MessageCacheFlushOut(mdlMessageCacheFlushOut *message);
    friend class mdlMessageCacheFlushLocal;
    void MessageCacheFlushLocal(mdlMessageCacheFlushLocal *message);
    friend class mdlMessageGridShare;
    void MessageGridShare(mdlMessageGridShare *message);
    friend class mdlMessageDFT_R2C;
    void MessageDFT_R2C(mdlMessageDFT_R2C *message);
    friend class mdlMessageDFT_C2R;
    void MessageDFT_C2R(mdlMessageDFT_C2R *message);
    friend class mdlMessageFFT_Sizes;
    void MessageFFT_Sizes(mdlMessageFFT_Sizes *message);
    friend class mdlMessageFFT_Plans;
    void MessageFFT_Plans(mdlMessageFFT_Plans *message);
    friend class mdlMessageAlltoallv;
    void MessageAlltoallv(mdlMessageAlltoallv *message);
    friend class mdlMessageSend;
    void MessageSend(mdlMessageSend *message);
    friend class mdlMessageReceive;
    void MessageReceive(mdlMessageReceive *message);
    friend class mdlMessageReceiveReply;
    void MessageReceiveReply(mdlMessageReceiveReply *message);
    void FinishReceiveReply(mdlMessageReceiveReply *message);
    friend class mdlMessageSendRequest;
    void MessageSendRequest(mdlMessageSendRequest *message);
    friend class mdlMessageSendReply;
    void MessageSendReply(mdlMessageSendReply *message);
    friend class mdlMessageCacheRequest;
    void MessageCacheRequest(mdlMessageCacheRequest *message);

protected:
    void flush_element(CacheHeader *pHdr,int iLineSize);
    void queue_local_flush(CacheHeader *ph);
    virtual int checkMPI();
    void processMessages();
    void finishRequests();
    int swapall(const char *buffer,int count,int datasize,/*const*/ int *counts);
public:
    explicit mpiClass(int (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),
    	     int argc=0, char **argv=0);
    virtual ~mpiClass();
    int Launch(int (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *));
    void KillAll(int signo);
#ifdef USE_CUDA
    void enqueue(const cudaMessage &M, basicQueue &replyTo);
    bool isCudaActive() {return cuda.isActive(); }
#else
    bool isCudaActive() {return false; }
#endif
    void enqueue(mdlMessage &M);
    void enqueue(const mdlMessage &M, basicQueue &replyTo, bool bWait=false);
    void pthreadBarrierWait();
    };
} // namespace mdl
#endif

#ifdef __cplusplus
extern "C" {
#endif

int mdlProcToThread(MDL mdl, int iProc);
int mdlThreadToProc(MDL mdl, int iThread);

/*
 ** General Functions
 */
int mdlLaunch(int,char **,int (*)(MDL,void *),void * (*)(MDL),void (*)(MDL,void *));
MDL   mdlMDL(void);
void *mdlWORKER(void);

void mdlAbort(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx);
void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx);
void mdlAddService(MDL,int,void *,fcnService_t *fcnService,int,int);
int  mdlReqService(MDL, int, int, void *, int);
void mdlGetReply(MDL,int,void *,int *);

static inline int mdlGridCoordCompare(const mdlGridCoord *a,const mdlGridCoord *b) {
    return a->x==b->x && a->y==b->y && a->z==b->z; 
    }

static inline mdlGridCoord *mdlGridCoordIncrement(mdlGridCoord *a) {
    ++a->i;
    ++a->II;
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

void mdlGridCoordFirstLast(MDL mdl,MDLGRID grid,mdlGridCoord *f,mdlGridCoord *l,int bCacheALign);

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

size_t mdlFFTlocalCount(MDL mdl,int n1,int n2,int n3,int *nz,int *sz,int *ny,int*sy);
MDLFFT mdlFFTNodeInitialize(MDL mdl,int nx,int ny,int nz,int bMeasure,FFTW3(real) *data);
MDLFFT mdlFFTInitialize(MDL mdl,int nx,int ny,int nz,int bMeasure,FFTW3(real) *data);
void mdlFFTNodeFinish( MDL mdl, MDLFFT fft );
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

void mdlAlltoallv(MDL mdl,int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps);

/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
void *mdlMallocArray(MDL mdl,size_t nmemb,size_t size,size_t minSize);
void *mdlSetArray(MDL mdl,size_t nmemb,size_t size,void *vdata);
void mdlFreeArray(MDL,void *);
void mdlSetCacheSize(MDL,int);
int mdlCacheStatus(MDL mdl,int cid); /* zero means not open */
void mdlROcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData);
void mdlCOcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData,
		void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,const void *));
void mdlAdvancedCacheRO(MDL mdl,int cid,void *pHash,int iDataSize);
void mdlAdvancedCacheCO(MDL mdl,int cid,void *pHash,int iDataSize,
	void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,const void *));
void mdlFinishCache(MDL,int);
void mdlCacheCheck(MDL);
void mdlCacheBarrier(MDL,int);
void mdlPrefetch(MDL mdl,int cid,int iIndex, int id);
void *mdlAcquire(MDL mdl,int cid,int iIndex,int id);
void *mdlFetch(MDL mdl,int cid,int iIndex,int id);
void *mdlVirtualFetch(MDL mdl,int cid,int iIndex,int id);
const void *mdlKeyFetch(MDL mdl,int cid,uint32_t uHash, void *pKey,int lock,int modify,int virt);
void *mdlKeyAcquire(MDL mdl,int cid,uint32_t uHash, void *pKey);
void mdlRelease(MDL,int,void *);
void mdlFlushCache(MDL,int);
void mdlThreadBarrier(MDL);
void mdlCompleteAllWork(MDL);
void *mdlAcquireWrite(MDL mdl, int cid, int iIndex);
void mdlReleaseWrite(MDL mdl,int cid,void *p);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);

void mdlSetWorkQueueSize(MDL,int,int);
void mdlSetCudaBufferSize(MDL,int,int);
int mdlCudaActive(MDL mdl);
void mdlAddWork(MDL mdl, void *ctx,
    int (*initWork)(void *ctx,void *vwork),
    int (*checkWork)(void *ctx,void *vwork),
    mdlWorkFunction doWork,
    mdlWorkFunction doneWork);

void mdlTimeReset(MDL mdl);
double mdlTimeComputing(MDL mdl);
double mdlTimeSynchronizing(MDL mdl);
double mdlTimeWaiting(MDL mdl);
double mdlWallTime(MDL mdl);
void mdlprintf(MDL mdl, const char *format, ...);

#ifdef __cplusplus
}
#endif
#endif
