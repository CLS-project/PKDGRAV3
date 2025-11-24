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
#include <boost/circular_buffer.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/qvm/vec.hpp>
#include "mdlstubs.h"
#include "mdlbase.h"
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
#endif
#include <stdint.h>
#include "mpi.h"
#include "mdlfft.h"



#ifdef INSTRUMENT_MPI
#include "WelfordStats.h"
using StatsAccum = WelfordStats<uint64_t>;
#else
struct StatsAccum {
    template <typename T>
    StatsAccum &operator+=(T) { return *this; }
};
#endif
class Instrument {
public:
    Instrument(StatsAccum &accumulator)
        : acc_(accumulator), start_(getticks()) {}

    ~Instrument() {
        auto end = getticks();
        acc_ += elapsed(end, start_);
    }
    void reset() {}

private:
    StatsAccum &acc_;
    ticks start_;
};
#include "arc.h"
#include "mpimessages.h"
#include "gpu/mdlgpu.h"
#ifdef USE_CUDA
    #include "mdlcuda.h"
#endif
#ifdef USE_METAL
    #include "metal/mdlmetal.h"
#endif
#include "rwlock.h"
#include <chrono>
#include <memory>
#include <tuple>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <optional>

#ifndef MPI_VERSION
    #define MPI_VERSION 1
#endif

#define SRV_STOP        0

/*                                   MQ   MPI */
#define MDL_TAG_BARRIER         1 /* Yes  Yes */
#define MDL_TAG_SWAPINIT    2 /* NO   Yes */
#define MDL_TAG_CACHE_FLUSH     2 /* Yes  NO  */
#define MDL_TAG_SWAP        3 /* NO   Yes */
#define MDL_TAG_REQ     4 /* NO   Yes */
#define MDL_TAG_RPL     5 /* This is treated specially */
#define MDL_TAG_SEND            6 /* NO   Yes */
#define MDL_TAG_CACHECOM    7 /* Yes  Yes */

#define MDL_TAG_MAX             8

typedef int (*mdlWorkFunction)(void *ctx);
typedef int (*mdlPack)(void *,int *,size_t,void *);

#ifdef __cplusplus
namespace mdl {
static const auto mdl_cache_size = 15'000'000;
static const auto mdl_check_mask = 0x7f;

// typedef struct cacheHeaderNew {
//     uint32_t mid     :  2; /*  2: Message ID: request, flush, reply */
//     uint32_t cid     :  6; /*  6: Cache for this entry */
//     uint32_t tidFrom : 12; /* 12: thread id that made the request */
//     uint32_t tidTo   : 12; /* 12: thread id that gets the reply */
//     uint32_t iLine;        /* Index of the cache line */
//     } CacheHeadernew;

#define MDL_CACHE_DATA_SIZE (4096)

enum class MPI_RECV_MODE {
    TOTAL = 0,
    PER_RANK = 1
};


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
    uint32_t nElements;  // Total number of elements in the cache
    uint32_t nDataSize;  // Size (bytes) of an element in the source container
    uint32_t nPackSize;  // Size (bytes) of a packed element
    uint32_t nFlushSize; // Size (bytes) of a flushed element
    uint32_t nCacheSize; // Size (bytes) of an element in the cache
    // NOTE: nCacheSize >= nDataSize. If larger this means there are decorations after the data.
    bool     bModify;    // If this cache will be modified (and thus flushed)
    // These overrides only work if nDataSize == nPackSize == nFlushSize == nCacheSize
    // This is the most common case and the "classic" cache
    virtual void    pack(void *dst, const void *src)                  { memcpy(dst,src,nDataSize); }
    virtual void  unpack(void *dst, const void *src, const void *key) { memcpy(dst,src,nDataSize); }
    virtual void    init(void *dst)                                   {}
    virtual void   flush(void *dst, const void *src)                  { memcpy(dst,src,nDataSize); }
    virtual void combine(void *dst, const void *src, const void *key) {}
    virtual void *create(uint32_t size, const void *pKey) { return nullptr; }
    virtual uint32_t getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) {return uId;}
    virtual bool  modify() {return bModify; }
    virtual void *get_element(uint32_t i) {return nullptr;}
    virtual void *get_element(uint32_t uHash, void *pKey) {return nullptr;}
public:
    explicit CACHEhelper(uint32_t nDataSize, uint32_t nElements, bool bModify=false)
        : nElements(nElements), nDataSize(nDataSize), nPackSize(nDataSize), nFlushSize(nDataSize), nCacheSize(nDataSize), bModify(bModify) {}
    virtual ~CACHEhelper() = default;
    auto data_size()  const {return nDataSize;}
    auto pack_size()  const {return nPackSize;}
    auto cache_size() const {return nCacheSize;}
    auto flush_size() const {return nFlushSize;}
    auto elements()   const {return nElements;}
};

class arrayCACHEhelper : public CACHEhelper {
protected:
    char *pArray;
public:
    explicit arrayCACHEhelper(uint32_t nDataSize, void *pArray, uint32_t nArray, bool bModify=false)
        : CACHEhelper(nDataSize,nArray,bModify), pArray(static_cast<char *>(pArray)) {}
    template<typename T>
    explicit arrayCACHEhelper(T *pArray, uint32_t nArray, bool bModify=false)
        : CACHEhelper(sizeof(T),nArray,bModify), pArray(reinterpret_cast<char *>(pArray)) {}

    virtual void *get_element(uint32_t i) override {
        return pArray + static_cast<size_t>(i) * data_size();
    }
};

class hashCACHEhelper : public CACHEhelper {
protected:
    hash::GHASH *hash_table;
public:
    explicit hashCACHEhelper(hash::GHASH *hash_table, uint32_t nDataSize, bool bModify=false)
        : CACHEhelper(nDataSize,0,bModify), hash_table(hash_table) {}
    virtual void *get_element(uint32_t uHash, void *pKey) override final {
        return hash_table->lookup(uHash,pKey);
    }
};

class CACHE : private ARChelper {
    friend class mdlClass;
    friend class mpiClass;
    std::shared_ptr<CACHEhelper> cache_helper;
    std::unique_ptr<GARC> arc_cache;

    virtual uint32_t  getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey)                            override;
    virtual void *getLocalData(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey)                            override;
    virtual void   flushElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) override;
    virtual void *invokeRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual)             override;
    virtual void *finishRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual, void *dst, const void *src) override;
    virtual void combineElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) override;
    hash::GHASH *hash_table = nullptr;
    class mdlClass *const mdl;  // MDL is needed for cache operations
    mdlMessageCacheRequest CacheRequest;

    auto lookup(uint32_t uHash, const void *pKey) {
        return hash_table->lookup(uHash,pKey);
    }

    auto fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual)
    {return arc_cache->fetch(uIndex,uId,bLock,bModify,bVirtual);}
    void *fetch(uint32_t uHash, void *pKey, int bLock,int bModify,bool bVirtual);
    void release(void *p) {return arc_cache->release(p);}
    void clear() {arc_cache->clear();}

    void initialize(uint32_t cacheSize,std::shared_ptr<CACHEhelper> helper);

    void initialize_advanced(uint32_t cacheSize,hash::GHASH *hash,std::shared_ptr<CACHEhelper> helper);

    uint16_t iCID;
    // Fast lock to allow local updates while cache is active
    RWLock rwlock;
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

    void close();
    bool isActive() {return cache_helper.get() != nullptr; }
    bool modify() {assert(cache_helper); return cache_helper->modify();}
    uint32_t key_size() {return arc_cache->key_size();}
    void *getElement(int i) {return cache_helper->get_element(i);}
    void *WriteLock(int iIndex) {rwlock.lock_write(); return getElement(iIndex);}
    void WriteUnlock(const void *p) {rwlock.unlock_write();}
    void *ReadLock(int iIndex) {rwlock.lock_read(); return getElement(iIndex);}
    void ReadUnlock(const void *p) {rwlock.unlock_read();}
public:
    explicit CACHE(mdlClass *mdl,uint16_t iCID);
    virtual ~CACHE() = default;
};

//*****************************************************************************

class mdlClass : public mdlBASE {
    int exit_code; // Only used by the worker thread
    friend class CACHE;
    friend class mpiClass;
protected:
    friend class mdlMessageFlushToCore;
    void MessageFlushToCore(mdlMessageFlushToCore *message);
    // Domain decomposition
public:
    using dd_offset_type = uint64_t;
protected:
    char *ddBuffer;
    std::vector<dd_offset_type> ddOffset;
    std::vector<dd_offset_type> ddCounts;
    dd_offset_type ddFreeOffset[2], ddFreeCounts[2]; // For global Swap
public:
    /*
     ** Special MDL memory allocation functions for allocating memory
     ** which must be visible to other processors thru the MDL cache
     ** functions.
     ** Malloc() is defined to return a pointer to AT LEAST iSize bytes
     ** of memory. This pointer will be passed to cache functions.
     ** For PVM and most machines these functions are trivial, but on the
     ** T3D and perhaps some future machines these functions are required.
     */
    template<typename T=void>
    auto Malloc(size_t iSize) const {
        return static_cast<T *>(malloc(iSize));
    }
    template<typename T>
    void Free(T *p) const {
        free(p);
    }

protected:
    void *exchangeArray(void *vdata);
public:
    template<typename T>
    auto SetArray(void *vdata,size_t nmemb=0) {
        if constexpr (std::is_void_v<T>) nMessageData = 0;
        else nMessageData = nmemb * sizeof(T);
        return static_cast<T *>(exchangeArray(vdata));
    }
    template<typename T>
    auto SetArray(T *vdata,size_t nmemb=0) {
        if constexpr (std::is_void_v<T>) nMessageData = 0;
        else nMessageData = nmemb * sizeof(T);
        return static_cast<T *>(exchangeArray(vdata));
    }

    class mpiClass *mpi;
protected:
    class mdlClass **pmdl;
    void *(*fcnWorkerInit)(MDL mdl);
    void   (*fcnWorkerDone)(MDL mdl, void *ctx);
    int    (*fcnMaster)(    MDL mdl, void *ctx);

    void *worker_ctx;
    mdlMessageQueue threadBarrierQueue;

    std::unique_ptr<mdlMessageQueue[]> queueReceive; // Receive "Send/Ssend"
    void *pvMessageData; /* These two are for the collective malloc */
    size_t nMessageData;
    void **pSegments;
    boost::container::small_vector<uint64_t, 8> nBytesSegment;
    int iCoreMPI;             /* Core that handles MPI requests */
    int cacheSize;

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
    std::vector< std::unique_ptr<CACHE>>cache;
    mdlMessageQueue wqCacheFlush;

    int nFlushOutBytes;

protected:
    void CommitServices();
    void Handler();
    int run_master();
    void mdl_MPI_Ssend(void *buf, int count, int dest, int tag);
    int mdl_MPI_Sendrecv(
        void *sendbuf, int sendcount,
        int dest, int sendtag, void *recvbuf, int recvcount,
        int source, int recvtag, int *nReceived);
    void mdl_start_MPI_Ssend(mdlMessageSend &M, mdlMessageQueue &replyTo);
    int mdl_MPI_Recv(void *buf, int count, int source, int tag, int *nBytes);
    virtual int checkMPI();
    void drainMPI();
    void *WorkerThread();
    void combine_all_incoming();

    void bookkeeping();
    void *finishCacheRequest(int cid,uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual, void *dst, const void *src);

protected:
    void flush_core_buffer();
    basicMessage &waitQueue(basicQueue &wait);
    void enqueue(mdlMessage &M);
    void enqueue(const mdlMessage &M, basicQueue &replyTo, bool bWait=false);
    void enqueueAndWait(const mdlMessage &M);
public:
    gpu::Client gpu;
#ifdef USE_METAL
    void enqueue(const metal::metalMessage &M);
    void enqueue(const metal::metalMessage &M, basicQueue &replyTo);
    void enqueueAndWait(const metal::metalMessage &M);
#endif
private:
    void init(bool bDiag = false);

public:
    explicit mdlClass(class mpiClass *mpi, int iMDL);
    explicit mdlClass(class mpiClass *mpi,
                      int (*fcnMaster)(MDL,void *),void *(*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),
                      int argc=0, char **argv=0);
    virtual ~mdlClass();
    void SetCacheMaxInflight(int iMax);
    void SetMpiRecvCount(int count, MPI_RECV_MODE mode=MPI_RECV_MODE::TOTAL);

    void SetCacheSize(int cacheSize) {
        this->cacheSize = cacheSize;
    }
    CACHE *CacheInitialize(int cid,std::shared_ptr<CACHEhelper> helper);
    CACHE *AdvancedCacheInitialize(int cid,hash::GHASH *hash,std::shared_ptr<CACHEhelper> helper);

    // Cache access
public:
    template<typename T=void>
    auto Fetch(int cid,int iIndex,int id) {
        const bool lock   = false; // we never lock in fetch
        const bool modify = false; // fetch can never modify
        const bool virt   = false; // really fetch the element
        return static_cast<T *>(Access(cid, iIndex, id, lock, modify, virt));
    }
    template<typename T=void>
    auto KeyFetch(int cid,uint32_t uHash, void *pKey,int lock,int modify,int virt) {
        return static_cast<T *>(Access(cid, uHash, pKey, lock, modify, virt));
    }
    template<typename T=void>
    auto Acquire(int cid,int iIndex,int id) {
        const bool lock   = true;  // we always lock in acquire
        const bool modify = cache[cid]->modify();
        const bool virt   = false; // really fetch the element
        return static_cast<T *>(Access(cid, iIndex, id, lock, modify, virt));
    }
    template<typename T=void>
    auto KeyAcquire(int cid,uint32_t uHash, void *pKey) {
        const bool lock   = true;  // we always lock in acquire
        const bool modify = cache[cid]->modify();
        const bool virt   = false; // really fetch the element
        return static_cast<T *>(Access(cid, uHash, pKey, lock, modify, virt));
    }
    /* Locks the element, but does not fetch or initialize */
    template<typename T=void>
    auto VirtualFetch(int cid,int iIndex,int id) {
        const int lock   = false; // fetch never locks
        const int modify = true;  // virtual always modifies
        const bool virt  = true;  // do not fetch the element
        return static_cast<T *>(Access(cid, iIndex, id, lock, modify, virt));
    }
    void Release(int cid,void *p) {
        cache[cid]->release(p);
    }
    template<typename T=void>
    auto AcquireWrite(int cid,int iIndex) {
        return static_cast<T *>(cache[cid]->WriteLock(iIndex));
    }
    void ReleaseWrite(int cid,void *p) {
        cache[cid]->WriteUnlock(p);
    }

    void FlushCache(int cid);
    void FinishCache(int cid);
    int ReqService(int id,int sid,void *vin,int nInBytes);
    int ReqService(int id,int sid) {
        return ReqService(id,sid,nullptr,0);
    }
    template<typename T>
    int ReqService(int id,int sid,T &in) {
        return ReqService(id,sid,&in,sizeof(in));
    }
    int GetReply(int rID,int nMaxBytes=0,void *vout=nullptr);
    template<typename T, typename std::enable_if<!std::is_pointer<T>::value, int>::type = 0>
    int GetReply(int rID,T &out) {
        return GetReply(rID,sizeof(out),&out);
    }
    void Send(int id,mdlPack pack, void *ctx);
    void Recv(int id,mdlPack unpack, void *ctx);
    int Swap(int id,size_t nBufBytes,void *vBuf,size_t nOutBytes, size_t *pnSndBytes,size_t *pnRcvBytes);

    void CacheCheck();
    int  CacheStatus(int cid);
    void CacheBarrier(int cid);
    int  ThreadBarrier(bool bGlobal=false,int iVote=0);
    void CompleteAllWork();
    bool CudaActive();
    bool isMetalActive();
    int numGPUs();
    void Backtrace() {show_backtrace();}

    void *Access(int cid, uint32_t uIndex,  int uId, bool bLock,bool bModify,bool bVirtual);
    void *Access(int cid, uint32_t uHash,void *pKey, bool bLock,bool bModify,bool bVirtual);

    void GridShare(MDLGRID grid);
#ifdef MDL_FFTW
    size_t FFTlocalCount(int n1,int n2,int n3,int *nz,int *sz,int *ny,int *sy);
    MDLFFT FFTNodeInitialize(int n1,int n2,int n3,int bMeasure,FFTW3(real) *data);
    MDLFFT FFTInitialize(int nx,int ny,int nz,int bMeasure,FFTW3(real) *data);
    void FFTNodeFinish( MDLFFT fft );
    void FFTFinish( MDLFFT fft );

    void FFT( MDLFFT fft, FFTW3(real) *data );
    void IFFT( MDLFFT fft, FFTW3(complex) *kdata );
#endif
    void Alltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps);

    int swaplocal( void *buffer,dd_offset_type count,dd_offset_type datasize,/*const*/ dd_offset_type *counts);
    int swapglobal(void *buffer,dd_offset_type count,dd_offset_type datasize,/*const*/ dd_offset_type *counts);
    size_t new_shared_array(void **p, int nSegments,  size_t *nElements,const size_t *nBytesPerElement,size_t nMinTotalStore=0);
    void delete_shared_array(void *p,size_t nBytes);

    double NumAccess(int cid) {
        return cache[cid]->nAccess;
    }
    double MissRatio(int cid) {
        auto &c = cache[cid];
        double dAccess = c->nAccess;
        return dAccess > 0.0 ? (c->nMiss/dAccess) : 0.0;
    }
};

class mpiClass : public mdlClass {
public:
    mdlMessageQueue queueWORK;
    mdlMessageQueue queueREGISTER;
#ifdef USE_CUDA
    CUDA cuda;
#endif
#ifdef USE_METAL
    metal::METAL metal;
#endif
protected:
    MPI_Comm commMDL;             /* Current active communicator */
    std::optional<boost::barrier> barrier;
    boost::thread_group worker_threads; // Worker threads
    mdlMessageQueue queueMPInew;     // Queue of work sent to MPI task
    std::vector<mdlMessageCacheRequest *> CacheRequestMessages;
    std::vector<bool> CacheRequestRendezvous;

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
    // By setting this to a positive value, it limits the number of cache messages inflight
    // to each rank to the specified number. When enabled (not zero), all REQUEST and RESPONSE
    // messages are added to the flush buffer. The flush buffer is sent when below this limit.
    // To avoid deadlock, the flush buffer is sent if there is at least one REPLY message after
    // processing a batch of incoming REQUEST/REPLY/FLUSH messages.
    int iCacheMaxInflight = 0;
    int iCacheBufSize;  /* Cache input buffer size */
    int iReplyBufSize;  /* Cache reply buffer size */
    MPI_RECV_MODE mpi_recv_mode;    // TOTAL or PER_RANK
    int mpi_recv_count;             // Number of outstanding receives to post

    std::vector<MPI_Request>    SendReceiveRequests;
    std::vector<MPI_Status>     SendReceiveStatuses;
    std::vector<int>            SendReceiveAvailable;
    std::vector<int>            SendReceiveIndices;
    std::vector<mdlMessageMPI *> SendReceiveMessages;
    struct BufferedCacheRequest {
        CacheHeader header;
        const void *data;
        BufferedCacheRequest() = default;
        BufferedCacheRequest(const CacheHeader *header,const void *data=nullptr) : header(*header), data(data) {}
    };
    boost::circular_buffer<BufferedCacheRequest> PendingRequests;
    void BufferCacheResponse(FlushBuffer *flush,BufferedCacheRequest &request);

#ifndef NDEBUG
    uint64_t nRequestsCreated, nRequestsReaped;
#endif
    MPI_Request *newRequest(mdlMessageMPI *message,MPI_Request request=MPI_REQUEST_NULL);

    std::vector<std::unique_ptr<mdlMessageCacheReceive>> msgCacheReceive;
    std::vector<uint32_t> countCacheInflight;
#ifdef DEBUG_COUNT_CACHE
    std::vector<uint64_t> countCacheSend, countCacheRecv;
#endif
    std::list<mdlMessageCacheReply *> freeCacheReplies;
#ifdef MDL_FFTW
    // Cached FFTW plans
    struct fft_plan_information {
        ptrdiff_t nz, sz, ny, sy, nLocal;
        FFTW3(plan) fplan, iplan;
    };
    typedef std::tuple<ptrdiff_t,ptrdiff_t,ptrdiff_t> fft_plan_key;
    std::map<fft_plan_key,fft_plan_information> fft_plans;
#endif

protected:
    // These are functions that are called as a result of a worker thread sending
    // us a message to perform a certain task.
    friend class mdlMessageSTOP;
    void MessageSTOP(mdlMessageSTOP *message);
    friend class mdlMessageBarrierMPI;
    void MessageBarrierMPI(mdlMessageBarrierMPI *message);
    friend class mdlMessageSwapGlobal;
    void MessageSwapGlobal(mdlMessageSwapGlobal *message);
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
    void FinishCacheReceive(mdlMessageCacheReceive *message, MPI_Request request, MPI_Status status);
    void CacheReceive(int bytes, CacheHeader *ph, int iProcFrom);
    int CacheReceiveRequest(int count, CacheHeader *ph, int iProcFrom);
    int CacheReceiveReply(int count, CacheHeader *ph, int iProcFrom);
    int CacheReceiveFlush(int count, CacheHeader *ph, int iProcFrom);
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
#ifdef MDL_FFTW
    friend class mdlMessageDFT_R2C;
    void MessageDFT_R2C(mdlMessageDFT_R2C *message);
    friend class mdlMessageDFT_C2R;
    void MessageDFT_C2R(mdlMessageDFT_C2R *message);
    friend class mdlMessageFFT_Sizes;
    void MessageFFT_Sizes(mdlMessageFFT_Sizes *message);
    friend class mdlMessageFFT_Plans;
    void MessageFFT_Plans(mdlMessageFFT_Plans *message);
#endif
    friend class mdlMessageAlltoallv;
    void MessageAlltoallv(mdlMessageAlltoallv *message);
    friend class mdlMessageSend;
    void MessageSend(mdlMessageSend *message);
    friend class mdlMessageReceive;
    void MessageReceive(mdlMessageReceive *message);
    friend class mdlMessageReceiveReply;
    void MessageReceiveReply(mdlMessageReceiveReply *message);
    void FinishReceiveReply(mdlMessageReceiveReply *message, MPI_Request request, MPI_Status status);
    friend class mdlMessageSendRequest;
    void MessageSendRequest(mdlMessageSendRequest *message);
    friend class mdlMessageSendReply;
    void MessageSendReply(mdlMessageSendReply *message);
    friend class mdlMessageCacheRequest;
    void MessageCacheRequest(mdlMessageCacheRequest *message);
    void FinishCacheRequest(mdlMessageCacheRequest *message, MPI_Request request, MPI_Status status);

protected:
    void expedite_flush(int iProc);
    mdlMessageFlushToRank *get_flush_buffer(int iProc,int iSize,bool bWait=true);
    void flush_element(CacheHeader *pHdr,int iLineSize);
    void queue_local_flush(CacheHeader *ph);
    virtual int checkMPI();
    int processMessages();
    void finishRequests();
public:
    explicit mpiClass(int (*fcnMaster)(MDL,void *),void *(*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),
                      int argc=0, char **argv=0);
    virtual ~mpiClass();
    void SetCacheMaxInflight(int iMax) {iCacheMaxInflight = iMax;}
    void SetMpiRecvCount(int count, MPI_RECV_MODE mode=MPI_RECV_MODE::TOTAL) {
        assert(count>=1 && (mode==MPI_RECV_MODE::TOTAL || mode==MPI_RECV_MODE::PER_RANK));
        this->mpi_recv_count = count;
        this->mpi_recv_mode  = mode;
    }
    int Launch();
#ifdef USE_CUDA
    bool CudaActive() {return cuda.isActive(); }
    int numGPUs() {return cuda.numDevices();}
#else
    bool CudaActive() {return false; }
    int numGPUs() {return 0;}
#endif
#ifdef USE_METAL
    void enqueue(const metal::metalMessage &M, basicQueue &replyTo);
    bool isMetalActive() {return metal.isActive(); }
#else
    bool isMetalActive() {return false; }
#endif
    void enqueue(mdlMessage &M);
    void enqueue(const mdlMessage &M, basicQueue &replyTo, bool bWait=false);
    void ThreadBarrierWait();
public:
    // The following are only used in MPI instrumentation is enabled (INSTRUMENT_MPI)
    struct {
        ticks ticks_cache_open;         // Tick counter when the first cache was opened
        ticks ticks_cache_close;        // Tick counter when the last cache was closed
        std::chrono::time_point<std::chrono::high_resolution_clock> time_cache_open;  // Time point when the first cache was opened
        std::chrono::time_point<std::chrono::high_resolution_clock> time_cache_close; // Time point when the last cache was closed
        StatsAccum mpi_loop;  // Statistics on the time spent in an MPI loop iteration
        StatsAccum mpi_test;  // Statistics on the time spent in an MPI test iteration
        StatsAccum gpu_test;  // Statistics on the time spent in a GPU test iteration
        StatsAccum msg_test;  // Statistics on the time spent in a message test iteration
        StatsAccum cache_request; // Statistics on the time spent processing a cache request
        StatsAccum cache_reply;   // Statistics on the time spent between the send and recv
    } cache_stats;
#ifdef INSTRUMENT_MPI
public:
    std::tuple<WelfordStats<double>, WelfordStats<double>, WelfordStats<double>, WelfordStats<double>>
    get_cache_stats() {
        auto total_ticks = elapsed(getticks(), cache_stats.ticks_cache_open);
        auto total_seconds = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - cache_stats.time_cache_open).count();
        auto seconds_per_tick = total_seconds / total_ticks;
        return {
            cache_stats.mpi_loop * seconds_per_tick,
            cache_stats.mpi_test * seconds_per_tick,
            cache_stats.gpu_test * seconds_per_tick,
            cache_stats.msg_test *seconds_per_tick
        };
    }
#endif
};
int Launch(int (*fcnMaster)(MDL,void *),void *(*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),int argc=0,char **argv=nullptr);
} // namespace mdl
#endif // __cplusplus

#ifdef __cplusplus
extern "C" {
#endif

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
** This gives the processor on which the given slab can be found.
*/
static inline int mdlGridId(MDL vmdl,MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    auto mdl = static_cast<mdl::mdlClass *>(vmdl);
    assert(z<grid->n3);
    return mdl->ProcToThread(grid->id[z]);
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

/* Grid accessors: r-space */
#define mdlFFTrId(mdl,fft,x,y,z) mdlGridId(mdl,(fft)->rgrid,x,y,z)
#define mdlFFTrIdx(mdl,fft,x,y,z) mdlGridIdx(mdl,(fft)->rgrid,x,y,z)

/* Grid accessors: k-space (note permuted indices) */
#define mdlFFTkId(mdl,fft,x,y,z) mdlGridId(mdl,(fft)->kgrid,x,z,y)
#define mdlFFTkIdx(mdl,fft,x,y,z) mdlGridIdx(mdl,(fft)->kgrid,x,z,y)

#endif

void mdlprintf(MDL mdl, const char *format, ...);

#ifdef __cplusplus
}
#endif
#endif
