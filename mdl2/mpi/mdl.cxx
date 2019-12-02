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

#include "mdl.h"

static inline int size_t_to_int(size_t v) {
    return (int)v;
    }


#ifdef USE_HWLOC
#include "hwloc.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>


#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef __linux__
#include <sys/resource.h>
#endif
#include "mpi.h"
#ifdef USE_CUDA
#include "cudautil.h"
#endif
#ifdef USE_CL
#include "clutil.h"
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif

#ifdef __cplusplus
#define CAST(T,V) reinterpret_cast<T>(V)
#else
#define CAST(T,V) ((T)(V))
#endif

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_DEFAULT_CACHEIDS	12

#define MDL_TRANS_SIZE		5000000
#define MDL_FLUSH_DATA_SIZE	32000
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

/*
** The MPI specification allows tags to at least 32767 (but most implementations go higher).
** We multiply this offset by the destination core, so 32768/8 = 4096 cores max per process,
** or more if the MPI implementation allows it.
*/
#define MDL_TAG_THREAD_OFFSET   MDL_TAG_MAX


/* Core to core queueing. */
#define MDL_SE_STOP             0
#define MDL_SE_SEND_REQUEST     1
#define MDL_SE_SEND_REPLY       2
#define MDL_SE_RECV_REQUEST     3
#define MDL_SE_RECV_REPLY       4
#define MDL_SE_BARRIER_REQUEST  5
#define MDL_SE_BARRIER_REPLY    6
#define MDL_SE_CACHE_OPEN       7
#define MDL_SE_CACHE_CLOSE      8
#define MDL_SE_CACHE_REQUEST    9
#define MDL_SE_CACHE_FLUSH     10
#define MDL_SE_CACHE_FLUSH_OUT 11
#define MDL_SE_CACHE_FLUSH_LCL 12
#define MDL_SE_MPI_SEND        13
#define MDL_SE_MPI_SSEND       14
#define MDL_SE_MPI_RECV        15

#define MDL_SE_FFT_SIZES       20
#define MDL_SE_FFT_PLANS       21
#define MDL_SE_FFT_DFT_R2C     22
#define MDL_SE_FFT_DFT_C2R     23

#define MDL_SE_GRID_SHARE      30

#define MDL_SE_ALLTOALLV       40

/*
** The bit values of these flags are NOT arbitrary!
** The implementation depends on their values being set this way.
*/
#define _P1_     0x10000000
#define _T1_     0x00000000
#define _B1_     0x20000000
#define _T2_     0x40000000
#define _B2_     0x60000000
#define _WHERE_  0x70000000
#define _IDMASK_ 0x0fffffff
#define _DIRTY_  0x80000000
#define _ARC_MAGIC_ 0xa6c3a91c00000000
#define _ARC_MASK_  0xa6c3a91cffffffff

static uint32_t swar32(register uint32_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return(x);
}

static inline uint32_t murmur2(const uint32_t *key, int len) {
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    uint32_t h = 0xdeadbeef /*^ len : len will be the same */;
    while(len--) {
	uint32_t k = *key++;
	k *= m;
	k ^= k >> r;
	k *= m;
	h *= m;
	h ^= k;
	}
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
    } 

/*
** This makes the following assumptions (changed from regular hash)
**   1. keys lengths are a multiple of four bytes (and at least four bytes)
**   2. length is always identical (so don't need to mix in the length)
**   3. We will always use the same seed
*/
static inline uint32_t murmur3(const uint32_t* key, size_t len) {
    uint32_t h = 0xdeadbeef;
    do {
	uint32_t k = *key++;
	k *= 0xcc9e2d51;
	k = (k << 15) | (k >> 17);
	k *= 0x1b873593;
	h ^= k;
	h = (h << 13) | (h >> 19);
	h = h * 5 + 0xe6546b64;
	} while (--len);
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
    }

/*
** MurmurHash2, by Austin Appleby
** adapted for hashing 2 uint32_t variables for mdl2
*/
static inline uint32_t MurmurHash2(uint32_t a,uint32_t b) {
    /* 
    ** 'm' and 'r' are mixing constants generated offline.
    ** They're not really 'magic', they just happen to work well.
    */
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    uint32_t h = 0xdeadbeef;

    /* Mix the 2 32-bit words into the hash */
    a *= m; 
    b *= m; 
    a ^= a >> r; 
    b ^= b >> r; 
    a *= m; 	
    b *= m; 	
    /* now work on the hash */
    h ^= a;
    h *= m; 
    h ^= b;	
    /* Do a few final mixes of the hash to ensure the last few
    ** bytes are well-incorporated. */
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
    }

CDB::CDB() {
    next = prev = this;
    uId = 0xdeadbeef;
    }

CDB *CDB::remove_from_list() {
    prev->next = next;
    next->prev = prev;
    return this;
    }

CDB *CDB::lru_remove() {
    return prev->remove_from_list();
    }

void CDB::lru_insert(CDB *p) {
    p->prev = prev;
    p->next = this;
    prev->next = p;
    prev = p;
    }

void CDB::mru_insert(CDB *p) {
    next->lru_insert(p);
    }

CDB * ARC::remove_from_hash(CDB *p) {
    ARC *arc = this;
    uint32_t uPage = p->uPage;
    uint32_t uId = p->uId&_IDMASK_;
    uint32_t uHash = (MurmurHash2(uPage,uId)&arc->uHashMask);
    CDB **pt;

    for(pt = &arc->Hash[uHash]; *pt != NULL; pt = &((*pt)->coll)) {
	if ( *pt == p) {
	    *pt = (*pt)->coll;
	    return p;
	    }
	}
    printf("Tried to remove uPage %d uId %x\n", uPage, uId);
    abort();  /* should never get here, the element should always be found in the hash */
    }

#ifdef HAVE_MEMKIND
#include <hbwmalloc.h>
#define ARC_malloc hbw_malloc
#define ARC_free hbw_free
#else
#define ARC_malloc malloc
#define ARC_free free
#endif

ARC::ARC(mdlClass * mdlIn,uint32_t nCacheIn,uint32_t uDataSizeIn,CACHE *c)
    : mdl(mdlIn), nCache(nCacheIn), uDataSize((uDataSizeIn+7)>>3), cache(c) {
    uint32_t i;

    /*
    ** Allocate stuff.
    */
    cdbBase = new CDB[2*nCache];
    /*
    ** Make sure we have sufficient alignment of data.
    ** In this case to nearest long word (8 bytes).
    ** We add one long word at the start to store the 
    ** magic number and lock count.
    */
    dataBase = new uint64_t[nCache*(uDataSize+1)];
    dataLast = dataBase + nCache*(uDataSize+1);
    /*
    ** Determine nHash.
    */
    uHashMask = swar32(3*nCache-1);
    nHash = uHashMask+1; 
    Hash = new CDB*[nHash];
    for (i=0;i<nHash;++i) Hash[i] = NULL;
    /*
    ** Create all lists as circular lists, with a sentinel
    ** CDB at the head/tail of the list.
    ** Initialize the lengths of the various lists.
    */
    T1 = new CDB;    T1Length = 0;
    B1 = new CDB;    B1Length = 0;
    T2 = new CDB;    T2Length = 0;
    B2 = new CDB;    B2Length = 0;
    Free = new CDB;
    /*
    ** Initialize target T1 length.
    */
    target_T1 = nCache/2;   /* is this ok? */
    /*
    ** Insert CDBs with data into the Free list first.
    */
    for (i=0;i<nCache;++i) {
	/*dataBase[i*(uDataSize+1)] = _ARC_MAGIC_;*/ /* Defer this until we use it. */
	cdbBase[i].data = &dataBase[i*(uDataSize+1)+1];
	cdbBase[i].coll = NULL;
	Free->lru_insert(&cdbBase[i]);
    }
    /*
    ** Finally insert CDBs without data pages into the Free list last.
    */
    for (i=nCache;i<2*nCache;++i) {
	cdbBase[i].data = 0;
	cdbBase[i].coll = NULL;
	Free->mru_insert(&cdbBase[i]);
    }
}


ARC::~ARC() {
    /*
    ** Free the sentinels.
    */
    delete Free;
    delete B2;
    delete T2;
    delete B1;
    delete T1;
    /*
    ** Free the hash table.
    */
    delete Hash;
    /*
    ** Free the data pages.
    */
    delete dataBase;
    /*
    ** Free the CDBs.
    */
    delete cdbBase;
    }

/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct srvHeader {
    int32_t idFrom;      /* Global thread ID */
    int16_t replyTag;    /* Which MPI tag to send the reply */
    int16_t sid;
    int32_t nInBytes;
    int32_t nOutBytes;
    } SRVHEAD;

typedef struct {
    MDLserviceElement se;
    uint32_t iDataSize;
    } cacheOpenClose;


typedef struct {
    MDLserviceElement se;
    ptrdiff_t nz, sz, ny, sy, nLocal;
    int n1,n2,n3;
    } fftSizes;

#ifdef MDL_FFTW
typedef struct {
    fftSizes sizes;
    FFTW3(plan) fplan, iplan;
    FFTW3(real) *data;
    FFTW3(complex) *kdata;
    } fftPlans;

typedef struct {
    MDLserviceElement se;
    MDLFFT fft;
    FFTW3(real) *data;
    FFTW3(complex) *kdata;
    } fftTrans;
#endif

typedef struct {
    MDLserviceElement se;
    MDLGRID grid;
    } gridShare;

typedef struct {
    MDLserviceElement se;
    void *sbuff;
    int *scount;
    int *sdisps;
    void *rbuff;
    int *rcount;
    int *rdisps;
    int dataSize;
    } alltoallv;

/*****************************************************************************\

The following is the MPI thread only functions

\*****************************************************************************/

mpiClass::mpiClass(void (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),int argc, char **argv)
    : mdlClass(fcnMaster,fcnWorkerInit,fcnWorkerDone,argc,argv) {
    }

mpiClass::~mpiClass() {
    }

MDLflushBuffer * mpiClass::get_local_flush_buffer() {
    MDLflushBuffer *flush;

    /* This is bad. We just have to wait and hope for the best. */
    while (OPA_Queue_is_empty(&mpi->localFlushBuffers)) {
	/*
	** We have no local flush buffers, but threads should wakeup soon and give them back.
	** In the mean time we need to flush buffers ourselves (if we are not a dedicated MPI).
	*/
	bookkeeping();
	}
    OPA_Queue_dequeue(&mpi->localFlushBuffers, flush, MDLflushBuffer, hdr.hdr);
    flush->nBytes = 0;
    return flush;
    }

void mpiClass::queue_local_flush(CAHEAD *ph) {
    int iCore = ph->idTo;
    mdlClass * mdl1 = pmdl[iCore];
    CACHE *c = &mdl1->cache[ph->cid];
    MDLflushBuffer *flush;
    int nLeft, nToAdd;

    if ((flush=mpi->flushBuffersByCore[iCore]) == NULL) {
	flush = mpi->flushBuffersByCore[iCore] = get_local_flush_buffer();
	}

    nLeft = flush->nBufferSize - flush->nBytes;
    nToAdd = sizeof(CAHEAD) + c->iLineSize;

    /* No room left? Send this buffer and grab a new one for this rank */
    if (nToAdd > nLeft ) {
	OPA_Queue_enqueue(&mdl1->wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	flush = mpi->flushBuffersByCore[iCore] = get_local_flush_buffer();
	}

    /* Buffer this element */
    CAHEAD *ca = (CAHEAD *)( (char *)(flush+1) + flush->nBytes);
    char *pData = (char *)(ca+1);
    memcpy(ca,ph,sizeof(CAHEAD) + c->iLineSize);
    flush->nBytes += sizeof(CAHEAD) + c->iLineSize;
    }

#define MDL_MID_CACHEREQ	1
#define MDL_MID_CACHERPL	2
#define MDL_MID_CACHEFLUSH	3

int mpiClass::CacheReceive(MPI_Status *status) {
    CACHE *c;
    CAHEAD *ph = (CAHEAD *)(mpi->pReqRcv+1);
    char *pszRcv = (char *)(ph+1);
    MDLserviceCacheReq *creq;
    MDLcacheReplyData *pdata;
    CAHEAD *phRpl;
    char *pszRpl;
    char *t;
    int tag;
    int s,n,i;
    int ret, count;
    int iLineSize, iProc, iCore;
    uint32_t iIndex;
    char *pLine;

    mpi->pReqRcv->mpiRequest = MPI_REQUEST_NULL;
    MPI_Get_count(status, MPI_BYTE, &count);
    int iRankFrom = status->MPI_SOURCE;
    assert(count>=sizeof(CAHEAD));

    /* Well, could be any threads cache */
    iCore = ph->idTo - pmdl[0]->Self();
    assert(iCore>=0 && iCore<Cores());
    c = &pmdl[iCore]->cache[ph->cid];
    if (c->iType == MDL_NOCACHE) printf("Cache %d\n",ph->cid);

    assert(c->iType != MDL_NOCACHE);

    switch (ph->mid) {
	/* A different process wants local cache data - we can grab it directly */
    case MDL_MID_CACHEREQ:
	assert( count == sizeof(CAHEAD) );
	iProc = ThreadToProc(ph->idFrom); /* Can use iRankFrom */
	assert(iProc==iRankFrom);
	while (mpi->freeCacheReplies == NULL ) { /* Spin waiting for a buffer */
	    MDLcacheReplyData **plast, *pnext;
	    plast = &mpi->busyCacheReplies;
	    for(pdata = mpi->busyCacheReplies; pdata != NULL; pdata = pnext) {
		int flag;
		pnext = pdata->next;
		assert(pdata->mpiRequest!=MPI_REQUEST_NULL);
		MPI_Test(&pdata->mpiRequest,&flag,status);
		if (flag) {
		    pdata->mpiRequest = MPI_REQUEST_NULL;
		    if (pdata->next==NULL) mpi->busyCacheRepliesTail = plast;
		    *plast = pdata->next;
		    pdata->next = mpi->freeCacheReplies;
		    mpi->freeCacheReplies = pdata;
		    break;
		    }
		else plast = &pdata->next;
		}
	    /*mpi->busyCacheRepliesTail = plast;*/
	    }
	if (mpi->freeCacheReplies) {
	    pdata = mpi->freeCacheReplies;
	    mpi->freeCacheReplies = pdata->next;
	    }
	else {
	    pdata = CAST(MDLcacheReplyData *,malloc(sizeof(MDLcacheReplyData) + mpi->iReplyBufSize));
	    assert(pdata != NULL);
	    }
	assert(*mpi->busyCacheRepliesTail == NULL);
	pdata->next = NULL;
	*mpi->busyCacheRepliesTail = pdata;
	mpi->busyCacheRepliesTail = &pdata->next;

	phRpl = (CAHEAD *)(pdata + 1);
	pszRpl = (char *)(phRpl + 1);
	phRpl->cid = ph->cid;
	phRpl->mid = MDL_MID_CACHERPL;
	phRpl->idFrom = Self();

	/* Caclulation based on threads */
	phRpl->idTo = ph->idFrom;
	phRpl->iLine = ph->iLine;
	phRpl->nItems = 1;

	s = ph->iLine << c->nLineBits;
	n = s + c->nLineElements;
	iLineSize = c->iLineSize;
	for(i=s; i<n; i++ ) {
	    if (i<c->nData) {
		t = CAST(char *,(*c->getElt)(c->pData,i,c->iDataSize));
		memcpy(pszRpl,t,c->iDataSize);
		}
	    else memset(pszRpl,0,c->iDataSize);
	    pszRpl += c->iDataSize;
	    }
	tag = MDL_TAG_CACHECOM /*+ MDL_TAG_THREAD_OFFSET * iCore*/;
	assert(pdata->mpiRequest == MPI_REQUEST_NULL);
	MPI_Isend(phRpl,(int)sizeof(CAHEAD)+iLineSize,MPI_BYTE,
		  iRankFrom, tag, mpi->commMDL,&pdata->mpiRequest);
	/* This is found on the "busyCacheReplies" list */
	ret = 0;
	break;
    case MDL_MID_CACHEFLUSH:
	assert(c->iType == MDL_COCACHE);
	while(count>0) {
	    assert(count > sizeof(CAHEAD));
	    pszRcv = (char *)(ph+1);
	    iCore = ph->idTo - pmdl[0]->Self();
	    c = &pmdl[iCore]->cache[ph->cid];
	    iIndex = ph->iLine << c->nLineBits;
	    while(iIndex >= c->nData) {
		iIndex -= c->nData;
		assert(iCore+1<Cores());
		c = &pmdl[++iCore]->cache[ph->cid];
		}
	    ph->iLine = iIndex >> c->nLineBits;
	    ph->idTo = iCore;
	    queue_local_flush(ph);
	    pszRcv += c->iLineSize;
	    ph = (CAHEAD *)(pszRcv);
	    count -= sizeof(CAHEAD) + c->iLineSize;
	    }
	assert(count==0);
	ret = 0;
	break;
	/* A remote process has sent us cache data - we need to let that thread know */
    case MDL_MID_CACHERPL:
	iLineSize = c->iLineSize;
	assert( count == sizeof(CAHEAD) + iLineSize );
	creq = mpi->pThreadCacheReq[iCore];
	assert(creq!=NULL);
	mpi->pThreadCacheReq[iCore] = NULL;
	assert(creq->request!=MPI_REQUEST_NULL);
	MPI_Request_free(&creq->request);
	creq->request = MPI_REQUEST_NULL;
	pLine = CAST(char *,creq->pLine);
	creq->caReq.nItems = ph->nItems; /* Let caller know actual number */
	/*
	 ** For now assume no prefetching!
	 ** This means that this WILL be the reply to this Acquire
	 ** request.
	 */
	assert(pLine != NULL);
	memcpy(pLine,pszRcv,iLineSize);
	if (c->iType == MDL_COCACHE && c->init) {
	    /*
	     ** Call the initializer function for all elements in
	     ** the cache line.
	     */
	    for (i=0;i<iLineSize;i+=c->iDataSize) {
		    (*c->init)(c->ctx,&pLine[i]);
		}
	    }
	SendThreadMessage(MDL_TAG_CACHECOM,iCore,creq,MDL_SE_CACHE_REQUEST);
	ret = 1;
	break;
    default:
	assert(0);
	}

    /* Fire up next receive */
    assert(mpi->pReqRcv->mpiRequest == MPI_REQUEST_NULL);
    MPI_Irecv(mpi->pReqRcv+1,mpi->iCacheBufSize, MPI_BYTE, MPI_ANY_SOURCE,
	      MDL_TAG_CACHECOM, mpi->commMDL, &mpi->pReqRcv->mpiRequest);

    return ret;
    }

/*
** Perform a global swap: effectively a specialized in-place alltoallv
** - data to be sent starts at "buffer" and send items are contiguous and in order by target
** - total size of the buffer is "count" and must be at least as big as what we sent
** - free space immediately follows send data and there must be some free space.
*/
int mpiClass::swapall(const char *buffer,int count,int datasize,/*const*/ int *counts) {
    size_t size = datasize; /* NB!! must be a 64-bit type, hence size_t */
    char * const pBufferEnd = (char *)buffer + size*count;
    int *rcount, *scount;
    MPI_Datatype mpitype;
    MPI_Request *requests;
    MPI_Status *statuses;
    char *pSendBuff;
    int nTotalReceived = 0;
    int nSend, nRecv, iSend, iRecv, nMaxSend;
    int nRequests;
    int i;

    /* MPI buffers and datatype */
    MPI_Type_contiguous(datasize,MPI_BYTE,&mpitype);
    MPI_Type_commit(&mpitype);
    rcount = CAST(int *,malloc(sizeof(int) * Procs()));
    assert(rcount != NULL);
    scount = CAST(int *,malloc(sizeof(int) * Procs()));
    assert(scount != NULL);
    /* Note twice as many request/status pairs because of send + recv possibility */
    requests = CAST(MPI_Request *,malloc(sizeof(MPI_Request) * 2 * Procs()));
    assert(requests != NULL);
    statuses = CAST(MPI_Status *,malloc(sizeof(MPI_Status) * 2 * Procs()));
    assert(statuses != NULL);

    /* Counts of what to send and then move the "rejects" to the end of the buffer */
    for(nSend=0,i=0; i<Procs(); ++i) nSend += counts[i];
    assert(nSend <= count);
    nRecv = count - nSend;
    pSendBuff = (char *)buffer + size * nRecv;
    if (nRecv && nSend) memmove(pSendBuff,buffer,size * nSend);

    for(;;) {
	/*
	** At the start of this loop:
	**   nRecv: amount of element we have room to receive
	**   buffer: pointer to the receive buffer
	**   nSend: number of element (total) that we need to send
	**   pSendBuff: pointer to first element to send
	**   counts[]: number to send to each target: sum = nSend
	*/

	/* We are done when there is nothing globally left to send */
	MPI_Allreduce(&nSend,&nMaxSend,1,MPI_INT,MPI_MAX,mpi->commMDL);
	if (nMaxSend==0) break;

	/* Collect how many we need to receive */
	MPI_Alltoall(counts,1,MPI_INT,rcount,1,MPI_INT,mpi->commMDL);

	nRequests = 0;
	iRecv = 0;
	/* Calculate how many we can actually receive and post the receive */
	for(i=0; i<Procs(); ++i) {
	    if (iRecv < nRecv && rcount[i]) {
		if ( nRecv-iRecv < rcount[i] ) rcount[i] = nRecv-iRecv;
		MPI_Irecv((char *)buffer + size * iRecv, rcount[i], mpitype,
		    i, MDL_TAG_SWAP, mpi->commMDL, requests + nRequests++ );
		iRecv += rcount[i];
		}
	    else rcount[i] = 0;
	    }
	assert(iRecv <= nRecv);

	/* Now communicate how many we actually want to receive */
	MPI_Alltoall(rcount,1,MPI_INT,scount,1,MPI_INT,mpi->commMDL);

	/* Now we post the matching sends. We can use "rsend" here because the receives are posted. */
	iSend = 0;
	for(i=0; i<Procs(); ++i) {
	    if (scount[i]) {
		MPI_Irsend(pSendBuff + size * iSend, scount[i], mpitype,
		    i, MDL_TAG_SWAP, mpi->commMDL, requests + nRequests++ );
		}
	    iSend += counts[i]; /* note "counts" here, not "scount" */
	    }

	/* Wait for all communication (send and recv) to finish */
	MPI_Waitall(nRequests, requests, statuses);

	/* Now we compact the buffer for the next iteration */
	char *pSrc = pBufferEnd;
	pSendBuff = pBufferEnd;
	nRecv -= iRecv; /* We received "iRecv" element, so we have that must less room */
	for(i=Procs()-1; i>=0; --i) {
	    if (counts[i]) {
		size_t nLeft = counts[i] - scount[i];
		nRecv += scount[i]; /* We sent data so we have more receive room */
		nSend -= scount[i];
		if (nLeft) {
		    size_t nMove = size*nLeft;
		    if (pSrc != pSendBuff) memmove(pSendBuff - nMove, pSrc - nMove, nMove);
		    pSendBuff -= nMove;
		    }
		pSrc -= size * counts[i];
		counts[i] = nLeft;
		}
	    }
	nTotalReceived += iRecv;
	buffer += size * iRecv; /* New receive area */
	}

    MPI_Type_free(&mpitype);
    free(rcount);
    free(scount);
    free(requests);
    free(statuses);

    return nTotalReceived;
    }

/* Remove the given flush buffer from whichever list it is in */
static inline MDLflushBuffer *flush_remove(MDLflushBuffer *t) {
    t->hdr.mpi.prev->hdr.mpi.next = t->hdr.mpi.next;
    t->hdr.mpi.next->hdr.mpi.prev = t->hdr.mpi.prev;
    return t;
}

/* Insert element "e" after "list" (head or other element) */
static inline void flush_insert_after(MDLflushBuffer *list, MDLflushBuffer *e) {
    e->hdr.mpi.prev = list;
    e->hdr.mpi.next = list->hdr.mpi.next;
    list->hdr.mpi.next->hdr.mpi.prev = e;
    list->hdr.mpi.next = e;
}

/* Check all pending sends to see if they have completed so we can reuse the buffer */
int mpiClass::flush_check_completion() {
    MDLflushBuffer *pBuffer, *pNext;

    for(pBuffer=mpi->flushHeadSent.hdr.mpi.next; pBuffer != &mpi->flushHeadSent; pBuffer=pNext) {
	int flag;
	pNext = pBuffer->hdr.mpi.next;
	MPI_Test(&pBuffer->request, &flag, MPI_STATUS_IGNORE);
	if (flag) {
	    flush_insert_after(&mpi->flushHeadFree,flush_remove(pBuffer));
	    }
	}
    return mpi->flushHeadSent.hdr.mpi.next != &mpi->flushHeadSent; /* true == still waiting */
    }

void mpiClass::flush_element(CAHEAD *pHdr,int iLineSize) {
    int iProc = ThreadToProc(pHdr->idTo);
    MDLflushBuffer *pBuffer;
    int iCore;

    /* If we have a buffer already then add to it if we can */
    if ( (pBuffer=mpi->flushBuffersByRank[iProc]) != NULL ) {
	int nLeft = pBuffer->nBufferSize - pBuffer->nBytes;
	int nToAdd = sizeof(CAHEAD) + iLineSize;
	/* No room left? Send this buffer and grab a new one for this rank */
	if (nToAdd > nLeft ) {
	    flush_send(pBuffer);
	    pBuffer = mpi->flushBuffersByRank[iProc] = NULL;
	    }
	else {
	    /* We have used this buffer, so let's make it most recently used. */
	    flush_insert_after(&mpi->flushHeadBusy,flush_remove(pBuffer));
	    }
	}
    /* Grab a free buffer, or wait for one to become free */
    if (pBuffer==NULL) {
	while (mpi->flushHeadFree.hdr.mpi.next == &mpi->flushHeadFree) {
	    MPI_Status status;
	    int flag;
	    flush_check_completion();
    	    if (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL) {
            	while(1) {
                    MPI_Test(&mpi->pReqRcv->mpiRequest, &flag, &status);
                    if (flag) CacheReceive(&status);
                    else break;
                    }
                }
	    }
	pBuffer = flush_remove(mpi->flushHeadFree.hdr.mpi.next);
	flush_insert_after(&mpi->flushHeadBusy,pBuffer);
	mpi->flushBusyCount++;
	mpi->flushBuffersByRank[iProc] = pBuffer;
	pBuffer->nBytes = 0;
	pBuffer->iRankTo = iProc;
	/* If we have run out of free buffers then flush the oldest right away */
//	if (mpi->flushHeadFree.hdr.mpi.next == &mpi->flushHeadFree ) {
	if (mpi->flushBusyCount >= mpi->flushBuffCount-32) {
	    assert(mpi->flushBusyCount>=0);
	    assert(mpi->flushBusyCount<=mpi->flushBuffCount);
	    if (mpi->flushHeadBusy.hdr.mpi.next != &mpi->flushHeadBusy ) {
		MDLflushBuffer *flush = mpi->flushHeadBusy.hdr.mpi.prev;
		mpi->flushBuffersByRank[flush->iRankTo] = NULL;
		flush_send(flush);
		}
	    }
	}

    /* Now just buffer the data */
    CAHEAD *ca = (CAHEAD *)( (char *)(pBuffer+1) + pBuffer->nBytes);
    char *pData = (char *)(ca+1);
    iLineSize += sizeof(CAHEAD);
    memcpy(ca,pHdr,iLineSize);
    pBuffer->nBytes += iLineSize;
    }

void mpiClass::flush_elements(MDLflushBuffer *pFlush) {
    CAHEAD *ca;
    char *pData;
    int count;

    count = pFlush->nBytes;
    ca = (CAHEAD *)(pFlush + 1);
    while(count > 0) {
	int iLineSize = pmdl[0]->cache[ca->cid].iLineSize;
	pData = (char *)(ca+1);
	flush_element(ca,iLineSize);
	pData += iLineSize;
	ca = (CAHEAD *)pData;
	count -= sizeof(CAHEAD) + iLineSize;
	}
    assert(count == 0);
    }

/*
** At the end we need to flush any pending cache elements
*/
void mpiClass::flush_all_elements() {
    while (mpi->flushHeadBusy.hdr.mpi.next != &mpi->flushHeadBusy ) {
	MDLflushBuffer *flush = mpi->flushHeadBusy.hdr.mpi.prev;
	mpi->flushBuffersByRank[flush->iRankTo] = NULL;
	flush_send(flush);
	}
    while(flush_check_completion()) {
	checkMPI();
	}
    }

void mpiClass::finish_local_flush() {
    MDLflushBuffer *flush;
    int iCore;

    for(iCore=0; iCore<Cores(); ++iCore) {
	mdlClass * mdl1 = pmdl[iCore];
	if ((flush=mpi->flushBuffersByCore[iCore])) {
	    OPA_Queue_enqueue(&mdl1->wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	    mpi->flushBuffersByCore[iCore] = NULL;
	    }
	}
    }

/* Send the given buffer. Remove it from whichever list and put it on the "sent" list. */
void mpiClass::flush_send(MDLflushBuffer *pBuffer) {
    if (pBuffer->nBytes>0) {
	MPI_Issend(pBuffer+1,pBuffer->nBytes,MPI_BYTE,pBuffer->iRankTo,
	    MDL_TAG_CACHECOM,mpi->commMDL,&pBuffer->request);
	flush_insert_after(mpi->flushHeadSent.hdr.mpi.prev,flush_remove(pBuffer));
	}
    else {
	flush_insert_after(&mpi->flushHeadFree,flush_remove(pBuffer));
	}
    mpi->flushBusyCount--;
    }

/*
** This routine must be called often by the MPI thread. It will drain
** any requests from the thread queue, and track MPI completions.
*/
int mpiClass::checkMPI() {
    while(1) {
	MDLserviceSend *send;
	MDLserviceCacheReq *caReq;
	MDLserviceElement *qhdr;
	MDLflushBuffer *pFlush;
	cacheOpenClose *coc;
#ifdef MDL_FFTW
	fftSizes* sizes;
	fftPlans* plans;
	fftTrans* trans;
#endif
	alltoallv *a2a;
	gridShare *share;
	SRVHEAD *head;
	int iProc, iCore, tag, i;
	int flag,indx;
	MPI_Status status;
	MPI_Datatype mpitype;

	/* These are messages from other threads */
	if (!OPA_Queue_is_empty(&mpi->queueMPI)) {
	    OPA_Queue_dequeue(&mpi->queueMPI, qhdr, MDLserviceElement, hdr);
	    switch(qhdr->iServiceID) {
		/* A thread has finished */
	    case MDL_SE_STOP:
		--mpi->nActiveCores;
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_STOP);
		break;
		/* mdlReqService() on behalf of a thread */
	    case MDL_SE_SEND_REQUEST:
		send = (MDLserviceSend *)qhdr;
		iProc = ThreadToProc(send->target);
		iCore = send->target - ProcToThread(iProc);
		assert(iCore>=0);
		tag = send->tag + MDL_TAG_THREAD_OFFSET * iCore;

		/* Grab a free tag for the reply */
		head = CAST(SRVHEAD *,send->buf);
		i = mpi->iRequestTarget;
		do {
		    if (++mpi->iRequestTarget == mpi->nRequestTargets) mpi->iRequestTarget = 0;
		    } while(i != mpi->iRequestTarget && mpi->pRequestTargets[mpi->iRequestTarget] >= 0);
		assert(mpi->pRequestTargets[mpi->iRequestTarget] < 0);
		head->replyTag = mpi->iRequestTarget;
		mpi->pRequestTargets[mpi->iRequestTarget] = iProc;
		MPI_Isend(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;
		/* A service has finished and wants to send back the reply */
	    case MDL_SE_SEND_REPLY:
		send = (MDLserviceSend *)qhdr;
		iProc = ThreadToProc(send->target);
		/* tag is really the request ID */
		head = CAST(SRVHEAD *,send->buf);
		tag = MDL_TAG_RPL + MDL_TAG_THREAD_OFFSET * head->replyTag;
		MPI_Isend(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;
		/* mdlGetReply() on behalf of a thread */
	    case MDL_SE_RECV_REPLY:
		send = (MDLserviceSend *)qhdr;
		/* Target is really the request ID */
		assert(send->target < mpi->nRequestTargets);
		iProc = mpi->pRequestTargets[send->target];
		tag = MDL_TAG_RPL + MDL_TAG_THREAD_OFFSET * send->target;
		head = CAST(SRVHEAD *,send->buf);
		MPI_Irecv(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;

	    case MDL_SE_MPI_SEND:
	    case MDL_SE_MPI_SSEND:
		send = (MDLserviceSend *)qhdr;
		iProc = ThreadToProc(send->target);
		iCore = send->target - ProcToThread(iProc);
		assert(iCore>=0);
		tag = send->tag + MDL_TAG_THREAD_OFFSET * iCore;
		if (qhdr->iServiceID==MDL_SE_MPI_SSEND)
		    MPI_Issend(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		else
		    MPI_Isend(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;
	    case MDL_SE_MPI_RECV:
	    case MDL_SE_RECV_REQUEST:
		send = (MDLserviceSend *)qhdr;
		if (send->target==MPI_ANY_SOURCE) iProc = MPI_ANY_SOURCE;
		else iProc = ThreadToProc(send->target);
		iCore = send->svc.iCoreFrom;
		assert(iCore>=0);
		tag = send->tag + MDL_TAG_THREAD_OFFSET * iCore;
		MPI_Irecv(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;
		/* A thread has opened a cache -- we may need to resize our buffers */
	    case MDL_SE_CACHE_OPEN:
		coc = (cacheOpenClose *)qhdr;
		if (mpi->pReqRcv->mpiRequest == MPI_REQUEST_NULL) {
		    MPI_Irecv(mpi->pReqRcv+1,mpi->iCacheBufSize, MPI_BYTE,
			MPI_ANY_SOURCE, MDL_TAG_CACHECOM,
			mpi->commMDL, &mpi->pReqRcv->mpiRequest);
		    }
		++mpi->nOpenCaches;
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_OPEN);
		break;
		/* A thread has closed a cache -- we can cancel the receive if no threads have a cache open */
	    case MDL_SE_CACHE_CLOSE:
		assert(mpi->nOpenCaches > 0);
		--mpi->nOpenCaches;
		assert (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL);
		if (mpi->nOpenCaches == 0) {
		    flush_all_elements();
		    MPI_Cancel(&mpi->pReqRcv->mpiRequest);
		    MPI_Wait(&mpi->pReqRcv->mpiRequest, &status);
		    mpi->pReqRcv->mpiRequest = MPI_REQUEST_NULL;
		    }
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
		/* A thread needs to request a missing element */
	    case MDL_SE_CACHE_REQUEST:
		caReq = (MDLserviceCacheReq *)qhdr;
		assert(mpi->pThreadCacheReq[caReq->svc.iCoreFrom]==NULL);
		iProc = ThreadToProc(caReq->caReq.idTo);
		mpi->pThreadCacheReq[caReq->svc.iCoreFrom] = caReq;
		assert(caReq->request==MPI_REQUEST_NULL);
		MPI_Isend(&caReq->caReq,sizeof(CAHEAD),MPI_BYTE,iProc,MDL_TAG_CACHECOM, mpi->commMDL,&caReq->request);
		break;
	    case MDL_SE_CACHE_FLUSH:
		pFlush = (MDLflushBuffer *)qhdr;
		flush_elements(pFlush);
		OPA_Queue_enqueue(&pmdl[qhdr->iCoreFrom]->coreFlushBuffers, pFlush, MDLflushBuffer, hdr.hdr);
		break;
	    case MDL_SE_CACHE_FLUSH_OUT:
		flush_all_elements();
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_FLUSH_OUT);
		break;
	    case MDL_SE_CACHE_FLUSH_LCL:
		finish_local_flush();
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_FLUSH_LCL);
		break;
#ifdef MDL_FFTW
	    case MDL_SE_FFT_SIZES:
		sizes = (fftSizes *)qhdr;
		sizes->nLocal = 2*FFTW3(mpi_local_size_3d_transposed)
		    (sizes->n3,sizes->n2,sizes->n1/2+1,mpi->commMDL,
		    &sizes->nz,&sizes->sz,&sizes->ny,&sizes->sy);
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
	    case MDL_SE_FFT_PLANS:
		plans = (fftPlans *)qhdr;
		plans->sizes.nLocal = 2*FFTW3(mpi_local_size_3d_transposed)
		    (plans->sizes.n3,plans->sizes.n2,plans->sizes.n1/2+1,mpi->commMDL,
		    &plans->sizes.nz,&plans->sizes.sz,&plans->sizes.ny,&plans->sizes.sy);
		plans->fplan = FFTW3(mpi_plan_dft_r2c_3d)(
		    plans->sizes.n3,plans->sizes.n2,plans->sizes.n1,plans->data,plans->kdata,
		    mpi->commMDL,FFTW_MPI_TRANSPOSED_OUT | (plans->data==NULL?FFTW_ESTIMATE:FFTW_MEASURE) );
		plans->iplan = FFTW3(mpi_plan_dft_c2r_3d)(
		    plans->sizes.n3,plans->sizes.n2,plans->sizes.n1,plans->kdata,plans->data,
		    mpi->commMDL,FFTW_MPI_TRANSPOSED_IN  | (plans->kdata==NULL?FFTW_ESTIMATE:FFTW_MEASURE) );
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
	    case MDL_SE_FFT_DFT_R2C:
		trans = (fftTrans *)qhdr;
		FFTW3(execute_dft_r2c)(trans->fft->fplan,trans->data,trans->kdata);
		pthread_barrier_wait(&pmdl[0]->barrier);
		break;
	    case MDL_SE_FFT_DFT_C2R:
		trans = (fftTrans *)qhdr;
		FFTW3(execute_dft_c2r)(trans->fft->iplan,trans->kdata,trans->data);
		pthread_barrier_wait(&pmdl[0]->barrier);
		break;
#endif
	    case MDL_SE_GRID_SHARE:
		share = (gridShare *)qhdr;
		MPI_Allgather(&share->grid->sSlab,sizeof(*share->grid->rs),MPI_BYTE,
		    share->grid->rs,sizeof(*share->grid->rs),MPI_BYTE,
		    mpi->commMDL);
		MPI_Allgather(&share->grid->nSlab,sizeof(*share->grid->rn),MPI_BYTE,
		    share->grid->rn,sizeof(*share->grid->rn),MPI_BYTE,
		    mpi->commMDL);
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
	    case MDL_SE_BARRIER_REQUEST:
		MPI_Barrier(mpi->commMDL);
		SendThreadMessage(MDL_TAG_BARRIER,0,qhdr,MDL_SE_BARRIER_REPLY);
		break;
	    case MDL_SE_ALLTOALLV:
		a2a = (alltoallv *)qhdr;
		MPI_Type_contiguous(a2a->dataSize,MPI_BYTE,&mpitype);
		MPI_Type_commit(&mpitype);
		MPI_Alltoallv(
		    a2a->sbuff, a2a->scount, a2a->sdisps, mpitype,
		    a2a->rbuff, a2a->rcount, a2a->rdisps, mpitype,
		    mpi->commMDL);
		MPI_Type_free(&mpitype);
		SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_ALLTOALLV);
		break;
	    default:
		assert(0);
		}
	    continue;
	    }

	/* Start any CUDA work packages */
#ifdef USE_CUDA
	if (cudaCtx) {
	    CUDA_registerBuffers(mpi->cudaCtx, &mpi->queueREGISTER);
            CUDA_startWork(mpi->cudaCtx, &mpi->queueWORK);
	    CUDA_flushDone(mpi->cudaCtx);
	    }
#endif

	/* These are messages/completions from/to other MPI processes. */
	if (mpi->nSendRecvReq) {
	    MPI_Testany(mpi->nSendRecvReq,mpi->pSendRecvReq, &indx, &flag, &status);
	    if (flag) {
		assert(indx>=0 && indx<mpi->nSendRecvReq);
		send = mpi->pSendRecvBuf[indx];
		for(i=indx+1;i<mpi->nSendRecvReq;++i) {
		    mpi->pSendRecvReq[i-1] = mpi->pSendRecvReq[i];
		    mpi->pSendRecvBuf[i-1] = mpi->pSendRecvBuf[i];
		    }
		--mpi->nSendRecvReq;
		switch(send->svc.iServiceID) {
		case MDL_SE_MPI_SEND:
		case MDL_SE_MPI_SSEND:
		case MDL_SE_SEND_REQUEST:
		case MDL_SE_SEND_REPLY:
		    assert(send == &pmdl[send->svc.iCoreFrom]->sendRequest);
		    tag = 0;
		    break;
		case MDL_SE_RECV_REPLY:
		    head = CAST(SRVHEAD *,send->buf);
		    assert(mpi->pRequestTargets[send->target]>=0);
		    mpi->pRequestTargets[send->target] = -1;
		case MDL_SE_MPI_RECV:
		case MDL_SE_RECV_REQUEST:
		    assert(send == &pmdl[send->svc.iCoreFrom]->recvRequest);
		    tag = send->tag % MDL_TAG_THREAD_OFFSET;
		    MPI_Get_count(&status, MPI_BYTE, &send->count);
		    send->target = status.MPI_SOURCE;
		    break;
		default:
		    assert(0);
		    }
		SendThreadMessage(tag,send->svc.iCoreFrom,send,MDL_SE_MPI_SSEND);
		continue;
		}
	    }
	if (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL) {
	    while(1) {
		MPI_Test(&mpi->pReqRcv->mpiRequest, &flag, &status);
		if (flag) CacheReceive(&status);
		else break;
		}
	    }
	break;
	}
    return mpi->nActiveCores;
    }

static void TERM_handler(int signo) {
    MPI_Abort(MPI_COMM_WORLD,130);
    }

extern "C"
void mdlLaunch(int argc,char **argv,void (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *)) {
    mpiClass *mdl = new mpiClass(fcnMaster,fcnWorkerInit,fcnWorkerDone,argc,argv);
    mdl->Launch(argc,argv,fcnMaster,fcnWorkerInit,fcnWorkerDone);
    delete mdl;
    }
void mpiClass::Launch(int argc,char **argv,void (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *)) {
    int i,n,bDiag,bThreads,bDedicated,thread_support,rc,flag,*piTagUB;
    char *p, ach[256];
#ifdef USE_HWLOC
    hwloc_topology_t topology;
    hwloc_cpuset_t set_proc, set_thread;
    hwloc_obj_t t = NULL;
#endif

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("MDL Startup");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    mpi = new mdlContextMPI;

#ifdef _SC_NPROCESSORS_CONF /* from unistd.h */
    nCores = sysconf(_SC_NPROCESSORS_CONF);
#endif

    /*
    ** Do some low level argument parsing for number of threads, and
    ** diagnostic flag!
    */
    bDiag = 0;
    bThreads = 0;
    bDedicated = -1;
    if(argv) {
	for (argc = 0; argv[argc]; argc++);
	i = 1;
	while (argv[i]) {
	    if (!strcmp(argv[i], "-sz") && !bThreads) {
		++i;
		nCores = atoi(argv[i]);
		if (argv[i]) bThreads = 1;
		}
	    if (!strcmp(argv[i], "-dedicated")) {
		if (bDedicated<1) bDedicated = 0;
		}
	    if (!strcmp(argv[i], "+dedicated")) {
		if (bDedicated<1) bDedicated = 1;
		}
	    if (!strcmp(argv[i], "+sharedmpi")) {
		bDedicated = 2;
		}
	    if (!strcmp(argv[i], "+d") && !bDiag) {
		p = getenv("MDL_DIAGNOSTIC");
		if (!p) p = getenv("HOME");
		if (!p) sprintf(ach, "/tmp");
		else sprintf(ach, "%s", p);
		bDiag = 1;
		}
	    ++i;
	    }
	argc = i;
	}
    if (!bThreads) {
	if ( (p=getenv("SLURM_CPUS_PER_TASK")) != NULL ) nCores = atoi(p);
	else if ( (p=getenv("OMP_NUM_THREADS")) != NULL ) nCores = atoi(p);
	}
    assert(Cores()>0);

    /* MPI Initialization */
#ifdef USE_ITT
    __itt_string_handle* shMPITask = __itt_string_handle_create("MPI");
    __itt_task_begin(domain, __itt_null, __itt_null, shMPITask);
#endif
    mpi->commMDL = MPI_COMM_WORLD;
    rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,&thread_support);
    if (rc!=MPI_SUCCESS) {
	MPI_Error_string(rc, ach, &i);
	perror(ach);
	MPI_Abort(mpi->commMDL,rc);
	}
#ifdef HAVE_SIGNAL_H
    signal(SIGINT,TERM_handler);
#endif
#ifdef MDL_FFTW
    if (Cores()>1) FFTW3(init_threads)();
    FFTW3(mpi_init)();
    if (Cores()>1) FFTW3(plan_with_nthreads)(Cores());
#endif

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
    MPI_Comm_size(mpi->commMDL, &nProcs);
    MPI_Comm_rank(mpi->commMDL, &iProc);

    /* Dedicate one of the threads for MPI, unless it would be senseless to do so */
    if (bDedicated == -1) {
	if (nProcs>0 && Cores()>3) bDedicated = 1;
	else bDedicated = 0;
	}
    if (bDedicated == 1) {
#ifndef USE_HWLOC
	/*if (Cores()==1 || nProcs==1) bDedicated=0;
	  else*/ --Cores();
#else
	hwloc_topology_init(&topology);
	hwloc_topology_load(topology);
	set_proc = hwloc_bitmap_alloc();
	set_thread = hwloc_bitmap_alloc();
	n = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
	if (n != Cores()) n = 1;
	else if (hwloc_get_cpubind(topology,set_proc,HWLOC_CPUBIND_PROCESS) != 0) n = 1;
	else {
	    // Get the first core and reserve it for our MDL use
	    t = hwloc_get_obj_inside_cpuset_by_type(topology,set_proc,HWLOC_OBJ_CORE,0);
	    if ( t != NULL ) {
		hwloc_bitmap_andnot(set_proc,set_proc,t->cpuset);
		n = hwloc_bitmap_weight(t->cpuset);
		}
	    else n = 1;
	    }
	nCores -= n;
#endif
	}

    /* Construct the thread/processor map */
    iProcToThread = CAST(int *,malloc((Procs() + 1) * sizeof(int)));
    assert(iProcToThread != NULL);
    iProcToThread[0] = 0;
    MPI_Allgather(&nCores, 1, MPI_INT, iProcToThread + 1, 1, MPI_INT, mpi->commMDL);
    for (i = 1; i < Procs(); ++i) iProcToThread[i + 1] += iProcToThread[i];
    nThreads = iProcToThread[Procs()];
    idSelf = iProcToThread[Proc()];

    iCoreMPI = bDedicated ? -1 : 0;
    iCore = iCoreMPI;

    /* We have an extra MDL context in case we want a dedicated MPI thread */
    pmdl = new mdlClass *[Cores()+1];
    assert(pmdl!=NULL);
    pmdl++;
    pmdl[iCoreMPI] = this;
    threadid = CAST(pthread_t *,malloc(Cores() * sizeof(pthread_t)));

#ifdef USE_CL
    void * clContext = CL_create_context();
#endif

    /* Allocate the other MDL structures for any threads. */
    for (i = iCoreMPI+1; i < Cores(); ++i) {
	//pmdl[i] = new mdlClass(fcnMaster,fcnWorkerInit,fcnWorkerDone,argc,argv);
	pmdl[i] = new mdlClass(this,i);
	}

/* All GPU work is funneled through the MPI thread */
#ifdef USE_CUDA
    mpi->inCudaBufSize = mpi->outCudaBufSize = 0;
    mpi->cudaCtx = CUDA_initialize(0,-1,NULL,NULL);
    OPA_Queue_init(&mpi->queueCUDA);
#endif

    OPA_Queue_init(&mpi->queueMPI);
    OPA_Queue_init(&mpi->queueREGISTER);
    OPA_Queue_init(&mpi->queueWORK);
    mpi->nSendRecvReq = 0;
    mpi->pSendRecvReq = NULL;
    mpi->pSendRecvBuf = NULL;
    mpi->pThreadCacheReq = NULL;
    mpi->pRequestTargets = NULL;
    mpi->nRequestTargets = 0;
    mpi->iRequestTarget = 0;

    assert(sizeof(CAHEAD) == 16); /* Well, should be a multiple of 8 at least. */
    mpi->nOpenCaches = 0;
    mpi->iReplyBufSize = sizeof(CAHEAD) + MDL_CACHE_DATA_SIZE;
    mpi->iCacheBufSize = sizeof(CAHEAD) + MDL_FLUSH_DATA_SIZE;
    mpi->pReqRcv = CAST(MDLcacheReplyData *,malloc(sizeof(MDLcacheReplyData) + mpi->iCacheBufSize));
    assert(mpi->pReqRcv != NULL);
    mpi->pReqRcv->mpiRequest = MPI_REQUEST_NULL;

    /* Start with a sensible number of cache buffers */
    mpi->freeCacheReplies = NULL;
    mpi->busyCacheReplies = NULL;
    mpi->busyCacheRepliesTail = &mpi->busyCacheReplies;
    n = Procs();
    if (n > 256) n = 256;
    else if (n < 64) n = 64;
    for (i = 0; i<n; ++i) {
	MDLcacheReplyData *pdata = CAST(MDLcacheReplyData *,malloc(sizeof(MDLcacheReplyData) + mpi->iReplyBufSize));
	assert( pdata != NULL );
	pdata->mpiRequest = MPI_REQUEST_NULL;
	pdata->next = mpi->freeCacheReplies;
	mpi->freeCacheReplies = pdata;
	}

    /* For flushing to remote processes/nodes */
    mpi->flushBuffersByRank = CAST(MDLflushBuffer **,malloc(sizeof(MDLflushBuffer *) * Procs()));
    assert(mpi->flushBuffersByRank);
    for(i=0; i<Procs(); ++i) mpi->flushBuffersByRank[i] = NULL;

    /* For flushing to threads on this processor */
    OPA_Queue_init(&mpi->localFlushBuffers);
    mpi->flushBuffersByCore = CAST(MDLflushBuffer **,malloc(sizeof(MDLflushBuffer *) * Cores()));
    assert(mpi->flushBuffersByCore);
    for(i=0; i<Cores(); ++i) mpi->flushBuffersByCore[i] = NULL;
    for(i=0; i<2*Cores(); ++i) {
	MDLflushBuffer *pBuffer = CAST(MDLflushBuffer *,malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE));
	assert(pBuffer!=NULL);
	pBuffer->nBufferSize = MDL_FLUSH_DATA_SIZE;
	pBuffer->nBytes = 0;
	OPA_Queue_enqueue(&mpi->localFlushBuffers, pBuffer, MDLflushBuffer, hdr.hdr);
	}

    mpi->flushHeadBusy.hdr.mpi.prev = mpi->flushHeadBusy.hdr.mpi.next = &mpi->flushHeadBusy;
    mpi->flushHeadFree.hdr.mpi.prev = mpi->flushHeadFree.hdr.mpi.next = &mpi->flushHeadFree;
    mpi->flushHeadSent.hdr.mpi.prev = mpi->flushHeadSent.hdr.mpi.next = &mpi->flushHeadSent;

    mpi->flushBuffCount = n;
    mpi->flushBusyCount = 0;
    for(i=0; i<n; ++i) {
	MDLflushBuffer *pBuffer = CAST(MDLflushBuffer *,malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE));
	assert(pBuffer!=NULL);
	pBuffer->nBufferSize = MDL_FLUSH_DATA_SIZE;
	pBuffer->nBytes = 0;
	flush_insert_after(mpi->flushHeadFree.hdr.mpi.prev,pBuffer);
	}

    /* Some bookeeping for the send/recv - 1 of each per thread */
    mpi->nSendRecvReq = 0;
    mpi->pSendRecvReq = CAST(MPI_Request *,malloc(Cores()*2*sizeof(MPI_Request)));
    assert(mpi->pSendRecvReq!=NULL);
    mpi->pSendRecvBuf = CAST(MDLserviceSend **,malloc(Cores()*2*sizeof(MDLserviceSend *)));
    assert(mpi->pSendRecvBuf!=NULL);
    mpi->pThreadCacheReq = CAST(MDLserviceCacheReq **,malloc(Cores()*sizeof(MDLserviceCacheReq *)));
    assert(mpi->pThreadCacheReq!=NULL);
    for (i = 0; i < Cores(); ++i) mpi->pThreadCacheReq[i] = NULL;

    /* Ring buffer of requests */
    mpi->iRequestTarget = 0;
    mpi->nRequestTargets = 2 * log2(1.0 * Threads()) + Cores();
    mpi->pRequestTargets = CAST(int *,malloc(mpi->nRequestTargets * sizeof(*mpi->pRequestTargets)));
    assert(mpi->pRequestTargets!=NULL);
    for(i=0; i<mpi->nRequestTargets; ++i) mpi->pRequestTargets[i] = -1;

    /* Make sure that MPI supports enough tags */
    rc = MPI_Comm_get_attr(mpi->commMDL,MPI_TAG_UB,&piTagUB,&flag);
    if (rc==MPI_SUCCESS && flag) {
	assert(Cores()*MDL_TAG_THREAD_OFFSET < *piTagUB);
	assert(mpi->nRequestTargets*MDL_TAG_THREAD_OFFSET < *piTagUB);
	}

#ifdef USE_ITT
    __itt_thread_set_name("MPI");
#endif

    /* Launch threads: if dedicated MPI thread then launch all worker threads. */
#ifdef USE_ITT
    __itt_string_handle* shPthreadTask = __itt_string_handle_create("pthread");
    __itt_task_begin(domain, __itt_null, __itt_null, shPthreadTask);
#endif
    mpi->nActiveCores = 0;
    if (Cores() > 1 || bDedicated) {
#ifdef USE_HWLOC
	int icpu = -1;
#endif
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_barrier_init(&pmdl[0]->barrier,NULL,Cores()+(bDedicated?1:0));
	for (i = iCoreMPI+1; i < Cores(); ++i) {
	    pthread_create(&threadid[i], &attr,
		mdlWorkerThread,
		pmdl[i]);
#ifdef USE_HWLOC
	    if (t) { // If requested, bind the new thread to the specific core
		icpu = hwloc_bitmap_next(set_proc,icpu);
		hwloc_bitmap_only(set_thread,icpu);
		hwloc_set_thread_cpubind(topology,threadid[i],set_thread,HWLOC_CPUBIND_THREAD);
		}
#endif
	    ++mpi->nActiveCores;
	    }
	pthread_attr_destroy(&attr);
	}
#ifdef USE_ITT
    __itt_task_end(domain);
    __itt_task_end(domain);
#endif
#ifdef USE_HWLOC
    hwloc_bitmap_free(set_thread);
    hwloc_bitmap_free(set_proc);
    hwloc_topology_destroy(topology);
#endif
    if (!bDedicated) {
	worker_ctx = (*fcnWorkerInit)(reinterpret_cast<MDL>(static_cast<mdlClass *>(this)));
	CommitServices();
	if (Self()) Handler();
	else run_master();
	(*fcnWorkerDone)(reinterpret_cast<MDL>(static_cast<mdlClass *>(this)),worker_ctx);
	}
    drainMPI();
    pthread_barrier_destroy(&pmdl[0]->barrier);
    }

void mdlAbort(MDL mdl) {
    abort();
    }

mdlClass::~mdlClass() {
    MDLcacheReplyData *pdata, *pnext;
    int i;

    for (i = iCoreMPI+1; i < Cores(); ++i) {
	pthread_join(threadid[i],0);
	pmdl[i]->cleanupMDL();
	free(pmdl[i]);
	}
    free(threadid);
    delete (pmdl-1);

    /* Finish any outstanding cache sends */
    for(pdata = mpi->busyCacheReplies; pdata != NULL; pdata=pnext) {
	MPI_Status status;
	assert(pdata->mpiRequest!=MPI_REQUEST_NULL);
	pnext = pdata->next;
	MPI_Wait(&pdata->mpiRequest,&status);
	pdata->next = mpi->freeCacheReplies;
	mpi->freeCacheReplies = pdata;
	}
    mpi->busyCacheReplies = NULL;
    mpi->busyCacheRepliesTail = &mpi->busyCacheReplies;
    MPI_Barrier(mpi->commMDL);
    MPI_Finalize();

#ifdef USE_CUDA
    if (cudaCtx) CUDA_finish(cudaCtx);
#endif

    cleanupMDL();

    /*
     ** Close Diagnostic file.
     */
    if (bDiag) {
	fclose(fpDiag);
	}

    /*
     ** Deallocate storage.
     */
    for(pdata = mpi->freeCacheReplies; pdata != NULL; pdata=pnext) {
	pnext = pdata->next;
	free(pdata);
	}
    mpi->freeCacheReplies = NULL;

#ifdef USE_CUDA
    if (mpi->cudaCtx) CUDA_finish(mpi->cudaCtx);
#endif
    free(mpi->pReqRcv);
    free(iProcToThread);
    free(mpi->pSendRecvReq);
    free(mpi->pSendRecvBuf);
    free(mpi->pThreadCacheReq);
    free(mpi->pRequestTargets);
    free(mpi);
    }

/*****************************************************************************\
\*****************************************************************************/

// This routine is overridden for the MPI thread.
int mdlClass::checkMPI() { return 0; }

ARC *mdlClass::arcReinitialize(ARC *arc,uint32_t nCache,uint32_t uDataSize,CACHE *c) {
    if (arc!=NULL) {
	assert(arc->cache == c);
	if ( arc->nCache == nCache && arc->uDataSize == ((uDataSize+7)>>3) )
	    return arc;
	delete arc;
	}
    return new ARC(this,nCache,uDataSize,c);
    }

void mdlClass::SendThreadMessage(int iQueue,int iCore, void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = Core();
    OPA_Queue_enqueue(pmdl[iCore]->inQueue+iQueue, qhdr, MDLserviceElement, hdr);
    }

void mdlClass::SendToMPI(void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = Core();
//    OPA_Queue_enqueue(&pmdl[iCoreMPI]->mpi->queueMPI, qhdr, MDLserviceElement, hdr);
    OPA_Queue_enqueue(&mpi->queueMPI, qhdr, MDLserviceElement, hdr);
    }

/* Accept pending combine requests, and call the combine function for each. */
void mdlClass::combine_all_incoming() {
    //mdlContextMPI *mpi = pmdl[iCoreMPI]->mpi;
    if (Core()<0) return; /* Dedicated MPI thread combines nothing */
    while (!OPA_Queue_is_empty(&wqCacheFlush)) {
	MDLflushBuffer *flush;
	CACHE *c;
	CAHEAD *ca;
	char *pData;
	int count, i, s, n;
	uint32_t uIndex;

	OPA_Queue_dequeue(&wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	count = flush->nBytes;
    	ca = (CAHEAD *)(flush + 1);
	while(count > 0) {
	    pData = (char *)(ca+1);
	    c = &cache[ca->cid];
	    uIndex = ca->iLine << c->nLineBits;
	    s = uIndex;
	    n = s + c->nLineElements;
	    for(i=s; i<n; i++ ) {
		if (i<c->nData)
		    (*c->combine)(c->ctx,(*c->getElt)(c->pData,i,c->iDataSize),pData);
		pData += c->iDataSize;
		}
	    count -= sizeof(CAHEAD) + c->iLineSize;
	    ca = (CAHEAD *)pData;
	    }
	assert(count == 0);
	OPA_Queue_enqueue(&mpi->localFlushBuffers, flush, MDLflushBuffer, hdr.hdr);
	}
    }

void mdlClass::bookkeeping() {
#ifdef USE_CUDA
    if (cudaCtx) CUDA_checkForRecovery(cudaCtx);
#endif
    combine_all_incoming();
    }

/* Do one piece of work. Return 0 if there was no work. */
int mdlClass::DoSomeWork() {
    MDLwqNode *work;
    int rc = 0;
    CacheCheck();
    if (!OPA_Queue_is_empty(&wq)) {
	/* Grab a work package, and perform it */
	OPA_Queue_dequeue(&wq, work, MDLwqNode, q.hdr);
	OPA_decr_int(&wqCurSize);
	while ( (*work->doFcn)(work->ctx) != 0 ) {
	    CacheCheck();
	    }
	rc = 1;
	/* Send it back to the original thread for completion (could be us) */
	if (work->iCoreOwner == Core()) goto reQueue;
	OPA_Queue_enqueue(&pmdl[work->iCoreOwner]->wqDone, work, MDLwqNode, q.hdr);
	}
    while (!OPA_Queue_is_empty(&wqDone)) {
	OPA_Queue_dequeue(&wqDone, work, MDLwqNode, q.hdr);
    reQueue:
	(*work->doneFcn)(work->ctx);
	work->ctx = NULL;
	work->doFcn = NULL;
	work->doneFcn = NULL;
	OPA_Queue_enqueue(&wqFree, work, MDLwqNode, q.hdr);
	}
    return rc;
    }

extern "C" void mdlCompleteAllWork(MDL mdl) { reinterpret_cast<mdlClass *>(mdl)->CompleteAllWork(); }
void mdlClass::CompleteAllWork() {
    while(DoSomeWork()) {}
#ifdef USE_CL
    while(CL_flushDone(clCtx)) {}
#endif
#ifdef USE_CUDA
    if (cudaCtx) while(CUDA_flushDone(cudaCtx)) {}
#endif
    }

void /*MDLserviceElement*/ * mdlClass::WaitThreadQueue(int iQueue) {
    MDLserviceElement *qhdr;
    while (OPA_Queue_is_empty(inQueue+iQueue)) {
	checkMPI(); // Only does something on the MPI thread
	bookkeeping();
	if (DoSomeWork() == 0) {
#ifdef _MSC_VER
	    SwitchToThread();
#else
	    sched_yield();
#endif
	    }
	}
    OPA_Queue_dequeue(inQueue+iQueue, qhdr, MDLserviceElement, hdr);
    return qhdr;
    }

/* Synchronize threads */
extern "C" void mdlThreadBarrier(MDL mdl) { reinterpret_cast<mdlClass *>(mdl)->ThreadBarrier(); }
void mdlClass::ThreadBarrier() {
    MDLserviceElement *qhdr;
    int i;

    if (Core()) {
	SendThreadMessage(MDL_TAG_BARRIER,0,&inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = CAST(MDLserviceElement *,WaitThreadQueue(MDL_TAG_BARRIER));
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	for(i=1; i<Cores(); ++i) {
	    qhdr = CAST(MDLserviceElement *,WaitThreadQueue(MDL_TAG_BARRIER));
	    assert(qhdr->iServiceID == MDL_SE_BARRIER_REQUEST);
	    }
	for(i=1; i<Cores(); ++i) {
	    SendThreadMessage(MDL_TAG_BARRIER,i,&pmdl[i]->inMessage,MDL_SE_BARRIER_REPLY);
	    }
	}
    }

int mdlClass::mdl_start_MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, int stype) {
    int iCore = dest - pmdl[0]->Self();
    int bOnNode = (iCore >= 0 && iCore < Cores());
    MDLserviceSend *send = &sendRequest;

    assert(dest != Self());
    send->buf = buf;
    send->count = count;
    send->datatype = datatype;
    send->target = dest;
    send->tag = tag;

    /* To send on node, we give the other thread the message. It will consume it after posting a receive */
    if (bOnNode && stype!=MDL_SE_SEND_REQUEST && stype!=MDL_SE_SEND_REPLY ) {
	SendThreadMessage(MDL_TAG_MAX + Core(),iCore,send,stype); /* From me */
	}
    /* Pass this off to the MPI thread (which could be ourselves) */
    else {
	SendToMPI(send,stype);
	}
    return MPI_SUCCESS;
    }

int mdlClass::mdl_MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag) {
    int rc = mdl_start_MPI_Ssend(buf,count,datatype,dest,tag,MDL_SE_MPI_SSEND);
    WaitThreadQueue(0); /* Wait for Send to complete */
    return rc;
    }

int mdlClass::mdl_MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag) {
    int rc = mdl_start_MPI_Ssend(buf,count,datatype,dest,tag,MDL_SE_MPI_SEND);
    WaitThreadQueue(0); /* Wait for Send to complete */
    return rc;
    }

int mdlClass::mdl_remote_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, int *nBytes, int iServiceID) {
    MDLserviceSend *send = &recvRequest;
    send->buf = buf;
    send->count = count;
    send->datatype = datatype;
    send->target = source;
    send->tag = tag;
    SendToMPI(send,iServiceID);
    send = CAST(MDLserviceSend *,WaitThreadQueue(tag)); /* Wait for the "send" to come back to us. */
    *nBytes = send->count;
    return MPI_SUCCESS;
    }

/*
**
*/
int mdlClass::mdl_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, int *nBytes) {
    int iCore = source - pmdl[0]->Self();
    int bOnNode = (iCore >= 0 && iCore < Cores());
    MDLserviceSend *send = &recvRequest;

    assert(source != Self());

    /* If the sender is on-node, we just wait for the send to appear and make a copy */
    if (bOnNode) {
	MDLserviceElement *qhdr;
	qhdr = CAST(MDLserviceElement *,WaitThreadQueue(MDL_TAG_MAX + iCore)); /* From me */
	assert(qhdr->iServiceID == MDL_SE_MPI_SSEND);
	assert(qhdr->iCoreFrom == iCore);
	send = (MDLserviceSend *)qhdr;
	assert(send->tag == tag);
	memcpy(buf,send->buf,send->count);
	*nBytes = send->count;
	SendThreadMessage(0,qhdr->iCoreFrom,qhdr,MDL_SE_MPI_SSEND);
	}

    /* Off-node: We ask the MPI thread to post a receive for us. */
    else {
	return mdl_remote_MPI_Recv(buf,count,datatype,source,tag,nBytes,MDL_SE_MPI_RECV);
	}
    return MPI_SUCCESS;
    }

/*
** Here we do a bidirectional message exchange between two threads.
*/
int mdlClass::mdl_MPI_Sendrecv(
    void *sendbuf, int sendcount, MPI_Datatype sendtype,
    int dest, int sendtag, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int source, int recvtag, int *nReceived) {
    mdl_start_MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag,MDL_SE_MPI_SSEND);
    mdl_MPI_Recv(recvbuf,recvcount,recvtype,source,recvtag,nReceived);
    WaitThreadQueue(0); /* Wait for Send to complete */
    return MPI_SUCCESS;
    }

int mdlClass::mdl_MPI_Barrier() {
    MDLserviceElement *qhdr;
    int i, id;
    int nBytes;

    CompleteAllWork();
    TimeAddComputing();
    wqAccepting = 1;
    if (Core()) {
	SendThreadMessage(MDL_TAG_BARRIER,0,&inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = CAST(MDLserviceElement *,WaitThreadQueue(MDL_TAG_BARRIER));
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	for(i=1; i<Cores(); ++i) {
	    qhdr = CAST(MDLserviceElement *,WaitThreadQueue(MDL_TAG_BARRIER));
	    assert(qhdr->iServiceID == MDL_SE_BARRIER_REQUEST);
	    }

	if (Self()==0) {
	    for(i=1; i<Procs(); ++i) {
		mdl_MPI_Recv(&id, sizeof(id), MPI_BYTE, MPI_ANY_SOURCE, MDL_TAG_BARRIER, &nBytes);
		}
	    for(i=1; i<Procs(); ++i) {
		id = i;
		mdl_MPI_Send(&id, sizeof(id),
		    MPI_BYTE, ProcToThread(i), MDL_TAG_BARRIER);
		}
	    }
	else {
	    id = Proc();
	    mdl_MPI_Ssend(&id,sizeof(id),MPI_BYTE,0,MDL_TAG_BARRIER);
	    mdl_MPI_Recv(&id, sizeof(id),MPI_BYTE,0,MDL_TAG_BARRIER,&nBytes);
	    }
#if 0 /* Hard barrier is not required here */
	SendToMPI(&inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = WaitThreadQueue(MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
#endif
	for(i=1; i<Cores(); ++i) {
	    SendThreadMessage(MDL_TAG_BARRIER,i,&pmdl[i]->inMessage,MDL_SE_BARRIER_REPLY);
	    }
	}
    wqAccepting = 0;
    TimeAddSynchronizing();
    return MPI_SUCCESS;
    }

void mdlClass::init(bool bDiag) {
    mpi = NULL;

    /*
    ** Our work queue. We can defer work for when we are waiting, or get work
    ** from other threads if we are idle.
    */
    OPA_Queue_init(&wq);
    OPA_Queue_init(&wqDone);
    OPA_Queue_init(&wqFree);
    wqMaxSize = 0;
    wqAccepting = 0;
    wqLastHelper = 0;
    OPA_store_int(&wqCurSize,0);

    OPA_Queue_header_init(&inMessage.hdr);
    OPA_Queue_header_init(&sendRequest.svc.hdr);
    OPA_Queue_header_init(&recvRequest.svc.hdr);

    /*
    ** Set default "maximums" for structures. These are NOT hard
    ** maximums, as the structures will be realloc'd when these
    ** values are exceeded.
    */
    nMaxSrvBytes = -1;
    /*
    ** Allocate service buffers.
    */
    pszIn = NULL;
    pszOut = NULL;
    pszBuf = NULL;
    /*
     ** Allocate swapping transfer buffer. This buffer remains fixed.
     */
    pszTrans = CAST(char *,malloc(MDL_TRANS_SIZE));
    assert(pszTrans != NULL);
    /*
    ** Allocate initial cache spaces.
    */
    nMaxCacheIds = MDL_DEFAULT_CACHEIDS;
    cache = CAST(CACHE *,malloc(nMaxCacheIds*sizeof(CACHE)));
    assert(cache != NULL);
    OPA_Queue_init(&wqCacheFlush);
    nFlushOutBytes = 0;

    OPA_Queue_init(&coreFlushBuffers);
    for(int i=0; i<9; ++i) {
	MDLflushBuffer *pFlush = CAST(MDLflushBuffer *,malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE));
	assert(pFlush!=NULL);
	pFlush->nBufferSize = MDL_FLUSH_DATA_SIZE;
	pFlush->nBytes = 0;
	OPA_Queue_enqueue(&coreFlushBuffers, pFlush, MDLflushBuffer, hdr.hdr);
	}
    coreFlushBuffer = CAST(MDLflushBuffer *,malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE));
    assert(coreFlushBuffer!=NULL);
    coreFlushBuffer->nBufferSize = MDL_FLUSH_DATA_SIZE;
    coreFlushBuffer->nBytes = 0;

    /*
    ** Initialize caching spaces.
    */
    cacheSize = MDL_CACHE_SIZE;
    for (int i = 0; i<nMaxCacheIds; ++i) {
	cache[i].iType = MDL_NOCACHE;
	cache[i].iCID = i;
	cache[i].arc = NULL;
	}

    /* We need a queue for each TAG, and a receive queue from each thread. */
    inQueue = CAST(OPA_Queue_info_t *,malloc((MDL_TAG_MAX+Cores()) * sizeof(*inQueue)));
    for(int i=0; i<(MDL_TAG_MAX+Cores()); ++i) OPA_Queue_init(inQueue+i);

#ifdef USE_CUDA
    inCudaBufSize = outCudaBufSize = 0;
    cudaCtx = CUDA_initialize(Cores(), Core(), &mpi->queueWORK, &mpi->queueREGISTER);
#endif
#ifdef USE_CL
    clCtx = CL_initialize(clContext,Cores(),Core());
#endif

    bDiag = bDiag;
    if (bDiag) {
	char achDiag[256], ach[256];
	const char *tmp = strrchr(argv[0], '/');
	if (!tmp) tmp = argv[0];
	else ++tmp;
	sprintf(achDiag, "%s/%s.%d", ach, tmp, Self());
	fpDiag = fopen(achDiag, "w");
	assert(fpDiag != NULL);
	}
    }

mdlClass::mdlClass(mdlClass *mdl, int iMDL)
	: mdlBASE(mdl->argc,mdl->argv)  {
    iCore = iMDL;
    idSelf = mdl->Self() + iMDL;
    nThreads = mdl->Threads();
    nCores = mdl->Cores();
    nProcs = mdl->nProcs;
    iProc = mdl->iProc;
    iProcToThread = mdl->iProcToThread;

    iCoreMPI = mdl->iCoreMPI;
    pmdl = mdl->pmdl;

    fcnWorkerInit = mdl->fcnWorkerInit;
    fcnWorkerDone = mdl->fcnWorkerDone;
    fcnMaster = mdl->fcnMaster;
    init();
    mpi = mdl->mpi;
    }

mdlClass::mdlClass(void (*fcnMaster)(MDL,void *),void * (*fcnWorkerInit)(MDL),void (*fcnWorkerDone)(MDL,void *),int argc, char **argv)
	: mdlBASE(argc,argv) {
    this->fcnWorkerInit = fcnWorkerInit;
    this->fcnWorkerDone = fcnWorkerDone;
    this->fcnMaster = fcnMaster;
    init();
    }

void mdlClass::drainMPI() {
    while(checkMPI()) {
#ifdef _MSC_VER
	SwitchToThread();
#else
	sched_yield();
#endif
	}
    }

void mdlClass::run_master() {
    (*fcnMaster)(reinterpret_cast<MDL>(this),worker_ctx);
    int id;
    for (id=1;id<Threads();++id) {
	int rID = ReqService(id,SRV_STOP,NULL,0);
	GetReply(rID,NULL,NULL);
	}
    }

void *mdlClass::WorkerThread() {
    void *result;
#ifdef USE_ITT
    char szName[20];
    sprintf(szName,"ID %d", Core());
    __itt_thread_set_name(szName);
#endif
    worker_ctx = (*fcnWorkerInit)(reinterpret_cast<MDL>(this));
    CommitServices();
    if (Self()) Handler();
    else run_master();
    (*fcnWorkerDone)(reinterpret_cast<MDL>(this),worker_ctx);

    if (Core() != iCoreMPI) {
	SendToMPI(&inMessage,MDL_SE_STOP);
	WaitThreadQueue(0); /* Wait for Send to complete */
	}
    else {
	drainMPI();
	}
    return result;
    }
void *mdlClass::mdlWorkerThread(void *vmdl) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(vmdl);
    return mdl->WorkerThread();
    }

void mdlClass::cleanupMDL() {
    int i;
    for (i = 0; i<nMaxCacheIds; ++i) {
	if (cache[i].arc) delete cache[i].arc;
	}
    free(pszIn);
    free(pszOut);
    free(pszBuf);
    free(inQueue);
    free(pszTrans);
    free(cache);
    }

/*
** This function will transfer a block of data using a pack function.
** The corresponding node must call mdlRecv.
*/

#define SEND_BUFFER_SIZE (1*1024*1024)

extern "C" void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx) { reinterpret_cast<mdlClass *>(mdl)->Send(id,pack,ctx); }
void mdlClass::Send(int id,mdlPack pack, void *ctx) {
    size_t nBuff;
    char *vOut;

    vOut = CAST(char *,malloc(SEND_BUFFER_SIZE));
    mdlassert(this,vOut!=NULL);

    do {
	nBuff = (*pack)(ctx,&id,SEND_BUFFER_SIZE,vOut);
	mdl_MPI_Ssend(vOut,nBuff,MPI_BYTE,id,MDL_TAG_SEND);
	}
    while ( nBuff != 0 );

    free(vOut);
    }

extern "C" void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx) { reinterpret_cast<mdlClass *>(mdl)->Recv(id,unpack,ctx); }
void mdlClass::Recv(int id,mdlPack unpack, void *ctx) {
    void *vIn;
    size_t nUnpack;
    int nBytes;
    int inid;

    if ( id < 0 ) id = MPI_ANY_SOURCE;

    vIn = malloc(SEND_BUFFER_SIZE);
    mdlassert(this,vIn!=NULL);

    do {
	mdl_MPI_Recv(vIn,SEND_BUFFER_SIZE,MPI_BYTE,id,MDL_TAG_SEND,&nBytes);
	inid = id; //status.MPI_SOURCE;
	nUnpack = (*unpack)(ctx,&inid,nBytes,vIn);
	}
    while (nUnpack>0 && nBytes>0);

    free(vIn);
    }

/*
 ** This is a tricky function. It initiates a bilateral transfer between
 ** two threads. Both threads MUST be expecting this transfer. The transfer
 ** occurs between idSelf <---> 'id' or 'id' <---> idSelf as seen from the
 ** opposing thread. It is designed as a high performance non-local memory
 ** swapping primitive and implementation will vary in non-trivial ways
 ** between differing architectures and parallel paradigms (eg. message
 ** passing and shared address space). A buffer is specified by 'pszBuf'
 ** which is 'nBufBytes' in size. Of this buffer the LAST 'nOutBytes' are
 ** transfered to the opponent, in turn, the opponent thread transfers his
 ** nBufBytes to this thread's buffer starting at 'pszBuf'.
 ** If the transfer completes with no problems the function returns 1.
 ** If the function returns 0 then one of the players has not received all
 ** of the others memory, however he will have successfully transfered all
 ** of his memory.
 */
extern "C" int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes, size_t *pnSndBytes,size_t *pnRcvBytes) {
    return reinterpret_cast<mdlClass *>(mdl)->Swap(id,nBufBytes,vBuf,nOutBytes,pnSndBytes,pnRcvBytes);
    }
int mdlClass::Swap(int id,size_t nBufBytes,void *vBuf,size_t nOutBytes, size_t *pnSndBytes,size_t *pnRcvBytes) {
    size_t nInBytes,nOutBufBytes;
    int nInMax,nOutMax,nBytes;
    char *pszBuf = CAST(char *,vBuf);
    char *pszIn,*pszOut;
    struct swapInit {
	size_t nOutBytes;
	size_t nBufBytes;
	} swi,swo;

    *pnRcvBytes = 0;
    *pnSndBytes = 0;
    /*
     **	Send number of rejects to target thread amount of free space
     */
    swi.nOutBytes = nOutBytes;
    swi.nBufBytes = nBufBytes;
    mdl_MPI_Sendrecv(&swi, sizeof(swi), MPI_BYTE, id, MDL_TAG_SWAPINIT,
	&swo, sizeof(swo), MPI_BYTE, id, MDL_TAG_SWAPINIT, &nBytes);
    assert(nBytes == sizeof(swo));
    nInBytes = swo.nOutBytes;
    nOutBufBytes = swo.nBufBytes;
    /*
     ** Start bilateral transfers. Note: One processor is GUARANTEED to
     ** complete all its transfers.
     */
    assert(nBufBytes >= nOutBytes);
    pszOut = &pszBuf[nBufBytes-nOutBytes];
    pszIn = pszBuf;
    while (nOutBytes && nInBytes) {
	/*
	 ** nOutMax is the maximum number of bytes allowed to be sent
	 ** nInMax is the number of bytes which will be received.
	 */
	nOutMax = size_t_to_int((nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE);
	nOutMax = size_t_to_int((nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes);
	nInMax = size_t_to_int((nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE);
	nInMax = size_t_to_int((nInMax < nBufBytes)?nInMax:nBufBytes);
	/*
	 ** Copy to a temp buffer to be safe.
	 */
	memcpy(pszTrans,pszOut,nOutMax);
	mdl_MPI_Sendrecv(pszTrans,nOutMax, MPI_BYTE, id, MDL_TAG_SWAP,
	    pszIn,nInMax, MPI_BYTE, id, MDL_TAG_SWAP, &nBytes);
	assert(nBytes == nInMax);
	/*
	 ** Adjust pointers and counts for next itteration.
	 */
	pszOut = &pszOut[nOutMax];
	nOutBytes -= nOutMax;
	nOutBufBytes -= nOutMax;
	*pnSndBytes += nOutMax;
	pszIn = &pszIn[nInMax];
	nInBytes -= nInMax;
	nBufBytes -= nInMax;
	*pnRcvBytes += nInMax;
	}

    /*
     ** At this stage we perform only unilateral transfers, and here we
     ** could exceed the opponent's storage capacity.
     ** Note: use of Ssend is mandatory here, also because of this we
     ** don't need to use the intermediate buffer mdl->pszTrans.
     */
    while (nOutBytes && nOutBufBytes) {
	nOutMax = size_t_to_int((nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE);
	nOutMax = size_t_to_int((nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes);
	mdl_MPI_Ssend(pszOut,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP);
	pszOut = &pszOut[nOutMax];
	nOutBytes -= nOutMax;
	nOutBufBytes -= nOutMax;
	*pnSndBytes += nOutMax;
	}
    while (nInBytes && nBufBytes) {
	nInMax = size_t_to_int((nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE);
	nInMax = size_t_to_int((nInMax < nBufBytes)?nInMax:nBufBytes);
	mdl_MPI_Recv(pszIn,nInMax,MPI_BYTE,id,MDL_TAG_SWAP,&nBytes);
	assert(nBytes == nInMax);
	pszIn = &pszIn[nInMax];
	nInBytes -= nInMax;
	nBufBytes -= nInMax;
	*pnRcvBytes += nInMax;
	}
    if (nOutBytes) return(0);
    else if (nInBytes) return(0);
    else return(1);
    }

void mdlClass::CommitServices() {
    int nMaxBytes;
    nMaxBytes = (nMaxInBytes > nMaxOutBytes) ? nMaxInBytes : nMaxOutBytes;
    if (nMaxBytes > nMaxSrvBytes) {
        pszIn = CAST(char *,realloc(pszIn, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement)));
        assert(pszIn != NULL);
	pszOut = CAST(char *,realloc(pszOut, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement)));
        assert(pszOut != NULL);
	pszBuf = CAST(char *,realloc(pszBuf, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement)));
        assert(pszBuf != NULL);
        nMaxSrvBytes = nMaxBytes;
        }
    /* We need a thread barrier here because we share these buffers */
    ThreadBarrier();
    }

void mdlAddService(MDL cmdl,int sid,void *p1,
		   fcnService_t *fcnService,
		   int nInBytes,int nOutBytes) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    mdl->AddService(sid, p1, fcnService, nInBytes, nOutBytes);
    }

extern "C" int mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes) { return reinterpret_cast<mdlClass *>(mdl)->ReqService(id,sid,vin,nInBytes); }
int mdlClass::ReqService(int id,int sid,void *vin,int nInBytes) {
    char *pszIn = CAST(char *,vin);
    SRVHEAD *ph = (SRVHEAD *)(pszBuf);
    char *pszOut = CAST(char *,ph + 1);
    ph->idFrom = Self();
    ph->sid = sid;
    if (!pszIn) ph->nInBytes = 0;
    else ph->nInBytes = nInBytes;
    ph->nOutBytes = 0;
    if (nInBytes > 0) {
	assert(pszIn != NULL);
	memcpy(pszOut, pszIn, nInBytes);
	}
    mdl_start_MPI_Ssend(ph, nInBytes + (int)sizeof(SRVHEAD), MPI_BYTE, id, MDL_TAG_REQ, MDL_SE_SEND_REQUEST);
    WaitThreadQueue(0); /* Wait for Send to complete */
    return ph->replyTag;
    }

extern "C" void mdlGetReply(MDL mdl,int rID,void *vout,int *pnOutBytes) { reinterpret_cast<mdlClass *>(mdl)->GetReply(rID,vout,pnOutBytes); }
void mdlClass::GetReply(int rID,void *vout,int *pnOutBytes) {
    char *pszOut = CAST(char *,vout);
    SRVHEAD *ph = (SRVHEAD *)pszBuf;
    char *pszIn = &pszBuf[sizeof(SRVHEAD)];
    int nBytes,id;
    id = rID;
    mdl_remote_MPI_Recv(pszBuf, nMaxSrvBytes + (int)sizeof(SRVHEAD), MPI_BYTE,
	id, MDL_TAG_RPL, &nBytes, MDL_SE_RECV_REPLY);
    assert(nBytes == ph->nOutBytes + sizeof(SRVHEAD));
    if (ph->nOutBytes > 0 && pszOut != NULL)
	memcpy(pszOut, pszIn, ph->nOutBytes);
    if (pnOutBytes) *pnOutBytes = ph->nOutBytes;
    }

void mdlClass::Handler() {
    MDLserviceElement *qhi = (MDLserviceElement *)(pszIn);
    MDLserviceElement *qho = (MDLserviceElement *)(pszOut);
    SRVHEAD *phi = (SRVHEAD *)(qhi + 1);
    SRVHEAD *pho = (SRVHEAD *)(qho + 1);
    char *pszIn = (char *)(phi + 1);
    char *pszOut = (char *)(pho + 1);
    int sid,id,tag,nOutBytes,nBytes;

    do {
	/* We ALWAYS use MPI to send requests. */
	mdl_remote_MPI_Recv(phi, nMaxSrvBytes + sizeof(SRVHEAD),
	    MPI_BYTE, MPI_ANY_SOURCE, MDL_TAG_REQ, &nBytes, MDL_SE_RECV_REQUEST);
	assert(nBytes == phi->nInBytes + sizeof(SRVHEAD));
	id = phi->idFrom;
	sid = phi->sid;
	assert(sid < nMaxServices);
	if (phi->nInBytes > psrv[sid].nInBytes) {
	    printf("ERROR: pid=%d, sid=%d, nInBytes=%d, sid.nInBytes=%d\n",
		Self(), sid, phi->nInBytes, psrv[sid].nInBytes);
	    }
	assert(phi->nInBytes <= psrv[sid].nInBytes);
	nOutBytes = 0;
	assert(psrv[sid].fcnService != NULL);

	nOutBytes = (*psrv[sid].fcnService)(psrv[sid].p1, pszIn, phi->nInBytes,
	    pszOut, psrv[sid].nOutBytes);
	if (nOutBytes > psrv[sid].nOutBytes) {
	    fprintf(stderr,"%d > %d: sid=%d\n",
		nOutBytes, psrv[sid].nOutBytes, sid);
	    assert(nOutBytes <= psrv[sid].nOutBytes);
	    }
	pho->idFrom = Self();
	pho->replyTag = phi->replyTag;
	pho->sid = sid;
	pho->nInBytes = phi->nInBytes;
	pho->nOutBytes = nOutBytes;
	tag = phi->replyTag;
	mdl_start_MPI_Ssend(pho, nOutBytes + sizeof(SRVHEAD), MPI_BYTE, id, tag, MDL_SE_SEND_REPLY);
	WaitThreadQueue(0); /* Wait for Send to complete */
	} while (sid != SRV_STOP);
    }

/*
 ** Special MDL memory allocation functions for allocating memory
 ** which must be visible to other processors thru the MDL cache
 ** functions.
 ** mdlMalloc() is defined to return a pointer to AT LEAST iSize bytes
 ** of memory. This pointer will be passed to either mdlROcache or
 ** mdlCOcache as the pData parameter.
 ** For PVM and most machines these functions are trivial, but on the
 ** T3D and perhaps some future machines these functions are required.
 */
void *mdlMalloc(MDL mdl,size_t iSize) {
/*    void *ptr;
    int rc = MPI_Alloc_mem(iSize,MPI_INFO_NULL,&ptr);
    if (rc) return NULL;
    return ptr;*/
    return(malloc(iSize));
    }

void mdlFree(MDL mdl,void *p) {
/*    MPI_Free_mem(p);*/
    free(p);
    }

/* This is a "thread collective" call. */
void *mdlSetArray(MDL cmdl,size_t nmemb,size_t size,void *vdata) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    char *data = CAST(char *,vdata);
    mdl->nMessageData = nmemb * size;
    mdl->ThreadBarrier();
    if (mdlCore(mdl)==0) {
	int i;
	for(i=0; i<mdlCores(mdl); ++i) {
	    mdl->pmdl[i]->pvMessageData = data;
	    data += mdl->pmdl[i]->nMessageData;
	    }
	}
    mdl->ThreadBarrier();
    return mdl->pvMessageData;
    }

/*
** New collective form. All cores must call this at the same time.
** Normally "size" is identical on all cores, while nmemb may differ
** but this is not strictly a requirement.
*/
void *mdlMallocArray(MDL cmdl,size_t nmemb,size_t size,size_t minSize) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    char *data;
    size_t iSize;
    mdl->nMessageData = nmemb * size;
    mdl->ThreadBarrier();
    if (mdlCore(mdl)==0) {
	int i;
	iSize = 0;
	for(i=0; i<mdlCores(mdl); ++i) iSize += mdl->pmdl[i]->nMessageData;
	if (iSize < minSize) iSize = minSize;
	data = CAST(char *,mdlMalloc(cmdl,iSize));
	for(i=0; i<mdlCores(mdl); ++i) {
	    mdl->pmdl[i]->pvMessageData = data;
	    data += mdl->pmdl[i]->nMessageData;
	    }
	}
    mdl->ThreadBarrier();
    data = CAST(char *,mdl->pvMessageData);
    iSize = nmemb * size;
    if (iSize > 4096) iSize -= 4096; /* Cheesy page hack */
    memset(data,0,iSize); /* First touch */
    mdl->ThreadBarrier();
    return data;
    }
void mdlFreeArray(MDL mdl,void *p) {
    if (mdlCore(mdl)==0) free(p);
    }

/*
** This is the default element fetch routine.  It impliments the old behaviour
** of a single large array.  New data structures need to be more clever.
*/
static void *getArrayElement(void *vData,int i,int iDataSize) {
    char *pData = CAST(char *,vData);
    return pData + (size_t)i*(size_t)iDataSize;
    }

void mdlSetCacheSize(MDL cmdl,int cacheSize) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    mdl->cacheSize = cacheSize;
    }

extern "C" void mdlCacheCheck(MDL mdl) { reinterpret_cast<mdlClass *>(mdl)->CacheCheck(); }
void mdlClass::CacheCheck() {
    checkMPI(); // Only does something on the MPI thread
    bookkeeping();
    }

int mdlCacheStatus(MDL cmdl,int cid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    assert(cid >= 0);
    if (cid >= mdl->nMaxCacheIds) return MDL_NOCACHE;
    CACHE *c = &mdl->cache[cid];
    return c->iType;
    }

/*
 ** Common initialization for all types of caches.
 */
CACHE *mdlClass::CacheInitialize(
    int cid,
    void * (*getElt)(void *pData,int i,int iDataSize),
    void *pData,int iDataSize,int nData,
    void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *)) {

    CACHE *c;
    int i;
    cacheOpenClose coc;

    /*
     ** We cannot reallocate this structure because there may be other threads accessing it.
     */
    assert(cid >= 0 && cid <nMaxCacheIds);
    if (cid<0 || cid >= nMaxCacheIds) abort();
    c = &cache[cid];
    assert(c->iType == MDL_NOCACHE);
    c->getElt = getElt==NULL ? getArrayElement : getElt;
    c->pData = pData;
    c->iDataSize = iDataSize;
    c->nData = nData;
    if (c->iDataSize > MDL_CACHE_DATA_SIZE) c->nLineBits = 0;
    else c->nLineBits = log2(MDL_CACHE_DATA_SIZE / c->iDataSize);
    /*
    ** Let's not have too much of a good thing. We want to fetch in cache lines both for
    ** performance, and to save memory (ARC cache has a high overhead), but too much has
    ** a negative effect on performance. Empirically, 16 elements maximal is optimal.
    ** This was tested with the particle, cell and P(k) caches.
    */
    if (c->nLineBits > 4) c->nLineBits = 4;
    c->nLineElements = 1 << c->nLineBits;
    c->nLineMask = (1 << c->nLineBits) - 1;
    c->iLineSize = c->nLineElements*c->iDataSize;

    c->nAccess = 0;
    c->nMiss = 0;				/* !!!, not NB */
    c->nColl = 0;				/* !!!, not NB */
    c->arc = arcReinitialize(c->arc,cacheSize/c->iLineSize,c->iLineSize,c);
    c->pOneLine = CAST(char *,malloc(c->iLineSize));

    /*
     ** Set up the request message as much as possible!
     */
    OPA_Queue_header_init(&c->cacheRequest.svc.hdr);
    c->cacheRequest.svc.iServiceID = MDL_SE_CACHE_REQUEST;
    c->cacheRequest.svc.iCoreFrom = Core();
    c->cacheRequest.caReq.cid = cid;
    c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
    c->cacheRequest.caReq.idFrom = Self();
    c->cacheRequest.caReq.idTo = -1;
    c->cacheRequest.caReq.iLine = -1;
    c->cacheRequest.caReq.nItems = 1;
    c->cacheRequest.request = MPI_REQUEST_NULL;

    /* Read-only or combiner caches */
    assert( init==NULL && combine==NULL || init!=NULL && combine!=NULL );
    c->iType = (init==NULL ? MDL_ROCACHE : MDL_COCACHE);
    c->init = init;
    c->combine = combine;
    c->ctx = ctx;

    /* Nobody should start using this cache until all threads have started it! */
    mdl_MPI_Barrier();

    /* We might need to resize the cache buffer */
    coc.iDataSize = c->iLineSize;
    SendToMPI(&coc,MDL_SE_CACHE_OPEN);
    WaitThreadQueue(0);

    return(c);
    }

/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData) {
    reinterpret_cast<mdlClass *>(mdl)->CacheInitialize(cid,getElt,pData,iDataSize,nData,NULL,NULL,NULL);
    }

/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData,
		void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *)) {
    assert(init);
    assert(combine);
    reinterpret_cast<mdlClass *>(mdl)->CacheInitialize(cid,getElt,pData,iDataSize,nData,ctx,init,combine);
    }

extern "C" void mdlCacheBarrier(MDL mdl,int cid) { reinterpret_cast<mdlClass *>(mdl)->CacheBarrier(cid); }
void mdlClass::CacheBarrier(int cid) {
    mdl_MPI_Barrier();
    }

MDLflushBuffer *mdlClass::flush_core_buffer() {
    MDLflushBuffer *pBuffer = coreFlushBuffer;
    if (pBuffer->nBytes) {
	SendToMPI(&pBuffer->hdr.svc,MDL_SE_CACHE_FLUSH);
	while (OPA_Queue_is_empty(&coreFlushBuffers)) {
	    bookkeeping();
	    checkMPI(); // Only does something on the MPI thread
	    }
	OPA_Queue_dequeue(&coreFlushBuffers, pBuffer, MDLflushBuffer, hdr.hdr);
	coreFlushBuffer = pBuffer;
	pBuffer->nBytes = 0;
	pBuffer->iRankTo = -1;
	}
    return pBuffer;
    }

/* This CDB must not be on any list when destage is called */
CDB *ARC::destage(CDB *temp) {
    assert(temp->data != NULL);
    assert(temp->data[-1] == _ARC_MAGIC_);
    if (temp->uId & _DIRTY_) {     /* if dirty, evict before free */
	CACHE *c = cache;
	int iSize = sizeof(CAHEAD) + c->iLineSize;
	MDLflushBuffer *pBuffer = mdl->coreFlushBuffer;
	if (pBuffer->nBytes + iSize > pBuffer->nBufferSize) {
	    pBuffer = mdl->flush_core_buffer();
	    }
	CAHEAD *ca = (CAHEAD *)( (char *)(pBuffer+1) + pBuffer->nBytes);
	char *pData = (char *)(ca+1);
	ca->cid = c->iCID;
	ca->mid = MDL_MID_CACHEFLUSH;
	ca->nItems = 1;
	ca->idFrom = mdlSelf(mdl);
	ca->idTo = temp->uId&_IDMASK_;
	ca->iLine = temp->uPage;
	memcpy(pData,temp->data,c->iLineSize);
	pBuffer->nBytes += iSize;
	temp->uId &= ~ _DIRTY_;    /* No longer dirty */
	}
    return temp;
    }

uint64_t *ARC::replace(bool bInB2) {
    CDB *temp;
    uint64_t *data;
    uint32_t max = (target_T1 > 1)?(target_T1+(bInB2?0:1)):1;
    if (T1Length >= max) { /* T1’s size exceeds target? */
                                        /* yes: T1 is too big */
	temp = T1->prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->prev;
	    if (temp == T1) {           /* all pages in T1 are currently locked, try T2 in this case */
		temp = T2->prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->prev;
		    if (temp == T2) return(NULL); /* all pages are currently locked, give up! */
		}
		goto replace_T2;
	    }
	}
	goto replace_T1;
    } else {
	/* no: T1 is not too big */
	temp = T2->prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->prev;
	    if (temp == T2) {           /* all pages in T2 are currently locked, try T1 in this case */
		temp = T1->prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->prev;
		    if (temp == T1) return(NULL); /* all pages are currently locked, give up! */
		}
		goto replace_T1;
	    }
	}
	goto replace_T2;
    }
    assert(0);
    if (0) {  /* using a Duff's device to handle the replacement */
    replace_T1: 
        temp->remove_from_list();          /* grab LRU unlocked page from T1 */
	destage(temp);     /* if dirty, evict before overwrite */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
	B1->mru_insert(temp);           /* put it on B1 */
	temp->uId = (temp->uId&_IDMASK_)|_B1_;  /* need to be careful here because it could have been _P1_ */
/* 	assert( (temp->uId&_WHERE_) == _B1_ ); */
        T1Length--; B1Length++;          /* bookkeep */
	/*assert(B1Length<=nCache);*/
    }
    if (0) {
    replace_T2:
        temp->remove_from_list();          /* grab LRU unlocked page from T2 */
	destage(temp);     /* if dirty, evict before overwrite */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
        B2->mru_insert(temp);           /* put it on B2 */
        temp->uId |= _B2_;          /* note that fact */
/* 	assert( (temp->uId&_WHERE_) == _B2_ ); */
        T2Length--; B2Length++;          /* bookkeep */
    }
    /*assert(data!=NULL);*/
    return data;
}

void ARC::RemoveAll() {
    CDB *temp;
    for(;T1Length;--T1Length) {
	temp = remove_from_hash(T1->lru_remove());
	destage(temp);     /* if dirty, evict before free */
	Free->lru_insert(temp);
	}
    for(;T2Length;--T2Length) {
	temp = remove_from_hash(T2->lru_remove());
	destage(temp);     /* if dirty, evict before free */
	Free->lru_insert(temp);
	}
    for(;B1Length;--B1Length) {
	temp = remove_from_hash(B1->lru_remove());
	assert(temp->data == NULL);
	Free->mru_insert(temp);
	}
    for(;B2Length;--B2Length) {
	temp = remove_from_hash(B2->lru_remove());
	assert(temp->data == NULL);
	Free->mru_insert(temp);
	}
    assert(T1->next == T1);
    assert(T1->prev == T1);
    assert(B1->next == B1);
    assert(B1->prev == B1);
    assert(T2->next == T2);
    assert(T2->prev == T2);
    assert(B2->next == B2);
    assert(B2->prev == B2);
    }

extern "C" void mdlFlushCache(MDL mdl,int cid) { reinterpret_cast<mdlClass *>(mdl)->FlushCache(cid); }
void mdlClass::FlushCache(int cid) {
    CACHE *c = &cache[cid];

    TimeAddComputing();
    wqAccepting = 1;
    c->arc->RemoveAll();
    flush_core_buffer();
    ThreadBarrier();
    if (Core()==0) { /* This flushes all buffered data, not just our thread */
	SendToMPI(&inMessage,MDL_SE_CACHE_FLUSH_OUT);
	WaitThreadQueue(0);
	}
    /* We must wait for all threads to finish with this cache */
    mdl_MPI_Barrier();
    if (Core()==0) { /* This flushes all buffered data, not just our thread */
	SendToMPI(&inMessage,MDL_SE_CACHE_FLUSH_LCL);
	WaitThreadQueue(0);
	}
    ThreadBarrier();
    wqAccepting = 0;
    TimeAddSynchronizing();
    }

extern "C" void mdlFinishCache(MDL mdl,int cid) { reinterpret_cast<mdlClass *>(mdl)->FinishCache(cid); }
void mdlClass::FinishCache(int cid) {
    CACHE *c = &cache[cid];
    cacheOpenClose coc;

    TimeAddComputing();
    wqAccepting = 1;
    FlushCache(cid);

    coc.iDataSize = c->iDataSize;
    SendToMPI(&coc,MDL_SE_CACHE_CLOSE);
    WaitThreadQueue(0);

    /*
     ** Free up storage and finish.
     */
    free(c->pOneLine);
#ifdef FINISH_ARC
    delete c->arc;
    c->arc = NULL;
#endif
    c->iType = MDL_NOCACHE;
    wqAccepting = 0;
    TimeAddSynchronizing();
    }

void mdlPrefetch(MDL mdl,int cid,int iIndex, int id) {
    }

/*
** The highest order bit (1<<31) of uId encodes the dirty bit.
** The next 3 highest order bits of uId ((1<<30), (1<<29) and (1<<28)) are reserved 
** for list location and should be zero! The maximum legal uId is then (1<<28)-1.
*/
static inline CDB *arcSetPrefetchDataByHash(MDL cmdl,ARC *arc,uint32_t uPage,uint32_t uId,void *data,uint32_t uHash) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    CDB *temp;
    uint32_t tuId = uId&_IDMASK_;
    uint32_t L1Length;
    bool inB2=false;

    assert(data);
    for( temp = arc->Hash[uHash]; temp; temp = temp->coll ) {
	if (temp->uPage == uPage && (temp->uId&_IDMASK_) == tuId) break;
	}
    if (temp != NULL) {                       /* found in cache directory? */
	switch (temp->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	case _T1_:
	case _T2_:
	    return temp;              /* do nothing */
	case _B1_:                            /* B1 "hit": don't favour anything */
	    arc->B1Length--;                                           /* bookkeep */
	    goto doBcase;
	case _B2_:                            /* B2 "hit": don't favour anything */
	    arc->B2Length--;                                           /* bookkeep*/
	    inB2=true;
	doBcase:
	    /* Better would be to put this back on P1, but L1 may be full. */
	    temp->remove_from_list();                                   /* take off whichever list */
	    temp->data = arc->replace(inB2);                                /* find a place to put new page */
	    temp->uPage = uPage;                          /* bookkeep */
	    if (inB2) {
		temp->uId = _T2_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 		assert( (temp->uId&_WHERE_) == _T2_ ); */
		arc->T2->mru_insert(temp);                                     /* not seen yet, put on T1 */
		arc->T2Length++;
		}
	    else {
		temp->uId = _P1_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 		assert( (temp->uId&_WHERE_) == _P1_ ); */
		arc->T1->mru_insert(temp);                                     /* not seen yet, put on T1 */
		arc->T1Length++;
		/*assert(arc->T1Length+arc->B1Length<=arc->nCache);*/
		}
	    break;
	}
    } else {                                                              /* page is not in cache directory */
	L1Length = arc->T1Length + arc->B1Length;
	/*assert(L1Length<=arc->nCache);*/
	if (L1Length == arc->nCache) {                                   /* B1 + T1 full? */
	    if (arc->T1Length < arc->nCache) {                                           /* Still room in T1? */
		temp = arc->remove_from_hash(arc->B1->lru_remove());        /* yes: take page off B1 */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = arc->replace();                                /* find new place to put page */
	    } else {                                                      /* no: B1 must be empty */
		temp = arc->remove_from_hash(arc->T1->lru_remove());         /* take page off T1 */
		arc->destage(temp);     /* if dirty, evict before overwrite */
		arc->T1Length--;                                               /* bookkeep that */
	    }
	} else {                                                          /* B1 + T1 have less than c pages */
	    uint32_t nCache = arc->T1Length + arc->T2Length + arc->B1Length + arc->B2Length;
	    /*assert(L1Length < arc->nCache);*/
	    if (nCache >= arc->nCache) {         /* cache full? */
		/* Yes, cache full: */
		if (arc->T1Length + arc->T2Length + arc->B1Length + arc->B2Length == 2*arc->nCache) {
		    /* directory is full: */
		    /*assert(arc->B2Length>0);*/
		    temp = arc->remove_from_hash(arc->B2->lru_remove());/* here we lose memory of what was in lru B2 */
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=true;
		} else {                                                   /* cache directory not full, easy case */
		    temp = arc->Free->lru_remove();
		}
		temp->data = arc->replace(inB2);                                /* new place for page */
	    } else {                                                      /* cache not full, easy case */
		temp = arc->Free->lru_remove();
		assert(temp->data != NULL);               /* This CDB should have an unused page associated with it */
		temp->data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
	    }
	}
	arc->T1->mru_insert(temp);                                             /* not been seen yet, but put on T1 */
	arc->T1Length++;                                                       /* bookkeep: */
	/*assert(arc->T1Length+arc->B1Length<=arc->nCache);*/
	temp->uId = _P1_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 	assert( (temp->uId&_WHERE_) == _P1_ ); */
	temp->uPage = uPage;
	temp->coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */
    }
    memcpy(temp->data,data,arc->cache->iLineSize); /* Copy actual cache data amount */
    return temp;
}

static inline CDB *arcSetPrefetchData(MDL mdl,ARC *arc,uint32_t uPage,uint32_t uId,void *data) {
    return arcSetPrefetchDataByHash(mdl,arc,uPage,uId,data,MurmurHash2(uPage,uId&_IDMASK_)&arc->uHashMask);
    }

/*
** This releases a lock on a page by decrementing the lock counter for that page. It also checks
** that the magic number matches, and that the lock count is greater than 0. This is to prevent
** strange errors if the user tries to release the wrong pointer (could silently modify memory 
** that it should not if this check was not there).
*/
static inline void arcRelease(ARC *arc,uint64_t *p) {
    if (p>arc->dataBase && p<arc->dataLast) { /* Might have been a fast, read-only grab */
	/* We will be given an element, but this needs to be turned into a cache line */
	p = arc->dataBase + (p - arc->dataBase) / (arc->uDataSize+1) * (arc->uDataSize+1) + 1;
	uint64_t t = p[-1]-1;
	assert((t^_ARC_MAGIC_) < 0x00000000ffffffff);
	p[-1] = t;
	}
}

void mdlClass::queueCacheRequest(int cid, int iLine, int id) {
    int iCore = id - pmdl[0]->Self();
    /* Local requests for combiner cache are handled in finishCacheRequest() */
    if (iCore < 0 || iCore >= Cores() ) {
	CACHE *c = &cache[cid];
	c->cacheRequest.caReq.cid = cid;
	c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
	c->cacheRequest.caReq.idFrom = Self();
	c->cacheRequest.caReq.idTo = id;
	c->cacheRequest.caReq.iLine = iLine;
	c->cacheRequest.caReq.nItems = c->nLineElements;
	c->cacheRequest.pLine = c->pOneLine;
	SendToMPI(&c->cacheRequest,MDL_SE_CACHE_REQUEST);
	}
    }

void mdlClass::finishCacheRequest(int cid, int iLine, int id, CDB *temp,int bVirtual) {
    int iCore = id - pmdl[0]->Self();
    CACHE *c = &cache[cid];
    int i,s,n;
    char *pData = (char *)temp->data;
    uint32_t iIndex = iLine << c->nLineBits;

    s = iIndex;
    n = s + c->nLineElements;

    if (bVirtual) {
	memset(temp->data,0,c->iLineSize);
	}
    /* Local requests must be from a combiner cache if we get here */
    else if (iCore >= 0 && iCore < Cores() ) {
	mdlClass * omdl = pmdl[iCore];
	CACHE *oc = &omdl->cache[cid];
	if ( n > oc->nData ) n = oc->nData;
	for(i=s; i<n; i++ ) {
	    memcpy(pData,(*oc->getElt)(oc->pData,i,oc->iDataSize),oc->iDataSize);
	    pData += oc->iDataSize;
	    }
	}
    else {
	assert(iLine == temp->uPage);
	WaitThreadQueue(MDL_TAG_CACHECOM);
	assert(iLine == c->cacheRequest.caReq.iLine );
	memcpy(temp->data,c->pOneLine, c->iLineSize);
	}

    if (c->init) {
	pData = (char *)temp->data;
	for(i=s; i<n; i++ ) {
	    (*c->init)(c->ctx,pData);
	    pData += c->iDataSize;
	    }
	}

    }

/*
** The highest order bit (1<<31) of uId encodes the dirty bit.
** The next 3 highest order bits of uId ((1<<30), (1<<29) and (1<<28)) are reserved 
** for list location and should be zero! The maximum legal uId is then (1<<28)-1.
*/
void *mdlClass::Access(int cid, uint32_t uIndex, int uId, int bLock,int bModify,int bVirtual) {
    CACHE *c = &cache[cid];
    ARC *arc = c->arc;
    CDB *temp;
    uint32_t L1Length;
    uint32_t uHash,rat;
    uint32_t tuId = uId&_IDMASK_;
    uint32_t uLine;
    uint32_t iInLine;
    bool inB2=false;

    if (!(++c->nAccess & MDL_CHECK_MASK)) CacheCheck();

    /* Short circuit the cache if this belongs to another thread (or ourselves) */
    uint32_t uCore = uId - pmdl[0]->Self();
    if (uCore < Cores() && c->iType == MDL_ROCACHE ) {
	mdlClass * omdl = pmdl[uCore];
	c = &omdl->cache[cid];
	return (*c->getElt)(c->pData,uIndex,c->iDataSize);
	}
    iInLine = uIndex & c->nLineMask;
    uLine = uIndex >> c->nLineBits;

    /* First check our own cache */
    uHash = (MurmurHash2(uLine,tuId)&arc->uHashMask);
    for ( temp = arc->Hash[uHash]; temp; temp = temp->coll ) {
	if (temp->uPage == uLine && (temp->uId&_IDMASK_) == tuId) break;
	}
    if (temp != NULL) {                       /* found in cache directory? */
	switch (temp->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	    temp->uId = uId;     /* clears prefetch flag and sets WHERE = _T1_ and dirty bit */
/* 	    assert( (temp->uId&_WHERE_) == _T1_ ); */
	    temp->remove_from_list();                           /* take off T1 list */
	    arc->T1->mru_insert(temp);                         /* first access but recently prefetched, put on T1 */
	    goto cachehit;
	case _T1_:
	    arc->T1Length--; arc->T2Length++;
	    temp->uId |= _T2_;
/* 	    assert( (temp->uId&_WHERE_) == _T2_ ); */
	    /* fall through */
	case _T2_:
	    temp->remove_from_list();                                   /* take off whichever list */
	    arc->T2->mru_insert(temp);                                     /* seen twice recently, put on T2 */
	    temp->uId |= uId;          /* if the dirty bit is now set we need to record this */
	cachehit:
	    if (bLock) {
		/*
		** We don't have to check if the lock counter rolls over, since it will increment the  
		** magic number first. This in turn will cause this page to be locked in a way that it 
		** cannot be unlocked without an error condition.
		*/
		++temp->data[-1];       /* increase lock count */
		}
	    /*
	    ** Get me outa here.
	    */
	    return (char *)temp->data + c->iDataSize*iInLine;
	case _B1_:                            /* B1 hit: favor recency */
	    /*
	    ** Can initiate the data request right here, and do the rest while waiting...
	    */
	    if (!bVirtual) queueCacheRequest(cid,uLine,uId);
/* 	    assert(arc->B1->next != arc->B1); */
/* 	    assert(arc->B1Length>0); */
	    rat = arc->B2Length/arc->B1Length;
	    if (rat < 1) rat = 1;
	    arc->target_T1 += rat;
	    if (arc->target_T1 > arc->nCache) arc->target_T1 = arc->nCache;
	    /* adapt the target size */
	    arc->B1Length--;                                           /* bookkeep */
	    goto doBcase;
	case _B2_:                            /* B2 hit: favor frequency */
	    /*
	    ** Can initiate the data request right here, and do the rest while waiting...
	    */
	    if (!bVirtual) queueCacheRequest(cid,uLine,uId);
/* 	    assert(arc->B2->next != arc->B2); */
/* 	    assert(arc->B2Length>0); */

	    rat = arc->B1Length/arc->B2Length;
	    if (rat < 1) rat = 1;
	    if (rat > arc->target_T1) arc->target_T1 = 0;
	    else arc->target_T1 = arc->target_T1 - rat;
	    /* adapt the target size */
	    arc->B2Length--;                                           /* bookkeep */
	    inB2=true;
	doBcase:
	    temp->remove_from_list();                                   /* take off whichever list */
	    temp->data = arc->replace(inB2);                                /* find a place to put new page */
	    temp->uId = _T2_|uId;     /* temp->ARC_where = _T2_; and set the dirty bit for this page */
/* 	    assert( (temp->uId&_WHERE_) == _T2_ ); */
	    temp->uPage = uLine;                          /* bookkeep */
	    arc->T2->mru_insert(temp);                                     /* seen twice recently, put on T2 */
	    arc->T2Length++;                 /* JS: this was not in the original code. Should it be? bookkeep */
	    /*assert(temp->data!=NULL);*/
	    finishCacheRequest(cid,uLine,uId,temp,bVirtual);
	    break;
	    }
	}

    else {                                                              /* page is not in cache directory */
	++c->nMiss;
	TimeAddComputing();
	/*
	** Can initiate the data request right here, and do the rest while waiting...
	*/
	if (!bVirtual) queueCacheRequest(cid,uLine,uId);
	L1Length = arc->T1Length + arc->B1Length;
	/*assert(L1Length<=arc->nCache);*/
	if (L1Length == arc->nCache) {                                   /* B1 + T1 full? */
	    if (arc->T1Length < arc->nCache) {                                           /* Still room in T1? */
		temp = arc->remove_from_hash(arc->B1->lru_remove());        /* yes: take page off B1 */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = arc->replace();                                /* find new place to put page */
		}
	    else {                                                      /* no: B1 must be empty */
		temp = arc->remove_from_hash(arc->T1->lru_remove());       /* take page off T1 */
		arc->destage(temp);     /* if dirty, evict before overwrite */
		arc->T1Length--;                                               /* bookkeep that */
		}
	    /*assert(temp->data!=NULL);*/
	    }
	else {                                                          /* B1 + T1 have less than c pages */
	    uint32_t nCache = arc->T1Length + arc->T2Length + arc->B1Length + arc->B2Length;
	    /*assert(L1Length < arc->nCache);*/
	    if (nCache >= arc->nCache) {         /* cache full? */
		/* Yes, cache full: */
		if (nCache == 2*arc->nCache) {
		    /* directory is full: */
		    temp = arc->remove_from_hash(arc->B2->lru_remove());
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=true;
		} else {                                                   /* cache directory not full, easy case */
		    temp = arc->Free->lru_remove();
		    assert(temp->data == NULL);            /* This CDB should not be associated with data */
		}
		temp->data = arc->replace(inB2);                                /* new place for page */
		/*assert(temp->data!=NULL);*/
	    } else {                                                      /* cache not full, easy case */
		temp = arc->Free->lru_remove();
		assert(temp->data != NULL);               /* This CDB should have an unused page associated with it */
		temp->data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
		}
	    }
	arc->T1->mru_insert(temp);                                             /* seen once recently, put on T1 */
	arc->T1Length++;                                                       /* bookkeep: */
	/*assert(arc->T1Length+arc->B1Length<=arc->nCache);*/
	temp->uId = uId;                  /* temp->dirty = dirty;  p->ARC_where = _T1_; as well! */
/* 	assert( (temp->uId&_WHERE_) == _T1_ ); */
	temp->uPage = uLine;
	finishCacheRequest(cid,uLine,uId,temp,bVirtual);
	temp->coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */
	TimeAddWaiting();
    }
    /*assert(temp!=NULL);*/
    if (bLock) {
	/*
	** We don't have to check if the lock counter rolls over, since it will increment the  
	** magic number first. This in turn will cause this page to be locked in a way that it 
	** cannot be unlocked without an error condition.
	*/
	++temp->data[-1];       /* increase lock count */
    }
    /* If we will modify the element, then it must eventually be flushed. */
    if (bModify) temp->uId |= _DIRTY_;

    return (char *)temp->data + c->iDataSize*iInLine;
}

/* Does not lock the element */
void *mdlFetch(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 0;  /* we never lock in fetch */
    const int modify = 0; /* fetch can never modify */
    const int virt = 0; /* really fetch the element */
    return reinterpret_cast<mdlClass *>(mdl)->Access(cid, iIndex, id, lock, modify, virt);
    }

/* Locks, so mdlRelease must be called eventually */
void *mdlAcquire(MDL cmdl,int cid,int iIndex,int id) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    const int lock = 1;  /* we always lock in acquire */
    const int modify = (mdl->cache[cid].iType == MDL_COCACHE);
    const int virt = 0; /* really fetch the element */
    return mdl->Access(cid, iIndex, id, lock, modify, virt);
    }

/* Locks the element, but does not fetch or initialize */
void *mdlVirtualFetch(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 0; /* fetch never locks */
    const int modify = 1; /* virtual always modifies */
    const int virt = 1; /* do not fetch the element */
    return reinterpret_cast<mdlClass *>(mdl)->Access(cid, iIndex, id, lock, modify, virt);
    }

void mdlRelease(MDL cmdl,int cid,void *p) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    CACHE *c = &mdl->cache[cid];
    arcRelease(c->arc,CAST(uint64_t *,p));
    }


double mdlNumAccess(MDL cmdl,int cid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    CACHE *c = &mdl->cache[cid];
    return(c->nAccess);
    }


double mdlMissRatio(MDL cmdl,int cid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccess;

    if (dAccess > 0.0) return(c->nMiss/dAccess);
    else return(0.0);
    }


double mdlCollRatio(MDL cmdl,int cid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccess;

    if (dAccess > 0.0) return(c->nColl/dAccess);
    else return(0.0);
    }

/*
** GRID Geometry information.  The basic process is as follows:
** - Initialize: Create a MDLGRID giving the global geometry information (total grid size)
** - SetLocal:   Set the local grid geometry (which slabs are on this processor)
** - GridShare:  Share this information between processors
** - Malloc:     Allocate one or more grid instances
** - Free:       Free the memory for all grid instances
** - Finish:     Free the GRID geometry information.
*/
void mdlGridInitialize(MDL cmdl,MDLGRID *pgrid,int n1,int n2,int n3,int a1) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    MDLGRID grid;
    assert(n1>0&&n2>0&&n3>0);
    assert(n1<=a1);
    *pgrid = grid = CAST(mdlGridContext *,malloc(sizeof(struct mdlGridContext))); assert(grid!=NULL);
    grid->n1 = n1;
    grid->n2 = n2;
    grid->n3 = n3;
    grid->a1 = a1;

    /* This will be shared later (see mdlGridShare) */
    grid->id = CAST(uint32_t *,malloc(sizeof(*grid->id)*(grid->n3)));    assert(grid->id!=NULL);
    grid->rs = CAST(uint32_t *,mdlMalloc(cmdl,sizeof(*grid->rs)*mdl->Procs())); assert(grid->rs!=NULL);
    grid->rn = CAST(uint32_t *,mdlMalloc(cmdl,sizeof(*grid->rn)*mdl->Procs())); assert(grid->rn!=NULL);

    /* The following need to be set to appropriate values still. */
    grid->sSlab = grid->nSlab = 0;
    grid->nLocal = 0;
    }

void mdlGridFinish(MDL cmdl, MDLGRID grid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    if (grid->rs) free(grid->rs);
    if (grid->rn) free(grid->rn);
    if (grid->id) free(grid->id);
    free(grid);
    }

void mdlGridSetLocal(MDL cmdl,MDLGRID grid,int s, int n, uint64_t nLocal) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    assert( s>=0 && s<grid->n3);
    assert( n>=0 && s+n<=grid->n3);
    grid->sSlab = s;
    grid->nSlab = n;
    grid->nLocal = nLocal;
    }

/*
** Share the local GRID information with other processors by,
**   - finding the starting slab and number of slabs on each processor
**   - building a mapping from slab to processor id.
*/
void mdlGridShare(MDL cmdl,MDLGRID grid) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    int i, id;
    gridShare share;

    share.grid = grid;
    mdl->SendToMPI(&share,MDL_SE_GRID_SHARE);
    mdl->WaitThreadQueue(0);

    /* Calculate on which processor each slab can be found. */
    for(id=0; id<mdl->Procs(); id++ ) {
	for( i=grid->rs[id]; i<grid->rs[id]+grid->rn[id]; i++ ) grid->id[i] = id;
	}
    }

/*
** Retrieve the first and last element for the calling thread.
*/
void mdlGridCoordFirstLast(MDL cmdl,MDLGRID grid,mdlGridCoord *f,mdlGridCoord *l,int bCacheAlign) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    uint64_t nPerCore, nThisCore;
    uint64_t nLocal = (uint64_t)(grid->a1) * grid->n2 * grid->nSlab;

    /* Number on each core with multiples of "a1" elments (complete pencils). */
    /* This needs to change to be complete MDL "cache lines" at some point. */
    uint64_t nAlign;
    if (bCacheAlign) nAlign = 1 << (int)(log2(MDL_CACHE_DATA_SIZE / sizeof(float)));
    else nAlign = grid->a1;

    /*nPerCore = nLocal / mdlCores(mdl) + nAlign - 1;*/
    nPerCore = (nLocal-1) / mdlCores(mdl) + nAlign;
    nPerCore -= nPerCore % nAlign;

    if ( mdlCore(mdl)*nPerCore >= nLocal) nThisCore = 0;
    else if (mdlCore(mdl) == mdlCores(mdl)-1) nThisCore = nLocal - (mdlCores(mdl)-1)*nPerCore;
    else if ( (1+mdlCore(mdl))*nPerCore < nLocal) nThisCore = nPerCore;
    else nThisCore = nLocal - mdlCore(mdl)*nPerCore;

    /* Calculate global x,y,z coordinates, and local "i" coordinate. */
    f->II = nPerCore * mdlCore(mdl);
    if (f->II > nLocal) f->II = nLocal;
    l->II = f->II + nThisCore;
    f->i = 0;
    l->i = nThisCore;
    f->z = f->II/(grid->a1*grid->n2);
    f->y = f->II/grid->a1 - f->z * grid->n2;
    f->x = f->II - grid->a1 * (f->z*grid->n2 + f->y);
    assert(bCacheAlign || f->x == 0); // MDL depends on this at the moment
    f->z += grid->sSlab;
    f->grid = grid;

    l->z = l->II/(grid->a1*grid->n2);
    l->y = l->II/grid->a1 - l->z * grid->n2;
    l->x = l->II - grid->a1 * (l->z*grid->n2 + l->y);
    assert(bCacheAlign || l->x == 0); // MDL depends on this at the moment
    l->z += grid->sSlab;
    l->grid = grid;
    }

/*
** Allocate the local elements.  The size of a single element is
** given and the local GRID information is consulted to determine
** how many to allocate.
*/
void *mdlGridMalloc(MDL mdl,MDLGRID grid,int nEntrySize) {
    return mdlMallocArray(mdl,mdlCore(mdl)?0:grid->nLocal,nEntrySize,0);
    }

void mdlGridFree( MDL mdl, MDLGRID grid, void *p ) {
    mdlFree(mdl,p);
    }

#ifdef MDL_FFTW
size_t mdlFFTlocalCount(MDL cmdl,int n1,int n2,int n3,int *nz,int *sz,int *ny,int*sy) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    fftSizes sizes;
    sizes.n1 = n1;
    sizes.n2 = n2;
    sizes.n3 = n3;
    mdl->SendToMPI(&sizes,MDL_SE_FFT_SIZES);
    mdl->WaitThreadQueue(0);
    if (nz!=NULL) *nz = sizes.nz;
    if (sz!=NULL) *sz = sizes.sz;
    if (ny!=NULL) *ny = sizes.ny;
    if (sy!=NULL) *sy = sizes.sy;
    return sizes.nLocal;
    }

MDLFFT mdlFFTNodeInitialize(MDL cmdl,int n1,int n2,int n3,int bMeasure,FFTW3(real) *data) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    fftPlans plans;
    MDLFFT fft = CAST(mdlFFTContext *,malloc(sizeof(struct mdlFFTContext)));
    assert( fft != NULL);
    assert(mdlCore(mdl) == 0);
    plans.sizes.n1 = n1;
    plans.sizes.n2 = n2;
    plans.sizes.n3 = n3;
    plans.data = 0;/*Estimate is faster? data;*/
    plans.kdata = 0;/*(FFTW3(complex) *)data;*/
    mdl->SendToMPI(&plans,MDL_SE_FFT_PLANS);
    mdl->WaitThreadQueue(0);
    fft->fplan = plans.fplan;
    fft->iplan = plans.iplan;

    /*
    ** Dimensions of k-space and r-space grid.  Note transposed order.
    ** Note also that the "actual" dimension 1 side of the r-space array
    ** can be (and usually is) larger than "n1" because of the inplace FFT.
    */
    mdlGridInitialize(cmdl,&fft->rgrid,n1,n2,n3,2*(n1/2+1));
    mdlGridInitialize(cmdl,&fft->kgrid,n1/2+1,n3,n2,n1/2+1);
    mdlGridSetLocal(cmdl,fft->rgrid,plans.sizes.sz,plans.sizes.nz,plans.sizes.nLocal);
    mdlGridSetLocal(cmdl,fft->kgrid,plans.sizes.sy,plans.sizes.ny,plans.sizes.nLocal/2);
    mdlGridShare(cmdl,fft->rgrid);
    mdlGridShare(cmdl,fft->kgrid);
    return fft;
    }

MDLFFT mdlFFTInitialize(MDL cmdl,int n1,int n2,int n3,int bMeasure,FFTW3(real) *data) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    MDLFFT fft;
    if (mdlCore(mdl) == 0) {
	mdl->pvMessageData = mdlFFTNodeInitialize(cmdl,n1,n2,n3,bMeasure,data);
	}
    mdl->ThreadBarrier(); /* To synchronize pvMessageData */
    fft = CAST(MDLFFT,mdl->pmdl[0]->pvMessageData);
    mdl->ThreadBarrier(); /* In case we reuse pvMessageData */
    return fft;
    }

void mdlFFTNodeFinish( MDL cmdl, MDLFFT fft ) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    FFTW3(destroy_plan)(fft->fplan);
    FFTW3(destroy_plan)(fft->iplan);
    mdlGridFinish(cmdl,fft->kgrid);
    mdlGridFinish(cmdl,fft->rgrid);
    free(fft);
    }
void mdlFFTFinish( MDL cmdl, MDLFFT fft ) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    mdl->ThreadBarrier();
    if (mdlCore(mdl) == 0) {
	mdlFFTNodeFinish(cmdl,fft);
	}
    }

FFTW3(real) *mdlFFTMalloc( MDL mdl, MDLFFT fft ) {
    return CAST(FFTW3(real) *,mdlGridMalloc(mdl,fft->rgrid,sizeof(FFTW3(real))));
    }

void mdlFFTFree( MDL mdl, MDLFFT fft, void *p ) {
    mdlGridFree(mdl,fft->rgrid,p);
    }

void mdlFFT( MDL cmdl, MDLFFT fft, FFTW3(real) *data ) { reinterpret_cast<mdlClass *>(cmdl)->FFT(fft,data); }
void mdlClass::FFT( MDLFFT fft, FFTW3(real) *data ) {
    fftTrans trans;
    ThreadBarrier();
    if (Core() == iCoreMPI) {
	FFTW3(execute_dft_r2c)(fft->fplan,data,(FFTW3(complex) *)(data));
	}
    else if (Core() == 0) {
	trans.fft = fft;
	trans.data = data;
	trans.kdata = (FFTW3(complex) *)data;
	SendToMPI(&trans,MDL_SE_FFT_DFT_R2C);
	}
    if (Cores()>1) pthread_barrier_wait(&pmdl[0]->barrier);
    ThreadBarrier();
    }

void mdlIFFT( MDL cmdl, MDLFFT fft, FFTW3(complex) *kdata ) { reinterpret_cast<mdlClass *>(cmdl)->IFFT(fft,kdata); }
void mdlClass::IFFT( MDLFFT fft, FFTW3(complex) *kdata ) {
    fftTrans trans;
    ThreadBarrier();
    if (Core() == iCoreMPI) {
	FFTW3(execute_dft_c2r)(fft->iplan,kdata,(FFTW3(real) *)(kdata));
	}
    else if (Core() == 0) {
	trans.fft = fft;
	trans.kdata = kdata;
	trans.data = (FFTW3(real) *)kdata;
	SendToMPI(&trans,MDL_SE_FFT_DFT_C2R);
	}
    if (Cores()>1) pthread_barrier_wait(&pmdl[0]->barrier);
    ThreadBarrier();
    }
#endif

void mdlAlltoallv(MDL cmdl,int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    mdl->Alltoallv(dataSize,sbuff,scount,sdisps,rbuff,rcount,rdisps);
    }

void mdlClass::Alltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps) {
    alltoallv a2a;
    assert(Core() == 0);
    a2a.sbuff = sbuff;
    a2a.scount = scount;
    a2a.sdisps = sdisps;
    a2a.rbuff = rbuff;
    a2a.rcount = rcount;
    a2a.rdisps = rdisps;
    a2a.dataSize = dataSize;
    SendToMPI(&a2a,MDL_SE_ALLTOALLV);
    WaitThreadQueue(0);
    }
#if defined(USE_CUDA) || defined(USE_CL)
void mdlSetCudaBufferSize(MDL cmdl,int inBufSize, int outBufSize) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    if (mdl->inCudaBufSize < inBufSize) mdl->inCudaBufSize = inBufSize;
    if (mdl->outCudaBufSize < outBufSize) mdl->outCudaBufSize = outBufSize;
    }
#endif

void mdlSetWorkQueueSize(MDL cmdl,int wqMaxSize,int cudaSize) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    MDLwqNode *work;
    int i;

#ifdef USE_CUDA
    if (mdl->cudaCtx) CUDA_SetQueueSize(mdl->cudaCtx,cudaSize,mdl->inCudaBufSize,mdl->outCudaBufSize);
#endif
#ifdef USE_CL
    CL_SetQueueSize(mdl->clCtx,cudaSize,mdl->inCudaBufSize,mdl->outCudaBufSize);
#endif
    while (wqMaxSize > mdl->wqMaxSize) {
	for(i=0; i<mdl->Cores(); ++i) {
	    work = CAST(MDLwqNode *,malloc(sizeof(MDLwqNode)));
	    OPA_Queue_header_init(&work->q.hdr);
	    work->iCoreOwner = mdl->Core();
	    work->ctx = NULL;
	    work->doFcn = NULL;
	    work->doneFcn = NULL;
	    OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	    }
	++mdl->wqMaxSize;
	}
    while (wqMaxSize < mdl->wqMaxSize) {
	for(i=0; i<mdl->Cores(); ++i) {
	    if (!OPA_Queue_is_empty(&mdl->wqFree)) 
		OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	    else assert(0);
	    free(work);
	    }
	--mdl->wqMaxSize;
	} 
    }

void mdlAddWork(MDL cmdl, void *ctx,
    int (*initWork)(void *ctx,void *vwork),
    int (*checkWork)(void *ctx,void *vwork),
    mdlWorkFunction doWork,
    mdlWorkFunction doneWork) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    mdlClass *Qmdl = NULL;
    MDLwqNode *work;
    int i;

    /* We prefer to let CUDA do the work */
#ifdef USE_CUDA
    assert(initWork==NULL && checkWork==NULL); // obsolete
//    if (CUDA_queue(mdl->cudaCtx,ctx,initWork,checkWork)) return;
#endif
    /* Obviously, we can only queue work if we have a free queue element */
    if (!OPA_Queue_is_empty(&mdl->wqFree)) {
	/* We have some room, so save work for later */
	if (OPA_load_int(&mdl->wqCurSize) < mdl->wqMaxSize) Qmdl = mdl;
	/* See if anyone else is accepting work */
	else {
	    i = mdl->wqLastHelper;
	    do {
		mdlClass * rMDL = mdl->pmdl[i];
		if (rMDL->wqAccepting && OPA_load_int(&rMDL->wqCurSize)<rMDL->wqMaxSize) {
		    Qmdl = rMDL;
		    mdl->wqLastHelper = i;
		    break;
		    }
		if (++i == mdl->Cores()) i = 0;
		} while(i!=mdl->wqLastHelper);
	    }
	if (Qmdl) {
	    OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	    work->ctx = ctx;
	    work->doFcn = doWork;
	    work->doneFcn = doneWork;
	    OPA_incr_int(&Qmdl->wqCurSize);
	    OPA_Queue_enqueue(&Qmdl->wq, work, MDLwqNode, q.hdr);
	    return;
	    }
	}

    /* Just handle it ourselves */
    while( doWork(ctx) != 0 ) {}
    doneWork(ctx);
    }

int mdlProcToThread(MDL cmdl, int iProc) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    return mdl->ProcToThread(iProc);
    }
int mdlThreadToProc(MDL cmdl, int iThread) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    return mdl->ThreadToProc(iThread);
    }

int mdlThreads(void *mdl) {  return reinterpret_cast<mdlClass *>(mdl)->Threads(); }
int mdlSelf(void *mdl)    {  return reinterpret_cast<mdlClass *>(mdl)->Self(); }
int mdlCore(void *mdl)    {  return reinterpret_cast<mdlClass *>(mdl)->Core(); }
int mdlCores(void *mdl)   {  return reinterpret_cast<mdlClass *>(mdl)->Cores(); }
int mdlProc(void *mdl)    {  return reinterpret_cast<mdlClass *>(mdl)->Proc(); }
int mdlProcs(void *mdl)   {  return reinterpret_cast<mdlClass *>(mdl)->Procs(); }
int mdlGetArgc(void *mdl) {  return reinterpret_cast<mdlClass *>(mdl)->argc; }
char **mdlGetArgv(void *mdl) {  return reinterpret_cast<mdlClass *>(mdl)->argv; }


void mdlTimeReset(MDL mdl)             {        reinterpret_cast<mdlClass *>(mdl)->TimeReset(); }
double mdlTimeComputing(MDL mdl)       { return reinterpret_cast<mdlClass *>(mdl)->TimeComputing(); }
double mdlTimeSynchronizing(MDL mdl)   { return reinterpret_cast<mdlClass *>(mdl)->TimeSynchronizing(); }
double mdlTimeWaiting(MDL mdl)         { return reinterpret_cast<mdlClass *>(mdl)->TimeWaiting(); }

void mdlprintf(MDL cmdl, const char *format, ...) {
    mdlClass *mdl = reinterpret_cast<mdlClass *>(cmdl);
    va_list args;
    va_start(args, format);
    mdl->mdl_vprintf(format,args);
    va_end(args);
    }