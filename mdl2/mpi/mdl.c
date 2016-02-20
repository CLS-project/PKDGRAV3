/*
 ** MPI version of MDL.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#if !defined(HAVE_CONFIG_H) || defined(HAVE_MALLOC_H)
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
#include "mdl.h"
#ifdef USE_CUDA
#include "cudautil.h"
#endif
#ifdef USE_ITT
#include "ittnotify.h"
#endif

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_DEFAULT_CACHEIDS	10

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

/*
** MurmurHash2, by Austin Appleby
** adapted for hashing 2 uint32_t variables for mdl2
*/
static inline uint32_t MurmurHash2(uint32_t a,uint32_t b)
{
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

static inline CDB *remove_from_list(CDB *t) {
    t->hdr.links.prev->hdr.links.next = t->hdr.links.next;
    t->hdr.links.next->hdr.links.prev = t->hdr.links.prev;
    return t;
}

static inline CDB *lru_remove(CDB *list) {
    return remove_from_list(list->hdr.links.prev);
}

static inline void lru_insert(CDB *p,CDB *list) {
    p->hdr.links.prev = list->hdr.links.prev;
    p->hdr.links.next = list;
    list->hdr.links.prev->hdr.links.next = p;
    list->hdr.links.prev = p;
}

static inline void mru_insert(CDB *p,CDB *list) {
    lru_insert(p,list->hdr.links.next);
    }

static inline CDB * remove_from_hash(ARC arc,CDB *p) {
    uint32_t uPage = p->uPage;
    uint32_t uId = p->uId&_IDMASK_;
    uint32_t uHash = (MurmurHash2(uPage,uId)&arc->uHashMask);
    CDB **pt;

    for(pt = &arc->Hash[uHash]; *pt != NULL; pt = &((*pt)->extra.coll)) {
	if ( *pt == p) {
	    *pt = (*pt)->extra.coll;
	    return p;
	    }
	}
    printf("Tried to remove uPage %d uId %x\n", uPage, uId);
    abort();  /* should never get here, the element should always be found in the hash */
    }

ARC arcInitialize(uint32_t nCache,uint32_t uDataSize,CACHE *c) {
    ARC arc;
    uint32_t i;
    
    /*
    ** Allocate stuff.
    */
    arc = malloc(sizeof(struct ArcContext));
    assert(arc != NULL);
    arc->cache = c;
    arc->nCache = nCache;
    arc->cdbBase = malloc(2*nCache*sizeof(CDB));
    assert(arc->cdbBase != NULL);
    /*
    ** Make sure we have sufficient alignment of data.
    ** In this case to nearest long word (8 bytes).
    ** We add one long word at the start to store the 
    ** magic number and lock count.
    */
    arc->uDataSize = (uDataSize+7)>>3;
    arc->dataBase = malloc(nCache*sizeof(uint64_t)*(arc->uDataSize+1));
    assert(arc->dataBase != NULL);
    arc->dataLast = arc->dataBase + nCache*(arc->uDataSize+1);
    /*
    ** Determine nHash.
    */
    arc->uHashMask = swar32(3*nCache-1);
    arc->nHash = arc->uHashMask+1; 
    arc->Hash = malloc(arc->nHash*sizeof(CDB *));
    assert(arc->Hash != NULL);
    for (i=0;i<arc->nHash;++i) {
	arc->Hash[i] = NULL;
    }
    /*
    ** Create all lists as circular lists, with a sentinel
    ** CDB at the head/tail of the list.
    ** Initialize the lengths of the various lists.
    */
    arc->T1 = malloc(sizeof(CDB));
    assert(arc->T1 != NULL);
    arc->T1->hdr.links.next = arc->T1;
    arc->T1->hdr.links.prev = arc->T1;
    arc->T1->uId = 0xdeadbeef;
    arc->T1Length = 0;
    arc->B1 = malloc(sizeof(CDB));
    assert(arc->B1 != NULL);
    arc->B1->hdr.links.next = arc->B1;
    arc->B1->hdr.links.prev = arc->B1;
    arc->B1->uId = 0xdeadbeef;
    arc->B1Length = 0;
    arc->T2 = malloc(sizeof(CDB));
    assert(arc->T2 != NULL);
    arc->T2->hdr.links.next = arc->T2;
    arc->T2->hdr.links.prev = arc->T2;
    arc->T2->uId = 0xdeadbeef;
    arc->T2Length = 0;
    arc->B2 = malloc(sizeof(CDB));
    assert(arc->B2 != NULL);
    arc->B2->hdr.links.next = arc->B2;
    arc->B2->hdr.links.prev = arc->B2;
    arc->B2->uId = 0xdeadbeef;
    arc->B2Length = 0;
    arc->Free = malloc(sizeof(CDB));
    assert(arc->Free != NULL);
    arc->Free->hdr.links.next = arc->Free;
    arc->Free->hdr.links.prev = arc->Free;
    /*
    ** Initialize target T1 length.
    */
    arc->target_T1 = nCache/2;   /* is this ok? */
    /*
    ** Insert CDBs with data into the Free list first.
    */
    for (i=0;i<nCache;++i) {
	/*arc->dataBase[i*(arc->uDataSize+1)] = _ARC_MAGIC_;*/ /* Defer this until we use it. */
	arc->cdbBase[i].data = &arc->dataBase[i*(arc->uDataSize+1)+1];
	arc->cdbBase[i].extra.coll = NULL;
	lru_insert(&arc->cdbBase[i],arc->Free);
    }
    /*
    ** Finally insert CDBs without data pages into the Free list last.
    */
    for (i=nCache;i<2*nCache;++i) {
	arc->cdbBase[i].data = 0;
	arc->cdbBase[i].extra.coll = NULL;
	mru_insert(&arc->cdbBase[i],arc->Free);
    }
    return(arc);
}


void arcFinish(ARC arc) {
    /*
    ** Free the sentinels.
    */
    free(arc->Free);
    free(arc->B2);
    free(arc->T2);
    free(arc->B1);
    free(arc->T1);
    /*
    ** Free the hash table.
    */
    free(arc->Hash);
    /*
    ** Free the data pages.
    */
    free(arc->dataBase);
    /*
    ** Free the CDBs.
    */
    free(arc->cdbBase);
    /*
    ** Free context.
    */
    free(arc);
}

ARC arcReinitialize(ARC arc,uint32_t nCache,uint32_t uDataSize,CACHE *c) {
    if (arc!=NULL) {
	assert(arc->cache == c);
	if ( arc->nCache == nCache && arc->uDataSize == ((uDataSize+7)>>3) )
	    return arc;
	arcFinish(arc);
	}
    return arcInitialize(nCache,uDataSize,c);
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


static void mdlSendThreadMessage(MDL mdl,int iQueue,int iCore, void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = mdl->base.iCore;
    OPA_Queue_enqueue(mdl->pmdl[iCore]->inQueue+iQueue, qhdr, MDLserviceElement, hdr);
    }

static void mdlSendToMPI(MDL mdl,void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = mdl->base.iCore;
    OPA_Queue_enqueue(&mdl->pmdl[mdl->iCoreMPI]->mpi->queueMPI, qhdr, MDLserviceElement, hdr);
    }

/* Accept pending combine requests, and call the combine function for each. */
static void combine_all_incoming(MDL mdl) {
    mdlContextMPI *mpi = mdl->pmdl[mdl->iCoreMPI]->mpi;
    if (mdl->base.iCore<0) return; /* Dedicated MPI thread combines nothing */
    while (!OPA_Queue_is_empty(&mdl->wqCacheFlush)) {
	MDLflushBuffer *flush;
	CACHE *c;
	CAHEAD *ca;
	char *pData;
	int count, i, s, n;
	uint32_t uIndex;

	OPA_Queue_dequeue(&mdl->wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	count = flush->nBytes;
    	ca = (CAHEAD *)(flush + 1);
	while(count > 0) {
	    pData = (char *)(ca+1);
	    c = &mdl->cache[ca->cid];
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

static void mpi_finish_local_flush(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
    MDLflushBuffer *flush;
    int iCore;

    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */
    for(iCore=0; iCore<mdlCores(mdl); ++iCore) {
	MDL mdl1 = mdl->pmdl[iCore];
	if ((flush=mpi->flushBuffersByCore[iCore])) {
	    OPA_Queue_enqueue(&mdl1->wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	    mpi->flushBuffersByCore[iCore] = NULL;
	    }
	}
    }

static MDLflushBuffer *get_local_flush_buffer(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
    MDLflushBuffer *flush;

    /* This is bad. We just have to wait and hope for the best. */
    while (OPA_Queue_is_empty(&mpi->localFlushBuffers)) {
	/*
	** We have no local flush buffers, but threads should wakeup soon and give them back.
	** In the mean time we need to flush buffers ourselves (if we are not a dedicated MPI).
	*/
	combine_all_incoming(mdl);
	}
    OPA_Queue_dequeue(&mpi->localFlushBuffers, flush, MDLflushBuffer, hdr.hdr);
    flush->nBytes = 0;
    return flush;
    }

static void mpi_queue_local_flush(MDL mdl,CAHEAD *ph) {
    mdlContextMPI *mpi = mdl->mpi;
    int iCore = ph->idTo;
    MDL mdl1 = mdl->pmdl[iCore];
    CACHE *c = &mdl1->cache[ph->cid];
    MDLflushBuffer *flush;
    int nLeft, nToAdd;

    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */
    assert(iCore < mdlCores(mdl));

    if ((flush=mpi->flushBuffersByCore[iCore]) == NULL) {
	flush = mpi->flushBuffersByCore[iCore] = get_local_flush_buffer(mdl);
	}

    nLeft = flush->nBufferSize - flush->nBytes;
    nToAdd = sizeof(CAHEAD) + c->iLineSize;

    /* No room left? Send this buffer and grab a new one for this rank */
    if (nToAdd > nLeft ) {
	OPA_Queue_enqueue(&mdl1->wqCacheFlush, flush, MDLflushBuffer, hdr.hdr);
	flush = mpi->flushBuffersByCore[iCore] = get_local_flush_buffer(mdl);
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

int mdlCacheReceive(MDL mdl,MPI_Status *status) {
    mdlContextMPI *mpi = mdl->mpi;
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
    iCore = ph->idTo - mdl->pmdl[0]->base.idSelf;
    assert(iCore>=0 && iCore<mdl->base.nCores);
    c = &mdl->pmdl[iCore]->cache[ph->cid];
    if (c->iType == MDL_NOCACHE) printf("Cache %d\n",ph->cid);

    assert(c->iType != MDL_NOCACHE);

    switch (ph->mid) {
	/* A different process wants local cache data - we can grab it directly */
    case MDL_MID_CACHEREQ:
	assert( count == sizeof(CAHEAD) );
	iProc = mdlThreadToProc(mdl,ph->idFrom); /* Can use iRankFrom */
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
	    pdata = malloc(sizeof(MDLcacheReplyData) + mpi->iReplyBufSize);
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
	phRpl->idFrom = mdl->base.idSelf;

	/* Caclulation based on threads */
	phRpl->idTo = ph->idFrom;
	phRpl->iLine = ph->iLine;
	phRpl->nItems = 1;

	s = ph->iLine << c->nLineBits;
	n = s + c->nLineElements;
	iLineSize = c->iLineSize;
	for(i=s; i<n; i++ ) {
	    if (i<c->nData) {
		t = (*c->getElt)(c->pData,i,c->iDataSize);
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
	    iCore = ph->idTo - mdl->pmdl[0]->base.idSelf;
	    c = &mdl->pmdl[iCore]->cache[ph->cid];
	    iIndex = ph->iLine << c->nLineBits;
	    while(iIndex >= c->nData) {
		iIndex -= c->nData;
		assert(iCore+1<mdlCores(mdl));
		c = &mdl->pmdl[++iCore]->cache[ph->cid];
		}
	    ph->iLine = iIndex >> c->nLineBits;
	    ph->idTo = iCore;
	    mpi_queue_local_flush(mdl,ph);
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
	pLine = creq->pLine;
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
	mdlSendThreadMessage(mdl,MDL_TAG_CACHECOM,iCore,creq,MDL_SE_CACHE_REQUEST);
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
static int flush_check_completion(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
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

/* Send the given buffer. Remove it from whichever list and put it on the "sent" list. */
static void flush_send(MDL mdl,MDLflushBuffer *pBuffer) {
    mdlContextMPI *mpi = mdl->mpi;
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

static void mpi_flush_element(MDL mdl,CAHEAD *pHdr,int iLineSize) {
    mdlContextMPI *mpi = mdl->mpi;
    int iProc = mdlThreadToProc(mdl,pHdr->idTo);
    MDLflushBuffer *pBuffer;
    int iCore;

    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */

    /* If we have a buffer already then add to it if we can */
    if ( (pBuffer=mpi->flushBuffersByRank[iProc]) != NULL ) {
	int nLeft = pBuffer->nBufferSize - pBuffer->nBytes;
	int nToAdd = sizeof(CAHEAD) + iLineSize;
	/* No room left? Send this buffer and grab a new one for this rank */
	if (nToAdd > nLeft ) {
	    flush_send(mdl,pBuffer);
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
	    flush_check_completion(mdl);
    	    if (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL) {
            	while(1) {
                    MPI_Test(&mpi->pReqRcv->mpiRequest, &flag, &status);
                    if (flag) mdlCacheReceive(mdl,&status);
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
		flush_send(mdl,flush);
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


static void mpi_flush_elements(MDL mdl,MDLflushBuffer *pFlush) {
    mdlContextMPI *mpi = mdl->mpi;
    CAHEAD *ca;
    char *pData;
    int count;

    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */
    count = pFlush->nBytes;
    ca = (CAHEAD *)(pFlush + 1);
    while(count > 0) {
	int iLineSize = mdl->pmdl[0]->cache[ca->cid].iLineSize;
	pData = (char *)(ca+1);
	mpi_flush_element(mdl,ca,iLineSize);
	pData += iLineSize;
	ca = (CAHEAD *)pData;
	count -= sizeof(CAHEAD) + iLineSize;
	}
    assert(count == 0);
    }

/*
** At the end we need to flush any pending cache elements
*/
static int checkMPI(MDL mdl);
void mpi_flush_all_elements(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */
    while (mpi->flushHeadBusy.hdr.mpi.next != &mpi->flushHeadBusy ) {
	MDLflushBuffer *flush = mpi->flushHeadBusy.hdr.mpi.prev;
	mpi->flushBuffersByRank[flush->iRankTo] = NULL;
	flush_send(mdl,flush);
	}
    while(flush_check_completion(mdl)) {
	checkMPI(mdl);
	}
    }

/*
** This routine must be called often by the MPI thread. It will drain
** any requests from the thread queue, and track MPI completions.
*/
static int checkMPI(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
    assert(mdl->base.iCore == mdl->iCoreMPI); /* MPI thread only */
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
	gridShare *share;
	SRVHEAD *head;
	int iProc, iCore, tag, i;
	int flag,indx;
	MPI_Status status;

	/* These are messages from other threads */
	if (!OPA_Queue_is_empty(&mpi->queueMPI)) {
	    OPA_Queue_dequeue(&mpi->queueMPI, qhdr, MDLserviceElement, hdr);
	    switch(qhdr->iServiceID) {
		/* A thread has finished */
	    case MDL_SE_STOP:
		--mpi->nActiveCores;
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_STOP);
		break;
		/* mdlReqService() on behalf of a thread */
	    case MDL_SE_SEND_REQUEST:
		send = (MDLserviceSend *)qhdr;
		iProc = mdlThreadToProc(mdl,send->target);
		iCore = send->target - mdlProcToThread(mdl,iProc);
		assert(iCore>=0);
		tag = send->tag + MDL_TAG_THREAD_OFFSET * iCore;

		/* Grab a free tag for the reply */
		head = send->buf;
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
		iProc = mdlThreadToProc(mdl,send->target);
		/* tag is really the request ID */
		head = send->buf;
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
		head = send->buf;
		MPI_Irecv(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;

	    case MDL_SE_MPI_SEND:
	    case MDL_SE_MPI_SSEND:
		send = (MDLserviceSend *)qhdr;
		iProc = mdlThreadToProc(mdl,send->target);
		iCore = send->target - mdlProcToThread(mdl,iProc);
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
		else iProc = mdlThreadToProc(mdl,send->target);
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
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_OPEN);
		break;
		/* A thread has closed a cache -- we can cancel the receive if no threads have a cache open */
	    case MDL_SE_CACHE_CLOSE:
		assert(mpi->nOpenCaches > 0);
		--mpi->nOpenCaches;
		assert (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL);
		if (mpi->nOpenCaches == 0) {
		    mpi_flush_all_elements(mdl);
		    MPI_Cancel(&mpi->pReqRcv->mpiRequest);
		    MPI_Wait(&mpi->pReqRcv->mpiRequest, &status);
		    mpi->pReqRcv->mpiRequest = MPI_REQUEST_NULL;
		    }
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
		/* A thread needs to request a missing element */
	    case MDL_SE_CACHE_REQUEST:
		caReq = (MDLserviceCacheReq *)qhdr;
		assert(mpi->pThreadCacheReq[caReq->svc.iCoreFrom]==NULL);
		iProc = mdlThreadToProc(mdl,caReq->caReq.idTo);
		mpi->pThreadCacheReq[caReq->svc.iCoreFrom] = caReq;
		assert(caReq->request==MPI_REQUEST_NULL);
		MPI_Isend(&caReq->caReq,sizeof(CAHEAD),MPI_BYTE,iProc,MDL_TAG_CACHECOM, mpi->commMDL,&caReq->request);
		break;
	    case MDL_SE_CACHE_FLUSH:
		pFlush = (MDLflushBuffer *)qhdr;
		mpi_flush_elements(mdl,pFlush);
		OPA_Queue_enqueue(&mdl->pmdl[qhdr->iCoreFrom]->coreFlushBuffers, pFlush, MDLflushBuffer, hdr.hdr);
		break;
	    case MDL_SE_CACHE_FLUSH_OUT:
		mpi_flush_all_elements(mdl);
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_FLUSH_OUT);
		break;
	    case MDL_SE_CACHE_FLUSH_LCL:
		mpi_finish_local_flush(mdl);
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_FLUSH_LCL);
		break;
#ifdef MDL_FFTW
	    case MDL_SE_FFT_SIZES:
		sizes = (fftSizes *)qhdr;
		sizes->nLocal = 2*FFTW3(mpi_local_size_3d_transposed)
		    (sizes->n3,sizes->n2,sizes->n1/2+1,mpi->commMDL,
		    &sizes->nz,&sizes->sz,&sizes->ny,&sizes->sy);
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
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
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
	    case MDL_SE_FFT_DFT_R2C:
		trans = (fftTrans *)qhdr;
		FFTW3(execute_dft_r2c)(trans->fft->fplan,trans->data,trans->kdata);
		pthread_barrier_wait(&mdl->pmdl[0]->barrier);
		break;
	    case MDL_SE_FFT_DFT_C2R:
		trans = (fftTrans *)qhdr;
		FFTW3(execute_dft_c2r)(trans->fft->iplan,trans->kdata,trans->data);
		pthread_barrier_wait(&mdl->pmdl[0]->barrier);
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
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
	    case MDL_SE_BARRIER_REQUEST:
		MPI_Barrier(mpi->commMDL);
		mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,0,qhdr,MDL_SE_BARRIER_REPLY);
		break;
	    default:
		assert(0);
		}
	    continue;
	    }
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
		    assert(send == &mdl->pmdl[send->svc.iCoreFrom]->sendRequest);
		    tag = 0;
		    break;
		case MDL_SE_RECV_REPLY:
		    head = send->buf;
		    assert(mpi->pRequestTargets[send->target]>=0);
		    mpi->pRequestTargets[send->target] = -1;
		case MDL_SE_MPI_RECV:
		case MDL_SE_RECV_REQUEST:
		    assert(send == &mdl->pmdl[send->svc.iCoreFrom]->recvRequest);
		    tag = send->tag % MDL_TAG_THREAD_OFFSET;
		    MPI_Get_count(&status, MPI_BYTE, &send->count);
		    send->target = status.MPI_SOURCE;
		    break;
		default:
		    assert(0);
		    }
		mdlSendThreadMessage(mdl,tag,send->svc.iCoreFrom,send,MDL_SE_MPI_SSEND);
		continue;
		}
	    }
	if (mpi->pReqRcv->mpiRequest != MPI_REQUEST_NULL) {
	    while(1) {
		MPI_Test(&mpi->pReqRcv->mpiRequest, &flag, &status);
		if (flag) mdlCacheReceive(mdl,&status);
		else break;
		}
	    }
	break;
	}
    return mpi->nActiveCores;
    }

/* Do one piece of work. Return 0 if there was no work. */
static int mdlDoSomeWork(MDL mdl) {
    MDLwqNode *work;
    int rc = 0;
#ifdef NO_TOO_OFTEN_USE_CUDA
    rc = CUDA_flushDone(mdl->cudaCtx);
#endif
    mdlCacheCheck(mdl);
    if (!OPA_Queue_is_empty(&mdl->wq)) {
	/* Grab a work package, and perform it */
	OPA_Queue_dequeue(&mdl->wq, work, MDLwqNode, q.hdr);
	OPA_decr_int(&mdl->wqCurSize);
	while ( (*work->doFcn)(work->ctx) != 0 ) {
	    mdlCacheCheck(mdl);
	    }
	rc = 1;
	/* Send it back to the original thread for completion (could be us) */
	if (work->iCoreOwner == mdl->base.iCore) goto reQueue;
	OPA_Queue_enqueue(&mdl->pmdl[work->iCoreOwner]->wqDone, work, MDLwqNode, q.hdr);
	}
    while (!OPA_Queue_is_empty(&mdl->wqDone)) {
	OPA_Queue_dequeue(&mdl->wqDone, work, MDLwqNode, q.hdr);
    reQueue:
	(*work->doneFcn)(work->ctx);
	work->ctx = NULL;
	work->doFcn = NULL;
	work->doneFcn = NULL;
	OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	}
    return rc;
    }

void mdlCompleteAllWork(MDL mdl) {
#ifdef USE_CUDA
    CUDA_sendWork(mdl->cudaCtx);
#endif
    while(mdlDoSomeWork(mdl)) {}
#ifdef USE_CUDA
    while(CUDA_flushDone(mdl->cudaCtx)) {}
#endif
    }

static void /*MDLserviceElement*/ *mdlWaitThreadQueue(MDL mdl,int iQueue) {
    MDLserviceElement *qhdr;
    while (OPA_Queue_is_empty(mdl->inQueue+iQueue)) {
	if (mdl->base.iCore == mdl->iCoreMPI) checkMPI(mdl);
	combine_all_incoming(mdl);
	if (mdlDoSomeWork(mdl) == 0) {
#ifdef _MSC_VER
	    SwitchToThread();
#else
	    sched_yield();
#endif
	    }
	}
    OPA_Queue_dequeue(mdl->inQueue+iQueue, qhdr, MDLserviceElement, hdr);
    return qhdr;
    }

/* Synchronize threads */
void mdlThreadBarrier(MDL mdl) {
    MDLserviceElement *qhdr;
    int i;

    if (mdl->base.iCore) {
	mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,0,&mdl->inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	for(i=1; i<mdl->base.nCores; ++i) {
	    qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	    assert(qhdr->iServiceID == MDL_SE_BARRIER_REQUEST);
	    }
	for(i=1; i<mdl->base.nCores; ++i) {
	    mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,i,&mdl->pmdl[i]->inMessage,MDL_SE_BARRIER_REPLY);
	    }
	}
    }

static int mdl_start_MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MDL mdl,int stype) {
    int iCore = dest - mdl->pmdl[0]->base.idSelf;
    int bOnNode = (iCore >= 0 && iCore < mdl->base.nCores);
    MDLserviceSend *send = &mdl->sendRequest;

    assert(dest != mdlSelf(mdl));
    send->buf = buf;
    send->count = count;
    send->datatype = datatype;
    send->target = dest;
    send->tag = tag;

    /* To send on node, we give the other thread the message. It will consume it after posting a receive */
    if (bOnNode && stype!=MDL_SE_SEND_REQUEST && stype!=MDL_SE_SEND_REPLY ) {
	mdlSendThreadMessage(mdl,MDL_TAG_MAX + mdl->base.iCore,iCore,send,stype); /* From me */
	}
    /* Pass this off to the MPI thread (which could be ourselves) */
    else {
	mdlSendToMPI(mdl,send,stype);
	}
    return MPI_SUCCESS;
    }

static int mdl_MPI_Ssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MDL mdl) {
    int rc = mdl_start_MPI_Ssend(buf,count,datatype,dest,tag,mdl,MDL_SE_MPI_SSEND);
    mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
    return rc;
    }

static int mdl_MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MDL mdl) {
    int rc = mdl_start_MPI_Ssend(buf,count,datatype,dest,tag,mdl,MDL_SE_MPI_SEND);
    mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
    return rc;
    }

static int mdl_remote_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,
    MDL mdl, int *nBytes, int iServiceID) {
    MDLserviceSend *send = &mdl->recvRequest;
    send->buf = buf;
    send->count = count;
    send->datatype = datatype;
    send->target = source;
    send->tag = tag;
    mdlSendToMPI(mdl,send,iServiceID);
    send = mdlWaitThreadQueue(mdl,tag); /* Wait for the "send" to come back to us. */
    *nBytes = send->count;
    return MPI_SUCCESS;
    }

/*
**
*/
static int mdl_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MDL mdl, int *nBytes) {
    int iCore = source - mdl->pmdl[0]->base.idSelf;
    int bOnNode = (iCore >= 0 && iCore < mdl->base.nCores);
    MDLserviceSend *send = &mdl->recvRequest;

    assert(source != mdlSelf(mdl));

    /* If the sender is on-node, we just wait for the send to appear and make a copy */
    if (bOnNode) {
	MDLserviceElement *qhdr;
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_MAX + iCore); /* From me */
	assert(qhdr->iServiceID == MDL_SE_MPI_SSEND);
	assert(qhdr->iCoreFrom == iCore);
	send = (MDLserviceSend *)qhdr;
	assert(send->tag == tag);
	memcpy(buf,send->buf,send->count);
	*nBytes = send->count;
	mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_MPI_SSEND);
	}

    /* Off-node: We ask the MPI thread to post a receive for us. */
    else {
	return mdl_remote_MPI_Recv(buf,count,datatype,source,tag,mdl,nBytes,MDL_SE_MPI_RECV);
	}
    return MPI_SUCCESS;
    }

/*
** Here we do a bidirectional message exchange between two threads.
*/
static int mdl_MPI_Sendrecv(
    void *sendbuf, int sendcount, MPI_Datatype sendtype,
    int dest, int sendtag, void *recvbuf, int recvcount,
    MPI_Datatype recvtype, int source, int recvtag,
    MDL mdl,  int *nReceived) {
    mdl_start_MPI_Ssend(sendbuf, sendcount, sendtype, dest, sendtag, mdl,MDL_SE_MPI_SSEND);
    mdl_MPI_Recv(recvbuf,recvcount,recvtype,source,recvtag,mdl,nReceived);
    mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
    return MPI_SUCCESS;
    }

static int mdl_MPI_Barrier(MDL mdl) {
    MDLserviceElement *qhdr;
    int i, id;
    int nBytes;

    mdlCompleteAllWork(mdl);
    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    if (mdl->base.iCore) {
	mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,0,&mdl->inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	for(i=1; i<mdl->base.nCores; ++i) {
	    qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	    assert(qhdr->iServiceID == MDL_SE_BARRIER_REQUEST);
	    }

	if (mdl->base.idSelf==0) {
	    for(i=1; i<mdl->base.nProcs; ++i) {
		mdl_MPI_Recv(&id, sizeof(id),
		    MPI_BYTE, MPI_ANY_SOURCE, MDL_TAG_BARRIER, mdl, &nBytes);
		}
	    for(i=1; i<mdl->base.nProcs; ++i) {
		id = i;
		mdl_MPI_Send(&id, sizeof(id),
		    MPI_BYTE, mdlProcToThread(mdl,i), MDL_TAG_BARRIER, mdl);
		}
	    }
	else {
	    id = mdl->base.iProc;
	    mdl_MPI_Ssend(&id,sizeof(id),MPI_BYTE,0,MDL_TAG_BARRIER,mdl);
	    mdl_MPI_Recv(&id, sizeof(id),
		MPI_BYTE, 0, MDL_TAG_BARRIER, mdl, &nBytes);
	    }
#if 0 /* Hard barrier is not required here */
	mdlSendToMPI(mdl,&mdl->inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
#endif
	for(i=1; i<mdl->base.nCores; ++i) {
	    mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,i,&mdl->pmdl[i]->inMessage,MDL_SE_BARRIER_REPLY);
	    }
	}
    mdl->wqAccepting = 0;
    mdlTimeAddSynchronizing(mdl);
    return MPI_SUCCESS;
    }

void mdlInitCommon(MDL mdl0, int iMDL,int bDiag,int argc, char **argv,
    void * (*fcnMaster)(MDL),void * (*fcnChild)(MDL)) {
    MDL mdl = mdl0->pmdl[iMDL];
    int i;

    if (mdl!=mdl0) {
	mdlBaseInitialize(&mdl->base,argc,argv);
	mdl->base.iCore = iMDL;
	mdl->base.idSelf = mdl0->base.idSelf + iMDL;
	mdl->base.nThreads = mdl0->base.nThreads;
	mdl->base.nCores = mdl0->base.nCores;
	mdl->base.nProcs = mdl0->base.nProcs;
	mdl->base.iProc = mdl0->base.iProc;
	mdl->base.iProcToThread = mdl0->base.iProcToThread;
	mdl->iCoreMPI = mdl0->iCoreMPI;
	mdl->pmdl = mdl0->pmdl;
	mdl->mpi = NULL;    /* Only used by the MPI thread */
	}
    if (mdl->base.idSelf) mdl->fcnWorker = fcnChild;
    else mdl->fcnWorker = fcnMaster;
    /* We need a queue for each TAG, and a receive queue from each thread. */
    mdl->inQueue = malloc((MDL_TAG_MAX+mdl->base.nCores) * sizeof(*mdl->inQueue));
    for(i=0; i<(MDL_TAG_MAX+mdl->base.nCores); ++i) OPA_Queue_init(mdl->inQueue+i);

    /*
    ** Our work queue. We can defer work for when we are waiting, or get work
    ** from other threads if we are idle.
    */
    OPA_Queue_init(&mdl->wq);
    OPA_Queue_init(&mdl->wqDone);
    OPA_Queue_init(&mdl->wqFree);
    mdl->wqMaxSize = 0;
    mdl->wqAccepting = 0;
    mdl->wqLastHelper = 0;
    OPA_store_int(&mdl->wqCurSize,0);

#ifdef USE_CUDA
    mdl->inCudaBufSize = mdl->outCudaBufSize = 0;
    mdl->cudaCtx = CUDA_initialize(mdl->base.iCore);
#endif

    /*
    ** Set default "maximums" for structures. These are NOT hard
    ** maximums, as the structures will be realloc'd when these
    ** values are exceeded.
    */
    mdl->nMaxSrvBytes = 0;
    /*
    ** Allocate service buffers.
    */
    mdl->pszIn = NULL;
    mdl->pszOut = NULL;
    mdl->pszBuf = NULL;
    /*
     ** Allocate swapping transfer buffer. This buffer remains fixed.
     */
    mdl->pszTrans = malloc(MDL_TRANS_SIZE);
    assert(mdl->pszTrans != NULL);
    /*
    ** Allocate initial cache spaces.
    */
    mdl->nMaxCacheIds = MDL_DEFAULT_CACHEIDS;
    mdl->cache = malloc(mdl->nMaxCacheIds*sizeof(CACHE));
    assert(mdl->cache != NULL);
    OPA_Queue_init(&mdl->wqCacheFlush);
    mdl->nFlushOutBytes = 0;

    OPA_Queue_init(&mdl->coreFlushBuffers);
    for(i=0; i<9; ++i) {
	MDLflushBuffer *pFlush = malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE);
	assert(pFlush!=NULL);
	pFlush->nBufferSize = MDL_FLUSH_DATA_SIZE;
	pFlush->nBytes = 0;
	OPA_Queue_enqueue(&mdl->coreFlushBuffers, pFlush, MDLflushBuffer, hdr.hdr);
	}
    mdl->coreFlushBuffer = malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE);
    assert(mdl->coreFlushBuffer!=NULL);
    mdl->coreFlushBuffer->nBufferSize = MDL_FLUSH_DATA_SIZE;
    mdl->coreFlushBuffer->nBytes = 0;

    /*
    ** Initialize caching spaces.
    */
    mdl->cacheSize = MDL_CACHE_SIZE;
    for (i = 0; i<mdl->nMaxCacheIds; ++i) {
	mdl->cache[i].iType = MDL_NOCACHE;
	mdl->cache[i].iCID = i;
	mdl->cache[i].arc = NULL;
	}

    mdl->base.bDiag = bDiag;
    if (mdl->base.bDiag) {
	char achDiag[256], ach[256];
	const char *tmp = strrchr(argv[0], '/');
	if (!tmp) tmp = argv[0];
	else ++tmp;
	sprintf(achDiag, "%s/%s.%d", ach, tmp, mdlSelf(mdl));
	mdl->base.fpDiag = fopen(achDiag, "w");
	assert(mdl->base.fpDiag != NULL);
	}
    }

static void drainMPI(MDL mdl) {
    while(checkMPI(mdl)) {
#ifdef _MSC_VER
	SwitchToThread();
#else
	sched_yield();
#endif
	}
    }

static void *mdlWorkerThread(void *vmdl) {
    MDL mdl = vmdl;
    void *result;
#ifdef USE_ITT
    char szName[20];
    sprintf(szName,"ID %d", mdl->base.iCore);
    __itt_thread_set_name(szName);
#endif
    result = (*mdl->fcnWorker)(mdl);
    if (mdl->base.iCore != mdl->iCoreMPI) {
	mdlSendToMPI(mdl,&mdl->inMessage,MDL_SE_STOP);
	mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
	}
    else {
	drainMPI(mdl);
	}
    return result;
    }

static void cleanupMDL(MDL mdl) {
    int i;
    for (i = 0; i<mdl->nMaxCacheIds; ++i) {
	if (mdl->cache[i].arc) arcFinish(mdl->cache[i].arc);
	}
    free(mdl->pszIn);
    free(mdl->pszOut);
    free(mdl->pszBuf);
    free(mdl->inQueue);
    free(mdl->pszTrans);
    free(mdl->cache);
    mdlBaseFinish(&mdl->base);
    }

#if defined(SIGRTMAX) && defined(HAVE_MALLOC_STATS)
static void SIGRTMAX0_handler(int signo) {
    malloc_stats();
    }
#endif


void mdlLaunch(int argc,char **argv,void * (*fcnMaster)(MDL),void * (*fcnChild)(MDL)) {
    MDL mdl;
    int i,n,bDiag,bThreads,bDedicated,thread_support,rc,flag,*piTagUB;
    char *p, ach[256];
    mdlContextMPI *mpi;

#if defined(SIGRTMAX) && defined(HAVE_MALLOC_STATS)
    signal(SIGRTMAX-1,SIGRTMAX0_handler);
#endif

#ifdef USE_ITT
    __itt_domain* domain = __itt_domain_create("MyTraces.MyDomain");
    __itt_string_handle* shMyTask = __itt_string_handle_create("MDL Startup");
    __itt_task_begin(domain, __itt_null, __itt_null, shMyTask);
#endif
    mpi = malloc(sizeof(mdlContextMPI));
    assert(mpi!=NULL);
    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    mdlBaseInitialize(&mdl->base,argc,argv);
    mdl->mpi = mpi;

#ifdef _SC_NPROCESSORS_CONF /* from unistd.h */
    mdl->base.nCores = sysconf(_SC_NPROCESSORS_CONF);
#endif

    /*
    ** Do some low level argument parsing for number of threads, and
    ** diagnostic flag!
    */
    for (argc = 0; argv[argc]; argc++);
    bDiag = 0;
    bThreads = 0;
    bDedicated = -1;
    i = 1;
    while (argv[i]) {
	if (!strcmp(argv[i], "-sz") && !bThreads) {
	    ++i;
	    mdl->base.nCores = atoi(argv[i]);
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

    if (!bThreads) {
	if ( (p=getenv("SLURM_CPUS_PER_TASK")) != NULL ) mdl->base.nCores = atoi(p);
	else if ( (p=getenv("OMP_NUM_THREADS")) != NULL ) mdl->base.nCores = atoi(p);
	}
    assert(mdl->base.nCores>0);

    /* MPI Initialization */
#ifdef USE_ITT
    __itt_string_handle* shMPITask = __itt_string_handle_create("MPI");
    __itt_task_begin(domain, __itt_null, __itt_null, shMPITask);
#endif
    rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,&thread_support);
    if (rc!=MPI_SUCCESS) {
	MPI_Error_string(rc, ach, &i);
	perror(ach);
	abort();
	}
#ifdef MDL_FFTW
    if (mdlCores(mdl)>1) FFTW3(init_threads)();
    FFTW3(mpi_init)();
    if (mdlCores(mdl)>1) FFTW3(plan_with_nthreads)(mdlCores(mdl));
#endif

#ifdef USE_ITT
    __itt_task_end(domain);
#endif
    mpi->commMDL = MPI_COMM_WORLD;
    MPI_Comm_size(mpi->commMDL, &mdl->base.nProcs);
    MPI_Comm_rank(mpi->commMDL, &mdl->base.iProc);

    /* Dedicate one of the threads for MPI, unless it would be senseless to do so */
    if (bDedicated == -1) {
	if (mdl->base.nProcs>0 && mdl->base.nCores>3) bDedicated = 1;
	else bDedicated = 0;
	}
    if (bDedicated == 1) {
	/*if (mdl->base.nCores==1 || mdl->base.nProcs==1) bDedicated=0;
	  else*/ --mdl->base.nCores;
	}

    /* Construct the thread/processor map */
    mdl->base.iProcToThread = malloc((mdl->base.nProcs + 1) * sizeof(int));
    assert(mdl->base.iProcToThread != NULL);
    mdl->base.iProcToThread[0] = 0;
    MPI_Allgather(&mdl->base.nCores, 1, MPI_INT, mdl->base.iProcToThread + 1, 1, MPI_INT, mpi->commMDL);
    for (i = 1; i < mdl->base.nProcs; ++i) mdl->base.iProcToThread[i + 1] += mdl->base.iProcToThread[i];
    mdl->base.nThreads = mdl->base.iProcToThread[mdl->base.nProcs];
    mdl->base.idSelf = mdl->base.iProcToThread[mdl->base.iProc];

    mdl->iCoreMPI = bDedicated ? -1 : 0;
    mdl->base.iCore = mdl->iCoreMPI;

    /* We have an extra MDL context in case we want a dedicated MPI thread */
    mdl->pmdl = malloc((mdl->base.nCores+1) * sizeof(struct mdlContext *));
    assert(mdl->pmdl!=NULL);
    mdl->pmdl++;
    mdl->pmdl[mdl->iCoreMPI] = mdl;
    mdl->threadid = malloc(mdl->base.nCores * sizeof(pthread_t));

    /* Allocate the other MDL structures for any threads. */
    for (i = mdl->iCoreMPI+1; i < mdl->base.nCores; ++i) {
	mdl->pmdl[i] = malloc(sizeof(struct mdlContext));
	assert(mdl->pmdl[i] != NULL);
	}
    for (i = mdl->iCoreMPI; i < mdl->base.nCores; ++i)
	mdlInitCommon(mdl, i, bDiag, argc, argv, fcnMaster, fcnChild);


    OPA_Queue_init(&mpi->queueMPI);
    OPA_Queue_header_init(&mdl->inMessage.hdr);
    OPA_Queue_header_init(&mdl->sendRequest.svc.hdr);
    OPA_Queue_header_init(&mdl->recvRequest.svc.hdr);
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
    mpi->pReqRcv = malloc(sizeof(MDLcacheReplyData) + mpi->iCacheBufSize);
    assert(mpi->pReqRcv != NULL);
    mpi->pReqRcv->mpiRequest = MPI_REQUEST_NULL;

    /* Start with a sensible number of cache buffers */
    mpi->freeCacheReplies = NULL;
    mpi->busyCacheReplies = NULL;
    mpi->busyCacheRepliesTail = &mpi->busyCacheReplies;
    n = mdl->base.nProcs;
    if (n > 256) n = 256;
    else if (n < 32) n = 64;
    for (i = 0; i<n; ++i) {
	MDLcacheReplyData *pdata = malloc(sizeof(MDLcacheReplyData) + mpi->iReplyBufSize);
	assert( pdata != NULL );
	pdata->mpiRequest = MPI_REQUEST_NULL;
	pdata->next = mpi->freeCacheReplies;
	mpi->freeCacheReplies = pdata;
	}

    /* For flushing to remote processes/nodes */
    mpi->flushBuffersByRank = malloc(sizeof(MDLflushBuffer *) * mdl->base.nProcs);
    assert(mpi->flushBuffersByRank);
    for(i=0; i<mdlProcs(mdl); ++i) mpi->flushBuffersByRank[i] = NULL;

    /* For flushing to threads on this processor */
    OPA_Queue_init(&mpi->localFlushBuffers);
    mpi->flushBuffersByCore = malloc(sizeof(MDLflushBuffer *) * mdlCores(mdl));
    assert(mpi->flushBuffersByCore);
    for(i=0; i<mdlCores(mdl); ++i) mpi->flushBuffersByCore[i] = NULL;
    for(i=0; i<2*mdlCores(mdl); ++i) {
	MDLflushBuffer *pBuffer = malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE);
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
	MDLflushBuffer *pBuffer = malloc(sizeof(MDLflushBuffer) + MDL_FLUSH_DATA_SIZE);
	assert(pBuffer!=NULL);
	pBuffer->nBufferSize = MDL_FLUSH_DATA_SIZE;
	pBuffer->nBytes = 0;
	flush_insert_after(mpi->flushHeadFree.hdr.mpi.prev,pBuffer);
	}

    /* Some bookeeping for the send/recv - 1 of each per thread */
    mpi->nSendRecvReq = 0;
    mpi->pSendRecvReq = malloc(mdl->base.nCores*2*sizeof(MPI_Request));
    assert(mpi->pSendRecvReq!=NULL);
    mpi->pSendRecvBuf = malloc(mdl->base.nCores*2*sizeof(MDLserviceSend *));
    assert(mpi->pSendRecvBuf!=NULL);
    mpi->pThreadCacheReq = malloc(mdl->base.nCores*sizeof(MDLserviceCacheReq *));
    assert(mpi->pThreadCacheReq!=NULL);
    for (i = 0; i < mdl->base.nCores; ++i) mpi->pThreadCacheReq[i] = NULL;

    /* Ring buffer of requests */
    mpi->iRequestTarget = 0;
    mpi->nRequestTargets = 2 * log2(1.0 * mdl->base.nThreads) + mdl->base.nCores;
    mpi->pRequestTargets = malloc(mpi->nRequestTargets * sizeof(*mpi->pRequestTargets));
    assert(mpi->pRequestTargets!=NULL);
    for(i=0; i<mpi->nRequestTargets; ++i) mpi->pRequestTargets[i] = -1;

    /* Make sure that MPI supports enough tags */
    rc = MPI_Comm_get_attr(mpi->commMDL,MPI_TAG_UB,&piTagUB,&flag);
    if (rc==MPI_SUCCESS && flag) {
	assert(mdl->base.nCores*MDL_TAG_THREAD_OFFSET < *piTagUB);
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
    if (mdl->base.nCores > 1 || bDedicated) {
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_barrier_init(&mdl->pmdl[0]->barrier,NULL,mdlCores(mdl)+(bDedicated?1:0));
	for (i = mdl->iCoreMPI+1; i < mdl->base.nCores; ++i) {
	    pthread_create(&mdl->threadid[i], &attr,
		mdlWorkerThread,
		mdl->pmdl[i]);
	    ++mpi->nActiveCores;
	    }
	pthread_attr_destroy(&attr);
	}
#ifdef USE_ITT
    __itt_task_end(domain);
    __itt_task_end(domain);
#endif
    if (!bDedicated) {
	if (mdl->base.idSelf) (*fcnChild)(mdl);
	else fcnMaster(mdl);
	}
    drainMPI(mdl);
    pthread_barrier_destroy(&mdl->pmdl[0]->barrier);
    mdlFinish(mdl);
    }

void mdlFinish(MDL mdl) {
    mdlContextMPI *mpi = mdl->mpi;
    MDLcacheReplyData *pdata, *pnext;
    int i;

    for (i = mdl->iCoreMPI+1; i < mdl->base.nCores; ++i) {
	pthread_join(mdl->threadid[i],0);
	cleanupMDL(mdl->pmdl[i]);
	free(mdl->pmdl[i]);
	}
    free(mdl->threadid);
    free(mdl->pmdl-1);

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
    MPI_Barrier(mdl->mpi->commMDL);
    MPI_Finalize();

#ifdef USE_CUDA
    CUDA_finish(mdl->cudaCtx);
#endif

    cleanupMDL(mdl);

    /*
     ** Close Diagnostic file.
     */
    if (mdl->base.bDiag) {
	fclose(mdl->base.fpDiag);
	}

    /*
     ** Deallocate storage.
     */
    for(pdata = mpi->freeCacheReplies; pdata != NULL; pdata=pnext) {
	pnext = pdata->next;
	free(pdata);
	}
    mpi->freeCacheReplies = NULL;

    free(mdl->mpi->pReqRcv);
    free(mdl->base.iProcToThread);
    free(mpi->pSendRecvReq);
    free(mpi->pSendRecvBuf);
    free(mpi->pThreadCacheReq);
    free(mpi->pRequestTargets);
    free(mdl->mpi);
    free(mdl);
    }

/*
** This function will transfer a block of data using a pack function.
** The corresponding node must call mdlRecv.
*/

#define SEND_BUFFER_SIZE (1*1024*1024)

void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx) {
    size_t nBuff;
    char *vOut;

    vOut = malloc(SEND_BUFFER_SIZE);
    mdlassert(mdl,vOut!=NULL);

    do {
	nBuff = (*pack)(ctx,&id,SEND_BUFFER_SIZE,vOut);
	mdl_MPI_Ssend(vOut,nBuff,MPI_BYTE,id,MDL_TAG_SEND,mdl);
	}
    while ( nBuff != 0 );

    free(vOut);
    }

void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx) {
    void *vIn;
    size_t nUnpack;
    int nBytes;
    int inid;

    if ( id < 0 ) id = MPI_ANY_SOURCE;

    vIn = malloc(SEND_BUFFER_SIZE);
    mdlassert(mdl,vIn!=NULL);

    do {
	mdl_MPI_Recv(vIn,SEND_BUFFER_SIZE,MPI_BYTE,id,MDL_TAG_SEND,mdl,&nBytes);
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
int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes,
	    size_t *pnSndBytes,size_t *pnRcvBytes) {
    size_t nInBytes,nOutBufBytes;
    int nInMax,nOutMax,nBytes;
    char *pszBuf = vBuf;
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
	&swo, sizeof(swo), MPI_BYTE, id, MDL_TAG_SWAPINIT,
	mdl, &nBytes);
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
	memcpy(mdl->pszTrans,pszOut,nOutMax);
	mdl_MPI_Sendrecv(mdl->pszTrans,nOutMax, MPI_BYTE, id, MDL_TAG_SWAP,
	    pszIn,nInMax, MPI_BYTE, id, MDL_TAG_SWAP,
	    mdl, &nBytes);
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
	mdl_MPI_Ssend(pszOut,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP,mdl);
	pszOut = &pszOut[nOutMax];
	nOutBytes -= nOutMax;
	nOutBufBytes -= nOutMax;
	*pnSndBytes += nOutMax;
	}
    while (nInBytes && nBufBytes) {
	nInMax = size_t_to_int((nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE);
	nInMax = size_t_to_int((nInMax < nBufBytes)?nInMax:nBufBytes);
	mdl_MPI_Recv(pszIn,nInMax,MPI_BYTE,id,MDL_TAG_SWAP,mdl,&nBytes);
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

void mdlCommitServices(MDL mdl) {
    int nMaxBytes;
    nMaxBytes = (mdl->base.nMaxInBytes > mdl->base.nMaxOutBytes) ? mdl->base.nMaxInBytes : mdl->base.nMaxOutBytes;
    if (nMaxBytes > mdl->nMaxSrvBytes) {
        mdl->pszIn = realloc(mdl->pszIn, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement));
        assert(mdl->pszIn != NULL);
	mdl->pszOut = realloc(mdl->pszOut, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement));
        assert(mdl->pszOut != NULL);
	mdl->pszBuf = realloc(mdl->pszBuf, nMaxBytes + sizeof(SRVHEAD) + sizeof(MDLserviceElement));
        assert(mdl->pszBuf != NULL);
        mdl->nMaxSrvBytes = nMaxBytes;
        }
    /* We need a thread barrier here because we share these buffers */
    mdlThreadBarrier(mdl);
    }

void mdlAddService(MDL mdl,int sid,void *p1,
		   void (*fcnService)(void *,void *,int,void *,int *),
		   int nInBytes,int nOutBytes) {
    mdlBaseAddService(&mdl->base, sid, p1, fcnService, nInBytes, nOutBytes);
    }


int mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes) {
    char *pszIn = vin;
    SRVHEAD *ph = (SRVHEAD *)(mdl->pszBuf);
    char *pszOut = (char *)(ph + 1);
    ph->idFrom = mdl->base.idSelf;
    ph->sid = sid;
    if (!pszIn) ph->nInBytes = 0;
    else ph->nInBytes = nInBytes;
    ph->nOutBytes = 0;
    if (nInBytes > 0) {
	assert(pszIn != NULL);
	memcpy(pszOut, pszIn, nInBytes);
	}
    mdl_start_MPI_Ssend(ph, nInBytes + (int)sizeof(SRVHEAD), MPI_BYTE, id, MDL_TAG_REQ, mdl, MDL_SE_SEND_REQUEST);
    mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
    return ph->replyTag;
    }


void mdlGetReply(MDL mdl,int rID,void *vout,int *pnOutBytes) {
    char *pszOut = vout;
    SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
    char *pszIn = &mdl->pszBuf[sizeof(SRVHEAD)];
    int nBytes,id;
    id = rID;
    mdl_remote_MPI_Recv(mdl->pszBuf, mdl->nMaxSrvBytes + (int)sizeof(SRVHEAD), MPI_BYTE,
	id, MDL_TAG_RPL, mdl, &nBytes, MDL_SE_RECV_REPLY);
    assert(nBytes == ph->nOutBytes + sizeof(SRVHEAD));
    if (ph->nOutBytes > 0 && pszOut != NULL)
	memcpy(pszOut, pszIn, ph->nOutBytes);
    if (pnOutBytes) *pnOutBytes = ph->nOutBytes;
    }

void mdlHandler(MDL mdl) {
    MDLserviceElement *qhi = (MDLserviceElement *)(mdl->pszIn);
    MDLserviceElement *qho = (MDLserviceElement *)(mdl->pszOut);
    SRVHEAD *phi = (SRVHEAD *)(qhi + 1);
    SRVHEAD *pho = (SRVHEAD *)(qho + 1);
    char *pszIn = (char *)(phi + 1);
    char *pszOut = (char *)(pho + 1);
    int sid,id,tag,nOutBytes,nBytes;

    do {
	/* We ALWAYS use MPI to send requests. */
	mdl_remote_MPI_Recv(phi, mdl->nMaxSrvBytes + sizeof(SRVHEAD),
	    MPI_BYTE, MPI_ANY_SOURCE, MDL_TAG_REQ, mdl, &nBytes, MDL_SE_RECV_REQUEST);
	assert(nBytes == phi->nInBytes + sizeof(SRVHEAD));
	id = phi->idFrom;
	sid = phi->sid;
	assert(sid < mdl->base.nMaxServices);
	if (phi->nInBytes > mdl->base.psrv[sid].nInBytes) {
	    printf("ERROR: pid=%d, sid=%d, nInBytes=%d, sid.nInBytes=%d\n",
		mdlSelf(mdl), sid, phi->nInBytes, mdl->base.psrv[sid].nInBytes);
	    }
	assert(phi->nInBytes <= mdl->base.psrv[sid].nInBytes);
	nOutBytes = 0;
	assert(mdl->base.psrv[sid].fcnService != NULL);
	(*mdl->base.psrv[sid].fcnService)(mdl->base.psrv[sid].p1, pszIn, phi->nInBytes,
	    pszOut, &nOutBytes);
	assert(nOutBytes <= mdl->base.psrv[sid].nOutBytes);
	pho->idFrom = mdl->base.idSelf;
	pho->replyTag = phi->replyTag;
	pho->sid = sid;
	pho->nInBytes = phi->nInBytes;
	pho->nOutBytes = nOutBytes;
	tag = phi->replyTag;
	mdl_start_MPI_Ssend(pho, nOutBytes + sizeof(SRVHEAD), MPI_BYTE, id, tag, mdl,MDL_SE_SEND_REPLY);
	mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
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
void *mdlSetArray(MDL mdl,size_t nmemb,size_t size,void *vdata) {
    char *data = vdata;
    mdl->nMessageData = nmemb * size;
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl)==0) {
	int i;
	for(i=0; i<mdlCores(mdl); ++i) {
	    mdl->pmdl[i]->pvMessageData = data;
	    data += mdl->pmdl[i]->nMessageData;
	    }
	}
    mdlThreadBarrier(mdl);
    return mdl->pvMessageData;
    }

/*
** New collective form. All cores must call this at the same time.
** Normally "size" is identical on all cores, while nmemb may differ
** but this is not strictly a requirement.
*/
void *mdlMallocArray(MDL mdl,size_t nmemb,size_t size,size_t minSize) {
    char *data;
    size_t iSize;
    mdl->nMessageData = nmemb * size;
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl)==0) {
	int i;
	iSize = 0;
	for(i=0; i<mdlCores(mdl); ++i) iSize += mdl->pmdl[i]->nMessageData;
	if (iSize < minSize) iSize = minSize;
	data = mdlMalloc(mdl,iSize);
	for(i=0; i<mdlCores(mdl); ++i) {
	    mdl->pmdl[i]->pvMessageData = data;
	    data += mdl->pmdl[i]->nMessageData;
	    }
	}
    mdlThreadBarrier(mdl);
    data = mdl->pvMessageData;
    iSize = nmemb * size;
    if (iSize > 4096) iSize -= 4096; /* Cheesy page hack */
    memset(data,0,iSize); /* First touch */
    mdlThreadBarrier(mdl);
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
    char *pData = vData;
    return pData + (size_t)i*(size_t)iDataSize;
    }

void mdlSetCacheSize(MDL mdl,int cacheSize) {
    mdl->cacheSize = cacheSize;
    }

void mdlCacheCheck(MDL mdl) {
    if (mdl->base.iCore == mdl->iCoreMPI) checkMPI(mdl);
    combine_all_incoming(mdl);
    }

int mdlCacheStatus(MDL mdl,int cid) {
    assert(cid >= 0);
    if (cid >= mdl->nMaxCacheIds) return MDL_NOCACHE;
    CACHE *c = &mdl->cache[cid];
    return c->iType;
    }

/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(
    MDL mdl,int cid,
    void * (*getElt)(void *pData,int i,int iDataSize),
    void *pData,int iDataSize,int nData,
    void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *)) {

    CACHE *c;
    int i,nMaxCacheIds;
    cacheOpenClose coc;

    /*
     ** We cannot reallocate this structure because there may be other threads accessing it.
     */
    assert(cid >= 0 && cid <mdl->nMaxCacheIds);
    if (cid >= mdl->nMaxCacheIds) {
	abort();
	/*
	 ** reallocate cache spaces, adding space for 2 new cache spaces
	 ** including the one just defined.
	 */
	nMaxCacheIds = cid + 3;
	mdl->cache = realloc(mdl->cache,nMaxCacheIds*sizeof(CACHE));
	assert(mdl->cache != NULL);
	/*
	 ** Initialize the new cache slots.
	 */
	for (i=mdl->nMaxCacheIds;i<nMaxCacheIds;++i) {
	    mdl->cache[i].iType = MDL_NOCACHE;
	    mdl->cache[i].iCID = i;
	    mdl->cache[i].arc = NULL;
	    }
	mdl->nMaxCacheIds = nMaxCacheIds;
	}
    c = &mdl->cache[cid];
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
    c->arc = arcReinitialize(c->arc,mdl->cacheSize/c->iLineSize,c->iLineSize,c);
    c->pOneLine = malloc(c->iLineSize);

    /*
     ** Set up the request message as much as possible!
     */
    OPA_Queue_header_init(&c->cacheRequest.svc.hdr);
    c->cacheRequest.svc.iServiceID = MDL_SE_CACHE_REQUEST;
    c->cacheRequest.svc.iCoreFrom = mdl->base.iCore;
    c->cacheRequest.caReq.cid = cid;
    c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
    c->cacheRequest.caReq.idFrom = mdl->base.idSelf;
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
    mdl_MPI_Barrier(mdl);

    /* We might need to resize the cache buffer */
    coc.iDataSize = c->iLineSize;
    mdlSendToMPI(mdl,&coc,MDL_SE_CACHE_OPEN);
    mdlWaitThreadQueue(mdl,0);

    return(c);
    }

/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData) {
    CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData,NULL,NULL,NULL);
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
    CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData,ctx,init,combine);
    }

void mdlCacheBarrier(MDL mdl,int cid) {
    mdl_MPI_Barrier(mdl);
    }

MDLflushBuffer *flush_core_buffer(MDL mdl) {
    MDLflushBuffer *pBuffer = mdl->coreFlushBuffer;
    if (pBuffer->nBytes) {
	mdlSendToMPI(mdl,&pBuffer->hdr.svc,MDL_SE_CACHE_FLUSH);
	while (OPA_Queue_is_empty(&mdl->coreFlushBuffers)) {
	    combine_all_incoming(mdl);
	    if (mdl->base.iCore == mdl->iCoreMPI) checkMPI(mdl);
	    }
	OPA_Queue_dequeue(&mdl->coreFlushBuffers, pBuffer, MDLflushBuffer, hdr.hdr);
	mdl->coreFlushBuffer = pBuffer;
	pBuffer->nBytes = 0;
	pBuffer->iRankTo = -1;
	}
    return pBuffer;
    }

/* This CDB must not be on any list when destage is called */
static CDB *destage(MDL mdl, ARC arc, CDB *temp) {
    assert(temp->data != NULL);
    assert(temp->data[-1] == _ARC_MAGIC_);
    if (temp->uId & _DIRTY_) {     /* if dirty, evict before free */
	CACHE *c = arc->cache;
	int iSize = sizeof(CAHEAD) + c->iLineSize;
	MDLflushBuffer *pBuffer = mdl->coreFlushBuffer;
	if (pBuffer->nBytes + iSize > pBuffer->nBufferSize) {
	    pBuffer = flush_core_buffer(mdl);
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

static void arcRemoveAll(MDL mdl,ARC arc) {
    CDB *temp;
    for(;arc->T1Length;--arc->T1Length) {
	temp = remove_from_hash(arc,lru_remove(arc->T1));
	destage(mdl,arc,temp);     /* if dirty, evict before free */
	lru_insert(temp,arc->Free);
	}
    for(;arc->T2Length;--arc->T2Length) {
	temp = remove_from_hash(arc,lru_remove(arc->T2));
	destage(mdl,arc,temp);     /* if dirty, evict before free */
	lru_insert(temp,arc->Free);
	}
    for(;arc->B1Length;--arc->B1Length) {
	temp = remove_from_hash(arc,lru_remove(arc->B1));
	assert(temp->data == NULL);
	mru_insert(temp,arc->Free);
	}
    for(;arc->B2Length;--arc->B2Length) {
	temp = remove_from_hash(arc,lru_remove(arc->B2));
	assert(temp->data == NULL);
	mru_insert(temp,arc->Free);
	}
    assert(arc->T1->hdr.links.next == arc->T1);
    assert(arc->T1->hdr.links.prev == arc->T1);
    assert(arc->B1->hdr.links.next == arc->B1);
    assert(arc->B1->hdr.links.prev == arc->B1);
    assert(arc->T2->hdr.links.next == arc->T2);
    assert(arc->T2->hdr.links.prev == arc->T2);
    assert(arc->B2->hdr.links.next == arc->B2);
    assert(arc->B2->hdr.links.prev == arc->B2);
    }

void mdlFlushCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];

    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    arcRemoveAll(mdl,c->arc);
    flush_core_buffer(mdl);
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl)==0) { /* This flushes all buffered data, not just our thread */
	mdlSendToMPI(mdl,&mdl->inMessage,MDL_SE_CACHE_FLUSH_OUT);
	mdlWaitThreadQueue(mdl,0);
	}
    /* We must wait for all threads to finish with this cache */
    mdl_MPI_Barrier(mdl);
    if (mdlCore(mdl)==0) { /* This flushes all buffered data, not just our thread */
	mdlSendToMPI(mdl,&mdl->inMessage,MDL_SE_CACHE_FLUSH_LCL);
	mdlWaitThreadQueue(mdl,0);
	}
    mdlThreadBarrier(mdl);
    mdl->wqAccepting = 0;
    mdlTimeAddSynchronizing(mdl);
    }

void mdlFinishCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    cacheOpenClose coc;

    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    mdlFlushCache(mdl,cid);

    coc.iDataSize = c->iDataSize;
    mdlSendToMPI(mdl,&coc,MDL_SE_CACHE_CLOSE);
    mdlWaitThreadQueue(mdl,0);

    /*
     ** Free up storage and finish.
     */
    free(c->pOneLine);
#ifdef FINISH_ARC
    arcFinish(c->arc);
    c->arc = NULL;
#endif
    c->iType = MDL_NOCACHE;
    mdl->wqAccepting = 0;
    mdlTimeAddSynchronizing(mdl);
    }

void mdlPrefetch(MDL mdl,int cid,int iIndex, int id) {
    }

static inline uint64_t *replace(MDL mdl,ARC arc, int iInB2) {
    CDB *temp;
    uint64_t *data;
    uint32_t max = (arc->target_T1 > 1)?(arc->target_T1+1-iInB2):1;
    if (arc->T1Length >= max) { /* T1’s size exceeds target? */
                                        /* yes: T1 is too big */
	temp = arc->T1->hdr.links.prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->hdr.links.prev;
	    if (temp == arc->T1) {           /* all pages in T1 are currently locked, try T2 in this case */
		temp = arc->T2->hdr.links.prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->hdr.links.prev;
		    if (temp == arc->T2) return(NULL); /* all pages are currently locked, give up! */
		}
		goto replace_T2;
	    }
	}
	goto replace_T1;
    } else {
	/* no: T1 is not too big */
	temp = arc->T2->hdr.links.prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->hdr.links.prev;
	    if (temp == arc->T2) {           /* all pages in T2 are currently locked, try T1 in this case */
		temp = arc->T1->hdr.links.prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->hdr.links.prev;
		    if (temp == arc->T1) return(NULL); /* all pages are currently locked, give up! */
		}
		goto replace_T1;
	    }
	}
	goto replace_T2;
    }
    assert(0);
    if (0) {  /* using a Duff's device to handle the replacement */
    replace_T1: 
        remove_from_list(temp);          /* grab LRU unlocked page from T1 */
	destage(mdl,arc,temp);     /* if dirty, evict before overwrite */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
	mru_insert(temp,arc->B1);           /* put it on B1 */
	temp->uId = (temp->uId&_IDMASK_)|_B1_;  /* need to be careful here because it could have been _P1_ */
/* 	assert( (temp->uId&_WHERE_) == _B1_ ); */
        arc->T1Length--; arc->B1Length++;          /* bookkeep */
	/*assert(arc->B1Length<=arc->nCache);*/
    }
    if (0) {
    replace_T2:
        remove_from_list(temp);          /* grab LRU unlocked page from T2 */
	destage(mdl,arc,temp);     /* if dirty, evict before overwrite */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
        mru_insert(temp,arc->B2);           /* put it on B2 */
        temp->uId |= _B2_;          /* note that fact */
/* 	assert( (temp->uId&_WHERE_) == _B2_ ); */
        arc->T2Length--; arc->B2Length++;          /* bookkeep */
    }
    /*assert(data!=NULL);*/
    return data;
}

/*
** The highest order bit (1<<31) of uId encodes the dirty bit.
** The next 3 highest order bits of uId ((1<<30), (1<<29) and (1<<28)) are reserved 
** for list location and should be zero! The maximum legal uId is then (1<<28)-1.
*/
static inline CDB *arcSetPrefetchDataByHash(MDL mdl,ARC arc,uint32_t uPage,uint32_t uId,void *data,uint32_t uHash) {
    CDB *temp;
    uint32_t tuId = uId&_IDMASK_;
    uint32_t L1Length;
    int inB2=0;

    assert(data);
    for( temp = arc->Hash[uHash]; temp; temp = temp->extra.coll ) {
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
	    inB2=1;
	doBcase:
	    /* Better would be to put this back on P1, but L1 may be full. */
	    remove_from_list(temp);                                   /* take off whichever list */
	    temp->data = replace(mdl,arc,inB2);                                /* find a place to put new page */
	    temp->uPage = uPage;                          /* bookkeep */
	    if (inB2) {
		temp->uId = _T2_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 		assert( (temp->uId&_WHERE_) == _T2_ ); */
		mru_insert(temp,arc->T2);                                     /* not seen yet, put on T1 */
		arc->T2Length++;
		}
	    else {
		temp->uId = _P1_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 		assert( (temp->uId&_WHERE_) == _P1_ ); */
		mru_insert(temp,arc->T1);                                     /* not seen yet, put on T1 */
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
		temp = remove_from_hash(arc,lru_remove(arc->B1));        /* yes: take page off B1 */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = replace(mdl,arc,0);                                /* find new place to put page */
	    } else {                                                      /* no: B1 must be empty */
		temp = remove_from_hash(arc,lru_remove(arc->T1));         /* take page off T1 */
		destage(mdl,arc,temp);     /* if dirty, evict before overwrite */
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
		    temp = remove_from_hash(arc,lru_remove(arc->B2));/* here we lose memory of what was in lru B2 */
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=1;
		} else {                                                   /* cache directory not full, easy case */
		    temp = lru_remove(arc->Free);
		}
		temp->data = replace(mdl,arc,inB2);                                /* new place for page */
	    } else {                                                      /* cache not full, easy case */
		temp = lru_remove(arc->Free);
		assert(temp->data != NULL);               /* This CDB should have an unused page associated with it */
		temp->data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
	    }
	}
	mru_insert(temp,arc->T1);                                             /* not been seen yet, but put on T1 */
	arc->T1Length++;                                                       /* bookkeep: */
	/*assert(arc->T1Length+arc->B1Length<=arc->nCache);*/
	temp->uId = _P1_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
/* 	assert( (temp->uId&_WHERE_) == _P1_ ); */
	temp->uPage = uPage;
	temp->extra.coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */
    }
    memcpy(temp->data,data,arc->cache->iLineSize); /* Copy actual cache data amount */
    return temp;
}

static inline CDB *arcSetPrefetchData(MDL mdl,ARC arc,uint32_t uPage,uint32_t uId,void *data) {
    return arcSetPrefetchDataByHash(mdl,arc,uPage,uId,data,MurmurHash2(uPage,uId&_IDMASK_)&arc->uHashMask);
    }

/*
** This releases a lock on a page by decrementing the lock counter for that page. It also checks
** that the magic number matches, and that the lock count is greater than 0. This is to prevent
** strange errors if the user tries to release the wrong pointer (could silently modify memory 
** that it should not if this check was not there).
*/
static inline void arcRelease(ARC arc,uint64_t *p) {
    if (p>arc->dataBase && p<arc->dataLast) { /* Might have been a fast, read-only grab */
	/* We will be given an element, but this needs to be turned into a cache line */
	p = arc->dataBase + (p - arc->dataBase) / (arc->uDataSize+1) * (arc->uDataSize+1) + 1;
	uint64_t t = p[-1]-1;
	assert((t^_ARC_MAGIC_) < 0x00000000ffffffff);
	p[-1] = t;
	}
}

static void queueCacheRequest(MDL mdl, int cid, int iLine, int id) {
    int iCore = id - mdl->pmdl[0]->base.idSelf;
    /* Local requests for combiner cache are handled in finishCacheRequest() */
    if (iCore < 0 || iCore >= mdl->base.nCores ) {
	CACHE *c = &mdl->cache[cid];
	c->cacheRequest.caReq.cid = cid;
	c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
	c->cacheRequest.caReq.idFrom = mdl->base.idSelf;
	c->cacheRequest.caReq.idTo = id;
	c->cacheRequest.caReq.iLine = iLine;
	c->cacheRequest.caReq.nItems = c->nLineElements;
	c->cacheRequest.pLine = c->pOneLine;
	mdlSendToMPI(mdl,&c->cacheRequest,MDL_SE_CACHE_REQUEST);
	}
    }

static void finishCacheRequest(MDL mdl, int cid, int iLine, int id, CDB *temp,int bVirtual) {
    int iCore = id - mdl->pmdl[0]->base.idSelf;
    CACHE *c = &mdl->cache[cid];
    int i,s,n;
    char *pData = (char *)temp->data;
    uint32_t iIndex = iLine << c->nLineBits;

    s = iIndex;
    n = s + c->nLineElements;

    if (bVirtual) {
	memset(temp->data,0,c->iLineSize);
	}
    /* Local requests must be from a combiner cache if we get here */
    else if (iCore >= 0 && iCore < mdl->base.nCores ) {
	MDL omdl = mdl->pmdl[iCore];
	CACHE *oc = &omdl->cache[cid];
	if ( n > oc->nData ) n = oc->nData;
	for(i=s; i<n; i++ ) {
	    memcpy(pData,(*oc->getElt)(oc->pData,i,oc->iDataSize),oc->iDataSize);
	    pData += oc->iDataSize;
	    }
	}
    else {
	assert(iLine == temp->uPage);
	mdlWaitThreadQueue(mdl,MDL_TAG_CACHECOM);
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
void *mdlAccess(MDL mdl, int cid, uint32_t uIndex, int uId, int bLock,int bModify,int bVirtual) {
    CACHE *c = &mdl->cache[cid];
    ARC arc = c->arc;
    CDB *temp;
    uint32_t L1Length;
    uint32_t uHash,rat;
    uint32_t tuId = uId&_IDMASK_;
    uint32_t uLine;
    uint32_t iInLine;
    int inB2=0;

    if (!(++c->nAccess & MDL_CHECK_MASK)) mdlCacheCheck(mdl);

    /* Short circuit the cache if this belongs to another thread (or ourselves) */
    uint32_t uCore = uId - mdl->pmdl[0]->base.idSelf;
    if (uCore < mdl->base.nCores && c->iType == MDL_ROCACHE ) {
	MDL omdl = mdl->pmdl[uCore];
	c = &omdl->cache[cid];
	return (*c->getElt)(c->pData,uIndex,c->iDataSize);
	}
    iInLine = uIndex & c->nLineMask;
    uLine = uIndex >> c->nLineBits;

    /* First check our own cache */
    uHash = (MurmurHash2(uLine,tuId)&arc->uHashMask);
    for ( temp = arc->Hash[uHash]; temp; temp = temp->extra.coll ) {
	if (temp->uPage == uLine && (temp->uId&_IDMASK_) == tuId) break;
	}
    if (temp != NULL) {                       /* found in cache directory? */
	switch (temp->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	    temp->uId = uId;     /* clears prefetch flag and sets WHERE = _T1_ and dirty bit */
/* 	    assert( (temp->uId&_WHERE_) == _T1_ ); */
	    remove_from_list(temp);                           /* take off T1 list */
	    mru_insert(temp,arc->T1);                         /* first access but recently prefetched, put on T1 */
	    goto cachehit;
	case _T1_:
	    arc->T1Length--; arc->T2Length++;
	    temp->uId |= _T2_;
/* 	    assert( (temp->uId&_WHERE_) == _T2_ ); */
	    /* fall through */
	case _T2_:
	    remove_from_list(temp);                                   /* take off whichever list */
	    mru_insert(temp,arc->T2);                                     /* seen twice recently, put on T2 */
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
	    if (!bVirtual) queueCacheRequest(mdl,cid,uLine,uId);
/* 	    assert(arc->B1->hdr.links.next != arc->B1); */
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
	    if (!bVirtual) queueCacheRequest(mdl,cid,uLine,uId);
/* 	    assert(arc->B2->hdr.links.next != arc->B2); */
/* 	    assert(arc->B2Length>0); */

	    rat = arc->B1Length/arc->B2Length;
	    if (rat < 1) rat = 1;
	    if (rat > arc->target_T1) arc->target_T1 = 0;
	    else arc->target_T1 = arc->target_T1 - rat;
	    /* adapt the target size */
	    arc->B2Length--;                                           /* bookkeep */
	    inB2=1;
	doBcase:
	    remove_from_list(temp);                                   /* take off whichever list */
	    temp->data = replace(mdl,arc,inB2);                                /* find a place to put new page */
	    temp->uId = _T2_|uId;     /* temp->ARC_where = _T2_; and set the dirty bit for this page */
/* 	    assert( (temp->uId&_WHERE_) == _T2_ ); */
	    temp->uPage = uLine;                          /* bookkeep */
	    mru_insert(temp,arc->T2);                                     /* seen twice recently, put on T2 */
	    arc->T2Length++;                 /* JS: this was not in the original code. Should it be? bookkeep */
	    /*assert(temp->data!=NULL);*/
	    finishCacheRequest(mdl,cid,uLine,uId,temp,bVirtual);
	    break;
	    }
	}

    else {                                                              /* page is not in cache directory */
	++c->nMiss;
	mdlTimeAddComputing(mdl);
	/*
	** Can initiate the data request right here, and do the rest while waiting...
	*/
	if (!bVirtual) queueCacheRequest(mdl,cid,uLine,uId);
	L1Length = arc->T1Length + arc->B1Length;
	/*assert(L1Length<=arc->nCache);*/
	if (L1Length == arc->nCache) {                                   /* B1 + T1 full? */
	    if (arc->T1Length < arc->nCache) {                                           /* Still room in T1? */
		temp = remove_from_hash(arc,lru_remove(arc->B1));        /* yes: take page off B1 */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = replace(mdl,arc,0);                                /* find new place to put page */
		}
	    else {                                                      /* no: B1 must be empty */
		temp = remove_from_hash(arc,lru_remove(arc->T1));       /* take page off T1 */
		destage(mdl,arc,temp);     /* if dirty, evict before overwrite */
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
		    temp = remove_from_hash(arc,lru_remove(arc->B2));
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=1;
		} else {                                                   /* cache directory not full, easy case */
		    temp = lru_remove(arc->Free);
		    assert(temp->data == NULL);            /* This CDB should not be associated with data */
		}
		temp->data = replace(mdl,arc,inB2);                                /* new place for page */
		/*assert(temp->data!=NULL);*/
	    } else {                                                      /* cache not full, easy case */
		temp = lru_remove(arc->Free);
		assert(temp->data != NULL);               /* This CDB should have an unused page associated with it */
		temp->data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
		}
	    }
	mru_insert(temp,arc->T1);                                             /* seen once recently, put on T1 */
	arc->T1Length++;                                                       /* bookkeep: */
	/*assert(arc->T1Length+arc->B1Length<=arc->nCache);*/
	temp->uId = uId;                  /* temp->dirty = dirty;  p->ARC_where = _T1_; as well! */
/* 	assert( (temp->uId&_WHERE_) == _T1_ ); */
	temp->uPage = uLine;
	finishCacheRequest(mdl,cid,uLine,uId,temp,bVirtual);
	temp->extra.coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */
	mdlTimeAddWaiting(mdl);
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
    const int virtual = 0; /* really fetch the element */
    return mdlAccess(mdl, cid, iIndex, id, lock, modify, virtual);
    }

/* Locks, so mdlRelease must be called eventually */
void *mdlAcquire(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 1;  /* we always lock in acquire */
    const int modify = (mdl->cache[cid].iType == MDL_COCACHE);
    const int virtual = 0; /* really fetch the element */
    return mdlAccess(mdl, cid, iIndex, id, lock, modify, virtual);
    }

/* Locks the element, but does not fetch or initialize */
void *mdlVirtualFetch(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 0; /* fetch never locks */
    const int modify = 1; /* virtual always modifies */
    const int virtual = 1; /* do not fetch the element */
    return mdlAccess(mdl, cid, iIndex, id, lock, modify, virtual);
    }

void mdlRelease(MDL mdl,int cid,void *p) {
    CACHE *c = &mdl->cache[cid];
    arcRelease(c->arc,p);
    }


double mdlNumAccess(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];

    return(c->nAccess);
    }


double mdlMissRatio(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccess;

    if (dAccess > 0.0) return(c->nMiss/dAccess);
    else return(0.0);
    }


double mdlCollRatio(MDL mdl,int cid) {
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
void mdlGridInitialize(MDL mdl,MDLGRID *pgrid,int n1,int n2,int n3,int a1) {
    MDLGRID grid;
    assert(n1>0&&n2>0&&n3>0);
    assert(n1<=a1);
    *pgrid = grid = malloc(sizeof(struct mdlGridContext)); assert(grid!=NULL);
    grid->n1 = n1;
    grid->n2 = n2;
    grid->n3 = n3;
    grid->a1 = a1;

    /* This will be shared later (see mdlGridShare) */
    grid->id = malloc(sizeof(*grid->id)*(grid->n3));    assert(grid->id!=NULL);
    grid->rs = mdlMalloc(mdl,sizeof(*grid->rs)*mdl->base.nProcs); assert(grid->rs!=NULL);
    grid->rn = mdlMalloc(mdl,sizeof(*grid->rn)*mdl->base.nProcs); assert(grid->rn!=NULL);

    /* The following need to be set to appropriate values still. */
    grid->sSlab = grid->nSlab = 0;
    grid->nLocal = 0;
    }

void mdlGridFinish(MDL mdl, MDLGRID grid) {
    if (grid->rs) free(grid->rs);
    if (grid->rn) free(grid->rn);
    if (grid->id) free(grid->id);
    free(grid);
    }

void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, uint64_t nLocal) {
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
void mdlGridShare(MDL mdl,MDLGRID grid) {
    int i, id;
    gridShare share;

    share.grid = grid;
    mdlSendToMPI(mdl,&share,MDL_SE_GRID_SHARE);
    mdlWaitThreadQueue(mdl,0);

    /* Calculate on which processor each slab can be found. */
    for(id=0; id<mdl->base.nProcs; id++ ) {
	for( i=grid->rs[id]; i<grid->rs[id]+grid->rn[id]; i++ ) grid->id[i] = id;
	}
    }

/*
** Retrieve the first and last element for the calling thread.
*/
void mdlGridCoordFirstLast(MDL mdl,MDLGRID grid,mdlGridCoord *f,mdlGridCoord *l,int bCacheAlign) {
    uint64_t nPerCore, nThisCore;
    uint64_t nLocal = (uint64_t)(grid->a1) * grid->n2 * grid->nSlab;

    /* Number on each core with multiples of "a1" elments (complete pencils). */
    /* This needs to change to be complete MDL "cache lines" at some point. */
    uint64_t nAlign;
    if (bCacheAlign) nAlign = 1 << (int)(log2(MDL_CACHE_DATA_SIZE / sizeof(float)));
    else nAlign = grid->a1;

    nPerCore = nLocal / mdlCores(mdl) + nAlign - 1;
    nPerCore -= nPerCore % nAlign;
    nThisCore = nPerCore;
    if (mdlCore(mdl) == mdlCores(mdl)-1) {
	nThisCore = nLocal - (mdlCores(mdl)-1)*nThisCore;
	}
    /* Calculate global x,y,z coordinates, and local "i" coordinate. */
    f->I = nPerCore * mdlCore(mdl);
    l->I = f->I + nThisCore;
    f->i = 0;
    l->i = nThisCore;

    f->z = f->I/(grid->a1*grid->n2);
    f->y = f->I/grid->a1 - f->z * grid->n2;
    f->x = f->I - grid->a1 * (f->z*grid->n2 + f->y);
    assert(bCacheAlign || f->x == 0); // MDL depends on this at the moment
    f->z += grid->sSlab;
    f->grid = grid;

    l->z = l->I/(grid->a1*grid->n2);
    l->y = l->I/grid->a1 - l->z * grid->n2;
    l->x = l->I - grid->a1 * (l->z*grid->n2 + l->y);
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
size_t mdlFFTlocalCount(MDL mdl,int n1,int n2,int n3,int *nz,int *sz,int *ny,int*sy) {
    fftSizes sizes;
    sizes.n1 = n1;
    sizes.n2 = n2;
    sizes.n3 = n3;
    mdlSendToMPI(mdl,&sizes,MDL_SE_FFT_SIZES);
    mdlWaitThreadQueue(mdl,0);
    if (nz!=NULL) *nz = sizes.nz;
    if (sz!=NULL) *sz = sizes.sz;
    if (ny!=NULL) *ny = sizes.ny;
    if (sy!=NULL) *sy = sizes.sy;
    return sizes.nLocal;
    }

MDLFFT mdlFFTNodeInitialize(MDL mdl,int n1,int n2,int n3,int bMeasure,FFTW3(real) *data) {
    fftPlans plans;
    MDLFFT fft = malloc(sizeof(struct mdlFFTContext));
    assert( fft != NULL);
    assert(mdlCore(mdl) == 0);
    plans.sizes.n1 = n1;
    plans.sizes.n2 = n2;
    plans.sizes.n3 = n3;
    plans.data = 0;/*Estimate is faster? data;*/
    plans.kdata = 0;/*(FFTW3(complex) *)data;*/
    mdlSendToMPI(mdl,&plans,MDL_SE_FFT_PLANS);
    mdlWaitThreadQueue(mdl,0);
    fft->fplan = plans.fplan;
    fft->iplan = plans.iplan;

    /*
    ** Dimensions of k-space and r-space grid.  Note transposed order.
    ** Note also that the "actual" dimension 1 side of the r-space array
    ** can be (and usually is) larger than "n1" because of the inplace FFT.
    */
    mdlGridInitialize(mdl,&fft->rgrid,n1,n2,n3,2*(n1/2+1));
    mdlGridInitialize(mdl,&fft->kgrid,n1/2+1,n3,n2,n1/2+1);
    mdlGridSetLocal(mdl,fft->rgrid,plans.sizes.sz,plans.sizes.nz,plans.sizes.nLocal);
    mdlGridSetLocal(mdl,fft->kgrid,plans.sizes.sy,plans.sizes.ny,plans.sizes.nLocal/2);
    mdlGridShare(mdl,fft->rgrid);
    mdlGridShare(mdl,fft->kgrid);
    return fft;
    }

MDLFFT mdlFFTInitialize(MDL mdl,int n1,int n2,int n3,int bMeasure,FFTW3(real) *data) {
    MDLFFT fft;
    if (mdlCore(mdl) == 0) {
	mdl->pvMessageData = mdlFFTNodeInitialize(mdl,n1,n2,n3,bMeasure,data);
	}
    mdlThreadBarrier(mdl); /* To synchronize pvMessageData */
    fft = mdl->pmdl[0]->pvMessageData;
    mdlThreadBarrier(mdl); /* In case we reuse pvMessageData */
    return fft;
    }

void mdlFFTNodeFinish( MDL mdl, MDLFFT fft ) {
    FFTW3(destroy_plan)(fft->fplan);
    FFTW3(destroy_plan)(fft->iplan);
    mdlGridFinish(mdl,fft->kgrid);
    mdlGridFinish(mdl,fft->rgrid);
    free(fft);
    }
void mdlFFTFinish( MDL mdl, MDLFFT fft ) {
    mdlThreadBarrier(mdl);
    if (mdlCore(mdl) == 0) {
	mdlFFTNodeFinish(mdl,fft);
	}
    }

FFTW3(real) *mdlFFTMalloc( MDL mdl, MDLFFT fft ) {
    return mdlGridMalloc(mdl,fft->rgrid,sizeof(FFTW3(real)));
    }

void mdlFFTFree( MDL mdl, MDLFFT fft, void *p ) {
    mdlGridFree(mdl,fft->rgrid,p);
    }

void mdlFFT( MDL mdl, MDLFFT fft, FFTW3(real) *data ) {
    fftTrans trans;
    mdlThreadBarrier(mdl);
    if (mdl->base.iCore == mdl->iCoreMPI) {
	FFTW3(execute_dft_r2c)(fft->fplan,data,(FFTW3(complex) *)(data));
	}
    else if (mdlCore(mdl) == 0) {
	trans.fft = fft;
	trans.data = data;
	trans.kdata = (FFTW3(complex) *)data;
	mdlSendToMPI(mdl,&trans,MDL_SE_FFT_DFT_R2C);
	}
    pthread_barrier_wait(&mdl->pmdl[0]->barrier);
    mdlThreadBarrier(mdl);
    }
void mdlIFFT( MDL mdl, MDLFFT fft, FFTW3(complex) *kdata ) {
    fftTrans trans;
    mdlThreadBarrier(mdl);
    if (mdl->base.iCore == mdl->iCoreMPI) {
	FFTW3(execute_dft_c2r)(fft->iplan,kdata,(FFTW3(real) *)(kdata));
	}
    else if (mdlCore(mdl) == 0) {
	trans.fft = fft;
	trans.kdata = kdata;
	trans.data = (FFTW3(real) *)kdata;
	mdlSendToMPI(mdl,&trans,MDL_SE_FFT_DFT_C2R);
	}
    pthread_barrier_wait(&mdl->pmdl[0]->barrier);
    mdlThreadBarrier(mdl);
    }
#endif

#ifdef USE_CUDA
void mdlSetCudaBufferSize(MDL mdl,int inBufSize, int outBufSize) {
    if (mdl->inCudaBufSize < inBufSize) mdl->inCudaBufSize = inBufSize;
    if (mdl->outCudaBufSize < outBufSize) mdl->outCudaBufSize = outBufSize;
    }
#endif

void mdlSetWorkQueueSize(MDL mdl,int wqMaxSize,int cudaSize) {
    MDLwqNode *work;
    int i;

#ifdef USE_CUDA
    CUDA_SetQueueSize(mdl->cudaCtx,cudaSize,mdl->inCudaBufSize,mdl->outCudaBufSize);
#endif
    while (wqMaxSize > mdl->wqMaxSize) {
	for(i=0; i<mdl->base.nCores; ++i) {
	    work = malloc(sizeof(MDLwqNode));
	    OPA_Queue_header_init(&work->q.hdr);
	    work->iCoreOwner = mdl->base.iCore;
	    work->ctx = NULL;
	    work->doFcn = NULL;
	    work->doneFcn = NULL;
	    OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	    }
	++mdl->wqMaxSize;
	}
    while (wqMaxSize < mdl->wqMaxSize) {
	for(i=0; i<mdl->base.nCores; ++i) {
	    if (!OPA_Queue_is_empty(&mdl->wqFree)) 
		OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, q.hdr);
	    else assert(0);
	    free(work);
	    }
	--mdl->wqMaxSize;
	} 
    }

void mdlAddWork(MDL mdl, void *ctx,
    int (*initWork)(void *ctx,void *vwork),
    int (*checkWork)(void *ctx,void *vwork),
    mdlWorkFunction doWork,
    mdlWorkFunction doneWork) {
    MDL Qmdl = NULL;
    MDLwqNode *work;
    int i;

    /* We prefer to let CUDA do the work */
#ifdef USE_CUDA
    if (CUDA_queue(mdl->cudaCtx,ctx,initWork,checkWork,doneWork)) return;
#endif
    /* Obviously, we can only queue work if we have a free queue element */
    if (!OPA_Queue_is_empty(&mdl->wqFree)) {
	/* We have some room, so save work for later */
	if (OPA_load_int(&mdl->wqCurSize) < mdl->wqMaxSize) Qmdl = mdl;
	/* See if anyone else is accepting work */
	else {
	    i = mdl->wqLastHelper;
	    do {
		MDL rMDL = mdl->pmdl[i];
		if (rMDL->wqAccepting && OPA_load_int(&rMDL->wqCurSize)<rMDL->wqMaxSize) {
		    Qmdl = rMDL;
		    mdl->wqLastHelper = i;
		    break;
		    }
		if (++i == mdl->base.nCores) i = 0;
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
