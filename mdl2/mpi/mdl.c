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
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef __linux__
#include <sys/resource.h>
#endif
#include "mpi.h"
#include "mdl.h"

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_DEFAULT_BYTES		80000
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE		5000000

#define MDL_TAG_BARRIER        	1
#define MDL_TAG_SWAPINIT 	2
#define MDL_TAG_SWAP		3
#define MDL_TAG_REQ	   	4
#define MDL_TAG_RPL		5 /* This is treated specially */
#define MDL_TAG_SEND            6
#define MDL_TAG_CACHECOM	7

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

#define MDL_SE_MPI_SEND        11
#define MDL_SE_MPI_SSEND       12
#define MDL_SE_MPI_RECV        13


#ifdef USE_ARC
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


static inline CDB *lru_remove(CDB *list) {
    CDB *p = list->prev;
    p->prev->next = list;
    list->prev = p->prev;
    return(p);
}

static inline void remove_from_list(CDB *t) {
    t->prev->next = t->next;
    t->next->prev = t->prev;
}

static inline void mru_insert(CDB *p,CDB *list) {
    p->prev = list;
    p->next = list->next;
    list->next->prev = p;
    list->next = p;
}

static inline void lru_insert(CDB *p,CDB *list) {
    p->prev = list->prev;
    p->next = list;
    list->prev->next = p;
    list->prev = p;
}

static inline void remove_from_hash(ARC arc,CDB *p) {
    uint32_t uIndex = p->uIndex;
    uint32_t uId = p->uId&_IDMASK_;
    uint32_t uHash = (MurmurHash2(uIndex,uId)&arc->uHashMask);
    CDB **pt = &arc->Hash[uHash];
    while (*pt) {
	if ((*pt)->uIndex == uIndex) {
	    if (((*pt)->uId&_IDMASK_) == uId) {
		/*
		** Unlink.
		*/
		*pt = (*pt)->coll;
		return;
		}
	    }
	pt = &((*pt)->coll);
	}
    printf("Tried to remove uIndex %d uId %x\n", uIndex, uId);
    assert(0);  /* should never get here, the element should always be found in the hash */
    }

ARC arcInitialize(uint32_t nCache,uint32_t uDataSize) {
    ARC arc;
    uint32_t i;
    
    /*
    ** Allocate stuff.
    */
    arc = malloc(sizeof(struct ArcContext));
    assert(arc != NULL);
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
    arc->T1->next = arc->T1;
    arc->T1->prev = arc->T1;
    arc->T1->uId = 0xdeadbeef;
    arc->T1Length = 0;
    arc->B1 = malloc(sizeof(CDB));
    assert(arc->B1 != NULL);
    arc->B1->next = arc->B1;
    arc->B1->prev = arc->B1;
    arc->B1->uId = 0xdeadbeef;
    arc->B1Length = 0;
    arc->T2 = malloc(sizeof(CDB));
    assert(arc->T2 != NULL);
    arc->T2->next = arc->T2;
    arc->T2->prev = arc->T2;
    arc->T2->uId = 0xdeadbeef;
    arc->T2Length = 0;
    arc->B2 = malloc(sizeof(CDB));
    assert(arc->B2 != NULL);
    arc->B2->next = arc->B2;
    arc->B2->prev = arc->B2;
    arc->B2->uId = 0xdeadbeef;
    arc->B2Length = 0;
    arc->Free = malloc(sizeof(CDB));
    assert(arc->Free != NULL);
    arc->Free->next = arc->Free;
    arc->Free->prev = arc->Free;
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
	arc->cdbBase[i].coll = NULL;
	lru_insert(&arc->cdbBase[i],arc->Free);
    }
    /*
    ** Finally insert CDBs without data pages into the Free list last.
    */
    for (i=nCache;i<2*nCache;++i) {
	arc->cdbBase[i].data = 0;
	arc->cdbBase[i].coll = NULL;
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
#endif


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
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    uint32_t iDataSize;
    } cacheOpenClose;

static void mdlSendThreadMessage(MDL mdl,int iQueue,int iCore, void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    /*OPA_Queue_header_init(&qhdr->hdr); DONE ONCE AT START */
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = mdl->base.iCore;
    OPA_Queue_enqueue(mdl->pmdl[iCore]->inQueue+iQueue, qhdr, MDLserviceElement, hdr);
    }

static void mdlSendToMPI(MDL mdl,void *vhdr, uint16_t iServiceID ) {
    MDLserviceElement *qhdr = (MDLserviceElement *)vhdr;
    /*OPA_Queue_header_init(&qhdr->hdr); DONE ONCE AT START */
    qhdr->iServiceID = iServiceID;
    qhdr->iCoreFrom = mdl->base.iCore;
    OPA_Queue_enqueue(&mdl->pmdl[mdl->iCoreMPI]->mpi.queueMPI, qhdr, MDLserviceElement, hdr);
    }

#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEFLSH	5

int mdlCacheReceive(MDL mdl) {
    mdlContextMPI *mpi = &mdl->mpi;
    CACHE *c;
    CAHEAD *ph = (CAHEAD *)mpi->pszRcv;
    MDLserviceCacheReq *creq;
    char *pszRcv = &mpi->pszRcv[sizeof(CAHEAD)];
    CAHEAD *phRpl;
    char *pszRpl;
    char *t;
    int id,tag;
    int s,n,i;
    MPI_Status status;
    int ret;
    int iLineSize, iProc, iCore;
    char *pLine;

    ret = MPI_Wait(&mpi->ReqRcv, &status); /* This will never block by design */
    assert(ret == MPI_SUCCESS);
    mpi->ReqRcv = MPI_REQUEST_NULL;

    /* Well, could be any threads cache */
    iCore = ph->idTo - mdl->pmdl[0]->base.idSelf;
    assert(iCore>=0 && iCore<mdl->base.nCores);
    c = &mdl->pmdl[iCore]->cache[ph->cid];
    assert(c->iType != MDL_NOCACHE);

    switch (ph->mid) {
	/* A different process wants local cache data - we can grab it directly */
    case MDL_MID_CACHEREQ:
	iProc = mdlThreadToProc(mdl,ph->idFrom);
	/*
	 ** This is the tricky part! Here is where the real deadlock
	 ** difficulties surface. Making sure to have one buffer per
	 ** thread solves those problems here.
	 */
	pszRpl = &mpi->ppszRpl[ph->idFrom][sizeof(CAHEAD)];
	phRpl = (CAHEAD *)mpi->ppszRpl[ph->idFrom];
	phRpl->cid = ph->cid;
	phRpl->mid = MDL_MID_CACHERPL;
	phRpl->idFrom = mdl->base.idSelf;
	phRpl->idTo = ph->idFrom;
	assert(ph->iLine>=0);

	s = ph->iLine*MDL_CACHELINE_ELTS;
	n = s + MDL_CACHELINE_ELTS;
	if ( n > c->nData ) n = c->nData;
	iLineSize = (n-s) * c->iDataSize;
	for(i=s; i<n; i++ ) {
	    t = (*c->getElt)(c->pData,i,c->iDataSize);
	    memcpy(pszRpl+(i-s)*c->iDataSize,t,c->iDataSize);
	    }
	if (mpi->pmidRpl[ph->idFrom] != -1) {
	    MPI_Wait(&mpi->pReqRpl[ph->idFrom], &status);
	    }
	mpi->pmidRpl[ph->idFrom] = 0;
	/*iCore = ph->idFrom - mdlProcToThread(mdl,iProc)*/;
	tag = MDL_TAG_CACHECOM /*+ MDL_TAG_THREAD_OFFSET * iCore*/;
	MPI_Isend(phRpl,(int)sizeof(CAHEAD)+iLineSize,MPI_BYTE,
		  iProc, tag, mpi->commMDL,
		  &mpi->pReqRpl[ph->idFrom]);
	ret = 0;
	break;
    case MDL_MID_CACHEFLSH:
	assert(c->iType == MDL_COCACHE);
	s = ph->iLine*MDL_CACHELINE_ELTS;
	n = s + MDL_CACHELINE_ELTS;
	if (n > c->nData) n = c->nData;
	for(i=s;i<n;i++) {
		(*c->combine)(c->ctx,(*c->getElt)(c->pData,i,c->iDataSize),
			      &pszRcv[(i-s)*c->iDataSize]);
	    }
	ret = 0;
	break;
	/* A remote process has sent us cache data - we need to let that thread know */
    case MDL_MID_CACHERPL:
	creq = mpi->pThreadCacheReq[iCore];
	assert(creq!=NULL);
	mpi->pThreadCacheReq[iCore] = NULL;
	pLine = creq->pLine;

	/*
	 ** For now assume no prefetching!
	 ** This means that this WILL be the reply to this Aquire
	 ** request.
	 */
	assert(pLine != NULL);
	iLineSize = c->iLineSize;
	for (i=0;i<iLineSize;++i) pLine[i] = pszRcv[i];
	if (c->iType == MDL_COCACHE && c->init) {
	    /*
	     ** Call the initializer function for all elements in
	     ** the cache line.
	     */
	    for (i=0;i<c->iLineSize;i+=c->iDataSize) {
		    (*c->init)(c->ctx,&pLine[i]);
		}
	    }
	mdlSendThreadMessage(mdl,MDL_TAG_CACHECOM,iCore,creq,MDL_SE_CACHE_REQUEST);
	ret = 1;
	break;
    default:
	assert(0);
	}

    /*
     * Fire up next receive
     */
    id = MPI_ANY_SOURCE;
    MPI_Irecv(mpi->pszRcv,mpi->iCaBufSize, MPI_BYTE, id,
	      MDL_TAG_CACHECOM, mpi->commMDL, &mpi->ReqRcv);

    return ret;
    }



/*
** This routine must be called often by the MPI thread. It will drain
** any requests from the thread queue, and track MPI completions.
*/
static int checkMPI(MDL mdl) {
    mdlContextMPI *mpi = &mdl->mpi;
    while(1) {
	MDLserviceSend *send;
	MDLserviceCacheReq *caReq;
	MDLserviceElement *qhdr;
	cacheOpenClose *coc;
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
		iCore = send->iCoreFrom;
		assert(iCore>=0);
		tag = send->tag + MDL_TAG_THREAD_OFFSET * iCore;
		MPI_Irecv(send->buf,send->count,send->datatype,iProc,tag,mpi->commMDL,mpi->pSendRecvReq+mpi->nSendRecvReq);
		mpi->pSendRecvBuf[mpi->nSendRecvReq] = send;
		++mpi->nSendRecvReq;
		break;
		/* A thread has opened a cache -- we may need to resize our buffers */
	    case MDL_SE_CACHE_OPEN:
		coc = (cacheOpenClose *)qhdr;
		/* If this cache is larger than the largest, we have to reallocate buffers */
		if (coc->iDataSize > mpi->iMaxDataSize) {
		    if (mpi->ReqRcv != MPI_REQUEST_NULL) MPI_Cancel(&mpi->ReqRcv);
		    mpi->ReqRcv = MPI_REQUEST_NULL;
		    mpi->iMaxDataSize = coc->iDataSize;
		    mpi->iCaBufSize = (int)sizeof(CAHEAD) +
			  mpi->iMaxDataSize*(1 << MDL_CACHELINE_BITS);
		    mpi->pszRcv = realloc(mpi->pszRcv,mpi->iCaBufSize);
		    assert(mpi->pszRcv != NULL);
		    for (i=0;i<mdl->base.nThreads;++i) {
			mpi->ppszRpl[i] = realloc(mpi->ppszRpl[i],mdl->mpi.iCaBufSize);
			assert(mpi->ppszRpl[i] != NULL);
			}
		    mdl->pszFlsh = realloc(mdl->pszFlsh,mdl->mpi.iCaBufSize);
		    assert(mdl->pszFlsh != NULL);
		    }
		if (mpi->ReqRcv == MPI_REQUEST_NULL) {
		    MPI_Irecv(mpi->pszRcv,mpi->iCaBufSize, MPI_BYTE,
			MPI_ANY_SOURCE, MDL_TAG_CACHECOM,
			mpi->commMDL, &mpi->ReqRcv);
		    }
		++mpi->nOpenCaches;
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_OPEN);
		break;
		/* A thread has closed a cache -- we can cancel the receive if no threads have a cache open */
	    case MDL_SE_CACHE_CLOSE:
		assert(mpi->nOpenCaches > 0);
		--mpi->nOpenCaches;
		assert (mpi->ReqRcv != MPI_REQUEST_NULL);
		if (mpi->nOpenCaches == 0) {
		    MPI_Cancel(&mpi->ReqRcv);
		    mpi->ReqRcv = MPI_REQUEST_NULL;
		    }
		mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_CACHE_CLOSE);
		break;
		/* A thread needs to request a missing element */
	    case MDL_SE_CACHE_REQUEST:
		caReq = (MDLserviceCacheReq *)qhdr;
		assert(mpi->pThreadCacheReq[caReq->iCoreFrom]==NULL);
		iProc = mdlThreadToProc(mdl,caReq->caReq.idTo);
		mpi->pThreadCacheReq[caReq->iCoreFrom] = caReq;
		MPI_Send(&caReq->caReq,sizeof(CAHEAD),MPI_BYTE,iProc,MDL_TAG_CACHECOM, mpi->commMDL);
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
		switch(send->iServiceID) {
		case MDL_SE_MPI_SEND:
		case MDL_SE_MPI_SSEND:
		case MDL_SE_SEND_REQUEST:
		case MDL_SE_SEND_REPLY:
		    assert(send == &mdl->pmdl[send->iCoreFrom]->sendRequest);
		    tag = 0;
		    break;
		case MDL_SE_RECV_REPLY:
		    head = send->buf;
		    assert(mpi->pRequestTargets[send->target]>=0);
		    mpi->pRequestTargets[send->target] = -1;
		case MDL_SE_MPI_RECV:
		case MDL_SE_RECV_REQUEST:
		    assert(send == &mdl->pmdl[send->iCoreFrom]->recvRequest);
		    tag = send->tag % MDL_TAG_THREAD_OFFSET;
		    MPI_Get_count(&status, MPI_BYTE, &send->count);
		    break;
		default:
		    assert(0);
		    }
		mdlSendThreadMessage(mdl,tag,send->iCoreFrom,send,MDL_SE_MPI_SSEND);
		continue;
		}
	    }
	if (mpi->ReqRcv != MPI_REQUEST_NULL) {
	    while(1) {
		MPI_Test(&mdl->mpi.ReqRcv, &flag, &status);
		if (flag) mdlCacheReceive(mdl);                    
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
    if (!OPA_Queue_is_empty(&mdl->wq)) {
	/* Grab a work package, and perform it */
	OPA_Queue_dequeue(&mdl->wq, work, MDLwqNode, hdr);
	OPA_decr_int(&mdl->wqCurSize);
	while ( (*work->doFcn)(work->ctx) != 0 ) {
	    mdlCacheCheck(mdl);
	    }
	/* Send it back to the original thread for completion (could be us) */
	OPA_Queue_enqueue(&mdl->pmdl[work->iCoreOwner]->wqDone, work, MDLwqNode, hdr);
	rc = 1;
	}
    while (!OPA_Queue_is_empty(&mdl->wqDone)) {
	OPA_Queue_dequeue(&mdl->wqDone, work, MDLwqNode, hdr);
	(*work->doneFcn)(work->ctx);
	OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, hdr);
	}
    return rc;
    }

static void mdlCompleteAllWork(MDL mdl) {
    while(mdlDoSomeWork(mdl)) {}
    }

static void /*MDLserviceElement*/ *mdlWaitThreadQueue(MDL mdl,int iQueue) {
    MDLserviceElement *qhdr;
    while (OPA_Queue_is_empty(mdl->inQueue+iQueue)) {
	if (mdl->base.iCore == mdl->iCoreMPI) checkMPI(mdl);
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
static void mdlThreadBarrier(MDL mdl, void (*handler)(MDL mdl, MDLserviceElement *qhdr) ) {
    MDLserviceElement *qhdr;
    int i;

    if (mdl->base.iCore) {
	mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,0,&mdl->inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	int nExpected = mdl->base.nCores-1;
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
	mdlSendThreadMessage(mdl,0,qhdr->iCoreFrom,qhdr,MDL_SE_MPI_SSEND);
	}

    /* Off-node: We ask the MPI thread to post a receive for us. */
    else {
	return mdl_remote_MPI_Recv(buf,count,datatype,source,tag,mdl,nBytes,MDL_SE_MPI_RECV);
	}
    *nBytes = send->count;
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

    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    if (mdl->base.iCore) {
	mdlSendThreadMessage(mdl,MDL_TAG_BARRIER,0,&mdl->inMessage,MDL_SE_BARRIER_REQUEST);
	qhdr = mdlWaitThreadQueue(mdl,MDL_TAG_BARRIER);
	assert(qhdr->iServiceID == MDL_SE_BARRIER_REPLY);
	}
    else {
	int nExpected = mdl->base.nCores-1;
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

    /* Only used by the MPI thread */
    OPA_Queue_init(&mdl->mpi.queueMPI);
    OPA_Queue_header_init(&mdl->inMessage.hdr);
    OPA_Queue_header_init(&mdl->sendRequest.hdr);
    OPA_Queue_header_init(&mdl->recvRequest.hdr);
    mdl->mpi.nSendRecvReq = 0;
    mdl->mpi.pSendRecvReq = NULL;
    mdl->mpi.pSendRecvBuf = NULL;
    mdl->mpi.pThreadCacheReq = NULL;
    mdl->mpi.pRequestTargets = NULL;
    mdl->mpi.nRequestTargets = 0;
    mdl->mpi.iRequestTarget = 0;

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
    /*
    ** Initialize caching spaces.
    */
    mdl->cacheSize = MDL_CACHE_SIZE;
    for (i = 0; i<mdl->nMaxCacheIds; ++i) {
	mdl->cache[i].iType = MDL_NOCACHE;
#ifdef USE_ARC
	mdl->cache[i].arc = NULL;
#endif
	}

    /*
    ** Allocate caching buffers, with initial data size of 0.
    ** We need one reply buffer for each thread, to deadlock situations.
    */
    mdl->mpi.nOpenCaches = 0;
    mdl->mpi.iMaxDataSize = 0;
    mdl->mpi.iCaBufSize = sizeof(CAHEAD);
    mdl->mpi.pszRcv = malloc(mdl->mpi.iCaBufSize);
    assert(mdl->mpi.pszRcv != NULL);
    mdl->mpi.ppszRpl = malloc(mdl->base.nThreads*sizeof(char *));
    assert(mdl->mpi.ppszRpl != NULL);
    mdl->mpi.pmidRpl = malloc(mdl->base.nThreads*sizeof(int));
    assert(mdl->mpi.pmidRpl != NULL);
    for (i = 0; i<mdl->base.nThreads; ++i)
	mdl->mpi.pmidRpl[i] = -1;
    mdl->mpi.pReqRpl = malloc(mdl->base.nThreads*sizeof(MPI_Request));
    assert(mdl->mpi.pReqRpl != NULL);
    for (i = 0; i<mdl->base.nThreads; ++i) {
	mdl->mpi.ppszRpl[i] = malloc(mdl->mpi.iCaBufSize);
	assert(mdl->mpi.ppszRpl[i] != NULL);
	}
    mdl->pszFlsh = malloc(mdl->mpi.iCaBufSize);
    assert(mdl->pszFlsh != NULL);

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
    void *result = (*mdl->fcnWorker)(mdl);
//    MDLserviceElement stop;
    if (mdl->base.iCore != mdl->iCoreMPI) {
	mdlSendToMPI(mdl,&mdl->inMessage,MDL_SE_STOP);
	mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
	}
    else {
	drainMPI(mdl);
	}
    return result;
    }

int mdlLaunch(int argc,char **argv,void * (*fcnMaster)(MDL),void * (*fcnChild)(MDL)) {
    MDL mdl;
    int i,j,bDiag,bThreads,bDedicated,thread_support,rc,flag,*piTagUB;
    char *p, ach[256];

    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    mdlBaseInitialize(&mdl->base,argc,argv);

    /*
    ** Do some low level argument parsing for number of threads, and
    ** diagnostic flag!
    */
    for (argc = 0; argv[argc]; argc++);
    bDiag = 0;
    bThreads = 0;
    bDedicated = 0;
    i = 1;
    while (argv[i]) {
	if (!strcmp(argv[i], "-sz") && !bThreads) {
	    ++i;
	    mdl->base.nCores = atoi(argv[i]);
	    if (argv[i]) bThreads = 1;
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
    /* Dedicate one of the threads for MPI */
    if (bDedicated == 1) {
	if (mdl->base.nCores==1) bDedicated=0;
	else --mdl->base.nCores;
	}

    /* MPI Initialization */
    rc = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED,&thread_support);
    if (rc!=MPI_SUCCESS) {
	MPI_Error_string(rc, ach, &i);
	perror(ach);
	abort();
	}
    mdl->mpi.commMDL = MPI_COMM_WORLD;
    mdl->mpi.ReqRcv = MPI_REQUEST_NULL;
    MPI_Comm_size(mdl->mpi.commMDL, &mdl->base.nProcs);
    MPI_Comm_rank(mdl->mpi.commMDL, &mdl->base.iProc);

    /* Construct the thread/processor map */
    mdl->base.iProcToThread = malloc((mdl->base.nProcs + 1) * sizeof(int));
    assert(mdl->base.iProcToThread != NULL);
    mdl->base.iProcToThread[0] = 0;
    MPI_Allgather(&mdl->base.nCores, 1, MPI_INT, mdl->base.iProcToThread + 1, 1, MPI_INT, mdl->mpi.commMDL);
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

    /* Some bookeeping for the send/recv - 1 of each per thread */
    mdl->mpi.nSendRecvReq = 0;
    mdl->mpi.pSendRecvReq = malloc(mdl->base.nCores*2*sizeof(MPI_Request));
    assert(mdl->mpi.pSendRecvReq!=NULL);
    mdl->mpi.pSendRecvBuf = malloc(mdl->base.nCores*2*sizeof(MDLserviceSend *));
    assert(mdl->mpi.pSendRecvBuf!=NULL);
    mdl->mpi.pThreadCacheReq = malloc(mdl->base.nCores*sizeof(MDLserviceCacheReq *));
    assert(mdl->mpi.pThreadCacheReq!=NULL);
    for (i = 0; i < mdl->base.nCores; ++i) mdl->mpi.pThreadCacheReq[i] = NULL;

    /* Ring buffer of requests */
    mdl->mpi.iRequestTarget = 0;
    mdl->mpi.nRequestTargets = 2 * log2(1.0 * mdl->base.nThreads) + mdl->base.nCores;
    mdl->mpi.pRequestTargets = malloc(mdl->mpi.nRequestTargets * sizeof(*mdl->mpi.pRequestTargets));
    assert(mdl->mpi.pRequestTargets!=NULL);
    for(i=0; i<mdl->mpi.nRequestTargets; ++i) mdl->mpi.pRequestTargets[i] = -1;

    /* Make sure that MPI supports enough tags */
    rc = MPI_Comm_get_attr(mdl->mpi.commMDL,MPI_TAG_UB,&piTagUB,&flag);
    if (rc==MPI_SUCCESS && flag) {
	assert(mdl->base.nCores*MDL_TAG_THREAD_OFFSET < *piTagUB);
	assert(mdl->mpi.nRequestTargets*MDL_TAG_THREAD_OFFSET < *piTagUB);
	}

    /* Launch threads: if dedicated MPI thread then launch all worker threads. */
    mdl->mpi.nActiveCores = 0;
    if (mdl->base.nCores > 1 || bDedicated) {
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	for (i = mdl->iCoreMPI+1; i < mdl->base.nCores; ++i) {
	    pthread_create(&mdl->threadid[i], &attr,
		mdlWorkerThread,
		mdl->pmdl[i]);
	    ++mdl->mpi.nActiveCores;
	    }
	pthread_attr_destroy(&attr);
	}
    if (!bDedicated) {
	if (mdl->base.idSelf) (*fcnChild)(mdl);
	else fcnMaster(mdl);
	}
    drainMPI(mdl);
    mdlFinish(mdl);
    return(mdl->base.nThreads);
    }

void mdlFinish(MDL mdl) {
    int i;

    for (i = mdl->iCoreMPI+1; i < mdl->base.nCores; ++i) {
	pthread_join(mdl->threadid[i],0);
	}
    MPI_Barrier(mdl->mpi.commMDL);
    MPI_Finalize();


    for (i = 0; i<mdl->nMaxCacheIds; ++i) {
#ifdef USE_ARC
	if (mdl->cache[i].arc) arcFinish(mdl->cache[i].arc);
#endif
	}


    /*
     ** Close Diagnostic file.
     */
    if (mdl->base.bDiag) {
	fclose(mdl->base.fpDiag);
	}
    /*
     ** Deallocate storage.
     */
    free(mdl->pszIn);
    free(mdl->pszOut);
    free(mdl->pszBuf);
    free(mdl->pszTrans);
    free(mdl->cache);
    free(mdl->mpi.pszRcv);
    free(mdl->pszFlsh);
    for (i=0;i<mdl->base.nThreads;++i) free(mdl->mpi.ppszRpl[i]);
    free(mdl->mpi.ppszRpl);
    free(mdl->mpi.pmidRpl);
    free(mdl->mpi.pReqRpl);
    free(mdl->base.iProcToThread);
    mdlBaseFinish(&mdl->base);
    free(mdl);
    }


/*
** This needs to be improved by abstracting away more of the MPI functionality
*/

int mdlBcast ( MDL mdl, void *buf, int count, MDL_Datatype datatype, int root ) {
    return MPI_Bcast( buf, count, datatype, root, mdl->mpi.commMDL );
    }

int mdlScan ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op ) {
    return MPI_Scan( sendbuf, recvbuf, count, datatype, op, mdl->mpi.commMDL );
    }

int mdlExscan ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op ) {
    return MPI_Exscan( sendbuf, recvbuf, count, datatype, op, mdl->mpi.commMDL );
    }

int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root ) {
    return MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, mdl->mpi.commMDL );
    }

int mdlAllreduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op ) {
    return MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, mdl->mpi.commMDL );
    }

int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
		 void *recvbuf, int rcount, MDL_Datatype rtype) {
    return MPI_Alltoall(sendbuf,scount,stype,
			recvbuf,rcount,rtype,mdl->mpi.commMDL);
    }

int mdlAlltoallv( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype sendtype,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype recvtype) {
    return MPI_Alltoallv( sendbuf, sendcnts, sdispls, sendtype, 
        recvbuf, recvcnts, rdispls, recvtype, mdl->mpi.commMDL );
    }

int mdlAlltoallw( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype *stypes,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype *rtypes) {
    return MPI_Alltoallw( sendbuf, sendcnts, sdispls, stypes,
        recvbuf, recvcnts, rdispls, rtypes, mdl->mpi.commMDL );
    }

int mdlAllGather( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype recvtype) {
    return MPI_Allgather(sendbuf, scount, stype, recvbuf, rcount, recvtype, mdl->mpi.commMDL);
    } 

int mdlAllGatherv( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int *recvcnts, int *rdisps, MDL_Datatype recvtype) {
    return MPI_Allgatherv(sendbuf, scount, stype, recvbuf, recvcnts, rdisps, recvtype, mdl->mpi.commMDL);
    } 

int mdlReduceScatter( MDL mdl, void* sendbuf, void* recvbuf, int *recvcounts,
    MDL_Datatype datatype, MDL_Op op) {
    return MPI_Reduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, mdl->mpi.commMDL );
    }

int mdlTypeContiguous(MDL mdl,int count, MDL_Datatype old_type, MDL_Datatype *newtype) {
    return MPI_Type_contiguous(count,old_type,newtype);
    }

int mdlTypeIndexed(MDL mdl, int count,
    int array_of_blocklengths[], int array_of_displacements[],
    MDL_Datatype oldtype, MDL_Datatype *newtype) {
    return MPI_Type_indexed(count,
	array_of_blocklengths,array_of_displacements,
	oldtype,newtype);
    }

int mdlTypeCommit(MDL mdl, MDL_Datatype *datatype) {
    return MPI_Type_commit(datatype);
    }

int mdlTypeFree (MDL mdl, MDL_Datatype *datatype ) {
    return MPI_Type_free(datatype);
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
    MPI_Status status;
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
    size_t i;
    int nInMax,nOutMax,nBytes;
    int pid;
    char *pszBuf = vBuf;
    char *pszIn,*pszOut;
    struct swapInit {
	size_t nOutBytes;
	size_t nBufBytes;
	} swi,swo;
    MPI_Status status;

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
    mdlThreadBarrier(mdl,NULL);
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
    int iProc = mdlThreadToProc(mdl, id);
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
    MDLserviceElement *qhdr;
    char *pszIn = &mdl->pszBuf[sizeof(SRVHEAD)];
    int i,nBytes,id;
    MPI_Status status;
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
    MDLserviceElement *qhdr;
    char *pszIn = (char *)(phi + 1);
    char *pszOut = (char *)(pho + 1);
    int sid,id,tag,nOutBytes,nBytes;
    MPI_Status status;
    MPI_Comm   comm;
    int iProc;

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
	iProc = mdlThreadToProc(mdl,id);
	pho->nInBytes = phi->nInBytes;
	pho->nOutBytes = nOutBytes;
	tag = phi->replyTag;
	mdl_start_MPI_Ssend(pho, nOutBytes + sizeof(SRVHEAD), MPI_BYTE, id, tag, mdl,MDL_SE_SEND_REPLY);
	mdlWaitThreadQueue(mdl,0); /* Wait for Send to complete */
	} while (sid != SRV_STOP);
    }

void AdjustDataSize(MDL mdl) {
    int i,iMaxDataSize;

    /*
     ** Change buffer size?
     */
    iMaxDataSize = 0;
    for (i=0;i<mdl->nMaxCacheIds;++i) {
	if (mdl->cache[i].iType == MDL_NOCACHE) continue;
	if (mdl->cache[i].iDataSize > iMaxDataSize) {
	    iMaxDataSize = mdl->cache[i].iDataSize;
	    }
	}
    if (iMaxDataSize != mdl->mpi.iMaxDataSize) {
	/*
	 ** Create new buffer with realloc?
	 ** Be very careful when reallocing buffers in other libraries
	 ** (not PVM) to be sure that the buffers are not in use!
	 ** A pending non-blocking receive on a buffer which is realloced
	 ** here will cause problems, make sure to take this into account!
	 ** This is certainly true in using the MPL library.
	 */
	if (mdl->base.iCore == mdl->iCoreMPI)
	    MPI_Cancel(&mdl->mpi.ReqRcv);
	mdl->mpi.ReqRcv = MPI_REQUEST_NULL;
	mdl->mpi.iMaxDataSize = iMaxDataSize;
	mdl->mpi.iCaBufSize = (int)sizeof(CAHEAD) +
			  iMaxDataSize*(1 << MDL_CACHELINE_BITS);
	mdl->mpi.pszRcv = realloc(mdl->mpi.pszRcv,mdl->mpi.iCaBufSize);
	assert(mdl->mpi.pszRcv != NULL);
	for (i=0;i<mdl->base.nThreads;++i) {
	    mdl->mpi.ppszRpl[i] = realloc(mdl->mpi.ppszRpl[i],mdl->mpi.iCaBufSize);
	    assert(mdl->mpi.ppszRpl[i] != NULL);
	    }
	mdl->pszFlsh = realloc(mdl->pszFlsh,mdl->mpi.iCaBufSize);
	assert(mdl->pszFlsh != NULL);

	/*
	 * Fire up receive again.
	 */
	if (mdl->base.iCore == mdl->iCoreMPI)
	    MPI_Irecv(mdl->mpi.pszRcv,mdl->mpi.iCaBufSize, MPI_BYTE,
		MPI_ANY_SOURCE, MDL_TAG_CACHECOM,
		mdl->mpi.commMDL, &mdl->mpi.ReqRcv);
	}
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
    int id;
    CAHEAD caIn;
    int i,nMaxCacheIds;
    cacheOpenClose coc;

    /*
     ** Allocate more cache spaces if required!
     */
    assert(cid >= 0);
    if (cid >= mdl->nMaxCacheIds) {
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
#ifdef USE_ARC
	    mdl->cache[i].arc = NULL;
#endif
	    }
	mdl->nMaxCacheIds = nMaxCacheIds;
	}
    c = &mdl->cache[cid];
    assert(c->iType == MDL_NOCACHE);
    c->getElt = getElt==NULL ? getArrayElement : getElt;
    c->pData = pData;
    c->iDataSize = iDataSize;
    c->nData = nData;
    c->iLineSize = MDL_CACHELINE_ELTS*c->iDataSize;

    c->iKeyShift = 0;
    while ((1 << c->iKeyShift) < mdl->base.nThreads) ++c->iKeyShift;
    c->iIdMask = (1 << c->iKeyShift) - 1;

    if (c->iKeyShift < MDL_CACHELINE_BITS) {
	/*
	 * Key will be (index & MDL_INDEX_MASK) | id.
	 */
	c->iInvKeyShift = MDL_CACHELINE_BITS;
	c->iKeyShift = 0;
	}
    else {
	/*
	 * Key will be (index & MDL_INDEX_MASK) << KeyShift | id.
	 */
	c->iInvKeyShift = c->iKeyShift;
	c->iKeyShift -= MDL_CACHELINE_BITS;
	}

    /*
     ** Determine the number of cache lines to be allocated.
     */
    assert( mdl->cacheSize >= (1<<MDL_CACHELINE_BITS)*10*c->iDataSize);
    c->nLines = (mdl->cacheSize/c->iDataSize) >> MDL_CACHELINE_BITS;
    c->nTrans = 1;
    while (c->nTrans < c->nLines) c->nTrans *= 2;
    c->nTrans *= 2;
    c->iTransMask = c->nTrans-1;
    /*
     **	Set up the translation table.
     */
    c->pTrans = malloc(c->nTrans*sizeof(int));
    assert(c->pTrans != NULL);
    for (i=0;i<c->nTrans;++i) c->pTrans[i] = 0;
    /*
     ** Set up the tags. Note pTag[0] is a Sentinel!
     */
    c->pTag = malloc(c->nLines*sizeof(CTAG));
    assert(c->pTag != NULL);
    for (i=0;i<c->nLines;++i) {
	c->pTag[i].iKey = MDL_INVALID_KEY;
	c->pTag[i].nLock = 0;
	c->pTag[i].iLink = 0;
	}
    c->pTag[0].nLock = 1;     /* always locked */
    c->iLastVictim = 0;       /* This makes sure we have the first iVictim = 1 */
    c->nAccess = 0;
    c->nMiss = 0;				/* !!!, not NB */
    c->nColl = 0;				/* !!!, not NB */
    /*
     ** Allocate cache data lines.
     */
    c->pLine = malloc(mdl->cacheSize);
    assert(c->pLine != NULL);
#ifdef USE_ARC
    if (c->arc==NULL) c->arc = arcInitialize(mdl->cacheSize/iDataSize,iDataSize);
    c->pOneLine = malloc(c->iLineSize);
#endif

    /*
     ** Set up the request message as much as possible!
     */
    OPA_Queue_header_init(&c->cacheRequest.hdr);
    c->cacheRequest.iServiceID = MDL_SE_CACHE_REQUEST;
    c->cacheRequest.iCoreFrom = mdl->base.iCore;
    c->cacheRequest.caReq.cid = cid;
    c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
    c->cacheRequest.caReq.idFrom = mdl->base.idSelf;
    c->cacheRequest.caReq.idTo = -1;

    /* Read-only or combiner caches */
    c->iType = (init==NULL ? MDL_ROCACHE : MDL_COCACHE);
    c->init = init;
    c->combine = combine;
    c->ctx = ctx;

    /* Nobody should start using this cache until all threads have started it! */
    mdl_MPI_Barrier(mdl);

    /* We might need to resize the cache buffer */
    coc.iDataSize = c->iDataSize;
    mdlSendToMPI(mdl,&coc,MDL_SE_CACHE_OPEN);
    mdlWaitThreadQueue(mdl,0);

    AdjustDataSize(mdl);

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

void mdlFlushCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    CAHEAD caOut;
    CAHEAD *caFlsh = (CAHEAD *)mdl->pszFlsh;
    char *pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
    mdlkey_t iKey;
    int i,id;
    char *t;
    int j;
    int last;
    MPI_Status status;
    MPI_Request reqFlsh;
    MPI_Request reqBoth[2];
    int index;

    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    if (c->iType == MDL_COCACHE) {
	/*
	 ** Must flush all valid data elements.
	 */
	caFlsh->cid = cid;
	caFlsh->mid = MDL_MID_CACHEFLSH;
	caFlsh->idFrom = mdl->base.idSelf;
	for (i=1;i<c->nLines;++i) {
	    iKey = c->pTag[i].iKey;
	    if (iKey != MDL_INVALID_KEY) {
		/*
		 ** Flush element since it is valid!
		 */
		id = mdlkey_t_to_int(iKey & c->iIdMask);
		caFlsh->iLine = mdlkey_t_to_int(iKey >> c->iInvKeyShift);
		t = &c->pLine[i*c->iLineSize];
		for (j = 0; j < c->iLineSize; ++j)
		    pszFlsh[j] = t[j];
		/*
		 * Use Synchronous send so as not to
		 * overwhelm the receiver.
		 */
		MPI_Issend(caFlsh, (int)sizeof(CAHEAD)+c->iLineSize,
			   MPI_BYTE, id, MDL_TAG_CACHECOM,
			   mdl->mpi.commMDL, &reqFlsh);
		/*
		 * Wait for the Flush to complete, but
		 * also service any incoming cache requests.
		*/
		reqBoth[0] = mdl->mpi.ReqRcv;
		reqBoth[1] = reqFlsh;

		while (1) {
		    MPI_Waitany(2, reqBoth, &index, &status);
		    assert(!(index != 0 && reqBoth[0] ==
			     MPI_REQUEST_NULL));
		    mdl->mpi.ReqRcv = reqBoth[0];
		    if (index == 1) /* Flush has completed */
			break;
		    else if (index == 0) {
			mdlCacheReceive(mdl);
			reqBoth[0] = mdl->mpi.ReqRcv;
			}
		    else
			assert(0);
		    }
		}
		c->pTag[i].iKey = MDL_INVALID_KEY;
	    }
	}
    mdl->wqAccepting = 0;
    mdlTimeAddSynchronizing(mdl);
    }

#ifdef USE_ARC
static void arcFlushAll(ARC arc) {
    int i;
    CDB *temp;
    for(;arc->T1Length;--arc->T1Length) {
	temp = lru_remove(arc->T1);
	assert(temp->data != NULL);
	assert(temp->data[-1] == _ARC_MAGIC_);
	remove_from_hash(arc,temp);
	lru_insert(temp,arc->Free);
	}
    for(;arc->T2Length;--arc->T2Length) {
	temp = lru_remove(arc->T2);
	assert(temp->data != NULL);
	assert(temp->data[-1] == _ARC_MAGIC_);
	remove_from_hash(arc,temp);
	lru_insert(temp,arc->Free);
	}
    for(;arc->B1Length;--arc->B1Length) {
	temp = lru_remove(arc->B1);
	assert(temp->data == NULL);
	remove_from_hash(arc,temp);
	mru_insert(temp,arc->Free);
	}
    for(;arc->B2Length;--arc->B2Length) {
	temp = lru_remove(arc->B2);
	assert(temp->data == NULL);
	remove_from_hash(arc,temp);
	mru_insert(temp,arc->Free);
	}
    assert(arc->T1->next == arc->T1);
    assert(arc->T1->prev == arc->T1);
    assert(arc->B1->next == arc->B1);
    assert(arc->B1->prev == arc->B1);
    assert(arc->T2->next == arc->T2);
    assert(arc->T2->prev == arc->T2);
    assert(arc->B2->next == arc->B2);
    assert(arc->B2->prev == arc->B2);
    }
#endif

void mdlFinishCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    CAHEAD caOut;
    CAHEAD *caFlsh = (CAHEAD *)mdl->pszFlsh;
    char *pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
    cacheOpenClose coc;
    mdlkey_t iKey;
    int i,id;
    char *t;
    int j;
    MPI_Status status;
    MPI_Request reqFlsh;
    MPI_Request reqBoth[2];
    int index;

    mdlTimeAddComputing(mdl);
    mdl->wqAccepting = 1;
    mdlFlushCache(mdl,cid);

    /* We must wait for all threads to finish with this cache */
    coc.iDataSize = c->iDataSize;
    mdlSendToMPI(mdl,&coc,MDL_SE_CACHE_CLOSE);
    mdlWaitThreadQueue(mdl,0);
    mdl_MPI_Barrier(mdl);

    /*
     ** Free up storage and finish.
     */
    free(c->pTrans);
    free(c->pTag);
    free(c->pLine);
#ifdef USE_ARC
    free(c->pOneLine);
    arcFlushAll(c->arc);
#ifdef FINISH_ARC
    arcFinish(c->arc);
    c->arc = NULL;
#endif
#endif
    c->iType = MDL_NOCACHE;
    mdl->wqAccepting = 0;
    mdlTimeAddSynchronizing(mdl);
    }

void mdlPrefetch(MDL mdl,int cid,int iIndex, int id) {
    }

#ifdef USE_ARC
static inline uint64_t *replace(ARC arc, int iInB2) {
    CDB *temp;
    uint64_t *data;
    uint32_t max = (arc->target_T1 > 1)?(arc->target_T1+1-iInB2):1;
    if (arc->T1Length >= max) { /* T1’s size exceeds target? */
                                        /* yes: T1 is too big */
	temp = arc->T1->prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->prev;
	    if (temp == arc->T1) {           /* all pages in T1 are currently locked, try T2 in this case */
		temp = arc->T2->prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->prev;
		    if (temp == arc->T2) return(NULL); /* all pages are currently locked, give up! */
		}
		goto replace_T2;
	    }
	}
	goto replace_T1;
    } else {
	/* no: T1 is not too big */
	temp = arc->T2->prev;                /* get LRU */
	while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
	    temp = temp->prev;
	    if (temp == arc->T2) {           /* all pages in T2 are currently locked, try T1 in this case */
		temp = arc->T1->prev;                /* get LRU */
		while (temp->data[-1] != _ARC_MAGIC_) {           /* is it a locked page? */
		    temp = temp->prev;
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
#if (0)
    if (temp->uId & _DIRTY_) destage(temp);     /* if dirty, evict before overwrite */
#endif
        remove_from_list(temp);          /* grab LRU unlocked page from T1 */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
	mru_insert(temp,arc->B1);           /* put it on B1 */
	temp->uId = (temp->uId&_IDMASK_)|_B1_;  /* need to be careful here because it could have been _P1_ */
        arc->T1Length--; arc->B1Length++;          /* bookkeep */
	/*assert(arc->B1Length<=arc->nCache);*/
    }
    if (0) {
    replace_T2:
#if (0)
    if (temp->uId & _DIRTY_) destage(temp);     /* if dirty, evict before overwrite */
#endif
        remove_from_list(temp);          /* grab LRU unlocked page from T2 */
	data = temp->data;
	temp->data = NULL; /*GHOST*/
        mru_insert(temp,arc->B2);           /* put it on B2 */
        temp->uId |= _B2_;          /* note that fact */
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
static inline void arcSetPrefetchData(ARC arc,uint32_t uIndex,uint32_t uId,void *data) {
    CDB *temp;
    uint32_t uHash,rat;
    uint32_t tuId = uId&_IDMASK_;
    uint32_t L1Length, L2Length;
    int inB2=0;

    uHash = (MurmurHash2(uIndex,tuId)&arc->uHashMask);
    temp = arc->Hash[uHash];
    while (temp) {
	if (temp->uIndex == uIndex) {
	    if ((temp->uId&_IDMASK_) == tuId) break;
	}
	temp = temp->coll;
    }
    if (temp != NULL) {                       /* found in cache directory? */
	switch (temp->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	case _T1_:
	case _T2_:
	    return;              /* do nothing */
	case _B1_:                            /* B1 "hit": don't favour anything */
	    arc->B1Length--;                                           /* bookkeep */
	    goto doBcase;
	case _B2_:                            /* B2 "hit": don't favour anything */
	    arc->B2Length--;                                           /* bookkeep*/
	    inB2=1;
	doBcase:
	    /* Better would be to put this back on P1, but L1 may be full. */
	    remove_from_list(temp);                                   /* take off whichever list */
	    temp->data = replace(arc,inB2);                                /* find a place to put new page */
	    temp->uIndex = uIndex;                          /* bookkeep */
	    if (inB2) {
		temp->uId = _T2_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
		mru_insert(temp,arc->T2);                                     /* not seen yet, put on T1 */
		arc->T2Length++;
		}
	    else {
		temp->uId = _P1_|(uId&_IDMASK_);     /* temp->ARC_where = _P1_; and clear the dirty bit for this page */
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
		temp = lru_remove(arc->B1);                                    /* yes: take page off B1 */
		remove_from_hash(arc,temp);            /* remove from hash table */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = replace(arc,0);                                /* find new place to put page */
	    } else {                                                      /* no: B1 must be empty */
		temp = lru_remove(arc->T1);                                    /* take page off T1 */
#if (0)
		if (temp->dirty) destage(temp);                           /* if dirty, evict before overwrite */
#endif
		remove_from_hash(arc,temp);            /* remove from hash table */
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
		    temp = lru_remove(arc->B2);            /* here we lose memory of what was in lru B2 */
		    remove_from_hash(arc,temp);            /* remove from hash table */
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=1;
		} else {                                                   /* cache directory not full, easy case */
		    temp = lru_remove(arc->Free);
		}
		temp->data = replace(arc,inB2);                                /* new place for page */
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
	temp->uIndex = uIndex;
	temp->coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */
    }
    memcpy(temp->data,data,arc->uDataSize*sizeof(uint64_t));
}


/*
** This releases a lock on a page by decrementing the lock counter for that page. It also checks
** that the magic number matches, and that the lock count is greater than 0. This is to prevent
** strange errors if the user tries to release the wrong pointer (could silently modify memory 
** that it should not if this check was not there).
*/
static inline void arcRelease(ARC arc,uint64_t *p) {
    if (p>arc->dataBase && p<arc->dataLast) {
	uint64_t t = p[-1]-1;
	assert((t^_ARC_MAGIC_) < 0x00000000ffffffff);
	p[-1] = t;
	}
}

static void queueCacheRequest(MDL mdl, int cid, int iIndex, int id) {
    CACHE *c = &mdl->cache[cid];
    int iLine = iIndex >> MDL_CACHELINE_BITS;
    int iElt = iIndex & MDL_CACHE_MASK;
    int s = iIndex & MDL_INDEX_MASK;
    c->cacheRequest.caReq.cid = cid;
    c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
    c->cacheRequest.caReq.idFrom = mdl->base.idSelf;
    c->cacheRequest.caReq.idTo = id;
    c->cacheRequest.caReq.iLine = iLine;
    c->cacheRequest.pLine = c->pOneLine;
    mdlSendToMPI(mdl,&c->cacheRequest,MDL_SE_CACHE_REQUEST);
    }

static void finishCacheRequest(MDL mdl, int cid, CDB *temp) {
    CACHE *c = &mdl->cache[cid];
    ARC arc = c->arc;
    uint32_t uIndex = temp->uIndex;
    uint32_t uId = temp->uId & _IDMASK_;
    int iLine = uIndex >> MDL_CACHELINE_BITS;
    int iElt = uIndex & MDL_CACHE_MASK;
    int s = uIndex & MDL_INDEX_MASK;

    assert(temp->data!=NULL);
    mdlWaitThreadQueue(mdl,MDL_TAG_CACHECOM);
    memcpy(temp->data,c->pOneLine + iElt*c->iDataSize, c->iDataSize);

#if 1
    ++temp->data[-1];       /* prefetch must never evict our data */
    int i;
    for(i=0; i<MDL_CACHELINE_ELTS; ++i) {
	int nid = s + i;
	if (nid != uIndex)
	    arcSetPrefetchData(arc,nid,uId,c->pOneLine + i*c->iDataSize);
	}
    --temp->data[-1];       /* may lock below */
#endif
    }

/*
** The highest order bit (1<<31) of uId encodes the dirty bit.
** The next 3 highest order bits of uId ((1<<30), (1<<29) and (1<<28)) are reserved 
** for list location and should be zero! The maximum legal uId is then (1<<28)-1.
*/
void *mdlDoMiss(MDL mdl, int cid, int iIndex, int id, mdlkey_t iKey, int bLock) {
    CACHE *c = &mdl->cache[cid];
    ARC arc = c->arc;
    uint32_t uIndex = iIndex;
    uint32_t uId = id;
    CDB *temp;
    uint32_t L1Length, L2Length;
    uint32_t uHash,rat;
    uint32_t tuId = uId&_IDMASK_;
    int inB2=0;

    uHash = (MurmurHash2(uIndex,tuId)&arc->uHashMask);
    temp = arc->Hash[uHash];
    while (temp) {
	if (temp->uIndex == uIndex) {
	    if ((temp->uId&_IDMASK_) == tuId) break;
	    }
	temp = temp->coll;
	}
    if (temp != NULL) {                       /* found in cache directory? */
	switch (temp->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	    temp->uId = uId;     /* clears prefetch flag and sets WHERE = _T1_ and dirty bit */
	    remove_from_list(temp);                           /* take off T1 list */
	    mru_insert(temp,arc->T1);                         /* first access but recently prefetched, put on T1 */
	    goto cachehit;
	case _T1_:
	    arc->T1Length--; arc->T2Length++;
	    temp->uId |= _T2_;
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
	    return((void *)temp->data);
	case _B1_:                            /* B1 hit: favor recency */
	    /*
	    ** Can initiate the data request right here, and do the rest while waiting...
	    */
	    queueCacheRequest(mdl,cid,iIndex,id);

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
	    queueCacheRequest(mdl,cid,iIndex,id);

	    rat = arc->B1Length/arc->B2Length;
	    if (rat < 1) rat = 1;
	    if (rat > arc->target_T1) arc->target_T1 = 0;
	    else arc->target_T1 = arc->target_T1 - rat;
	    /* adapt the target size */
	    arc->B2Length--;                                           /* bookkeep */
	    inB2=1;
	doBcase:
	    remove_from_list(temp);                                   /* take off whichever list */
	    temp->data = replace(arc,inB2);                                /* find a place to put new page */
	    temp->uId = _T2_|uId;     /* temp->ARC_where = _T2_; and set the dirty bit for this page */
	    temp->uIndex = uIndex;                          /* bookkeep */
	    mru_insert(temp,arc->T2);                                     /* seen twice recently, put on T2 */
	    arc->T2Length++;                 /* JS: this was not in the original code. Should it be? bookkeep */
	    /*assert(temp->data!=NULL);*/
	    finishCacheRequest(mdl,cid,temp);
	    break;
	    }
	}

    else {                                                              /* page is not in cache directory */
	++c->nMiss;
	mdlTimeAddComputing(mdl);
	/*
	** Can initiate the data request right here, and do the rest while waiting...
	*/
	queueCacheRequest(mdl,cid,iIndex,id);
	L1Length = arc->T1Length + arc->B1Length;
	/*assert(L1Length<=arc->nCache);*/
	if (L1Length == arc->nCache) {                                   /* B1 + T1 full? */
	    if (arc->T1Length < arc->nCache) {                                           /* Still room in T1? */
		temp = lru_remove(arc->B1);                                    /* yes: take page off B1 */
		remove_from_hash(arc,temp);            /* remove from hash table */
		arc->B1Length--;                                               /* bookkeep that */
		temp->data = replace(arc,0);                                /* find new place to put page */
		}
	    else {                                                      /* no: B1 must be empty */
		temp = lru_remove(arc->T1);                                    /* take page off T1 */
#if (0)
		if (temp->dirty) destage(temp);                           /* if dirty, evict before overwrite */
#endif
		remove_from_hash(arc,temp);            /* remove from hash table */
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
		    temp = lru_remove(arc->B2);
		    remove_from_hash(arc,temp);            /* remove from hash table */
		    arc->B2Length--;                                           /* find and reuse B2’s LRU */
		    inB2=1;
		} else {                                                   /* cache directory not full, easy case */
		    temp = lru_remove(arc->Free);
		    assert(temp->data == NULL);            /* This CDB should not be associated with data */
		}
		temp->data = replace(arc,inB2);                                /* new place for page */
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
	temp->uIndex = uIndex;
	temp->coll = arc->Hash[uHash];                  /* add to collision chain */
	arc->Hash[uHash] = temp;                               /* insert into hash table */

	finishCacheRequest(mdl,cid,temp);

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
    return((void *)temp->data);
}
#else
void *mdlDoMiss(MDL mdl, int cid, int iIndex, int id, mdlkey_t iKey, int lock) {
    CACHE *c = &mdl->cache[cid];
    char *pLine;
    int iElt,i,iLine;
    mdlkey_t iKeyVic;
    int idVic;
    int iVictim,*pi;
    char ach[80];
    CAHEAD *caFlsh;
    char *pszFlsh;
    MPI_Status status;
    MPI_Request reqFlsh;

    /* Short circuit the cache if this belongs to another thread */
    int iCore = id - mdl->pmdl[0]->base.idSelf;
    if (iCore >= 0 && iCore < mdl->base.nCores && c->iType == MDL_ROCACHE ) {
	MDL omdl = mdl->pmdl[iCore];
	c = &omdl->cache[cid];
	return (*c->getElt)(c->pData,iIndex,c->iDataSize);
	}

    mdlTimeAddComputing(mdl);

    /*
    ** Cache Miss.
    */
    iLine = iIndex >> MDL_CACHELINE_BITS;
    c->cacheRequest.caReq.cid = cid;
    c->cacheRequest.caReq.mid = MDL_MID_CACHEREQ;
    c->cacheRequest.caReq.idFrom = mdl->base.idSelf;
    c->cacheRequest.caReq.idTo = id;
    c->cacheRequest.caReq.iLine = iLine;
    /* We could send the request now, but we don't know the target cache line yet. */

    ++c->nMiss;

    if (++c->iLastVictim == c->nLines) c->iLastVictim = 1;
    iVictim = c->iLastVictim;

    iElt = iIndex & MDL_CACHE_MASK;
    for (i=1;i<c->nLines;++i) {
	if (!c->pTag[iVictim].nLock) {
	    /*
	     ** Found victim.
	     */
	    iKeyVic = c->pTag[iVictim].iKey;
	    /*
	     ** 'pLine' will point to the actual data line in the cache.
	     */
	    pLine = &c->pLine[iVictim*c->iLineSize];
	    caFlsh = NULL;
	    if (iKeyVic != MDL_INVALID_KEY) {
		if (c->iType == MDL_COCACHE) {
		    /*
		    ** Flush element since it is valid!
		    */
		    idVic = mdlkey_t_to_int(iKeyVic&c->iIdMask);
		    caFlsh = (CAHEAD *)mdl->pszFlsh;
		    pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
		    caFlsh->cid = cid;
		    caFlsh->mid = MDL_MID_CACHEFLSH;
		    caFlsh->idFrom = mdl->base.idSelf;
		    caFlsh->iLine = mdlkey_t_to_int(iKeyVic >> c->iInvKeyShift);
		    for (i = 0; i < c->iLineSize; ++i)
			pszFlsh[i] = pLine[i];
		    MPI_Isend(caFlsh, (int)sizeof(CAHEAD)+c->iLineSize,
			      MPI_BYTE, idVic,
			      MDL_TAG_CACHECOM, mdl->mpi.commMDL, &reqFlsh);
		    }
		/*
		 ** If valid iLine then "unlink" it from the cache.
		 */
		pi = &c->pTrans[iKeyVic & c->iTransMask];
		while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
		*pi = c->pTag[iVictim].iLink;
		}
	    c->pTag[iVictim].iKey = iKey;
	    if (lock) c->pTag[iVictim].nLock = 1;
	    /*
	     **	Add the modified victim tag back into the cache.
	     ** Note: the new element is placed at the head of the chain.
	     */
	    pi = &c->pTrans[iKey & c->iTransMask];
	    c->pTag[iVictim].iLink = *pi;
	    *pi = iVictim;
	    goto Await;
	    }
	if (++iVictim == c->nLines) iVictim = 1;
	}
    /*
     ** Cache Failure!
     */
    sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
    mdlDiag(mdl,ach);
    exit(1);
Await:
    c->iLastVictim = iVictim;
    /*
     ** At this point 'pLine' is the recipient cache line for the
     ** data requested from processor 'id'.
     */
    c->cacheRequest.pLine = pLine;
    mdlSendToMPI(mdl,&c->cacheRequest,MDL_SE_CACHE_REQUEST);
    mdlWaitThreadQueue(mdl,MDL_TAG_CACHECOM);

    mdlTimeAddWaiting(mdl);
    return(&pLine[iElt*c->iDataSize]);
    }
#endif

void *mdlAquire(MDL mdl,int cid,int iIndex,int id) {
    const int lock = 1;  /* we always lock in aquire */
    CACHE *c = &mdl->cache[cid];
    char *pLine;
    mdlkey_t iKey;
    int iElt;
    int i;

    ++c->nAccess;

    if (!(c->nAccess & MDL_CHECK_MASK))
	mdlCacheCheck(mdl);

    /* Short circuit the cache if this belongs to another thread (or ourselves) */
    int iCore = id - mdl->pmdl[0]->base.idSelf;
    if (iCore >= 0 && iCore < mdl->base.nCores && c->iType == MDL_ROCACHE ) {
	MDL omdl = mdl->pmdl[iCore];
	c = &omdl->cache[cid];
	return (*c->getElt)(c->pData,iIndex,c->iDataSize);
	}

#ifdef USE_ARC
#else
    /*
     ** Determine memory block key value and cache line.
     */
    iKey = iIndex & MDL_INDEX_MASK;
    iKey <<= c->iKeyShift;
    iKey |= id;
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
#endif
    return(mdlDoMiss(mdl, cid, iIndex, id, iKey, lock));
    }


void mdlRelease(MDL mdl,int cid,void *p) {
    CACHE *c = &mdl->cache[cid];
#ifdef USE_ARC
    arcRelease(c->arc,p);
#else
    int iLine;

    iLine = ((char *)p - c->pLine) / c->iLineSize;
    /*
     ** Check if the pointer fell in a cache line, otherwise it
     ** must have been a local pointer.
     */
    if (iLine > 0 && iLine < c->nLines) {
	--c->pTag[iLine].nLock;
	assert(c->pTag[iLine].nLock >= 0);
	}
#endif
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
    grid->rs = mdlMalloc(mdl,sizeof(*grid->rs)*mdl->base.nThreads); assert(grid->rs!=NULL);
    grid->rn = mdlMalloc(mdl,sizeof(*grid->rn)*mdl->base.nThreads); assert(grid->rn!=NULL);

    /* The following need to be set to appropriate values still. */
    grid->s = grid->n = grid->nlocal = 0;
    }

void mdlGridFinish(MDL mdl, MDLGRID grid) {
    if (grid->rs) free(grid->rs);
    if (grid->rn) free(grid->rn);
    if (grid->id) free(grid->id);
    free(grid);
    }

void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, int nlocal) {
    assert( s>=0 && s<grid->n3);
    assert( n>=0 && s+n<=grid->n3);
    grid->s = s;
    grid->n = n;
    grid->nlocal = nlocal;
    }

/*
** Share the local GRID information with other processors by,
**   - finding the starting slab and number of slabs on each processor
**   - building a mapping from slab to processor id.
*/
void mdlGridShare(MDL mdl,MDLGRID grid) {
    int i, id;
#if 0
    MPI_Allgather(&grid->s,sizeof(*grid->rs),MPI_BYTE,
		  grid->rs,sizeof(*grid->rs),MPI_BYTE,
		  mdl->commMDL);
    MPI_Allgather(&grid->n,sizeof(*grid->rn),MPI_BYTE,
		  grid->rn,sizeof(*grid->rn),MPI_BYTE,
		  mdl->commMDL);
#else
    assert(0);
#endif
    /* Calculate on which processor each slab can be found. */
    for(id=0; id<mdl->base.nThreads; id++ ) {
	for( i=grid->rs[id]; i<grid->rs[id]+grid->rn[id]; i++ ) grid->id[i] = id;
	}
    }

/*
** Allocate the local elements.  The size of a single element is
** given and the local GRID information is consulted to determine
** how many to allocate.
*/
void *mdlGridMalloc(MDL mdl,MDLGRID grid,int nEntrySize) {
    return mdlMalloc(mdl,nEntrySize*grid->nlocal);
    }

void mdlGridFree( MDL mdl, MDLGRID grid, void *p ) {
    mdlFree(mdl,p);
    }

#ifdef MDL_FFTW
size_t mdlFFTInitialize(MDL mdl,MDLFFT *pfft,
			int n1,int n2,int n3,int bMeasure) {
    MDLFFT fft;
    int sz,nz,sy,ny,nlocal;

    *pfft = NULL;
    fft = malloc(sizeof(struct mdlFFTContext));
    assert(fft != NULL);

    /* Contruct the FFTW plans */
    fft->fplan = rfftw3d_mpi_create_plan(mdl->commMDL,
					 n3, n2, n1,
					 FFTW_REAL_TO_COMPLEX,
					 (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE) );
    fft->iplan = rfftw3d_mpi_create_plan(mdl->commMDL,
					 /* dim.'s of REAL data --> */ n3, n2, n1,
					 FFTW_COMPLEX_TO_REAL,
					 (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE));
    rfftwnd_mpi_local_sizes( fft->fplan, &nz, &sz, &ny, &sy,&nlocal);

    /*
    ** Dimensions of k-space and r-space grid.  Note transposed order.
    ** Note also that the "actual" dimension 1 side of the r-space array
    ** can be (and usually is) larger than "n1" because of the inplace FFT.
    */
    mdlGridInitialize(mdl,&fft->rgrid,n1,n2,n3,2*(n1/2+1));
    mdlGridInitialize(mdl,&fft->kgrid,n1/2+1,n3,n2,n1/2+1);

    mdlGridSetLocal(mdl,fft->rgrid,sz,nz,nlocal);
    mdlGridSetLocal(mdl,fft->kgrid,sy,ny,nlocal/2);
    mdlGridShare(mdl,fft->rgrid);
    mdlGridShare(mdl,fft->kgrid);

    *pfft = fft;
    return nlocal;
    }

void mdlFFTFinish( MDL mdl, MDLFFT fft ) {
    rfftwnd_mpi_destroy_plan(fft->fplan);
    rfftwnd_mpi_destroy_plan(fft->iplan);
    mdlGridFinish(mdl,fft->kgrid);
    mdlGridFinish(mdl,fft->rgrid);
    free(fft);
    }

fftw_real *mdlFFTMAlloc( MDL mdl, MDLFFT fft ) {
    return mdlGridMalloc(mdl,fft->rgrid,sizeof(fftw_real));
    }

void mdlFFTFree( MDL mdl, MDLFFT fft, void *p ) {
    mdlGridFree(mdl,fft->rgrid,p);
    }

void mdlFFT( MDL mdl, MDLFFT fft, fftw_real *data, int bInverse ) {
    rfftwnd_mpi_plan plan = bInverse ? fft->iplan : fft->fplan;
    rfftwnd_mpi(plan,1,data,0,FFTW_TRANSPOSED_ORDER);
    }
#endif

void mdlSetWorkQueueSize(MDL mdl,int wqMaxSize,int cudaSize) {
    MDLwqNode *work;
    while (wqMaxSize > mdl->wqMaxSize) {
	work = malloc(sizeof(MDLwqNode));
	OPA_Queue_header_init(&work->hdr);
	work->iCoreOwner = mdl->base.iCore;
	OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, hdr);
	work = malloc(sizeof(MDLwqNode));
	OPA_Queue_header_init(&work->hdr);
	work->iCoreOwner = mdl->base.iCore;
	OPA_Queue_enqueue(&mdl->wqFree, work, MDLwqNode, hdr);
	++mdl->wqMaxSize;
	}
    while (wqMaxSize < mdl->wqMaxSize && !OPA_Queue_is_empty(&mdl->wqFree)) {
	OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, hdr);
	free(work);
	OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, hdr);
	free(work);
	--mdl->wqMaxSize;
	}
    /* This needs to be called when there is no active work so see if we removed enough */
    assert(mdl->wqMaxSize == wqMaxSize);
    }

void mdlAddWork(MDL mdl, void *ctx, mdlWorkFunction initWork, mdlWorkFunction checkWork, mdlWorkFunction doWork, mdlWorkFunction doneWork) {
    MDL Qmdl = NULL;
    int i;

    /* Obviously, we can only queue work if we have a free queue element */
    if (!OPA_Queue_is_empty(&mdl->wqFree)) {
	/* We have some room, so save work for later */
	if (OPA_load_int(&mdl->wqCurSize) < mdl->wqMaxSize) Qmdl = mdl;
	/* See if anyone else is accepting work */
	else {
	    i = mdl->wqLastHelper;
	    do {
		MDL rMDL = mdl->pmdl[i];
		if (rMDL->wqAccepting && OPA_load_int(&rMDL->wqCurSize)<=1) {
		    Qmdl = rMDL;
		    mdl->wqLastHelper = i;
		    break;
		    }
		if (++i == mdl->base.nCores) i = 0;
		} while(i!=mdl->wqLastHelper);
	    }
	if (Qmdl) {
	    MDLwqNode *work;
	    OPA_Queue_dequeue(&mdl->wqFree, work, MDLwqNode, hdr);
	    work->ctx = ctx;
	    work->checkFcn = checkWork;
	    work->doFcn = doWork;
	    work->doneFcn = doneWork;
	    OPA_incr_int(&Qmdl->wqCurSize);
	    OPA_Queue_enqueue(&Qmdl->wq, work, MDLwqNode, hdr);
	    return;
	    }

	}

    /* Just handle it ourselves */
    while( doWork(ctx) != 0 ) {}
    doneWork(ctx);

    }
