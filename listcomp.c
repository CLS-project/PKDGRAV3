#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include "listcomp.h"

/*
#define LISTCOMP_TEST
*/

#define LIST_GROW  100

uint32_t msb32(register uint32_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return(x & ~(x >> 1));
}

uint32_t swar32(register uint32_t x)
{
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return(x);
}

uint32_t ones32(register uint32_t x)
{
    /* 32-bit recursive reduction using SWAR...
       but first step is mapping 2-bit values
       into sum of 2 1-bit values in sneaky way
    */
    x -= ((x >> 1) & 0x55555555);
    x = (((x >> 2) & 0x33333333) + (x & 0x33333333));
    x = (((x >> 4) + x) & 0x0f0f0f0f);
    x += (x >> 8);
    x += (x >> 16);
    return(x & 0x0000003f);
}


void lcodePrintList(LIST *p,int nList) {
    int i;
    
    for (i=0;i<nList;) {
	printf("%d:%d ",p[i].iPid,p[i].iIndex);
	if ((++i)%20 == 0) printf("\n");
    }
    if (i%20 != 0) printf("\n");
}

void uint2bin(uint32_t s,char *bin)
{
    int i;

    for (i=0;i<32;++i) {
	bin[i] = '0' + (s&1);
	s = s>>1;
    }
    bin[32] = 0;
}


int lcodeCmpList(const void *v1,const void *v2) {
    LIST *p1 = (LIST *)v1;
    LIST *p2 = (LIST *)v2;
    int iRet;
    iRet = p1->iPid-p2->iPid;
    if (iRet) return(iRet);
    else return(p1->iIndex-p2->iIndex);
}


LCODE lcodeInit(uint32_t nThreads,uint32_t idSelf,uint32_t nLocal,uint32_t nSmooth) {
    LCODE ctx;

    ctx = malloc(sizeof(struct lcodeContext));
    assert(ctx != NULL);
    ctx->mPid = swar32(nThreads);
    ctx->mPrefix = swar32(nLocal);
    ctx->mSuffix = swar32(nSmooth);
    ctx->mSuffix >>= 0;  /* shift back to an interval with most likely runs */
    ctx->mPrefix ^= ctx->mSuffix;
    ctx->nPid = ones32(ctx->mPid);
    ctx->idSelf = idSelf;
    ctx->nPrefix = ones32(ctx->mPrefix);
    ctx->nSuffix = ones32(ctx->mSuffix);
    printf("\nnBitsPrefix:%d nBitsSuffix:%d\n",ctx->nPrefix,ctx->nSuffix);
    ctx->aCode = NULL;
    ctx->nCode = 0;
    ctx->uIndex = 0;
    ctx->uMask = 1;
    ctx->inCode = NULL;
    return ctx;
};


void lcodeFinish(LCODE ctx) {
    if (ctx->aCode) free(ctx->aCode);
    free(ctx);
}


void OutOne(LCODE ctx) {
    if (ctx->uIndex == ctx->nCode) {
	ctx->nCode += 256;
	ctx->aCode = realloc(ctx->aCode,ctx->nCode);
	assert(ctx->aCode != NULL);
    }
    ctx->aCode[ctx->uIndex] |= ctx->uMask;
    ctx->uMask <<= 1;
    if (!ctx->uMask) {
	++ctx->uIndex;
	ctx->uMask = 1;
    }
}

int InOne(LCODE ctx) {
    int iRet;
    iRet = (ctx->inCode[ctx->uIndex] & ctx->uMask);
    ctx->uMask <<= 1;
    if (!ctx->uMask) {
	++ctx->uIndex;
	ctx->uMask = 1;
    }
    return(iRet);
}

void OutZero(LCODE ctx) {
    if (ctx->uIndex == ctx->nCode) {
	ctx->nCode += 256;
	ctx->aCode = realloc(ctx->aCode,ctx->nCode);
	assert(ctx->aCode != NULL);
    }
    ctx->aCode[ctx->uIndex] &= ~ctx->uMask;
    ctx->uMask <<= 1;
    if (!ctx->uMask) {
	++ctx->uIndex;
	ctx->uMask = 1;
    }
}

int OutPid(LCODE ctx,uint32_t uPid) {
    int i;
    OutOne(ctx);
    for (i=ctx->nPid-1;i>=0;--i) {
	if (uPid & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    return(1+ctx->nPid);
}


void InPid(LCODE ctx,uint32_t *uPid) {
    int i;
    *uPid = 0;
    for (i=0;i<ctx->nPid;++i) {
	*uPid <<= 1;
	if (InOne(ctx)) *uPid |= 1;
    }
}


int OutPrefix(LCODE ctx,uint32_t uPrefix) {
    int i;
    uPrefix >>= ctx->nSuffix;
    for (i=ctx->nPrefix-1;i>=0;--i) {
	if (uPrefix & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    return(ctx->nPrefix);
}

void InPrefix(LCODE ctx,uint32_t *uPrefix) {
    int i;

    *uPrefix = 0;
    for (i=0;i<ctx->nPrefix;++i) {
	*uPrefix <<= 1;
	if (InOne(ctx)) *uPrefix |= 1;
    }
    *uPrefix <<= ctx->nSuffix;
}


void SkipPrefix(LCODE ctx) {
    int i;

    for (i=0;i<ctx->nPrefix;++i) {
	ctx->uMask <<= 1;
	if (!ctx->uMask) {
	    ++ctx->uIndex;
	    ctx->uMask = 1;
	}
    }
}

int OutRun(LCODE ctx,uint32_t uStart,uint32_t uEnd) {
    int i;
    OutOne(ctx);
    for (i=ctx->nSuffix-1;i>=0;--i) {
	if (uStart & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    for (i=ctx->nSuffix-1;i>=0;--i) {
	if (uEnd & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    return(1+2*ctx->nSuffix);
}

void InRun(LCODE ctx,uint32_t uPrefix,uint32_t *uStart,uint32_t *uEnd) {
    int i;
    *uStart = 0;
    *uEnd = 0;
    for (i=0;i<ctx->nSuffix;++i) {
	*uStart <<= 1;
	if (InOne(ctx)) *uStart |= 1;
    }
    for (i=0;i<ctx->nSuffix;++i) {
	*uEnd <<= 1;
	if (InOne(ctx)) *uEnd |= 1;
    }
    *uStart |= uPrefix;
    *uEnd |= uPrefix;
}

void SkipRun(LCODE ctx) {
    int i;

    for (i=0;i<(2*ctx->nSuffix);++i) {
	ctx->uMask <<= 1;
	if (!ctx->uMask) {
	    ++ctx->uIndex;
	    ctx->uMask = 1;
	}
    }
}

int OutSingle(LCODE ctx,uint32_t uStart) {
    int i;
    OutZero(ctx);
    for (i=ctx->nSuffix-1;i>=0;--i) {
	if (uStart & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    return(1+ctx->nSuffix);
}

void InSingle(LCODE ctx,uint32_t uPrefix,uint32_t *uStart) {
    int i;
    *uStart = 0;
    for (i=0;i<ctx->nSuffix;++i) {
	*uStart <<= 1;
	if (InOne(ctx)) *uStart |= 1;
    }
    *uStart |= uPrefix;
}

void SkipSingle(LCODE ctx) {
    int i;

    for (i=0;i<ctx->nSuffix;++i) {
	ctx->uMask <<= 1;
	if (!ctx->uMask) {
	    ++ctx->uIndex;
	    ctx->uMask = 1;
	}
    }
}

int OutSuffix(LCODE ctx,uint32_t uStart) {
    int i;
    for (i=ctx->nSuffix-1;i>=0;--i) {
	if (uStart & (1<<i)) OutOne(ctx);
	else OutZero(ctx);
    }
    return(ctx->nSuffix);
}

void BackSkipSuffix(LCODE ctx) {
    int i;

    for (i=0;i<ctx->nSuffix;++i) {
	ctx->uMask >>= 1;
	if (!ctx->uMask) {
	    --ctx->uIndex;
	    ctx->uMask = 1<<7;
	}
    }
}


int lcodeEncode(LCODE ctx,LIST *aList,uint32_t nList,char **ppOutput) {
    uint32_t uPrefix,iPid;
    uint32_t nOutBits;
    uint32_t ip,il,isStart,isEnd;
    int nBytes,i;
    /*
    ** Start of encoding. First local members, so skip forward until
    ** iPid in the list first equals ctx->idSelf.
    */
    ctx->uIndex = 0;
    ctx->uMask = 1;
    nOutBits = 0;
    for (ip=0;ip<nList;++ip) if (aList[ip].iPid == ctx->idSelf) break;
    for (il=ip;il<nList;++il) if (aList[il].iPid != ctx->idSelf) break;
    while (ip < il) {
	uPrefix = aList[ip].iIndex & ctx->mPrefix;
	OutOne(ctx);
        ++nOutBits;
	nOutBits += OutPrefix(ctx,uPrefix);
	/*
	** Determine whether there are runs to do.
	*/
	isStart = ip;
	isEnd = isStart+1;
	while (1) {
	    if (isEnd == il) {
		if (isEnd-1 > isStart) {
		    nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		}
		break;
	    }
	    else if ((aList[isEnd].iIndex & ctx->mPrefix) != uPrefix) {
		if (isEnd-1 > isStart) {
		    nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		}
		break;
	    }
	    else if (((aList[isEnd].iIndex-1)&ctx->mSuffix) != (aList[isEnd-1].iIndex&ctx->mSuffix)) {
		if (isEnd-1 > isStart) {
		    nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		}
		isStart = isEnd;
		isEnd = isStart+1;
	    }
	    else {
		++isEnd;
	    }
	}
	OutZero(ctx);
	++nOutBits;
	isStart = ip;
	isEnd = isStart+1;
	while (1) {
	    if (isEnd == il) {
		if (isEnd-1 == isStart) {
		    nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		}
		ip = isEnd;
		break;
	    }
	    else if ((aList[isEnd].iIndex & ctx->mPrefix) != uPrefix) {
		if (isEnd-1 == isStart) {
		    nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		}
		ip = isEnd;
		break;
	    }
	    else if (((aList[isEnd].iIndex-1)&ctx->mSuffix) != (aList[isEnd-1].iIndex&ctx->mSuffix)) {
		if (isEnd-1 == isStart) {
		    nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		}
		isStart = isEnd;
		isEnd = isStart+1;
	    }
	    else {
		++isEnd;
	    }
	}
	OutOne(ctx);
	++nOutBits;
    }
    OutZero(ctx);
    ++nOutBits;
    ip=0;
    il=nList;
    while (ip < il) {
	iPid = aList[ip].iPid;
	if (iPid == ctx->idSelf) {
	    ++ip;
	    continue;
	}
	nOutBits += OutPid(ctx,iPid);
	while (1) {
	    uPrefix = aList[ip].iIndex & ctx->mPrefix;
	    nOutBits += OutPrefix(ctx,uPrefix);
	    /*
	    ** Determine whether there are runs to do.
	    */
	    isStart = ip;
	    isEnd = isStart+1;
	    while (1) {
		if (isEnd == il) {
		    if (isEnd-1 > isStart) {
			nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		    }
		    break;
		}
		else if ((aList[isEnd].iIndex & ctx->mPrefix) != uPrefix || 
			 aList[isEnd].iPid != iPid) {
		    if (isEnd-1 > isStart) {
			nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		    }
		    break;
		}
		else if (((aList[isEnd].iIndex-1)&ctx->mSuffix) != (aList[isEnd-1].iIndex&ctx->mSuffix)) {
		    if (isEnd-1 > isStart) {
			nOutBits += OutRun(ctx,aList[isStart].iIndex,aList[isEnd-1].iIndex);
		    }
		    isStart = isEnd;
		    isEnd = isStart+1;
		}
		else {
		    ++isEnd;
		}
	    }
	    OutZero(ctx);
	    ++nOutBits;
	    isStart = ip;
	    isEnd = isStart+1;
	    while (1) {
		if (isEnd == il) {
		    if (isEnd-1 == isStart) {
			nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		    }
		    ip = isEnd;
		    break;
		}
		else if (aList[isEnd].iPid != iPid) {
		    if (isEnd-1 == isStart) {
			nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		    }
		    il = isEnd; /* we need this to generate a new iPid */
		    ip = isEnd;
		    break;
		}
		else if ((aList[isEnd].iIndex & ctx->mPrefix) != uPrefix) {
		    if (isEnd-1 == isStart) {
			nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		    }
		    ip = isEnd;
		    break;
		}
		else if (((aList[isEnd].iIndex-1)&ctx->mSuffix) != (aList[isEnd-1].iIndex&ctx->mSuffix)) {
		    if (isEnd-1 == isStart) {
			nOutBits += OutSingle(ctx,aList[isStart].iIndex);
		    }
		    isStart = isEnd;
		    isEnd = isStart+1;
		}
		else {
		    ++isEnd;
		}
	    }
	    OutOne(ctx);
	    ++nOutBits;
	    if (ip == il) {
		OutZero(ctx);
		++nOutBits;
		break;
	    }
	    else {
		OutOne(ctx);
		++nOutBits;
	    }
	}
	il = nList;
    }
    OutZero(ctx);
    ++nOutBits;
    printf("nOutBits:%d\n",nOutBits);
    nBytes = (ctx->uMask == 1)?ctx->uIndex:ctx->uIndex+1;
    assert(nBytes*8 >= nOutBits);
    *ppOutput = malloc(nBytes);
    assert(*ppOutput != NULL);
    /*
    ** Now copy the code to the output, but only the number of bytes actually needed.
    */
    for (i=0;i<nBytes;++i) {
	(*ppOutput)[i] = ctx->aCode[i];
    }
    return(nBytes);
}


int lcodeDecode(LCODE ctx,char *pInput,LIST **ppList,int *pnMaxList,int *pnList) {
    int nList = 0;
    int nBytes;
    uint32_t uPrefix;
    uint32_t uStart,uEnd,u,uPid;

    ctx->uIndex = 0;
    ctx->uMask = 1;
    ctx->inCode = pInput;
    while (InOne(ctx)) {
	InPrefix(ctx,&uPrefix);
	while (InOne(ctx)) {
	    InRun(ctx,uPrefix,&uStart,&uEnd);
	    if (nList+(uEnd-uStart+1) > *pnMaxList) {
		*pnMaxList = nList+(uEnd-uStart+1)+LIST_GROW;
		*ppList = realloc(*ppList,(*pnMaxList)*sizeof(LIST));
		assert(*ppList != NULL);
	    }
	    for (u=uStart;u<=uEnd;++u) {
		(*ppList)[nList].iIndex = u;
		(*ppList)[nList].iPid = ctx->idSelf;
		++nList;
	    }
	}
	while (!InOne(ctx)) {
	    InSingle(ctx,uPrefix,&u);
	    if (nList == *pnMaxList) {
		*pnMaxList = nList+1+LIST_GROW;
		*ppList = realloc(*ppList,(*pnMaxList)*sizeof(LIST));
		assert(*ppList != NULL);
	    }
	    (*ppList)[nList].iIndex = u;
	    (*ppList)[nList].iPid = ctx->idSelf;
	    ++nList;
	}
    }
    while (InOne(ctx)) {
	/*
	** Non-local members.
	*/
	InPid(ctx,&uPid);
	do {
	    InPrefix(ctx,&uPrefix);
	    while (InOne(ctx)) {
		InRun(ctx,uPrefix,&uStart,&uEnd);
		if (nList+(uEnd-uStart+1) > *pnMaxList) {
		    *pnMaxList = nList+(uEnd-uStart+1)+LIST_GROW;
		    *ppList = realloc(*ppList,(*pnMaxList)*sizeof(LIST));
		    assert(*ppList != NULL);
		}
		for (u=uStart;u<=uEnd;++u) {
		    (*ppList)[nList].iIndex = u;
		    (*ppList)[nList].iPid = uPid;
		    ++nList;
		}
	    }
	    while (!InOne(ctx)) {
		InSingle(ctx,uPrefix,&u);
		if (nList == *pnMaxList) {
		    *pnMaxList = nList+1+LIST_GROW;
		    *ppList = realloc(*ppList,(*pnMaxList)*sizeof(LIST));
		    assert(*ppList != NULL);
		}
		(*ppList)[nList].iIndex = u;
		(*ppList)[nList].iPid = uPid;
		++nList;
	    }
	} while (InOne(ctx));	
    }
    *pnList = nList;
    nBytes = (ctx->uMask == 1)?ctx->uIndex:ctx->uIndex+1;
    return(nBytes);
}


int bInListLocal(LCODE ctx,char *pInput,uint32_t iIndex) {
    int nList = 0;
    uint32_t uPrefix;
    uint32_t uStart,uEnd,u,uPid;

    ctx->uIndex = 0;
    ctx->uMask = 1;
    ctx->inCode = pInput;
    while (InOne(ctx)) {
	InPrefix(ctx,&uPrefix);
	while (InOne(ctx)) {
	    InRun(ctx,uPrefix,&uStart,&uEnd);
	    if (iIndex >= uStart && iIndex <= uEnd) return(1);
	}
	while (!InOne(ctx)) {
	    InSingle(ctx,uPrefix,&u);
	    if (iIndex == u) return(1);
	}
    }
    return(0);
}
    

#if (0)
/*
** This function is supposed to directly remove elements from a list, but the bookkeepping is 
** a bit of a pain, so will implement it completely later.
*/
int bListRemoveLocal(LCODE ctx,uint32_t iIndex) {
    int nList = 0;
    uint32_t uPrefix;
    uint32_t uStart,uEnd,u[2];
    int nOutSingle;

    ctx->uIndex = 0;
    ctx->uMask = 1;
    while (InOne(ctx)) {
	InPrefix(ctx,&uPrefix);
	while (InOne(ctx)) {
	    InRun(ctx,uPrefix,&uStart,&uEnd);
	    if (iIndex >= uStart && iIndex <= uEnd) {
		/*
		** Found in this run, start copy up to this point.
		*/
		CurrentToCopy(ctx);
		if (uEnd-uStart < 2) {
		    /*
		    ** Run disappears and one single is added under this
		    ** prefix. There were only 2 elements in the run.
		    ** This means that the code needs ctx->nSuffix less
		    ** bits in total.
		    */
		    if (uStart == iIndex) u[0] = uEnd;
		    else u[0] = uStart;
		    nOutSingle = 1;
		    }
		else if (uStart == iIndex) {
		    /*
		    ** Run remains and start is modified.
		    */
		    BackSkipSuffix(ctx);
		    BackSkipSuffix(ctx);
		    OutSuffix(ctx,iIndex);
		    SkipSingle(ctx);
		    nOutSingle = 0;
		}
		else if (uEnd == iIndex) {
		    /*
		    ** Run remains and end is modified.
		    */
		    BackSkipSuffix(ctx);
		    OutSuffix(ctx,iIndex);		    
		    nOutSingle = 0;
		}
		else if (uEnd-uStart == 2) {
		    /*
		    ** Run disappears and 2 singles, uStart and uEnd are
		    ** added under this prefix. We need to make 1-bit of
		    ** extra space for the 2 singles in this case.
		    */
		    u[0] = uStart;
		    u[1] = uEnd;
		    nOutSingle = 2;
		}
		else if (uStart == iIndex-1) {
		    /*
		    ** Run remains with modified start, and a single (uStart)
		    ** is added under this prefix. We need 1+ctx->nSuffix 
		    ** extra bits in this case.
		    */
		    BackSkipSuffix(ctx);
		    BackSkipSuffix(ctx);
		    OutSuffix(ctx,iIndex+1);		    
		    SkipSingle(ctx);
		    u[0] = uStart;
		    nOutSingle = 1;
		}
		else if (uEnd == iIndex+1) {
		    /*
		    ** Run remains with modified uEnd, and a single (uEnd)
		    ** is added under this prefix. We need 1+ctx->nSuffix 
		    ** extra bits in this case.
		    */
		    BackSkipSuffix(ctx);
		    OutSuffix(ctx,iIndex-1);		    
		    u[0] = uEnd;
		    nOutSingle = 1;
		}
		else {
		    assert(uEnd-uStart > 3);
		    /*
		    ** The run is end modified (iIndex-1) and a new run 
		    ** from (iIndex+1,uEnd) is added. We need 1+2*ctx->nSuffix 
		    ** extra bits in this case.
		    */
		    BackSkipSuffix(ctx);
		    OutSuffix(ctx,iIndex-1);
		    InsertRun(ctx,iIndex+1,uEnd);
		    nOutSingle = 0;
		}
		return(1);
	    }
	}
	while (!InOne(ctx)) {
	    InSingle(ctx,uPrefix,&u[0]);
	    if (iIndex == u[0]) {
		/*
		** Single is removed, but if this is the only single under
		** this prefix, then possibly the Prefix must also be 
		** removed (unless there was a run).
		*/
		return(1);
	    }
	}
    }
    return(0);
}
#endif 

int bInList(LCODE ctx,char *pInput,uint32_t iIndex,uint32_t iPid) {
    int nList = 0;
    uint32_t uPrefix;
    uint32_t uStart,uEnd,u,uPid;

    ctx->uIndex = 0;
    ctx->uMask = 1;
    ctx->inCode = pInput;
    if (iPid == ctx->idSelf) {
	while (InOne(ctx)) {
	    InPrefix(ctx,&uPrefix);
	    while (InOne(ctx)) {
		InRun(ctx,uPrefix,&uStart,&uEnd);
		if (iIndex >= uStart && iIndex <= uEnd) return(1);
	    }
	    while (!InOne(ctx)) {
		InSingle(ctx,uPrefix,&u);
		if (iIndex == u) return(1);
	    }
	}
	return(0);
    }
    else {
	while (InOne(ctx)) {
	    SkipPrefix(ctx);
	    while (InOne(ctx)) {
		SkipRun(ctx);
	    }
	    while (!InOne(ctx)) {
		SkipSingle(ctx);
	    }
	}
    }
    while (InOne(ctx)) {
	/*
	** Non-local members.
	*/
	InPid(ctx,&uPid);
	if (uPid == iPid) {
	    do {
		InPrefix(ctx,&uPrefix);
		while (InOne(ctx)) {
		    InRun(ctx,uPrefix,&uStart,&uEnd);
		    if (iIndex >= uStart && iIndex <= uEnd) return(1);
		}
		while (!InOne(ctx)) {
		    InSingle(ctx,uPrefix,&u);
		    if (iIndex == u) return(1);
		}
	    } while (InOne(ctx));
	    return(0);
	}
	else {
	    do {
		SkipPrefix(ctx);
		while (InOne(ctx)) {
		    SkipRun(ctx);
		}
		while (!InOne(ctx)) {
		    SkipSingle(ctx);
		}
	    } while (InOne(ctx));
	}
    }
    return(0);
}

#ifdef LISTCOMP_TEST
int main(void) {
    LCODE ctx;
    LIST aList[1000];
    LIST *bList;
    int nbListMax=0;
    int nList=0;
    uint32_t nLocal;
    uint32_t nSmooth=30;
    uint32_t iPid;
    int i,nOutBytes,nInBytes;
    char *code;
    
    while (scanf("id:%d",&iPid) == 1) {
	while (scanf("%d",&aList[nList].iIndex) == 1) {
	    aList[nList].iPid = iPid;
	    ++nList;
	}
    }
    qsort(aList,nList,sizeof(LIST),lcodeCmpList);

    PrintList(aList,nList);

    nLocal = 10*aList[nList-1].iIndex;
    ctx = lcodeInit(8,4,nLocal,nSmooth);
    nOutBytes = lcodeEncode(ctx,aList,nList,&code);
    printf("\nnOutBytes = %d = %f bytes/element = %5.2f%%\n",
	   nOutBytes,1.0*nOutBytes/nList,100.0*nOutBytes/4/nList);

    /*
    ** Start of decoding.
    */
    nInBytes = lcodeDecode(ctx,code,&bList,&nbListMax,&nList);
    assert(nInBytes == nOutBytes);

    PrintList(bList,nList);

    for (i=0;i<nList;++i) {
	if (bInList(ctx,code,aList[i].iIndex,aList[i].iPid)) {
	    printf("%d:%d found\n",aList[i].iPid,aList[i].iIndex);
	}
	else {
	    printf("%d:%d NOT FOUND\n",aList[i].iPid,aList[i].iIndex);
	}
    }
}
#endif
