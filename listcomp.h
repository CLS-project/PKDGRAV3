#ifndef LISTCOMP_INCLUDED
#define LISTCOMP_INCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <stdint.h>

typedef struct listElement {
    uint32_t iIndex;
    uint32_t iPid;
} LIST;

typedef struct lcodeContext {
    uint32_t nPid;
    uint32_t idSelf;
    uint32_t nPrefix;
    uint32_t nSuffix;
    uint32_t mPid;
    uint32_t mPrefix;
    uint32_t mSuffix;
/*
** the following variables are used for creating the encoded list
*/
    uint32_t uIndex;
    unsigned char uMask;
    uint32_t nCode;
    char *aCode;
    char *inCode; /* we set this pointer to the input string to decode */
} * LCODE;


static inline int InOne(LCODE ctx) {
    int iRet;
    iRet = (ctx->inCode[ctx->uIndex] & ctx->uMask);
    ctx->uMask <<= 1;
    if (!ctx->uMask) {
	++ctx->uIndex;
	ctx->uMask = 1;
    }
    return(iRet);
}

static inline void InPrefix(LCODE ctx,uint32_t *uPrefix) {
    int i;

    *uPrefix = 0;
    for (i=0;i<ctx->nPrefix;++i) {
	*uPrefix <<= 1;
	if (InOne(ctx)) *uPrefix |= 1;
    }
    *uPrefix <<= ctx->nSuffix;
}

static inline void InRun(LCODE ctx,uint32_t uPrefix,uint32_t *uStart,uint32_t *uEnd) {
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

static inline void InSingle(LCODE ctx,uint32_t uPrefix,uint32_t *uStart) {
    int i;
    *uStart = 0;
    for (i=0;i<ctx->nSuffix;++i) {
	*uStart <<= 1;
	if (InOne(ctx)) *uStart |= 1;
    }
    *uStart |= uPrefix;
}

static inline int bInListLocal(LCODE ctx,char *pInput,uint32_t iIndex) {
    uint32_t uPrefix;
    uint32_t uStart,uEnd,u;

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
    

void lcodePrintList(LIST *p,int nList);
int lcodeCmpList(const void *v1,const void *v2);
LCODE lcodeInit(uint32_t nThreads,uint32_t idSelf,uint32_t nLocal,uint32_t nSmooth);
void lcodeFinish(LCODE ctx);
int lcodeEncode(LCODE ctx,LIST *aList,uint32_t nList,char **ppOutput);
int lcodeDecode(LCODE ctx,char *pInput,LIST **ppList,int *pnMaxList,int *pnList);
int bInList(LCODE ctx,char *pInput,uint32_t iIndex,uint32_t iPid);

#endif
