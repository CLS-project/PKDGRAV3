#ifndef LISTCOMP_INCLUDED
#define LISTCOMP_INCLUDED

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
    char uMask;
    uint32_t nCode;
    char *aCode;
    char *inCode; /* we set this pointer to the input string to decode */
} * LCODE;


int lcodeCmpList(const void *v1,const void *v2);
LCODE lcodeInit(uint32_t nThreads,uint32_t idSelf,uint32_t nLocal,uint32_t nSmooth);
void lcodeFinish(LCODE ctx);
int lcodeEncode(LCODE ctx,LIST *aList,uint32_t nList,char **ppOutput);
int lcodeDecode(LCODE ctx,char *pInput,LIST **ppList,int *pnMaxList,int *pnList);
int bInListLocal(LCODE ctx,char *pInput,uint32_t iIndex);
int bInList(LCODE ctx,char *pInput,uint32_t iIndex,uint32_t iPid);

#endif
