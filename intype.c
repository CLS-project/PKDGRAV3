#include "intype.h"

static void *openMark(PKD pkd,FILE *fp) {
    uint64_t N, nDark, nGas, nStar, i;
    fscanf(fp,"%"PRIu64" %"PRIu64" %"PRIu64, &nDark, &nGas, &nStar);
    N = nDark + nGas + nStar;
    }

static int readMark(void *ctx,FILE *fp,uint64_t i,int iType,int iDim) {
    uint64_t id;
    return fscanf(fp,"%"PRIu64, &id );
    }

static void doneMark(void *ctx) {
    }

void pkdInASCII(PKD pkd,char *pszFileName,int iType,int iDim) {
    FILE *fp;
    uint64_t N, nDark, nGas, nStar, i;
    void *ctx;
    void * (*fnOpen)(PKD pkd,FILE *fp);
    int (*fnRead)(void *ctx,FILE *fp,uint64_t i,int iType,int iDim);
    void (*fnDone)(void *ctx);

    switch(iType) {
    case IN_SRC_MARK:
    case IN_DST_MARK:
	fnOpen = openMark;
	fnRead = readMark;
	fnDone = doneMark;
	break;

    default:
	assert(0);
	}

    fp = fopen (pszFileName,"r");
    assert(fp != NULL);

    ctx = (*fnOpen)(pkd,fp);
    while ((*fnRead)(ctx,fp,i,iType,iDim)>0) {}
    (*fnDone)(ctx);


    if (fclose(fp) != 0) {
	perror("pkdInASCII: could not close file");
	exit(1);
	}
    }
