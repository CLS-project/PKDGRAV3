#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "output.h"


/* Generic context: all things/stuff from iIndex (starts at zero) */
struct packCtx {
    PKD pkd;
    int iIndex;
    };

static int packGroupStats(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packCtx *ctx = (struct packCtx *)vctx;
    int nLeft = ctx->pkd->nLocalGroups - ctx->iIndex;
    int n = nSize / sizeof(TinyGroupTable);
    if ( n > nLeft ) n = nLeft;
    memcpy(vBuff,ctx->pkd->tinyGroupTable + 1 + ctx->iIndex, n*sizeof(TinyGroupTable) );
    ctx->iIndex += n;
    return n*sizeof(TinyGroupTable);
    }

static int unpackGroupStats(void *vctx, int *id, size_t nSize, void *vBuff) {
    asyncFileInfo *info = vctx;
    io_write(info,vBuff,nSize);
    return 1;
    }

void pkdOutput(PKD pkd, int eOutputType, int iProcessor,int nProcessor,
    int iPartner,int nPartner, const char *fname ) {
    struct packCtx ctx = {pkd,0};


    /* I do all of the writing */
    if (iPartner == pkd->idSelf) {
	asyncFileInfo info;
	char achOutFile[256];
	strcpy(achOutFile,fname);
	sprintf(achOutFile+strlen(achOutFile),".%d",iProcessor);
	io_init(&info);
	if (io_create(&info,achOutFile) < 0) { perror(fname); abort(); }
	io_write(&info,pkd->tinyGroupTable+1,sizeof(TinyGroupTable)*pkd->nLocalGroups);
	while(--nPartner) {
	    mdlRecv(pkd->mdl,++iPartner,unpackGroupStats,&info);
	    }
	io_close(&info);
	io_free(&info);
	}
    /* We just send all of our data onward */
    else {
	struct packCtx ctx;
	ctx.pkd = pkd;
	ctx.iIndex = 0;
	mdlSend(pkd->mdl,iPartner,packGroupStats, &ctx);
	}
    }

