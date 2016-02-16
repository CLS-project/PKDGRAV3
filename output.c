#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "output.h"
#include "pst.h"

/* Generic context: all things/stuff from iIndex (starts at zero) */
struct packCtx {
    PKD pkd;
    int iIndex;
    };

/*
** Tiny Group statistics
*/
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




/*
** We do not do the write, rather we send to another thread.
*/

void pkdOutputSend(PKD pkd, outType eOutputType, int iPartner) {
    struct packCtx ctx;
    ctx.pkd = pkd;
    ctx.iIndex = 0;
    switch(eOutputType) {
    case OUT_TINY_GROUP:
	mdlSend(pkd->mdl,iPartner, packGroupStats, &ctx);
	break;
    default:
	fprintf(stderr,"ERROR: invalid output type %d\n", eOutputType);
	abort();
	}
    }

/*
** We are the writer. We may need to receive as well.
*/

void pkdOutput(PKD pkd, outType eOutputType, int iProcessor,int nProcessor,
    int iPartner,int nPartner, const char *fname ) {
    struct packCtx ctx = {pkd,0};
    mdlPack unpack;

    /* I do all of the writing */
    if (iPartner == pkd->idSelf) {
	asyncFileInfo info;
	char achOutFile[256];
	strcpy(achOutFile,fname);
	sprintf(achOutFile+strlen(achOutFile),".%d",iProcessor);
	io_init(&info);
	if (io_create(&info,achOutFile) < 0) { perror(fname); abort(); }

	switch(eOutputType) {
	case OUT_TINY_GROUP:
	    io_write(&info,pkd->tinyGroupTable+1,sizeof(TinyGroupTable)*pkd->nLocalGroups);
	    unpack = unpackGroupStats;
	    break;
	default:
	    unpack = NULL;
	    fprintf(stderr,"ERROR: invalid output type %d\n", eOutputType);
	    abort();
	    }
	while(--nPartner) {
	    struct inOutputSend send;
	    send.iPartner = pkd->idSelf;
	    send.eOutputType = eOutputType;
	    ++iPartner;
	    int rID = mdlReqService(pkd->mdl,iPartner,PST_OUTPUT_SEND,&send,sizeof(send));
	    mdlRecv(pkd->mdl,iPartner,unpack,&info);
	    mdlGetReply(pkd->mdl,rID,NULL,NULL);
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

