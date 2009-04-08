#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include "pkd.h"
#include "outtype.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#ifdef LARGEF
#include <fcntl.h>
#endif

const char *outtype_c_module_id = "$Id$";
const char *outtype_h_module_id = OUTTYPE_H_MODULE_ID;

#define OUTTYPE_UNKNOWN 0
#define OUTTYPE_INTEGER 1
#define OUTTYPE_FLOAT 2

static int getType(int iType) {
    switch(iType) {
    case OUT_IORDER_ARRAY:
    case OUT_GROUP_ARRAY:
	return OUTTYPE_INTEGER;

    case OUT_DENSITY_ARRAY:
    case OUT_COLOR_ARRAY:
    case OUT_POT_ARRAY:
    case OUT_AMAG_ARRAY:
    case OUT_RUNG_ARRAY:
    case OUT_SOFT_ARRAY:
    case OUT_RELAX_ARRAY:
    case OUT_DIVV_ARRAY:
    case OUT_VELDISP2_ARRAY:
    case OUT_VELDISP_ARRAY:
    case OUT_PHASEDENS_ARRAY:

    case OUT_POS_VECTOR:
    case OUT_VEL_VECTOR:
    case OUT_MEANVEL_VECTOR:
    case OUT_ACCEL_VECTOR:
	return OUTTYPE_FLOAT;

    default:
	return OUTTYPE_UNKNOWN;
	}
    }

/* Write an integer */
static uint64_t fetchInteger(PKD pkd,PARTICLE *p,int iType,int iDim) {
    uint64_t v;

    switch (iType) {
    case OUT_IORDER_ARRAY:
	v = p->iOrder;
	break;
    case OUT_GROUP_ARRAY:
	assert(pkd->oGroup);
	v = *pkdInt32(p,pkd->oGroup);
	break;
    default:
	v = 0;
	}
    return v;
    }
static double fetchFloat(PKD pkd,PARTICLE *p,int iType,int iDim) {
    float *a;
    double v;
    VELSMOOTH *pvel;
    switch (iType) {
    case OUT_DENSITY_ARRAY:
	v = p->fDensity;
	break;
    case OUT_COLOR_ARRAY:
	assert(0);
    case OUT_POT_ARRAY:
	assert(pkd->oPotential);
	a = pkdPot(pkd,p);
	v = *a;
	break;
    case OUT_AMAG_ARRAY:
	assert(pkd->oAcceleration); /* Validate memory model */
	a = pkdAccel(pkd,p);
	v = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	break;
    case OUT_RUNG_ARRAY:
	v = p->uRung;
	break;
    case OUT_SOFT_ARRAY:
	v = pkdSoft(pkd,p);
	break;
    case OUT_RELAX_ARRAY:
	assert(pkd->oRelaxation);
	a = pkdField(p,pkd->oRelaxation);
	v = *a;
	break;
    case OUT_DIVV_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->divv;
	break;
    case OUT_VELDISP2_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->veldisp2;
    case OUT_VELDISP_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = sqrt(pvel->veldisp2);
	break;
    case OUT_PHASEDENS_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = p->fDensity*pow(pvel->veldisp2,-1.5);
	break;
    case OUT_POS_VECTOR:
	v = p->r[iDim];
	break;
    case OUT_VEL_VECTOR:
	assert(pkd->oVelocity); /* Validate memory model */
	v = pkdVel(pkd,p)[iDim];
	break;
    case OUT_MEANVEL_VECTOR:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->vmean[iDim];
	break;
    case OUT_ACCEL_VECTOR:
	assert(pkd->oAcceleration); /* Validate memory model */
	a = pkdAccel(pkd,p);
	v = a[iDim];
	break;
    default:
	v = 0.0;
	}
    return v;
    }

/******************************************************************************\
 * Generic buffered output - flush needs to be customized for the type of
 * output stream (ASCII,BZIP2,GZIP,etc.).
\******************************************************************************/
static void storeInteger(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int n = ctx->inOffset - ctx->inBuffer;
    if ( PKDOUT_BUFFER_SIZE - n < 20 ) {
	(*ctx->fnFlush)(pkd,ctx,0);
	}
    sprintf(ctx->inOffset,"%"PRIu64"\n",fetchInteger(pkd,p,iType,iDim));
    assert(strlen(ctx->inOffset) < 20 );
    while( *ctx->inOffset ) ++ctx->inOffset;
    }
static void storeFloat(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    if ( PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer) < 20 )
	(*ctx->fnFlush)(pkd,ctx,0);
    sprintf(ctx->inOffset,"%.8g\n",fetchFloat(pkd,p,iType,iDim));
    assert(strlen(ctx->inOffset) < 20 );
    while( *ctx->inOffset ) ++ctx->inOffset;
    }
static void storeHdr(PKD pkd,PKDOUT ctx,uint64_t N) {
    if ( PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer) < 20 )
	(*ctx->fnFlush)(pkd,ctx,0);
    sprintf(ctx->inOffset,"%"PRIu64"\n",N);
    while( *ctx->inOffset ) ++ctx->inOffset;
    }
static void finish(PKD pkd,PKDOUT ctx) {
    (*ctx->fnFlush)(pkd,ctx,0); /* Flush input buffer */
    (*ctx->fnFlush)(pkd,ctx,1); /* Finish output stream */
    }

/******************************************************************************\
 * Regular ASCII
\******************************************************************************/

/*
** Flush the input buffer by writing it directly to a file
*/
static void writeASCII(PKD pkd,PKDOUT ctx,int final) {
    int n = ctx->inOffset - ctx->inBuffer;
    fwrite(ctx->inBuffer,n,1,ctx->fp);
    ctx->inOffset = ctx->inBuffer;
    }

static void closeASCII(PKD pkd,PKDOUT ctx) {
    (*ctx->fnFlush)(pkd,ctx,0); /* Flush input buffer */
    }

/******************************************************************************\
 * bzip2 compression
\******************************************************************************/

#ifdef HAVE_LIBBZ2
static void closeBZ2(PKD pkd,PKDOUT ctx) {
    (*ctx->fnFlush)(pkd,ctx,0); /* Flush input buffer */
    (*ctx->fnFlush)(pkd,ctx,1); /* Flush output buffer */
    }
/*
** Write routine to write to an output file
*/
static void outputBZ2(PKD pkd,PKDOUT ctx) {
    int n = ctx->CTX.bzStream->next_out - ctx->outBuffer->data;
    fwrite(ctx->outBuffer->data,n,1,ctx->fp);
    ctx->CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
    ctx->CTX.bzStream->next_out = ctx->outBuffer->data;
    }
/*
** Write routine to allocate a new buffer.
*/
static void bufferBZ2(PKD pkd,PKDOUT ctx) {
    ctx->outBuffer->nBytes = PKDOUT_BUFFER_SIZE - ctx->CTX.bzStream->avail_out;
    if ( ctx->CTX.bzStream->avail_out == 0 ) {
	ctx->outBuffer = ctx->outBuffer->next = malloc(sizeof(PKDOUTBUFFER));
	ctx->outBuffer->next = NULL;
	ctx->outBuffer->nBytes = 0;
	ctx->CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
	ctx->CTX.bzStream->next_out = ctx->outBuffer->data;
	}
    }
/*
** Called when the input buffer is (nearly) full.  Empties the input buffer by
** compressing the data.  New output buffers are allocated as required.
*/
static void flushBZ2(PKD pkd,PKDOUT ctx,int final) {
    int bzerror;
    ctx->CTX.bzStream->next_in = ctx->inBuffer;
    ctx->CTX.bzStream->avail_in = ctx->inOffset - ctx->inBuffer;
    ctx->inOffset = ctx->inBuffer;
    while(ctx->CTX.bzStream->avail_in) {
	if (ctx->CTX.bzStream->avail_out==0) {
	    ctx->nBytes += (ctx->outBuffer->nBytes=PKDOUT_BUFFER_SIZE);
	    (*ctx->fnWrite)(pkd,ctx);
	    }
	bzerror = BZ2_bzCompress(ctx->CTX.bzStream,BZ_RUN);
	assert(bzerror>=0);
	}
    if (final) {
	do {
	    bzerror = BZ2_bzCompress(ctx->CTX.bzStream,BZ_FINISH);
	    (*ctx->fnWrite)(pkd,ctx);
	    } while(bzerror>0&&bzerror!=BZ_STREAM_END);
	(*ctx->fnWrite)(pkd,ctx);
	ctx->outBuffer->nBytes = PKDOUT_BUFFER_SIZE - ctx->CTX.bzStream->avail_out;
	ctx->nBytes += ctx->outBuffer->nBytes;
	BZ2_bzCompressEnd(ctx->CTX.bzStream);
	ctx->CTX.bzStream = NULL;
	}
    }

static void setupBZ2(PKD pkd,PKDOUT ctx) {
    ctx->fnClose = closeBZ2;
    ctx->fnFlush = flushBZ2;
    ctx->fnWrite = NULL; /* Still must be set */

    ctx->CTX.bzStream = malloc(sizeof(bz_stream));
    assert(ctx->CTX.bzStream!=NULL);
    ctx->CTX.bzStream->bzalloc = NULL;
    ctx->CTX.bzStream->bzfree = NULL;
    ctx->CTX.bzStream->opaque = NULL;
    ctx->CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
    ctx->CTX.bzStream->next_out = ctx->outBuffer->data;
    BZ2_bzCompressInit(ctx->CTX.bzStream,9,0,0);
    }
#endif

/******************************************************************************\
 * zlib compression
\******************************************************************************/

#ifdef HAVE_LIBZ
/*
** Flush the input buffer by writing it directly to a file
*/
static void closeZ(PKD pkd,PKDOUT ctx) {
    (*ctx->fnFlush)(pkd,ctx,0); /* Flush input buffer */
    (*ctx->fnFlush)(pkd,ctx,1); /* Flush output buffer */
    }
/*
** Write routine to write to an output file
*/
static void outputZ(PKD pkd,PKDOUT ctx) {
    int n = (char *)ctx->CTX.gzStream->next_out - ctx->outBuffer->data;
    fwrite(ctx->outBuffer->data,n,1,ctx->fp);
    ctx->CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
    ctx->CTX.gzStream->next_out = ctx->outBuffer->data;
    }
/*
** Write routine to allocate a new buffer.
*/
static void bufferZ(PKD pkd,PKDOUT ctx) {
    ctx->outBuffer->nBytes = PKDOUT_BUFFER_SIZE - ctx->CTX.gzStream->avail_out;
    if ( ctx->CTX.gzStream->avail_out == 0 ) {
	ctx->outBuffer = ctx->outBuffer->next = malloc(sizeof(PKDOUTBUFFER));
	ctx->outBuffer->next = NULL;
	ctx->outBuffer->nBytes = 0;
	ctx->CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
	ctx->CTX.gzStream->next_out = ctx->outBuffer->data;
	}
    }
/*
** Called when the input buffer is (nearly) full.  Empties the input buffer by
** compressing the data.  New output buffers are allocated as required.
*/
static void flushZ(PKD pkd,PKDOUT ctx,int final) {
    int gzerror;
    ctx->CTX.gzStream->next_in = ctx->inBuffer;
    ctx->CTX.gzStream->avail_in = ctx->inOffset - ctx->inBuffer;
    ctx->inOffset = ctx->inBuffer;

    while(ctx->CTX.gzStream->avail_in) {
	if (ctx->CTX.gzStream->avail_out==0) {
	    ctx->nBytes += (ctx->outBuffer->nBytes=PKDOUT_BUFFER_SIZE);
	    (*ctx->fnWrite)(pkd,ctx);
	    }
	gzerror = deflate(ctx->CTX.gzStream,0);
	assert(gzerror>=0);
	}
    if (final) {
	do {
	    gzerror = deflate(ctx->CTX.gzStream,Z_FINISH);
	    assert(gzerror>=0);
	    (*ctx->fnWrite)(pkd,ctx);
	    } while(gzerror!=Z_STREAM_END);
	ctx->outBuffer->nBytes = PKDOUT_BUFFER_SIZE - ctx->CTX.gzStream->avail_out;
	ctx->nBytes += ctx->outBuffer->nBytes;
	deflateEnd(ctx->CTX.gzStream);
	ctx->CTX.gzStream = NULL;
	}
    }

void setupZ(PKD pkd,PKDOUT ctx) {
    ctx->fnClose = closeZ;
    ctx->fnFlush = flushZ;
    ctx->fnWrite = NULL;

    ctx->CTX.gzStream = malloc(sizeof(z_stream));
    assert(ctx->CTX.gzStream!=NULL);
    ctx->CTX.gzStream->zalloc = Z_NULL;
    ctx->CTX.gzStream->zfree = Z_NULL;
    ctx->CTX.gzStream->opaque = Z_NULL;
    ctx->CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
    ctx->CTX.gzStream->next_out = ctx->outBuffer->data;
    deflateInit2(ctx->CTX.gzStream,Z_BEST_COMPRESSION,Z_DEFLATED,31,8,Z_DEFAULT_STRATEGY);
    }
#endif

/******************************************************************************\
 * Storage routines
\******************************************************************************/

PKDOUT pkdStartOutASCII(PKD pkd,int iFile, int iType) {
    PKDOUT ctx;

    /*
    ** Allocate the context, input buffer, and the first output buffer.
    */
    ctx = malloc(sizeof(struct pkdout)); assert(ctx!=NULL);
    ctx->fp = NULL;
    ctx->outBuffer = ctx->headBuffer = malloc(sizeof(PKDOUTBUFFER));
    assert(ctx->outBuffer!=NULL);
    ctx->outBuffer->next = NULL;
    ctx->inBuffer = ctx->inOffset = malloc(PKDOUT_BUFFER_SIZE);
    assert(ctx->inBuffer!=NULL);
    ctx->nBytes = 0;

    switch(getType(iType)) {
    case OUTTYPE_INTEGER:
	ctx->fnOut = storeInteger;
	break;
    case OUTTYPE_FLOAT:
	ctx->fnOut = storeFloat;
	break;
    default:
	assert(0);
	}
    ctx->fnHdr = storeHdr;
    ctx->fnClose = finish;

    switch(iFile) {
#ifdef HAVE_LIBBZ2
    case PKDOUT_TYPE_BZIP2:
	setupBZ2(pkd,ctx);
	ctx->fnWrite = bufferBZ2;
	break;
#endif
#ifdef HAVE_LIBZ
    case PKDOUT_TYPE_ZLIB:
	setupZ(pkd,ctx);
	ctx->fnWrite = bufferZ;
	break;
#endif
    default:
	assert(0);
	}
    return ctx;
    }

/*
** Finish the output stream by flushing any remaining characters in the output
** stream, and then finalizing it.  Compressed data is still available in the
** output buffers
*/
void pkdFinishOutASCII(PKD pkd,PKDOUT ctx) {
    (*ctx->fnClose)(pkd,ctx);
    }

/*
** Returns the total number of bytes available in the output buffers.  This will
** only return a valid number after Finish has been called.
*/
uint64_t pkdCountOutASCII(PKD pkd,PKDOUT ctx) {
    return ctx->nBytes;
    }

void pkdDumpOutASCII(PKD pkd,PKDOUT ctx,FILE *fp) {
    PKDOUTBUFFER *buf;
    for(buf=ctx->headBuffer;buf!=NULL;buf=buf->next) {
	fwrite(buf->data,buf->nBytes,1,fp);
	}
    }

/*
** Free up all buffers, and the PKDOUT context
*/
void pkdFreeOutASCII(PKD pkd,PKDOUT ctx) {
    PKDOUTBUFFER *buf, *nxt;

    free(ctx->inBuffer);
    for(buf=ctx->headBuffer;buf!=NULL;buf=nxt) {
	nxt = buf->next;
	free(buf);
	}
    free(ctx);
    }

/******************************************************************************\
 * File I/O routines
\******************************************************************************/

PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iFile,int iType) {
    PKDOUT ctx;

    ctx = malloc(sizeof(struct pkdout)); assert(ctx!=NULL);
    ctx->fp = NULL;
    ctx->outBuffer = ctx->headBuffer = malloc(sizeof(PKDOUTBUFFER));
    assert(ctx->outBuffer!=NULL);
    ctx->outBuffer->next = NULL;
    ctx->inBuffer = ctx->inOffset = malloc(PKDOUT_BUFFER_SIZE);
    assert(ctx->inBuffer!=NULL);
    ctx->nBytes = 0;

    /* Determine how to handle the header and data */
    switch(getType(iType)) {
    case OUTTYPE_INTEGER:
	ctx->fnOut = storeInteger;
	break;
    case OUTTYPE_FLOAT:
	ctx->fnOut = storeFloat;
	break;
    default:
	assert(0);
	}
    ctx->fnHdr = storeHdr;

    ctx->fp = fopen (pszFileName,mode);
    assert(ctx->fp != NULL);
    setvbuf(ctx->fp,NULL,_IOFBF,PKDOUT_BUFFER_SIZE);

    switch(iFile) {

#ifdef HAVE_LIBBZ2
    case PKDOUT_TYPE_BZIP2:
	setupBZ2(pkd,ctx);
	ctx->fnWrite = outputBZ2;
	break;
#endif
#ifdef HAVE_LIBZ
    case PKDOUT_TYPE_ZLIB:
	setupZ(pkd,ctx);
	ctx->fnWrite = outputZ;
	break;
#endif
    default:
	ctx->fnFlush = writeASCII;
	ctx->fnClose = closeASCII;
	}
    return ctx;
    }

void pkdCloseOutASCII(PKD pkd,PKDOUT ctx) {
    (*ctx->fnClose)(pkd,ctx);
    fclose(ctx->fp);
    free(ctx->inBuffer);
    free(ctx);
    }

void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N) {
    (*ctx->fnHdr)(pkd,ctx,N);
    }

void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim) {
    int i;

    /*
    ** Write Elements!
    */
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p = pkdParticle(pkd,i);
	if ( pkdIsSrcActive(p,0,MAX_RUNG) )
	    (*ctx->fnOut)(pkd,ctx,p,iType,iDim);
	}
    }

#ifdef USE_HDF5
void pkdOutHDF5(PKD pkd,char *pszFileName,int iType,int iDim) {
    assert(0);
    }
#endif

void pkdOutGroup(PKD pkd,char *pszFileName,int iType, int nStart,double dvFac) {
    FILE *fp;
    int i,j,nout,lStart;

    /*
     ** Write Group Data!
     */
    if (iType == OUT_GROUP_STATS) {
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	for (i=0;i<pkd->nGroups;++i) {
	    fprintf(fp,"%d ",pkd->groupData[i].iGlobalId);
	    fprintf(fp,"%d ",pkd->groupData[i].nTotal);
	    fprintf(fp,"%.8g ",pkd->groupData[i].fMass);
            fprintf(fp,"%.8g ",pkd->groupData[i].fRMSRadius);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[0]);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[1]);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[2]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[0]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[1]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[2]);
	    fprintf(fp,"\n");
	    }
	}
    else if (iType == OUT_GROUP_TIPSY_NAT) {

	struct star_particle sp;
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	/*
	 ** Seek past the header
	 */
	lStart = sizeof(struct dump);
	fseek(fp,lStart,SEEK_SET);

	for (i=0;i<pkd->nGroups;++i) {
	    if (pkd->groupData[i].bMyGroup) {
		for (j=0;j<3;++j) {
		  sp.pos[j] = pkd->groupData[i].r[j];
		  sp.vel[j] = dvFac*pkd->groupData[i].v[j];
		}
		sp.mass = pkd->groupData[i].fMass;
		sp.eps = pkd->groupData[i].fRMSRadius;
		sp.tform = 0.0;
		sp.metals = 0.0;
		nout = fwrite(&sp,sizeof(struct star_particle),1,fp);
		mdlassert(pkd->mdl,nout == 1);
		}
	    }
	}
    else if (iType == OUT_GROUP_TIPSY_STD) {

	float fTmp;
	XDR xdrs;

	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	/*
	 ** Seek past the header
	 */
	lStart = 32;
	lStart += nStart*44;

	fseek(fp,lStart,SEEK_SET);

	xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	for (i=0;i<pkd->nGroups;++i) {
	    if (pkd->groupData[i].bMyGroup) {
		fTmp = pkd->groupData[i].fMass;
		xdr_float(&xdrs,&fTmp);
		for (j=0;j<3;++j) {
		    fTmp = pkd->groupData[i].r[j];
		    xdr_float(&xdrs,&fTmp);
		    }
		for (j=0;j<3;++j) {
		    fTmp = dvFac*pkd->groupData[i].v[j];
		    xdr_float(&xdrs,&fTmp);
		    }
		fTmp = 0.0;
		xdr_float(&xdrs,&fTmp); /* metals */
		xdr_float(&xdrs,&fTmp); /* t form */
		fTmp = pkd->groupData[i].fRMSRadius;
		xdr_float(&xdrs,&fTmp); /* softening eps*/
		fTmp = 0.0;
		xdr_float(&xdrs,&fTmp); /* grav. potential phi*/
		}
	    }
	xdr_destroy(&xdrs);

	}
    else if (iType == OUT_GROUP_PROFILES) {
	fp = fopen(pszFileName,"at");
	assert(fp != NULL);
	for (i=0;i< pkd->nBins;++i) {
	    fprintf(fp,"%.8g ",pkd->groupBin[i].fRadius);
	    fprintf(fp,"%d ",pkd->groupBin[i].nMembers);
	    fprintf(fp,"%.8g ",pkd->groupBin[i].fMassInBin);
	    fprintf(fp,"\n");
	    }
	}
    else assert(0);
    i = fclose(fp);
    if (i != 0) {
	perror("pkdOutGroup: could not close file");
	exit(1);
	}
    }
