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

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cinttypes>
#include "pkd.h"
#include "outtype.h"
#ifdef LARGEF
    #include <fcntl.h>
#endif

#define OUTTYPE_UNKNOWN 0
#define OUTTYPE_INTEGER 1
#define OUTTYPE_FLOAT 2
#define OUTTYPE_RUNGDEST 3
#define OUTTYPE_PSGROUP 4

static int getType(int iType) {
    switch (iType) {
    case OUT_IORDER_ARRAY:
    case OUT_GROUP_ARRAY:
    case OUT_GLOBALGID_ARRAY:
    case OUT_MARKED_ARRAY:
        return OUTTYPE_INTEGER;
    case OUT_PSGROUP_ARRAY:
        return OUTTYPE_PSGROUP;

    case OUT_BALL_ARRAY:
    case OUT_DENSITY_ARRAY:
    case OUT_POT_ARRAY:
    case OUT_AMAG_ARRAY:
    case OUT_RUNG_ARRAY:
    case OUT_SOFT_ARRAY:
    case OUT_DIVV_ARRAY:
    case OUT_VELDISP2_ARRAY:
    case OUT_VELDISP_ARRAY:
    case OUT_PHASEDENS_ARRAY:
    case OUT_C_ARRAY:
    case OUT_HSPH_ARRAY:
    case OUT_POS_VECTOR:
    case OUT_VEL_VECTOR:
    case OUT_MEANVEL_VECTOR:
    case OUT_ACCEL_VECTOR:
        return OUTTYPE_FLOAT;

    case OUT_CACHEFLUX_ARRAY:
        return OUTTYPE_INTEGER;
    case OUT_CACHECOLL_ARRAY:
        return OUTTYPE_INTEGER;
    case OUT_AVOIDEDFLUXES_ARRAY:
        return OUTTYPE_INTEGER;
    case OUT_COMPUTEDFLUXES_ARRAY:
        return OUTTYPE_INTEGER;
    case OUT_RUNGDEST_ARRAY:
        return OUTTYPE_RUNGDEST;

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
        v = pkdGetGroup(pkd,p);
        break;
    case OUT_GLOBALGID_ARRAY:
        v = pkdGetGlobalGid(pkd,p);
        break;
    case OUT_MARKED_ARRAY:
        v = p->bMarked;
        break;
    case OUT_PSGROUP_ARRAY:
        assert(0);
        /*v = pkd->psGroupData[pkdGetGroup(pkd,p)].iGlobalId;*/
        v = pkdGetGroup(pkd,p);
        break;
#ifdef DEBUG_CACHED_FLUXES
    case OUT_CACHEFLUX_ARRAY:
        if (pkdIsGas(pkd,p)) {
            v = pkdSph(pkd,p)->flux_cache;
        }
        else {
            v = 0;
        }
        break;
    case OUT_CACHECOLL_ARRAY:
        if (pkdIsGas(pkd,p)) {
            v = pkdSph(pkd,p)->coll_cache;
        }
        else {
            v = 0;
        }
        break;
    case OUT_AVOIDEDFLUXES_ARRAY:
        if (pkdIsGas(pkd,p)) {
            v = (uint64_t)(pkdSph(pkd,p)->avoided_fluxes);
        }
        else {
            v = 0;
        }
        break;
    case OUT_COMPUTEDFLUXES_ARRAY:
        if (pkdIsGas(pkd,p)) {
            v = (uint64_t)(pkdSph(pkd,p)->computed_fluxes);
        }
        else {
            v = 0;
        }
        break;
#endif
    default:
        v = 0;
    }
    return v;
}
static double fetchFloat(PKD pkd,PARTICLE *p,int iType,int iDim) {
    double v;
//    VELSMOOTH *pvel;
    switch (iType) {
    case OUT_DENSITY_ARRAY:
        v = pkd->particles.density(p);
        break;
    case OUT_BALL_ARRAY:
        v = pkd->particles.ball(p);
        break;
    case OUT_POT_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oPotential));
        v = pkd->particles.potential(p);
        break;
    case OUT_AMAG_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oAcceleration)); /* Validate memory model */
        v = sqrt(blitz::dot(pkd->particles.acceleration(p),pkd->particles.acceleration(p)));
        break;
    case OUT_RUNG_ARRAY:
        v = p->uRung;
        break;
    case OUT_SOFT_ARRAY:
        v = pkd->particles.soft(p);
        break;
    case OUT_DIVV_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oVelSmooth)); /* Validate memory model */
        v = pkd->particles.get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth).divv;
        break;
    case OUT_VELDISP2_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oVelSmooth)); /* Validate memory model */
        v = pkd->particles.get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth).veldisp2;
        break;
    case OUT_VELDISP_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oVelSmooth)); /* Validate memory model */
        v = sqrt(pkd->particles.get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth).veldisp2);
        break;
    case OUT_PHASEDENS_ARRAY:
        assert(pkd->particles.present(PKD_FIELD::oVelSmooth)); /* Validate memory model */
        v = pkdDensity(pkd,p)*pow(pkd->particles.get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth).veldisp2,-1.5);
        break;
    case OUT_C_ARRAY:
        v = pkd->particles.sph(p).c;
        break;
    case OUT_HSPH_ARRAY:
        v = pkd->particles.ball(p) * 0.5;
        break;
    case OUT_POS_VECTOR:
        v = pkd->particles.position(p)[iDim];
        break;
    case OUT_VEL_VECTOR:
        assert(pkd->particles.present(PKD_FIELD::oVelocity)); /* Validate memory model */
        v = pkd->particles.velocity(p)[iDim];
        break;
    case OUT_MEANVEL_VECTOR:
        assert(pkd->particles.present(PKD_FIELD::oVelSmooth)); /* Validate memory model */
        v = pkd->particles.get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth).vmean[iDim];
        break;
    case OUT_ACCEL_VECTOR:
        assert(pkd->particles.present(PKD_FIELD::oAcceleration)); /* Validate memory model */
        v = pkd->particles.acceleration(p)[iDim];
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
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset - ctx->inBuffer);
    if ( n < 24 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset - ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%" PRIu64 "\n",fetchInteger(pkd,p,iType,iDim));
    assert(strlen(ctx->inOffset) < 24 );
    while ( *ctx->inOffset ) ++ctx->inOffset;
}
static void storeFloat(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 20 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%.8g\n",fetchFloat(pkd,p,iType,iDim));
    assert(strlen(ctx->inOffset) < 20 );
    while ( *ctx->inOffset ) ++ctx->inOffset;
}

extern uint64_t hilbert2d(float x,float y);
extern uint64_t hilbert3d(float x,float y,float z);
static void storeRungDest(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int iRung;
    float x,y,z;
    int64_t lKey;
    const auto &RungDest = pkd->particles.RungDest(p);


    x = pkdPos(pkd,p,0) + 1.5;
    if (x < 1.0) x = 1.0;
    else if (x >= 2.0) x = 2.0;
    y = pkdPos(pkd,p,1) + 1.5;
    if (y < 1.0) y = 1.0;
    else if (y >= 2.0) y = 2.0;
    z = pkdPos(pkd,p,2) + 1.5;
    if (z < 1.0) z = 1.0;
    else if (z >= 2.0) z = 2.0;

#if PEANO_HILBERT_KEY_MAX > 0x3ffffffffffll
    lKey = hilbert3d(x,y,z);
#else
    lKey = hilbert2d(x,y);
#endif
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 100 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%016" PRIx64 " %d",lKey,static_cast<int>(p->uRung));
    ctx->inOffset += strlen(ctx->inOffset);
    n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    for (iRung=0; iRung<8; iRung++) {
        snprintf(ctx->inOffset,n," %d", RungDest[iRung]);
        ctx->inOffset += strlen(ctx->inOffset);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    *ctx->inOffset++ = '\n';
}
static void storeHdr(PKD pkd,PKDOUT ctx,uint64_t N) {
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 20 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%" PRIu64 "\n",N);
    while ( *ctx->inOffset ) ++ctx->inOffset;
}
static void finish(PKD pkd,PKDOUT ctx) {
    (*ctx->fnFlush)(pkd,ctx,0); /* Flush input buffer */
    (*ctx->fnFlush)(pkd,ctx,1); /* Finish output stream */
}

static void storePsGroup(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 40 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%" PRIu64 " %i\n",(uint64_t)p->iOrder, pkdGetGroup(pkd,p));
    assert(strlen(ctx->inOffset) < 40 );
    while ( *ctx->inOffset ) ++ctx->inOffset;
}

/******************************************************************************\
 * Generic binary output - flush needs to be customized
\******************************************************************************/
static void storeIntegerBinary(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int n = ctx->inOffset - ctx->inBuffer;
    if ( PKDOUT_BUFFER_SIZE - n < sizeof(int) ) (*ctx->fnFlush)(pkd,ctx,0);
    *(int *)ctx->inOffset = fetchInteger(pkd,p,iType,iDim);
    ctx->inOffset += sizeof(int);
}
static void storeFloatBinary(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    int n = ctx->inOffset - ctx->inBuffer;
    if ( PKDOUT_BUFFER_SIZE - n < sizeof(float) ) (*ctx->fnFlush)(pkd,ctx,0);
    *(float *)ctx->inOffset = fetchFloat(pkd,p,iType,iDim);
    ctx->inOffset += sizeof(float);
}

extern uint64_t hilbert2d(float x,float y);
extern uint64_t hilbert3d(float x,float y,float z);
static void storeRungDestBinary(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    assert(0);
    int iRung;
    float x,y,z;
    int64_t lKey;
    const auto &RungDest = pkd->particles.RungDest(p);


    x = pkdPos(pkd,p,0) + 1.5;
    if (x < 1.0) x = 1.0;
    else if (x >= 2.0) x = 2.0;
    y = pkdPos(pkd,p,1) + 1.5;
    if (y < 1.0) y = 1.0;
    else if (y >= 2.0) y = 2.0;
    z = pkdPos(pkd,p,2) + 1.5;
    if (z < 1.0) z = 1.0;
    else if (z >= 2.0) z = 2.0;

#if PEANO_HILBERT_KEY_MAX > 0x3ffffffffffll
    lKey = hilbert3d(x,y,z);
#else
    lKey = hilbert2d(x,y);
#endif
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 100 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%016" PRIx64 " %d",lKey,static_cast<int>(p->uRung));
    ctx->inOffset += strlen(ctx->inOffset);
    n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    for (iRung=0; iRung<8; iRung++) {
        snprintf(ctx->inOffset,n," %d", RungDest[iRung]);
        ctx->inOffset += strlen(ctx->inOffset);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    *ctx->inOffset++ = '\n';
}
static void storeHdrBinary(PKD pkd,PKDOUT ctx,uint64_t N) {
    int n = ctx->inOffset - ctx->inBuffer;
    if ( PKDOUT_BUFFER_SIZE - n < sizeof(int) ) (*ctx->fnFlush)(pkd,ctx,0);
    *(int *)ctx->inOffset = N;
    ctx->inOffset += sizeof(int);
}
static void storePsGroupBinary(PKD pkd,PKDOUT ctx,PARTICLE *p,int iType,int iDim) {
    assert(0);
    int n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    if ( n < 40 ) {
        (*ctx->fnFlush)(pkd,ctx,0);
        n = PKDOUT_BUFFER_SIZE - (ctx->inOffset-ctx->inBuffer);
    }
    snprintf(ctx->inOffset,n,"%" PRIu64 " %i\n",(uint64_t)p->iOrder, pkdGetGroup(pkd,p));
    assert(strlen(ctx->inOffset) < 40 );
    while ( *ctx->inOffset ) ++ctx->inOffset;
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
        ctx->outBuffer = ctx->outBuffer->next = new PKDOUTBUFFER;
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
    while (ctx->CTX.bzStream->avail_in) {
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
        } while (bzerror>0&&bzerror!=BZ_STREAM_END);
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

    ctx->CTX.bzStream = new bz_stream;
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
        ctx->outBuffer = ctx->outBuffer->next = new PKDOUTBUFFER;
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

    while (ctx->CTX.gzStream->avail_in) {
        if (ctx->CTX.gzStream->avail_out==0) {
            ctx->nBytes += (ctx->outBuffer->nBytes=PKDOUT_BUFFER_SIZE);
            (*ctx->fnWrite)(pkd,ctx);
        }
        gzerror = deflate(ctx->CTX.gzStream,0);
        assert(gzerror>=0);
    }
    if (final) {
        do {
            if (ctx->CTX.gzStream->avail_out==0) {
                ctx->nBytes += (ctx->outBuffer->nBytes=PKDOUT_BUFFER_SIZE);
                (*ctx->fnWrite)(pkd,ctx);
            }
            gzerror = deflate(ctx->CTX.gzStream,Z_FINISH);
            assert(gzerror>=0);
            (*ctx->fnWrite)(pkd,ctx);
        } while (gzerror!=Z_STREAM_END);
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

    ctx->CTX.gzStream = new z_stream;
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
    ctx = new struct pkdout; assert(ctx!=NULL);
    ctx->fp = NULL;
    ctx->outBuffer = ctx->headBuffer = new PKDOUTBUFFER;
    assert(ctx->outBuffer!=NULL);
    ctx->outBuffer->next = NULL;
    ctx->inBuffer = ctx->inOffset = new char[PKDOUT_BUFFER_SIZE];
    assert(ctx->inBuffer!=NULL);
    ctx->nBytes = 0;

    if (iFile == PKDOUT_TYPE_BINARY) {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            ctx->fnOut = storeIntegerBinary;
            break;
        case OUTTYPE_FLOAT:
            ctx->fnOut = storeFloatBinary;
            break;
        case OUTTYPE_RUNGDEST:
            ctx->fnOut = storeRungDestBinary;
            break;
        case OUTTYPE_PSGROUP:
            ctx->fnOut = storePsGroupBinary;
            break;

        default:
            assert(0);
        }
        ctx->fnHdr = storeHdrBinary;
    }
    else {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            ctx->fnOut = storeInteger;
            break;
        case OUTTYPE_FLOAT:
            ctx->fnOut = storeFloat;
            break;
        case OUTTYPE_RUNGDEST:
            ctx->fnOut = storeRungDest;
            break;
        case OUTTYPE_PSGROUP:
            ctx->fnOut = storePsGroup;
            break;

        default:
            assert(0);
        }
        ctx->fnHdr = storeHdr;
    }
    ctx->fnClose = finish;

    switch (iFile) {
    case PKDOUT_TYPE_BINARY:
        break;
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
    for (buf=ctx->headBuffer; buf!=NULL; buf=buf->next) {
        fwrite(buf->data,buf->nBytes,1,fp);
    }
}

/*
** Free up all buffers, and the PKDOUT context
*/
void pkdFreeOutASCII(PKD pkd,PKDOUT ctx) {
    PKDOUTBUFFER *buf, *nxt;

    delete [] ctx->inBuffer;
    for (buf=ctx->headBuffer; buf!=NULL; buf=nxt) {
        nxt = buf->next;
        delete buf;
    }
    delete ctx;
}

/******************************************************************************\
 * File I/O routines
\******************************************************************************/

PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iFile,int iType) {
    PKDOUT ctx;

    ctx = new struct pkdout; assert(ctx!=NULL);
    ctx->fp = NULL;
    ctx->outBuffer = ctx->headBuffer = new PKDOUTBUFFER;
    assert(ctx->outBuffer!=NULL);
    ctx->outBuffer->next = NULL;
    ctx->inBuffer = ctx->inOffset = new char [PKDOUT_BUFFER_SIZE];
    assert(ctx->inBuffer!=NULL);
    ctx->nBytes = 0;

    /* Determine how to handle the header and data */
    if (iFile == PKDOUT_TYPE_BINARY) {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            ctx->fnOut = storeIntegerBinary;
            break;
        case OUTTYPE_FLOAT:
            ctx->fnOut = storeFloatBinary;
            break;
        case OUTTYPE_RUNGDEST:
            ctx->fnOut = storeRungDestBinary;
            break;
        case OUTTYPE_PSGROUP:
            ctx->fnOut = storePsGroupBinary;
            break;

        default:
            assert(0);
        }
        ctx->fnHdr = storeHdrBinary;
    }
    else {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            ctx->fnOut = storeInteger;
            break;
        case OUTTYPE_FLOAT:
            ctx->fnOut = storeFloat;
            break;
        case OUTTYPE_RUNGDEST:
            ctx->fnOut = storeRungDest;
            break;
        case OUTTYPE_PSGROUP:
            ctx->fnOut = storePsGroup;
            break;
        default:
            assert(0);
        }
        ctx->fnHdr = storeHdr;
    }
    ctx->fp = fopen (pszFileName,mode);
    assert(ctx->fp != NULL);
    /*WTF: corrupts!!! setvbuf(ctx->fp,NULL,_IOFBF,PKDOUT_BUFFER_SIZE);*/

    switch (iFile) {

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
    delete [] ctx->inBuffer;
    delete ctx;
}

void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N) {
    (*ctx->fnHdr)(pkd,ctx,N);
}

void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim) {
    /*
    ** Write Elements!
    */
    for (auto &p : pkd->particles) {
        (*ctx->fnOut)(pkd,ctx,&p,iType,iDim);
    }
}

#ifdef USE_HDF5
void pkdOutHDF5(PKD pkd,char *pszFileName,int iType,int iDim) {
    assert(0);
}
#endif
