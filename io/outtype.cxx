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
#include <cmath>
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
    case OUT_CACHECOLL_ARRAY:
    case OUT_AVOIDEDFLUXES_ARRAY:
    case OUT_COMPUTEDFLUXES_ARRAY:
        return OUTTYPE_INTEGER;

    case OUT_RUNGDEST_ARRAY:
        return OUTTYPE_RUNGDEST;

    default:
        return OUTTYPE_UNKNOWN;
    }
}

/* Write an integer */
auto PKDOUT::fetchInteger(particle P, int iType, int iDim) -> uint64_t {
    uint64_t v;

    switch (iType) {
    case OUT_IORDER_ARRAY:
        v = P.order();
        break;
    case OUT_GROUP_ARRAY:
        v = P.group();
        break;
    case OUT_GLOBALGID_ARRAY:
        v = P.global_gid();
        break;
    case OUT_MARKED_ARRAY:
        v = P.is_marked();
        break;
    case OUT_PSGROUP_ARRAY:
        assert(0);
        v = P.group();
        break;
#ifdef DEBUG_CACHED_FLUXES
    case OUT_CACHEFLUX_ARRAY:
        if (P.is_gas()) {
            v = P.sph().flux_cache;
        }
        else v = 0;
        break;
    case OUT_CACHECOLL_ARRAY:
        if (P.is_gas()) {
            v = P.sph().coll_cache;
        }
        else v = 0;
        break;
    case OUT_AVOIDEDFLUXES_ARRAY:
        if (P.is_gas()) {
            v = P.sph().avoided_fluxes;
        }
        else v = 0;
        break;
    case OUT_COMPUTEDFLUXES_ARRAY:
        if (P.is_gas()) {
            v = P.sph().computed_fluxes;
        }
        else v = 0;
        break;
#endif
    default:
        v = 0;
    }
    return v;
}
auto PKDOUT::fetchFloat(particle P, int iType, int iDim) -> double {
    double v;
    switch (iType) {
    case OUT_DENSITY_ARRAY:
        v = P.density();
        break;
    case OUT_BALL_ARRAY:
        v = P.ball();
        break;
    case OUT_POT_ARRAY:
        assert(P.have_potential());
        v = P.potential();
        break;
    case OUT_AMAG_ARRAY:
        assert(P.have_acceleration()); /* Validate memory model */
        v = sqrt(blitz::dot(P.acceleration(), P.acceleration()));
        break;
    case OUT_RUNG_ARRAY:
        v = P.rung();
        break;
    case OUT_SOFT_ARRAY:
        v = P.soft();
        break;
    case OUT_DIVV_ARRAY:
        assert(P.have_vel_smooth()); /* Validate memory model */
        v = P.VelSmooth().divv;
        break;
    case OUT_VELDISP2_ARRAY:
        assert(P.have_vel_smooth()); /* Validate memory model */
        v = P.VelSmooth().veldisp2;
        break;
    case OUT_VELDISP_ARRAY:
        assert(P.have_vel_smooth()); /* Validate memory model */
        v = sqrt(P.VelSmooth().veldisp2);
        break;
    case OUT_PHASEDENS_ARRAY:
        assert(P.have_vel_smooth()); /* Validate memory model */
        v = P.density() * pow(P.VelSmooth().veldisp2, -1.5);
        break;
    case OUT_C_ARRAY:
        v = P.sph().c;
        break;
    case OUT_HSPH_ARRAY:
        v = P.ball() * 0.5;
        break;
    case OUT_POS_VECTOR:
        v = P.position(iDim);
        break;
    case OUT_VEL_VECTOR:
        assert(P.have_velocity()); /* Validate memory model */
        v = P.velocity()[iDim];
        break;
    case OUT_MEANVEL_VECTOR:
        assert(P.have_vel_smooth()); /* Validate memory model */
        v = P.VelSmooth().vmean[iDim];
        break;
    case OUT_ACCEL_VECTOR:
        assert(P.have_acceleration()); /* Validate memory model */
        v = P.acceleration()[iDim];
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
void PKDOUT::storeInteger(particle P, int iType, int iDim) {
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 24) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%" PRIu64 "\n", fetchInteger(P, iType, iDim));
    assert(strlen(inOffset) < 24);
    while (*inOffset) ++inOffset;
}
void PKDOUT::storeFloat(particle P, int iType, int iDim) {
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 20) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%.8g\n", fetchFloat(P, iType, iDim));
    assert(strlen(inOffset) < 20);
    while (*inOffset) ++inOffset;
}

extern uint64_t hilbert2d(float x, float y);
extern uint64_t hilbert3d(float x, float y, float z);
void PKDOUT::storeRungDest(particle P, int iType, int iDim) {
    int iRung;
    int64_t lKey;
    const auto& RungDest = P.RungDest();
    auto r = blitz::min(blitz::max(P.position() + 1.5, 1.0), std::nextafter(2.0, 0.0));

#if PEANO_HILBERT_KEY_MAX > 0x3ffffffffffll
    lKey = hilbert3d(r[0], r[1], r[2]);
#else
    lKey = hilbert2d(r[0], r[1]);
#endif
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 100) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%016" PRIx64 " %d", lKey, static_cast<int>(P.rung()));
    inOffset += strlen(inOffset);
    n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    for (iRung = 0; iRung < 8; iRung++) {
        snprintf(inOffset, n, " %d", RungDest[iRung]);
        inOffset += strlen(inOffset);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    *inOffset++ = '\n';
}
void PKDOUT::storeHdr(uint64_t N) {
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 20) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%" PRIu64 "\n", N);
    while (*inOffset) ++inOffset;
}
void PKDOUT::finish() {
    (this->*fnFlush)(0); /* Flush input buffer */
    (this->*fnFlush)(1); /* Finish output stream */
}

void PKDOUT::storePsGroup(particle P, int iType, int iDim) {
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 40) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%" PRIu64 " %i\n", P.order(), P.group());
    assert(strlen(inOffset) < 40);
    while (*inOffset) ++inOffset;
}

/******************************************************************************\
 * Generic binary output - flush needs to be customized
\******************************************************************************/
void PKDOUT::storeIntegerBinary(particle P, int iType, int iDim) {
    int n = inOffset - inBuffer;
    if (PKDOUT_BUFFER_SIZE - n < sizeof(int)) (this->*fnFlush)(0);
    *(int*)inOffset = fetchInteger(P, iType, iDim);
    inOffset += sizeof(int);
}
void PKDOUT::storeFloatBinary(particle P, int iType, int iDim) {
    int n = inOffset - inBuffer;
    if (PKDOUT_BUFFER_SIZE - n < sizeof(float)) (this->*fnFlush)(0);
    *(float*)inOffset = fetchFloat(P, iType, iDim);
    inOffset += sizeof(float);
}

extern uint64_t hilbert2d(float x, float y);
extern uint64_t hilbert3d(float x, float y, float z);
void PKDOUT::storeRungDestBinary(particle P, int iType, int iDim) {
    assert(0);
    int iRung;
    int64_t lKey;
    const auto& RungDest = P.RungDest();
    auto r = blitz::min(blitz::max(P.position() + 1.5, 1.0), std::nextafter(2.0, 0.0));

#if PEANO_HILBERT_KEY_MAX > 0x3ffffffffffll
    lKey = hilbert3d(r[0], r[1], r[2]);
#else
    lKey = hilbert2d(r[0], r[1]);
#endif
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 100) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%016" PRIx64 " %d", lKey, static_cast<int>(P.rung()));
    inOffset += strlen(inOffset);
    n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    for (iRung = 0; iRung < 8; iRung++) {
        snprintf(inOffset, n, " %d", RungDest[iRung]);
        inOffset += strlen(inOffset);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    *inOffset++ = '\n';
}
void PKDOUT::storeHdrBinary(uint64_t N) {
    int n = inOffset - inBuffer;
    if (PKDOUT_BUFFER_SIZE - n < sizeof(int)) (this->*fnFlush)(0);
    *(int*)inOffset = N;
    inOffset += sizeof(int);
}
void PKDOUT::storePsGroupBinary(particle P, int iType, int iDim) {
    assert(0);
    int n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    if (n < 40) {
        (this->*fnFlush)(0);
        n = PKDOUT_BUFFER_SIZE - (inOffset - inBuffer);
    }
    snprintf(inOffset, n, "%" PRIu64 " %i\n", P.order(), P.group());
    assert(strlen(inOffset) < 40);
    while (*inOffset) ++inOffset;
}

/******************************************************************************\
 * Regular ASCII
\******************************************************************************/

/*
** Flush the input buffer by writing it directly to a file
*/
void PKDOUT::writeASCII(int final) {
    int n = inOffset - inBuffer;
    fwrite(inBuffer, n, 1, fp);
    inOffset = inBuffer;
}

void PKDOUT::closeASCII() {
    (this->*fnFlush)(0); /* Flush input buffer */
}

/******************************************************************************\
 * bzip2 compression
\******************************************************************************/

#ifdef HAVE_LIBBZ2
void PKDOUT::closeBZ2() {
    (this->*fnFlush)(0); /* Flush input buffer */
    (this->*fnFlush)(1); /* Flush output buffer */
}
/*
** Write routine to write to an output file
*/
void PKDOUT::outputBZ2() {
    int n = CTX.bzStream->next_out - outBuffer->data;
    fwrite(outBuffer->data, n, 1, fp);
    CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
    CTX.bzStream->next_out = outBuffer->data;
}
/*
** Write routine to allocate a new buffer.
*/
void PKDOUT::bufferBZ2() {
    outBuffer->nBytes = PKDOUT_BUFFER_SIZE - CTX.bzStream->avail_out;
    if (CTX.bzStream->avail_out == 0) {
        outBuffer = outBuffer->next = new PKDOUTBUFFER;
        outBuffer->next = NULL;
        outBuffer->nBytes = 0;
        CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
        CTX.bzStream->next_out = outBuffer->data;
    }
}
/*
** Called when the input buffer is (nearly) full.  Empties the input buffer by
** compressing the data.  New output buffers are allocated as required.
*/
void PKDOUT::flushBZ2(int final) {
    int bzerror;
    CTX.bzStream->next_in = inBuffer;
    CTX.bzStream->avail_in = inOffset - inBuffer;
    inOffset = inBuffer;
    while (CTX.bzStream->avail_in) {
        if (CTX.bzStream->avail_out == 0) {
            nBytes += (outBuffer->nBytes = PKDOUT_BUFFER_SIZE);
            (*fnWrite)(pkd, ctx);
        }
        bzerror = BZ2_bzCompress(CTX.bzStream, BZ_RUN);
        assert(bzerror >= 0);
    }
    if (final) {
        do {
            bzerror = BZ2_bzCompress(CTX.bzStream, BZ_FINISH);
            (*fnWrite)(pkd, ctx);
        } while (bzerror > 0 && bzerror != BZ_STREAM_END);
        (*fnWrite)(pkd, ctx);
        outBuffer->nBytes = PKDOUT_BUFFER_SIZE - CTX.bzStream->avail_out;
        nBytes += outBuffer->nBytes;
        BZ2_bzCompressEnd(CTX.bzStream);
        CTX.bzStream = NULL;
    }
}

void PKDOUT::setupBZ2() {
    fnClose = closeBZ2;
    fnFlush = flushBZ2;
    fnWrite = NULL; /* Still must be set */

    CTX.bzStream = new bz_stream;
    assert(CTX.bzStream != NULL);
    CTX.bzStream->bzalloc = NULL;
    CTX.bzStream->bzfree = NULL;
    CTX.bzStream->opaque = NULL;
    CTX.bzStream->avail_out = PKDOUT_BUFFER_SIZE;
    CTX.bzStream->next_out = outBuffer->data;
    BZ2_bzCompressInit(CTX.bzStream, 9, 0, 0);
}
#endif

/******************************************************************************\
 * zlib compression
\******************************************************************************/

#ifdef HAVE_LIBZ
/*
** Flush the input buffer by writing it directly to a file
*/
void PKDOUT::closeZ() {
    (this->*fnFlush)(0); /* Flush input buffer */
    (this->*fnFlush)(1); /* Flush output buffer */
}
/*
** Write routine to write to an output file
*/
void PKDOUT::outputZ() {
    int n = (char*)CTX.gzStream->next_out - outBuffer->data;
    fwrite(outBuffer->data, n, 1, fp);
    CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
    CTX.gzStream->next_out = outBuffer->data;
}
/*
** Write routine to allocate a new buffer.
*/
void PKDOUT::bufferZ() {
    outBuffer->nBytes = PKDOUT_BUFFER_SIZE - CTX.gzStream->avail_out;
    if (CTX.gzStream->avail_out == 0) {
        outBuffer = outBuffer->next = new PKDOUTBUFFER;
        outBuffer->next = NULL;
        outBuffer->nBytes = 0;
        CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
        CTX.gzStream->next_out = outBuffer->data;
    }
}
/*
** Called when the input buffer is (nearly) full.  Empties the input buffer by
** compressing the data.  New output buffers are allocated as required.
*/
void PKDOUT::flushZ(int final) {
    int gzerror;
    CTX.gzStream->next_in = inBuffer;
    CTX.gzStream->avail_in = inOffset - inBuffer;
    inOffset = inBuffer;
    while (CTX.gzStream->avail_in) {
        if (CTX.gzStream->avail_out == 0) {
            nBytes += (outBuffer->nBytes = PKDOUT_BUFFER_SIZE);
            (*fnWrite)(pkd, ctx);
        }
        gzerror = deflate(ctx->CTX.gzStream, 0);
        assert(gzerror >= 0);
    }
    if (final) {
        do {
            if (CTX.gzStream->avail_out == 0) {
                nBytes += (outBuffer->nBytes = PKDOUT_BUFFER_SIZE);
                (*fnWrite)(pkd, ctx);
            }
            gzerror = deflate(CTX.gzStream, Z_FINISH);
            assert(gzerror >= 0);
            (*fnWrite)(pkd, ctx);
        } while (gzerror != Z_STREAM_END);
        outBuffer->nBytes = PKDOUT_BUFFER_SIZE - CTX.gzStream->avail_out;
        nBytes += outBuffer->nBytes;
        deflateEnd(CTX.gzStream);
        CTX.gzStream = NULL;
    }
}

void PKDOUT::setupZ() {
    fnClose = &closeZ;
    fnFlush = flushZ;
    fnWrite = NULL;

    CTX.gzStream = new z_stream;
    assert(CTX.gzStream != NULL);
    CTX.gzStream->zalloc = Z_NULL;
    CTX.gzStream->zfree = Z_NULL;
    CTX.gzStream->opaque = Z_NULL;
    CTX.gzStream->avail_out = PKDOUT_BUFFER_SIZE;
    CTX.gzStream->next_out = outBuffer->data;
    deflateInit2(CTX.gzStream, Z_BEST_COMPRESSION, Z_DEFLATED, 31, 8, Z_DEFAULT_STRATEGY);
}
#endif

/******************************************************************************\
 * Storage routines
\******************************************************************************/

PKDOUT::PKDOUT(PKD pkd, std::string_view filename, int iFile, int iType) : PKDOUT(pkd, iFile, iType) {
    fp = fopen(filename.data(), "wb");
    assert(fp != NULL);
    /*WTF: corrupts!!! setvbuf(fp,NULL,_IOFBF,8192);*/
}

PKDOUT::PKDOUT(PKD pkd, int iFile, int iType) : pkd(pkd), fnWrite(nullptr) {

    /*
    ** Allocate the context, input buffer, and the first output buffer.
    */
    fp = nullptr;
    outBuffer = headBuffer = new PKDOUTBUFFER;
    assert(outBuffer != NULL);
    outBuffer->next = NULL;
    inBuffer = inOffset = new char[PKDOUT_BUFFER_SIZE];
    assert(inBuffer != NULL);
    nBytes = 0;

    if (iFile == PKDOUT_TYPE_BINARY) {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            fnOut = &PKDOUT::storeIntegerBinary;
            break;
        case OUTTYPE_FLOAT:
            fnOut = &PKDOUT::storeFloatBinary;
            break;
        case OUTTYPE_RUNGDEST:
            fnOut = &PKDOUT::storeRungDestBinary;
            break;
        case OUTTYPE_PSGROUP:
            fnOut = &PKDOUT::storePsGroupBinary;
            break;

        default:
            assert(0);
        }
        fnHdr = &PKDOUT::storeHdrBinary;
    }
    else {
        switch (getType(iType)) {
        case OUTTYPE_INTEGER:
            fnOut = &PKDOUT::storeInteger;
            break;
        case OUTTYPE_FLOAT:
            fnOut = &PKDOUT::storeFloat;
            break;
        case OUTTYPE_RUNGDEST:
            fnOut = &PKDOUT::storeRungDest;
            break;
        case OUTTYPE_PSGROUP:
            fnOut = &PKDOUT::storePsGroup;
            break;

        default:
            assert(0);
        }
        fnHdr = &PKDOUT::storeHdr;
    }
    fnClose = &PKDOUT::finish;

    switch (iFile) {
    case PKDOUT_TYPE_BINARY:
        break;
#ifdef HAVE_LIBBZ2
    case PKDOUT_TYPE_BZIP2:
        setupBZ2();
        fnWrite = bufferBZ2;
        break;
#endif
#ifdef HAVE_LIBZ
    case PKDOUT_TYPE_ZLIB:
        setupZ();
        fnWrite = bufferZ;
        break;
#endif
    default:
        fnFlush = &PKDOUT::writeASCII;
        fnClose = &PKDOUT::closeASCII;
    }
}

/*
** Finish the output stream by flushing any remaining characters in the output
** stream, and then finalizing it.  Compressed data is still available in the
** output buffers
*/
void PKDOUT::FinishOutASCII() {
    (this->*fnClose)();
}

/*
** Returns the total number of bytes available in the output buffers.  This will
** only return a valid number after Finish has been called.
*/
auto PKDOUT::CountOutASCII() -> uint64_t {
    return nBytes;
}

void PKDOUT::DumpOutASCII(FILE* fp) {
    PKDOUTBUFFER* buf;
    for (buf = headBuffer; buf != NULL; buf = buf->next) {
        fwrite(buf->data, buf->nBytes, 1, fp);
    }
}

/*
** Free up all buffers, and the PKDOUT context
*/
PKDOUT::~PKDOUT() {
    PKDOUTBUFFER *buf, *nxt;

    delete[] inBuffer;
    for (buf = headBuffer; buf != nullptr; buf = nxt) {
        nxt = buf->next;
        delete buf;
    }
}

void PKDOUT::CloseOutASCII() {
    (this->*fnClose)();
    fclose(fp);
}

void PKDOUT::OutHdr(uint64_t N) {
    (this->*fnHdr)(N);
}

void PKDOUT::OutASCII(int iType, int iDim) {
    /*
    ** Write Elements!
    */
    for (auto& P : pkd->particles) {
        (this->*fnOut)(P, iType, iDim);
    }
}
