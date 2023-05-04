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

#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

enum OUT_TYPE {
    OUT_DENSITY_ARRAY        = 11,
    OUT_POT_ARRAY            = 12,
    OUT_AMAG_ARRAY           = 13,
    OUT_IMASS_ARRAY          = 14,
    OUT_RUNG_ARRAY           = 15,
    OUT_DIVV_ARRAY           = 16,
    OUT_VELDISP2_ARRAY       = 17,
    OUT_VELDISP_ARRAY        = 18,
    OUT_PHASEDENS_ARRAY      = 19,
    OUT_SOFT_ARRAY           = 20,
    OUT_POS_VECTOR           = 21,
    OUT_VEL_VECTOR           = 22,
    OUT_ACCEL_VECTOR         = 23,
    OUT_MEANVEL_VECTOR       = 24,

    OUT_IORDER_ARRAY         = 25,

    OUT_C_ARRAY              = 28,
    OUT_HSPH_ARRAY           = 29,

    OUT_RUNGDEST_ARRAY       = 30,
    OUT_MARKED_ARRAY         = 31,

    OUT_CACHEFLUX_ARRAY      = 32,
    OUT_CACHECOLL_ARRAY      = 33,
    OUT_AVOIDEDFLUXES_ARRAY  = 34,
    OUT_COMPUTEDFLUXES_ARRAY = 35,

    OUT_HOP_STATS            = 100,

    OUT_GROUP_ARRAY          = 112,
    OUT_GLOBALGID_ARRAY      = 113,
    OUT_BALL_ARRAY           = 121,
    OUT_PSGROUP_ARRAY        = 122,
    OUT_PSGROUP_STATS        = 123
};

#define PKDOUT_TYPE_ASCII 0
#define PKDOUT_TYPE_BINARY 1
#define PKDOUT_TYPE_ZLIB   2
#define PKDOUT_TYPE_BZIP2  3

#ifdef HAVE_LIBBZ2
    #include <bzlib.h>
#endif
#ifdef HAVE_LIBZ
    #include <zlib.h>
#endif

#define PKDOUT_BUFFER_SIZE (1024*1024)
typedef struct pkdOutBuffer {
    struct pkdOutBuffer *next;
    uint32_t nBytes;
    char data[PKDOUT_BUFFER_SIZE];
} PKDOUTBUFFER;

typedef struct pkdout {
    FILE *fp;
#if defined(HAVE_LIBBZ2) || defined(HAVE_LIBZ)
    union {
#ifdef HAVE_LIBBZ2
        bz_stream *bzStream;
#endif
#ifdef HAVE_LIBZ
        z_stream *gzStream;
#endif
    } CTX;
#endif
    uint64_t nBytes;
    char *inBuffer;
    char *inOffset;
    PKDOUTBUFFER *outBuffer;
    PKDOUTBUFFER *headBuffer;
    void (*fnOut)(PKD pkd,struct pkdout *ctx,PARTICLE *p,int iType,int iDim);
    void (*fnHdr)(PKD pkd,struct pkdout *ctx,uint64_t N);
    void (*fnWrite)(PKD pkd,struct pkdout *ctx);
    void (*fnFlush)(PKD pkd,struct pkdout *ctx,int final);
    void (*fnClose)(PKD pkd,struct pkdout *ctx);
} *PKDOUT;

#ifdef __cplusplus
extern "C" {
#endif
PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iFile,int iType);
void pkdCloseOutASCII(PKD pkd,PKDOUT ctx);
void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N);
void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim);

PKDOUT pkdStartOutASCII(PKD pkd,int iFile, int iType);
void pkdFinishOutASCII(PKD pkd,PKDOUT ctx);
uint64_t pkdCountOutASCII(PKD pkd,PKDOUT ctx);
void pkdDumpOutASCII(PKD pkd,PKDOUT ctx,FILE *fp);
void pkdFreeOutASCII(PKD pkd,PKDOUT ctx);
#ifdef __cplusplus
}
#endif
#endif
