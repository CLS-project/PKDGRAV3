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

#define OUT_TIPSY_STD            0
#define OUT_TIPSY_DBL            1

#define OUT_DENSITY_ARRAY   11
#define OUT_POT_ARRAY       12
#define OUT_AMAG_ARRAY      13
#define OUT_IMASS_ARRAY     14
#define OUT_RUNG_ARRAY      15
#define OUT_DIVV_ARRAY          16
#define OUT_VELDISP2_ARRAY      17
#define OUT_VELDISP_ARRAY       18
#define OUT_PHASEDENS_ARRAY     19
#define OUT_SOFT_ARRAY          20
#define OUT_POS_VECTOR      21
#define OUT_VEL_VECTOR      22
#define OUT_ACCEL_VECTOR    23
#define OUT_MEANVEL_VECTOR      24

#define OUT_IORDER_ARRAY        25

#define OUT_UDOT_ARRAY          26
#define OUT_U_ARRAY             27
#define OUT_C_ARRAY             28
#define OUT_HSPH_ARRAY          29

#define OUT_RUNGDEST_ARRAY      30
#define OUT_MARKED_ARRAY        31

#define OUT_CACHEFLUX_ARRAY   32
#define OUT_CACHECOLL_ARRAY   33
#define OUT_AVOIDEDFLUXES_ARRAY   34
#define OUT_COMPUTEDFLUXES_ARRAY   35

#define OUT_HOP_STATS          100

#define OUT_GROUP_ARRAY        112
#define OUT_RELAX_ARRAY        120
#define OUT_BALL_ARRAY         121
#define OUT_PSGROUP_ARRAY      122
#define OUT_PSGROUP_STATS      123

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
