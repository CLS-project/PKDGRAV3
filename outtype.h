#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define OUT_TIPSY_STD            0
#define OUT_TIPSY_DBL            1

#define OUT_COLOR_ARRAY		10
#define OUT_DENSITY_ARRAY	11
#define OUT_POT_ARRAY		12
#define OUT_AMAG_ARRAY		13
#define OUT_IMASS_ARRAY		14
#define OUT_RUNG_ARRAY		15
#define OUT_DIVV_ARRAY          16
#define OUT_VELDISP2_ARRAY      17
#define OUT_VELDISP_ARRAY       18
#define OUT_PHASEDENS_ARRAY     19
#define OUT_SOFT_ARRAY          20
#define OUT_POS_VECTOR		21
#define OUT_VEL_VECTOR		22
#define OUT_ACCEL_VECTOR	23
#define OUT_MEANVEL_VECTOR      24

#define OUT_IORDER_ARRAY        25

#define OUT_UDOT_ARRAY          26
#define OUT_U_ARRAY             27
#define OUT_C_ARRAY             28
#define OUT_HSPH_ARRAY          29

#define OUT_RUNGDEST_ARRAY      30

#define OUT_GROUP_ARRAY	112
#define OUT_GROUP_TIPSY_NAT	113
#define OUT_GROUP_TIPSY_STD	114
#define OUT_GROUP_STATS 115
#define OUT_GROUP_PROFILES 116
#define OUT_RELAX_ARRAY 120
#define OUT_BALL_ARRAY 121
#define OUT_PSGROUP_ARRAY	122
#define OUT_PSGROUP_STATS	123

#define PKDOUT_TYPE_ASCII 0
#define PKDOUT_TYPE_ZLIB  1
#define PKDOUT_TYPE_BZIP2 2

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

PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iFile,int iType);
void pkdCloseOutASCII(PKD pkd,PKDOUT ctx);
void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N);
void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim);

PKDOUT pkdStartOutASCII(PKD pkd,int iFile, int iType);
void pkdFinishOutASCII(PKD pkd,PKDOUT ctx);
uint64_t pkdCountOutASCII(PKD pkd,PKDOUT ctx);
void pkdDumpOutASCII(PKD pkd,PKDOUT ctx,FILE *fp);
void pkdFreeOutASCII(PKD pkd,PKDOUT ctx);


void pkdOutGroup(PKD,char *,int,int,double);

#endif
