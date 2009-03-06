#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED
#define OUTTYPE_H_MODULE_ID "$Id$"

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

#define OUT_GROUP_ARRAY	112
#define OUT_GROUP_TIPSY_NAT	113
#define OUT_GROUP_TIPSY_STD	114
#define OUT_GROUP_STATS 115
#define OUT_GROUP_PROFILES 116
#define OUT_RELAX_ARRAY 120

#ifdef HAVE_LIBBZ2
#include <bzlib.h>
#endif
#ifdef HAVE_LIBZ
#include <zlib.h>
#endif

typedef struct pkdout {
    FILE *fp;
    union {
#ifdef HAVE_LIBBZ2
	BZFILE *bz;
#endif
#ifdef HAVE_LIBZ
	gzFile gz;
#endif
	} CTX;
    void (*fnOut)(PKD pkd,void *fp,PARTICLE *p,int iType,int iDim);
    void (*fnHdr)(PKD pkd,void *ctx,uint64_t N);
    void (*fnClose)(PKD pkd,void *ctx);
    } *PKDOUT;

PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iType);
void pkdCloseOutASCII(PKD pkd,PKDOUT ctx);
void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N);
void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim);

void pkdOutGroup(PKD,char *,int,int,double);

#endif
