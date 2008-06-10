#ifndef OUTTYPE_HINCLUDED
#define OUTTYPE_HINCLUDED

#include "pkd.h"

#define OUT_COLOR_ARRAY		1
#define OUT_DENSITY_ARRAY	2
#define OUT_POT_ARRAY		3
#define OUT_AMAG_ARRAY		4
#define OUT_IMASS_ARRAY		5
#define OUT_RUNG_ARRAY		6

#define OUT_DIVV_ARRAY          7
#define OUT_VELDISP2_ARRAY      8
#define OUT_PHASEDENS_ARRAY     9

#define OUT_SOFT_ARRAY          32

#define OUT_GROUP_ARRAY	112
#define OUT_GROUP_TIPSY_NAT	113
#define OUT_GROUP_TIPSY_STD	114
#define OUT_GROUP_STATS 115
#define OUT_GROUP_PROFILES 116
void pkdOutGroup(PKD,char *,int,int,double);

#define OUT_RELAX_ARRAY 120
#define OUT_POS_VECTOR		1
#define OUT_VEL_VECTOR		2
#define OUT_ACCEL_VECTOR	3
#define OUT_MEANVEL_VECTOR      4

void pkdOutArray(PKD,char *,int);
void pkdOutVector(PKD,char *,int,int);

#endif
