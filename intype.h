#ifndef INTYPE_HINCLUDED
#define INTYPE_HINCLUDED
#define INTYPE_H_MODULE_ID "$Id$"

#include "pkd.h"

#define IN_TIPSY_STD            0
#define IN_TIPSY_DBL            1
#define IN_TIPSY_NAT            2

#define IN_SRC_MARK             10
#define IN_DST_MARK             11

void pkdInASCII(PKD pkd,char *pszFileName,int iType,int iDim);

#endif
