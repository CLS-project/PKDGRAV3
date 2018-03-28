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
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#else
#define PRIu64 "llu"
#endif
#include "intype.h"

static void *openMark(PKD pkd,FILE *fp) {
    uint64_t N, nDark, nGas, nStar;
    fscanf(fp,"%"PRIu64" %"PRIu64" %"PRIu64, &nDark, &nGas, &nStar);
    N = nDark + nGas + nStar;
    return NULL;
    }

static int readMark(void *ctx,FILE *fp,uint64_t i,int iType,int iDim) {
    uint64_t id;
    return fscanf(fp,"%"PRIu64, &id );
    }

static void doneMark(void *ctx) {
    }

void pkdInASCII(PKD pkd,char *pszFileName,int iType,int iDim) {
    FILE *fp;
    uint64_t i;
    void *ctx;
    void * (*fnOpen)(PKD pkd,FILE *fp);
    int (*fnRead)(void *ctx,FILE *fp,uint64_t i,int iType,int iDim);
    void (*fnDone)(void *ctx);

    switch(iType) {
    case IN_SRC_MARK:
    case IN_DST_MARK:
	fnOpen = openMark;
	fnRead = readMark;
	fnDone = doneMark;
	break;

    default:
	assert(0);
	}

    fp = fopen (pszFileName,"r");
    assert(fp != NULL);

    ctx = (*fnOpen)(pkd,fp);
    i = 0;
    while ((*fnRead)(ctx,fp,i++,iType,iDim)>0) {}
    (*fnDone)(ctx);


    if (fclose(fp) != 0) {
	perror("pkdInASCII: could not close file");
	exit(1);
	}
    }
