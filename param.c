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
#ifdef USE_PYTHON
    #include <Python.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#ifdef HAVE_INTTYPES_H
    #include <inttypes.h>
#else
    #define PRIu64 "llu"
#endif
#ifdef HAVE_MALLOC_H
    #include <malloc.h>
#endif
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include "param.h"

void prmInitialize(PRM *pprm,void (*fcnLeader)(void),void (*fcnTrailer)(void)) {
    PRM prm;

    prm = (PRM)malloc(sizeof(struct prmContext));
    assert(prm != NULL);
    *pprm = prm;
    prm->pnHead = NULL;
    prm->pszFilename = NULL;
    prm->fcnLeader = fcnLeader;
    prm->fcnTrailer = fcnTrailer;
}


void prmFinish(PRM prm) {
    PRM_NODE *pn,*pnKill;

    pn = prm->pnHead;
    while (pn) {
        pnKill = pn;
        pn = pn->pnNext;
        free(pnKill->pszName);
        if (pnKill->pszArg) free(pnKill->pszArg);
        if (pnKill->pszArgUsage) free(pnKill->pszArgUsage);
        free(pnKill);
    }
    free(prm);
}

void prmAddArray(PRM prm,const char *pszName,int iType,void *pValue,int iSize,int *pCount) {
    PRM_NODE *pn,*pnTail;
    pn = (PRM_NODE *)malloc(sizeof(PRM_NODE));
    assert(pn != NULL);
    pn->pszName = (char *)malloc(strlen(pszName)+1);
    assert(pn->pszName != NULL);
    strcpy(pn->pszName,pszName);
    pn->iType = iType;
    pn->iSize = iSize;
    pn->bArg = 0;
    pn->bFile = 0;
    pn->pValue = pValue;
    pn->pCount = pCount;
    pn->pszArg = NULL;
    pn->pszArgUsage = NULL;
    pn->pnNext = NULL;
    if (!prm->pnHead) prm->pnHead = pn;
    else {
        pnTail = prm->pnHead;
        while (pnTail->pnNext) pnTail = pnTail->pnNext;
        pnTail->pnNext = pn;
    }
}

void prmAddParam(PRM prm,const char *pszName,int iType,void *pValue,
                 int iSize,const char *pszArg,const char *pszArgUsage) {
    PRM_NODE *pn,*pnTail;

    pn = (PRM_NODE *)malloc(sizeof(PRM_NODE));
    assert(pn != NULL);
    pn->pszName = (char *)malloc(strlen(pszName)+1);
    assert(pn->pszName != NULL);
    strcpy(pn->pszName,pszName);
    pn->iType = iType;
    pn->iSize = iSize;
    pn->bArg = 0;
    pn->bFile = 0;
    pn->pValue = pValue;
    pn->pCount = NULL;
    if (pszArg) {
        pn->pszArg = (char *)malloc(strlen(pszArg)+1);
        assert(pn->pszArg != NULL);
        strcpy(pn->pszArg,pszArg);
    }
    else pn->pszArg = NULL;
    if (pszArgUsage) {
        pn->pszArgUsage = (char *)malloc(strlen(pszArgUsage)+1);
        assert(pn->pszArgUsage != NULL);
        strcpy(pn->pszArgUsage,pszArgUsage);
    }
    else pn->pszArgUsage = NULL;
    pn->pnNext = NULL;
    if (!prm->pnHead) prm->pnHead = pn;
    else {
        pnTail = prm->pnHead;
        while (pnTail->pnNext) pnTail = pnTail->pnNext;
        pnTail->pnNext = pn;
    }
}


int prmArgSpecified(PRM prm,const char *pszName) {
    PRM_NODE *pn;

    pn = prm->pnHead;
    while (pn) {
        if (pn->pszArg)
            if (!strcmp(pn->pszName,pszName)) break;
        pn = pn->pnNext;
    }
    if (!pn) return (0);
    return (pn->bArg);
}


int prmFileSpecified(PRM prm,const char *pszName) {
    PRM_NODE *pn;

    pn = prm->pnHead;
    while (pn) {
        if (!strcmp(pn->pszName,pszName)) break;
        pn = pn->pnNext;
    }
    if (!pn) return (0);
    return (pn->bFile);
}


int prmSpecified(PRM prm,const char *pszName) {
    return (prmArgSpecified(prm,pszName) || prmFileSpecified(prm,pszName));
}
