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

#ifndef PARAM_HINCLUDED
#define PARAM_HINCLUDED

#include "fio.h"

typedef struct prmNode {
    struct prmNode *pnNext;
    char *pszName;
    int iType;
    int bArg;
    int bFile;
    int iSize;
    void *pValue;
    char *pszArg;
    char *pszArgUsage;
    } PRM_NODE;

typedef struct prmContext {
    PRM_NODE *pnHead;
    const char *pszFilename;
    int script_argc;
    char **script_argv;
    void (*fcnLeader)(void);
    void (*fcnTrailer)(void);
    } * PRM;

#define PRM_LINE_SIZE	128

void prmInitialize(PRM *,void (*)(void),void (*)(void));
void prmFinish(PRM);
void prmAddParam(PRM,const char *,int,void *,int,const char *,const char *);
void prmArgUsage(PRM prm);
void prmSave(PRM prm, FIO fio);
int prmParseParam(PRM);
int prmArgProc(PRM,int,char **);
int prmSpecified(PRM,const char *);
int prmArgSpecified(PRM,const char *);
int prmFileSpecified(PRM,const char *);

#endif







