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







