#ifndef PARAM_HINCLUDED
#define PARAM_HINCLUDED

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
void prmAddParam(PRM,char *,int,void *,int,char *,char *);
void prmArgUsage(PRM prm);
int prmParseParam(PRM);
int prmArgProc(PRM,int,char **);
int prmSpecified(PRM,char *);
int prmArgSpecified(PRM,char *);
int prmFileSpecified(PRM,char *);

#endif







