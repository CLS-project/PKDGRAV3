#ifdef HAVE_CONFIG_H
#include "config.h"
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
    prm->script_argc = 0;
    prm->script_argv = NULL;
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


void prmAddParam(PRM prm,char *pszName,int iType,void *pValue,
		 int iSize,char *pszArg,char *pszArgUsage) {
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


void prmArgUsage(PRM prm) {
    PRM_NODE *pn;

    if (prm->fcnLeader) (*prm->fcnLeader)();
    pn = prm->pnHead;
    while (pn) {
	if (pn->pszArg && pn->pszArgUsage) {
	    if (pn->iType == 0) {
		printf("[+%s][-%s] %s\n",pn->pszArg,pn->pszArg,
		       pn->pszArgUsage);
		}
	    else {
		printf("[-%s %s]\n",pn->pszArg,pn->pszArgUsage);
		}
	    }
	pn = pn->pnNext;
	}
    if (prm->fcnTrailer) (*prm->fcnTrailer)();
    }

/*
** This function saves all of the input parameters, as well as single-variable
** state information.
*/
void prmSave(PRM prm, FIO fio) {
    PRM_NODE *pn;

    /* We really shouldn't know about this structure, but what can you do? */
    for( pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	assert(pn->pValue);
	switch (pn->iType) {
	case 0:
	case 1:
	    assert(pn->iSize == sizeof(int));
	    fioSetAttr(fio,pn->pszName,FIO_TYPE_INT,pn->pValue);
	    break;
	case 2:
	    assert(pn->iSize == sizeof(double));
	    fioSetAttr(fio,pn->pszName,FIO_TYPE_DOUBLE,pn->pValue);
	    break;
	case 3:
	    fioSetAttr(fio,pn->pszName,FIO_TYPE_STRING,pn->pValue);
	    break;
	case 4:
	    assert(pn->iSize == sizeof(uint64_t));
	    fioSetAttr(fio,pn->pszName,FIO_TYPE_UINT64,pn->pValue);
	    break;
	    }
	}
    }

/*
** If there was a parameter file specified, then parse it the "old"
** way without using python.  We keep this around in case python is
** not available or uses too much memory.
*/
int prmParseParam(PRM prm) {
    FILE *fpParam;
    PRM_NODE *pn;
    char achBuf[PRM_LINE_SIZE];
    char *p,*q,*pszCmd,t;
    int iLine,ret;
    int bWarnQuote=0;

    if (prm->pszFilename==NULL) return(1);
    fpParam = fopen(prm->pszFilename,"r");
    if (!fpParam) {
	printf("Could not open file:%s\n",prm->pszFilename);
	return(0);
	}
    p = fgets(achBuf,PRM_LINE_SIZE,fpParam);
    iLine = 1;
    while (p) {
	if (*p == 0) goto new_line;
	if (*p == '#') goto new_line;
	while (isspace((int) *p)) {
	    ++p;
	    if (*p == 0) goto new_line;
	    }
	if (isalpha((int) *p)) {
	    pszCmd = p;
	    ++p;
	    if (*p == 0) goto lookup_cmd;
	    }
	else goto syntax_error;
	while (isalnum((int) *p)||strchr("_$",*p)) {
	    ++p;
	    if (*p == 0) goto lookup_cmd;
	    }
    lookup_cmd:
	t = *p;
	*p = 0;
	pn = prm->pnHead;
	while (pn) {
	    if (!strcmp(pszCmd,pn->pszName)) break;
	    pn = pn->pnNext;
	    }
	if (!pn) goto cmd_error;
	*p = t;
	if (*p == 0) goto syntax_error;
	while (isspace((int) *p)) {
	    ++p;
	    if (*p == 0) goto syntax_error;
	    }
	if (*p != '=') goto syntax_error;
	++p;
	if (*p == 0) goto syntax_error;
	while (isspace((int) *p)) {
	    ++p;
	    if (*p == 0) goto syntax_error;
	    }
	pn->bFile = 1;
	/* The command line is authoritative */
	if (pn->bArg) goto new_line;
	switch (pn->iType) {
	case 0:
	    assert(pn->iSize == sizeof(int));
	    ret = sscanf(p,"%d",(int *)pn->pValue);
	    if (ret != 1) goto syntax_error;
	    break;
	case 1:
	    assert(pn->iSize == sizeof(int));
	    ret = sscanf(p,"%d",(int *)pn->pValue);
	    if (ret != 1) goto syntax_error;
	    break;
	case 2:
	    assert(pn->iSize == sizeof(double));
	    ret = sscanf(p,"%lf",(double *)pn->pValue);
	    if (ret != 1) goto syntax_error;
	    break;
	case 3:
	    /*
	    ** Make sure there is enough space to handle the string.
	    ** This is a CONSERVATIVE test.
	    */
	    assert((size_t)pn->iSize > strlen(p));
	    ret = sscanf(p,"%[^\n#]",(char *)pn->pValue);
	    if (ret != 1) goto syntax_error;
	    /*
	     ** Strip trailing whitespace. OKAY!
	     */
	    p = pn->pValue;
	    q = &p[strlen(p)];
	    while (--q >= p) if (!isspace((int) *q)) break;
	    if ( *q=='"' && *p=='"' && p != q) {
		*q = 0;
		strcpy(p,p+1);
		}
	    else {
		*++q = 0;
		if ( !bWarnQuote ) {
		    fprintf(stderr,"WARNING: strings in parameter file should now be enclosed in quotes\n");
		    bWarnQuote=1;
		    }
		}
	    break;
	case 4:
	    assert(pn->iSize == sizeof(uint64_t));
	    ret = sscanf(p,"%"PRIu64,(uint64_t *)pn->pValue);
	    if (ret != 1) goto syntax_error;
	    break;
	default:
	    goto cmd_error;
	    }
    new_line:
	p = fgets(achBuf,PRM_LINE_SIZE,fpParam);
	++iLine;
	}
    fclose(fpParam);
    return(1);
syntax_error:
    q = achBuf;
    while (*q) {
	if (*q == '\n') *q = 0;
	else ++q;
	}
    printf("Syntax error in %s(%d):\n%s",prm->pszFilename,iLine,achBuf);
    fclose(fpParam);
    return(0);
cmd_error:
    q = achBuf;
    while (*q) {
	if (*q == '\n') *q = 0;
	else ++q;
	}
    printf("Unrecognized command in %s(%d):%s\n",prm->pszFilename,iLine,achBuf);
    fclose(fpParam);
    return(0);
    }

int prmArgProc(PRM prm,int argc,char **argv) {
    int i,ret;
    PRM_NODE *pn;

    if (argc < 2) return(1);

    /*
    ** Run through all command line arguments, setting the correct value
    ** for each.  If the very last item is not an option or option value,
    ** then it must be a parameter file, or python script if python is
    ** supported.
    */
    for (i=1;i<argc;++i) {
	if ((*argv[i] == '-' || *argv[i] == '+') && argv[i][1] ) {
	    pn = prm->pnHead;
	    while (pn) {
		if (pn->pszArg)
		    if (!strcmp(&argv[i][1],pn->pszArg)) break;
		pn = pn->pnNext;
		}
	    if (!pn) {
		printf("Unrecognized command line argument:%s\n",argv[i]);
		prmArgUsage(prm);
		return(0);
		}
	    }
#ifdef USE_PYTHON
	else {
	    int j;
	    if ( strcmp(argv[i],"-") != 0 ) {
		prm->pszFilename = argv[i];
		}
	    i++;
	    prm->script_argc = argc - i + 1;
	    prm->script_argv = malloc(prm->script_argc * sizeof(char *));
	    assert(prm->script_argv != NULL);
	    for( j=1; j<prm->script_argc; j++)
		prm->script_argv[j] = argv[i++];
	    prm->script_argv[0] = NULL;
	    return(1);
	    }
#else
	else if (i == argc-1) {
	    prm->pszFilename = argv[i];
	    return(1);
	    }
	else {
	    printf("Unrecognized command line argument:%s\n",argv[i]);
	    prmArgUsage(prm);
	    return(0);
	    }
#endif
	switch (pn->iType) {
	case 0:
	    /*
	     ** It's a boolean.
	     */
	    if (argv[i][0] == '-') *((int *)pn->pValue) = 0;
	    else *((int *)pn->pValue) = 1;
	    break;
	case 1:
	    /*
	     ** It's an int
	     */
	    ++i;
	    if (i == argc) {
		printf("Missing integer value after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    ret = sscanf(argv[i],"%d",(int *) pn->pValue);
	    if (ret != 1) {
		printf("Expected integer after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    break;
	case 2:
	    /*
	     ** It's a DOUBLE
	     */
	    ++i;
	    if (i == argc) {
		printf("Missing double value after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    ret = sscanf(argv[i],"%lf",(double *)pn->pValue);
	    if (ret != 1) {
		printf("Expected double after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    break;
	case 3:
	    /*
	     ** It's a string
	     */
	    ++i;
	    if (i == argc) {
		printf("Missing string after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    assert(pn->iSize > strlen(argv[i]));
	    strcpy((char *)pn->pValue,argv[i]);
	    break;
	case 4:
	    /*
	     ** It's an uint64_t
	     */
	    ++i;
	    if (i == argc) {
		printf("Missing integer value after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    ret = sscanf(argv[i],"%"PRIu64,(uint64_t *) pn->pValue);
	    if (ret != 1) {
		printf("Expected integer after command line ");
		printf("argument:%s\n",argv[i-1]);
		prmArgUsage(prm);
		return(0);
		}
	    break;
	default:
	    assert(0);
	    }
	pn->bArg = 1;
	}

#ifdef USE_PYTHON
    prm->script_argc = 1;
    prm->script_argv = malloc(prm->script_argc * sizeof(char *));
    assert(prm->script_argv != NULL);
    prm->script_argv[0] = NULL;
#endif

    return(1);
    }


int prmArgSpecified(PRM prm,char *pszName) {
    PRM_NODE *pn;

    pn = prm->pnHead;
    while (pn) {
	if (pn->pszArg)
	    if (!strcmp(pn->pszName,pszName)) break;
	pn = pn->pnNext;
	}
    if (!pn) return(0);
    return(pn->bArg);
    }


int prmFileSpecified(PRM prm,char *pszName) {
    PRM_NODE *pn;

    pn = prm->pnHead;
    while (pn) {
	if (!strcmp(pn->pszName,pszName)) break;
	pn = pn->pnNext;
	}
    if (!pn) return(0);
    return(pn->bFile);
    }


int prmSpecified(PRM prm,char *pszName) {
    return(prmArgSpecified(prm,pszName) || prmFileSpecified(prm,pszName));
    }
