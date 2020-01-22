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
#include "pkdtinypy.h"

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
#ifdef USE_PYTHON
    if (prm->script_argv) free(prm->script_argv);
#endif
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

#ifdef USE_PYTHON
static void setNode(PRM_NODE *pn,int i,PyObject *v) {
    const char *s;
    switch(pn->iType) {
    case 0:
    case 1:
	assert(pn->iSize == sizeof(int));
#if PY_MAJOR_VERSION >= 3
	if (PyLong_Check(v)) ((int *)pn->pValue)[i] = PyLong_AsLong(v);
#else
	if (PyInt_Check(v)) ((int *)pn->pValue)[i] = PyInt_AsLong(v);
#endif		
	else if (PyFloat_Check(v)) ((int *)pn->pValue)[i] = (int)PyFloat_AsDouble(v);
	else fprintf(stderr,"Invalid type for %s\n",pn->pszName);
	break;
    case 2:
	assert(pn->iSize == sizeof(double));
	if (PyFloat_Check(v)) ((double *)pn->pValue)[i] = PyFloat_AsDouble(v);
#if PY_MAJOR_VERSION >= 3
	else if (PyLong_Check(v)) ((double *)pn->pValue)[i] = PyLong_AsLong(v);
#else
	else if (PyInt_Check(v)) ((double *)pn->pValue)[i] = PyInt_AsLong(v);
#endif		
	else fprintf(stderr,"Invalid type for %s\n",pn->pszName);
	break;
    case 3:
#if PY_MAJOR_VERSION >= 3
	if (PyUnicode_Check(v)) {
	    PyObject *ascii = PyUnicode_AsASCIIString(v);
	    s = PyBytes_AsString(ascii);
	    Py_DECREF(ascii);
	    }
#else 
	if (PyString_Check(v)) s = PyString_AsString(v);
#endif
	else {
	    fprintf(stderr,"Invalid type for %s\n",pn->pszName);
	    s = NULL;
	    }
	if (s!=NULL) {
	    assert(pn->iSize > strlen(s));
	    strcpy((char *)pn->pValue,s);
	    }
	else *(char *)pn->pValue = 0;
	break;
    case 4:
	assert(pn->iSize == sizeof(uint64_t));
#if PY_MAJOR_VERSION >= 3
	((uint64_t *)pn->pValue)[i] = PyLong_AsLong(v);
#else
	((uint64_t *)pn->pValue)[i] = PyInt_AsLong(v);
#endif
	break;
	}
    }

/* Copy parameters from python dictionary back into parameters. */
static int ppy2prm(PRM prm,PyObject *global) {
    PyObject *v;
    PRM_NODE *pn;
    int bOK = 1;

    for( pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	v = PyDict_GetItemString(global, pn->pszName);
	if (v!=NULL) {
	    if (v == Py_None) continue;
	    pn->bArg = 1;
	    if (PyList_Check(v)) {
		if (pn->pCount==NULL) {
	            fprintf(stderr,"The parameter %s cannot be a list!\n",pn->pszName);
	            bOK = 0;
		    }
		else {
		    int i, n = PyList_Size(v);
		    for(i=0; i<n; ++i) setNode(pn,i,PyList_GetItem(v,i));
		    *pn->pCount = n;
		    }
		}
	    else setNode(pn,0,v);
	    }
	}
    return bOK;
    }
#else
/*
** If there was a parameter file specified, then parse it the "old"
** way without using python.  We keep this around in case python is
** not available or uses too much memory.
*/
static const char *type_names[] = {
    "Boolean", "Integer", "Real", "String", "Long Integer"
};

static int setNode(PRM_NODE *pn,int i,tp_obj o) {
    switch(pn->iType) {
    	case 0:
        case 1:
        case 2:
        case 4:
            if (o.type == TP_NUMBER) {
                if (pn->iType == 4) ((uint64_t *)pn->pValue)[i] = o.number.val;
                else if (pn->iType < 2)  ((int *)pn->pValue)[i] = o.number.val;
                else                  ((double *)pn->pValue)[i] = o.number.val;
		if (pn->pCount) *pn->pCount = 1;
        	}
            else {
        	fprintf(stderr, "ERROR: %s has invalid type, expected %s\n",
        		pn->pszName, type_names[pn->iType]);
        	return 0;
        	}
            break;
        case 3:
            if (o.type == TP_STRING) {
                if (pn->iSize <= o.string.len) {
                    fprintf(stderr, "Parameter %s too long\n",pn->pszName);
                    return 0;
                    }
		memcpy((char *)pn->pValue,o.string.val,o.string.len);
		((char *)pn->pValue)[o.string.len] = 0;
		}
	    break;
	}
    return 1;
    }
#endif

int prmParseParam(PRM prm,void *msr) {
    FILE *fpParam;
    PRM_NODE *pn;
    char achBuf[PRM_LINE_SIZE];
    char *p,*q,*pszCmd,t;
    int iLine,ret;
    int bWarnQuote=0;
    extern void prm2ppy();

    if (prm->pszFilename==NULL) return(1);
    prm->script_argv[0] = prm->pszFilename;

#ifdef USE_PYTHON
    Py_Initialize();
    PyObject *mainModule = PyImport_AddModule("__main__"); 

    FILE *fp;
    PyObject *globals;
#if PY_MAJOR_VERSION > 2
    wchar_t **wargv;
    int i;
    wargv = malloc(sizeof(*wargv)*prm->script_argc);
    for(i=0; i<prm->script_argc; ++i) wargv[i] = Py_DecodeLocale(prm->script_argv[i],NULL);
#endif
    assert(Py_IsInitialized());
    assert(prm->script_argc>0);

    globals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__", PyEval_GetBuiltins());

#if PY_MAJOR_VERSION > 2
    PySys_SetArgv(prm->script_argc, wargv);
#else
    PySys_SetArgv(prm->script_argc, prm->script_argv);
#endif
    fp = fopen(prm->script_argv[0],"r");
    if (fp==NULL) {perror(prm->script_argv[0]); abort();}

    PyRun_SimpleFile(fp,prm->script_argv[0]);
    fclose(fp);

#if PY_MAJOR_VERSION > 2
    for(i=0; i<prm->script_argc; ++i) PyMem_RawFree(wargv[i]);
    free(wargv);
#endif
    return ppy2prm(prm,PyModule_GetDict(mainModule));
#else
    tp_vm *tp = tp_init(0, NULL);
    tpyInitialize(tp,msr);

    /* Look for, and retrieve all valid parameters */
    tp_obj dict = tp_main(tp,(char *)prm->pszFilename,0,0);
    for( pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
    	tp_obj o;
	if (tp_iget(tp,&o,dict,tp_string(pn->pszName))) {
	    pn->bFile = 1; // This came from the file
	    if (!pn->bArg) { // but command line takes precedence
		if (o.type == TP_LIST) {
		    int i;
		    for (i=0; i<o.list.val->len; ++i) {
			if (!setNode(pn,i,o.list.val->items[i])) return 0;
			}
		    *pn->pCount = o.list.val->len;
		    }
		else if (setNode(pn,0,o)) {
		    if (pn->pCount) *pn->pCount = 1;
		    }
		else return 0;
	    }
	    tp_del(tp,dict,tp_string(pn->pszName)); /* Delete as we processed it */
	    }
	}
    /* Now look for unprocessed parameters: this would be an error */
    int i,n = tp_len(tp,dict).number.val;
    int bOK = 1;
    for(i=0; i<n; ++i) {
	tp_obj o, k = tp_iter(tp,dict,tp_number(i));
	tp_iget(tp,&o,dict,k);
	const char *key = k.string.val;
	if (islower(*key) && o.type!=TP_NONE) {
	    int j;
	    for(j=0; j<k.string.len; ++j){
		if (isupper(key[j])) {
		    printf("Unrecognized parameter: %s\n", key);
		    bOK = 0;
		    }
		}
            }
	}

    tp_deinit(tp);
    return bOK;
#endif
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


int prmArgSpecified(PRM prm,const char *pszName) {
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


int prmFileSpecified(PRM prm,const char *pszName) {
    PRM_NODE *pn;

    pn = prm->pnHead;
    while (pn) {
	if (!strcmp(pn->pszName,pszName)) break;
	pn = pn->pnNext;
	}
    if (!pn) return(0);
    return(pn->bFile);
    }


int prmSpecified(PRM prm,const char *pszName) {
    return(prmArgSpecified(prm,pszName) || prmFileSpecified(prm,pszName));
    }
