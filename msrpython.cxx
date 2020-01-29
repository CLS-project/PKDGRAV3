/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2020 Douglas Potter & Joachim Stadel
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
#include <Python.h>
#include <structmember.h> // for PyMemberDef
#include "m_parse.h"

#include "master.h"
#include "csmpython.h"

#define MASTER_MODULE_NAME "MASTER"

/******************************************************************************\
*   Copy parameters from Python to pkdgrav3 (may go away eventually)
\******************************************************************************/

static void setNode(PRM_NODE *pn,int i,PyObject *v) {
    const char *s;
    switch(pn->iType) {
    case 0:
    case 1:
	assert(pn->iSize == sizeof(int));
	if (PyLong_Check(v)) ((int *)pn->pValue)[i] = PyLong_AsLong(v);
	else if (PyFloat_Check(v)) ((int *)pn->pValue)[i] = (int)PyFloat_AsDouble(v);
	else fprintf(stderr,"Invalid type for %s\n",pn->pszName);
	break;
    case 2:
	assert(pn->iSize == sizeof(double));
	if (PyFloat_Check(v)) ((double *)pn->pValue)[i] = PyFloat_AsDouble(v);
	else if (PyLong_Check(v)) ((double *)pn->pValue)[i] = PyLong_AsLong(v);
	else fprintf(stderr,"Invalid type for %s\n",pn->pszName);
	break;
    case 3:
	if (PyUnicode_Check(v)) {
	    PyObject *ascii = PyUnicode_AsASCIIString(v);
	    s = PyBytes_AsString(ascii);
	    Py_DECREF(ascii);
	    }
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
	((uint64_t *)pn->pValue)[i] = PyLong_AsLong(v);
	break;
	}
    }

static int ppy2prm(PRM prm,PyObject *arguments, PyObject *specified) {
    int bOK = 1;

    for( auto pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	//auto v = PyDict_GetItemString(arguments, pn->pszName);
	//auto f = PyDict_GetItemString(specified, pn->pszName);
	auto v = PyObject_GetAttrString(arguments, pn->pszName); // A Namespace
	if (v!=NULL) {
	    if (v != Py_None) {
		auto f = PyObject_GetAttrString(specified, pn->pszName); // A Namespace
		if (f) { pn->bArg = PyObject_IsTrue(f)>0; Py_DECREF(v); }
		else pn->bArg = 0;
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
	    Py_DECREF(v);
	    }
	}
    return bOK;
    }

/******************************************************************************\
*   Copy parameters from pkdgrav3 to Python (may go away eventually)
\******************************************************************************/

static void prm2ppy(PRM prm,PyObject *arguments, PyObject *specified) {
    for( auto pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	if (pn->pCount!=NULL) continue; // Lists are read-only for now
	PyObject *v;

	switch(pn->iType) {
	case 0:
	case 1:
	    v = PyLong_FromLong(*(int *)pn->pValue);
	    break;
	case 2:
	    v = PyFloat_FromDouble(*(double *)pn->pValue);
	    break;
	default:
	    v = NULL;
	    }
	if (v) PyObject_SetAttrString(arguments,pn->pszName,v);
	}
    }

void msrSaveParameters(MSR msr) {
    prm2ppy(msr->prm,msr->arguments,msr->specified);
    }

/******************************************************************************\
*   MSR Module methods
\******************************************************************************/

typedef struct {
    PyObject_HEAD
    MSR msr;
    } MSRINSTANCE;

static PyObject *
ppy_msr_GenerateIC(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"csm","z","grid","seed","L",NULL};
    int nGrid, iSeed;
    double dBoxSize = 1.0;
    CSMINSTANCE *cosmo = NULL;
    self->msr->param.bPeriodic = 1;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "Odi|id:GenerateIC", const_cast<char **>(kwlist),
	     &cosmo,&self->msr->param.dRedFrom,&self->msr->param.nGrid,
	     &self->msr->param.iSeed, &self->msr->param.dBoxSize) )
	return NULL;
    if (!PyObject_TypeCheck(cosmo,&csmType)) return PyErr_Format(PyExc_TypeError,"Expected CSM object");
    self->msr->csm->val = cosmo->csm->val;
    if (self->msr->csm->val.classData.bClass)
        csmClassGslInitialize(self->msr->csm);
    double dTime = msrGenerateIC(self->msr);
    return Py_BuildValue("d", dTime );
    }

static PyObject *
ppy_msr_DomainDecomp(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"rung",NULL};
    int iRung    = 0;
    int bOthers  = 0;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|i:DomainDecomp", const_cast<char **>(kwlist),
	     &iRung ) )
	return NULL;
    msrDomainDecomp(self->msr,iRung,bOthers);
    Py_RETURN_NONE;
    }

static PyObject *
ppy_msr_Checkpoint(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={/*"name",*/"step","time",NULL};
//    const char *fname;
    int iStep = 0;
    double dTime = 1.0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|id:Checkpoint", const_cast<char **>(kwlist),
	     /*&fname,*/&iStep,&dTime ) )
	return NULL;
    msrCheckpoint(self->msr,iStep,dTime);
    Py_RETURN_NONE;
    }

static PyObject *
ppy_msr_Restart(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    MSR msr = self->msr;
    static char const *kwlist[]={"arguments","specified","species","classes","n","name","step","time","E","U", "Utime", NULL};
    PyObject *species, *classes;
    int n, iStep;
    const char *name;
    double dTime;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "OOOOisidddd:Restart", const_cast<char **>(kwlist),
	     &msr->arguments,&msr->specified,&species,&classes,&n,&name,
	     &iStep,&dTime,&msr->dEcosmo,&msr->dUOld, &msr->dTimeOld ) )
	return NULL;

    // Create a vector of number of species
    species = PySequence_Fast(species,"species must be a list");
    int nSpecies = PySequence_Fast_GET_SIZE(species);
    std::vector<uint64_t> vecSpecies;
    vecSpecies.reserve(nSpecies);
    for(auto i=0; i < nSpecies; ++i) {
        PyObject *item = PySequence_Fast_GET_ITEM(species, i);
        vecSpecies.push_back(PyNumber_AsSsize_t(item,NULL));
	}
    Py_DECREF(species); // PySequence_Fast creates a new reference
    assert(vecSpecies.size()==4);
    msr->N     = vecSpecies[0];
    msr->nDark = vecSpecies[1];
    msr->nGas  = vecSpecies[2];
    msr->nStar = vecSpecies[3];

    // Process the array of class information
    classes = PySequence_Fast(classes,"species must be a list");
    int nClasses = PySequence_Fast_GET_SIZE(classes);
    msr->nCheckpointClasses = nClasses;
    for(auto i=0; i < nClasses; ++i) {
        PyObject *item = PySequence_Fast_GET_ITEM(classes, i);
	auto cls = PySequence_Fast(item,"class entry must be a list");
	assert(PySequence_Fast_GET_SIZE(cls)==3);
	PyObject *itemSpecies = PySequence_Fast_GET_ITEM(cls, 0);
	PyObject *itemMass    = PySequence_Fast_GET_ITEM(cls, 1);
	PyObject *itemSoft    = PySequence_Fast_GET_ITEM(cls, 2);
	msr->aCheckpointClasses[i].eSpecies = (FIO_SPECIES)PyNumber_AsSsize_t(itemSpecies,NULL);
	msr->aCheckpointClasses[i].fMass = PyFloat_AsDouble(itemMass);
	msr->aCheckpointClasses[i].fSoft = PyFloat_AsDouble(itemSoft);
	Py_DECREF(cls); // PySequence_Fast creates a new reference
	}
    Py_DECREF(classes); // PySequence_Fast creates a new reference

    ppy2prm(msr->prm,msr->arguments,msr->specified);

    msrRestart(msr, n, name, iStep, dTime);

    Py_RETURN_NONE;
    }

static PyObject *
ppy_msr_Write(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"name","time",NULL};
    const char *fname;
    double dTime = 1.0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|d:Write", const_cast<char **>(kwlist),
	     &fname,&dTime ) )
	return NULL;
    msrWrite(self->msr,fname,dTime,0);
    Py_RETURN_NONE;
    }

static PyMethodDef msr_methods[] = {
    {"GenerateIC", (PyCFunction)ppy_msr_GenerateIC, METH_VARARGS|METH_KEYWORDS,
     "Generate Initial Condition"},
    {"DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
    {"Checkpoint", (PyCFunction)ppy_msr_Checkpoint, METH_VARARGS|METH_KEYWORDS,
     "Write a checkpoint"},
    {"Restart", (PyCFunction)ppy_msr_Restart, METH_VARARGS|METH_KEYWORDS,
     "Restart from a checkpoint"},
    {"Write", (PyCFunction)ppy_msr_Write, METH_VARARGS|METH_KEYWORDS,
     "Write a particle output"},
    {NULL, NULL, 0, NULL}
};

/******************************************************************************\
*   MSR Module Definition
\******************************************************************************/

struct msrModuleState {
    MSR msr;
    bool bImported;
    };

static struct PyModuleDef msrModule = { 
    PyModuleDef_HEAD_INIT,
    MASTER_MODULE_NAME,
    "pkdgrav3 orchestration module",
    sizeof(struct msrModuleState),
    NULL, // Methods
    NULL, // Slots
    NULL, NULL, NULL
};

/******************************************************************************\
*   MSR OBJECT Boilerplate code
\******************************************************************************/

static void msr_dealloc(MSRINSTANCE *self) {
    Py_TYPE(self)->tp_free((PyObject*)self);
    }

static PyObject *msr_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    MSRINSTANCE *self = reinterpret_cast<MSRINSTANCE *>(type->tp_alloc(type, 0));
    if (self == NULL) { return NULL; }
    // Make a copy of the MSR object
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState*>(PyModule_GetState(msr_module));
    self->msr = moduleState->msr;
    return reinterpret_cast<PyObject *>(self);
    }

static int msr_init(MSRINSTANCE *self, PyObject *args, PyObject *kwds) {
    return 0;
    }

static PyMemberDef msr_members[] = {
	{NULL} /* Sentinel */
    };

static PyGetSetDef msr_getseters[] = {
    /*{"son", (getter)Father_getson, (setter)Father_setson, "son", NULL},*/
	{NULL} /* Sentinel */
    };

static PyTypeObject msrType = {
    PyVarObject_HEAD_INIT(NULL,0)
    "MSR", /*tp_name*/
    sizeof(MSRINSTANCE), /*tp_basicsize */
    0, /*tp_itemsize*/
    (destructor)msr_dealloc, /*tp_dealloc*/
    0, /*tp_print*/
    0, /*tp_getattr*/
    0, /*tp_setattr*/
    0, /*tp_compare*/
    0, /*tp_repr*/
    0, /*tp_as_number*/
    0, /*tp_as_sequence*/
    0, /*tp_as_mapping*/
    0, /*tp_hash */
    0, /*tp_call*/
    0, /*tp_str*/
    0, /*tp_getattro*/
    0, /*tp_setattro*/
    0, /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "MSR objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    msr_methods, /* tp_methods */
    msr_members, /* tp_members */
    msr_getseters, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc)msr_init, /* tp_init */
    0, /* tp_alloc */
    msr_new, /* tp_new */
    };

/******************************************************************************\
*   MSR Module Boilerplate code
\******************************************************************************/

// This routine is called if the user does an "import msr" from the script.
// We initialize the MSR module which will later cause "main" to be skipped.
static PyObject * initModuleMSR(void) {
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState*>(PyModule_GetState(msr_module));
    moduleState->bImported = true;
    //auto mainModule = PyImport_AddModule("__main__"); // Borrowed reference
    //PyObject *mainDict = PyModule_GetDict(mainModule);
    //setConstants(mainDict);
    if (PyType_Ready(&msrType) >= 0) { // Initialize "MSR" object as well.
	Py_INCREF(&msrType);
	PyModule_AddObject(msr_module, "MSR", (PyObject *)&msrType);
	}
    return msr_module;
    }

/******************************************************************************\
*   Setup MSR using Python to parse parameters / enter analysis mode
\******************************************************************************/

extern "C" PyObject *PyInit_CSM(void);

int msrPython(MSR *msr, int argc, char *argv[]) {
    PyImport_AppendInittab(MASTER_MODULE_NAME,initModuleMSR);
    PyImport_AppendInittab("CSM",PyInit_CSM);
    Py_InitializeEx(0);

    // Contruct the "MSR" context and module
    auto msr_module = PyModule_Create(&msrModule);
    PyState_AddModule(msr_module,&msrModule);
    msrInitialize(msr,mdlMDL(),mdlWORKER(),0,NULL);
    auto moduleState = reinterpret_cast<struct msrModuleState*>(PyModule_GetState(msr_module));
    moduleState->msr = *msr;
    moduleState->bImported = false; // If imported then we enter script mode

    // Convert program arguments to unicode
    auto wargv = new wchar_t *[argc];
    for(int i=0; i<argc; ++i) wargv[i] = Py_DecodeLocale(argv[i],NULL);
    PySys_SetArgv(argc, wargv);

    auto globals = PyDict_New();
    auto locals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__",PyEval_GetBuiltins());

    // Parse the command line
    auto PARSE = PyModule_New("parse");
    PyModule_AddStringConstant(PARSE, "__file__", "parse.py");
    PyObject *localDict = PyModule_GetDict(PARSE);
    PyDict_SetItemString(localDict, "__builtins__", PyEval_GetBuiltins());
    PyObject *pyValue = PyRun_String(parse_py, Py_file_input, localDict, localDict);
    PyObject *parse = PyObject_GetAttrString(PARSE, "parse");
    if (!PyCallable_Check(parse)) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST be callable\n"); abort(); }
    PyObject *update = PyObject_GetAttrString(PARSE, "update");
    if (!PyCallable_Check(update)){ fprintf(stderr,"INTERNAL ERROR: parse.update() MUST be callable\n"); abort(); }
    PyObject *result= PyObject_CallObject(parse,NULL);
    if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
    // Retrieve the results
    int n = PyTuple_Size(result);
    if (n!=2) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST return a tuple\n"); abort();	}
    PyObject *arguments = PyTuple_GetItem(result,0); /* Values of each parameter */
    PyObject *specified = PyTuple_GetItem(result,1); /* If it was explicitely specified */
    PyObject *script = PyObject_GetAttrString(arguments,"script");
    (*msr)->arguments = arguments;
    (*msr)->specified = specified;

    // If a script was specified then we run it.
    if (script != Py_None) {
    	char *filename;
	PyObject *ascii;
	if (PyUnicode_Check(script)) {
	    ascii = PyUnicode_AsASCIIString(script);
	    filename = PyBytes_AsString(ascii);
	    }
	else { fprintf(stderr,"INTERNAL ERROR: script filename is invalid\n"); abort();	}
	FILE *fp = fopen(filename,"r");
	if (fp == NULL) { perror(filename); exit(errno); }
	PyRun_FileEx(fp,filename,Py_file_input,globals,locals,1); // fp is closed on return
	if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	Py_DECREF(ascii);
	}

    // If "MASTER" was imported then we are done -- the script should have done its job
    if (!moduleState->bImported) { // We must prepare for a normal legacy execution
	PyObject *args = PyTuple_New(3);
	PyTuple_SetItem(args,0,locals);
	PyTuple_SetItem(args,1,arguments);
	PyTuple_SetItem(args,2,specified);
	PyObject_CallObject(update,args); // Copy and variables into the arguments Namespace
	if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	ppy2prm((*msr)->prm,arguments,specified); // Update the pkdgrav parameter state
	}

    return moduleState->bImported;
    }