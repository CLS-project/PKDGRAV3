#include <Python.h>
#include <structmember.h> // for PyMemberDef
#include "m_parse.h"

#include "master.h"
#include "csmpython.h"

#define MASTER_MODULE_NAME "MASTER"

/******************************************************************************\
*   MSR Module methods
\******************************************************************************/

typedef struct {
    PyObject_HEAD
    MSR msr;
    } MSRINSTANCE;

static PyObject *
ppy_msr_GenerateIC(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"csm","grid","seed",NULL};
    int nGrid, iSeed;
    CSMINSTANCE *cosmo = NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "O|ii:GenerateIC", const_cast<char **>(kwlist),
	     &cosmo,&self->msr->param.nGrid, &self->msr->param.iSeed ) )
	return NULL;
    if (!PyObject_TypeCheck(cosmo,&csmType)) return PyErr_Format(PyExc_TypeError,"Expected CSM object");
    self->msr->csm->val = cosmo->csm->val;
    if (self->msr->csm->val.classData.bClass)
        csmClassGslInitialize(self->msr->csm);
    msrGenerateIC(self->msr);
    Py_RETURN_NONE;
    }

static PyMethodDef msr_methods[] = {
    {"GenerateIC", (PyCFunction)ppy_msr_GenerateIC, METH_VARARGS|METH_KEYWORDS,
     "Generate Initial Condition"},
    {NULL, NULL, 0, NULL}
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
    PyObject * strMSR = Py_BuildValue("s",MASTER_MODULE_NAME);
    auto msr_module = PyImport_GetModule(strMSR);
    Py_DECREF(strMSR);
    self->msr = * reinterpret_cast<MSR*>(PyModule_GetState(msr_module));
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

static struct PyModuleDef msrModule = { 
    PyModuleDef_HEAD_INIT,
    MASTER_MODULE_NAME,
    "pkdgrav3 orchestration module",
    sizeof(MSR),
    NULL, // Methods
    NULL, // Slots
    NULL, NULL, NULL
};

// This routine is called if the user does an "import msr" from the script.
// We initialize the MSR module which will later cause "main" to be skipped.
static PyObject * initModuleMSR(void) {
    auto msr_module = PyModule_Create(&msrModule);
    auto pmsr = reinterpret_cast<MSR*>(PyModule_GetState(msr_module));
    msrInitialize(pmsr,mdlMDL(),mdlWORKER(),0,NULL);
    auto mainModule = PyImport_AddModule("__main__"); // Borrowed reference
    //PyObject *mainDict = PyModule_GetDict(mainModule);
    //setConstants(mainDict);
    if (PyType_Ready(&msrType) >= 0) { // Initialize "MSR" object as well.
	Py_INCREF(&msrType);
	PyModule_AddObject(msr_module, "MSR", (PyObject *)&msrType);
	}
    return msr_module;
    }

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

/* Copy parameters from python dictionary back into parameters. */
static int ppy2prm(PRM prm,PyObject *arguments, PyObject *specified) {
    PRM_NODE *pn;
    int bOK = 1;

    for( pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	//auto v = PyDict_GetItemString(arguments, pn->pszName);
	//auto f = PyDict_GetItemString(specified, pn->pszName);
	auto v = PyObject_GetAttrString(arguments, pn->pszName); // These are a Namespace
	auto f = PyObject_GetAttrString(specified, pn->pszName); // These are a Namespace
	if (v!=NULL) {
	    if (v == Py_None) continue;
	    pn->bArg = (f && PyObject_IsTrue(f)>0);
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

/******************************************************************************\
*   Setup MSR using Python to parse parameters / enter analysis mode
\******************************************************************************/

extern "C"
int msrPython(MSR *msr, int argc, char *argv[]) {
    PyImport_AppendInittab(MASTER_MODULE_NAME,initModuleMSR);
    { PyObject *PyInit_CSM(void); PyImport_AppendInittab("CSM",PyInit_CSM); }
    Py_InitializeEx(0);

    auto wargv = new wchar_t *[argc];
    for(int i=0; i<argc; ++i) wargv[i] = Py_DecodeLocale(argv[i],NULL);
    PySys_SetArgv(argc, wargv);

    auto globals = PyDict_New();
    auto locals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__",PyEval_GetBuiltins());
//    PyObject *PARSE = PyImport_ImportModule("parse");
//    if (PyErr_Occurred()) { /* Use the internal parser */
    auto PARSE = PyModule_New("parse");
    PyModule_AddStringConstant(PARSE, "__file__", "parse.py");
    PyObject *localDict = PyModule_GetDict(PARSE);
    PyDict_SetItemString(localDict, "__builtins__", PyEval_GetBuiltins());
    PyObject *pyValue = PyRun_String(parse_py, Py_file_input, localDict, localDict);

    PyObject *parse = PyObject_GetAttrString(PARSE, "parse");
    if (!PyCallable_Check(parse)) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST be callable\n"); abort(); }
    PyObject *update = PyObject_GetAttrString(PARSE, "update");
    if (!PyCallable_Check(update)){ fprintf(stderr,"INTERNAL ERROR: parse.update() MUST be callable\n"); abort(); }
//    PyObject *master = PyObject_GetAttrString(PARSE, "main");
    PyObject *result= PyObject_CallObject(parse,NULL);
    if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
    int n = PyTuple_Size(result);
    if (n!=2) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST return a tuple\n"); abort();	}
    PyObject *arguments = PyTuple_GetItem(result,0); /* Values of each parameter */
    PyObject *specified = PyTuple_GetItem(result,1); /* If it was explicitely specified */
    PyObject *script = PyObject_GetAttrString(arguments,"script");
    if (script != Py_None) {
    	char *filename;
	if (PyUnicode_Check(script)) {
	    PyObject *ascii = PyUnicode_AsASCIIString(script);
	    filename = PyBytes_AsString(ascii);
	    Py_DECREF(ascii);
	    }
	else { fprintf(stderr,"INTERNAL ERROR: script filename is invalid\n"); abort();	}
	FILE *fp = fopen(filename,"r");
	if (fp == NULL) { perror(filename); exit(errno); }
	PyRun_FileEx(fp,filename,Py_file_input,globals,locals,1);
	fclose(fp);
	if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	}

    PyObject * strMSR = Py_BuildValue("s",MASTER_MODULE_NAME);
    auto msr_module = PyImport_GetModule(strMSR);
    int bImported = msr_module != NULL;
    Py_DECREF(strMSR);

    if (msr_module == NULL) { // We must prepare for a normal legacy execution
    	msrInitialize(msr,mdlMDL(),mdlWORKER(),0,NULL);
	PyObject *args = PyTuple_New(3);
	PyTuple_SetItem(args,0,locals);
	PyTuple_SetItem(args,1,arguments);
	PyTuple_SetItem(args,2,specified);
	PyObject_CallObject(update,args);
	if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	ppy2prm((*msr)->prm,arguments,specified);
	}
    else *msr = * reinterpret_cast<MSR*>(PyModule_GetState(msr_module));

 //    if (!rc) {
 //    	printf("main() will be automatically run\n");
	// PyObject *args = PyTuple_New(3);
	// PyTuple_SetItem(args,0,locals);
	// PyTuple_SetItem(args,1,arguments);
	// PyTuple_SetItem(args,2,specified);
	// PyObject_CallObject(update,args);
	// if (PyErr_Occurred()) { PyErr_Print(); exit(1); }

	// args = PyTuple_New(2);
	// PyTuple_SetItem(args,0,arguments);
	// PyTuple_SetItem(args,1,specified);
	// PyObject_CallObject(master,args);
	// if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	// }
    return msr_module != NULL;
    }