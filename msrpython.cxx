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
#include <string>
#include <Python.h>
#define USE_NUMPY
#ifdef USE_NUMPY
    #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
    #include <numpy/arrayobject.h>
#endif
#include <structmember.h> // for PyMemberDef
#include "m_parse.h"

#include "master.h"
#include "csmpython.h"

#define MASTER_MODULE_NAME "MASTER"
#define MASTER_TYPE_NAME "MSR"


/******************************************************************************\
*   Flush stdout and stderr
\******************************************************************************/
static void flush_std_files(void) {
    auto out = PySys_GetObject("stdout");
    auto err = PySys_GetObject("stderr");
    if (out != NULL && out != Py_None) {
        auto tmp = PyObject_CallMethod(out, "flush", "");
        if (tmp == NULL) PyErr_WriteUnraisable(out);
        else Py_DECREF(tmp);
    }

    if (err != NULL && err != Py_None) {
        auto tmp = PyObject_CallMethod(err, "flush", "");
        if (tmp == NULL) PyErr_Clear();
        else Py_DECREF(tmp);
    }
}

/******************************************************************************\
*   Modern, safe, direct from Python parameter fetching
\******************************************************************************/

int64_t MSR::getScalarInteger(const char *name, PyObject *v) {
    if (PyLong_Check(v)) return PyLong_AsLong(v);
    else if (PyFloat_Check(v)) return static_cast<int64_t>(PyFloat_AsDouble(v));
    else throw std::domain_error(name);
}
int64_t MSR::getScalarInteger(const char *name) {
    int64_t result;
    auto v = PyObject_GetAttrString(arguments, name);
    result = getScalarInteger(name,v);
    Py_DECREF(v);
    return result;
}
std::vector<int64_t> MSR::getVectorInteger(const char *name) {
    std::vector<int64_t> result;
    auto v = PyObject_GetAttrString(arguments, name);
    if (PyList_Check(v)) {
        auto n = PyList_Size(v);
        result.reserve(n);
        for (auto i=0; i<n; ++i) result.push_back(getScalarInteger(name,PyList_GetItem(v,i)));
    }
    else result.push_back(getScalarInteger(name));
    Py_DECREF(v);
    return result;
}

double MSR::getScalarNumber(const char *name, PyObject *v) {
    if (PyFloat_Check(v)) return PyFloat_AsDouble(v);
    else if (PyLong_Check(v)) return PyLong_AsLong(v);
    else throw std::domain_error(name);
}
double MSR::getScalarNumber(const char *name) {
    double result;
    auto v = PyObject_GetAttrString(arguments, name);
    result = getScalarNumber(name,v);
    Py_DECREF(v);
    return result;
}
std::vector<double> MSR::getVectorNumber(const char *name) {
    std::vector<double> result;
    auto v = PyObject_GetAttrString(arguments, name);
    if (PyList_Check(v)) {
        auto n = PyList_Size(v);
        result.reserve(n);
        for (auto i=0; i<n; ++i) result.push_back(getScalarNumber(name,PyList_GetItem(v,i)));
    }
    else result.push_back(getScalarNumber(name));
    Py_DECREF(v);
    return result;
}

std::string MSR::getScalarString(const char *name, PyObject *v) {
    std::string result;
    if (PyUnicode_Check(v)) {
        auto ascii = PyUnicode_AsASCIIString(v);
        result = PyBytes_AsString(ascii);
        Py_DECREF(ascii);
    }
    else throw std::domain_error(name);
    return result;
}
std::string MSR::getScalarString(const char *name) {
    std::string result;
    auto v = PyObject_GetAttrString(arguments, name);
    result = getScalarString(name,v);
    Py_DECREF(v);
    return result;
}
std::vector<std::string> MSR::getVectorString(const char *name) {
    std::vector<std::string> result;
    auto v = PyObject_GetAttrString(arguments, name);
    if (PyList_Check(v)) {
        auto n = PyList_Size(v);
        result.reserve(n);
        for (auto i=0; i<n; ++i) result.push_back(getScalarString(name,PyList_GetItem(v,i)));
    }
    else result.push_back(getScalarString(name));
    Py_DECREF(v);
    return result;
}


/******************************************************************************\
*   Copy parameters from Python to pkdgrav3 (may go away eventually)
\******************************************************************************/

static void setNode(PRM_NODE *pn,int i,PyObject *v) {
    const char *s;
    switch (pn->iType) {
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

    for ( auto pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
        //auto v = PyDict_GetItemString(arguments, pn->pszName);
        //auto f = PyDict_GetItemString(specified, pn->pszName);
        auto v = PyObject_GetAttrString(arguments, pn->pszName); // A Namespace
        if (v!=NULL) {
            if (v != Py_None) {
                auto f = PyObject_GetAttrString(specified, pn->pszName); // A Namespace
                if (f) { pn->bArg = PyObject_IsTrue(f)>0; Py_DECREF(f); }
                else pn->bArg = 0;
                if (PyList_Check(v)) {
                    if (pn->pCount==NULL) {
                        fprintf(stderr,"The parameter %s cannot be a list!\n",pn->pszName);
                        bOK = 0;
                    }
                    else {
                        int i, n = PyList_Size(v);
                        for (i=0; i<n; ++i) setNode(pn,i,PyList_GetItem(v,i));
                        *pn->pCount = n;
                    }
                }
                else setNode(pn,0,v);
            }
            Py_DECREF(v);
        }
        else {PyErr_Clear();}
    }
    return bOK;
}

/******************************************************************************\
*   Copy parameters from pkdgrav3 to Python (may go away eventually)
\******************************************************************************/

static void prm2ppy(PRM prm,PyObject *arguments, PyObject *specified) {
    for ( auto pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
        if (pn->pCount!=NULL) continue; // Lists are read-only for now
        PyObject *v;

        switch (pn->iType) {
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
        if (v) {
            PyObject_SetAttrString(arguments,pn->pszName,v);
            Py_DECREF(v);
        }
    }
}

void MSR::SaveParameters() {
    prm2ppy(prm,arguments,specified);
}

/******************************************************************************\
*   MSR Module Definition
\******************************************************************************/

struct msrModuleState {
    MSR *msr;
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
*   Ephemeral OBJECT Instance
\******************************************************************************/

struct EPHEMERALINSTANCE {
    PyObject_HEAD
    class MSR *msr;
    size_t nBytesPerNode;
    size_t nBytesPerParticle;
};

/******************************************************************************\
*   Ephemeral Object methods
\******************************************************************************/
static PyObject *
ephemeral_enter(EPHEMERALINSTANCE *self, PyObject *args, PyObject *kwobj) {
    if (!PyArg_ParseTuple(args, "")) return NULL;
    Py_INCREF(self);
    return reinterpret_cast<PyObject *>(self);
}

static PyObject *
ephemeral_exit(EPHEMERALINSTANCE *self, PyObject *args, PyObject *kwobj) {
    Py_RETURN_FALSE;
}

static PyObject *
ppy_ephemeral_add_grid(EPHEMERALINSTANCE *self, PyObject *args, PyObject *kwobj) {
    MSR *msr = self->msr;
    static char const *kwlist[]= {"grid","count",NULL};
    int grid, count=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "i|i:add_grid", const_cast<char **>(kwlist),
                &grid,&count) )
        return NULL;
    auto nBytes = msr->getLocalGridMemory(grid) * count;
    self->nBytesPerNode += nBytes;
    return Py_BuildValue("L",nBytes);
}

/******************************************************************************\
*   Ephemeral object methods list
\******************************************************************************/

static PyMethodDef ephemeral_methods[] = {
    {"__enter__",(PyCFunction)ephemeral_enter,METH_VARARGS,"__enter__"},
    {"__exit__", (PyCFunction)ephemeral_exit, METH_VARARGS,"__exit__"},

    {
        "add_grid", (PyCFunction)ppy_ephemeral_add_grid, METH_VARARGS|METH_KEYWORDS,
        "Add Ephemeral memory for one or more grids"
    },
    {NULL, NULL, 0, NULL}
};
/******************************************************************************\
*   Ephemeral OBJECT Boilerplate code - Keeps track of home much Ephemeral
\******************************************************************************/

static void ephemeral_dealloc(EPHEMERALINSTANCE *self) {
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *ephemeral_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    auto self = reinterpret_cast<EPHEMERALINSTANCE *>(type->tp_alloc(type, 0));
    if (self == NULL) { return NULL; }
    // Make a copy of the MSR object
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState *>(PyModule_GetState(msr_module));
    self->msr = moduleState->msr;
    return reinterpret_cast<PyObject *>(self);
}

static int ephemeral_init(EPHEMERALINSTANCE *self, PyObject *args, PyObject *kwds) {
    MSR *msr = self->msr;
    self->nBytesPerNode = 0;
    self->nBytesPerParticle = 0;
    static char const *kwlist[]= {"grid","count","per_particle",NULL};
    int grid=0, count=1, per_particle=0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwds, "|iii:ephemeral", const_cast<char **>(kwlist),
                &grid, &count, &per_particle) )
        return -1;
    if (grid) self->nBytesPerNode += msr->getLocalGridMemory(grid) * count;
    self->nBytesPerParticle = per_particle;
    return 0;
}

static PyMemberDef ephemeral_members[] = {
    {NULL} /* Sentinel */
};

static PyObject *ephemeral_get_bytes_per_node(EPHEMERALINSTANCE *self, void *) {
    return PyLong_FromSize_t(self->nBytesPerNode);
}

static PyObject *ephemeral_get_bytes_per_particle(EPHEMERALINSTANCE *self, void *) {
    return PyLong_FromSize_t(self->nBytesPerParticle);
}

// This warning should be fixed in newer Python versions
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwrite-strings"
static PyGetSetDef ephemeral_getseters[] = {
    {
        "bytes_per_node", (getter)ephemeral_get_bytes_per_node, NULL,
        "Number of bytes of ephemeral needed for each node", NULL
    },
    {
        "bytes_per_particle", (getter)ephemeral_get_bytes_per_particle, NULL,
        "Number of bytes of ephemeral needed for each particle", NULL
    },
    {NULL} /* Sentinel */
};
#pragma GCC diagnostic pop

static PyTypeObject ephemeralType = {
    PyVarObject_HEAD_INIT(NULL,0)
    MASTER_MODULE_NAME ".ephemeral", /*tp_name*/
    sizeof(EPHEMERALINSTANCE), /*tp_basicsize */
    0, /*tp_itemsize*/
    (destructor)ephemeral_dealloc, /*tp_dealloc*/
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
    "Ephemeral objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    ephemeral_methods, /* tp_methods */
    ephemeral_members, /* tp_members */
    ephemeral_getseters, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc)ephemeral_init, /* tp_init */
    0, /* tp_alloc */
    ephemeral_new, /* tp_new */
};

/******************************************************************************\
*   MSR Object method
\******************************************************************************/

/********** Set Internal parameters **********/

bool MSR::setParameters(PyObject *kwobj,bool bIgnoreUnknown) {
    bool bSuccess = true;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    auto allow = PyDict_GetItemString(kwobj,"bIgnoreUnknown");
    if (allow) bIgnoreUnknown = PyObject_IsTrue(allow)>0;

    while (PyDict_Next(kwobj, &pos, &key, &value)) {
        const char *keyString;
        if (PyUnicode_Check(key)) {
            PyObject *ascii = PyUnicode_AsASCIIString(key);
            keyString = PyBytes_AsString(ascii);
            Py_DECREF(ascii);
            if (keyString[0]=='_') continue;
        }
        if (PyObject_HasAttr(arguments,key)) {
            PyObject_SetAttr(arguments,key,value);
            PyObject_SetAttr(specified,key,Py_True);
        }
        else if (!bIgnoreUnknown) {
            PyErr_Format(PyExc_AttributeError,"invalid parameter %A",key);
            PyErr_Print();
            bSuccess=false;
        }
    }
    ppy2prm(prm,arguments,specified);
    ValidateParameters();
    return bSuccess;
}

static PyObject *
ppy_msr_setParameters(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    MSR *msr = self->msr;
    int n = PyTuple_Size(args);
    if (n) return PyErr_Format(PyExc_ValueError,"setParameters() takes 0 positional arguments but %d %s given",n,n==1?"was":"were");
    if (!msr->setParameters(kwobj)) return PyErr_Format(PyExc_AttributeError,"invalid parameters specified");
    Py_RETURN_NONE;
}

/********** Analysis Hooks **********/

static PyObject *
ppy_msr_ephemeral(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    MSR *msr = self->msr;
    auto ephemeral = reinterpret_cast<EPHEMERALINSTANCE *>(PyObject_Call(reinterpret_cast<PyObject *>(&ephemeralType),args,kwobj));
    if (ephemeral) ephemeral->msr = msr;
    return reinterpret_cast<PyObject *>(ephemeral);
}

static PyObject *
ppy_msr_add_analysis(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    MSR *msr = self->msr;
    static char const *kwlist[]= {"callback","memory",NULL};
    PyObject *callback, *memory;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "OO:add_analysis", const_cast<char **>(kwlist),
                &callback, &memory) )
        return NULL;
    if (!PyCallable_Check(callback)) return PyErr_Format(PyExc_AttributeError,"callback must be callable");
    msr->addAnalysis(callback,self,memory);
    Py_RETURN_NONE;
}

/********** Start normal simulation **********/

static PyObject *
ppy_msr_simulate(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    MSR *msr = self->msr;
    static char const *kwlist[]= {NULL};
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, ":simulate", const_cast<char **>(kwlist)) )
        return NULL;
    msr->Hostname(); // List all host names
    auto dTime = msr->LoadOrGenerateIC();
    if (dTime != -HUGE_VAL) msr->Simulate(dTime);
    Py_RETURN_NONE;
}

/********** Initial Condition Generation **********/

static PyObject *
ppy_msr_GenerateIC(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"csm","z","grid","seed","L",NULL};
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
    double dTime = self->msr->GenerateIC();
    return Py_BuildValue("d", dTime );
}

/********** File I/O **********/

static PyObject *
ppy_msr_Load(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"filename",NULL};
    const char *fname;
    double dTime;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "s:Load", const_cast<char **>(kwlist),
                &fname ) )
        return NULL;
    dTime = self->msr->Read(fname);
    return Py_BuildValue("d", dTime );
}

static PyObject *
ppy_msr_Save(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"name","time",NULL};
    const char *fname;
    double dTime = 1.0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "s|d:Save", const_cast<char **>(kwlist),
                &fname,&dTime ) )
        return NULL;
    self->msr->Write(fname,dTime,0);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_WriteArray(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"name","field","binary",NULL};
    const char *fname;
    int iField = 0;
    int bBinary = 0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "si|p:WriteArray", const_cast<char **>(kwlist),
                &fname,&iField,&bBinary ) )
        return NULL;
    self->msr->OutArray(fname,iField,bBinary);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_WriteVector(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"name","field","binary",NULL};
    const char *fname;
    int iField = 0;
    int bBinary = 0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "si|p:WriteVector", const_cast<char **>(kwlist),
                &fname,&iField,&bBinary ) )
        return NULL;
    self->msr->OutVector(fname,iField,bBinary);
    Py_RETURN_NONE;
}

/********** Checkpoint I/O **********/

static PyObject *
ppy_msr_Checkpoint(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {/*"name",*/"step","time",NULL};
//    const char *fname;
    int iStep = 0;
    int nSteps = 0;
    double dTime = 1.0;
    double dDelta = 0.0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|iidd:Checkpoint", const_cast<char **>(kwlist),
                /*&fname,*/&iStep,&nSteps,&dTime,&dDelta ) )
        return NULL;
    self->msr->Checkpoint(iStep,nSteps,dTime,dDelta);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_Restart(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    MSR *msr = self->msr;
    static char const *kwlist[]= {"arguments","specified","species","classes","n","name","step","steps","time","delta","E","U", "Utime", NULL};
    PyObject *species, *classes;
    int n, iStep, nSteps;
    const char *name;
    double dTime, dDelta;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "OOOOisiiddddd:Restart", const_cast<char **>(kwlist),
                &msr->arguments,&msr->specified,&species,&classes,&n,&name,
                &iStep,&nSteps,&dTime,&dDelta,&msr->dEcosmo,&msr->dUOld, &msr->dTimeOld ) )
        return NULL;

    // Create a vector of number of species
    species = PySequence_Fast(species,"species must be a list");
    int nSpecies = PySequence_Fast_GET_SIZE(species);
    std::vector<uint64_t> vecSpecies;
    vecSpecies.reserve(nSpecies);
    for (auto i=0; i < nSpecies; ++i) {
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
    for (auto i=0; i < nClasses; ++i) {
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

    msr->Restart(n, name, iStep, nSteps, dTime, dDelta);

    Py_RETURN_NONE;
}

/********** Data Structures: Domain, Trees **********/

static PyObject *
ppy_msr_DomainDecomp(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"rung",NULL};
    int iRung    = 0;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|i:DomainDecomp", const_cast<char **>(kwlist),
                &iRung ) )
        return NULL;
    self->msr->DomainDecomp(iRung);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_BuildTree(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"ewald","rung","active","marked",NULL};
    int bEwald = self->msr->param.bEwald;
    uint8_t uRungDT = 0; /* Zero rung means build a single tree */
    int bActive = 0;
    int bMarked = 0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|pBpp:BuildTree", const_cast<char **>(kwlist),
                &bEwald, &uRungDT, &bActive, &bMarked ) )
        return NULL;
    if (bMarked) self->msr->BuildTreeMarked();
    else if (uRungDT) {
        if (bActive) self->msr->BuildTreeActive(bEwald,uRungDT);
        else self->msr->BuildTreeFixed(bEwald,uRungDT);
    }
    else if (bActive) return PyErr_Format(PyExc_ValueError,"Building an active tree requires a valid rung");
    else self->msr->BuildTree(bEwald);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_Reorder(MSRINSTANCE *self, PyObject *args) {
    flush_std_files();
    self->msr->Reorder();
    Py_RETURN_NONE;
}

/********** Algorithms: Gravity **********/

static PyObject *
ppy_msr_Gravity(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    MSR *msr = self->msr;
    static char const *kwlist[]= {"time","delta","theta","rung","ewald","step","KickClose","KickOpen", "onlyMarked", NULL};
    double dTime = 0.0;
    double dDelta = 0.0;
    double dTheta = msr->param.dTheta;

    int bEwald = msr->param.bEwald;
    int iRungLo    = 0;
    int iRungHi    = MAX_RUNG;
    int iRoot1 = ROOT;
    int iRoot2 = 0;
    double dStep = 0.0;
    //double dTheta = msr->getParameterDouble("dTheta");
    int bKickOpen = 1;
    int bKickClose = 1;
    int onlyMarked = 0;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "d|ddipdppp:Gravity", const_cast<char **>(kwlist),
                &dTime, &dDelta, &dTheta, &iRungLo, &bEwald, &dStep, &bKickClose, &bKickOpen, &onlyMarked ) )
        return NULL;
    if (onlyMarked) iRoot2 = FIXROOT;

    SPHOptions SPHoptions = initializeSPHOptions(msr->param,msr->csm,dTime);
    SPHoptions.doGravity = 1;
    uint8_t uRungMax = msr->Gravity(iRungLo,iRungHi,iRoot1,iRoot2,dTime,dDelta,dStep,dTheta,bKickClose,bKickOpen,bEwald,
                                    msr->param.bGravStep, msr->param.nPartRhoLoc, msr->param.iTimeStepCrit, msr->param.nGroup, SPHoptions);
    return Py_BuildValue("i", uRungMax);
}

/********** Algorithms: Smooth **********/

static PyObject *
ppy_msr_Smooth(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"type","n","symmetric","time","delta","resmooth",NULL};
    int iSmoothType;
    int bSymmetric = 0;
    int bResmooth = 0;
    double dTime = 0.0;
    double dDelta = 0.0;
    int nSmooth = self->msr->param.nSmooth;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "i|ipddp:Smooth", const_cast<char **>(kwlist),
                &iSmoothType,&nSmooth,&bSymmetric, &dTime, &dDelta, &bResmooth ) )
        return NULL;
    if (bResmooth) self->msr->ReSmooth(dTime,dDelta,iSmoothType,bSymmetric);
    else self->msr->Smooth(dTime,dDelta,iSmoothType,bSymmetric,nSmooth);
    Py_RETURN_NONE;
}

/********** Analysis: grid management **********/

static PyObject *
ppy_msr_grid_prepare(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"grid",NULL};
    int iGrid;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "i:grid_prepare", const_cast<char **>(kwlist),
                &iGrid ) )
        return NULL;
    self->msr->GridCreateFFT(iGrid);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_grid_release(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {NULL};
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, ":grid_release", const_cast<char **>(kwlist)) )
        return NULL;
    self->msr->GridDeleteFFT();
    Py_RETURN_NONE;
}

/********** Analysis: Measure Bispectrum **********/

static PyObject *
ppy_msr_bispectrum_select(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","source","kmin","kmax",NULL};
    int target=0, source=1;
    double kmin, kmax;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "iidd:bispectrum_select", const_cast<char **>(kwlist),
                &target, &source, &kmin, &kmax ) )
        return NULL;
    self->msr->BispectrumSelect(target,source,kmin,kmax);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_bispectrum_normalize(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","kmin","kmax",NULL};
    int target=0, source=-1;
    double kmin, kmax;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "idd:bispectrum_select", const_cast<char **>(kwlist),
                &target, &kmin, &kmax ) )
        return NULL;
    // If source is negative, then we fill with 0 or 1.
    self->msr->BispectrumSelect(target,source,kmin,kmax);
    Py_RETURN_NONE;
}
static PyObject *
ppy_msr_bispectrum_calculate(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"grid1","grid2","grid3",NULL};
    int grid1, grid2, grid3;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "iii:bispectrum_calculate", const_cast<char **>(kwlist),
                &grid1, &grid2, &grid3 ) )
        return NULL;
    double result = self->msr->BispectrumCalculate(grid1,grid2,grid3);
    return Py_BuildValue("d", result);
}

/********** Analysis: Measure P(k) **********/

static PyObject *
ppy_msr_assign_mass(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","order","delta",NULL};
    int order=4, target=0;
    double delta = 0.0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|iid:assign_mass", const_cast<char **>(kwlist),
                &target, &order, &delta ) )
        return NULL;
    self->msr->AssignMass(order,target,delta);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_density_contrast(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","k",NULL};
    int target=0, k=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|ip:density_contrast", const_cast<char **>(kwlist),
                &target, &k ) )
        return NULL;
    self->msr->DensityContrast(target,k);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_interlace(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","source",NULL};
    int target=0, source=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "ii:interlace", const_cast<char **>(kwlist),
                &target, &source ) )
        return NULL;
    self->msr->Interlace(target,source);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_window_correction(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","order",NULL};
    int order=4, target=0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|ii:interlace", const_cast<char **>(kwlist),
                &target, &order ) )
        return NULL;
    self->msr->WindowCorrection(order,target);
    Py_RETURN_NONE;
}

static PyObject *
ppy_msr_add_linear_signal(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"target","seed","L","a","fixed","phase",NULL};
    int target, seed, fixed=0;
    float phase = 0;
    double L, a;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "iidd|pf:add_linear_signal", const_cast<char **>(kwlist),
                &target, &seed, &L, &a, &fixed, &phase ) )
        return NULL;
    self->msr->AddLinearSignal(target,seed,L,a,fixed,phase);
    Py_RETURN_NONE;
}


static PyObject *
ppy_msr_grid_bin_k(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"source", "bins",NULL};
    int nBins = -1;
    int iGrid;
    std::vector<float> fK,fPk;
    std::vector<uint64_t> nK;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "ii:grid_bin_k", const_cast<char **>(kwlist),
                &iGrid, &nBins ) )
        return NULL;
    fPk.resize(nBins+1);
    fK.resize(nBins+1);
    nK.resize(nBins+1);
    self->msr->GridBinK(nBins,iGrid,nK.data(),fK.data(),fPk.data());
    auto ListK = PyList_New( nBins+1 );
    auto ListPk = PyList_New( nBins+1 );
    auto ListNk = PyList_New( nBins+1 );
    for ( auto i=0; i<=nBins; i++ ) {
        PyList_SetItem(ListK,i,Py_BuildValue("f",fK[i]));
        PyList_SetItem(ListPk,i,Py_BuildValue("f",fPk[i]));
        PyList_SetItem(ListNk,i,Py_BuildValue("L",nK[i]));
    }
    return Py_BuildValue("(NNN)",ListK,ListPk,ListNk);
}

static PyObject *
ppy_msr_MeasurePk(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"grid","bins","a",NULL};
    double a = 1.0;
    int nBins = -1;
    int nGrid, i;
    std::vector<float> fK,fPk,fPkAll;
    std::vector<uint64_t> nK;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "i|id:MeasurePk", const_cast<char **>(kwlist),
                &nGrid, &nBins, &a ) )
        return NULL;
    if (nBins<0) nBins = nGrid/2;

    fPk.resize(nBins+1);
    fPkAll.resize(nBins+1);
    fK.resize(nBins+1);
    nK.resize(nBins+1);
    self->msr->MeasurePk(4,1,nGrid,a,nBins,nK.data(),fK.data(),fPk.data(),fPkAll.data());
    auto ListK = PyList_New( nBins+1 );
    auto ListPk = PyList_New( nBins+1 );
    auto ListPkAll = PyList_New( nBins+1 );
    auto ListNk = PyList_New( nBins+1 );
    for ( i=0; i<=nBins; i++ ) {
        PyList_SetItem(ListK,i,Py_BuildValue("f",fK[i]));
        PyList_SetItem(ListPk,i,Py_BuildValue("f",fPk[i]));
        PyList_SetItem(ListNk,i,Py_BuildValue("L",nK[i]));
        PyList_SetItem(ListPkAll,i,Py_BuildValue("f",fPkAll[i]));
    }
    return Py_BuildValue("(NNNN)",ListK,ListPk,ListNk,ListPkAll);
}

static PyObject *
ppy_msr_write_grid(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"name","source","k",NULL};
    const char *name;
    int source=0, k=0;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "s|ip:write_grid", const_cast<char **>(kwlist),
                &name, &source, &k ) )
        return NULL;
    self->msr->OutputGrid(name, k, source, 1);
    Py_RETURN_NONE;
}

/********** Analysis: Mark particles **********/

static PyObject *
ppy_msr_MarkSpecies(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"species","setIfTrue","clearIfFalse",NULL};
    PyObject *species;
    int setIfTrue=1, clearIfFalse=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "O|pp:MarkSpecies", const_cast<char **>(kwlist),
                &species, &setIfTrue, &clearIfFalse ) )
        return NULL;
    uint64_t mSpecies = 0;

    if (PyNumber_Check(species)) {
        if (auto o = PyNumber_Long(species)) {
            mSpecies |= 1 << PyLong_AsLongLong(o);
            Py_DECREF(o);
        }
    }
    else if (PySequence_Check(species)) {
        species = PySequence_Fast(species,"Expected a single species or a list of species");
        if (species==NULL) return NULL;
        int nSpecies = PySequence_Fast_GET_SIZE(species);
        for (auto i=0; i < nSpecies; ++i) {
            PyObject *item = PySequence_Fast_GET_ITEM(species, i);
            mSpecies |= 1 << PyNumber_AsSsize_t(item,NULL);
        }
        Py_DECREF(species); // PySequence_Fast creates a new reference
    }
    else return PyErr_Format(PyExc_TypeError,"Expected a single species or a list of species");
    auto n = self->msr->SelSpecies(mSpecies,setIfTrue,clearIfFalse);
    return Py_BuildValue("L",n);
}

static PyObject *
ppy_msr_MarkBox(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"center","size","setIfTrue","clearIfFalse",NULL};
    double center[3], size[3];
    int setIfTrue=1, clearIfFalse=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "(ddd)(ddd)|pp:MarkBox", const_cast<char **>(kwlist),
                &center[0], &center[1], &center[2], &size[0], &size[1], &size[2], &setIfTrue, &clearIfFalse ) )
        return NULL;
    auto n = self->msr->SelBox(center,size,setIfTrue,clearIfFalse);
    return Py_BuildValue("L",n);
}

static PyObject *
ppy_msr_MarkSphere(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"center","radius","setIfTrue","clearIfFalse",NULL};
    double center[3], radius;
    int setIfTrue=1, clearIfFalse=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "(ddd)d|pp:MarkSphere", const_cast<char **>(kwlist),
                &center[0], &center[1], &center[2], &radius, &setIfTrue, &clearIfFalse ) )
        return NULL;
    auto n = self->msr->SelSphere(center,radius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L",n);
}

static PyObject *
ppy_msr_MarkCylinder(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"point1","point2","radius","setIfTrue","clearIfFalse",NULL};
    double point1[3], point2[3], radius;
    int setIfTrue=1, clearIfFalse=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "(ddd)(ddd)d|pp:MarkCylinder", const_cast<char **>(kwlist),
                &point1[0], &point1[1], &point1[2], &point2[0], &point2[1], &point2[2], &radius, &setIfTrue, &clearIfFalse ) )
        return NULL;
    auto n = self->msr->SelCylinder(point1,point2,radius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L",n);
}

static PyObject *
ppy_msr_MarkBlackholes(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"setIfTrue","clearIfFalse",NULL};
    int setIfTrue=1, clearIfFalse=1;
    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|pp:MarkBlackholes", const_cast<char **>(kwlist),
                &setIfTrue, &clearIfFalse ) )
        return NULL;
    auto n = self->msr->SelBlackholes(setIfTrue,clearIfFalse);
    return Py_BuildValue("L",n);
}

/********** Analysis: Retrieve the values for all particles (DANGER! not memory friendly) **********/

#ifdef USE_NUMPY
static PyObject *
ppy_msr_GetArray(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    flush_std_files();
    static char const *kwlist[]= {"field","time","marked",NULL};
    PKD_FIELD field;
    double dTime = 1.0;
    int bMarked = 0;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "i|dp:GetArray", const_cast<char **>(kwlist),
                &field, &dTime, &bMarked ) )
        return NULL;
    npy_intp N[2];
    int typenum = NPY_FLOAT32;
    int iUnitSize = sizeof(float);
    N[0] = self->msr->N;
    N[1] = 1;
    switch (field) {
    case oPosition:     N[1]=3; typenum = NPY_FLOAT64; iUnitSize = sizeof(double);  break;
    case oAcceleration: N[1]=3; break;
    case oVelocity:     N[1]=3; break;
    case oPotential:    break;
    case oGroup:                 typenum = NPY_UINT32; iUnitSize = sizeof(uint32_t); break;
    case oMass:         break;
    case oSoft:         break;
    case oDensity:      break;
    case oBall:         break;
    case oParticleID:            typenum = NPY_UINT64; iUnitSize = sizeof(uint64_t); break;
    default: abort();
//  oSph, /* Sph structure */
//  oStar, /* Star structure */
//  oRelaxation,
//  oVelSmooth,
//  oRungDest, /* Destination processor for each rung */
    }
    if (bMarked) N[0] = self->msr->CountSelected();
    auto array = PyArray_SimpleNew(N[1]>1?2:1, N, typenum);
    auto data = reinterpret_cast<double *>(PyArray_DATA(reinterpret_cast<PyArrayObject *>(array)));
    self->msr->RecvArray(data,field,N[1]*iUnitSize,dTime,bMarked);
    return array;
}
#if 0
static PyObject *
ppy_msr_GetParticles(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]= {"time","marked",NULL};
    double dTime = 1.0;
    int bMarked = 0;

    if ( !PyArg_ParseTupleAndKeywords(
                args, kwobj, "|dp:GetArray", const_cast<char **>(kwlist),
                &dTime, &bMarked ) )
        return NULL;
    npy_intp N[2];
    int typenum = NPY_FLOAT32;
    int iUnitSize = sizeof(float);
    N[0] = self->msr->N;
    N[1] = 1;
    if (bMarked) N[0] = self->msr->CountSelected();
//    auto array = PyArray_SimpleNew(N[1]>1?2:1, N, typenum);
//    auto data = reinterpret_cast<double*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(array)));
//    self->msr->RecvArray(data,field,N[1]*iUnitSize,dTime,bMarked);
//    return array;
    Py_RETURN_NONE; // FOR NOW, FIX
}
#endif

#endif

/******************************************************************************\
*   MSR Object methods list
\******************************************************************************/

static PyMethodDef msr_methods[] = {
    {
        "setParameters", (PyCFunction)ppy_msr_setParameters, METH_VARARGS|METH_KEYWORDS,
        "Set global MSR parameters for future function calls"
    },

    {
        "add_analysis", (PyCFunction)ppy_msr_add_analysis, METH_VARARGS|METH_KEYWORDS,
        "Add an analysis callback hook"
    },
    {
        "ephemeral", (PyCFunction)ppy_msr_ephemeral, METH_VARARGS|METH_KEYWORDS,
        "Return an Ephemeral class for memory management"
    },

    {
        "simulate", (PyCFunction)ppy_msr_simulate, METH_VARARGS|METH_KEYWORDS,
        "Start a regular simulation"
    },

    {
        "GenerateIC", (PyCFunction)ppy_msr_GenerateIC, METH_VARARGS|METH_KEYWORDS,
        "Generate Initial Condition"
    },

    {
        "Load", (PyCFunction)ppy_msr_Load, METH_VARARGS|METH_KEYWORDS,
        "Load an input file"
    },
    {
        "Save", (PyCFunction)ppy_msr_Save, METH_VARARGS|METH_KEYWORDS,
        "Write a particle output"
    },
    {
        "WriteArray", (PyCFunction)ppy_msr_WriteArray, METH_VARARGS|METH_KEYWORDS,
        "Write an array of some feld"
    },
    {
        "WriteVector", (PyCFunction)ppy_msr_WriteVector, METH_VARARGS|METH_KEYWORDS,
        "Write a vector output"
    },
    {
        "write_grid", (PyCFunction)ppy_msr_write_grid, METH_VARARGS|METH_KEYWORDS,
        "Write a grid output"
    },

    {
        "Checkpoint", (PyCFunction)ppy_msr_Checkpoint, METH_VARARGS|METH_KEYWORDS,
        "Write a checkpoint"
    },
    {
        "Restart", (PyCFunction)ppy_msr_Restart, METH_VARARGS|METH_KEYWORDS,
        "Restart from a checkpoint"
    },

    {
        "DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
        "Reorder the particles by position"
    },
    {
        "BuildTree", (PyCFunction)ppy_msr_BuildTree, METH_VARARGS|METH_KEYWORDS,
        "Build the spatial tree"
    },
    {
        "Reorder", (PyCFunction)ppy_msr_Reorder, METH_NOARGS,
        "Reorders the particles by iOrder"
    },

    {
        "Gravity", (PyCFunction)ppy_msr_Gravity, METH_VARARGS|METH_KEYWORDS,
        "Calculate gravity"
    },
    {
        "Smooth", (PyCFunction)ppy_msr_Smooth, METH_VARARGS|METH_KEYWORDS,
        "Smooth"
    },

    {
        "grid_prepare", (PyCFunction)ppy_msr_grid_prepare, METH_VARARGS|METH_KEYWORDS,
        "Prepare for grid operations"
    },
    {
        "grid_release", (PyCFunction)ppy_msr_grid_release, METH_VARARGS|METH_KEYWORDS,
        "Release grid resources"
    },
    {
        "assign_mass", (PyCFunction)ppy_msr_assign_mass, METH_VARARGS|METH_KEYWORDS,
        "Assign particle mass to a grid"
    },
    {
        "density_contrast", (PyCFunction)ppy_msr_density_contrast, METH_VARARGS|METH_KEYWORDS,
        "Calculate density contrast (delta)"
    },
    {
        "interlace", (PyCFunction)ppy_msr_interlace, METH_VARARGS|METH_KEYWORDS,
        "Interlace two k-space delta fields"
    },
    {
        "window_correction", (PyCFunction)ppy_msr_window_correction, METH_VARARGS|METH_KEYWORDS,
        "Correct k-space grid with mass assignment window function"
    },
    {
        "grid_bin_k", (PyCFunction)ppy_msr_grid_bin_k, METH_VARARGS|METH_KEYWORDS,
        "Bin the Grid in k space"
    },
    {
        "add_linear_signal", (PyCFunction)ppy_msr_add_linear_signal, METH_VARARGS|METH_KEYWORDS,
        "Add the linear signal to the existing grid"
    },
    {
        "MeasurePk", (PyCFunction)ppy_msr_MeasurePk, METH_VARARGS|METH_KEYWORDS,
        "Measure the power spectrum"
    },
    {
        "bispectrum_select", (PyCFunction)ppy_msr_bispectrum_select, METH_VARARGS|METH_KEYWORDS,
        "Select a k-space shell for the bispectrum and copy it to a target grid"
    },
    {
        "bispectrum_normalize", (PyCFunction)ppy_msr_bispectrum_normalize, METH_VARARGS|METH_KEYWORDS,
        "Fill a k-space shell for the bispectrum with 1 inside and 0 outside"
    },
    {
        "bispectrum_calculate", (PyCFunction)ppy_msr_bispectrum_calculate, METH_VARARGS|METH_KEYWORDS,
        "Calculate the bispectrum for three grids"
    },

    {
        "MarkSpecies", (PyCFunction)ppy_msr_MarkSpecies, METH_VARARGS|METH_KEYWORDS,
        "Mark one or more species of particles"
    },
    {
        "MarkBox", (PyCFunction)ppy_msr_MarkBox, METH_VARARGS|METH_KEYWORDS,
        "Mark particles in a box"
    },
    {
        "MarkSphere", (PyCFunction)ppy_msr_MarkSphere, METH_VARARGS|METH_KEYWORDS,
        "Mark particles in a sphere"
    },
    {
        "MarkCylinder", (PyCFunction)ppy_msr_MarkCylinder, METH_VARARGS|METH_KEYWORDS,
        "Mark particles in a cylinder"
    },
    {
        "MarkBlackholes", (PyCFunction)ppy_msr_MarkBlackholes, METH_VARARGS|METH_KEYWORDS,
        "Marks all blackholes"
    },
#ifdef USE_NUMPY
    {
        "GetArray", (PyCFunction)ppy_msr_GetArray, METH_VARARGS|METH_KEYWORDS,
        "Get a complete array of the given value"
    },
#endif
    {NULL, NULL, 0, NULL}
};

/******************************************************************************\
*   MSR OBJECT Boilerplate code
\******************************************************************************/

static void msr_dealloc(MSRINSTANCE *self) {
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *msr_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    MSRINSTANCE *self = reinterpret_cast<MSRINSTANCE *>(type->tp_alloc(type, 0));
    if (self == NULL) { return NULL; }
    // Make a copy of the MSR object
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState *>(PyModule_GetState(msr_module));
    self->msr = moduleState->msr;
    return reinterpret_cast<PyObject *>(self);
}

static int msr_init(MSRINSTANCE *self, PyObject *args, PyObject *kwds) {
    return 0;
}

static PyMemberDef msr_members[] = {
    {NULL} /* Sentinel */
};

static PyObject *msr_get_parm(MSRINSTANCE *self, void *) {
    Py_INCREF(self->msr->arguments);
    return self->msr->arguments;
}

static PyObject *msr_get_spec(MSRINSTANCE *self, void *) {
    Py_INCREF(self->msr->specified);
    return self->msr->specified;
}

// This warning should be fixed in newer Python versions
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wwrite-strings"
static PyGetSetDef msr_getseters[] = {
    {"parm", (getter)msr_get_parm, NULL, "All parameters", NULL},
    {"spec", (getter)msr_get_spec, NULL, "If parameters were specified", NULL},
    {NULL} /* Sentinel */
};
#pragma GCC diagnostic pop

static PyTypeObject msrType = {
    PyVarObject_HEAD_INIT(NULL,0)
    MASTER_MODULE_NAME "." MASTER_TYPE_NAME, /*tp_name*/
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
static PyObject *initModuleMSR(void) {
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState *>(PyModule_GetState(msr_module));
    moduleState->bImported = true;
#ifdef USE_NUMPY
    import_array();
#endif
    //auto mainModule = PyImport_AddModule("__main__"); // Borrowed reference
    //PyObject *mainDict = PyModule_GetDict(mainModule);
    //setConstants(mainDict);
    if (PyType_Ready(&msrType) >= 0) { // Initialize "MSR" object as well.
        Py_INCREF(&msrType);
        PyModule_AddObject(msr_module, MASTER_TYPE_NAME, (PyObject *)&msrType);
    }
    PyType_Ready(&ephemeralType);
    return msr_module;
}


/******************************************************************************\
*   Analysis callback
\******************************************************************************/

void MSR::addAnalysis(PyObject *callback,MSRINSTANCE *msr, PyObject *memory) {
    analysis_callbacks.emplace_back(callback,msr,memory);
}

void MSR::runAnalysis(int iStep,double dTime) {
    PyObject *call_args = PyTuple_New(1);
    PyObject *kwargs = PyDict_New();
    PyObject *step, *time, *a = NULL;

    PyDict_SetItemString(kwargs, "step", step=PyLong_FromLong(iStep));
    PyDict_SetItemString(kwargs, "time", time=PyFloat_FromDouble(dTime));
    if (csm->val.bComove) PyDict_SetItemString(kwargs, "a", a=PyFloat_FromDouble(csmTime2Exp(csm,dTime)));

    for ( msr_analysis_callback &i : analysis_callbacks) {
        PyTuple_SetItem(call_args,0,reinterpret_cast<PyObject *>(i.msr));
        PyObject_Call(i.callback,call_args,kwargs);
        if (PyErr_Occurred()) PyErr_Print();
    }
    PyDict_Clear(kwargs);
    Py_DECREF(step);
    Py_DECREF(time);
    Py_XDECREF(a);
    Py_DECREF(kwargs);
    Py_DECREF(call_args);
}

/******************************************************************************\
*   Setup MSR using Python to parse parameters / enter analysis mode
\******************************************************************************/

bool MSR::wasParameterSpecified(const char *name) const {
    bool bSpecified = false;
    if (auto f = PyObject_GetAttrString(specified,name)) {
        bSpecified = PyObject_IsTrue(f)>0;
        Py_DECREF(f);
    }
    return bSpecified;
}

bool MSR::getParameterBoolean(const char *name) const {
    bool v = false;
    if (auto o = PyObject_GetAttrString(arguments,name)) {
        v = PyObject_IsTrue(o)>0;
        Py_DECREF(o);
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
    return v;
}
void MSR::setParameter(const char *name,bool v,int bSpecified) {
    auto o = v ? Py_True : Py_False;
    PyObject_SetAttrString(arguments,name,o);
    if (bSpecified) {
        Py_INCREF(Py_True);
        PyObject_SetAttrString(specified,name,Py_True);
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
}


double MSR::getParameterDouble(const char *name) const {
    double v = 0.0;
    if (auto n = PyObject_GetAttrString(arguments,name)) {
        if (auto o = PyNumber_Float(n)) {
            v = PyFloat_AsDouble(o);
            Py_DECREF(o);
        }
        Py_DECREF(n);
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
    return v;
}
void MSR::setParameter(const char *name,double v,int bSpecified) {
    auto o = PyFloat_FromDouble(v);
    //v = PyLong_FromLong(*(int *)pn->pValue);
    if (o) {
        PyObject_SetAttrString(arguments,name,o);
        Py_DECREF(o);
        if (bSpecified) {
            Py_INCREF(Py_True);
            PyObject_SetAttrString(specified,name,Py_True);
        }
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
}

long long MSR::getParameterLongLong(const char *name) const {
    long long v = 0;
    if (auto n = PyObject_GetAttrString(arguments,name)) {
        if (auto o = PyNumber_Long(n)) {
            v = PyLong_AsLongLong(o);
            Py_DECREF(o);
        }
        Py_DECREF(n);
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
    return v;
}
void MSR::setParameter(const char *name,long long v,int bSpecified) {
    auto o = PyLong_FromLongLong(v);
    if (o) {
        PyObject_SetAttrString(arguments,name,o);
        Py_DECREF(o);
        if (bSpecified) {
            Py_INCREF(Py_True);
            PyObject_SetAttrString(specified,name,Py_True);
        }
    }
    if (PyErr_Occurred()) { PyErr_Print(); abort(); }
}

extern "C" PyObject *PyInit_CSM(void);

int MSR::Python(int argc, char *argv[]) {
    PyImport_AppendInittab(MASTER_MODULE_NAME,initModuleMSR);
    PyImport_AppendInittab("CSM",PyInit_CSM);

    // I don't like this, but it works for pyenv. See:
    //   https://bugs.python.org/issue22213
    //   https://www.python.org/dev/peps/pep-0432/
    auto PYENV_VIRTUAL_ENV = getenv("PYENV_VIRTUAL_ENV");
    if (PYENV_VIRTUAL_ENV) {
        std::string path = PYENV_VIRTUAL_ENV;
        path += "/bin/python";
        Py_SetProgramName(Py_DecodeLocale(path.c_str(),NULL));
    }
    Py_InitializeEx(0);

    // Contruct the "MSR" context and module
    auto msr_module = PyModule_Create(&msrModule);
    PyState_AddModule(msr_module,&msrModule);
    Initialize();
    auto moduleState = reinterpret_cast<struct msrModuleState *>(PyModule_GetState(msr_module));
    moduleState->msr = this;
    moduleState->bImported = false; // If imported then we enter script mode

    // Convert program arguments to unicode
    auto wargv = new wchar_t *[argc];
    for (int i=0; i<argc; ++i) wargv[i] = Py_DecodeLocale(argv[i],NULL);
    PySys_SetArgv(argc, wargv);
    delete[] wargv;

    PyObject *main_module = PyImport_ImportModule("__main__");
    auto globals = PyModule_GetDict(main_module);
    auto locals = globals;
    PyDict_SetItemString(globals, "__builtins__",PyEval_GetBuiltins());

    // Parse the command line
    auto PARSE = PyModule_New("parse");
    PyModule_AddStringConstant(PARSE, "__file__", "parse.py");
    PyObject *localDict = PyModule_GetDict(PARSE);
    PyDict_SetItemString(localDict, "__builtins__", PyEval_GetBuiltins());
    PyObject *pyValue = PyRun_String(parse_py, Py_file_input, localDict, localDict);
    Py_XDECREF(pyValue);
    PyObject *parse = PyObject_GetAttrString(PARSE, "parse");
    if (!PyCallable_Check(parse)) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST be callable\n"); abort(); }
    PyObject *update = PyObject_GetAttrString(PARSE, "update");
    if (!PyCallable_Check(update)) { fprintf(stderr,"INTERNAL ERROR: parse.update() MUST be callable\n"); abort(); }
    PyObject *result= PyObject_CallObject(parse,NULL);
    if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
    // Retrieve the results
    int n = PyTuple_Size(result);
    if (n!=2) { fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST return a tuple\n"); abort(); }
    Py_XDECREF(arguments);
    arguments = PyTuple_GetItem(result,0); /* Values of each parameter */
    Py_INCREF(arguments);
    Py_XDECREF(specified);
    specified = PyTuple_GetItem(result,1); /* If it was explicitely specified */
    Py_INCREF(specified);
    Py_DECREF(result);
    PyObject *script = PyObject_GetAttrString(arguments,"script");

    ppy2prm(prm,arguments,specified); // Update the pkdgrav parameter state
    bVDetails = getParameterBoolean("bVDetails");

    // If a script was specified then we run it.
    if (script != Py_None) {
        char *filename;
        if (PyUnicode_Check(script)) {
            PyObject *ascii = PyUnicode_AsASCIIString(script);
            filename = PyBytes_AsString(ascii);
            Py_DECREF(ascii);
        }
        else { fprintf(stderr,"INTERNAL ERROR: script filename is invalid\n"); abort(); }
        FILE *fp = fopen(filename,"r");
        if (fp == NULL) { perror(filename); exit(errno); }
        auto s = PyRun_FileEx(fp,filename,Py_file_input,globals,locals,1); // fp is closed on return
        Py_XDECREF(s);
        if (PyErr_Occurred()) {
            int rc = 1;
            if (PyErr_ExceptionMatches(PyExc_SystemExit)) {
                PyObject *etype, *evalue, *etrace;
                PyErr_Fetch(&etype, &evalue, &etrace);
                if (auto o = PyNumber_Long(evalue)) {
                    rc = PyLong_AsLong(o);
                    Py_DECREF(o);
                }
                Py_DECREF(etype);
                Py_DECREF(evalue);
                Py_DECREF(etrace);
            }
            else PyErr_Print();
            return rc;
        }
    }
    Py_DECREF(script);

    // If "MASTER" was imported then we are done -- the script should have done its job
    if (!moduleState->bImported) { // We must prepare for a normal legacy execution
        PyObject *args = PyTuple_New(3);
        PyTuple_SetItem(args,0,locals);
        PyTuple_SetItem(args,1,arguments);
        PyTuple_SetItem(args,2,specified);
        PyObject_CallObject(update,args); // Copy and variables into the arguments Namespace
        if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
        ppy2prm(prm,arguments,specified); // Update the pkdgrav parameter state
        bVDetails = getParameterBoolean("bVDetails");
    }

    return moduleState->bImported ? 0 : -1;
}