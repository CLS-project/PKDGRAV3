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
		if (f) { pn->bArg = PyObject_IsTrue(f)>0; Py_DECREF(f); }
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
*   MSR Module methods
\******************************************************************************/

struct MSRINSTANCE {
    PyObject_HEAD
    MSR *msr;
    };

/********** Initial Condition Generation **********/

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
    double dTime = self->msr->GenerateIC();
    return Py_BuildValue("d", dTime );
    }

/********** File I/O **********/

static PyObject *
ppy_msr_Load(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"filename",NULL};
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
    static char const *kwlist[]={"name","time",NULL};
    const char *fname;
    double dTime = 1.0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|d:Save", const_cast<char **>(kwlist),
	     &fname,&dTime ) )
	return NULL;
    self->msr->Write(fname,dTime,0);
    Py_RETURN_NONE;
    }

/********** Checkpoint I/O **********/

static PyObject *
ppy_msr_Checkpoint(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={/*"name",*/"step","time",NULL};
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
    MSR *msr = self->msr;
    static char const *kwlist[]={"arguments","specified","species","classes","n","name","step","steps","time","delta","E","U", "Utime", NULL};
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

    msr->Restart(n, name, iStep, nSteps, dTime, dDelta);

    Py_RETURN_NONE;
    }

/********** Data Structures: Domain, Trees **********/

static PyObject *
ppy_msr_DomainDecomp(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"rung",NULL};
    int iRung    = 0;
    int bOthers  = 0;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|i:DomainDecomp", const_cast<char **>(kwlist),
	     &iRung ) )
	return NULL;
    self->msr->DomainDecomp(iRung);
    Py_RETURN_NONE;
    }

static PyObject *
ppy_msr_BuildTree(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"ewald","rung","active",NULL};
    int bEwald = self->msr->param.bEwald;
    uint8_t uRungDT = 0; /* Zero rung means build a single tree */
    int bActive = 0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|iBi:BuildTree", const_cast<char **>(kwlist),
	     &bEwald, &uRungDT, &bActive ) )
	return NULL;
    if (uRungDT) {
	if (bActive) self->msr->BuildTreeActive(bEwald,uRungDT);
	else self->msr->BuildTreeFixed(bEwald,uRungDT);
	}
    else if (bActive) return PyErr_Format(PyExc_TypeError,"Building an active tree requires a valid rung");
    else self->msr->BuildTree(bEwald);
    Py_RETURN_NONE;
    }

static PyObject *
ppy_msr_Reorder(MSRINSTANCE *self, PyObject *args) {
    self->msr->Reorder();
    Py_RETURN_NONE;
    }

/********** Algorithms: Gravity **********/

static PyObject *
ppy_msr_Gravity(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    MSR *msr = self->msr;
    static char const *kwlist[]={"time","delta","rung","ewald","step","KickClose","KickOpen", NULL};
    double dTime = 0.0;
    double dDelta = 0.0;
    uint64_t nActive;

    int bEwald = msr->param.bEwald;
    int iRungLo    = 0;
    int iRungHi    = MAX_RUNG;
    int iRoot1 = ROOT;
    int iRoot2 = 0;
    double dStep = 0.0;
    //double dTheta = msr->getParameterDouble("dTheta");
    int bKickOpen = 1;
    int bKickClose = 1;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "d|dipdpp:Gravity", const_cast<char **>(kwlist),
	     &dTime, &dDelta, &iRungLo, &bEwald, &dStep, &bKickClose, &bKickOpen ) )
	return NULL;
    uint8_t uRungMax = msr->Gravity(iRungLo,iRungHi,iRoot1,iRoot2,dTime,dDelta,dStep,bKickClose,bKickOpen,bEwald,
	msr->param.bGravStep, msr->param.nPartRhoLoc, msr->param.iTimeStepCrit, msr->param.nGroup);
    return Py_BuildValue("i", uRungMax);
    }

/********** Algorithms: Smooth **********/

static PyObject *
ppy_msr_Smooth(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"type","n","symmetric","time","delta","resmooth",NULL};
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

/********** Analysis: Measure P(k) **********/

static PyObject *
ppy_msr_MeasurePk(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"grid","bins","a",NULL};
    double a = 1.0;
    int nBins = -1;
    int nGrid, i;
    std::vector<float> fK,fPk,fPkAll;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|ia:MeasurePk", const_cast<char **>(kwlist),
	     &nGrid ) )
	return NULL;
    if (nBins<0) nBins = nGrid/2;

    fPk.resize(nBins+1);
    fPkAll.resize(nBins+1);
    fK.resize(nBins+1);
    self->msr->MeasurePk(4,1,nGrid,a,nGrid/2,NULL,fK.data(),fPk.data(),fPkAll.data());
    auto ListK = PyList_New( nBins+1 );
    auto ListPk = PyList_New( nBins+1 );
    auto ListPkAll = PyList_New( nBins+1 );
    for( i=0; i<=nBins; i++ ) {
	PyList_SetItem(ListK,i,Py_BuildValue("f",fK[i]));
	PyList_SetItem(ListPk,i,Py_BuildValue("f",fPk[i]));
	PyList_SetItem(ListPkAll,i,Py_BuildValue("f",fPkAll[i]));
	}
    return Py_BuildValue("(NNN)",ListK,ListPk,ListPkAll);
    }


/********** Analysis: Retrieve the values for all particles (DANGER! not memory friendly) **********/

#ifdef USE_NUMPY
static PyObject *
ppy_msr_GetArray(MSRINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char const *kwlist[]={"field","time",NULL};
    PKD_FIELD field;
    double dTime = 1.0;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|d:GetArray", const_cast<char **>(kwlist),
	     &field, &dTime ) )
	return NULL;
    npy_intp N[2];
    int typenum = NPY_FLOAT32;
    int iUnitSize = sizeof(float);
    N[0] = self->msr->N;
    N[1] = 1;
    switch(field) {
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
//	oSph, /* Sph structure */
//	oStar, /* Star structure */
//	oRelaxation,
//	oVelSmooth,
//	oRungDest, /* Destination processor for each rung */
	}
    auto array = PyArray_SimpleNew(N[1]>1?2:1, N, typenum);
    auto data = reinterpret_cast<double*>(PyArray_DATA(reinterpret_cast<PyArrayObject*>(array)));
    self->msr->RecvArray(data,field,N[1]*iUnitSize,dTime);
    return array;
    }
#endif

/******************************************************************************\
*   MSR Module methods list
\******************************************************************************/

static PyMethodDef msr_methods[] = {
    {"GenerateIC", (PyCFunction)ppy_msr_GenerateIC, METH_VARARGS|METH_KEYWORDS,
     "Generate Initial Condition"},

    {"Load", (PyCFunction)ppy_msr_Load, METH_VARARGS|METH_KEYWORDS,
     "Load an input file"},
    {"Save", (PyCFunction)ppy_msr_Save, METH_VARARGS|METH_KEYWORDS,
     "Write a particle output"},

    {"Checkpoint", (PyCFunction)ppy_msr_Checkpoint, METH_VARARGS|METH_KEYWORDS,
     "Write a checkpoint"},
    {"Restart", (PyCFunction)ppy_msr_Restart, METH_VARARGS|METH_KEYWORDS,
     "Restart from a checkpoint"},

    {"DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
    {"BuildTree", (PyCFunction)ppy_msr_BuildTree, METH_VARARGS|METH_KEYWORDS,
     "Build the spatial tree"},
    {"Reorder", (PyCFunction)ppy_msr_Reorder, METH_NOARGS,
     "Reorders the particles by iOrder"},

    {"Gravity", (PyCFunction)ppy_msr_Gravity, METH_VARARGS|METH_KEYWORDS,
     "Calculate gravity"},
    {"Smooth", (PyCFunction)ppy_msr_Smooth, METH_VARARGS|METH_KEYWORDS,
     "Smooth"},

    {"MeasurePk", (PyCFunction)ppy_msr_MeasurePk, METH_VARARGS|METH_KEYWORDS,
     "Measure the power spectrum"},

#ifdef USE_NUMPY
    {"GetArray", (PyCFunction)ppy_msr_GetArray, METH_VARARGS|METH_KEYWORDS,
     "Get a complete array of the given value"},
#endif
    {NULL, NULL, 0, NULL}
};

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

static PyObject *msr_get_parm(MSRINSTANCE *self, void *) {
    Py_INCREF(self->msr->arguments);
    return self->msr->arguments;
    }

static PyObject *msr_get_spec(MSRINSTANCE *self, void *) {
    Py_INCREF(self->msr->specified);
    return self->msr->specified;
    }

static PyGetSetDef msr_getseters[] = {
	{"parm", (getter)msr_get_parm, NULL, "All parameters", NULL},
	{"spec", (getter)msr_get_spec, NULL, "If parameters were specified", NULL},
	{NULL} /* Sentinel */
    };

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
static PyObject * initModuleMSR(void) {
    auto msr_module = PyState_FindModule(&msrModule); // We created this already
    auto moduleState = reinterpret_cast<struct msrModuleState*>(PyModule_GetState(msr_module));
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
    return msr_module;
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
    auto moduleState = reinterpret_cast<struct msrModuleState*>(PyModule_GetState(msr_module));
    moduleState->msr = this;
    moduleState->bImported = false; // If imported then we enter script mode

    // Convert program arguments to unicode
    auto wargv = new wchar_t *[argc];
    for(int i=0; i<argc; ++i) wargv[i] = Py_DecodeLocale(argv[i],NULL);
    PySys_SetArgv(argc, wargv);

    PyObject * main_module = PyImport_ImportModule("__main__");
    auto globals = PyModule_GetDict(main_module);
    auto locals = globals;
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
    arguments = PyTuple_GetItem(result,0); /* Values of each parameter */
    specified = PyTuple_GetItem(result,1); /* If it was explicitely specified */
    PyObject *script = PyObject_GetAttrString(arguments,"script");

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
	auto s = PyRun_FileEx(fp,filename,Py_file_input,globals,locals,1); // fp is closed on return
	if (PyErr_Occurred()) { PyErr_Print(); exit(1); }
	Py_DECREF(s);
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
	ppy2prm(prm,arguments,specified); // Update the pkdgrav parameter state
	}

    return moduleState->bImported;
    }