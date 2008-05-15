#include <Python.h>
#include <marshal.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include "master.h"
#include "python.h"

/**********************************************************************\
 ** GLOBAL variables...  Sadly Python is not reentrant
\**********************************************************************/
static MSR ppy_msr;

typedef struct {
    PyObject *module;
    } ppyCtx;

static ppyCtx *global_ppy = NULL;

/**********************************************************************\
 ** MASTER Interface (msr functions)
\**********************************************************************/

static PyObject *
ppy_msr_SelSrcAll(PyObject *self, PyObject *args) {
    msrSelSrcAll(ppy_msr);
    return Py_BuildValue("L", ppy_msr->N);
}

static PyObject *
ppy_msr_SelDstAll(PyObject *self, PyObject *args) {
    msrSelDstAll(ppy_msr);
    return Py_BuildValue("L", ppy_msr->N);
}

static PyObject *
ppy_msr_SelSrcMass(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinMass", "MaxMass", "setIfTrue", "clearIfFalse", NULL };
    double dMinMass, dMaxMass;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelSrcMass", kwlist,
	     &dMinMass, &dMaxMass, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcMass(ppy_msr,dMinMass,dMaxMass,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstMass(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinMass", "MaxMass", "setIfTrue", "clearIfFalse", NULL };
    double dMinMass, dMaxMass;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelDstMass", kwlist,
	     &dMinMass, &dMaxMass, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstMass(ppy_msr,dMinMass,dMaxMass,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelSrcBox(PyObject *self, PyObject *args) {
    double dCenter[3], dSize[3];
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)|ii:SelSrcBox",
	     dCenter+0, dCenter+1, dCenter+2, dSize+0, dSize+1, dSize+2, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcBox(ppy_msr,dCenter,dSize,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}
static PyObject *
ppy_msr_SelDstBox(PyObject *self, PyObject *args) {
    double dCenter[3], dSize[3];
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)|ii:SelDstBox",
	     dCenter+0, dCenter+1, dCenter+2, dSize+0, dSize+1, dSize+2, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstBox(ppy_msr,dCenter,dSize,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelSrcSphere(PyObject *self, PyObject *args) {
    double r[3], dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)d|ii:SelSrcSphere",
	     r+0, r+1, r+2, &dRadius, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcSphere(ppy_msr,r,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstSphere(PyObject *self, PyObject *args) {
    double r[3], dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)d|ii:SelDstSphere",
	     r+0, r+1, r+2, &dRadius, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstSphere(ppy_msr,r,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelSrcCylinder(PyObject *self, PyObject *args) {
    double dP1[3],dP2[3],dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)d|ii:SelSrcSphere",
	     dP1+0, dP1+1, dP1+2, dP2+0, dP2+1, dP2+2, &dRadius,
	     &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcCylinder(ppy_msr,dP1,dP2,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstCylinder(PyObject *self, PyObject *args) {
    double dP1[3],dP2[3],dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)(ddd)d|ii:SelDstSphere",
	     dP1+0, dP1+1, dP1+2, dP2+0, dP2+1, dP2+2, &dRadius,
	     &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstCylinder(ppy_msr,dP1,dP2,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

#if 0
static PyObject *
ppy_msr_SelSrc(PyObject *self, PyObject *args) {
    uint64_t nSelected = 0;
    PyObject *object;
    if ( !PyArg_ParseTuple(
	     args, "O:SelSrc",
	     &object) )
	return NULL;
    if ( !PyCallable_Check(object) ) {
	PyErr_SetString(PyExc_TypeError, "SelSrc argument is not callable");
	return NULL;
	}
#if 0
    PyObject *o = PyMarshal_WriteObjectToString(object,Py_MARSHAL_VERSION);
    assert( o!=NULL );

    if ( !PyString_Check(o) ) {
	PyErr_SetString(PyExc_TypeError, "PyMarshal did not return a string");
	return NULL;
	}
    int iSize = PyString_GET_SIZE(o);
    const char *szData = PyString_AS_STRING(o);

    printf( "Object size: %d\n", iSize );
    Py_DECREF(o);
#endif
    return Py_BuildValue("L", nSelected);
}
#endif

static PyObject *
ppy_msr_DeepestPotential(PyObject *self, PyObject *args) {
    double r[3];
    float fPot;
    msrDeepestPot(ppy_msr,r,&fPot);
    return Py_BuildValue("((ddd)f)", r[0], r[1], r[2], fPot);
}

static PyObject *
ppy_msr_Profile(PyObject *self, PyObject *args) {
    double r[3], dMinR, dMaxR, dMassEncl, dRho;
    int nBins, i;
    const PROFILEBIN *pBins;
    PyObject *List;

    nBins = 100;
    if ( !PyArg_ParseTuple(
	     args, "(ddd)dd|i:Profile",
	     r+0, r+1, r+2, &dMinR, &dMaxR, &nBins ) )
	return NULL;
    pBins = msrProfile(ppy_msr,r,dMinR,dMaxR,nBins);
    List = PyList_New( nBins );
    dMassEncl = 0.0;
    for(i=0; i<nBins; i++ ) {
	dMassEncl += pBins[i].dMassInBin;
	dRho = pBins[i].dMassInBin/pBins[i].dVolume;
	PyList_SetItem(List,i,Py_BuildValue("(dLdd)",pBins[i].dRadius, pBins[i].nParticles, dRho, dMassEncl));
	}
    msrDeleteProfile(ppy_msr);
    return List;
}

static PyObject *
ppy_msr_Reorder(PyObject *self, PyObject *args) {
    msrReorder(ppy_msr);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_DomainDecomp(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Rung","SplitVA",NULL};
    int iRung    = 0;
    int bSplitVA = 0;
    int bGreater = 1;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|ii:DomainDecomp", kwlist,
	     &iRung, &bSplitVA ) )
	return NULL;
    msrDomainDecomp(ppy_msr,iRung,bGreater,bSplitVA);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_BuildTree(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time","Ewald",NULL};
    double dTime = 0.0;
    int bEwald = ppy_msr->param.bEwald;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|di:BuildTree", kwlist,
	     &dTime, &bEwald ) )
	return NULL;
    msrBuildTree(ppy_msr,dTime,bEwald);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Fof(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"FOFsDone","Symmetric","Time",NULL};
    double dExp;
    double dTime = 0.0;
    int nFOFsDone = 0;
    int bSymmetric = 0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|iid:BuildTree", kwlist,
	     &nFOFsDone,& bSymmetric, &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrFof(ppy_msr,nFOFsDone,SMX_FOF,bSymmetric,dExp);
    msrGroupMerge(ppy_msr,dExp);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
ppy_msr_GroupProfiles(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"FOFsDone","Symmetric","Time",NULL};
    double dExp;
    double dTime = 0.0;
    int nFOFsDone = 0;
    int bSymmetric = 0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|iid:BuildTree", kwlist,
	     &nFOFsDone,& bSymmetric, &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrGroupProfiles(ppy_msr,nFOFsDone,SMX_FOF,bSymmetric,dExp);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Write(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Checkpoint","Time",NULL};
    double dTime = 0.0;
    int bCheckpoint = 0;
    const char *fname;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|id:Write", kwlist,
	     &fname,&bCheckpoint, &dTime ) )
	return NULL;
    msrWrite(ppy_msr,fname,dTime,bCheckpoint);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyMethodDef ppy_msr_methods[] = {
/*
    {"SelSrc", ppy_msr_SelSrc, METH_VARARGS,
     "Selects source particles based on a supplied function"},
*/
    {"SelSrcAll", ppy_msr_SelSrcAll, METH_NOARGS,
     "Selects all particles as operation source."},
    {"SelDstAll", ppy_msr_SelDstAll, METH_NOARGS,
     "Selects all particles as operation destination."},
    {"SelSrcMass", (PyCFunction)ppy_msr_SelSrcMass, METH_VARARGS|METH_KEYWORDS,
     "Selects source particles with a specific mass range."},
    {"SelDstMass", (PyCFunction)ppy_msr_SelDstMass, METH_VARARGS|METH_KEYWORDS,
     "Selects destination particles with a specific mass range."},
    {"SelSrcBox", ppy_msr_SelSrcBox, METH_VARARGS,
     "Selects source particles inside a given box."},
    {"SelDstBox", ppy_msr_SelDstBox, METH_VARARGS,
     "Selects destination particles inside a given box."},
    {"SelSrcSphere", ppy_msr_SelSrcSphere, METH_VARARGS,
     "Selects source particles inside a given sphere."},
    {"SelDstSphere", ppy_msr_SelDstSphere, METH_VARARGS,
     "Selects destination particles inside a given sphere."},
    {"SelSrcCylinder", ppy_msr_SelSrcCylinder, METH_VARARGS,
     "Selects source particles inside a given cylinder."},
    {"SelDstCylinder", ppy_msr_SelDstCylinder, METH_VARARGS,
     "Selects destination particles inside a given cylinder."},

    {"DeepestPotential", ppy_msr_DeepestPotential, METH_NOARGS,
     "Finds the most bound particle (deepest potential)"},
    {"Profile", ppy_msr_Profile, METH_VARARGS,
     "Generate a density profile"},

    {"Reorder", ppy_msr_Reorder, METH_NOARGS,
     "Reorders the particles by iOrder"},
    {"DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
    {"BuildTree", (PyCFunction)ppy_msr_BuildTree, METH_VARARGS|METH_KEYWORDS,
     "Build the spatial tree"},
    {"Fof", (PyCFunction)ppy_msr_Fof, METH_VARARGS|METH_KEYWORDS,
     "Friends of Friends"},
    {"GroupProfiles", (PyCFunction)ppy_msr_GroupProfiles, METH_VARARGS|METH_KEYWORDS,
     "Group Profiles"},
    {"Write", (PyCFunction)ppy_msr_Write, METH_VARARGS|METH_KEYWORDS,
     "Write source particles to a file"},

    {NULL, NULL, 0, NULL}
};

/**********************************************************************\
 ** Parallel Python (ppy) setup
\**********************************************************************/

void ppyInitialize(PPY *pvppy, MSR msr, double dTime) {
    ppyCtx *ppy;
    PyObject *dict;

    ppy = malloc(sizeof(ppyCtx));
    assert(ppy!=NULL);
    *pvppy = global_ppy = ppy;

    ppy_msr = msr;
    Py_Initialize();
    ppy->module = Py_InitModule( "msr", ppy_msr_methods );

    dict = PyModule_GetDict(ppy->module);
    PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",dTime));
    }

void ppyFinish(PPY vppy) {
    ppyCtx *ppy = (ppyCtx *)vppy;
    Py_Finalize();
    free(ppy);
    global_ppy = NULL;
    }

void ppyRunScript(PPY vppy,const char *achFilename) {
    ppyCtx *ppy = (ppyCtx *)vppy;
    FILE *fp;
    PyObject *dict;

    assert(Py_IsInitialized());

    // Set parameters for easy access
    dict = PyModule_GetDict(ppy->module);

    PyDict_SetItemString(dict, "achOutName", Py_BuildValue("s",ppy_msr->param.achOutName));





    printf("---------------------------------------"
	   "---------------------------------------\n"
	   "Running Python Script %s\n"
	   "---------------------------------------"
	   "---------------------------------------\n",
	   achFilename );
    fp = fopen(achFilename,"r");
    PyRun_SimpleFile(fp,achFilename);
    fclose(fp);
    printf("---------------------------------------"
	   "---------------------------------------\n" );

    }
