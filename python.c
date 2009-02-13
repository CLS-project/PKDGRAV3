#include <Python.h>
#include <marshal.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include "master.h"
#include "outtype.h"
#include "intype.h"
#include "python.h"

const char *python_c_module_id = "$Id$";
const char *python_h_module_id = PYTHON_H_MODULE_ID;

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
ppy_msr_SelSrcPhaseDensity(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinDensity", "MaxDensity", "setIfTrue", "clearIfFalse", NULL };
    double dMinDensity, dMaxDensity;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelSrcPhaseDensity", kwlist,
	     &dMinDensity, &dMaxDensity, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcPhaseDensity(ppy_msr,dMinDensity,dMaxDensity,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstPhaseDensity(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinDensity", "MaxDensity", "setIfTrue", "clearIfFalse", NULL };
    double dMinDensity, dMaxDensity;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelDstPhaseDensity", kwlist,
	     &dMinDensity, &dMaxDensity, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstPhaseDensity(ppy_msr,dMinDensity,dMaxDensity,setIfTrue,clearIfFalse);
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

static PyObject *
ppy_msr_SelSrcById(PyObject *self, PyObject *args) {
    uint64_t idStart,idEnd;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "LL|ii:SelSrcById",
	     &idStart, &idEnd, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSrcById(ppy_msr,idStart,idEnd,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstById(PyObject *self, PyObject *args) {
    uint64_t idStart,idEnd;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    if ( !PyArg_ParseTuple(
	     args, "ll|ii:SelDstById",
	     &idStart, &idEnd, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelDstById(ppy_msr,idStart,idEnd,setIfTrue,clearIfFalse);
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
    double r[3], dMinR, dLogR, dMaxR, dMassEncl, dRho, dVel;
    int nBins, nPerBin, nAccuracy, i;
    const PROFILEBIN *pBins;
    PyObject *List;
    double radius_old;

    nAccuracy = 2; /* +/- 2 particles per bin */
    nBins = 200;
    nPerBin = 0;
    if ( !PyArg_ParseTuple(
	     args, "(ddd)dddi|ii:Profile",
	     r+0, r+1, r+2, &dMinR, &dLogR, &dMaxR, &nPerBin, &nBins, &nAccuracy ) )
	return NULL;
    msrProfile(ppy_msr,&pBins,&nBins,r,dMinR,dLogR,dMaxR,nPerBin,nBins,nAccuracy);

    List = PyList_New( nBins );
    assert( List !=NULL );
    dMassEncl = 0.0;

    radius_old = 0.0;
    for(i=0; i<nBins; i++ ) {
	PyObject *tuple;
	double ang, ang_theta, ang_phi, vel_circ;
	double radius_mean;

	dMassEncl += pBins[i].dMassInBin;
	assert(pBins[i].dVolume>0.0);
	dRho = pBins[i].dMassInBin/pBins[i].dVolume;
	dVel = sqrt(dMassEncl/pBins[i].dRadius);

	radius_mean = (pBins[i].dRadius + radius_old ) * 0.5;
	radius_old = pBins[i].dRadius;
	ang = sqrt(dot_product(pBins[i].L,pBins[i].L));
	if(ang > 0.0) ang_theta = 180.*acos(pBins[i].L[2]/ang)/M_PI;
	else ang_theta = 0.0;
	ang_phi = 180.*atan2(pBins[i].L[1],pBins[i].L[0])/M_PI ;
	vel_circ = ang / radius_mean ;
	tuple = Py_BuildValue("(dLdddddddddd)",
			      pBins[i].dRadius,
			      pBins[i].nParticles, dRho, dMassEncl,
			      dVel,pBins[i].vel_radial,pBins[i].vel_radial_sigma,
			      vel_circ, pBins[i].vel_tang_sigma, ang, ang_theta, ang_phi );
	assert(tuple != NULL);
	assert( PyList_SetItem(List,i,tuple) >= 0 );
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
ppy_msr_Gravity(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time","Ewald",NULL};
    double dTime = 0.0;
    int bEwald = ppy_msr->param.bEwald;
    PyObject *v, *dict;
    uint64_t nActive;
    int iSec = 0;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|di:Gravity", kwlist,
	     &dTime, &bEwald ) )
	return NULL;

    msrGravity(ppy_msr,0,MAX_RUNG,dTime,ppy_msr->param.iStartStep,bEwald,&iSec,&nActive);

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
ppy_msr_Smooth(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"iSmoothType","bSymmetric","dTime",NULL};
    PyObject *v, *dict;
    int iSmoothType;
    int bSymmetric = 0;
    double dTime;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|id:Smooth", kwlist,
	     &iSmoothType,&bSymmetric, &dTime ) )
	return NULL;

    msrSmooth(ppy_msr,dTime,iSmoothType,bSymmetric);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_ReSmooth(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"iSmoothType","bSymmetric","dTime",NULL};
    PyObject *v, *dict;
    int iSmoothType;
    int bSymmetric = 0;
    double dTime;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|id:ReSmooth", kwlist,
	     &iSmoothType,&bSymmetric, &dTime ) )
	return NULL;
    msrReSmooth(ppy_msr,dTime,iSmoothType,bSymmetric);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
ppy_msr_PeakVc(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Centers",NULL};
    PyObject *list, *item;
    int N, i;
    struct inPeakVc *in;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "O:PeakVc", kwlist,
	     &list ) )
	return NULL;

    N = PyList_Size(list);
    in = malloc(sizeof(struct inPeakVc)*N);
    assert(in!=NULL);

    for( i=0; i<N; i++ ) {
	item = PyList_GetItem(list,i);
	if ( !PyArg_ParseTuple(
		 item, "ddd:PeakVc",
		 &in[i].dCenter[0], &in[i].dCenter[1], &in[i].dCenter[2] ) )
	    return NULL;
	}
    msrPeakVc(ppy_msr,N,in);

    free(in);

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
	     args, kwobj, "|iid:GroupProfiles", kwlist,
	     &nFOFsDone,& bSymmetric, &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrGroupProfiles(ppy_msr,nFOFsDone,SMX_FOF,bSymmetric,dExp);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Load(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Type",NULL};
    const char *fname;
    double dTime;
    int iType = IN_TIPSY_STD;
    PyObject *dict, *v;

    dict = PyModule_GetDict(global_ppy->module);

    if ( (v = PyDict_GetItemString(dict, "bMemAcceleration")) != NULL )
	ppy_msr->param.bMemAcceleration = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemVelocity")) != NULL )
	ppy_msr->param.bMemVelocity = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemPotential")) != NULL )
	ppy_msr->param.bMemPotential = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemGroups")) != NULL )
	ppy_msr->param.bMemGroups = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemMass")) != NULL )
	ppy_msr->param.bMemMass = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemSoft")) != NULL )
	ppy_msr->param.bMemSoft = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemHermite")) != NULL )
	ppy_msr->param.bMemHermite = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemRelaxation")) != NULL )
	ppy_msr->param.bMemRelaxation = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemVelSmooth")) != NULL )
	ppy_msr->param.bMemVelSmooth = PyInt_AsLong(v);

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|i:Load", kwlist,
	     &fname, &iType ) )
	return NULL;
    switch(iType) {
    case IN_TIPSY_STD:
    case IN_TIPSY_DBL:
    case IN_TIPSY_NAT:
	dTime = msrRead(ppy_msr,fname);
	msrInitStep(ppy_msr);
	PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",dTime));
	break;

    case IN_SRC_MARK:
    case IN_DST_MARK:
	break;

    default:
	assert(0);
	}



    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Save(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Type","Time",NULL};
    const char *fname;
    int iType = 0;
    double dTime = 0.0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|id:Save", kwlist,
	     &fname,&iType, &dTime ) )
	return NULL;
    switch(iType) {
    case OUT_TIPSY_STD:
    case OUT_TIPSY_DBL:
	msrWrite(ppy_msr,fname,dTime,iType==OUT_TIPSY_DBL);
	break;

    case OUT_IORDER_ARRAY:
    case OUT_COLOR_ARRAY:
    case OUT_DENSITY_ARRAY:
    case OUT_POT_ARRAY:
    case OUT_AMAG_ARRAY:
    case OUT_IMASS_ARRAY:
    case OUT_RUNG_ARRAY:
    case OUT_DIVV_ARRAY:
    case OUT_VELDISP2_ARRAY:
    case OUT_VELDISP_ARRAY:
    case OUT_PHASEDENS_ARRAY:
    case OUT_SOFT_ARRAY:
    case OUT_GROUP_ARRAY:
    case OUT_RELAX_ARRAY:
	msrOutArray(ppy_msr,fname,iType);
	break;

    case OUT_POS_VECTOR:
    case OUT_VEL_VECTOR:
    case OUT_ACCEL_VECTOR:
    case OUT_MEANVEL_VECTOR:
	msrOutVector(ppy_msr,fname,iType);
	break;

    case OUT_GROUP_TIPSY_NAT:
    case OUT_GROUP_TIPSY_STD:
    case OUT_GROUP_STATS:
    case OUT_GROUP_PROFILES:
	msrOutGroups(ppy_msr,fname,iType,dTime);
	break;

    default:
	assert(0);
	}

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_SaveVector(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Type",NULL};
    const char *fname;
    int iType;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "si:SaveVector", kwlist,
	     &fname,&iType ) )
	return NULL;
    msrOutVector(ppy_msr,fname,iType);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_SaveArray(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Type",NULL};
    const char *fname;
    int iType;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "si:SaveArray", kwlist,
	     &fname,&iType ) )
	return NULL;
    msrOutArray(ppy_msr,fname,iType);
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
    {"SelSrcById", (PyCFunction)ppy_msr_SelSrcById, METH_VARARGS,
     "Selects source particles with a specific id range."},
    {"SelDstById", (PyCFunction)ppy_msr_SelDstById, METH_VARARGS,
     "Selects destination particles with a specific id range."},
    {"SelSrcPhaseDensity", (PyCFunction)ppy_msr_SelSrcPhaseDensity, METH_VARARGS|METH_KEYWORDS,
     "Selects source particles with a specific phase space density."},
    {"SelDstPhaseDensity", (PyCFunction)ppy_msr_SelDstPhaseDensity, METH_VARARGS|METH_KEYWORDS,
     "Selects destination particles with a specific phase space density."},
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
    {"Gravity", (PyCFunction)ppy_msr_Gravity, METH_VARARGS|METH_KEYWORDS,
     "Calculate gravity"},
    {"Smooth", (PyCFunction)ppy_msr_Smooth, METH_VARARGS|METH_KEYWORDS,
     "Smooth"},
    {"ReSmooth", (PyCFunction)ppy_msr_ReSmooth, METH_VARARGS|METH_KEYWORDS,
     "ReSmooth"},
    {"Fof", (PyCFunction)ppy_msr_Fof, METH_VARARGS|METH_KEYWORDS,
     "Friends of Friends"},
    {"GroupProfiles", (PyCFunction)ppy_msr_GroupProfiles, METH_VARARGS|METH_KEYWORDS,
     "Group Profiles"},
    {"PeakVc", (PyCFunction)ppy_msr_PeakVc, METH_VARARGS|METH_KEYWORDS,
     "Calculate peak circular velocities"},
    {"Load", (PyCFunction)ppy_msr_Load, METH_VARARGS|METH_KEYWORDS,
     "Load an input file"},
    {"Save", (PyCFunction)ppy_msr_Save, METH_VARARGS|METH_KEYWORDS,
     "Save particles to a file"},
    {"SaveVector", (PyCFunction)ppy_msr_SaveVector, METH_VARARGS|METH_KEYWORDS,
     "Save a vector to a file"},
    {"SaveArray", (PyCFunction)ppy_msr_SaveArray, METH_VARARGS|METH_KEYWORDS,
     "Save an array to a file"},

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

    PyDict_SetItemString(dict, "SMX_DENSITY", Py_BuildValue("i",SMX_DENSITY));
    PyDict_SetItemString(dict, "SMX_MEANVEL", Py_BuildValue("i",SMX_MEANVEL));
    PyDict_SetItemString(dict, "SMX_DIVV", Py_BuildValue("i",SMX_DIVV));
    PyDict_SetItemString(dict, "SMX_VELDISP2", Py_BuildValue("i",SMX_VELDISP2));
    PyDict_SetItemString(dict, "SMX_FOF", Py_BuildValue("i",SMX_FOF));
    PyDict_SetItemString(dict, "SMX_RELAXATION", Py_BuildValue("i",SMX_RELAXATION));

    PyDict_SetItemString(dict, "OUT_TIPSY_STD", Py_BuildValue("i",OUT_TIPSY_STD));
    PyDict_SetItemString(dict, "OUT_TIPSY_DBL", Py_BuildValue("i",OUT_TIPSY_DBL));
    PyDict_SetItemString(dict, "OUT_POS_VECTOR", Py_BuildValue("i",OUT_POS_VECTOR));
    PyDict_SetItemString(dict, "OUT_VEL_VECTOR", Py_BuildValue("i",OUT_VEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_ACCEL_VECTOR", Py_BuildValue("i",OUT_ACCEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_MEANVEL_VECTOR", Py_BuildValue("i",OUT_MEANVEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_IORDER_ARRAY", Py_BuildValue("i",OUT_IORDER_ARRAY));
    PyDict_SetItemString(dict, "OUT_COLOR_ARRAY", Py_BuildValue("i",OUT_COLOR_ARRAY));
    PyDict_SetItemString(dict, "OUT_DENSITY_ARRAY", Py_BuildValue("i",OUT_DENSITY_ARRAY));
    PyDict_SetItemString(dict, "OUT_POT_ARRAY", Py_BuildValue("i",OUT_POT_ARRAY));
    PyDict_SetItemString(dict, "OUT_AMAG_ARRAY", Py_BuildValue("i",OUT_AMAG_ARRAY));
    PyDict_SetItemString(dict, "OUT_IMASS_ARRAY", Py_BuildValue("i",OUT_IMASS_ARRAY));
    PyDict_SetItemString(dict, "OUT_RUNG_ARRAY", Py_BuildValue("i",OUT_RUNG_ARRAY));
    PyDict_SetItemString(dict, "OUT_SOFT_ARRAY", Py_BuildValue("i",OUT_SOFT_ARRAY));
    PyDict_SetItemString(dict, "OUT_DIVV_ARRAY", Py_BuildValue("i",OUT_DIVV_ARRAY));
    PyDict_SetItemString(dict, "OUT_VELDISP2_ARRAY", Py_BuildValue("i",OUT_VELDISP2_ARRAY));
    PyDict_SetItemString(dict, "OUT_VELDISP_ARRAY", Py_BuildValue("i",OUT_VELDISP_ARRAY));
    PyDict_SetItemString(dict, "OUT_PHASEDENS_ARRAY", Py_BuildValue("i",OUT_PHASEDENS_ARRAY));
    PyDict_SetItemString(dict, "OUT_SOFT_ARRAY", Py_BuildValue("i",OUT_SOFT_ARRAY));
    PyDict_SetItemString(dict, "OUT_GROUP_ARRAY", Py_BuildValue("i",OUT_GROUP_ARRAY));
    PyDict_SetItemString(dict, "OUT_RELAX_ARRAY", Py_BuildValue("i",OUT_RELAX_ARRAY));
    PyDict_SetItemString(dict, "OUT_GROUP_TIPSY_NAT", Py_BuildValue("i",OUT_GROUP_TIPSY_NAT));
    PyDict_SetItemString(dict, "OUT_GROUP_TIPSY_STD", Py_BuildValue("i",OUT_GROUP_TIPSY_STD));
    PyDict_SetItemString(dict, "OUT_GROUP_STATS", Py_BuildValue("i",OUT_GROUP_STATS));
    PyDict_SetItemString(dict, "OUT_GROUP_PROFILES", Py_BuildValue("i",OUT_GROUP_PROFILES));
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
    PyObject *dict, *globals;
    struct _node *node;

    assert(Py_IsInitialized());

    // Set parameters for easy access
    dict = PyModule_GetDict(ppy->module);

    PyDict_SetItemString(dict, "achOutName", Py_BuildValue("s",ppy_msr->param.achOutName));



    globals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__",
			 PyEval_GetBuiltins());



    printf("---------------------------------------"
	   "---------------------------------------\n"
	   "Running Python Script %s\n"
	   "---------------------------------------"
	   "---------------------------------------\n",
	   achFilename );
    fp = fopen(achFilename,"r");
#if 1
    PyRun_SimpleFile(fp,achFilename);
    fclose(fp);
#else
    node = PyParser_SimpleParseFile(fp,achFilename,Py_file_input);
    fclose(fp);
    if ( node ) {
	PyCodeObject *code = PyNode_Compile(node,achFilename);
	if ( code ) {
//	    PyObject *pFunc = PyObject_GetAttrString(code, "Nada");
//	    if ( pFunc ) {
//		printf("Yep\n");
//		}
	    PyObject *result = PyEval_EvalCode(code,globals,dict);
	    Py_DECREF(result);
	    }
	Py_DECREF(code);
	}
    if(PyErr_Occurred()) PyErr_Print();

#endif
    printf("---------------------------------------"
	   "---------------------------------------\n" );

    }
