#include <Python.h>
#include <structmember.h>
#include <marshal.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdint.h>
#include "master.h"
#include "outtype.h"
#include "intype.h"
#include "python.h"

/**********************************************************************\
 ** GLOBAL variables...  Sadly Python is not reentrant
\**********************************************************************/
static MSR ppy_msr;

typedef struct {
    PyObject *module;
    PyObject *mainModule;
    double dTime;
    int bImported;
    } ppyCtx;

static ppyCtx *global_ppy = NULL;

/**********************************************************************\
 ** Interface to pkdgrav parameters
\**********************************************************************/

static void addToDict(PyObject *dict,PRM_NODE *pn) {
    switch (pn->iType) {
    case 0:
    case 1:
	assert(pn->iSize == sizeof(int));
	PyDict_SetItemString(dict, pn->pszName, Py_BuildValue("i",*(int *)pn->pValue));
	break;
    case 2:
	assert(pn->iSize == sizeof(double));
	PyDict_SetItemString(dict, pn->pszName, Py_BuildValue("d",*(double *)pn->pValue));
	break;
    case 3:
	PyDict_SetItemString(dict, pn->pszName, Py_BuildValue("s",pn->pValue));
	break;
    case 4:
	assert(pn->iSize == sizeof(uint64_t));
	PyDict_SetItemString(dict, pn->pszName, Py_BuildValue("L",*(uint64_t *)pn->pValue));
	break;
	}
    }

/* Copy parameters into python dictionary. */
static void prm2ppy(void) {
    PRM_NODE *pn;
    PyObject *local, *global;

    global = PyModule_GetDict(global_ppy->mainModule);
    local = PyModule_GetDict(global_ppy->module);

    /* We really shouldn't know about this structure, but what can you do? */
    for( pn=ppy_msr->prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	if (pn->bArg) addToDict(global,pn);
	addToDict(local,pn);
	}
    }

/* Copy parameters from python dictionary back into parameters. */
static void ppy2prm(void) {
    PyObject *global;
    PyObject *v;
    const char *s;
    PRM_NODE *pn;

    global = PyModule_GetDict(global_ppy->mainModule);

    for( pn=ppy_msr->prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	v = PyDict_GetItemString(global, pn->pszName);
	if (v!=NULL) {
	    pn->bArg = 1;
	    switch(pn->iType) {
	    case 0:
	    case 1:
		assert(pn->iSize == sizeof(int));
		*(int *)pn->pValue = PyInt_AsLong(v);
		break;
	    case 2:
		assert(pn->iSize == sizeof(double));
		*(double *)pn->pValue = PyFloat_AsDouble(v);
		break;
	    case 3:
		s = PyString_AsString(v);
		assert(pn->iSize > strlen(s));
		strcpy((char *)pn->pValue,s);
		break;
	    case 4:
		assert(pn->iSize == sizeof(uint64_t));
		*(uint64_t *)pn->pValue = PyInt_AsLong(v);
		break;
		}
	    }
	}
    }


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
	     args, "(ddd)(ddd)d|ii:SelDstSphere",
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
ppy_msr_TotalMass(PyObject *self, PyObject *args) {
    double dMass;
    dMass = msrTotalMass(ppy_msr);
    return Py_BuildValue("d", dMass );
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

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|ii:DomainDecomp", kwlist,
	     &iRung, &bSplitVA ) )
	return NULL;
    msrDomainDecomp(ppy_msr,iRung,bSplitVA);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_UpdateRung(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Rung",NULL};
    int iRung    = 0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|i:UpdateRung", kwlist,
	     &iRung ) )
	return NULL;
    msrUpdateRung(ppy_msr,iRung);
    Py_INCREF(Py_None);
    return Py_None;
}

#ifdef MPI_VERSION

static PyObject *
ppy_msr_DomainDecompNew(PyObject *self, PyObject *args, PyObject *kwobj) {
    msrDomainDecompNew(ppy_msr);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_RungOrder(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Rung",NULL};
    int iRung    = 0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|i:RungOrder", kwlist,
	     &iRung ) )
	return NULL;
    msrRungOrder(ppy_msr,iRung);
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

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
    static char *kwlist[]={"Rung","Time","Ewald",NULL};
    double dTime = 0.0;
    int bEwald = ppy_msr->param.bEwald;
    PyObject *v, *dict;
    uint64_t nActive;
    int iSec = 0;
    int iRung    = 0;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|idi:Gravity", kwlist,
	     &iRung, &dTime, &bEwald ) )
	return NULL;

    msrGravity(ppy_msr,iRung,MAX_RUNG,dTime,ppy_msr->param.iStartStep,bEwald,ppy_msr->param.nGroup,&iSec,&nActive);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Fof(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time",NULL};
    double dExp;
    double dTime = 0.0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|d:Fof", kwlist,
	     &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrFof(ppy_msr,dExp);
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
ppy_msr_AdjustTime(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"aNew",NULL};
    int N, i;
    double aOld, aNew;
    double dTime;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    aOld = csmTime2Exp(ppy_msr->param.csm,dTime);

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "d:AdjustTime", kwlist,
	     &aNew ) )
	return NULL;

    dTime = msrAdjustTime(ppy_msr,aOld,aNew);
    PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",dTime));

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_InitGrid(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"x","y","z",NULL};
    int x, y, z;


    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "iii:InitGrid", kwlist,
	     &x, &y, &z ) )
	return NULL;

    msrInitGrid(ppy_msr,x,y,z);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Project(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"x","y","z",NULL};
    double x=0.0, y=0.0, z=0.0;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|ddd:GridProject", kwlist,
	     &x, &y, &z ) )
	return NULL;
    msrGridProject(ppy_msr,x,y,z);
    Py_INCREF(Py_None);
    return Py_None;
}


#ifdef MDL_FFTW
static PyObject *
ppy_msr_MeasurePk(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"nGrid","x","y","z","r",NULL};
    double dCenter[3] = {0.0,0.0,0.0};
    double dRadius = 0.5;
    int nGrid, iNyquist, i;
    float *fPk;
    PyObject *List, *value;

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|dddd:MeasurePk", kwlist,
	     &nGrid, dCenter+0, dCenter+1, dCenter+2, &dRadius ) )
	return NULL;
    iNyquist = nGrid/2;


    fPk = malloc(sizeof(float)*(iNyquist+1));
    msrMeasurePk(ppy_msr,dCenter,dRadius,nGrid,fPk);

    List = PyList_New( iNyquist+1 );
    assert( List !=NULL );
    for( i=0; i<=iNyquist; i++ ) {
	value = Py_BuildValue("f",fPk[i]);
	assert( value != NULL );
	assert( PyList_SetItem(List,i,value) >= 0 );
	}
    return List;
}
#endif

static PyObject *
ppy_msr_GroupProfiles(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time",NULL};
    double dExp;
    double dTime = 0.0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|d:GroupProfiles", kwlist,
	     &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrGroupProfiles(ppy_msr,dExp);
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
    if ( (v = PyDict_GetItemString(dict, "bMemNodeVBnd")) != NULL )
	ppy_msr->param.bMemNodeVBnd = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "iDomainMethod")) != NULL )
	ppy_msr->param.iDomainMethod = PyInt_AsLong(v);

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
    case OUT_BALL_ARRAY:
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
    case OUT_PSGROUP_ARRAY:
    case OUT_RELAX_ARRAY:
    case OUT_RUNGDEST_ARRAY:
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

    case OUT_PSGROUP_STATS:
        msrOutPsGroups(ppy_msr,fname,iType,dTime);
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

#if 0
static PyObject *
ppy_msr_BuildPsdTree(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time","Ewald",NULL};
    double dTime = 0.0;
    int bEwald = ppy_msr->param.bEwald;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|di:BuildPsdTree", kwlist,
	     &dTime, &bEwald ) )
	return NULL;
    msrBuildPsdTree(ppy_msr,dTime,bEwald);
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

static PyObject *
ppy_msr_PSGroupFinder(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={NULL};
    PyObject *v, *dict;

    msrPSGroupFinder(ppy_msr);

    Py_INCREF(Py_None);
    return Py_None;
}

#if 0
static PyObject *
ppy_msr_PsFof(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time",NULL};
    double dExp;
    double dTime = 0.0;
    PyObject *v, *dict;

    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return NULL;
    dTime = PyFloat_AsDouble(v);
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|d:PsFof", kwlist,
	     &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrPsFof(ppy_msr,dExp);
    //msrGroupMerge(ppy_msr,dExp);
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

/**********************************************************************\
 * MSR methods.  These methods are shared by both the "msr" module,
 * and the "MSR" object.
\**********************************************************************/

static PyMethodDef msr_methods[] = {
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
    {"TotalMass", ppy_msr_TotalMass, METH_NOARGS,
     "Returns the total mass of the selected particles"},
    {"Profile", ppy_msr_Profile, METH_VARARGS,
     "Generate a density profile"},
    {"Reorder", ppy_msr_Reorder, METH_NOARGS,
     "Reorders the particles by iOrder"},
    {"DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
#ifdef MPI_VERSION
    {"DomainDecompNew", (PyCFunction)ppy_msr_DomainDecompNew, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
    {"RungOrder", (PyCFunction)ppy_msr_RungOrder, METH_VARARGS|METH_KEYWORDS,
     "Reorder particles to the domains for the given rung"},
#endif
    {"UpdateRung", (PyCFunction)ppy_msr_UpdateRung, METH_VARARGS|METH_KEYWORDS,
     "Update particle rungs"},
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
    {"AdjustTime", (PyCFunction)ppy_msr_AdjustTime, METH_VARARGS|METH_KEYWORDS,
     "Changing starting time for Zel'dovich ICs"},
    {"InitGrid", (PyCFunction)ppy_msr_InitGrid, METH_VARARGS|METH_KEYWORDS,
     "Initialize/allocate a GRID"},
    {"Project", (PyCFunction)ppy_msr_Project, METH_VARARGS|METH_KEYWORDS,
     "Project density onto the grid"},
#ifdef MDL_FFTW
    {"MeasurePk", (PyCFunction)ppy_msr_MeasurePk, METH_VARARGS|METH_KEYWORDS,
     "Measure the power spectrum"},
#endif
    {"Load", (PyCFunction)ppy_msr_Load, METH_VARARGS|METH_KEYWORDS,
     "Load an input file"},
    {"Save", (PyCFunction)ppy_msr_Save, METH_VARARGS|METH_KEYWORDS,
     "Save particles to a file"},
    {"SaveVector", (PyCFunction)ppy_msr_SaveVector, METH_VARARGS|METH_KEYWORDS,
     "Save a vector to a file"},
    {"SaveArray", (PyCFunction)ppy_msr_SaveArray, METH_VARARGS|METH_KEYWORDS,
     "Save an array to a file"},
    //{"BuildPsdTree", (PyCFunction)ppy_msr_BuildPsdTree, METH_VARARGS|METH_KEYWORDS,
     //"Build the phase-space tree"},
    {"PSGroupFinder", (PyCFunction)ppy_msr_PSGroupFinder, METH_NOARGS,
     "Calculate phase space density using EnBiD algorithm"},
#if 0
    {"PsFof", (PyCFunction)ppy_msr_PsFof, METH_VARARGS|METH_KEYWORDS,
     "Phase-space Friends of Friends"},
#endif

    {NULL, NULL, 0, NULL}
};

/**********************************************************************\
 * MSR object.  We use a Python object for the MSR so that it can be
 * accessed remotely via Pyro.  Here is what the pkdgrav2 python script
 * ("the server") would look like:
 *
 **********************************************************************
 * import msr
 * import Pyro.core
 *
 * class remoteMSR(msr.MSR,Pyro.core.ObjBase):
 *     def __init__(self):
 *         Pyro.core.ObjBase.__init__(self)
 *         msr.MSR.__init__(self)
 *
 * Pyro.core.initServer()
 * daemon=Pyro.core.Daemon()
 * uri=daemon.connect(remoteMSR(),"msr")
 *
 * print "The daemon runs on port:",daemon.port
 * print "The object's uri is:",uri
 *
 * daemon.requestLoop()
 **********************************************************************
 *
 * Here is how the client would connect to the server:
 *
 **********************************************************************
 * import Pyro.core
 * msr = Pyro.core.getProxyForURI("PYRO://ip:port/guid")
 * msr.Load("runc03.bin")
 **********************************************************************
 *
 * Obviously, the "Pyro" module needs to be available, and the URI
 * needs to be replaced by what the server prints out.
 *
\**********************************************************************/

typedef struct {
    PyObject_HEAD
    } MSRINSTANCE;

static void msr_dealloc(MSRINSTANCE *self) {
    self->ob_type->tp_free((PyObject*)self);
    }

static PyObject *msr_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
    MSRINSTANCE *self;
    self = (MSRINSTANCE *)type->tp_alloc(type, 0);
    if (self == NULL) { return NULL; }
    return (PyObject *)self;
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
    PyObject_HEAD_INIT(NULL)
    0, /*ob_size*/
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
/**********************************************************************\
\**********************************************************************/

static void setConstants( PyObject *dict ) {
    PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",global_ppy->dTime));

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
    PyDict_SetItemString(dict, "OUT_RUNGDEST_ARRAY", Py_BuildValue("i",OUT_RUNGDEST_ARRAY));
    PyDict_SetItemString(dict, "OUT_COLOR_ARRAY", Py_BuildValue("i",OUT_COLOR_ARRAY));
    PyDict_SetItemString(dict, "OUT_DENSITY_ARRAY", Py_BuildValue("i",OUT_DENSITY_ARRAY));
    PyDict_SetItemString(dict, "OUT_BALL_ARRAY", Py_BuildValue("i",OUT_BALL_ARRAY));
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
    PyDict_SetItemString(dict, "OUT_PSGROUP_ARRAY", Py_BuildValue("i",OUT_PSGROUP_ARRAY));
    PyDict_SetItemString(dict, "OUT_RELAX_ARRAY", Py_BuildValue("i",OUT_RELAX_ARRAY));
    PyDict_SetItemString(dict, "OUT_GROUP_TIPSY_NAT", Py_BuildValue("i",OUT_GROUP_TIPSY_NAT));
    PyDict_SetItemString(dict, "OUT_GROUP_TIPSY_STD", Py_BuildValue("i",OUT_GROUP_TIPSY_STD));
    PyDict_SetItemString(dict, "OUT_GROUP_STATS", Py_BuildValue("i",OUT_GROUP_STATS));
    PyDict_SetItemString(dict, "OUT_GROUP_PROFILES", Py_BuildValue("i",OUT_GROUP_PROFILES));
    PyDict_SetItemString(dict, "OUT_PSGROUP_STATS", Py_BuildValue("i",OUT_PSGROUP_STATS));
    }


/**********************************************************************\
 ** Parallel Python (ppy) setup
\**********************************************************************/
static void initModuleMSR(void) {
    PyObject *dict;

//    global_ppy->module = Py_InitModule("msr", msr_methods);
    global_ppy->module = global_ppy->mainModule;
    global_ppy->bImported = 1;

    prm2ppy();

    dict = PyModule_GetDict(global_ppy->module);
    setConstants(dict);

    /* Import Pyro.core. */
    /*
    PyObject *main = PyImport_ImportModule("__main__");
    PyObject *main_dict = PyModule_GetDict(main);
    PyObject *pyro = PyImport_ImportModule("Pyro");
    PyObject *pyro_dict = PyModule_GetDict(pyro);
    PyObject *core = PyImport_ImportModule("Pyro.core");
    PyObject *core_dict = PyModule_GetDict(core);
    PyObject *ObjBase = PyDict_GetItemString(core_dict,"ObjBase");
    PyDict_SetItemString(main_dict, "Pyro", pyro);
    PyDict_SetItemString(pyro_dict, "core", core);
    */
    /* Initialize "MSR" object as well. */
    if (PyType_Ready(&msrType) >= 0) {
	Py_INCREF(&msrType);
	PyModule_AddObject(global_ppy->module, "MSR", (PyObject *)&msrType);
	}
    }


void ppyInitialize(PPY *pvppy, MSR msr, double dTime) {
    ppyCtx *ppy;
    ppy = malloc(sizeof(ppyCtx));
    assert(ppy!=NULL);
    *pvppy = global_ppy = ppy;
    ppy_msr = msr;
    ppy->module = 0;
    ppy->dTime = dTime;
    ppy->bImported = 0;
    PyImport_AppendInittab("msr",initModuleMSR);
    Py_Initialize();
    ppy->mainModule = PyImport_AddModule("__main__"); 
    global_ppy->module = global_ppy->mainModule;
    if (PyType_Ready(&msrType) >= 0) {
	Py_INCREF(&msrType);
	PyObject *pymsr = PyObject_CallObject((PyObject *) &msrType, NULL);
	PyModule_AddObject(ppy->mainModule, "msr", pymsr);
	setConstants(PyModule_GetDict(global_ppy->mainModule));
	}
    }

void ppyFinish(PPY vppy) {
    ppyCtx *ppy = (ppyCtx *)vppy;
    Py_Finalize();
    free(ppy);
    global_ppy = NULL;
    }

void ppyRunScript(PPY vppy,int argc, char *argv[]) {
    ppyCtx *ppy = (ppyCtx *)vppy;
    FILE *fp;
    PyObject *dict, *globals;
    struct _node *node;

    assert(Py_IsInitialized());
    assert(argc>0);

    globals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__",
			 PyEval_GetBuiltins());

    PySys_SetArgv(argc, argv);

    printf("---------------------------------------"
	   "---------------------------------------\n"
	   "Running Python Script %s\n"
	   "---------------------------------------"
	   "---------------------------------------\n",
	   argv[0] );
    fp = fopen(argv[0],"r");
#if 1
    PyRun_SimpleFile(fp,argv[0]);
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
