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

#include <Python.h>
#include <structmember.h>
#include <marshal.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#ifndef Py_TYPE
    #define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif
#ifndef PyVarObject_HEAD_INIT
    #define PyVarObject_HEAD_INIT(type, size) \
        PyObject_HEAD_INIT(type) size,
#endif

#include <stdint.h>
#include "master.h"
#include "outtype.h"
#include "intype.h"
#include "pkdpython.h"

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

static void addToDict(PyObject *dict,PRM prm,PRM_NODE *pn) {
    if (!prmSpecified(prm,pn->pszName))
    {
	PyDict_SetItemString(dict, pn->pszName, Py_None);
	return;
    }

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
	if (pn->bArg) addToDict(global,ppy_msr->prm,pn);
	addToDict(local,ppy_msr->prm,pn);
	}

    PyDict_SetItemString(local, "nThreads", Py_BuildValue("i",ppy_msr->nThreads));
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
	    if (v == Py_None) continue;
	    pn->bArg = 1;
	    switch(pn->iType) {
	    case 0:
	    case 1:
		assert(pn->iSize == sizeof(int));
#if PY_MAJOR_VERSION >= 3
		if (PyLong_Check(v)) *(int *)pn->pValue = PyLong_AsLong(v);
#else
		if (PyInt_Check(v)) *(int *)pn->pValue = PyInt_AsLong(v);
#endif		
		else if (PyFloat_Check(v)) *(int *)pn->pValue = (int)PyFloat_AsDouble(v);
		else fprintf(stderr,"Invalid type for %s\n",pn->pszName);
		break;
	    case 2:
		assert(pn->iSize == sizeof(double));
		if (PyFloat_Check(v)) *(double *)pn->pValue = PyFloat_AsDouble(v);
#if PY_MAJOR_VERSION >= 3
		else if (PyLong_Check(v)) *(double *)pn->pValue = PyLong_AsLong(v);
#else
		else if (PyInt_Check(v)) *(double *)pn->pValue = PyInt_AsLong(v);
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
		*(uint64_t *)pn->pValue = PyLong_AsLong(v);
#else
		*(uint64_t *)pn->pValue = PyInt_AsLong(v);
#endif
		break;
		}
	    }
	}
    }

static int ppy_get_dTime(double *dTime)
{
    PyObject *v, *dict;
    dict = PyModule_GetDict(global_ppy->module);
    if ( (v = PyDict_GetItemString(dict, "dTime")) == NULL )
	return 0;
    *dTime = PyFloat_AsDouble(v);
    return 1;
}


/**********************************************************************\
 ** MASTER Interface (msr functions)
\**********************************************************************/

static PyObject *
ppy_msr_SelAll(PyObject *self, PyObject *args) {
    ppy2prm();
    msrSelAll(ppy_msr);
    return Py_BuildValue("L", ppy_msr->N);
}

static PyObject *
ppy_msr_SelMass(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinMass", "MaxMass", "setIfTrue", "clearIfFalse", NULL };
    double dMinMass, dMaxMass;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelMass", kwlist,
	     &dMinMass, &dMaxMass, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelMass(ppy_msr,dMinMass,dMaxMass,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelPhaseDensity(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[] = { "MinDensity", "MaxDensity", "setIfTrue", "clearIfFalse", NULL };
    double dMinDensity, dMaxDensity;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "dd|ii:SelPhaseDensity", kwlist,
	     &dMinDensity, &dMaxDensity, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelPhaseDensity(ppy_msr,dMinDensity,dMaxDensity,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelBox(PyObject *self, PyObject *args) {
    double dCenter[3], dSize[3];
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)|ii:SelBox",
	     dCenter+0, dCenter+1, dCenter+2, dSize+0, dSize+1, dSize+2, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelBox(ppy_msr,dCenter,dSize,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelSphere(PyObject *self, PyObject *args) {
    double r[3], dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTuple(
	     args, "(ddd)d|ii:SelSphere",
	     r+0, r+1, r+2, &dRadius, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelSphere(ppy_msr,r,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelCylinder(PyObject *self, PyObject *args) {
    double dP1[3],dP2[3],dRadius;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTuple(
	     args, "(ddd)(ddd)d|ii:SelSphere",
	     dP1+0, dP1+1, dP1+2, dP2+0, dP2+1, dP2+2, &dRadius,
	     &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelCylinder(ppy_msr,dP1,dP2,dRadius,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelById(PyObject *self, PyObject *args) {
    uint64_t idStart,idEnd;
    int setIfTrue=1, clearIfFalse=1;
    uint64_t nSelected;

    ppy2prm();
    if ( !PyArg_ParseTuple(
	     args, "LL|ii:SelById",
	     &idStart, &idEnd, &setIfTrue, &clearIfFalse) )
	return NULL;
    nSelected = msrSelById(ppy_msr,idStart,idEnd,setIfTrue,clearIfFalse);
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_TotalMass(PyObject *self, PyObject *args) {
    double dMass;
    ppy2prm();
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

    ppy2prm();
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
    ppy2prm();
    msrReorder(ppy_msr);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_DomainDecomp(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Rung","SplitVA",NULL};
    int iRung    = 0;
    int bOthers  = 0;
    int bSplitVA = 0;

    ppy2prm();
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|iii:DomainDecomp", kwlist,
	     &iRung, &bOthers, &bSplitVA ) )
	return NULL;
    msrDomainDecomp(ppy_msr,iRung,bOthers,bSplitVA);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_UpdateRung(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Rung",NULL};
    int iRung    = 0;
    ppy2prm();
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|i:UpdateRung", kwlist,
	     &iRung ) )
	return NULL;
    msrUpdateRung(ppy_msr,iRung);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_BuildTree(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time","Ewald",NULL};
    double dTime = 0.0;
    int bEwald = ppy_msr->param.bEwald;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;

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
    uint64_t nActive;
    int iSec = 0;
    int iRung    = 0;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|idi:Gravity", kwlist,
	     &iRung, &dTime, &bEwald ) )
	return NULL;
    msrGravity(ppy_msr,iRung,MAX_RUNG,ROOT,0,dTime,ppy_msr->param.iStartStep,0,0,bEwald,ppy_msr->param.nGroup,&iSec,&nActive);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_Hop(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Time",NULL};
    double dExp;
    double dTime = 0.0;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|d:Hop", kwlist,
	     &dTime ) )
	return NULL;
    dExp = csmTime2Exp(ppy_msr->param.csm,dTime);
    msrHop(ppy_msr,dExp);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
ppy_msr_Smooth(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"iSmoothType","bSymmetric","dTime","nSmooth",NULL};
    int iSmoothType;
    int bSymmetric = 0;
    double dTime;
    int nSmooth = 64;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|idi:Smooth", kwlist,
	     &iSmoothType,&bSymmetric, &dTime, &nSmooth ) )
	return NULL;

    msrSmooth(ppy_msr,dTime,iSmoothType,bSymmetric,nSmooth);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_ReSmooth(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"iSmoothType","bSymmetric","dTime",NULL};
    int iSmoothType;
    int bSymmetric = 0;
    double dTime;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|id:ReSmooth", kwlist,
	     &iSmoothType,&bSymmetric, &dTime ) )
	return NULL;
    msrReSmooth(ppy_msr,dTime,iSmoothType,bSymmetric);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
ppy_msr_AdjustTime(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"aNew",NULL};
    int N, i;
    double aOld, aNew;
    double dTime;
    PyObject *dict;

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    aOld = csmTime2Exp(ppy_msr->param.csm,dTime);

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "d:AdjustTime", kwlist,
	     &aNew ) )
	return NULL;

    dTime = msrAdjustTime(ppy_msr,aOld,aNew);
    dict = PyModule_GetDict(global_ppy->module);
    PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",dTime));

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
ppy_msr_InitGrid(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"x","y","z",NULL};
    int x, y, z;


    ppy2prm();
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
    ppy2prm();
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
    float *fK,*fPk;
    PyObject *List, *value;

    ppy2prm();

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "i|dddd:MeasurePk", kwlist,
	     &nGrid, dCenter+0, dCenter+1, dCenter+2, &dRadius ) )
	return NULL;
    iNyquist = nGrid/2;


    fPk = malloc(sizeof(float)*(iNyquist+1));
    fK = malloc(sizeof(float)*(iNyquist+1));
    msrMeasurePk(ppy_msr,nGrid,nGrid/2,NULL,fK,fPk);

    List = PyList_New( iNyquist+1 );
    assert( List !=NULL );
    for( i=0; i<=iNyquist; i++ ) {
	value = Py_BuildValue("f",fPk[i]);
	assert( value != NULL );
	assert( PyList_SetItem(List,i,value) >= 0 );
	}
    free(fK);
    free(fPk);

    return List;
}
#endif

static PyObject *
ppy_msr_GenerateIC(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Grid","Seed",NULL};
    int nGrid, iSeed;

    ppy2prm();

    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|ii:GenerateIC", kwlist,
	     &ppy_msr->param.nGrid, &ppy_msr->param.iSeed ) )
	return NULL;

    msrGenerateIC(ppy_msr);    


    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
ppy_msr_Load(PyObject *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"Name","Type",NULL};
    const char *fname;
    double dTime;
    int iType = IN_TIPSY_STD;
    PyObject *dict;

    ppy2prm();
#if 0
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
    if ( (v = PyDict_GetItemString(dict, "bMemRelaxation")) != NULL )
	ppy_msr->param.bMemRelaxation = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemVelSmooth")) != NULL )
	ppy_msr->param.bMemVelSmooth = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemNodeVBnd")) != NULL )
	ppy_msr->param.bMemNodeVBnd = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "iDomainMethod")) != NULL )
	ppy_msr->param.iDomainMethod = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bParaRead")) != NULL )
	ppy_msr->param.bParaRead = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bParaWrite")) != NULL )
	ppy_msr->param.bParaWrite = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bMemNodeSphBounds")) != NULL )
	ppy_msr->param.bMemNodeSphBounds = PyInt_AsLong(v);
    if ( (v = PyDict_GetItemString(dict, "bStandard")) != NULL )
	ppy_msr->param.bStandard = PyInt_AsLong(v);
#endif

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
	dict = PyModule_GetDict(global_ppy->module);
	PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",dTime));
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

    ppy2prm();
    if (!ppy_get_dTime(&dTime))
	return NULL;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "s|id:Save", kwlist,
	     &fname,&iType, &dTime ) )
	return NULL;
    switch(iType) {
    case OUT_TIPSY_STD:
    case OUT_TIPSY_DBL:
	msrWrite(ppy_msr,fname,dTime,iType==OUT_TIPSY_DBL);
	break;

    case OUT_HOP_STATS:
	msrHopWrite(ppy_msr,fname);
	break;

    case OUT_IORDER_ARRAY:
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
    case OUT_MARKED_ARRAY:
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
    ppy2prm();
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
    ppy2prm();
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "si:SaveArray", kwlist,
	     &fname,&iType ) )
	return NULL;
    msrOutArray(ppy_msr,fname,iType);
    Py_INCREF(Py_None);
    return Py_None;
}

/**********************************************************************\
 * MSR methods.  These methods are shared by both the "msr" module,
 * and the "MSR" object.
\**********************************************************************/

static PyMethodDef msr_methods[] = {
    {"SelAll", ppy_msr_SelAll, METH_NOARGS,
     "Selects all particles as operation source."},
    {"SelMass", (PyCFunction)ppy_msr_SelMass, METH_VARARGS|METH_KEYWORDS,
     "Selects source particles with a specific mass range."},
    {"SelById", (PyCFunction)ppy_msr_SelById, METH_VARARGS,
     "Selects source particles with a specific id range."},
    {"SelPhaseDensity", (PyCFunction)ppy_msr_SelPhaseDensity, METH_VARARGS|METH_KEYWORDS,
     "Selects source particles with a specific phase space density."},
    {"SelBox", ppy_msr_SelBox, METH_VARARGS,
     "Selects source particles inside a given box."},
    {"SelSphere", ppy_msr_SelSphere, METH_VARARGS,
     "Selects source particles inside a given sphere."},
    {"SelCylinder", ppy_msr_SelCylinder, METH_VARARGS,
     "Selects source particles inside a given cylinder."},
    {"TotalMass", ppy_msr_TotalMass, METH_NOARGS,
     "Returns the total mass of the selected particles"},
    {"Profile", ppy_msr_Profile, METH_VARARGS,
     "Generate a density profile"},
    {"Reorder", ppy_msr_Reorder, METH_NOARGS,
     "Reorders the particles by iOrder"},
    {"DomainDecomp", (PyCFunction)ppy_msr_DomainDecomp, METH_VARARGS|METH_KEYWORDS,
     "Reorder the particles by position"},
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
    {"Hop", (PyCFunction)ppy_msr_Hop, METH_VARARGS|METH_KEYWORDS,
     "Grasshopper"},
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
    {"GenerateIC", (PyCFunction)ppy_msr_GenerateIC, METH_VARARGS|METH_KEYWORDS,
     "Generate Initial Condition"},
    {"Load", (PyCFunction)ppy_msr_Load, METH_VARARGS|METH_KEYWORDS,
     "Load an input file"},
    {"Save", (PyCFunction)ppy_msr_Save, METH_VARARGS|METH_KEYWORDS,
     "Save particles to a file"},
    {"SaveVector", (PyCFunction)ppy_msr_SaveVector, METH_VARARGS|METH_KEYWORDS,
     "Save a vector to a file"},
    {"SaveArray", (PyCFunction)ppy_msr_SaveArray, METH_VARARGS|METH_KEYWORDS,
     "Save an array to a file"},
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
    Py_TYPE(self)->tp_free((PyObject*)self);
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
/**********************************************************************\
\**********************************************************************/

static void setConstants( PyObject *dict ) {
    PyDict_SetItemString(dict, "dTime", Py_BuildValue("d",global_ppy->dTime));

    PyDict_SetItemString(dict, "SMX_DENSITY", Py_BuildValue("i",SMX_DENSITY));
    PyDict_SetItemString(dict, "SMX_DENSITY_M3", Py_BuildValue("i",SMX_DENSITY_M3));
    PyDict_SetItemString(dict, "SMX_GRADIENT_M3", Py_BuildValue("i",SMX_GRADIENT_M3));
    PyDict_SetItemString(dict, "SMX_MEANVEL", Py_BuildValue("i",SMX_MEANVEL));
    PyDict_SetItemString(dict, "SMX_DIVV", Py_BuildValue("i",SMX_DIVV));
    PyDict_SetItemString(dict, "SMX_VELDISP2", Py_BuildValue("i",SMX_VELDISP2));
    PyDict_SetItemString(dict, "SMX_FOF", Py_BuildValue("i",SMX_FOF));
    PyDict_SetItemString(dict, "SMX_RELAXATION", Py_BuildValue("i",SMX_RELAXATION));

    PyDict_SetItemString(dict, "OUT_TIPSY_STD", Py_BuildValue("i",OUT_TIPSY_STD));
    PyDict_SetItemString(dict, "OUT_TIPSY_DBL", Py_BuildValue("i",OUT_TIPSY_DBL));
    PyDict_SetItemString(dict, "OUT_HOP_STATS", Py_BuildValue("i",OUT_HOP_STATS));
    PyDict_SetItemString(dict, "OUT_POS_VECTOR", Py_BuildValue("i",OUT_POS_VECTOR));
    PyDict_SetItemString(dict, "OUT_VEL_VECTOR", Py_BuildValue("i",OUT_VEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_ACCEL_VECTOR", Py_BuildValue("i",OUT_ACCEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_MEANVEL_VECTOR", Py_BuildValue("i",OUT_MEANVEL_VECTOR));
    PyDict_SetItemString(dict, "OUT_IORDER_ARRAY", Py_BuildValue("i",OUT_IORDER_ARRAY));
    PyDict_SetItemString(dict, "OUT_RUNGDEST_ARRAY", Py_BuildValue("i",OUT_RUNGDEST_ARRAY));
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
    PyDict_SetItemString(dict, "OUT_MARKED_ARRAY", Py_BuildValue("i",OUT_MARKED_ARRAY));
    PyDict_SetItemString(dict, "OUT_PSGROUP_ARRAY", Py_BuildValue("i",OUT_PSGROUP_ARRAY));
    PyDict_SetItemString(dict, "OUT_RELAX_ARRAY", Py_BuildValue("i",OUT_RELAX_ARRAY));
    PyDict_SetItemString(dict, "OUT_PSGROUP_STATS", Py_BuildValue("i",OUT_PSGROUP_STATS));
    }


/**********************************************************************\
 ** Parallel Python (ppy) setup
\**********************************************************************/
static
#if PY_MAJOR_VERSION >= 3
PyObject *
#else
void
#endif
initModuleMSR(void) {
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
#if PY_MAJOR_VERSION >= 3
    return (PyObject *)&msrType;
#endif
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
	prm2ppy();
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
#if PY_MAJOR_VERSION > 2
    wchar_t **wargv;
    int i;
    wargv = malloc(sizeof(*wargv)*argc);
    for(i=0; i<argc; ++i) wargv[i] = Py_DecodeLocale(argv[i],NULL);
#endif
    assert(Py_IsInitialized());
    assert(argc>0);

    globals = PyDict_New();
    PyDict_SetItemString(globals, "__builtins__",
			 PyEval_GetBuiltins());

#if PY_MAJOR_VERSION > 2
    PySys_SetArgv(argc, wargv);
#else
    PySys_SetArgv(argc, argv);
#endif

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

#if PY_MAJOR_VERSION > 2
    for(i=0; i<argc; ++i) PyMem_RawFree(wargv[i]);
    free(wargv);
#endif
    }
