#include <Python.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include "master.h"

/**********************************************************************\
 ** GLOBAL variables...  Sadly Python is not reentrant
\**********************************************************************/
static MSR ppy_msr;

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
ppy_msr_SelSrcMass(PyObject *self, PyObject *args) {
    int ok;
    double dMinMass, dMaxMass;
    uint64_t nSelected;

    ok = PyArg_ParseTuple(args, "dd", &dMinMass, &dMaxMass);
    if ( ok ) {
	nSelected = msrSelSrcMass(ppy_msr,dMinMass,dMaxMass);
	}
    else {
	printf( "ERROR: Called SelSrcMass with invalid arguments\n");
	nSelected = 0;
	}
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_SelDstMass(PyObject *self, PyObject *args) {
    int ok;
    double dMinMass, dMaxMass;
    uint64_t nSelected;

    ok = PyArg_ParseTuple(args, "dd", &dMinMass, &dMaxMass);
    if ( ok ) {
	nSelected = msrSelDstMass(ppy_msr,dMinMass,dMaxMass);
	}
    else {
	printf( "ERROR: Called SelDstMass with invalid arguments\n");
	nSelected = 0;
	}
    return Py_BuildValue("L", nSelected);
}

static PyObject *
ppy_msr_DeepestPotential(PyObject *self, PyObject *args) {
    double r[3];
    float fPot;
    msrDeepestPot(ppy_msr,r,&fPot);
    return Py_BuildValue("((ddd)f)", r[0], r[1], r[2], fPot);
}

static PyObject *
ppy_msr_Profile(PyObject *self, PyObject *args) {
    int ok;
    double r[3], dMinR, dMaxR, dMassEncl, dRho;
    int nBins, i;
    const PROFILEBIN *pBins;
    PyObject *List = NULL;

    nBins = 100;
    ok = PyArg_ParseTuple(args, "(ddd)dd|i:Profile",
			  r+0, r+1, r+2,
			  &dMinR, &dMaxR, &nBins );
    if ( ok ) {
	pBins = msrProfile(ppy_msr,r,dMinR,dMaxR,nBins);
	List = PyList_New( nBins );
	dMassEncl = 0.0;
	for(i=0; i<nBins; i++ ) {
	    dMassEncl += pBins[i].dMassInBin;
	    dRho = pBins[i].dMassInBin/pBins[i].dVolume;
	    PyList_SetItem(List,i,Py_BuildValue("(dLdd)",pBins[i].dRadius, pBins[i].nParticles, dRho, dMassEncl));
	    }
	msrDeleteProfile(ppy_msr);
	}
    return List;
}

static PyMethodDef ppy_msr_methods[] = {
    {"SelSrcAll", ppy_msr_SelSrcAll, METH_NOARGS,
     "Selects all particles as operation source."},
    {"SelDstAll", ppy_msr_SelDstAll, METH_NOARGS,
     "Selects all particles as operation destination."},
    {"SelSrcMass", ppy_msr_SelSrcMass, METH_VARARGS,
     "Selects source particles with a specific mass range."},
    {"SelDstMass", ppy_msr_SelDstMass, METH_VARARGS,
     "Selects destination particles with a specific mass range."},
    {"DeepestPotential", ppy_msr_DeepestPotential, METH_NOARGS,
     "Finds the most bound particles (deepest potential)"},
    {"Profile", ppy_msr_Profile, METH_VARARGS,
     "Generate a density profile"},
    {NULL, NULL, 0, NULL}
};

/**********************************************************************\
 ** Parallel Python (ppy) setup
\**********************************************************************/

void ppyInitialize(MSR msr) {
    ppy_msr = msr;
    Py_Initialize();
    Py_InitModule( "msr", ppy_msr_methods );
    }

void ppyFinish() {
    Py_Finalize();
    }

void ppyRunScript(const char *achFilename) {
    FILE *fp;
    assert(Py_IsInitialized());
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
