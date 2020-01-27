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
#include "cosmo.h"
#include "csmpython.h"

#ifndef Py_TYPE
    #define Py_TYPE(ob) (((PyObject*)(ob))->ob_type)
#endif
#ifndef PyVarObject_HEAD_INIT
    #define PyVarObject_HEAD_INIT(type, size) \
        PyObject_HEAD_INIT(type) size,
#endif


/**********************************************************************\
 ** Converter functions
\**********************************************************************/

static int double_or_none(PyObject* o, void* d) {
    double value = *(double *)d;
    if (o==Py_None) {
    	*(double *)d = NAN;
    	return 1;
	}
    else if (PyFloat_Check(o)) value = PyFloat_AsDouble(o);
    else if (PyLong_Check(o)) value = PyLong_AsDouble(o);
    else {
	PyErr_SetString(PyExc_TypeError, "must be a number or None");
	return 0; // conversion failed
	}
    if (isnan(value)) PyErr_SetString(PyExc_TypeError, "cannot be nan");
    *(double *)d = value;
    return !PyErr_Occurred();
    }

/**********************************************************************\
 ** CSM Object methods
\**********************************************************************/

static PyObject *
ppy_csmClassRead(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"","Lbox","LinSpecies","kpivot","alphas","As","ns",NULL};
    double Lbox;
    double kpivot = 0.05;
    double alphas = 0.0;
    int i, nSpecies = 0;
    const char **aSpecies = NULL;
    const char *fname;
    PyObject *species = Py_None;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "sd|Odddd:ClassRead", kwlist,
	     &fname, &Lbox, &species, &kpivot, &alphas,
	     &self->csm->val.dNormalization,
	     &self->csm->val.dSpectral ) )
	return NULL;
    if (species!=Py_None) {
	if (!PyList_Check(species)) return PyErr_Format(PyExc_TypeError,"ClassRead() argument LinSpecies must be list or None");
	nSpecies = PyList_Size(species);
	aSpecies = malloc(nSpecies * sizeof(*aSpecies));    	assert(aSpecies!=NULL);
	for( i=0; i<nSpecies; ++i) {
            PyObject *o = PyList_GetItem(species,i);
#if PY_MAJOR_VERSION >= 3
            if (PyUnicode_Check(o)) aSpecies[i] = PyUnicode_AsUTF8(o);
#else
            if (PyString_Check(o)) aSpecies[i] = PyString_AsString(o);
#endif
	    else return PyErr_Format(PyExc_TypeError,"ClassRead() argument LinSpecies list must contain string elements");
            }
	}
    self->csm->val.classData.bClass = 1;
    self->csm->val.dPivot = kpivot;
    self->csm->val.dRunning = alphas;
    csmClassRead(self->csm,fname,Lbox,0.0,nSpecies,aSpecies,0,NULL);
    if (aSpecies) free(aSpecies);
    Py_INCREF(Py_None);
    return Py_None;
    }

static PyObject *
ppy_csmSetCosmology(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={ "dHubble0", "dOmega0", "dLambda", "dOmegaRad", "dOmegab", "dOmegaDE",
			    "w0", "wa", "dSigma8", "As", "ns", "running", "pivot", NULL};
    if ( !PyArg_ParseTupleAndKeywords(
	    args, kwobj, "|O&O&O&O&O&O&O&O&O&O&O&O&O&:SetCosmology", kwlist,
	    double_or_none, &self->csm->val.dHubble0,
	    double_or_none, &self->csm->val.dOmega0,
	    double_or_none, &self->csm->val.dLambda,
	    double_or_none, &self->csm->val.dOmegaRad,
	    double_or_none, &self->csm->val.dOmegab,
	    double_or_none, &self->csm->val.dOmegaDE,
	    double_or_none, &self->csm->val.w0,
	    double_or_none, &self->csm->val.wa,
	    double_or_none, &self->csm->val.dSigma8,
	    double_or_none, &self->csm->val.dNormalization,
	    double_or_none, &self->csm->val.dSpectral,
	    double_or_none, &self->csm->val.dRunning,
	    double_or_none, &self->csm->val.dPivot ) )
	return NULL;
    if (isnan(self->csm->val.dLambda))
    	self->csm->val.dLambda = 1.0 - self->csm->val.dOmega0
				     - self->csm->val.dOmegaRad
				     - self->csm->val.dOmegaDE;


    Py_INCREF(Py_None);
    return Py_None;
    }

static PyObject *
ppy_csmRhoBar_m(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a",NULL};
    double a;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:RhoBar_m", kwlist, &a) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmRhoBar_m(self->csm,a) );
    }

static PyObject *
ppy_csmRhoBar_lin(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a",NULL};
    double a;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:RhoBar_lin", kwlist, &a) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmRhoBar_lin(self->csm,a) );
    }

static PyObject *
ppy_csmDelta_m(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a","k",NULL};
    double a,k;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "dd:Delta_m", kwlist, &a, &k) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmDelta_m(self->csm,a,k));
    }

static PyObject *
ppy_csmTheta_m(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a","k",NULL};
    double a,k;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "dd:Theta_m", kwlist, &a, &k) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmTheta_m(self->csm,a,k) );
    }

static PyObject *
ppy_csmDelta_lin(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a","k",NULL};
    double a,k;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "dd:Delta_lin", kwlist, &a, &k) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmDelta_lin(self->csm,a,k) );
    }

static PyObject *
ppy_csmDeltaRho_lin(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a0","a1","k",NULL};
    double a0,a1,k;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "ddd:DeltaRho_lin", kwlist, &a0, &a1, &k) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmDeltaRho_lin(self->csm,a0,a1,k) );
    }

static PyObject *
ppy_csmZeta(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"k",NULL};
    double k;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:Zeta", kwlist, &k) ) return NULL;
    if (!self->csm->val.classData.bClass) return PyErr_Format(PyExc_AssertionError,"Class data not initialized (use ClassRead)");
    csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmZeta(self->csm,k) );
    }

static PyObject *
ppy_csmExp2Time(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a",NULL};
    double a;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:Exp2Time", kwlist, &a) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmExp2Time(self->csm,a) );
    }

static PyObject *
ppy_csmTime2Exp(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"t",NULL};
    double t;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:Time2Exp", kwlist, &t) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmTime2Exp(self->csm,t) );
    }

static PyObject *
ppy_csmTime2Hub(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"t",NULL};
    double t;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:Time2Hub", kwlist, &t) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmTime2Hub(self->csm,t) );
    }

static PyObject *
ppy_csmExp2Hub(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a",NULL};
    double a;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:Exp2Hub", kwlist, &a) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmExp2Hub(self->csm,a) );
    }

static PyObject *
ppy_csmComoveDriftFac(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"time","delta",NULL};
    double time,delta;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "dd:ComoveDriftFac", kwlist, &time, &delta) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmComoveDriftFac(self->csm,time,delta) );
    }

static PyObject *
ppy_csmComoveKickFac(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"time","delta",NULL};
    double time,delta;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "dd:ComoveKickFac", kwlist, &time, &delta) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    return Py_BuildValue("d", csmComoveKickFac(self->csm,time,delta) );
    }

static PyObject *
ppy_csmComoveGrowth(CSMINSTANCE *self, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"a",NULL};
    double a, D1LPT, D2LPT, f1LPT, f2LPT;
    if ( !PyArg_ParseTupleAndKeywords(args, kwobj, "d:ComoveGrowth", kwlist, &a) ) return NULL;
    if (self->csm->val.classData.bClass) csmClassGslInitialize(self->csm);
    csmComoveGrowth(self->csm, a, &D1LPT, &D2LPT, &f1LPT, &f2LPT);
    return Py_BuildValue("(dddd)", D1LPT, D2LPT, f1LPT, f2LPT );
    }

/**********************************************************************\
 ** CSM Object definition
\**********************************************************************/

static PyMethodDef csm_methods[] = {
    {"ClassRead", (PyCFunction)ppy_csmClassRead, METH_VARARGS|METH_KEYWORDS,
     "Read Cosmology from Class HDF5 file"},
    {"SetCosmology", (PyCFunction)ppy_csmSetCosmology, METH_VARARGS|METH_KEYWORDS,
     "Set Cosmology directly with parameters"},
    {"RhoBar_m", (PyCFunction)ppy_csmRhoBar_m, METH_VARARGS|METH_KEYWORDS,
     "Return RhoBar_m"},
    {"RhoBar_lin", (PyCFunction)ppy_csmRhoBar_lin, METH_VARARGS|METH_KEYWORDS,
     "Return RhoBar_lin"},
    {"Delta_m", (PyCFunction)ppy_csmDelta_m, METH_VARARGS|METH_KEYWORDS,
     "Return Delta_m"},
    {"Theta_m", (PyCFunction)ppy_csmTheta_m, METH_VARARGS|METH_KEYWORDS,
     "Return Theta_m"},
    {"Delta_lin", (PyCFunction)ppy_csmDelta_lin, METH_VARARGS|METH_KEYWORDS,
     "Return Delta_lin"},
    {"DeltaRho_lin", (PyCFunction)ppy_csmDeltaRho_lin, METH_VARARGS|METH_KEYWORDS,
     "Return DeltaRho_lin"},
    {"Zeta", (PyCFunction)ppy_csmZeta, METH_VARARGS|METH_KEYWORDS,
     "Return Zeta"},
    {"Exp2Time", (PyCFunction)ppy_csmExp2Time, METH_VARARGS|METH_KEYWORDS,
     "Return simulation time given expansion factor"},
    {"Time2Exp", (PyCFunction)ppy_csmTime2Exp, METH_VARARGS|METH_KEYWORDS,
     "Return expansion factor given simulation time"},
    {"Exp2Hub", (PyCFunction)ppy_csmExp2Hub, METH_VARARGS|METH_KEYWORDS,
     "Return Hubble parameter given expansion factor"},
    {"Time2Hub", (PyCFunction)ppy_csmTime2Hub, METH_VARARGS|METH_KEYWORDS,
     "Return Hubble parameter given simulation time"},
    {"ComoveDriftFac", (PyCFunction)ppy_csmComoveDriftFac, METH_VARARGS|METH_KEYWORDS,
     "Returns the drift factor"},
    {"ComoveKickFac", (PyCFunction)ppy_csmComoveKickFac, METH_VARARGS|METH_KEYWORDS,
     "Returns the kick factor"},
    {"ComoveGrowth", (PyCFunction)ppy_csmComoveGrowth, METH_VARARGS|METH_KEYWORDS,
     "Returns the growth factors"},
    {NULL, NULL, 0, NULL}
    };

static void csm_dealloc(CSMINSTANCE *self) {
    csmFinish(self->csm);
    Py_TYPE(self)->tp_free((PyObject*)self);
    }

static PyObject *csm_new(PyTypeObject *type, PyObject *args, PyObject *kwobj) {
    static char *kwlist[]={"As","ns","comoving",NULL};
    double As=0.0, ns=0.0;
    int bComove = 1;
    if ( !PyArg_ParseTupleAndKeywords(
	     args, kwobj, "|ddi:CSM", kwlist,
	     &As, &ns, &bComove ) )
	return NULL;
    CSMINSTANCE *self;
    self = (CSMINSTANCE *)type->tp_alloc(type, 0);
    if (self == NULL) { return NULL; }
    csmInitialize(&self->csm);
    self->csm->val.classData.bClass = 0;
    self->csm->val.dNormalization = As;
    self->csm->val.dSpectral = ns;
    self->csm->val.bComove = bComove;
    self->csm->val.dHubble0 = sqrt(8.0 * M_PI / 3.0);
    self->csm->val.dPivot = 0.05;
    return (PyObject *)self;
    }

static int csm_init(CSMINSTANCE *self, PyObject *args, PyObject *kwds) {
    return 0;
    }

static PyMemberDef csm_members[] = {
	{NULL} /* Sentinel */
    };

static int csm_setboolean(CSMINSTANCE *self, PyObject *v, void *offset) {
    int *pv = (int *)((char *)&self->csm->val + (long)(offset));
    *pv = PyObject_IsTrue(v);
    return 0;
    }

static PyObject *csm_getboolean(CSMINSTANCE *self, void *offset) {
    int *pv = (int *)((char *)&self->csm->val + (long)(offset));
    if (*pv) Py_RETURN_TRUE;
    else Py_RETURN_FALSE;
    }

static int csm_setdouble(CSMINSTANCE *self, PyObject *v, void *offset) {
    double *pv = (double *)((char *)&self->csm->val + (long)(offset));
    *pv = PyFloat_AsDouble(v);
    return 0;
    }

static PyObject *csm_getdouble(CSMINSTANCE *self, void *offset) {
    double *pv = (double *)((char *)&self->csm->val + (long)(offset));
    return Py_BuildValue("d", *pv);
    }

#define MAKE_BOOLEAN(a) {#a, (getter)csm_getboolean, (setter)csm_setboolean, #a, (void *)(offsetof(struct csmVariables,a))}
#define MAKE_DOUBLE(a) {#a, (getter)csm_getdouble, (setter)csm_setdouble, #a, (void *)(offsetof(struct csmVariables,a))}

static PyGetSetDef csm_getseters[] = {
    MAKE_BOOLEAN(bComove),
    MAKE_DOUBLE(dHubble0),
    MAKE_DOUBLE(dOmega0),
    MAKE_DOUBLE(dLambda),
    MAKE_DOUBLE(dOmegaRad),
    MAKE_DOUBLE(dOmegab),
    MAKE_DOUBLE(dOmegaDE),
    MAKE_DOUBLE(w0),
    MAKE_DOUBLE(wa),
    MAKE_DOUBLE(dSigma8),
    MAKE_DOUBLE(dNormalization),
    MAKE_DOUBLE(dSpectral),
    MAKE_DOUBLE(dRunning),
    MAKE_DOUBLE(dPivot),
    {NULL} /* Sentinel */
    };

PyTypeObject csmType = {
    PyVarObject_HEAD_INIT(NULL,0)
    "CSM", /*tp_name*/
    sizeof(CSMINSTANCE), /*tp_basicsize */
    0, /*tp_itemsize*/
    (destructor)csm_dealloc, /*tp_dealloc*/
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
    "CSM objects", /* tp_doc */
    0, /* tp_traverse */
    0, /* tp_clear */
    0, /* tp_richcompare */
    0, /* tp_weaklistoffset */
    0, /* tp_iter */
    0, /* tp_iternext */
    csm_methods, /* tp_methods */
    csm_members, /* tp_members */
    csm_getseters, /* tp_getset */
    0, /* tp_base */
    0, /* tp_dict */
    0, /* tp_descr_get */
    0, /* tp_descr_set */
    0, /* tp_dictoffset */
    (initproc)csm_init, /* tp_init */
    0, /* tp_alloc */
    csm_new, /* tp_new */
    };
/**********************************************************************\
\**********************************************************************/

struct csm_module_state {
    };

static PyMethodDef csm_module_methods[] = {
    {NULL, NULL}
    };

/**********************************************************************\
 ** Module definition. Note that the class will be defined after import.
\**********************************************************************/
#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef csm_moduledef = {
        PyModuleDef_HEAD_INIT,
        "CSM",
        NULL,
        sizeof(struct csm_module_state),
        csm_module_methods,
        NULL,
        NULL,
        NULL,
        NULL
    };
#endif

/**********************************************************************\
 ** Parallel Python (ppy) setup
\**********************************************************************/
#if PY_MAJOR_VERSION >= 3
PyObject * PyInit_CSM(void)
#else
void initCSM(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&csm_moduledef);
#else
    PyObject *module = Py_InitModule("CSM", csm_module_methods);
#endif
    /* Initialize "CSM" object as well. */
    if (PyType_Ready(&csmType) >= 0) {
	Py_INCREF(&csmType);
	PyModule_AddObject(module, "CSM", (PyObject *)&csmType);
	}
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
    }
