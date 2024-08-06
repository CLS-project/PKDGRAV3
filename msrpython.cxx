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
#include <numpy/arrayobject.h>
#include <structmember.h> // for PyMemberDef
#include "parse.h"
#include "master.h"
#include "csmpython.h"

#define CYTHON_EXTERN_C extern "C++"
#include "modules/checkpoint.h"
#include "modules/PKDGRAV.h"

/******************************************************************************\
*   Analysis callback
\******************************************************************************/

void MSR::addAnalysis(PyObject *callback,uint64_t per_particle, uint64_t per_process) {
    analysis_callbacks.emplace_back(callback,per_particle,per_process);
}

void MSR::runAnalysis(int iStep,double dTime) {
    set_dynamic(iStep,dTime);
    for ( msr_analysis_callback &i : analysis_callbacks) {
        auto v = parameters.call_dynamic(i.callback);
        Py_XDECREF(v);
        if (PyErr_Occurred()) PyErr_Print();
    }
}

/******************************************************************************\
*   Setup MSR using Python to parse parameters / enter analysis mode
\******************************************************************************/

extern "C" PyObject *PyInit_CSM(void);
extern "C" PyObject *PyInit_accuracy(void);
extern "C" PyObject *PyInit_cosmology(void);
extern "C" PyObject *PyInit_ephemeral(void);

int MSR::Python(int argc, char *argv[]) {
    PyImport_AppendInittab("PKDGRAV",PyInit_PKDGRAV);
    PyImport_AppendInittab("cosmology",PyInit_cosmology);
    PyImport_AppendInittab("CSM",PyInit_CSM);
    PyImport_AppendInittab("parse", PyInit_parse);
    PyImport_AppendInittab("checkpoint", PyInit_checkpoint);
    PyImport_AppendInittab("accuracy", PyInit_accuracy);
    PyImport_AppendInittab("ephemeral", PyInit_ephemeral);

    PKDGRAV_msr0 = this;
    PyStatus status;
    PyConfig config;
    PyConfig_InitPythonConfig(&config);
    PyConfig_Read(&config);
    config.site_import = 1;
    config.user_site_directory = 1;
    config.parse_argv = 0;
    status = PyConfig_SetBytesArgv(&config, argc, argv);
    // I don't like this, but it works for pyenv. See:
    //   https://github.com/python/cpython/issues/66409
    //   https://www.python.org/dev/peps/pep-0432/
    auto PYENV_VIRTUAL_ENV = getenv("VIRTUAL_ENV");
    if (!PYENV_VIRTUAL_ENV) PYENV_VIRTUAL_ENV = getenv("PYENV_VIRTUAL_ENV");
    if (PYENV_VIRTUAL_ENV) {
        std::string exec = PYENV_VIRTUAL_ENV;
        exec += "/bin/python";
        status = PyConfig_SetBytesString(&config,&config.program_name,exec.c_str());
    }
    status = Py_InitializeFromConfig(&config);
    PyConfig_Clear(&config);
    Initialize();

    PyObject *main_module = PyImport_ImportModule("__main__");
    auto globals = PyModule_GetDict(main_module);
    auto locals = globals;
    PyDict_SetItemString(globals, "__builtins__",PyEval_GetBuiltins());

    if (!PyImport_ImportModule("checkpoint")) {
        PyErr_Print();
        abort();
    }

    // This module is used for checkpointing and restoring the state
    pDill = PyImport_ImportModule("dill");
    if (!pDill) {
        PyErr_Print();
        abort();
    }
    pDill_load = PyObject_GetAttrString(pDill, "load");
    if (!pDill_load || !PyCallable_Check(pDill_load)) {
        PyErr_Print();
        abort();
    }
    pDill_dump = PyObject_GetAttrString(pDill, "dump");
    if (!pDill_dump || !PyCallable_Check(pDill_dump)) {
        PyErr_Print();
        abort();
    }
    // load_session/dump_session were "renamed" to load_module/dump_module in dill 0.3.6
    // the old name is still available in 0.3.8, but it is deprecated
    // we use the old name, and fallback to the new name if the old name is not available
    // presumably, the old name will be removed in a future version
    pDill_load_module = PyObject_GetAttrString(pDill, "load_session");
    if (!pDill_load_module || !PyCallable_Check(pDill_load_module)) {
        PyErr_Clear();
        pDill_load_module = PyObject_GetAttrString(pDill, "load_module");
    }
    if (!pDill_load_module || !PyCallable_Check(pDill_load_module)) {
        PyErr_Print();
        abort();
    }
    pDill_dump_module = PyObject_GetAttrString(pDill, "dump_session");
    if (!pDill_dump_module || !PyCallable_Check(pDill_dump_module)) {
        PyErr_Clear();
        pDill_dump_module = PyObject_GetAttrString(pDill, "dump_module");
    }
    if (!pDill_dump_module || !PyCallable_Check(pDill_dump_module)) {
        PyErr_Print();
        abort();
    }

    // Parse the command line
    auto PARSE = PyImport_ImportModule("parse");
    if (!PARSE) {
        PyErr_Print();
        abort();
    }
    auto result = parse();
    if (!result) {
        PyErr_Print();
        abort();
    }
    // Retrieve the results
    int n = PyTuple_Size(result);
    if (n!=2) {
        fprintf(stderr,"INTERNAL ERROR: parse.parse() MUST return a tuple\n");
        abort();
    }
    auto arguments = PyTuple_GetItem(result,0);         // Borrowed: Values of each parameter
    auto specified = PyTuple_GetItem(result,1);         // Borrowed: If it was explicitely specified
    parameters = pkd_parameters(arguments,specified);   // This will take ownership
    Py_DECREF(result);
    PyObject *script = PyObject_GetAttrString(arguments,"script");

    // If a script was specified then we run it.
    if (script != Py_None) {
        char *filename;
        if (PyUnicode_Check(script)) {
            PyObject *ascii = PyUnicode_AsASCIIString(script);
            filename = PyBytes_AsString(ascii);
            Py_DECREF(ascii);
        }
        else {
            fprintf(stderr,"INTERNAL ERROR: script filename is invalid\n");
            abort();
        }
        FILE *fp = fopen(filename,"r");
        if (fp == NULL) {
            perror(filename);
            exit(errno);
        }
        PyDict_SetItemString(locals, "__file__", script);
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

    auto imported = is_PKDGRAV_imported(locals);

    // If "MASTER" was imported then we are done -- the script should have done its job
    if (!imported) { // We must prepare for a normal legacy execution
        if (!parameters.verify(locals)) {
            PyErr_Print();
            fprintf(stderr,
                    "To avoid accidentially mistyping parameter names, you must prefix any additional\n"
                    "variables with an underscore. Verify the above listed variables/parameters.\n");
            exit(1);
        }
        update(locals,arguments,specified);
        if (PyErr_Occurred()) {
            PyErr_Print();
            exit(1);
        }
        parameters.set_dynamic("msr",PyImport_ImportModule("PKDGRAV"));
    }

    return imported ? 0 : -1;
}
