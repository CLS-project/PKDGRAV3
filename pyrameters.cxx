/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2023 Douglas Potter & Joachim Stadel
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
#include "pyrameters.h"

bool pyrameters::has(const char *name) const {
    bool bSpecified = false;
    if (auto f = PyObject_GetAttrString(specified,name)) {
        bSpecified = PyObject_IsTrue(f)>0;
        Py_DECREF(f);
    }
    return bSpecified;
}

void pyrameters::merge(PyObject *o1, PyObject *o2) {
    auto d1 = PyObject_GenericGetDict(o1,nullptr);
    auto d2 = PyObject_GenericGetDict(o2,nullptr);
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(d2, &pos, &key, &value)) {
        PyDict_SetItem(d1,key,value);
    }
}

void pyrameters::merge(const pyrameters &other) {
    merge(this->arguments,other.arguments);
    merge(this->specified,other.specified);
}

bool pyrameters::verify(PyObject *kwobj) {
    bool bSuccess = true;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
    while (PyDict_Next(kwobj, &pos, &key, &value)) {
        const char *keyString;
        if (PyUnicode_Check(key)) {
            PyObject *ascii = PyUnicode_AsASCIIString(key);
            keyString = PyBytes_AsString(ascii);
            Py_DECREF(ascii);
            if (keyString[0]=='_') continue; // skip things that start with underscore
        }
        // If this key is not a valid argument, then print an error,
        // unless it is a module, or a callable imported from another module.
        if (!PyObject_HasAttr(arguments,key) && !PyModule_Check(value)) {
            if (PyCallable_Check(value)) {
                auto module = PyObject_GetAttrString(value,"__module__");
                if (module && PyUnicode_Check(module)) {
                    auto ascii = PyUnicode_AsASCIIString(module);
                    auto result = PyBytes_AsString(ascii);
                    Py_DECREF(ascii);
                    Py_DECREF(module);
                    if (result && strcmp(result,"__main__")!=0) continue;
                }
            }
            PyErr_Format(PyExc_AttributeError,"invalid parameter %A",key);
            PyErr_Print();
            bSuccess=false;
        }
    }
    return bSuccess;
}

bool pyrameters::update(PyObject *kwobj,bool bIgnoreUnknown) {
    bool bSuccess = true;
    PyObject *key, *value;
    Py_ssize_t pos = 0;
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
    return bSuccess;
}

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
        assert(pn->iSize == sizeof(std::uint64_t));
        ((std::uint64_t *)pn->pValue)[i] = PyLong_AsLong(v);
        break;
    }
}

void pyrameters::prm2ppy(PRM prm) {
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

bool pyrameters::ppy2prm(PRM prm) {
    bool bOK = true;

    for ( auto pn=prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
        //auto v = PyDict_GetItemString(arguments, pn->pszName);
        //auto f = PyDict_GetItemString(specified, pn->pszName);
        auto v = PyObject_GetAttrString(arguments, pn->pszName); // A Namespace
        if (v!=NULL) {
            if (v != Py_None) {
                auto f = PyObject_GetAttrString(specified, pn->pszName); // A Namespace
                if (f) {
                    pn->bArg = PyObject_IsTrue(f)>0;
                    Py_DECREF(f);
                }
                else pn->bArg = 0;
                if (PyList_Check(v)) {
                    if (pn->pCount==NULL) {
                        fprintf(stderr,"The parameter %s cannot be a list!\n",pn->pszName);
                        bOK = false;
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
        else {
            PyErr_Clear();
        }
    }
    return bOK;
}

template<> PyObject *pyrameters::get<PyObject *>(const char *name) {
    auto v = PyObject_GetAttrString(arguments, name);
    if (!v) throw std::domain_error(name);
    if (PyCallable_Check(v)) {
        auto callback = v;
        auto call_args = PyTuple_New(0);
        v = PyObject_Call(callback,call_args,dynamic);
        Py_DECREF(call_args);
        Py_DECREF(callback);
    }
    if (PyErr_Occurred()) PyErr_Print();
    return v;
}

template<> double pyrameters::get<double>(const char *name, PyObject *v) {
    if (PyFloat_Check(v)) return PyFloat_AsDouble(v);
    else if (PyLong_Check(v)) return PyLong_AsLong(v);
    else throw std::domain_error(name);
}
template<> double pyrameters::get<double>(const char *name) {
    double result;
    auto v = get<PyObject *>(name);
    result = get<double>(name,v);
    Py_DECREF(v);
    return result;
}

template<> std::int64_t pyrameters::get<std::int64_t>(const char *name, PyObject *v) {
    if (PyLong_Check(v)) return PyLong_AsLong(v);
    else if (PyFloat_Check(v)) return static_cast<std::int64_t>(PyFloat_AsDouble(v));
    else throw std::domain_error(name);
}
template<> std::int64_t pyrameters::get<std::int64_t>(const char *name) {
    std::int64_t result;
    auto v = get<PyObject *>(name);
    result = get<std::int64_t>(name,v);
    Py_DECREF(v);
    return result;
}

template<> std::string pyrameters::get<std::string>(const char *name, PyObject *v) {
    std::string result;
    if (PyUnicode_Check(v)) {
        auto ascii = PyUnicode_AsASCIIString(v);
        result = PyBytes_AsString(ascii);
        Py_DECREF(ascii);
    }
    else throw std::domain_error(name);
    return result;
}
template<> std::string pyrameters::get<std::string>(const char *name) {
    std::string result;
    auto v = get<PyObject *>(name);
    result = get<std::string>(name,v);
    Py_DECREF(v);
    return result;
}

template<> bool pyrameters::get<bool>(const char *name) {
    bool result = false;
    auto v = get<PyObject *>(name);
    result = PyObject_IsTrue(v)>0;
    Py_DECREF(v);
    return result;
}

template<> void pyrameters::set_dynamic(const char *name, double value) {
    auto o = PyFloat_FromDouble(value);
    if (o) {
        PyDict_SetItemString(dynamic,name,o);
        Py_DECREF(o);
    }
    if (PyErr_Occurred()) {
        PyErr_Print();
        abort();
    }
}

template<> void pyrameters::set_dynamic(const char *name, std::int64_t value) {
    auto o = PyLong_FromSsize_t(value);
    if (o) {
        PyDict_SetItemString(dynamic,name,o);
        Py_DECREF(o);
    }
    if (PyErr_Occurred()) {
        PyErr_Print();
        abort();
    }
}

template<> void pyrameters::set_dynamic(const char *name, std::uint64_t value) {
    auto o = PyLong_FromSize_t(value);
    if (o) {
        PyDict_SetItemString(dynamic,name,o);
        Py_DECREF(o);
    }
    if (PyErr_Occurred()) {
        PyErr_Print();
        abort();
    }
}

template<> void pyrameters::set_dynamic(const char *name, float value) {
    set_dynamic<double>(name,value);
}

template<> void pyrameters::set_dynamic(const char *name, std::int32_t value) {
    set_dynamic<std::int64_t>(name,value);
}

template<> void pyrameters::set_dynamic(const char *name, std::uint32_t value) {
    set_dynamic<std::uint64_t>(name,value);
}
