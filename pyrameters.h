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

#ifndef PYRAMETERS_H
#define PYRAMETERS_H
#include <Python.h>
#include <stdexcept>
#include <string>

#include "param.h"  // This should go away at some point

class pyrameters {
protected:
    PyObject *arguments=nullptr, *specified=nullptr;

public:
    auto get_arguments() {return arguments; }
    auto get_specified() {return specified; }

protected:

    bool has(const char *name) const {
        bool bSpecified = false;
        if (auto f = PyObject_GetAttrString(specified,name)) {
            bSpecified = PyObject_IsTrue(f)>0;
            Py_DECREF(f);
        }
        return bSpecified;
    }

    template<typename T> T get(const char *name);
    template<typename T> T get(const char *name, PyObject *v);

public:

    virtual ~pyrameters() {
        // If we have arguments and specified then we need to decrement the reference count
        if (Py_IsInitialized()) {
            Py_XDECREF(arguments);
            Py_XDECREF(specified);
        }
    }

    pyrameters() : arguments(nullptr), specified(nullptr) {}

    pyrameters(PyObject *arguments, PyObject *specified) : arguments(arguments), specified(specified) {
        Py_INCREF(this->arguments);
        Py_INCREF(this->specified);
    }

    pyrameters(const pyrameters &other) {           // copy constructor
        this->arguments = other.arguments;
        this->specified = other.specified;
        Py_INCREF(this->arguments);
        Py_INCREF(this->specified);
    }

    pyrameters &operator=(const pyrameters &rhs) {  // copy assignment
        if (this != &rhs) {
            Py_XDECREF(this->arguments);
            Py_XDECREF(this->specified);

            this->arguments = rhs.arguments;
            this->specified = rhs.specified;
            Py_INCREF(this->arguments);
            Py_INCREF(this->specified);
        }
        return *this;
    }

    pyrameters(pyrameters &&other) {                // move constructor
        this->arguments = other.arguments;
        this->specified = other.specified;
        other.arguments = nullptr;
        other.specified = nullptr;
    }

    // move assignment
    pyrameters &operator=(pyrameters &&rhs) {       // move assignment
        if (this != &rhs) {
            this->arguments = rhs.arguments;
            this->specified = rhs.specified;
            rhs.arguments = nullptr;
            rhs.specified = nullptr;
        }
        return *this;
    }

protected:
    static void merge(PyObject *o1, PyObject *o2) {
        auto d1 = PyObject_GenericGetDict(o1,nullptr);
        auto d2 = PyObject_GenericGetDict(o2,nullptr);
        PyObject *key, *value;
        Py_ssize_t pos = 0;
        while (PyDict_Next(d2, &pos, &key, &value)) {
            PyDict_SetItem(d1,key,value);
        }
    }

public:
    /// @brief Merge the parameters from another parameter set into this one
    /// @param other the other parameter set
    void merge(const pyrameters &other) {
        merge(this->arguments,other.arguments);
        merge(this->specified,other.specified);
    }

    /// @brief Given a dictionary like object make sure that the keys exist in our parameters
    /// @param kwobj a dictionary like list of keys
    /// @return true if there are no rogue variables, false otherwise
    bool verify(PyObject *kwobj) {
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
            if (!PyObject_HasAttr(arguments,key) && !PyModule_Check(value)) {
                PyErr_Format(PyExc_AttributeError,"invalid parameter %A",key);
                PyErr_Print();
                bSuccess=false;
            }
        }
        return bSuccess;
    }

    bool update(PyObject *kwobj,bool bIgnoreUnknown=false) {
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

// This interface will go away at some point
protected:
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

public:
    void prm2ppy(PRM prm) {
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

    bool ppy2prm(PRM prm) {
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

};

template<> inline PyObject *pyrameters::get<PyObject *>(const char *name) {
    auto v = PyObject_GetAttrString(arguments, name);
    if (!v) throw std::domain_error(name);
    return v;
}

template<> inline double pyrameters::get<double>(const char *name, PyObject *v) {
    if (PyFloat_Check(v)) return PyFloat_AsDouble(v);
    else if (PyLong_Check(v)) return PyLong_AsLong(v);
    else throw std::domain_error(name);
}
template<> inline double pyrameters::get<double>(const char *name) {
    double result;
    auto v = get<PyObject *>(name);
    result = get<double>(name,v);
    Py_DECREF(v);
    return result;
}

template<> inline int64_t pyrameters::get<int64_t>(const char *name, PyObject *v) {
    if (PyLong_Check(v)) return PyLong_AsLong(v);
    else if (PyFloat_Check(v)) return static_cast<int64_t>(PyFloat_AsDouble(v));
    else throw std::domain_error(name);
}
template<> inline int64_t pyrameters::get<int64_t>(const char *name) {
    int64_t result;
    auto v = get<PyObject *>(name);
    result = get<int64_t>(name,v);
    Py_DECREF(v);
    return result;
}

template<> inline std::string pyrameters::get<std::string>(const char *name, PyObject *v) {
    std::string result;
    if (PyUnicode_Check(v)) {
        auto ascii = PyUnicode_AsASCIIString(v);
        result = PyBytes_AsString(ascii);
        Py_DECREF(ascii);
    }
    else throw std::domain_error(name);
    return result;
}
template<> inline std::string pyrameters::get<std::string>(const char *name) {
    std::string result;
    auto v = get<PyObject *>(name);
    result = get<std::string>(name,v);
    Py_DECREF(v);
    return result;
}

template<> inline bool pyrameters::get<bool>(const char *name) {
    bool v = false;
    if (auto o = PyObject_GetAttrString(arguments,name)) {
        v = PyObject_IsTrue(o)>0;
        Py_DECREF(o);
    }
    if (PyErr_Occurred()) {
        PyErr_Print();
        abort();
    }
    return v;
}

#endif
