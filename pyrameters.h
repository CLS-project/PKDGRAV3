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

class pyrameters {
protected:
    PyObject *arguments=nullptr, *specified=nullptr;
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

    template<> PyObject *get<PyObject *>(const char *name) {
        auto v = PyObject_GetAttrString(arguments, name);
        if (!v) throw std::domain_error(name);
        return v;
    }
    template<> double get<double>(const char *name, PyObject *v) {
        if (PyFloat_Check(v)) return PyFloat_AsDouble(v);
        else if (PyLong_Check(v)) return PyLong_AsLong(v);
        else throw std::domain_error(name);
    }
    template<> double get<double>(const char *name) {
        double result;
        auto v = get<PyObject *>(name);
        result = get<double>(name,v);
        Py_DECREF(v);
        return result;
    }

    template<> int64_t get<int64_t>(const char *name, PyObject *v) {
        if (PyLong_Check(v)) return PyLong_AsLong(v);
        else if (PyFloat_Check(v)) return static_cast<int64_t>(PyFloat_AsDouble(v));
        else throw std::domain_error(name);
    }
    template<> int64_t get<int64_t>(const char *name) {
        int64_t result;
        auto v = get<PyObject *>(name);
        result = get<int64_t>(name,v);
        Py_DECREF(v);
        return result;
    }

    template<> std::string get<std::string>(const char *name, PyObject *v) {
        std::string result;
        if (PyUnicode_Check(v)) {
            auto ascii = PyUnicode_AsASCIIString(v);
            result = PyBytes_AsString(ascii);
            Py_DECREF(ascii);
        }
        else throw std::domain_error(name);
        return result;
    }
    template<> std::string get<std::string>(const char *name) {
        std::string result;
        auto v = get<PyObject *>(name);
        result = get<std::string>(name,v);
        Py_DECREF(v);
        return result;
    }

    template<> bool get<bool>(const char *name) {
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

public:

    virtual ~pyrameters() {
        // If we have arguments and specified then we need to decrement the reference count
        Py_XDECREF(arguments);
        Py_XDECREF(specified);
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

};

#endif
