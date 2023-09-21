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
#include <cstdint>
#include "blitz/array.h"

#include "param.h"  // This should go away at some point

class pyrameters {
protected:
    PyObject *arguments_=nullptr, *specified_=nullptr;
    PyObject *dynamic_=nullptr;

public:
    auto arguments() {Py_INCREF(arguments_); return arguments_; }
    auto specified() {Py_INCREF(specified_); return specified_; }

protected:

    bool has(const char *name) const;

    template<typename T> T get(const char *name);
    template<typename T> T get(const char *name, PyObject *v);
    template<typename T,int N> blitz::TinyVector<T,N> get(const char *name);

public:
    template<typename T, typename std::enable_if<std::is_integral<T>::value &&std::is_signed<T>::value&& !std::is_same<T, bool>::value, int>::type = 0>
    void set(const char *name, T value) {
        PyObject_SetAttrString(arguments_,name,PyLong_FromSsize_t(value));
    }
    template<typename T, typename std::enable_if<std::is_integral<T>::value &&std::is_unsigned<T>::value&& !std::is_same<T, bool>::value, int>::type = 0>
    void set(const char *name, T value) {
        PyObject_SetAttrString(arguments_,name,PyLong_FromSize_t(value));
    }
    template<typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 0>
    void set(const char *name, T value) {
        PyObject_SetAttrString(arguments_,name,PyFloat_FromDouble(value));
    }
    template<typename T, typename std::enable_if<std::is_same<T, bool>::value, int>::type = 0>
    void set(const char *name, T value) {
        PyObject_SetAttrString(arguments_,name,value ? Py_True : Py_False);
    }

public:

    virtual ~pyrameters() {
        // If we have arguments and specified then we need to decrement the reference count
        if (Py_IsInitialized()) {
            Py_XDECREF(arguments_);
            Py_XDECREF(specified_);
            Py_XDECREF(dynamic_);
        }
    }

    pyrameters() : arguments_(nullptr), specified_(nullptr), dynamic_(nullptr) {}

    pyrameters(PyObject *arguments, PyObject *specified) : arguments_(arguments), specified_(specified), dynamic_(PyDict_New()) {
        Py_INCREF(this->arguments_);
        Py_INCREF(this->specified_);
    }

    pyrameters(const pyrameters &other) {           // copy constructor
        this->arguments_ = other.arguments_;
        this->specified_ = other.specified_;
        this->dynamic_ = PyDict_Copy(other.dynamic_);
        Py_INCREF(this->arguments_);
        Py_INCREF(this->specified_);
    }

    pyrameters &operator=(const pyrameters &rhs) {  // copy assignment
        if (this != &rhs) {
            Py_XDECREF(this->arguments_);
            Py_XDECREF(this->specified_);
            Py_XDECREF(this->dynamic_);

            this->arguments_ = rhs.arguments_;
            this->specified_ = rhs.specified_;
            this->dynamic_ = rhs.dynamic_;
            Py_INCREF(this->arguments_);
            Py_INCREF(this->specified_);
            Py_INCREF(this->dynamic_);
        }
        return *this;
    }

    pyrameters(pyrameters &&other) {                // move constructor
        this->arguments_ = other.arguments_;
        this->specified_ = other.specified_;
        this->dynamic_ = other.dynamic_;
        other.arguments_ = nullptr;
        other.specified_ = nullptr;
        other.dynamic_ = nullptr;
    }

    // move assignment
    pyrameters &operator=(pyrameters &&rhs) {       // move assignment
        if (this != &rhs) {
            this->arguments_ = rhs.arguments_;
            this->specified_ = rhs.specified_;
            this->dynamic_ = rhs.dynamic_;
            rhs.arguments_ = nullptr;
            rhs.specified_ = nullptr;
            rhs.dynamic_ = nullptr;
        }
        return *this;
    }

protected:
    static void merge(PyObject *o1, PyObject *o2);
    PyObject *call_or_return(PyObject *value);

public:
    /// @brief Merge the parameters from another parameter set into this one
    /// @param other the other parameter set
    void merge(const pyrameters &other);

    /// @brief Given a dictionary like object make sure that the keys exist in our parameters
    /// @param kwobj a dictionary like list of keys
    /// @return true if there are no rogue variables, false otherwise
    bool verify(PyObject *kwobj);

    bool update(PyObject *kwobj,bool bIgnoreUnknown=false);

public:
    void prm2ppy(PRM prm);
    bool ppy2prm(PRM prm);

public:
    template<typename T> void set_dynamic(const char *name, T value);

};

template<> PyObject    *pyrameters::get<PyObject *>(const char *name);
template<> double       pyrameters::get<double>(const char *name, PyObject *v);
template<> double       pyrameters::get<double>(const char *name);
template<> std::int64_t pyrameters::get<std::int64_t>(const char *name, PyObject *v);
template<> std::int64_t pyrameters::get<std::int64_t>(const char *name);
template<> std::string  pyrameters::get<std::string>(const char *name, PyObject *v);
template<> std::string  pyrameters::get<std::string>(const char *name);
template<> bool         pyrameters::get<bool>(const char *name);

template<> void pyrameters::set_dynamic(const char *name, float         value);
template<> void pyrameters::set_dynamic(const char *name, double        value);
template<> void pyrameters::set_dynamic(const char *name, std::int32_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::int64_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::uint32_t value);
template<> void pyrameters::set_dynamic(const char *name, std::uint64_t value);

#endif
