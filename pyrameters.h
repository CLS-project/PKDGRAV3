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
#include <string_view>
#include <cstdint>
#include <vector>
#include "blitz/array.h"

#include "param.h"  // This should go away at some point

class pyrameters {
    template<typename T, typename _ = void>
    struct is_container : std::false_type {};

    template<typename T>
    struct is_container<T, std::void_t<typename T::value_type>> : std::true_type {};

    // A helper template to determine if a type is a std::vector of some kind.
    template<typename>
    struct is_std_vector : std::false_type {};

    template<typename U>
    struct is_std_vector<std::vector<U>> : std::true_type {};

    template<typename>
    struct is_blitz_tinyvector : std::false_type {};

    template<typename U, int N>
    struct is_blitz_tinyvector<blitz::TinyVector<U, N>> : std::true_type {};

protected:
    PyObject *arguments_=nullptr, *specified_=nullptr;
    PyObject *dynamic_=nullptr;

public:
    auto arguments() {Py_INCREF(arguments_); return arguments_; }
    auto specified() {Py_INCREF(specified_); return specified_; }

protected:

    bool has(const char *name) const;

    // Generic get function
    template<typename T>
    T get(const char *name) const {
        auto obj = get<PyObject *>(name);  // Look up the attribute (and call if required); specialized
        T result = get<T>(name, obj);      // Convert the PyObject to the requested type
        Py_XDECREF(obj);
        return result;
    }

    template<typename T>
    T get(const char *name, PyObject *v) const {
        if constexpr (is_blitz_tinyvector<T>::value) {
            if (PyList_Check(v)) {
                Py_ssize_t size = PyList_Size(v);
                T result;
                if (result.length() != size) throw std::domain_error(name);
                for (Py_ssize_t i = 0; i < size; ++i) {
                    PyObject *item = PyList_GetItem(v, i);
                    result[i] = get<typename T::T_numtype>(name, item);
                }
                return result;
            }
            else return T(get<typename T::T_numtype>(name,v));
        }
        else if constexpr (std::is_same<T, std::string_view>::value) {
            if (PyUnicode_Check(v)) {
                Py_ssize_t length;
                auto c_string = PyUnicode_AsUTF8AndSize(v,&length);         // convert to a UTF8 string
                if (c_string == nullptr) throw std::domain_error(name);     // must be a string
                return std::string_view(c_string,length);                   // return as a string_view
            }
            else if (v == Py_None) return std::string_view();
            else throw std::domain_error(name);
        }
        else if constexpr (is_std_vector<T>::value) {
            if (PyList_Check(v)) {
                Py_ssize_t size = PyList_Size(v);
                T result(size);
                for (Py_ssize_t i = 0; i < size; ++i) {
                    PyObject *item = PyList_GetItem(v, i);
                    result[i] = get<typename T::value_type>(name, item);
                }
                return result;
            }
        }
        else if constexpr (std::is_same<T, bool>::value) {
            return PyObject_IsTrue(v)>0;
        }
        else if constexpr (std::is_integral<T>::value) {
            if (PyLong_Check(v)) return PyLong_AsLong(v);
            else if (PyFloat_Check(v)) return static_cast<std::int64_t>(PyFloat_AsDouble(v));
        }
        else if constexpr (std::is_floating_point<T>::value) {
            if (PyFloat_Check(v)) return PyFloat_AsDouble(v);
            else if (PyLong_Check(v)) return PyLong_AsLong(v);
        }
        else {
            throw std::runtime_error("Unsupported type for get");
        }
        throw std::runtime_error("Failed to convert PyObject to requested type.");
    }

protected:
    template<typename T>
    auto get_value(const char *name, const T &value) {
        if constexpr (std::is_same<T, bool>::value) {
            return Py_NewRef(value ? Py_True : Py_False);
        }
        else if constexpr (std::is_integral<T>::value) {
            if (std::is_signed<T>::value) {
                return PyLong_FromSsize_t(value);
            }
            else if constexpr (std::is_unsigned<T>::value) {
                return PyLong_FromSize_t(value);
            }
        }
        else if constexpr (std::is_floating_point<T>::value) {
            return PyFloat_FromDouble(value);
        }
        else if constexpr (std::is_same<T, const char *>::value) {
            return PyUnicode_FromString(value);
        }
        else if constexpr (std::is_same<T, std::string_view>::value) {
            return PyUnicode_FromString(value.data());
        }
        else if constexpr (std::is_same<T, PyObject *>::value) {
            return Py_NewRef(value);
        }
        else if constexpr (is_blitz_tinyvector<T>::value) {
            auto py_object = PyList_New(value.length());
            for (int i = 0; i < value.length(); ++i) {
                auto item = get_value(name,value[i]);
                PyList_SetItem(py_object,i,item);
            }
            return py_object;
        }
        else if constexpr (is_std_vector<T>::value) {
            auto py_object = PyList_New(value.size());
            for (int i = 0; i < value.size(); ++i) {
                auto item = get_value(name,value[i]);
                PyList_SetItem(py_object,i,item);
            }
            return py_object;
        }
        else {
            static_assert(std::is_same_v<T, void>, "Unsupported type for set");
        }
    }

public:
    template<typename T>
    void set(const char *name, const T &value) {
        auto py_object = get_value(name,value);
        PyObject_SetAttrString(arguments_,name,py_object);
        Py_DECREF(py_object);
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

template<> PyObject    *pyrameters::get<PyObject *>(const char *name) const;

template<> void pyrameters::set_dynamic(const char *name, float         value);
template<> void pyrameters::set_dynamic(const char *name, double        value);
template<> void pyrameters::set_dynamic(const char *name, std::int32_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::int64_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::uint32_t value);
template<> void pyrameters::set_dynamic(const char *name, std::uint64_t value);

#endif
