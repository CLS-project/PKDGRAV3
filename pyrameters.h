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

namespace pyrameters {
struct arg {
    template <typename T>
    PyObject *convert_to_pyobject(const T& val) {
        using DecayedT = std::decay_t<T>; // Removes references and converts arrays to pointers
        if constexpr (std::is_integral_v<DecayedT>) {
            if constexpr (std::is_unsigned_v<DecayedT>) return PyLong_FromUnsignedLongLong(static_cast<unsigned long long>(val));
            else                                        return PyLong_FromLongLong(static_cast<long long>(val));
        }
        else if constexpr (std::is_floating_point<DecayedT>::value)      return PyFloat_FromDouble(val);
        else if constexpr (std::is_same_v<DecayedT, std::string_view> || std::is_same_v<DecayedT, std::string>) {
            if (val.data() == nullptr) return Py_None;
            return PyUnicode_FromStringAndSize(val.data(), val.size());
        }
        else if constexpr (std::is_convertible_v<DecayedT, const char*>) {
            return (val != nullptr) ? PyUnicode_FromString(val) : Py_None;
        }
        else {
            static_assert(std::is_same_v<DecayedT, void>, "Unsupported type in convert_to_pyobject");
            return nullptr; // This line will never be reached
        }
    }
    const char *key;
    PyObject *value;

    arg(const char* k, PyObject* v) : key(k), value(v) {}
    template<typename T>
    arg(const char* k, const T &v) : key(k), value(convert_to_pyobject(v)) {}
};

struct ArgKey {
    const char *key;
    constexpr explicit ArgKey(const char* k) : key(k) {}
    template <typename T>
    arg operator=(const T& val) const {
        return arg(key, val);
    }
};

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

    PyObject *get_and_call(const char *name, PyObject *dynamic) const;

    // Generic get function
    template<typename T>
    T get(const char *name) const {
        auto obj = get_and_call(name,dynamic_); // Look up the attribute (and call if required)
        T result = convert_to_C<T>(name, obj);  // Convert the PyObject to the requested type
        Py_XDECREF(obj);
        return result;
    }

    template <typename T, typename... Args>
    T get(const char *name, arg first, Args... rest) const {
        PyObject* dynamic = PyDict_Copy(dynamic_);  // Start with a copy of dynamic_
        if (dynamic == nullptr) throw std::runtime_error("Failed to create dynamic dictionary.");
        add_to_dict(dynamic, first, rest...);       // Add all key-value pairs to temp_dict
        auto obj = get_and_call(name,dynamic);      // Look up the attribute (and call if required)
        T result = convert_to_C<T>(name, obj);      // Convert the PyObject to the requested type
        Py_XDECREF(obj);
        Py_DECREF(dynamic);
        return result;
    }

private:
    template<typename T>
    T convert_to_C(const char *name, PyObject *v) const {
        if constexpr (is_blitz_tinyvector<T>::value) {
            if (PySequence_Check(v)) { // Must be a sequence
                PyObject *seq = PySequence_Fast(v, name);
                if (!seq) throw std::domain_error(name);
                Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
                T result;
                if (result.length() != size) { // of a fixed size (like TinyVector)
                    Py_DECREF(seq);
                    throw std::domain_error(name);
                }
                PyObject **items = PySequence_Fast_ITEMS(seq);
                for (Py_ssize_t i = 0; i < size; ++i) {
                    result[i] = convert_to_C<typename T::T_numtype>(name, items[i]);
                }
                Py_DECREF(seq);
                return result;
            }
            else return T(convert_to_C<typename T::T_numtype>(name,v));
        }
//      The following works just fine, but returning a string_view is not a good idea because the underlying data may be destroyed.
//      else if constexpr (std::is_same<T, std::string_view>::value || std::is_same<T, std::string>::value) {
        else if constexpr (std::is_same<T, std::string>::value) {
            if (PyUnicode_Check(v)) {
                Py_ssize_t length;
                auto c_string = PyUnicode_AsUTF8AndSize(v, &length);         // convert to a UTF8 string
                if (c_string == nullptr) throw std::domain_error(name);     // must be a string
                return T(c_string,length);                   // return as a string_view
            }
            else if (v == Py_None) return T("");
            else throw std::domain_error(name);
        }
        else if constexpr (is_std_vector<T>::value) {
            T result;
            if (PySequence_Check(v)) { // Must be a sequence
                PyObject *seq = PySequence_Fast(v, name);
                if (!seq) throw std::domain_error(name);
                Py_ssize_t size = PySequence_Fast_GET_SIZE(seq);
                result.reserve(size);
                PyObject **items = PySequence_Fast_ITEMS(seq);
                for (Py_ssize_t i = 0; i < size; ++i) {
                    result.emplace_back(convert_to_C<typename T::value_type>(name, items[i]));
                }
                Py_DECREF(seq);
                return result;
            }
            else if (v != Py_None) result.push_back(convert_to_C<typename T::value_type>(name,v));
            return result; // uses NVRO or move semantics
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
#if PY_VERSION_HEX < 0x030a0000
    auto Py_NewRef(PyObject *object) {
        Py_INCREF(object);
        return object;
    }
#endif
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
        else if constexpr (std::is_same<T, std::string_view>::value || std::is_same<T, std::string>::value) {
            return PyUnicode_FromStringAndSize(value.data(), value.size());
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
    template<typename T> void set_dynamic(const char *name, T value);

    PyObject *call_dynamic(PyObject *callback) const;
private:
    PyObject *call_dynamic(PyObject *callback, PyObject *dynamic) const;
private:
    // Base case: Add a single argument to the dictionary
    void add_to_dict(PyObject* dict, arg arg) const {
        if (arg.value) {
            PyDict_SetItemString(dict, arg.key, arg.value);
            Py_DECREF(arg.value);  // Safe to decref since PyDict_SetItemString increments ref
        }
    }

    // Recursive case: Process multiple arguments
    template <typename... Args>
    void add_to_dict(PyObject* dict, arg first, Args... rest) const {
        add_to_dict(dict, first);
        add_to_dict(dict, rest...);
    }
};

template<> void pyrameters::set_dynamic(const char *name, PyObject     *value);
template<> void pyrameters::set_dynamic(const char *name, float         value);
template<> void pyrameters::set_dynamic(const char *name, double        value);
template<> void pyrameters::set_dynamic(const char *name, std::int32_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::int64_t  value);
template<> void pyrameters::set_dynamic(const char *name, std::uint32_t value);
template<> void pyrameters::set_dynamic(const char *name, std::uint64_t value);

namespace literals {

// Define a `_key` suffix operator (like fmt::literals)
constexpr ArgKey operator"" _key(const char* key, size_t) {
    return ArgKey(key);
}
} // namespace literals
} // namespace pyrameters

#endif
