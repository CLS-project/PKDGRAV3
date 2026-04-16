// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) Douglas Potter and Joachim Stadel
#pragma once

#include "datafields.h" // datafields::is_tinyvector_v, MAKE_STRUCTURE schema
#include <Python.h>
#include <numpy/arrayobject.h>
#include <cstdint>     // int8_t … uint64_t
#include <tuple>       // std::tuple_size_v, std::get
#include <type_traits> // std::remove_cvref_t, std::void_t, std::integral_constant
#include <utility>     // std::index_sequence, std::make_index_sequence

// ---------------------------------------------------------------------------
// Type helpers
// ---------------------------------------------------------------------------

template <class T>
using bare_t = std::remove_cvref_t<T>;

// ---------------------------------------------------------------------------
// numpy_typenum — maps a C scalar type to the corresponding NumPy type number.
// ---------------------------------------------------------------------------

template <class T>
struct numpy_typenum;

template <> struct numpy_typenum<float> : std::integral_constant<int, NPY_FLOAT32> {};
template <> struct numpy_typenum<double> : std::integral_constant<int, NPY_FLOAT64> {};
template <> struct numpy_typenum<int8_t> : std::integral_constant<int, NPY_INT8> {};
template <> struct numpy_typenum<uint8_t> : std::integral_constant<int, NPY_UINT8> {};
template <> struct numpy_typenum<int16_t> : std::integral_constant<int, NPY_INT16> {};
template <> struct numpy_typenum<uint16_t> : std::integral_constant<int, NPY_UINT16> {};
template <> struct numpy_typenum<int32_t> : std::integral_constant<int, NPY_INT32> {};
template <> struct numpy_typenum<uint32_t> : std::integral_constant<int, NPY_UINT32> {};
template <> struct numpy_typenum<int64_t> : std::integral_constant<int, NPY_INT64> {};
template <> struct numpy_typenum<uint64_t> : std::integral_constant<int, NPY_UINT64> {};
template <> struct numpy_typenum<bool> : std::integral_constant<int, NPY_BOOL> {};

// ---------------------------------------------------------------------------
// has_schema — detects types produced by MAKE_STRUCTURE (they carry a schema
// member tuple describing their fields).
// ---------------------------------------------------------------------------

template <class, class = void>
struct has_schema : std::false_type {};
template <class T>
struct has_schema<T, std::void_t<decltype(T::schema)>> : std::true_type {};
template <class T>
inline constexpr bool has_schema_v = has_schema<T>::value;

// ---------------------------------------------------------------------------
// tinyvector_traits — extracts element type and extent from a TinyVector.
// ---------------------------------------------------------------------------

template <class>
struct tinyvector_traits;
template <class T, int N>
struct tinyvector_traits<blitz::TinyVector<T, N>> {
    using value_type = T;
    static constexpr int extent = N;
};

// ---------------------------------------------------------------------------
// get_numpy_dtype<T> — returns a new NumPy dtype descriptor for T.
// Handles scalars, TinyVectors (sub-arrays), and MAKE_STRUCTURE structs.
// ---------------------------------------------------------------------------

template <class T>
auto get_numpy_dtype() -> PyArray_Descr*; // forward declaration

template <class T>
auto get_numpy_scalar_dtype() -> PyArray_Descr* {
    return PyArray_DescrFromType(numpy_typenum<T>::value);
}

template <class TV>
auto get_numpy_tinyvector_dtype() -> PyArray_Descr* {
    using traits = tinyvector_traits<bare_t<TV>>;
    using value_type = typename traits::value_type;
    constexpr int extent = traits::extent;

    PyArray_Descr* base = get_numpy_dtype<value_type>();
    if (!base) {
        return nullptr;
    }

    // Build the Python equivalent of: np.dtype((base, (extent,)))
    PyObject* shape = PyTuple_New(1);
    if (!shape) {
        Py_DECREF(reinterpret_cast<PyObject*>(base));
        return nullptr;
    }

    PyObject* n_obj = PyLong_FromLong(extent);
    if (!n_obj) {
        Py_DECREF(shape);
        Py_DECREF(reinterpret_cast<PyObject*>(base));
        return nullptr;
    }
    PyTuple_SET_ITEM(shape, 0, n_obj); // steals n_obj

    PyObject* spec = PyTuple_New(2);
    if (!spec) {
        Py_DECREF(shape);
        Py_DECREF(reinterpret_cast<PyObject*>(base));
        return nullptr;
    }
    PyTuple_SET_ITEM(spec, 0, reinterpret_cast<PyObject*>(base)); // steals base
    PyTuple_SET_ITEM(spec, 1, shape);                             // steals shape

    PyArray_Descr* out = nullptr;
    if (!PyArray_DescrConverter(spec, &out)) {
        Py_DECREF(spec);
        return nullptr;
    }
    Py_DECREF(spec);
    return out;
}

template <class T, std::size_t... I>
auto get_numpy_struct_dtype_impl(std::index_sequence<I...> /*indices*/) -> PyArray_Descr* {
    constexpr auto const& schema = T::schema;
    constexpr std::size_t field_count = std::tuple_size_v<bare_t<decltype(schema)>>;

    PyObject* names = PyList_New(field_count);
    PyObject* formats = PyList_New(field_count);
    PyObject* offsets = PyList_New(field_count);
    PyObject* spec = PyDict_New();

    if (!names || !formats || !offsets || !spec) {
        Py_XDECREF(names);
        Py_XDECREF(formats);
        Py_XDECREF(offsets);
        Py_XDECREF(spec);
        return nullptr;
    }

    bool success = true;

    auto add_field = [&](auto index_c) -> void {
        if (!success) {
            return;
        }

        constexpr std::size_t idx = decltype(index_c)::value;
        constexpr auto field_info = std::get<idx>(schema);
        using field_type = typename decltype(field_info)::type;

        PyObject* py_name = PyUnicode_FromStringAndSize(field_info.name.data(), static_cast<Py_ssize_t>(field_info.name.size()));
        if (!py_name) {
            success = false;
            return;
        }

        PyArray_Descr* subdescr = get_numpy_dtype<field_type>();
        if (!subdescr) {
            Py_DECREF(py_name);
            success = false;
            return;
        }

        PyObject* py_offset = PyLong_FromSize_t(field_info.offset);
        if (!py_offset) {
            Py_DECREF(py_name);
            Py_DECREF(reinterpret_cast<PyObject*>(subdescr));
            success = false;
            return;
        }

        PyList_SET_ITEM(names, idx, py_name);                                 // steals
        PyList_SET_ITEM(formats, idx, reinterpret_cast<PyObject*>(subdescr)); // steals
        PyList_SET_ITEM(offsets, idx, py_offset);                             // steals
    };

    (add_field(std::integral_constant<std::size_t, I>{}), ...);

    if (!success) {
        Py_DECREF(names);
        Py_DECREF(formats);
        Py_DECREF(offsets);
        Py_DECREF(spec);
        return nullptr;
    }

    PyObject* itemsize = PyLong_FromSize_t(sizeof(T));
    if (!itemsize) {
        Py_DECREF(names);
        Py_DECREF(formats);
        Py_DECREF(offsets);
        Py_DECREF(spec);
        return nullptr;
    }

    bool dict_success = PyDict_SetItemString(spec, "names", names) >= 0 && PyDict_SetItemString(spec, "formats", formats) >= 0 &&
                        PyDict_SetItemString(spec, "offsets", offsets) >= 0 && PyDict_SetItemString(spec, "itemsize", itemsize) >= 0;

    Py_DECREF(names);
    Py_DECREF(formats);
    Py_DECREF(offsets);
    Py_DECREF(itemsize);

    if (!dict_success) {
        Py_DECREF(spec);
        return nullptr;
    }

    PyArray_Descr* out = nullptr;
    if (!PyArray_DescrConverter(spec, &out)) {
        Py_DECREF(spec);
        return nullptr;
    }
    Py_DECREF(spec);
    return out;
}

template <class T>
auto get_numpy_dtype() -> PyArray_Descr* {
    using U = bare_t<T>;

    if constexpr (std::rank_v<U> == 1) {
        // C array — treated as a fixed-length sub-array, same as a TinyVector.
        using elem_t = std::remove_extent_t<U>;
        constexpr int extent = static_cast<int>(std::extent_v<U>);
        return get_numpy_tinyvector_dtype<blitz::TinyVector<elem_t, extent>>();
    }
    else if constexpr (datafields::is_tinyvector_v<U>) {
        return get_numpy_tinyvector_dtype<U>();
    }
    else if constexpr (has_schema_v<U>) {
        constexpr std::size_t N = std::tuple_size_v<bare_t<decltype(U::schema)>>;
        return get_numpy_struct_dtype_impl<U>(std::make_index_sequence<N>{});
    }
    else if constexpr (requires { numpy_typenum<U>::value; }) {
        return get_numpy_scalar_dtype<U>();
    }
    else {
        return nullptr; // Type not representable as a numpy dtype — field will be skipped.
    }
}

// ---------------------------------------------------------------------------
// get_python_type<T> — top-level entry point.
// Returns a new reference to the NumPy dtype object for T.
// ---------------------------------------------------------------------------

template <class T>
auto get_python_type() -> PyObject* {
    return reinterpret_cast<PyObject*>(get_numpy_dtype<T>());
}
