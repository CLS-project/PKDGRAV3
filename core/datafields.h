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

#ifndef CORE_DATAFIELDS_H
#define CORE_DATAFIELDS_H

#include <cstddef>
#include <cstring> // std::memcpy
#include <limits>
#include <tuple>
#include <vector>
#include <type_traits>
#include <concepts>
#include <string_view>
#include "blitz/array.h"

namespace datafields {

// Forward declaration — FieldBuilder (defined in databuilder.h) is friended by Fields.
template <typename DATA, typename FIELD> class FieldBuilder;

// ---- Fields<DATA, FIELD> -------------------------------------------------
//! \brief Lean access class — the runtime half of the field-layout system.
//!
//! Holds only the computed layout array and the swap buffer.  The layout is a
//! single blitz::TinyVector<offset_type, field_count+1> where indices
//! 0..field_count-1 are per-field byte offsets (NO_OFFSET if not present) and
//! index field_count stores the total element size.
//!
//! No build-time state; no add() or commit().  This is the base class
//! inherited by dataStore and splitStore.
//!
//! A Fields object can be initialised in two ways:
//!
//! (a) Implicit conversion from a completed FieldBuilder (see databuilder.h):
//! @code
//!   FieldBuilder<PARTICLE, PKD_FIELD> b;
//!   b.add<double[3]>(PKD_FIELD::oPosition, "r");
//!   b.commit();
//!   particleStore ps;
//!   ps.set_layout(b);
//! @endcode
//!
//! (b) From a serialised layout (MPI broadcast / worker-thread init):
//! @code
//!   MPI_Bcast(b.build().get_layout().data(),
//!             (PKD_FIELD::MAX_FIELD + 1) * sizeof(size_t), MPI_BYTE, 0, comm);
//!   // on worker:
//!   worker_store.set_layout(recv_layout);
//! @endcode
template <typename DATA, typename FIELD>
class Fields {
    template <typename D, typename F> friend class FieldBuilder;

  public:
    using offset_type = size_t;
    static constexpr int field_count = static_cast<int>(FIELD::MAX_FIELD);
    static constexpr offset_type NO_OFFSET = std::numeric_limits<offset_type>::max();
    using field = FIELD;
    //! Single flat array encoding the full layout: offsets[0..field_count-1]
    //! are per-field byte offsets; layout[field_count] is the element size.
    using layout_type = blitz::TinyVector<offset_type, field_count + 1>;

  protected:
    layout_type layout;            // layout[0..field_count-1] = offsets; layout[field_count] = element size
    std::vector<char> swap_buffer; // Buffer for element swapping

  public:
    Fields() : layout(NO_OFFSET) { layout[field_count] = sizeof(DATA); }

    //! Initialise from a serialised layout — used after MPI broadcast or to
    //! set up worker threads from a layout computed on the main thread.
    void set_layout(const layout_type& new_layout) {
        layout = new_layout;
        swap_buffer.resize(layout[field_count]);
    }

    //! Raw layout array for serialisation / MPI broadcast.
    auto get_layout() const -> const layout_type& { return layout; }

    //! Returns the offset (in bytes) of the specified field, or NO_OFFSET if not present.
    auto offset(FIELD f) const noexcept { return layout[static_cast<unsigned int>(f)]; }
    //! Returns true if the field is present, otherwise false.
    auto present(FIELD f) const noexcept { return layout[static_cast<unsigned int>(f)] != NO_OFFSET; }

    //! Returns a reference to the scalar field of type T within element *p.
    template <typename T, std::enable_if_t<!std::is_array_v<T>, bool> = true>
    const auto& get(const DATA* p, FIELD f) const noexcept {
        auto v = reinterpret_cast<const char*>(p);
        return *reinterpret_cast<const T*>(v + layout[static_cast<unsigned int>(f)]);
    }
    template <typename T, std::enable_if_t<!std::is_array_v<T>, bool> = true>
    auto& get(DATA* p, FIELD f) const noexcept {
        auto v = reinterpret_cast<char*>(p);
        return *reinterpret_cast<T*>(v + layout[static_cast<unsigned int>(f)]);
    }

    //! Array-type get — returns a blitz++ TinyVector reference.
    template <typename T, std::enable_if_t<std::rank_v<T> == 1, bool> = true>
    const auto& get(const DATA* p, FIELD f) const noexcept {
        return get<blitz::TinyVector<typename std::remove_extent<T>::type, std::extent_v<T>>>(p, f);
    }
    template <typename T, std::enable_if_t<std::rank_v<T> == 1, bool> = true>
    auto& get(DATA* p, FIELD f) const noexcept {
        return get<blitz::TinyVector<typename std::remove_extent<T>::type, std::extent_v<T>>>(p, f);
    }

    void swap(void* f1, void* f2) noexcept {
        auto elem_size = layout[field_count];
        std::memcpy(swap_buffer.data(), f1, elem_size);
        std::memcpy(f1, f2, elem_size);
        std::memcpy(f2, swap_buffer.data(), elem_size);
    }

    //! Returns the size (in bytes) of one element (base struct + all present fields).
    auto ElementSize() const noexcept { return layout[field_count]; }

    //! Returns a typed pointer to the i'th element starting at pBase.
    DATA* Element(void* pBase, int i) const noexcept {
        auto v = static_cast<char*>(pBase);
        return reinterpret_cast<DATA*>(v + (static_cast<uint64_t>(i) * layout[field_count]));
    }
};

// ---- Type trait helpers --------------------------------------------------

// Checks if a type is a std::tuple
template <class>
struct is_tuple : std::false_type {};
template <class... Ts>
struct is_tuple<std::tuple<Ts...>> : std::true_type {};
template <class T>
inline constexpr bool is_tuple_v = is_tuple<T>::value;

// Checks if a type is a TinyVector
template <class>
struct is_tinyvector : std::false_type {};
template <class T, int N>
struct is_tinyvector<blitz::TinyVector<T, N>> : std::true_type {};
template <class T>
inline constexpr bool is_tinyvector_v = is_tinyvector<std::remove_cvref_t<T>>::value;

// Concept that accepts a scalar type, a TinyVector, or a type with a tuple schema
// (e.g., a struct defined with MAKE_STRUCTURE)
template <class T>
concept ScalarOrHasTupleSchema = std::is_scalar_v<T> || is_tinyvector_v<T> || requires {
    { std::remove_cvref_t<T>::schema };
} && is_tuple_v<std::remove_cvref_t<decltype(std::remove_cvref_t<T>::schema)>>;

template <ScalarOrHasTupleSchema T>
struct Field {
    std::string_view name; // Name of the field (the structure member name as a string)
    std::size_t offset;    // Offset of the field within the structure (in bytes)
    using type = T;
};

} // namespace datafields

// For example:
// #define DATA_FIELDS(X) X(x,float) X(y,float) X(i,int)
// MAKE_STRUCTURE(DATA, DATA_FIELDS)
// Results in:
// namespace schema_details {
//     struct DATA {
//         float x;
//         float y;
//         int i;
//    };
// }
// struct DATA : public schema_details::NAME {
//     using base_type = schema_details::DATA;
//     static constexpr auto schema = std::tuple{datafields::Field<float>{"x"}, datafields::Field<float>{"y"}, datafields::Field<int>{"i"},};
// };
// It is not an accident that the type comes second and that we use __VA_ARGS__.
// This allows templated types (e.g., TinyVector) to be used as field types without needing to worry about commas (extra parameters technically)
#define STRUCTURE_FIELD_(name, ...) __VA_ARGS__ name;
#define STRUCTURE_FIELD_ELEM_(name, ...) datafields::Field<__VA_ARGS__>{#name, offsetof(base_type, name)},

#define MAKE_STRUCTURE(NAME, FIELDS)                                                                                                                           \
    namespace schema_details {                                                                                                                                 \
    struct NAME {                                                                                                                                              \
        FIELDS(STRUCTURE_FIELD_)                                                                                                                               \
    };                                                                                                                                                         \
    }                                                                                                                                                          \
    struct NAME : public schema_details::NAME {                                                                                                                \
        using base_type = schema_details::NAME;                                                                                                                \
        static constexpr auto schema = std::tuple{FIELDS(STRUCTURE_FIELD_ELEM_)};                                                                              \
    };                                                                                                                                                         \
    static_assert(std::is_standard_layout_v<NAME>);
#endif
