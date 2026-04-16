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
#include <limits>
#include <tuple>
#include <vector>
#include <type_traits>
#include <concepts>
#include <algorithm>
#include <string_view>
#include "blitz/array.h"

namespace datafields {

static constexpr size_t log2(size_t n) {
    return (n <= 1) ? 0 : 1 + log2(n / 2);
}

// Forward declarations — Overlay method bodies need FieldBuilder complete (defined below).
template <typename DATA, typename FIELD> class Fields;
template <typename DATA, typename FIELD> class FieldBuilder;

// ---- fields<FIELD> -------------------------------------------------------
// Private build-time base: per-field sizes/alignments and the offset-packing
// algorithm.  Only used by Overlay and FieldBuilder; never in the public
// inheritance chain seen by dataStore / splitStore.
template <typename FIELD>
class fields {
  protected:
    using offset_type = size_t;
    static constexpr offset_type NO_OFFSET = std::numeric_limits<offset_type>::max();
    static constexpr offset_type log_max_align = datafields::log2(alignof(std::max_align_t));
    static constexpr offset_type field_count = static_cast<offset_type>(FIELD::MAX_FIELD);
    blitz::TinyVector<offset_type, field_count> field_length; // Length of the field
    blitz::TinyVector<offset_type, field_count> field_ialign; // Alignment of the field in log2(align)
  public:
    using field = FIELD;
    fields() : field_length(offset_type(0)), field_ialign(offset_type(0)) {}

    // Enlarge an overlay (master) field to fit the biggest sub-field.
    void enlarge(FIELD f, offset_type size, offset_type alignment) {
        auto ifield = static_cast<offset_type>(f);
        field_length[ifield] = std::max(field_length[ifield], size);
        field_ialign[ifield] = std::max(field_ialign[ifield], alignment);
    }

    void clear() {
        field_length = 0;
        field_ialign = 0;
    }

    template <typename T>
    void add(FIELD f, std::string_view name) {
        static_assert(sizeof(T) > 0);               // Field should have non-zero size
        static_assert(sizeof(T) % alignof(T) == 0); // Field should be properly aligned
        auto ifield = static_cast<offset_type>(f);  // Index of the field
        assert(field_length[ifield] == 0);          // Field should not already exist
        field_length[ifield] = sizeof(T);           // Length of the field
        field_ialign[ifield] = log2(alignof(T));    // Alignment of the field in log2(align)
    }

    template <int N>
    std::tuple<offset_type, offset_type> assign_offsets(blitz::TinyVector<offset_type, N>& field_offset, offset_type start_offset = 0) {
        static_assert(N >= field_count, "field_offset vector must have at least field_count entries");
        // Find the maximum alignment (the required alignment for the element storage)
        auto alignment = blitz::max(field_ialign);
        auto alignment_mask = (1 << alignment) - 1;

        // If the current offset (the start offset) does not match the necessary alignment,
        // then keep padding with a field of the necessary alignment (or nothing) until it does.
        auto offset = start_offset;
        for (auto i = 0; i < alignment; ++i) {
            auto size = 1 << i;
            if ((2 * size - 1) & offset) {
                // Look for a field with the necessary alignment, e.g., a float or 3 floats
                // If we do not find a suitable field, just pad it with nothing
                // If the length of the field is suitable, then the alignment is satisfied
                int j = blitz::first((field_length & size) == size);
                if (j >= 0) {
                    field_offset[j] = offset;
                    offset += field_length[j];
                    field_length[j] = 0; // consume the field (do not assign it below)
                }
                else offset += size;
            }
        }
        assert((offset & alignment_mask) == 0);

        // Calculate how much space each alignment requires.
        blitz::TinyVector<offset_type, log_max_align + 1> align_length(offset_type(0));
        for (auto i = 0; i < field_count; ++i) {
            align_length[field_ialign[i]] += field_length[i];
        }
        // Now calculate the offsets to each alignment block (exclusive scan in reverse order)
        blitz::TinyVector<offset_type, log_max_align + 1> align_offset;
        align_offset[log_max_align] = offset;
        for (auto i = log_max_align; i > 0; --i) {
            align_offset[i - 1] = align_offset[i] + align_length[i];
        }
        offset = align_offset[0] + align_length[0];
        // Assign the offset to the field, ignoring fields that have no data
        for (auto i = 0; i < field_count; ++i) {
            if (field_length[i]) {
                assert(field_offset[i] == NO_OFFSET); // Field should not already exist
                field_offset[i] = align_offset[field_ialign[i]];
                align_offset[field_ialign[i]] += field_length[i];
            }
        }
        // Calculate the size of the element, including padding to the required alignment
        offset = (offset + alignment_mask) & ~alignment_mask;
        return {offset - start_offset, alignment};
    }
};

// ---- Overlay<DATA, FIELD> ------------------------------------------------
// RAII helper returned by FieldBuilder::overlay().  Owns its own fields<FIELD>
// instance to track sizes/alignments of the sub-fields within this overlay.
// On destruction (or explicit commit()) it:
//   1. Packs the sub-fields and stores their relative offsets in the builder.
//   2. Resizes the overlay's master field in the builder to fit the largest variant.
// Method bodies are defined after FieldBuilder because they read its private members.
template <typename DATA, typename FIELD>
class Overlay : protected fields<FIELD> {
    using fields<FIELD>::field_length;
    FieldBuilder<DATA, FIELD>& builder;
    unsigned int master_idx;

  public:
    using field = FIELD;
    Overlay(FieldBuilder<DATA, FIELD>& b, FIELD f) : builder(b), master_idx(static_cast<unsigned int>(f)) {}
    ~Overlay() { commit(); }

    void commit();

    // Add a sub-field; its offset within the overlay is computed at commit().
    template <typename T>
    void add(FIELD f, std::string_view name);

    // Add a sub-field at an explicit byte offset within the overlay.
    template <typename T>
    void add(FIELD f, std::string_view name, size_t offset);
};

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
//! (a) Implicit conversion from a completed FieldBuilder:
//! @code
//!   FieldBuilder<PARTICLE, PKD_FIELD> b;
//!   b.add<double[3]>(PKD_FIELD::oPosition, "r");
//!   b.commit();
//!   particleStore ps;
//!   static_cast<datafields::Fields<PARTICLE, PKD_FIELD>&>(ps) = b;
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

// ---- FieldBuilder<DATA, FIELD> -------------------------------------------
//! \brief Build class — computes the field layout and produces a Fields.
//!
//! Owns all build-time state (per-field sizes, overlay mapping).  It is NOT
//! part of the dataStore / splitStore inheritance chain.  After commit() the
//! result can be obtained as a Fields object:
//! @code
//!   FieldBuilder<PARTICLE, PKD_FIELD> b;
//!   b.set_header<PARTICLE>();          // optional: sets base element size
//!   b.add<double[3]>(PKD_FIELD::oPosition, "r");
//!   b.add<float>(PKD_FIELD::oMass, "m");
//!   b.commit();
//!
//!   // (a) implicit conversion — assign to the Fields base of a store:
//!   static_cast<datafields::Fields<PARTICLE,PKD_FIELD>&>(particleStore) = b;
//!
//!   // (b) MPI / thread distribution — layout is a single flat array:
//!   MPI_Bcast(b.build().get_layout().data(),
//!             (PKD_FIELD::MAX_FIELD + 1) * sizeof(size_t), MPI_BYTE, 0, comm);
//!   // on worker:
//!   worker_store.set_layout(recv_layout);
//! @endcode
template <typename DATA, typename FIELD>
class FieldBuilder : protected fields<FIELD> {
    friend class Overlay<DATA, FIELD>;

    using fields<FIELD>::field_length;
    using offset_type = size_t;
    using overlay_type = size_t;
    static constexpr int field_count = static_cast<int>(FIELD::MAX_FIELD);
    static constexpr size_t NO_OVERLAY = std::numeric_limits<overlay_type>::max();

    blitz::TinyVector<overlay_type, field_count> overlay_field; // Overlay master index, or NO_OVERLAY
    blitz::TinyVector<offset_type, field_count> overlay_offset; // Sub-field offset within overlay

    Fields<DATA, FIELD> result; // Accumulates computed layout; valid after commit()

  public:
    using field = FIELD;

    FieldBuilder() : overlay_field(NO_OVERLAY), overlay_offset(Fields<DATA, FIELD>::NO_OFFSET) {
        result.layout = Fields<DATA, FIELD>::NO_OFFSET;
        result.layout[field_count] = 0; // Caller must call set_header<H>() if there is a fixed header.
    }

    //! Override the base element size (default is sizeof(DATA)).
    //! Must be called before any add() calls when the header type differs from DATA.
    template <typename H>
    void set_header() {
        result.layout[field_count] = sizeof(H);
    }

    //! Register a field of type T.
    template <typename T>
    void add(FIELD f, std::string_view name) {
        fields<FIELD>::template add<T>(f, name);
    }

    //! Begin an overlay group rooted at master field f.
    //! Returns an Overlay RAII object; the overlay is committed when that object
    //! goes out of scope or when overlay.commit() is called explicitly.
    auto overlay(FIELD f) { return Overlay<DATA, FIELD>(*this, f); }

    //! Finalise the layout.  Must be called exactly once, after all add() /
    //! overlay() calls.  Populates the internal Fields result.
    void commit() {
        auto [size, alignment] = this->assign_offsets(result.layout, result.layout[field_count]);
        result.layout[field_count] += size;
        // Adjust sub-field offsets: absolute = master_absolute + relative_within_overlay
        for (auto i = 0; i < field_count; ++i) {
            if (overlay_field[i] != NO_OVERLAY && overlay_field[i] != static_cast<overlay_type>(i)) {
                result.layout[i] = result.layout[overlay_field[i]] + overlay_offset[i];
            }
        }
        // A field that is both the master and a sub-field gets its relative offset added.
        for (auto i = 0; i < field_count; ++i) {
            if (overlay_field[i] == static_cast<overlay_type>(i)) {
                result.layout[i] += overlay_offset[i];
            }
        }
        result.swap_buffer.resize(result.layout[field_count]);
    }

    using layout_type = typename Fields<DATA, FIELD>::layout_type;

    //! Access the completed layout (call after commit()).
    const Fields<DATA, FIELD>& build() const { return result; }

    //! Implicit conversion to layout_type — allows: store.set_layout(builder)
    operator layout_type() const { return result.layout; } // NOLINT: conversion operators cannot use trailing return type
};

// ---- Overlay method bodies (FieldBuilder is now complete) ----------------

template <typename DATA, typename FIELD>
void Overlay<DATA, FIELD>::commit() {
    // Pack sub-fields starting at offset 0 within the overlay.
    auto [size, alignment] = this->assign_offsets(builder.overlay_offset);
    // Ensure the master field in the builder is large enough for this overlay.
    builder.fields<FIELD>::enlarge(static_cast<FIELD>(master_idx), size, alignment);
    // Zero out so a redundant commit() (e.g., from the destructor) is a no-op.
    this->clear();
}

template <typename DATA, typename FIELD>
template <typename T>
void Overlay<DATA, FIELD>::add(FIELD f, std::string_view name) {
    fields<FIELD>::template add<T>(f, name);
    builder.overlay_field[static_cast<int>(f)] = master_idx;
}

template <typename DATA, typename FIELD>
template <typename T>
void Overlay<DATA, FIELD>::add(FIELD f, std::string_view name, size_t offset) {
    builder.overlay_field[static_cast<int>(f)] = master_idx;
    builder.overlay_offset[static_cast<int>(f)] = offset;
}

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

// Concept that accepts a scalar type, a TinyVector, or a type with a tuple schema (e.g., a struct defined with MAKE_STRUCTURE)
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
