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

#ifndef CORE_DATABUILDER_H
#define CORE_DATABUILDER_H

#include "datafields.h"
#include <algorithm> // std::max
#include <cassert>

namespace datafields {

static constexpr size_t log2(size_t n) {
    return (n <= 1) ? 0 : 1 + log2(n / 2);
}

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
        static_assert(N >= static_cast<int>(field_count), "field_offset vector must have at least field_count entries");
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
        for (auto i = 0; i < static_cast<int>(field_count); ++i) {
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
        for (auto i = 0; i < static_cast<int>(field_count); ++i) {
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
template <typename DATA, typename FIELD> class FieldBuilder; // forward

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

// ---- FieldBuilder<DATA, FIELD> -------------------------------------------
//! \brief Build class — computes the field layout and produces a Fields.
//!
//! Owns all build-time state (per-field sizes, overlay mapping).  It is NOT
//! part of the dataStore / splitStore inheritance chain.  After commit() the
//! result can be obtained via implicit conversion to layout_type, allowing
//! direct use with store.set_layout():
//! @code
//!   FieldBuilder<PARTICLE, PKD_FIELD> b;
//!   b.set_header<PARTICLE>();          // optional: sets base element size
//!   b.add<double[3]>(PKD_FIELD::oPosition, "r");
//!   b.add<float>(PKD_FIELD::oMass, "m");
//!   b.commit();
//!
//!   particleStore.set_layout(b);        // implicit conversion to layout_type
//!
//!   // MPI / thread distribution — layout is a single flat array:
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
    using layout_type = typename Fields<DATA, FIELD>::layout_type;

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

} // namespace datafields

#endif
