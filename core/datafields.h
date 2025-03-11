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
#include <vector>
#include <algorithm>
#include <string_view>
#include "blitz/array.h"

namespace datafields {
static constexpr size_t log2(size_t n) {
    return (n <= 1) ? 0 : 1 + log2(n / 2);
}

template<typename FIELD>
class fields {
protected:
    using offset_type = size_t;
    static constexpr offset_type log_max_align = datafields::log2(alignof(std::max_align_t));
    static constexpr offset_type field_count = static_cast<offset_type>(FIELD::MAX_FIELD);
    blitz::TinyVector<offset_type,field_count> field_length; // Length of the field
    blitz::TinyVector<offset_type,field_count> field_ialign; // Alignment of the field in log2(align)
public:
    using field = FIELD;
    fields() : field_length(offset_type(0)), field_ialign(offset_type(0)) {}

    // This will enlarge an overlay field if necessary
    void enlarge(FIELD f, offset_type size, offset_type alignment) {
        auto ifield = static_cast<offset_type>(f);
        field_length[ifield] = std::max(field_length[ifield],size);
        field_ialign[ifield] = std::max(field_ialign[ifield],alignment);
    }

    void clear() {
        field_length = 0;
        field_ialign = 0;
    }

    template<typename T>
    void add(FIELD f,std::string_view name) {
        static_assert(sizeof(T) % alignof(T) == 0); // Field should be properly aligned
        auto ifield = static_cast<offset_type>(f);      // Index of the field
        assert(field_length[ifield] == 0);          // Field should not already exist
        field_length[ifield] = sizeof(T);           // Length of the field
        field_ialign[ifield] = log2(alignof(T));    // Alignment of the field in log2(align)
    }

    std::tuple<offset_type,offset_type>
    assign_offsets(blitz::TinyVector<offset_type,field_count> &field_offset, offset_type start_offset=0) {
        // Find the maximum alignment (the required alignment for the element storage)
        auto alignment = blitz::max(field_ialign);
        auto alignment_mask = (1 << alignment) - 1;

        // If the current offset (the start offset) does not match the necessary alignment,
        // then keep padding with a field of the necessary alignment (or nothing) until it does.
        auto offset = start_offset;
        for (auto i=0; i<alignment; ++i) {
            auto size = 1 << i;
            if ( (2*size-1) & offset) {
                // Look for a field with the necessary alignment, e.g., a float or 3 floats
                // If we do not find a suitable field, just pad it with nothing
                // If the length of the field is suitable, then the alignment is satisfied
                int j = blitz::first((field_length&size) == size);
                if (j>=0) {
                    field_offset[j] = offset;
                    offset += field_length[j];
                    field_length[j] = 0;
                }
                else offset += size;
            }
        }
        assert((offset & alignment_mask)==0);

        // Calculate how much space each alignment requires.
        blitz::TinyVector<offset_type,log_max_align+1> align_length(offset_type(0));
        for (auto i=0; i<field_count; ++i) {
            align_length[field_ialign[i]] += field_length[i];
        }
        // Now calculate the offsets to each alignment block (exclusive scan in reverse order)
        blitz::TinyVector<offset_type,log_max_align+1> align_offset;
        align_offset[log_max_align] = offset;
        for (auto i=log_max_align; i>0; --i) {
            align_offset[i-1] = align_offset[i] + align_length[i];
        }
        offset = align_offset[0] + align_length[0];
        // Assign the offset to the field, ignoring fields that have no data
        for (auto i=0; i<field_count; ++i) {
            if (field_length[i]) {
                assert(field_offset[i] == 0); // Field should not already exist
                field_offset[i] = align_offset[field_ialign[i]];
                align_offset[field_ialign[i]] += field_length[i];
            }
        }
        // Calculate the size of the element, including padding to the required alignment
        offset = (offset + alignment_mask) & ~alignment_mask;
        return {offset - start_offset,alignment};
    }
};

template<typename DATA,typename FIELD>
class Fields;

template<typename DATA,typename FIELD>
class Overlay : protected fields<FIELD> {
    Fields<DATA,FIELD> &dfields;
    unsigned int overlay_field;

public:
    using field = FIELD;
    Overlay(Fields<DATA,FIELD> &dfields,FIELD f) : dfields(dfields), overlay_field(static_cast<unsigned int>(f)) {
    }
    ~Overlay() {
        commit();
    }

    void commit() {
        // This will assign offsets to the fields in the overlay starting at zero.
        // We have already recorded the overlay field so that we can later adjust the offset.
        auto [size,alignment] = this->assign_offsets(dfields.overlay_offset);
        // Make sure that the overlay field is large enough for the largest overlay
        dfields.enlarge(static_cast<FIELD>(overlay_field),size,alignment);
        // We set the lengths to zero so that the overlay fields are not processed again
        // We can "commit" the overlay, or it can be committed in the destructor (or both)
        this->clear();
    }

    template<typename T>
    void add(FIELD f,std::string_view name) {
        auto ifield = static_cast<int>(f);
        fields<FIELD>::template add<T>(f,name);
        dfields.overlay_field[ifield] = overlay_field;
    }
    template<typename T>
    void add(FIELD f,std::string_view name, size_t offset) {
        auto ifield = static_cast<int>(f);
        dfields.overlay_field[ifield] = overlay_field;
        dfields.overlay_offset[ifield] = offset;
    }
};

//! \brief Track offset of optional data fields
//!
//! The size of particles and tree nodes can vary depending on what
//! features are enabled at run-time. This class tracks which fields
//! are present, and their offset in the structure.
template<typename DATA,typename FIELD>
class Fields : protected fields<FIELD> {
    friend class Overlay<DATA,FIELD>;
public:
    using offset_type = size_t;
    using overlay_type = size_t;
    static constexpr auto field_count = static_cast<int>(FIELD::MAX_FIELD);
    using field = FIELD;
    static constexpr size_t NO_OVERLAY = std::numeric_limits<overlay_type>::max();
protected:
    size_t max_align = log2(alignof(DATA));
    size_t element_size = sizeof(DATA);
    blitz::TinyVector<offset_type, field_count> field_offset;  // Offset to the field within the specified alignment
    blitz::TinyVector<overlay_type,field_count> overlay_field; // Overlay index or NO_OVERLAY
    blitz::TinyVector<offset_type, field_count> overlay_offset;// Offset within the overlay

    // blitz::TinyVector<size_t,field_count> oFieldOffset;  //!< The offset (in bytes) to each field
    // size_t iElementSize = sizeof(DATA);      //!< Size (in bytes) of a single element
    // size_t nElementAlign = alignof(DATA);    //!< Requirement alignment of the element (usually 4 or 8)
    // size_t iElement32 = 0;
public:
    Fields() :
        field_offset(offset_type(0)),
        overlay_field(NO_OVERLAY),
        overlay_offset(offset_type(0))
        // oFieldOffset(size_t(0))
    {}
    //! Setup the field offsets for adding fields.
    //! Here we set the basic size of the element to a fixed value (the header size).
    //! Following this call, "add" is called for each element that should be present.
    //! \param iBasicSize The size of the element header (the fixed part)
    // void initialize(size_t iBasicSize=0) {
    //     iElementSize = nElementAlign = iBasicSize;
    //     iElement32 = 0;
    //     // oFieldOffset = 0;
    // }
    template<typename H>
    void set_header() {
        max_align = log2(alignof(H));
        element_size = sizeof(H);
    }

    //! Make sure that elements are properly aligned.
    //! This is called after all fields have been added. The size is adjusted to respect alignment.
    void commit(void) {
        // iElementSize = (iElementSize + nElementAlign - 1 ) & ~(nElementAlign-1);

        auto [size,alignment] = this->assign_offsets(field_offset,element_size);
        element_size += size;
        // Adjust the offsets for overlay fields (but skip the "master" field)
        for (auto i=0; i<field_count; ++i) {
            if (overlay_field[i] != NO_OVERLAY && overlay_field[i] != i) {
                field_offset[i] = field_offset[overlay_field[i]] + overlay_offset[i];
            }
        }
        // Now correct the "master" field
        for (auto i=0; i<field_count; ++i) {
            if (overlay_field[i] == i) {
                field_offset[i] += overlay_offset[i];
            }
        }
        this->clear();
    }

    //! Add a field to the store. A field of type T is added with the specified offset.
    //! \param f The field id
    //! \param offset The offset to the field
    // template<typename T,std::enable_if_t<!std::is_void_v<T>,bool> = true>
    // void add(FIELD f,std::string_view name, size_t offset) {
    //     oFieldOffset[static_cast<unsigned int>(f)] = offset;
    // }

    //! Add a field to the store. A field of type T is added and the offset recorded. The size of the
    //! element is increased to account for this. For types with multiple values use array syntax.
    //! @code
    //!   particles.add<double[3]>(PKD_FIELD::oPosition)
    //!   particles.add<float>(PKD_FIELD::oMass)
    //! @endcode
    //! Overlaying fields is possible by specifying a type with the maximum size and alignment
    //! for S. Normally this is done with a union. For example, to overlay a double[3] and a float[4]
    //! @code
    //!   union u {
    //!       double d[3];
    //!       float f[4];
    //!   };
    //!   particles.add<double[3],u>(PKD_FIELD::oPosition,"r")
    //!   particles.add<float[4]>(PKD_FIELD::oWhatever,"q",particles.offset(PKD_FIELD::oPosition))
    //! @endcode
    //! \param f The field id
    template<typename T>
    void add(FIELD f,std::string_view name) {
        fields<FIELD>::template add<T>(f,name);
    }

    //! Overlay a field on top of another field. This is used to create a union of fields.
    //! \param f The field id
    auto overlay(FIELD f) {
        return Overlay<DATA,FIELD>(*this,f);
    }

public:
    //! Returns the offset (in bytes) of the specified field or zero if it not present
    //! \param f The field id
    auto offset(FIELD f) const noexcept {return field_offset[static_cast<unsigned int>(f)];}
    //! Returns true if the field is present, otherwise false.
    //! \param f The field id
    auto present(FIELD f) const noexcept {return field_offset[static_cast<unsigned int>(f)]!= 0;}
    //! Returns a pointer to the field if present or nullptr if not.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<!std::is_array_v<T>,bool> = true>
    const auto & get(const DATA *p,FIELD f) const noexcept {
        auto v = reinterpret_cast<const char *>(p);
        return * reinterpret_cast<const T *>(v + field_offset[static_cast<unsigned int>(f)]);
    }
    //! Returns a pointer to the field if present or nullptr if not.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<!std::is_array_v<T>,bool> = true>
    auto & get(DATA *p,FIELD f) const noexcept {
        auto v = reinterpret_cast<char *>(p);
        return * reinterpret_cast<T *>(v + field_offset[static_cast<unsigned int>(f)]);
    }

    //! Returns the field as a blitz++ TinyVector.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<std::rank_v<T> == 1,bool> = true>
    const auto &get(const DATA *p,FIELD f) const noexcept {
        return get<blitz::TinyVector<typename std::remove_extent<T>::type,std::extent_v<T>>>(p,f);
    }

    //! Returns the field as a blitz++ TinyVector.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<std::rank_v<T> == 1,bool> = true>
    auto &get(DATA *p,FIELD f) const noexcept {
        return get<blitz::TinyVector<typename std::remove_extent<T>::type,std::extent_v<T>>>(p,f);
    }

    //! Returns the size (in bytes) of an element (all fields).
    auto ElementSize() const noexcept {return element_size; }
    //! Returns a pointer to the i'th element at pBase.
    //! \param pBase Pointer to the first element
    //! \param i Index of the desired element
    DATA *Element(void *pBase, int i) const noexcept {
        auto v = static_cast<char *>(pBase);
        return reinterpret_cast<DATA *>(v + ((uint64_t)i)*ElementSize());
    }
};
} // namespace datafields
#endif
