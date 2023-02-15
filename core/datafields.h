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

#include <vector>
#include <algorithm>
#include "blitz/array.h"

//! \brief Track offset of optional data fields
//!
//! The size of particles and tree nodes can vary depending on what
//! features are enabled at run-time. This class tracks which fields
//! are present, and their offset in the structure.
template<typename DATA,typename FIELD>
class dataFields {
public:
    using field = FIELD;
protected:
    std::vector<int> oFieldOffset;  //!< The offset (in bytes) to each field of the element
    size_t iElementSize = sizeof(DATA);     //!< Size (in bytes) of a single element
    size_t nElementAlign = alignof(DATA);   //!< Requirement alignment of the element (usually 4 or 8)
    size_t iElement32 = 0;
public:
    dataFields() : oFieldOffset(static_cast<size_t>(FIELD::MAX_FIELD),0) {}
    //! Setup the field offsets for adding fields.
    //! Here we set the basic size of the element to a fixed value (the header size).
    //! Following this call, "add" is called for each element that should be present.
    //! \param iBasicSize The size of the element header (the fixed part)
    void initialize(int iBasicSize=0) {
        iElementSize = nElementAlign = iBasicSize;
        iElement32 = 0;
        oFieldOffset.clear();
        oFieldOffset.insert(oFieldOffset.end(),static_cast<size_t>(FIELD::MAX_FIELD),0);
    }
    //! Make sure that elements are properly aligned.
    //! This is called after all fields have been added. The size is adjusted to respect alignment.
    void align(void) {
        iElementSize = (iElementSize + nElementAlign - 1 ) & ~(nElementAlign-1);
    }
    //! Add a "void" type field. For some element the offset isn't relevant, only if it is present.
    //! \param f The field id
    template<typename T,std::enable_if_t<std::is_void_v<T>,bool> = true>
    void add(FIELD f) {
        oFieldOffset[static_cast<unsigned int>(f)] = 1;
    }
    //! Add a field to the store. A field of type T is added and the offset recorded. The size of the
    //! element is increased to account for this. For types with multiple values use array syntax.
    //! @code
    //!   particles.add<double[3]>(PKD_FIELD::oPosition)
    //!   particles.add<float>(PKD_FIELD::oMass)
    //! @endcode
    //! \param f The field id
    template<typename T,std::enable_if_t<!std::is_void_v<T>,bool> = true>
    void add(FIELD f) {
        static_assert(std::is_standard_layout<T>());
        static_assert(std::is_void_v<T> || alignof(T) <= 4 || alignof(T) == 8);
        int iOffset = iElementSize;
        int iAlign = std::max(4ul,alignof(T));
        if (iAlign==4 && iElement32 && (sizeof(T)%8)) {
            iOffset = iElement32;
            iElement32 = 0;
            auto iExtra = sizeof(T) - 4;
            assert(iExtra%8 == 0);
            for (auto &i : oFieldOffset) {
                if (i>iOffset) i += iExtra;
            }
            iElementSize += iExtra;
        }
        else {
            auto iMask = iOffset & (iAlign-1);
            if (iMask) {
                iElement32 = iOffset;
                iOffset += iAlign - iMask;
            }
            assert((iOffset & (iAlign-1)) == 0);
            iElementSize = iOffset + sizeof(T);
        }
        if (nElementAlign < iAlign) nElementAlign = iAlign;
        oFieldOffset[static_cast<unsigned int>(f)] = iOffset;
    }
public:
    //! Returns true if the field is present, otherwise false.
    //! \param f The field id
    auto present(FIELD f) const noexcept {return oFieldOffset[static_cast<unsigned int>(f)]!= 0;}
    //! Returns a pointer to the field if present or nullptr if not.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<!std::is_array_v<T>,bool> = true>
    const auto & get(const DATA *p,FIELD f) const noexcept {
        auto v = reinterpret_cast<const char *>(p);
        return * reinterpret_cast<const T *>(v + oFieldOffset[static_cast<unsigned int>(f)]);
    }
    //! Returns a pointer to the field if present or nullptr if not.
    //! \param p Pointer to an element
    //! \param f The field id
    template<typename T,std::enable_if_t<!std::is_array_v<T>,bool> = true>
    auto & get(DATA *p,FIELD f) const noexcept {
        auto v = reinterpret_cast<char *>(p);
        return * reinterpret_cast<T *>(v + oFieldOffset[static_cast<unsigned int>(f)]);
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
    auto ElementSize() const noexcept {return iElementSize; }
    //! Returns a pointer to the i'th element at pBase.
    //! \param pBase Pointer to the first element
    //! \param i Index of the desired element
    DATA *Element(void *pBase, int i) const noexcept {
        auto v = static_cast<char *>(pBase);
        return reinterpret_cast<DATA *>(v + ((uint64_t)i)*ElementSize());
    }
};
#endif
