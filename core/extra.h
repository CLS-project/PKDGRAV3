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

#ifndef CORE_extra_H
#define CORE_extra_H
#include "pkd_config.h"
#include <cstdint>
#include <limits>
#include <sys/types.h>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <iterator>
#include "datastore.h"
#include "element.h"
#include "meshless.h"

// This should be in a different header file
enum class EXTRA_FIELD {
    oMFM, /* Sph structure */
    oMFV, /* MFV structure */
    oStar, /* Star structure */
    oBH, /* BH structure */

    MAX_FIELD
};

class extraStore : public dataStore<meshless::EXTRAFIELDS, EXTRA_FIELD> {
public:
    void initialize() {
    }
    auto ExtraSize() const {return ElementSize(); }
    void Add(int i) { SetLocal(Local()+i); }
    uint64_t append() {
        SetLocal(Local()+1);
        return Local()-1; //NB: we return the index, not the number
    }
    // auto &MFM( PARTICLE *p ) const {
    //     return get<meshless::FIELDS>(p,EXTRA_FIELD::oMFM);
    // }
    // const auto &MFM( const PARTICLE *p ) const {
    //     return get<meshless::FIELDS>(p,EXTRA_FIELD::oMFM);
    // }
    // auto &MFV( PARTICLE *p ) const {
    //     return get<meshless::MFV>(p,EXTRA_FIELD::oMFV);
    // }
    // const auto &MFV( const PARTICLE *p ) const {
    //     return get<meshless::MFV>(p,EXTRA_FIELD::oMFV);
    // }
    // auto &star( PARTICLE *p ) const {
    //     return get<meshless::STAR>(p,EXTRA_FIELD::oStar);
    // }
    // const auto &star( const PARTICLE *p ) const {
    //     return get<meshless::STAR>(p,EXTRA_FIELD::oStar);
    // }
    // auto &bh( PARTICLE *p ) const {
    //     return get<meshless::STAR>(p,EXTRA_FIELD::oBH);
    // }
    // const auto &bh( const PARTICLE *p ) const {
    //     return get<meshless::STAR>(p,EXTRA_FIELD::oBH);
    // }

    class ExtraField : public SingleElement<meshless::EXTRAFIELDS,extraStore> {
    public:
        using SingleElement::SingleElement;
    public:
    };

    /// @brief A "pointer" reference to a "extra field"
    ///
    /// This uses the default copy/assignment semantics (the internal pointers are copied).
    /// The Particle class uses direct assignment.
    using ExtraFieldPointer = SinglePointer<ExtraField>;

    /// @brief A "reference" to a extra field
    /// This allows extra fields to be assigned and swapped
    using ExtraFieldReference = SingleReference<ExtraField>;

    /// @brief Return a pointer proxy to a extra fields
    /// @param i Index of the extra fields
    /// @return pointer proxy
    auto operator[](int i) {
        return ExtraFieldReference(*this,i);
    }
    auto operator[](meshless::EXTRAFIELDS *e) {
        return ExtraFieldReference(*this,e);
    }
    // auto operator[](SingleReference<Particle> p) {
    //     auto i = p.extra_index();
    //     return ExtraFieldReference(*this,i);
    // }
    // const auto operator[](const PARTICLE *p) const {
    //     auto i = p.extra_index();
    //     return ExtraFieldReference(*this,i);
    // }

};


#endif //CORE_extra_H

