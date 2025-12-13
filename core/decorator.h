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

#ifndef CORE_DECORATOR_H
#define CORE_DECORATOR_H
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
#include "decorator_fields.h"
#include "boost/dynamic_bitset.hpp"

// This should be in a different header file
enum class DECORATOR_FIELD {
    oMFM,  /* Sph structure */
    oMFV,  /* MFV structure */
    oStar, /* Star structure */
    oBH,   /* BH structure */

    MAX_FIELD
};

class decoratorStore : public dataStore<decorators::FIELDS, DECORATOR_FIELD> {
  public:
    auto DecoratorSize() const { return ElementSize(); }
    // Given a pointer, return the index in the store
    uint64_t DecoratorIndex(value_type* d) const { return (reinterpret_cast<char*>(d) - reinterpret_cast<char* const>(pStore)) / ElementSize(); }
    uint64_t DecoratorIndex(const value_type* d) const { return (reinterpret_cast<const char*>(d) - reinterpret_cast<char* const>(pStore)) / ElementSize(); }
    void Add(int i) { append(i); }
    int SetLocal(int n) {
        resize(n);
        return n;
    }

    const auto& mfm(const decorators::FIELDS* d) const { return get<decorators::MFM>(d, DECORATOR_FIELD::oMFM); }
    const auto& mfv(const decorators::FIELDS* d) const { return get<decorators::MFV>(d, DECORATOR_FIELD::oMFV); }
    const auto& star(const decorators::FIELDS* d) const { return get<decorators::STAR>(d, DECORATOR_FIELD::oStar); }
    const auto& bh(const decorators::FIELDS* d) const { return get<decorators::BLACKHOLE>(d, DECORATOR_FIELD::oBH); }

    class DecoratorField : public SingleElement<decorators::FIELDS, decoratorStore> {
      public:
        using SingleElement::SingleElement;

      public:
        auto& mfm() const { return store().mfm(p); }
        auto& mfv() const { return store().mfv(p); }
        auto& star() const { return store().star(p); }
        auto& bh() const { return store().bh(p); }
    };

    /// @brief A "pointer" reference to a decorator
    ///
    /// This uses the default copy/assignment semantics (the internal pointers are copied).
    /// The Particle class uses direct assignment.
    using DecoratorFieldPointer = SinglePointer<DecoratorField>;

    /// @brief A "reference" to a decorator
    /// This allows decorator fields to be assigned and swapped
    using DecoratorFieldReference = SingleReference<DecoratorField>;

    /// @brief Return a pointer proxy to a decorator fields
    /// @param i Index of the decorator fields
    /// @return pointer proxy
    auto operator[](int i) { return DecoratorFieldReference(*this, i); }
    auto operator[](decorators::FIELDS* e) { return DecoratorFieldReference(*this, e); }
};

#endif // CORE_DECORATOR_H
