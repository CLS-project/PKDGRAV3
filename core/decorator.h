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
    oMFM, /* Sph structure */
    oMFV, /* MFV structure */
    oStar, /* Star structure */
    oBH, /* BH structure */

    MAX_FIELD
};

class decoratorStore : public dataStore<decorators::FIELDS, DECORATOR_FIELD> {
private:
    int nCache = 0;     //!< Number of occupied elements in the cache
    bool readonly;      //!< Whether the cache is used as part of a RO particle cache
    boost::dynamic_bitset<> cacheUsed;
public:
    void initialize() {
        cacheUsed = boost::dynamic_bitset();
        cacheUsed.reserve(nStore/2);
    }
    auto DecoratorSize() const {return ElementSize();}
    // Given a pointer, return the index in the store
    uint64_t DecoratorIndex(value_type *d) const {
        return (reinterpret_cast<char *>(d)-reinterpret_cast<char *const>(pStore))/ElementSize();
    }
    uint64_t DecoratorIndex(const value_type *d) const {
        return (reinterpret_cast<const char *>(d)-reinterpret_cast<char *const>(pStore))/ElementSize();
    }
    void Add(int i) { SetLocal(Local()+i); }
    int SetLocal(int n) {
        nLocal = n;
        cacheUsed.resize(nStore-nLocal);
        return nLocal;
    }
    uint64_t append() {
        SetLocal(Local()+1);
        return Local()-1; //NB: we return the index, not the number
    }
    //! Prepare the cache depending on whether it will act as read-only or not
    void CacheStart(bool readonly) {
        CacheClean();
        nCache = 0;
        readonly = readonly;
    }
    //! Clean the cache
    void CacheClean() {
        nCache = 0;
        cacheUsed.reset();
    }
    //! Return the current total size of the cache
    auto CacheSize() {
        return nStore - nLocal;
    }
    //! Return the current number of available cache elements
    auto CacheAvail() {
        return CacheSize() - nCache;
    }
    //! Return a element of the cache that can be used to
    //! temporarily store a decorator
    value_type *CacheAdd() {
        //CacheDebug();
        if (readonly) {
            // We just take the next element in line
            nCache++;
            if (!CacheAvail())
                nCache = 1;
            return Element(nLocal+nCache-1);
        }
        else {
            // Try as if the cache was read-only, but if we find that that element
            // is occupied, we start looping forward.
            auto size = cacheUsed.size();
            auto i = (nCache==size)?0:nCache;
            //assert(cacheUsed.count()!=cacheUsed.size()); // Decorator cache is full
            while (cacheUsed.test(i)) {
                i++;
                if (i==size) i=0;
            }
            nCache++;
            cacheUsed.flip(i);
            return Element(nLocal+i);
        }
    }
    void CacheDebug() {
        //printf("\n\tCache size: %d\n\tCache usage: %d(%ld)\n",
        //       CacheSize(), nCache, cacheUsed.count());

        assert(nCache == cacheUsed.count());

    }
    const void CacheRelease(const value_type *e) {
        //CacheDebug();
        auto i = DecoratorIndex(e) - nLocal;
        // if (!cacheUsed.test(i-nLocal))
        //     printf("Releasing %d \t %e \n",i-nLocal, e->sph.Uint);
        //assert(cacheUsed[i]);
        if (cacheUsed.test(i)) {
            cacheUsed.flip(i);
            nCache--;
        }
        else {
            // This seems to happen when flushing..
            // assert();
        }
    }


    const auto &mfm( const decorators::FIELDS *d ) const {
        return get<decorators::MFM>(
                   d, DECORATOR_FIELD::oMFM
               );
    }
    const auto &mfv( const decorators::FIELDS *d ) const {
        return get<decorators::MFV>(
                   d, DECORATOR_FIELD::oMFV
               );
    }
    const auto &star( const decorators::FIELDS *d ) const {
        return get<decorators::STAR>(
                   d, DECORATOR_FIELD::oStar
               );
    }
    const auto &bh( const decorators::FIELDS *d ) const {
        return get<decorators::BLACKHOLE>(
                   d, DECORATOR_FIELD::oBH
               );
    }

    class DecoratorField : public SingleElement<decorators::FIELDS,decoratorStore> {
    public:
        using SingleElement::SingleElement;
    public:
        auto &mfm()         const { return store().mfm(p); }
        auto &mfv()         const { return store().mfv(p); }
        auto &star()         const { return store().star(p); }
        auto &bh()         const { return store().bh(p); }
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
    auto operator[](int i) {
        return DecoratorFieldReference(*this,i);
    }
    auto operator[](decorators::FIELDS *e) {
        return DecoratorFieldReference(*this,e);
    }

};


#endif //CORE_DECORATOR_H

