/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2020 Douglas Potter & Joachim Stadel
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
#include "arc.h"
#include <cassert>
#include <iostream>
#include <algorithm>

namespace mdl {

using std::uint32_t;
using std::uint64_t;
using std::size_t;

namespace murmur {
template<int len>
uint64_t murmur(const uint64_t *key) {
    const uint64_t m = 0xc6a4a7935bd1e995;
    const int r = 47;
    uint64_t h = 0xdeadbeefdeadbeef;
    for(auto i=0; i<len; ++i) {
	uint64_t k = *key++;
	k *= m;
	k ^= k >> r;
	k *= m;
	h ^= k;
	h *= m;
	}
    h ^= h >> r;
    h *= m;
    h ^= h >> r;
    return h;
    }

template<int len>
uint32_t murmur2(const uint32_t *key) {
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    uint32_t h = 0xdeadbeef /*^ len : len will be the same */;
    for(auto i=0; i<len; ++i) {
	uint32_t k = *key++;
	k *= m;
	k ^= k >> r;
	k *= m;
	h *= m;
	h ^= k;
	}
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
    return h;
    }

/*
** This makes the following assumptions (changed from regular hash)
**   1. keys lengths are a multiple of four bytes (and at least four bytes)
**   2. length is always identical (so don't need to mix in the length)
**   3. We will always use the same seed
*/
template<int len>
uint32_t murmur3(const uint32_t* key) {
    uint32_t h = 0xdeadbeef;
    for(auto i=0; i<len; ++i) {
	uint32_t k = *key++;
	k *= 0xcc9e2d51;
	k = (k << 15) | (k >> 17);
	k *= 0x1b873593;
	h ^= k;
	h = (h << 13) | (h >> 19);
	h = h * 5 + 0xe6546b64;
	}
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
    }

} // namespace murmur

/*****************************************************************************\
* HASH table
\*****************************************************************************/

namespace hash {

static uint32_t hash(uint32_t uLine,uint32_t uId) {
    uint32_t key[] = {uLine,uId};
    return murmur::murmur2<2>(key);
    }


} // namespace hash

/*****************************************************************************\
* ARC cache
\*****************************************************************************/

template<>
void *ARC<>::fetch(uint32_t uIndex, int uId, bool bLock,bool bModify,bool bVirtual) {
    auto uLine = uIndex >> nLineBits;
    auto iInLine = uIndex & nLineMask;
    uint32_t rat;

    /* First check our own cache */
    auto uHash = hash::hash(uLine,uId);
    CDBL::iterator item;
    auto &Hash = HashChains[uHash&uHashMask];
    auto pEntry = find_key(Hash,uLine,uId);
    if (pEntry) {   // Page is reference by the cache
	item = CDBL::s_iterator_to(*pEntry);
	update(item,bLock,bModify);
	auto &cdb = std::get<0>(*item);
	if (cdb.data == nullptr) {
	    if (cdb.absent()) return nullptr;
	    helper->invokeRequest(uLine,uId,nullptr,bVirtual); // Request the element be fetched
	    cdb.data = replace(cdb.where()==B2);
	    cdb.uId = uId;
	    cdb.uPage = uLine;
	    helper->finishRequest(uLine,uId,nullptr,bVirtual,cdb.data);
	    }
	}
    else {          // Page is not in the cache
	helper->invokeRequest(uLine,uId,nullptr,bVirtual);
	item = insert(Hash);
	auto &cdb = std::get<0>(*item);
	cdb.uId = uId;
	cdb.uPage = uLine;
	helper->finishRequest(uLine,uId,nullptr,bVirtual,cdb.data);
	}
    auto &cdb = std::get<0>(*item);
    if (bLock) ++cdb.data[-1];   // Increase the lock count if requested
    if (bModify) cdb.setDirty(); // Mark page as dirty

    return reinterpret_cast<char*>(cdb.data) + uDataSizeInBytes*iInLine;
    }

} // namespace mdl
