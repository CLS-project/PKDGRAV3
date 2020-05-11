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
//    while(len--) {
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

static uint32_t swar32(uint32_t x) {
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    return(x);
    }
} // namespace murmur

/*****************************************************************************\
* HASH table
\*****************************************************************************/

namespace hash {
template<>
uint32_t HASH<>::hash(uint32_t uLine,uint32_t uId) {
    uint32_t key[] = {uLine,uId,uLine,uId,uLine,uId,uLine,uId};
    return murmur::murmur2<2>(key) & uHashMask;
    }

// Remove the CDB entry (given by a list iterator) from the hash table
template<>
HASH<>::CDBL::iterator HASH<>::remove(const CDBL::iterator p) {
    auto &cdb = std::get<0>(*p);
    uint32_t uHash = hash(cdb.uPage,cdb.getId());
    auto & Hash = HashChains[uHash];
    assert(!Hash.empty());
    // Move the matching hash table entry to the start of the collision chain, then move it to the free list
    std::partition(Hash.begin(), Hash.end(), [&](const CDBL::iterator & i) { return i == p; });
    assert(Hash.front() == p);
    HashFree.splice_after(HashFree.before_begin(),Hash,Hash.before_begin());
    return p;
    }

// Typically we want to remove the LRU (front) element from a list, in anticipation of reusing the entry
template<>
HASH<>::CDBL::iterator HASH<>::remove(uint32_t iTarget) {
    return remove(L[iTarget].begin());
    }

/*****************************************************************************\
* Generic interface: Deals with elements that have a key
\*****************************************************************************/
template<>
void HASH<>::drop_all() {
    HASH::clear();
    }

template<>
void * HASH<>::lookup(uint32_t uHash, const void *vKey) {
    auto key = reinterpret_cast<const uint32_t*>(vKey);
    auto uLine = key[0];
    auto uId = key[1];
    auto &Hash = HashChains[uHash&uHashMask];
    auto iskey = [uLine,uId](CDBL::iterator &i) {return std::get<0>(*i).uPage == uLine && (std::get<0>(*i).getId()) == uId;};
    auto match = std::find_if(Hash.begin(),Hash.end(),iskey);
    if (match != Hash.end()) return std::get<0>(**match).data;
    else return nullptr;
    }

template<>
void HASH<>::insert(uint32_t uHash, const void *vKey, void *data) {
    auto key = reinterpret_cast<const uint32_t*>(vKey);
    auto uLine = key[0];
    auto uId = key[1];
    auto &Hash = HashChains[uHash&uHashMask];
    auto item = move(0);
    auto &cdb = std::get<0>(*item);
    cdb.uId = uId;
    cdb.uPage = uLine;
    cdb.data = reinterpret_cast<uint64_t*>(data);
    assert(!HashFree.empty());
    Hash.splice_after(Hash.before_begin(),HashFree,HashFree.before_begin());
    assert(!Hash.empty());
    Hash.front() = item;
    }

template<>
void HASH<>::remove(uint32_t uHash, const void *vKey) {
    auto key = reinterpret_cast<const uint32_t*>(vKey);
    auto uLine = key[0];
    auto uId = key[1];
    auto &Hash = HashChains[uHash&uHashMask];
    auto iskey = [uLine,uId](CDBL::iterator &i) {return std::get<0>(*i).uPage == uLine && (std::get<0>(*i).getId()) == uId;};
    auto match = std::find_if(Hash.begin(),Hash.end(),iskey);
    if (match != Hash.end()) {
	// Move the matching hash table entry to the start of the collision chain, then move it to the free list
	std::partition(Hash.begin(), Hash.end(), [match](const CDBL::iterator & i) { return i == *match; });
	HashFree.splice_after(HashFree.before_begin(),Hash,Hash.before_begin());
	free(*match);
	}
    }

/*****************************************************************************\
* Factory method: Contruct a suitable ARC cache for this HASH table.
\*****************************************************************************/

template<>
ARC * HASH<>::clone(int nMaxElements) {
    return nullptr;
    }

} // namespace hash

/*****************************************************************************\
* ARC cache
\*****************************************************************************/

template<>
void ARC<>::initialize(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
    this->uLineSizeInWords = (uLineSizeInBytes+7) >> 3;       // Size of a cache line (aligned properly)
    // Calculate nCache based on the number of cache lines that will fit in out cache buffer
    // Account for the "magic" number before the cache line.
    auto nCache = (uCacheSizeInBytes>>3) / (this->uLineSizeInWords+1);
    if (this->nCache == nCache && this->uLineSizeInBytes == uLineSizeInBytes) return; // Already setup
    this->uLineSizeInBytes = uLineSizeInBytes;
    this->uDataSizeInBytes = (uLineSizeInBytes) >> nLineBits; // Size of a single element
    this->nLineBits = nLineBits;
    this->nLineMask = (1 << nLineBits) - 1;
    this->nCache = nCache;

    HASH::resize(nCache*2);

    // Allocate the total possible amount of storage. If we change cache types we won't have to reallocate.
    dataBase.resize(uCacheSizeInBytes >> 3);

    //target_T1 = nCache/2;   /* is this ok? */
    target_T1 = 0;
    }

template<>
ARC<>::ARC()
    : uLineSizeInWords(0), uLineSizeInBytes(0), uDataSizeInBytes(0),
      nCache(0), nLineBits(0), nLineMask(0), target_T1(0) {}

template<>
ARC<>::ARC(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
    initialize(uCacheSizeInBytes,uLineSizeInBytes,nLineBits);
    }

template<>
ARC<>::~ARC() {
    }

template<>
void ARC<>::lock(void *vp) {
    uint64_t *p = reinterpret_cast<uint64_t*>(vp);
    if (p>&dataBase.front() && p<=&dataBase.back()) { // Might have been a fast, read-only grab. If so ignore it.
	/* We will be given an element, but this needs to be turned into a cache line */
	p = &dataBase[(p - &dataBase.front()) / (uLineSizeInWords+1) * (uLineSizeInWords+1) + 1];
	++p[-1]; // Increment the lock count
	}
    }

template<>
void ARC<>::release(void *vp) {
    uint64_t *p = reinterpret_cast<uint64_t*>(vp);
    if (p>&dataBase.front() && p<=&dataBase.back()) { // Might have been a fast, read-only grab. If so ignore it.
	/* We will be given an element, but this needs to be turned into a cache line */
	p = &dataBase[(p - &dataBase.front()) / (uLineSizeInWords+1) * (uLineSizeInWords+1) + 1];
	uint64_t t = p[-1]-1;
	assert((t^_ARC_MAGIC_) < 0x00000000ffffffff); // Not an element or too many unlocks
	p[-1] = t;
	}
    }

template<>
void ARC<>::destage(const char *data,uint32_t uIndex,uint32_t uId) {}

template<>
void ARC<>::evict(ENTRY &temp) {
    auto &cdb = std::get<0>(temp);
    if (cdb.dirty()) {
    	destage(cdb.getData(),cdb.getPage(),cdb.getId());
	cdb.setClean();    /* No longer dirty */
	}
    }

template<>
uint64_t *ARC<>::replace(WHERE iTarget, CDBL::iterator item) {
    evict(*item);     /* if dirty, evict before overwrite */
    auto &cdb = std::get<0>(*item);
    auto data = cdb.data;
    cdb.data = NULL; /*GHOST*/
    cdb.setClean();
    move(iTarget,item);
    return data;
    }

template<>
uint64_t *ARC<>::replace(bool bInB2) {
    uint64_t *data;
    uint32_t max = (target_T1 > 1)?(target_T1+(bInB2?0:1)):1;
    auto unlocked = [](ENTRY&i) {return std::get<0>(i).data[-1] == _ARC_MAGIC_;}; // True if there are no locks
    CDBL::iterator item;
    if (L[T1].size() >= max) { // T1’s size exceeds target?
	item = std::find_if(L[T1].begin(),L[T1].end(),unlocked); // Grab the oldest unlocked entry
	if (item != L[T1].end()) return replace(B1,item);
	item = std::find_if(L[T2].begin(),L[T2].end(),unlocked); // Try the same in T2 if all T1 entries are locked
	if (item != L[T2].end()) return replace(B2,item);
	}
    else { // no: T1 is not too big
	item = std::find_if(L[T2].begin(),L[T2].end(),unlocked); // Grab the oldest unlocked entry
	if (item != L[T2].end()) return replace(B2,item);
	item = std::find_if(L[T1].begin(),L[T1].end(),unlocked); // Try the same in T1 if all T2 entries are locked
	if (item != L[T1].end()) return replace(B1,item);
	}
    std::cerr << "ERROR: all ARC entries are locked, aborting" << std::endl;
    abort();
    }

template<>
void ARC<>::RemoveAll() {
    for(auto &i : L[T1]) { evict(i); } // Flush any dirty (modified) cache entries
    for(auto &i : L[T2]) { evict(i); }
    drop_all();
    }

template<>
void *ARC<>::fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual) {
    auto uLine = uIndex >> nLineBits;
    auto iInLine = uIndex & nLineMask;
    uint32_t rat;

    /* First check our own cache */
    auto uHash = hash(uLine,uId);
    CDBL::iterator item;
    auto &Hash = HashChains[uHash];
    auto iskey = [uLine,uId](CDBL::iterator &i) {return std::get<0>(*i).uPage == uLine && (std::get<0>(*i).getId()) == uId;};
    auto match = std::find_if(Hash.begin(),Hash.end(),iskey);
    if (match != Hash.end()) {                       /* found in cache directory? */
	item = *match;
	auto &cdb = std::get<0>(**match);
	switch (cdb.where()) {                   /* yes, which list? */
	// If the element is in P1, T1 or T2 then we have a cache hit
	case P1: // Prefetched (this is an extension to ARC)
	    move(T1,item);
	    if (bLock) ++cdb.data[-1]; // Increase the lock count if requested
	    return reinterpret_cast<char*>(cdb.data) + uDataSizeInBytes*iInLine;
	case T1:
	case T2:
	    move(T2,item);
	    if (bLock) ++cdb.data[-1]; // Increase the lock count if requested
	    return reinterpret_cast<char*>(cdb.data) + uDataSizeInBytes*iInLine;
	case B1:                            /* B1 hit: favor recency */
	    if (cdb.absent()) { // We return "not found" (NULL)
		move(B2,item);
		return NULL;
		}
	    rat = L[B2].size() / L[B1].size();
	    if (rat < 1) rat = 1;
	    target_T1 += rat;
	    if (target_T1 > nCache) target_T1 = nCache;
	    break;
	case B2:                            /* B2 hit: favor frequency */
	    if (cdb.absent()) { // We return "not found" (NULL)
		move(B2,item);
		return NULL;
		}
	    rat = L[B1].size() / L[B2].size();
	    if (rat < 1) rat = 1;
	    if (rat > target_T1) target_T1 = 0;
	    else target_T1 = target_T1 - rat;
	    break;
	default:
	    std::cerr << "FATAL: found element in list " << cdb.where() << " which is invalid" << std::endl;
	    abort();
	    }
	// We only get here if the element was in B1 or B2
	move(T2,item);
	invokeRequest(uLine,uId,bVirtual); // Request the element be fetched
	cdb.data = replace(cdb.where()==B2);                                /* find a place to put new page */
	cdb.uId = uId;     /* temp->ARC_where = _T2_; and set the dirty bit for this page */
	cdb.uPage = uLine;                          /* bookkeep */
	finishRequest(uLine,uId,bVirtual,cdb.data);
	}

    else {                                                              /* page is not in cache directory */
	/*
	** Can initiate the data request right here, and do the rest while waiting...
	*/
	invokeRequest(uLine,uId,bVirtual);
	auto L1Length = L[T1].size() + L[B1].size();
	if (L1Length == nCache) {                                   /* B1 + T1 full? */
	    if (L[T1].size() < nCache) {                                           /* Still room in T1? */
		item = HASH::remove(B1);        /* yes: take page off B1 */
		auto &cdb = std::get<0>(*item);
		cdb.data = replace();                                /* find new place to put page */
		move(T1,item);
		}
	    else {                                                      /* no: B1 must be empty */
		item = HASH::remove(T1);       /* take page off T1 */
		evict(*item);     /* if dirty, evict before overwrite */
		move(T1,item);
		}
	    }
	else {                                                          /* B1 + T1 have less than c pages */
	    uint32_t nInCache = L[T1].size() + L[T2].size() + L[B1].size() + L[B2].size();
	    if (nInCache >= nCache) { // Cache is full
		auto data = replace(); // Careful: replace() can take from T1
		if (nInCache == 2*nCache) item = move(T1,remove(B2)); // Directory full so remove from B2
		else item = move(T1);                                 // cache directory not full, easy case
		auto &cdb = std::get<0>(*item);
		cdb.data = data;
		}
	    else { // Cache is not full; just grab an element and point to our internal storage
		item = move(T1);
		auto &cdb = std::get<0>(*item);
		cdb.data = &dataBase[nInCache*(uLineSizeInWords+1)+1];
		cdb.data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
		}
	    }
	assert(item == std::prev(L[T1].end()));
	auto &cdb = std::get<0>(*item);
	cdb.uId = uId;                  /* temp->dirty = dirty;  p->ARC_where = _T1_; as well! */
	cdb.uPage = uLine;
	finishRequest(uLine,uId,bVirtual,cdb.data);
	assert(!HashFree.empty());
	Hash.splice_after(Hash.before_begin(),HashFree,HashFree.before_begin());
	assert(!Hash.empty());
	Hash.front() = item;
    }
    auto &cdb = std::get<0>(*item);
    if (bLock) ++cdb.data[-1]; // Increase the lock count if requested
    // If we will modify the element, then it must eventually be flushed.
    // It is possible that we don't mark this dirty, but that it is already dirty.
    if (bModify) cdb.setDirty();

    return reinterpret_cast<char*>(cdb.data) + uDataSizeInBytes*iInLine;
    }

} // namespace mdl
