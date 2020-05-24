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
#ifndef ARC_H
#define ARC_H

#include <iostream>
#include <cstdint>
#include <algorithm>
#include <vector>
#include <tuple>
#include <boost/intrusive/list.hpp>

namespace mdl {

template<typename... KEYS> class ARC;

namespace hash {
using std::uint32_t;
using std::uint64_t;
using std::size_t;

template<typename... KEYS>
class HASH {
protected:
    struct CDB {
	static const int whereBits = 3;
	uint64_t *data;                 // page's location in cache
	uint32_t uId    : 31-whereBits; // Processor that owns this
	uint32_t iWhere :    whereBits; // Which list we are on
	bool     bDirty :            1; // True if we need to flush
	uint32_t uPage;                 // page's ID number (or hash id for advanced keys)

	explicit CDB(uint64_t *data=nullptr) : uId(0), iWhere(0), bDirty(0), uPage(0), data(data) {}
	uint32_t    where()   const {return iWhere;}
	uint32_t    getId()   const {return uId;}
	uint32_t    getPage() const {return uPage;}
	const char *getData() const {return reinterpret_cast<char*>(data);}
	bool        dirty()   const {return bDirty;}  // Data has been modified (data not null)
	bool        absent()  const {return dirty();} // ... or is absent (data is null)
	void setDirty() { bDirty = true; }
	void setClean() { bDirty = false; }
	};

protected:
    struct HashChain { struct HashChain *next; };
    typedef std::tuple<KEYS...> KEY;
    typedef std::tuple<CDB,KEY> PAIR;
    class ENTRY : public PAIR, public boost::intrusive::list_base_hook<>, public HashChain {};
    typedef boost::intrusive::list<ENTRY,boost::intrusive::constant_time_size<true> > CDBL;

    template<class T>
    constexpr T& get(HASH::ENTRY& t) noexcept
    { return std::get<T,CDB,KEY>(static_cast<PAIR&>(t)); }

protected:
    uint32_t uHashMask = 0;
private:
    std::vector<ENTRY> entries;
    CDBL freeList;
protected:
    CDBL L[1<<CDB::whereBits];
    std::vector<HashChain> HashChains; // Collision chain (one for each hash value)
protected:
    // Locate by advanced key, or by simple key
    ENTRY * find_key(HashChain &Hash,const KEYS&... keys, bool bRemove=false);
    ENTRY * find_key(HashChain &Hash,uint32_t uLine,uint32_t uId, bool bRemove=false);

    typename CDBL::iterator move(uint32_t iTarget);
    typename CDBL::iterator move(uint32_t iTarget, typename CDBL::iterator item);

    // Remove the CDB entry (given by a list iterator) or front of a list from the hash table
    typename CDBL::iterator take(const typename CDBL::iterator p);
    typename CDBL::iterator take(uint32_t iTarget);

    typename CDBL::iterator free(typename CDBL::iterator item);

public:
    void  insert(uint32_t uHash, const KEYS&... keys, void *data);
    void *lookup(uint32_t uHash, const KEYS&... keys);
    void  remove(uint32_t uHash, const KEYS&... keys);

public:
    class ARC<KEYS...> *clone();
    void print_statistics();
    void clear(); // Empty the hash table (by moving all elements to the free list)
    void resize(size_t count);
    virtual ~HASH() = default;

    HASH(const HASH&)=delete;
    HASH& operator=(const HASH&)=delete;
    HASH(HASH&&)=delete;
    HASH& operator=(HASH&&)=delete;

    explicit HASH() = default;
    explicit HASH(int nMaxElements) {
	resize(nMaxElements);
	}
    };

} // namespace hash

/*****************************************************************************\
* ARC CACHE
\*****************************************************************************/

// This class tells the ARC cache how to:
// 1. Find the remote processor (thread) based on the simple or advanced key
// 2. Request that a remote element be fetched
// 3. Finish the request and copy the data into ARC
// 4. Flush modified cache entries
class ARChelper {
public:
    virtual uint32_t  getThread(uint32_t uLine, uint32_t uId, const void *pKey) {return uId;}
    virtual void  invokeRequest(uint32_t uLine, uint32_t uId, const void *pKey, bool bVirtual)             = 0;
    virtual void  finishRequest(uint32_t uLine, uint32_t uId, const void *pKey, bool bVirtual, void *data) = 0;
    virtual void   flushElement(uint32_t uLine, uint32_t uId, const void *pKey,          const void *data) = 0;
    virtual void combineElement(uint32_t uLine, uint32_t uId, const void *pKey,          const void *data) = 0;
    };

template<typename ...KEYS>
class ARC : private hash::HASH<KEYS...> {
private:
    static const uint64_t _ARC_MAGIC_ = 0xa6c3a91c00000000;
    static const uint64_t _ARC_MASK_  = 0xa6c3a91cffffffff;
    enum WHERE {
	P1, T1, T2, B1, B2, // ARC cache lists
	};

private:
    ARChelper *helper = nullptr;
    std::vector<uint64_t> dataBase; // Contains all cached data (continguous)
    using HASH = typename hash::HASH<KEYS...>;
    using CDB =  typename HASH::CDB;
    using CDBL = typename HASH::CDBL;
    using KEY = typename HASH::KEY;
    using PAIR = typename HASH::PAIR;
    using ENTRY = typename HASH::ENTRY;
    using HashChain = typename HASH::HashChain;
    using HASH::L;
    using HASH::move;
    using HASH::take;
    using HASH::find_key;
    using HASH::uHashMask;
    using HASH::HashChains;
    template<class T>
    constexpr T& get(ENTRY& t) noexcept
    { return std::get<T,CDB,KEY>(static_cast<PAIR&>(t)); }

    typename std::vector<CDBL>::size_type target_T1;
    uint32_t nCache = 0;
    uint32_t nLineBits = 0;
    uint32_t nLineMask = 0;
    uint32_t uLineSizeInWords = 0;
    uint32_t uDataSizeInBytes = 0;
private:
    // If the entry is marked as dirty, then flush it (send it back to its owner)
    void flush(ENTRY &temp);

    // Find a victim to give us some cache space. Victim is evicted.
    auto replace(WHERE iTarget, typename CDBL::iterator item);
    auto replace(bool bInB2=false);

    // The element is already in the cache so we process a "hit"
    // On return, data will be valid, or it will be NULL and we have to fetch it.
    void update(typename CDBL::iterator item,bool bLock,bool bModify);

    // Here we want to insert a new key that is not in the ARC cache.
    // This routine will return a free entry that we can use.
    auto insert(HashChain &Hash);

public:
    ARC(const ARC&)=delete;
    ARC& operator=(const ARC&)=delete;
    ARC(ARC&&)=delete;
    ARC& operator=(ARC&&)=delete;

    explicit ARC() = default;
    virtual ~ARC() = default;

    void initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    void clear(); // Empty the cache. Evict/flush any dirty elements.

    // Decrements the lock count of an element. When the lock count is zero then it is eligible to be flushed.
    void release(void *vp);

    // Fetch using simple or advanced key
    void *fetch(uint32_t uIndex, int uId,bool bLock,bool bModify,bool bVirtual);
    void *fetch(uint32_t uHash,const KEYS&... keys,bool bLock,bool bModify,bool bVirtual);
    };

/*****************************************************************************\
* HASH implementation
\*****************************************************************************/
namespace hash {

    template<typename... KEYS>
    class ARC<KEYS...> *HASH<KEYS...>::clone() {
    	return new ARC<KEYS...>;
	}

    // Locate by advanced key
    template<typename... KEYS>
    typename HASH<KEYS...>::ENTRY * HASH<KEYS...>::find_key(HashChain &Hash,const KEYS&... keys, bool bRemove) {
	auto key = std::make_tuple(keys...);
	for(auto pHash = &Hash; pHash->next != &Hash; pHash = pHash->next) {
	    auto pEntry = static_cast<ENTRY*>(pHash->next);
	    if (get<KEY>(*pEntry) == key) { // We ignore the CDB() part (see operator==() in CDB)
	    	if (bRemove) pHash->next = pHash->next->next;
	    	return pEntry;
		}
	    }
	return nullptr;
	}

    // Locate by simple key (index and processor id)
    template<typename... KEYS>
    typename HASH<KEYS...>::ENTRY * HASH<KEYS...>::find_key(HashChain &Hash,uint32_t uLine,uint32_t uId, bool bRemove) {
	for(auto pHash = &Hash; pHash->next != &Hash; pHash = pHash->next) {
	    auto pEntry = static_cast<ENTRY*>(pHash->next);
	    auto const &cdb = get<CDB>(*pEntry);
	    if (cdb.uPage == uLine && (cdb.getId()) == uId) {
	    	if (bRemove) pHash->next = pHash->next->next;     	     
	    	return pEntry;
		}
	    }
	return nullptr;
	}

    template<typename... KEYS>
    void HASH<KEYS...>::print_statistics() {
	auto nInUse = entries.size() - freeList.size();
	std::cout << "Cache Size: " << entries.size()
	    << ", in use: " << nInUse
	    << std::endl;
	size_t nBucketsHit = 0;
	for(auto &i : HashChains) {
	    if (i.next != &i) ++nBucketsHit;
	    }
	std::cout << "Buckets: " << HashChains.size()
	    << ", buckets in use: " << nBucketsHit << std::endl;
	if (nInUse) {
	    std::cout << "Collisions: " << nInUse - nBucketsHit
		<< ", collision rate: " << 100.0 - 100.0 * nBucketsHit / nInUse << "%"
		<< std::endl;
	    }
	}

    template<typename... KEYS>
    typename HASH<KEYS...>::CDBL::iterator
    HASH<KEYS...>::move(uint32_t iTarget) {
	L[iTarget].splice(L[iTarget].end(),freeList,freeList.begin());
	get<CDB>(L[iTarget].back()).iWhere = iTarget;
	return std::prev(L[iTarget].end());
	}


    template<typename... KEYS>
    typename HASH<KEYS...>::CDBL::iterator
    HASH<KEYS...>::move(uint32_t iTarget, typename HASH<KEYS...>::CDBL::iterator item) {
	auto &cdb = get<CDB>(*item);
	L[iTarget].splice(L[iTarget].end(),L[cdb.iWhere],item);
	cdb.iWhere = iTarget;
	return item;
	}

    // Remove the CDB entry (given by a list iterator) from the hash table
    template<typename... KEYS>
    typename HASH<KEYS...>::CDBL::iterator
    HASH<KEYS...>::take(const typename HASH<KEYS...>::CDBL::iterator p) {
	// Find the hash value by trolling through the chain. The chain is supposed to be short.
	HashChain *pHash, *me = &*p;
	for (pHash=me->next; pHash!=me && pHash < &HashChains.front() && pHash > &HashChains.back(); pHash=pHash->next) {}
	uint32_t uHash = pHash - &HashChains.front();
	auto &Hash = HashChains[uHash&uHashMask];
	for(pHash = &Hash; pHash->next != me; pHash = pHash->next) {} // Find previous value
	pHash->next = pHash->next->next; // Remove this entry from the hash chain
	return p;
	}

    // Typically we want to take the LRU (front) element from a list, in anticipation of reusing the entry
    template<typename... KEYS>
    typename HASH<KEYS...>::CDBL::iterator
    HASH<KEYS...>::take(uint32_t iTarget) {
	return take(L[iTarget].begin());
	}

    template<typename... KEYS>
    typename HASH<KEYS...>::CDBL::iterator
    HASH<KEYS...>::free(typename HASH<KEYS...>::CDBL::iterator item) {
	auto &cdb = get<CDB>(*item);
	freeList.splice(freeList.end(),L[cdb.iWhere],item);
	return item;
	}

    template<typename... KEYS>
    void HASH<KEYS...>::clear() { // Empty the hash table (by moving all elements to the free list)
	for(auto &i : HashChains) i.next = &i; // Point to self = empty
	for(auto &list : L) freeList.splice(freeList.end(),list);
	}

    template<typename... KEYS>
    void HASH<KEYS...>::resize(size_t count) {
	clear();
	uHashMask = count + count/2;    // Reserve extra for collision resolution
	uHashMask |= (uHashMask >> 1);  // Calculate the bitmask for the next higher
	uHashMask |= (uHashMask >> 2);  //   power of two so we can mask for an index.
	uHashMask |= (uHashMask >> 4);
	uHashMask |= (uHashMask >> 8);
	uHashMask |= (uHashMask >> 16);
	HashChains.resize(uHashMask+1);        // Collision chains for each possible hashed value
	for(auto &i : HashChains) i.next = &i; // Point to self = empty
	entries.reserve(count);
	entries.resize(count); // Our complete list of ENTRY elements
	freeList.clear();      // Add them all to the free list
	for(auto &i : entries) freeList.push_back(i);
	}


    template<typename... KEYS>
    void HASH<KEYS...>::insert(uint32_t uHash, const KEYS&... keys,void *data) {
	auto item = move(0);                          // Grab a new item from the free list
	auto &cdb = get<CDB>(*item);             // CDB (cache data block)
	cdb.uId = 0;                                  // Should set this to processor id probably
	cdb.uPage = uHash;                            // Page is the hash value for KEY types
	cdb.data = reinterpret_cast<uint64_t*>(data); // Points to the data for this element
	get<KEY>(*item) = std::make_tuple(keys...);
	auto &Hash = HashChains[uHash&uHashMask];
	assert(find_key(Hash,keys...)==nullptr);         // Duplicate keys are not allowed
	item->next = Hash.next;
	Hash.next = &*item;
	}

    template<typename... KEYS>
    void *HASH<KEYS...>::lookup(uint32_t uHash, const KEYS&... keys) {
	auto &Hash = HashChains[uHash&uHashMask];
	auto pEntry = find_key(Hash,keys...);
	if (pEntry) return get<CDB>(*pEntry).data;
	return nullptr;
	}

    template<typename... KEYS>
    void HASH<KEYS...>::remove(uint32_t uHash, const KEYS&... keys) {
	auto &Hash = HashChains[uHash&uHashMask];
	auto pEntry = find_key(Hash,keys...,true);
	free(CDBL::s_iterator_to(*pEntry));
	}

    } // namespace hash

/*****************************************************************************\
* ARC implementation
\*****************************************************************************/

    // If the entry is marked as dirty, then flush it (send it back to its owner)
    template<typename ...KEYS>
    void ARC<KEYS...>::flush(ENTRY &item) {
	auto &cdb = get<CDB>(item);
	if (cdb.dirty()) {
	    helper->flushElement(cdb.getPage(),cdb.getId(),&cdb+1,cdb.getData());
	    cdb.setClean();    /* No longer dirty */
	    }
	}

    template<typename ...KEYS>
    auto ARC<KEYS...>::replace(WHERE iTarget, typename CDBL::iterator item) {
	flush(*item); // Flush this element if it is dirty
	auto &cdb = get<CDB>(*item);
	auto data = cdb.data;
	cdb.data = nullptr; /*GHOST*/
	cdb.setClean();
	move(iTarget,item);
	return data;
	}

    // Find a victim to give us some cache space. Victim is evicted.
    template<typename ...KEYS>
    auto ARC<KEYS...>::replace(bool bInB2) {
	uint64_t *data;
	uint32_t max = (target_T1 > 1)?(target_T1+(bInB2?0:1)):1;
	auto unlocked = [](PAIR&i) {return std::get<CDB>(i).data[-1] == _ARC_MAGIC_;}; // True if there are no locks
	typename CDBL::iterator item;
	if (L[T1].size() >= max) { // T1’s size exceeds target?
	    // Take the oldest unlocked entry in T1 or T2 if T1 is full
	    item = std::find_if(L[T1].begin(),L[T1].end(),unlocked);
	    if (item != L[T1].end()) return replace(B1,item);
	    item = std::find_if(L[T2].begin(),L[T2].end(),unlocked);
	    if (item != L[T2].end()) return replace(B2,item);
	    }
	else { // no: T1 is not too big
	    // Take the oldest unlocked entry in T2 or T1 if T2 is full
	    item = std::find_if(L[T2].begin(),L[T2].end(),unlocked);
	    if (item != L[T2].end()) return replace(B2,item);
	    item = std::find_if(L[T1].begin(),L[T1].end(),unlocked);
	    if (item != L[T1].end()) return replace(B1,item);
	    }
	std::cerr << "ERROR: all ARC entries are locked, aborting" << std::endl;
	abort();
	}

    // The element is already in the cache so we process a "hit"
    // On return, data will be valid, or it will be NULL and we have to fetch it.
    template<typename ...KEYS>
    void ARC<KEYS...>::update(typename CDBL::iterator item,bool bLock,bool bModify) {
	uint32_t rat;
	auto &cdb = get<CDB>(*item);
	switch (cdb.where()) {              // Which list is the element on
	// If the element is in P1, T1 or T2 then we have a cache hit
	case P1: // Prefetched (this is an extension to ARC)
	    move(T1,item);
	    break;
	case T1:
	case T2:
	    move(T2,item);
	    break;
	case B1:                            /* B1 hit: favor recency */
	    if (cdb.absent()) { // We return "not found" (NULL)
		move(B2,item);
		break;
		}
	    rat = L[B2].size() / L[B1].size();
	    if (rat < 1) rat = 1;
	    target_T1 += rat;
	    if (target_T1 > nCache) target_T1 = nCache;
	    move(T2,item);
	    break;
	case B2:                            /* B2 hit: favor frequency */
	    if (cdb.absent()) { // We return "not found" (NULL)
		move(B2,item);
		break;
		}
	    rat = L[B1].size() / L[B2].size();
	    if (rat < 1) rat = 1;
	    if (rat > target_T1) target_T1 = 0;
	    else target_T1 = target_T1 - rat;
	    move(T2,item);
	    break;
	default:
	    std::cerr << "FATAL: found element in list " << cdb.where() << " which is invalid" << std::endl;
	    abort();
	    }
	}

    // Here we want to insert a new key that is not in the ARC cache.
    // This routine will return a free entry that we can use.
    template<typename ...KEYS>
    auto ARC<KEYS...>::insert(HashChain &Hash) {
	typename CDBL::iterator item;
	auto L1Length = L[T1].size() + L[B1].size();
	// If we have used all of the available cache space then we need to evict an entry for reuse
	if (L1Length == nCache) {               // T1 + B1 contain entries with data.
	    if (L[T1].size() < nCache) {  // T1 still has room?
		item = take(B1);          // Yes, take from B1
		auto &cdb = get<CDB>(*item);
		cdb.data = replace();     // Find a some cache space. Old value is evicted.
		move(T1,item);            // Move to the end of T1 (MRU)
		}
	    else { // B1 must be empty
		item = take(T1);          // Take a value from T1
		flush(*item);             // flush if it was dirty
		move(T1,item);            // Move to the end of T1 (MRU)
		}
	    }
	else {
	    uint32_t nInCache = L[T1].size() + L[T2].size() + L[B1].size() + L[B2].size();
	    if (nInCache >= nCache) {           // Cache is full
		auto data = replace();          // Careful: replace() can take from T1
		if (nInCache == 2*nCache)
		    item = move(T1,take(B2));   // Directory full so take from B2
		else item = move(T1);           // cache directory not full, easy case
		auto &cdb = get<CDB>(*item);
		cdb.data = data;
		}
	    else { // Cache is not full; just grab an element and point to our internal storage
		item = move(T1);
		auto &cdb = get<CDB>(*item);
		cdb.data = &dataBase[nInCache*(uLineSizeInWords+1)+1];
		cdb.data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
		}
	    }
	item->next = Hash.next; // Add to hash table
	Hash.next = &*item;
	return item;
	}

    template<typename ...KEYS>
    void ARC<KEYS...>::initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
	this->helper = helper;
	this->uLineSizeInWords = (uLineSizeInBytes+7) >> 3;       // Size of a cache line (aligned properly)
	// Calculate nCache based on the number of cache lines that will fit in our cache buffer
	// Account for the "magic" number before the cache line.
	auto nCache = (uCacheSizeInBytes>>3) / (this->uLineSizeInWords+1);
	if (this->nCache == nCache) return; // Already setup
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

    // Empty the cache. Evict/flush any dirty elements.
    template<typename ...KEYS>
    void ARC<KEYS...>::clear() {
	for(auto &i : L[T1]) { flush(i); } // Flush any dirty (modified) cache entries
	for(auto &i : L[T2]) { flush(i); }
	HASH::clear(); // Now empty the hash table
	}

    // Decrements the lock count of an element. When the lock count is zero then it is eligible to be flushed.
    template<typename ...KEYS>
    void ARC<KEYS...>::release(void *vp) {
	uint64_t *p = reinterpret_cast<uint64_t*>(vp);
	if (p>&dataBase.front() && p<=&dataBase.back()) { // Might have been a fast, read-only grab. If so ignore it.
	    /* We will be given an element, but this needs to be turned into a cache line */
	    p = &dataBase[(p - &dataBase.front()) / (uLineSizeInWords+1) * (uLineSizeInWords+1) + 1];
	    uint64_t t = p[-1]-1;
	    assert((t^_ARC_MAGIC_) < 0x00000000ffffffff); // Not an element or too many unlocks
	    p[-1] = t;
	    }
	}

/*****************************************************************************\
* ARC fetch
\*****************************************************************************/

template<typename... KEYS>
void *ARC<KEYS...>::fetch(uint32_t uHash, const KEYS&... keys, bool bLock,bool bModify,bool bVirtual) {
    uint32_t rat;

    /* First check our own cache */
    typename CDBL::iterator item;
    auto &Hash = HashChains[uHash&uHashMask];
    auto pEntry = find_key(Hash,keys...);
    auto key = std::make_tuple(keys...);

    if (pEntry) {   // Page is reference by the cache
	item = CDBL::s_iterator_to(*pEntry);
	update(item,bLock,bModify);
	auto &cdb = get<CDB>(*item);
	if (cdb.data == nullptr) {
	    if (cdb.absent()) return nullptr;
	    auto uId = helper->getThread(uHash,0,&key);
	    helper->invokeRequest(uHash,uId,&key,bVirtual); // Request the element be fetched
	    cdb.data = replace(cdb.where()==B2);
	    cdb.uId = uId;
	    cdb.uPage = uHash;
	    helper->finishRequest(uHash,uId,&key,bVirtual,cdb.data);
	    }
	}
    else {          // Page is not in the cache
	auto uId = helper->getThread(uHash,0,&key);
	helper->invokeRequest(uHash,uId,&key,bVirtual);
	item = insert(Hash);

	auto &cdb = get<CDB>(*item);
	cdb.uId = uId;
	cdb.uPage = uHash;
	get<KEY>(*item) = std::make_tuple(keys...);


	helper->finishRequest(uHash,uId,&key,bVirtual,cdb.data);
	}
    auto &cdb = get<CDB>(*item);
    if (bLock) ++cdb.data[-1];   // Increase the lock count if requested
    if (bModify) cdb.setDirty(); // Mark page as dirty

    return cdb.data;
    }



} // namespace mdl

#endif
