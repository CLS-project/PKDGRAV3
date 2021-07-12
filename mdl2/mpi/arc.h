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
#include <typeinfo>
#include <cstdint>
#include <cstring> // memcpy
#include <algorithm>
#include <vector>
#include <tuple>
#include <boost/intrusive/list.hpp>

namespace mdl {

template<typename... KEYS> class ARC;

// This class tells the ARC cache how to:
// 1. Find the remote processor (thread) based on the simple or advanced key
// 2. Request that a remote element be fetched
// 3. Finish the request and copy the data into ARC
// 4. Flush modified cache entries
class ARChelper {
public:
    virtual uint32_t  getThread(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) {return uId;}
    virtual void * getLocalData(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey) {return nullptr;}
    virtual void *invokeRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual)             = 0;
    virtual void *finishRequest(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey, bool bVirtual, void *dst, const void *src) = 0;
    virtual void   flushElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) = 0;
    virtual void combineElement(uint32_t uLine, uint32_t uId, uint32_t size, const void *pKey,          const void *data) = 0;
    };

class GARC {
public:
    virtual ~GARC() = default;
    virtual void *fetch(uint32_t uIndex,    uint32_t uId,bool bLock,bool bModify,bool bVirtual) = 0;
    virtual void *fetch(uint32_t uHash, const void *pKey,bool bLock,bool bModify,bool bVirtual) = 0;
    virtual void *inject(uint32_t uHash, uint32_t uId, const void *pKey) = 0;
    virtual void initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0) = 0;
    virtual void release(void *vp) = 0;
    virtual void clear() = 0;
    virtual uint32_t key_size() = 0;
    };

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

namespace hash {
using std::uint32_t;
using std::uint64_t;
using std::size_t;

static inline uint32_t hash(uint32_t uLine,uint32_t uId) {
    uint32_t key[] = {uLine,uId};
    return murmur::murmur2<2>(key);
    }

class GHASH {
public:
    virtual ~GHASH() = default;
    virtual class GARC *clone(GARC *arc=nullptr) = 0;
    virtual void  insert(uint32_t uHash, const void *pKey, void *data) = 0;
    virtual void *lookup(uint32_t uHash, const void *pKey) = 0;
    virtual void  remove(uint32_t uHash, const void *pKey) = 0;
    };


template<typename... KEYS>
class HASH : public GHASH {
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
	const void *getData() const {return data;}
	bool        dirty()   const {return bDirty;}  // Data has been modified (data not null)
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
    ENTRY * find_key(HashChain &Hash,const KEY& key, bool bRemove=false);
    ENTRY * find_key(HashChain &Hash,const KEYS&... keys, bool bRemove=false);
    ENTRY * find_key(HashChain &Hash,uint32_t uLine,uint32_t uId, bool bRemove=false);

    typename CDBL::iterator move(uint32_t iTarget);
    typename CDBL::iterator move(uint32_t iTarget, typename CDBL::iterator item);

    // Remove the CDB entry (given by a list iterator) or front of a list from the hash table
    // and return the entry. It is ready to be moved to another list.
    typename CDBL::iterator take(const typename CDBL::iterator p);
    typename CDBL::iterator take(uint32_t iTarget);

    typename CDBL::iterator free(typename CDBL::iterator item);

public:
    void  insert(uint32_t uHash, const KEY& key,      void *data);
    void  insert(uint32_t uHash, const KEYS&... keys, void *data);
    void *lookup(uint32_t uHash, const KEY& key);
    void *lookup(uint32_t uHash, const KEYS&... keys);
    void  remove(uint32_t uHash, const KEY& key);
    void  remove(uint32_t uHash, const KEYS&... keys);
private:
    virtual void  insert(uint32_t uHash, const void *pKey, void *data) override;
    virtual void *lookup(uint32_t uHash, const void *pKey) override;
    virtual void  remove(uint32_t uHash, const void *pKey) override;

public:
    virtual class GARC *clone(GARC *arc=nullptr) override;
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

template<typename ...KEYS>
class ARC : private hash::HASH<KEYS...>, public GARC {
private:
    static const uint64_t _ARC_MAGIC_ = 0xa6c3a91c00000000;
    static const uint64_t _ARC_MASK_  = 0xa6c3a91cffffffff;
    enum WHERE {
	P1, A1, T1, T2, B1, B2, // ARC cache lists
	};

private:
    ARChelper *helper = nullptr;
    std::vector<uint64_t> dataBase; // Contains all cached data (continguous)
    std::vector<uint64_t> cacheLine;// Room for one cache line

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
    uint32_t nGhost = 0;
    uint32_t nAbsent= 0;
    uint32_t nLineBits = 0;
    uint32_t nLineMask = 0;
    uint32_t uLineSizeInWords = 0;
    uint32_t uDataSizeInBytes = 0;
private:
    // If the entry is marked as dirty, then flush it (send it back to its owner)
    void flush(ENTRY &temp);

    // Find a victim to give us some cache space. Victim is evicted.
    auto victim(bool bInB2=false);
    auto replace(WHERE iTarget, typename CDBL::iterator item);
    auto replace(typename CDBL::iterator item);
    auto replace(bool bInB2=false);

    // The element is already in the cache so we process a "hit"
    // On return, data will be valid, or it will be NULL and we have to fetch it.
    void update(typename CDBL::iterator item,bool bLock,bool bModify);

    // Here we want to insert a new key that is not in the ARC cache.
    // This routine will return a free entry that we can use.
    auto insert(HashChain &Hash);
    auto insert_absent(HashChain &Hash);

    void *fetch(uint32_t uHash,uint32_t uId,const KEY &key,bool bLock,bool bModify,bool bVirtual);
    void *inject(uint32_t uHash, uint32_t uId, const KEY &key);

public:
    ARC(const ARC&)=delete;
    ARC& operator=(const ARC&)=delete;
    ARC(ARC&&)=delete;
    ARC& operator=(ARC&&)=delete;

    explicit ARC() = default;
    virtual ~ARC() = default;

    virtual void initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0) override;
    virtual void clear() override; // Empty the cache. Evict/flush any dirty elements.
    virtual uint32_t key_size() override {return std::tuple_size<KEY>::value ? sizeof(KEY) : 0;}

    // Decrements the lock count of an element. When the lock count is zero then it is eligible to be flushed.
    virtual void release(void *vp) override;

    // Fetch using simple or advanced key
    void *fetch(uint32_t uHash,        const KEYS&... keys,bool bLock,bool bModify,bool bVirtual);
    virtual void *fetch(uint32_t uIndex,    uint32_t uId,bool bLock,bool bModify,bool bVirtual) override;
    virtual void *fetch(uint32_t uHash, const void *pKey,bool bLock,bool bModify,bool bVirtual) override;
    virtual void *inject(uint32_t uHash, uint32_t uId, const void *pKey) override;
    };

/*****************************************************************************\
* HASH implementation
\*****************************************************************************/
namespace hash {

    // conditionally clone
    template<typename... KEYS>
    class GARC *HASH<KEYS...>::clone(GARC *arc) {
	if (arc && typeid(*arc) == typeid(ARC<KEYS...>)) return nullptr;
	return new ARC<KEYS...>;
	}

    // Locate by advanced key
    template<typename... KEYS>
    typename HASH<KEYS...>::ENTRY * HASH<KEYS...>::find_key(HashChain &Hash,const KEY& key, bool bRemove) {
	for(auto pHash = &Hash; pHash->next != &Hash; pHash = pHash->next) {
	    auto pEntry = static_cast<ENTRY*>(pHash->next);
	    if (get<KEY>(*pEntry) == key) {
	    	if (bRemove) pHash->next = pHash->next->next;
	    	return pEntry;
		}
	    }
	return nullptr;
	}

    template<typename... KEYS>
    typename HASH<KEYS...>::ENTRY * HASH<KEYS...>::find_key(HashChain &Hash,const KEYS&... keys, bool bRemove) {
	return find_key(Hash,std::make_tuple(keys...),bRemove);
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
	for (pHash=me->next; pHash < &HashChains.front() || pHash > &HashChains.back(); pHash=pHash->next) {assert(pHash!=me);}
	uint32_t uHash = pHash - &HashChains.front();
	auto &Hash = HashChains[uHash];
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
//	clear();
	for(auto &list : L) list.clear();
	freeList.clear();
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
	for(auto &i : entries) freeList.push_back(i);
	}


    /*****************************************************************************/
    template<typename... KEYS>
    void HASH<KEYS...>::insert(uint32_t uHash, const KEY& key,void *data) {
	auto item = move(0);                          // Grab a new item from the free list
	auto &cdb = get<CDB>(*item);             // CDB (cache data block)
	cdb.uId = 0;                                  // Should set this to processor id probably
	cdb.uPage = uHash;                            // Page is the hash value for KEY types
	cdb.data = static_cast<uint64_t*>(data); // Points to the data for this element
	get<KEY>(*item) = key;
	auto &Hash = HashChains[uHash&uHashMask];
	assert(find_key(Hash,key)==nullptr);         // Duplicate keys are not allowed
	item->next = Hash.next;
	Hash.next = &*item;
	}
    template<typename... KEYS>
    void HASH<KEYS...>::insert(uint32_t uHash, const KEYS&... keys,void *data) {
    	insert(uHash,std::make_tuple(keys...),data);
	}
    template<typename... KEYS>
    void HASH<KEYS...>::insert(uint32_t uHash, const void *pKey,void *data) {
    	auto const &key = * static_cast<const KEY*>(pKey);
    	insert(uHash,key,data);
	}
    /*****************************************************************************/
    template<typename... KEYS>
    void *HASH<KEYS...>::lookup(uint32_t uHash, const KEY& key) {
	auto &Hash = HashChains[uHash&uHashMask];
	auto pEntry = find_key(Hash,key);
	if (pEntry) return get<CDB>(*pEntry).data;
	return nullptr;
	}
    template<typename... KEYS>
    void *HASH<KEYS...>::lookup(uint32_t uHash, const KEYS&... keys) {
    	return lookup(uHash,std::make_tuple(keys...));
	}
    template<typename... KEYS>
    void *HASH<KEYS...>::lookup(uint32_t uHash, const void *pKey) {
    	auto const &key = * static_cast<const KEY*>(pKey);
    	return lookup(uHash,key);
	}
    /*****************************************************************************/
    template<typename... KEYS>
    void HASH<KEYS...>::remove(uint32_t uHash, const KEY& key) {
	auto &Hash = HashChains[uHash&uHashMask];
	auto pEntry = find_key(Hash,key,true);
	assert(pEntry);
	free(CDBL::s_iterator_to(*pEntry));
	}
    template<typename... KEYS>
    void HASH<KEYS...>::remove(uint32_t uHash, const KEYS&... keys) {
    	remove(uHash,std::make_tuple(keys...));
	}
    template<typename... KEYS>
    void HASH<KEYS...>::remove(uint32_t uHash, const void *pKey) {
    	auto const &key = * static_cast<const KEY*>(pKey);
    	remove(uHash,key);
	}
    /*****************************************************************************/

    } // namespace hash

/*****************************************************************************\
* ARC implementation
\*****************************************************************************/

    // If the entry is marked as dirty, then flush it (send it back to its owner)
    template<typename ...KEYS>
    void ARC<KEYS...>::flush(ENTRY &item) {
	auto &cdb = get<CDB>(item);
	if (cdb.dirty()) {
	    auto &key = get<KEY>(item);
	    helper->flushElement(cdb.getPage(),cdb.getId(),
	    	std::tuple_size<KEY>::value ? sizeof(key) : 0,&key,cdb.getData());
	    cdb.setClean();    /* No longer dirty */
	    }
	}

    // We are going to reuse the data from "item" and move it from T -> B
    // The data pointer is returned so we can use it somewhere else.
    template<typename ...KEYS>
    auto ARC<KEYS...>::replace(WHERE iTarget, typename CDBL::iterator item) {
	assert(iTarget==B1 || iTarget==B2);
	flush(*item); // Flush this element if it is dirty
	auto &cdb = get<CDB>(*item);
	auto data = cdb.data;
	cdb.data = nullptr; /*GHOST*/
	move(iTarget,item); // Target here will be B (B1 or B2)
	return data;
	}

    // We know which item we want to use for replacement.
    // Move it to the right B cache and return the data.
    template<typename ...KEYS>
    auto ARC<KEYS...>::replace(typename CDBL::iterator item) {
	auto &cdb = get<CDB>(*item);
    	switch(cdb.where()) {
    	     case T1: return replace(B1,item); // T1 -> B1
    	     case T2: return replace(B2,item); // T2 -> B2
	    }
	std::cerr << "FATAL: found element in list " << cdb.where() << " which is invalid" << std::endl;
	abort();
	}

    // We don't have any cache data space left and potentially want to insert a new element.
    // This routine will find the appropriate "victim" in T1 or T2 and return it.
    // If there are no unlocked elements it is a fatal error and we abort.
    template<typename ...KEYS>
    auto ARC<KEYS...>::victim(bool bInB2) {
	uint32_t max = (target_T1 > 1)?(target_T1+(bInB2?0:1)):1;
	auto unlocked = [](PAIR&i) {return std::get<CDB>(i).data[-1] == _ARC_MAGIC_;}; // True if there are no locks
	typename CDBL::iterator item;
	if (L[T1].size() >= max && (item = std::find_if(L[T1].begin(),L[T1].end(),unlocked)) != L[T1].end()) return item;
	if ((item = std::find_if(L[T2].begin(),L[T2].end(),unlocked)) != L[T2].end()) return item;
	if (L[T1].size() < max && (item = std::find_if(L[T1].begin(),L[T1].end(),unlocked)) != L[T1].end()) return item;
	std::cerr << "ERROR: all ARC entries are locked, aborting" << std::endl;
	abort();
	}

    template<typename ...KEYS>
    auto ARC<KEYS...>::replace(bool bInB2) {
	return replace(victim(bInB2));
	}

    // The element is already in the cache so we process a "hit"
    // On return, data will be valid, or it will be NULL and we have to fetch it.
    template<typename ...KEYS>
    void ARC<KEYS...>::update(typename CDBL::iterator item,bool bLock,bool bModify) {
	auto delta = [this](int i,int j) { return std::max<size_t>(L[i].size() / L[j].size(),1); };
	auto &cdb = get<CDB>(*item);
	switch (cdb.where()) {              // Which list is the element on
	// If the element is in P1, T1 or T2 then we have a cache hit
	case P1: // Prefetched (this is an extension to ARC)
	    move(T1,item);
	    break;
	case A1: // Absent (this is an extension to ARC)
	    move(A1,item);
	    break;
	case T1:
	case T2:
	    move(T2,item);
	    break;
	case B1:                            /* B1 hit: favor recency */
	    target_T1 += std::min<size_t>(delta(B2,B1),nCache-target_T1);
	    move(T2,item);
	    break;
	case B2:                            /* B2 hit: favor frequency */
	    target_T1 -= std::min<size_t>(delta(B1,B2),target_T1);
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
	if (L1Length == nCache) {
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
		if (nInCache == nCache+nGhost)
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

    // Here we want to insert a new key that is not in the ARC cache.
    // This routine will return a free entry that we can use.
    template<typename ...KEYS>
    auto ARC<KEYS...>::insert_absent(HashChain &Hash) {
	assert(Hash.next);
	typename CDBL::iterator item;
	if (L[A1].size() < nAbsent) item = move(A1);
	else item = move(A1,take(A1));
	auto &cdb = get<CDB>(*item);
	cdb.data = nullptr;
	item->next = Hash.next; // Add to hash table
	Hash.next = &*item;
	return item;
	}
    template<typename ...KEYS>
    void ARC<KEYS...>::initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
	this->helper = helper;
	// Size of a cache line (aligned properly)
        this->uLineSizeInWords = (uLineSizeInBytes+sizeof(uint64_t)-1) / sizeof(uint64_t);
	cacheLine.resize(this->uLineSizeInWords);
	// Calculate nCache based on the number of cache lines that will fit in our cache buffer
	// Account for the "magic" number before the cache line.
	auto nCache = (uCacheSizeInBytes/sizeof(uint64_t)) / (this->uLineSizeInWords+1);
	if (this->nCache != nCache) {
	    this->uDataSizeInBytes = (uLineSizeInBytes) >> nLineBits; // Size of a single element
	    this->nLineBits = nLineBits;
	    this->nLineMask = (1 << nLineBits) - 1;
	    this->nCache = nCache;
	    this->nGhost = nCache;
	    this->nAbsent= nCache;
	    HASH::resize(nCache + nGhost + nAbsent);
	    // Allocate the total possible amount of storage. If we change cache types we won't have to reallocate.
	    dataBase.resize(uCacheSizeInBytes / sizeof(uint64_t));
	    }
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
	uint64_t *p = static_cast<uint64_t*>(vp);
	if (p>&dataBase.front() && p<=&dataBase.back()) { // Might have been a fast, read-only grab. If so ignore it.
	    /* We will be given an element, but this needs to be turned into a cache line */
	    p = &dataBase[(p - &dataBase.front()) / (uLineSizeInWords+1) * (uLineSizeInWords+1) + 1];
	    uint64_t t = p[-1]-1;
	    assert((t^_ARC_MAGIC_) < 0x00000000ffffffff); // Not an element or too many unlocks
	    p[-1] = t;
	    }
	}

/*****************************************************************************\
* ARC inject
\*****************************************************************************/

template<typename... KEYS>
void *ARC<KEYS...>::inject(uint32_t uHash, uint32_t uId, const KEY& key) {
    int uIndex = uHash;
    auto key_size = sizeof(PAIR) == sizeof(CDB) ? 0 : sizeof(KEY);
    auto &Hash = HashChains[uHash&uHashMask];
    auto pEntry = key_size ? find_key(Hash,key) : find_key(Hash,uIndex,uId);
    if (pEntry) return nullptr; // Already have a cached value
    auto item = insert(Hash);
    auto &cdb = get<CDB>(*item);
    cdb.uId   = uId;                        // Remote processor that owns this entry
    cdb.uPage = uIndex;                     // Remote array index (simple) or hash id (advanced)
    if (key_size) get<KEY>(*item) = key;    // Advanced key (or nothing)
    return cdb.data; // Data is not initialized here
    }

template<typename... KEYS>
void *ARC<KEYS...>::inject(uint32_t uHash, uint32_t uId, const void *pKey) {
    return inject(uHash,uId,*static_cast<const KEY*>(pKey));
    }

/*****************************************************************************\
* ARC fetch
\*****************************************************************************/

template<typename... KEYS>
void *ARC<KEYS...>::fetch(uint32_t uIndex, uint32_t uId, const KEY& key, bool bLock,bool bModify,bool bVirtual) {
    auto uHash = uIndex;
    auto key_size = sizeof(PAIR) == sizeof(CDB) ? 0 : sizeof(KEY);
    uint32_t iInLine = 0;
    void *data;

    if (key_size==0) {                  // SIMPLE KEY
	iInLine = uIndex & nLineMask;   // Which element in the cache line
	uIndex >>= nLineBits;           // uIndex is now uLine (part of the key)
	uHash = hash::hash(uIndex,uId);
	}

    /* First check our own cache */
    auto &Hash = HashChains[uHash&uHashMask];
    auto pEntry = key_size ? find_key(Hash,key) : find_key(Hash,uIndex,uId);

    typename CDBL::iterator item;
    if (pEntry) {   // Page is reference by the cache
	item = CDBL::s_iterator_to(*pEntry);
	update(item,bLock,bModify);
	auto &cdb = get<CDB>(*item);
	if (cdb.data == nullptr) {
	    if (cdb.where()==A1) return nullptr; // Absent
	    data = helper->invokeRequest(cdb.uPage,cdb.uId,key_size,&key,bVirtual); // Request the element be fetched
	    cdb.data = replace(cdb.where()==B2);
	    helper->finishRequest(cdb.uPage,cdb.uId,key_size,&key,bVirtual,cdb.data,data);
	    }
	}
    else {          // Page is not in the cache
	if (key_size) uId = helper->getThread(uIndex,uId,key_size,&key);
	data = helper->invokeRequest(uIndex,uId,key_size,&key,bVirtual);
	if (!bModify && !bVirtual && data && key_size==0) return data; // Simple local case
	data = helper->finishRequest(uIndex,uId,key_size,&key,bVirtual,cacheLine.data(),data);
	if (data) {
	    item = insert(Hash);
	    auto &cdb = get<CDB>(*item);
	    std::memcpy(cdb.data,data,uLineSizeInWords*sizeof(uint64_t));
	    }
	else item = insert_absent(Hash); // Whelp, no remote element
	auto &cdb = get<CDB>(*item);
	cdb.uId   = uId;                        // Remote processor that owns this entry
	cdb.uPage = uIndex;                     // Remote array index (simple) or hash id (advanced)
	if (key_size) get<KEY>(*item) = key;    // Advanced key (or nothing)
	if (cdb.where()==A1) return nullptr; // Absent
	}
    auto &cdb = get<CDB>(*item);
    if (bLock) ++cdb.data[-1];   // Increase the lock count if requested
    if (bModify) cdb.setDirty(); // Mark page as dirty
    return reinterpret_cast<char*>(cdb.data) + uDataSizeInBytes*iInLine;
    }

template<typename... KEYS>
void *ARC<KEYS...>::fetch(uint32_t uHash, const KEYS&... keys, bool bLock,bool bModify,bool bVirtual) {
    return fetch(uHash,0,std::make_tuple(keys...),bLock,bModify,bVirtual);
    }

template<typename... KEYS>
void *ARC<KEYS...>::fetch(uint32_t uHash, const void *pKey, bool bLock,bool bModify,bool bVirtual) {
    return fetch(uHash,0,*static_cast<const KEY*>(pKey),bLock,bModify,bVirtual);
    }

template<typename... KEYS>
void *ARC<KEYS...>::fetch(uint32_t uIndex, uint32_t uId, bool bLock,bool bModify,bool bVirtual) {
    return fetch(uIndex,uId,KEY(),bLock,bModify,bVirtual);
    }


} // namespace mdl

#endif
