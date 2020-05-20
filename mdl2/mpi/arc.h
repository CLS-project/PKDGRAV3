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
namespace hash {
using std::uint32_t;
using std::uint64_t;
using std::size_t;

class GHASH {
public:
    virtual ~GHASH() = default;
    virtual void   insert(uint32_t uHash, const void *pKey, void *data) = 0;
    virtual void  *lookup(uint32_t uHash, const void *pKey) = 0;
    virtual void   remove(uint32_t uHash, const void *pKey) = 0;
    virtual void   drop_all() = 0;
    virtual size_t key_length() = 0;
    virtual void   print_statistics() = 0;
    };

// Special case when there is no "key" in the hash table. Instead the index
// and id are used. The appropriate template specialization are provided in
// the module, hence we call abort() as this is never called.
struct no_key {
    bool operator==(const struct no_key &rhs) const {abort(); return false;}    
    };

    struct CDB {
	static const int whereBits = 3;
	uint64_t *data;                 // page's location in cache
	uint32_t uId    : 31-whereBits; // Processor that owns this
	uint32_t iWhere :    whereBits; // Which list we are on
	bool     bDirty :            1; // True if we need to flush
	uint32_t uPage;                 // page's ID number (or hash id for advanced keys)

	explicit CDB(uint64_t *data=nullptr) : uId(0), iWhere(0), data(data) {}
	uint32_t    where()   const {return iWhere;}
	uint32_t    getId()   const {return uId;}
	uint32_t    getPage() const {return uPage;}
	const char *getData() const {return reinterpret_cast<char*>(data);}
	bool        dirty()   const {return bDirty;}  // Data has been modified (data not null)
	bool        absent()  const {return dirty();} // ... or is absent (data is null)
	void setDirty() { this->bDirty = true; }
	void setClean() { this->bDirty = false; }
	};

template<typename KEY=no_key>
class HASH : public GHASH {
protected:
    struct HashChain { struct HashChain *next; };
    class ENTRY : public std::tuple<CDB,KEY>, public boost::intrusive::list_base_hook<>, public HashChain {};
    typedef boost::intrusive::list<ENTRY,boost::intrusive::constant_time_size<true> > CDBL;

protected:
    uint32_t uHashMask;
private:
    std::vector<ENTRY> entries;
    CDBL freeList;
protected:
    CDBL L[1<<CDB::whereBits];
    std::vector<HashChain> HashChains; // Collision chain (one for each hash value)
protected:
    // Search for the ENTRY matching the specified key, and return it if found,
    // or return nullptr if not found. Optionally remove it from the hash chain.
    ENTRY * find_key(HashChain &Hash,const KEY &key, bool bRemove=false) {
	for(auto pHash = &Hash; pHash->next != &Hash; pHash = pHash->next) {
	    auto pEntry = static_cast<ENTRY*>(pHash->next);
	    if (std::get<1>(*pEntry) == key) {
	    	if (bRemove) pHash->next = pHash->next->next;     	     
	    	return pEntry;
		}
	    }
	return nullptr;
	}
    ENTRY * find_key(HashChain &Hash,uint32_t uLine,uint32_t uId, bool bRemove=false) {
	for(auto pHash = &Hash; pHash->next != &Hash; pHash = pHash->next) {
	    auto pEntry = static_cast<ENTRY*>(pHash->next);
	    if (std::get<0>(*pEntry).uPage == uLine && (std::get<0>(*pEntry).getId()) == uId) {
	    	if (bRemove) pHash->next = pHash->next->next;     	     
	    	return pEntry;
		}
	    }
	return nullptr;
	}
    uint32_t hash(uint32_t uLine,uint32_t uId);
    typename CDBL::iterator move(uint32_t iTarget) {
	L[iTarget].splice(L[iTarget].end(),freeList,freeList.begin());
	std::get<0>(L[iTarget].back()).iWhere = iTarget;
	return std::prev(L[iTarget].end());
	}
    typename CDBL::iterator move(uint32_t iTarget, typename CDBL::iterator item) {
	auto &cdb = std::get<0>(*item);
	L[iTarget].splice(L[iTarget].end(),L[cdb.iWhere],item);
	cdb.iWhere = iTarget;
	return item;
	}
    typename CDBL::iterator free(typename CDBL::iterator item) {
	auto &cdb = std::get<0>(*item);
	freeList.splice(freeList.end(),L[cdb.iWhere],item);
	return item;
	}

public:
    virtual void  insert(uint32_t uHash, const void *pKey, void *data);
    virtual void *lookup(uint32_t uHash, const void *pKey);
    virtual void  remove(uint32_t uHash, const void *pKey);

    virtual void  drop_all() { clear(); }

    // Interface for generic key access and manipulation
    virtual size_t key_length() {return sizeof(KEY);}

public:
    class ARC *clone(int nMaxElements);

    virtual void print_statistics() {
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

    void clear() { // Empty the hash table (by moving all elements to the free list)
	for(auto &i : HashChains) i.next = &i; // Point to self = empty
	for(auto &list : L) freeList.splice(freeList.end(),list);
	}
    void resize(size_t count) {
	clear();
	uHashMask = count + count/2;    // Reserve extra for collision resolution
	uHashMask |= (uHashMask >> 1);  // Calculate the bitmask for the next higher
	uHashMask |= (uHashMask >> 2);  //   power of two so we can mask for an index.
	uHashMask |= (uHashMask >> 4);
	uHashMask |= (uHashMask >> 8);
	uHashMask |= (uHashMask >> 16);
	HashChains.resize(uHashMask+1);        // Collision chains for each possible hashed value
	for(auto &i : HashChains) i.next = &i; // Point to self = empty
	entries.resize(count); // Our complete list of ENTRY elements
	freeList.clear();      // Add them all to the free list
	for(auto &i : entries) freeList.push_back(i);
	}
    virtual ~HASH() = default;
    explicit HASH(int nMaxElements=0) : uHashMask(0) {
	resize(nMaxElements);
	}
    typename CDBL::iterator remove(const typename CDBL::iterator p);
    typename CDBL::iterator remove(uint32_t iTarget);

    };

template<typename KEY>
void *HASH<KEY>::lookup(uint32_t uHash, const void *vKey) {
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto &Hash2 = HashChains[uHash&uHashMask];
    auto pEntry = find_key(Hash2,key);
    if (pEntry) return std::get<0>(*pEntry).data;
    return nullptr;
    }

template<typename KEY>
void HASH<KEY>::insert(uint32_t uHash, const void *vKey, void *data) {
    assert(!freeList.empty()); // We cannot insert if the table is full
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto item = move(0);                          // Grab a new item from the free list
    auto &cdb = std::get<0>(*item);               // CDB (cache data block)
    cdb.uId = 0;                                  // Should set this to processor id probably
    cdb.uPage = uHash;                            // Page is the hash value for KEY types
    cdb.data = reinterpret_cast<uint64_t*>(data); // Points to the data for this element
    std::get<1>(*item) = key;                     // Record the key
    auto &Hash2 = HashChains[uHash&uHashMask];
    assert(find_key(Hash2,key)==nullptr);         // Duplicate keys are not allowed
    item->next = Hash2.next;
    Hash2.next = &*item;
    }

template<typename KEY>
void HASH<KEY>::remove(uint32_t uHash, const void *vKey) {
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto &Hash2 = HashChains[uHash&uHashMask];
    auto pEntry = find_key(Hash2,key,true);
    free(CDBL::s_iterator_to(*pEntry));
    }

} // namespace hash

// This class tells the ARC cache how to:
// 1. Request that a remote element be fetched
// 2. Finish the request and copy the data into ARC
// 3. Flush modified cache entries
class ARChelper {
public:
    virtual void invokeRequest(uint32_t uLine, uint32_t uId, void *pKey, bool bVirtual)             = 0;
    virtual void finishRequest(uint32_t uLine, uint32_t uId, void *pKey, bool bVirtual, void *data) = 0;
    virtual void flushElement( uint32_t uLine, uint32_t uId, void *pKey,          const void *data) = 0;
    };

template<typename ...Keys>
class ARC : private hash::HASH<Keys...> {
private:
    static const uint64_t _ARC_MAGIC_ = 0xa6c3a91c00000000;
    static const uint64_t _ARC_MASK_  = 0xa6c3a91cffffffff;
    enum WHERE {
	P1, T1, T2, B1, B2, // ARC cache lists
	};

private:
    ARChelper *helper;
    std::vector<uint64_t> dataBase; // Contains all cached data (continguous)
    typedef typename hash::HASH<Keys...>::CDBL CDBL;
    typedef typename hash::HASH<Keys...>::ENTRY ENTRY;

    typename std::vector<CDBL>::size_type target_T1;
    uint32_t nCache;
    uint32_t nLineBits, nLineMask;
    uint32_t uLineSizeInWords;
    uint32_t uLineSizeInBytes;
    uint32_t uDataSizeInBytes;
private:
    uint64_t *replace(WHERE iTarget, typename CDBL::iterator item);
    uint64_t *replace(bool iInB2=false);
    void evict(ENTRY &temp);
public:
    explicit ARC(ARChelper *helper=nullptr);
    virtual ~ARC();
    void initialize(ARChelper *helper, uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    void *fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual);
    void release(void *p);
    void RemoveAll();
    };

} // namespace mdl

#endif
