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

#include <cstdint>
#include <vector>
#include <list>
#include <tuple>
#include <forward_list>

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
    typedef std::tuple<CDB,KEY> ENTRY;
    typedef std::list<ENTRY> CDBL;
    typedef std::forward_list<typename CDBL::iterator> CDBIL; // Memory efficient linked-list

private:
    uint32_t uHashMask;
    CDBL freeList;
protected:
    std::vector<CDBL> L;            // All of our lists: P1, T1, etc.
    std::vector<CDBIL> HashChains;  // List of hash entries for each value
    CDBIL HashFree;                 // List of free hash entries (added to HashChains)
protected:
    typename CDBIL::iterator find_key(CDBIL &Hash,const KEY &key) {
	auto iskey = [key](const typename CDBL::iterator &i) {return std::get<1>(*i) == key;};
	auto match = std::find_if(Hash.begin(),Hash.end(),iskey);
	return match;    
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

    void clear() { // Empty the hash table (by moving all elements to the free list)
	for(auto &i : HashChains)
	    HashFree.splice_after(HashFree.before_begin(),i);
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
	HashFree.resize(count);         // Preallocate the maximum number of entries
	HashChains.resize(uHashMask+1); // Collision chains for each possible hashed value
	freeList.resize(count);
	}
    virtual ~HASH() = default;
    explicit HASH(int nMaxElements=0) : uHashMask(0) {
	L.resize(1<<CDB::whereBits); // Maximum number of lists
	resize(nMaxElements);
	}
    typename CDBL::iterator remove(const typename CDBL::iterator p);
    typename CDBL::iterator remove(uint32_t iTarget);
    };

template<typename KEY>
void *HASH<KEY>::lookup(uint32_t uHash, const void *vKey) {
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto &Hash = HashChains[uHash&uHashMask];     // Hash collision chain
    auto match = find_key(Hash,key);              // Look for our key here
    if (match != Hash.end())
    	return std::get<0>(**match).data;         // Return the data pointer
    else return nullptr;                          // .. or null if not found
    }

template<typename KEY>
void HASH<KEY>::insert(uint32_t uHash, const void *vKey, void *data) {
    assert(!freeList.empty()); // We cannot insert if the table is full
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto &Hash = HashChains[uHash&uHashMask];     // Hash collision chain
    auto match = find_key(Hash,key);              // Look for our key here
    assert(match == Hash.end());                  // Duplicate keys are not allowed
    auto item = move(0);                          // Grab a new item from the free list
    auto &cdb = std::get<0>(*item);               // CDB (cache data block)
    cdb.uId = 0;                                  // Should set this to processor id probably
    cdb.uPage = uHash;                            // Page is the hash value for KEY types
    cdb.data = reinterpret_cast<uint64_t*>(data); // Points to the data for this element
    std::get<1>(*item) = key;                     // Record the key
    assert(!HashFree.empty());                    // Add to the collision chain
    Hash.splice_after(Hash.before_begin(),HashFree,HashFree.before_begin());
    Hash.front() = item;
    }

template<typename KEY>
void HASH<KEY>::remove(uint32_t uHash, const void *vKey) {
    const KEY &key = * reinterpret_cast<const KEY*>(vKey);
    auto &Hash = HashChains[uHash&uHashMask];     // Hash collision chain
    auto match = find_key(Hash,key);              // Look for our key here
    assert(match != Hash.end());                  // Removing a missing key is INVALID
    if (match != Hash.end()) {
	// Move the matching hash table entry to the start of the collision chain, then move it to the free list
	std::partition(Hash.begin(), Hash.end(), [match](const typename CDBL::iterator & i) { return i == *match; });
	HashFree.splice_after(HashFree.before_begin(),Hash,Hash.before_begin());
	free(*match); // Return the entry to the free list
	}
    }

} // namespace hash

class ARChelper {
protected:
    virtual void flushElement2( uint32_t uLine, uint32_t uId, void *pKey,          const void *data) = 0;
    virtual void invokeRequest2(uint32_t uLine, uint32_t uId, void *pKey, bool bVirtual)             = 0;
    virtual void finishRequest2(uint32_t uLine, uint32_t uId, void *pKey, bool bVirtual, void *data) = 0;
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
protected:
    virtual void destage(const char *data,uint32_t uIndex,uint32_t uId);
    virtual void invokeRequest(uint32_t uLine, uint32_t uId, bool bVirtual) = 0;
    virtual void finishRequest(uint32_t uLine, uint32_t uId, bool bVirtual, void *data) = 0;
public:
    explicit ARC();
    explicit ARC(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    virtual ~ARC();
    void initialize(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    void *fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual);
    void lock(void *p);
    void release(void *p);
    void RemoveAll();
    };

} // namespace mdl

#endif
