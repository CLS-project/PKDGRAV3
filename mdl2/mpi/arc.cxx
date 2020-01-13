#include "arc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

static inline uint32_t murmur2(const uint32_t *key, int len) {
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    uint32_t h = 0xdeadbeef /*^ len : len will be the same */;
    while(len--) {
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
static inline uint32_t murmur3(const uint32_t* key, size_t len) {
    uint32_t h = 0xdeadbeef;
    do {
	uint32_t k = *key++;
	k *= 0xcc9e2d51;
	k = (k << 15) | (k >> 17);
	k *= 0x1b873593;
	h ^= k;
	h = (h << 13) | (h >> 19);
	h = h * 5 + 0xe6546b64;
	} while (--len);
    h ^= h >> 16;
    h *= 0x85ebca6b;
    h ^= h >> 13;
    h *= 0xc2b2ae35;
    h ^= h >> 16;
    return h;
    }

/*
** MurmurHash2, by Austin Appleby
** adapted for hashing 2 uint32_t variables for mdl2
*/
static inline uint32_t MurmurHash2(uint32_t a,uint32_t b) {
    /* 
    ** 'm' and 'r' are mixing constants generated offline.
    ** They're not really 'magic', they just happen to work well.
    */
    const uint32_t m = 0x5bd1e995;
    const int r = 24;
    uint32_t h = 0xdeadbeef;

    /* Mix the 2 32-bit words into the hash */
    a *= m; 
    b *= m; 
    a ^= a >> r; 
    b ^= b >> r; 
    a *= m; 	
    b *= m; 	
    /* now work on the hash */
    h ^= a;
    h *= m; 
    h ^= b;	
    /* Do a few final mixes of the hash to ensure the last few
    ** bytes are well-incorporated. */
    h ^= h >> 13;
    h *= m;
    h ^= h >> 15;
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

ARC::ARC()
    : nCache(0), uLineSizeInWords(0), uLineSizeInBytes(0), uDataSizeInBytes(0),
      nLineBits(0), nLineMask(0), uHashMask(0), target_T1(0) {}

ARC::ARC(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
    initialize(uCacheSizeInBytes,uLineSizeInBytes,nLineBits);
    }

void ARC::initialize(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits) {
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

    // Allocate the total possible amount of storage. If we change cache types we won't have to reallocate.
    dataBase.resize(uCacheSizeInBytes >> 3);

    // Calculate the size of the hash table
    uHashMask = swar32(3*nCache-1);
    HashFree.resize(nCache);
    HashChains.resize(uHashMask+1);

    target_T1 = nCache/2;   /* is this ok? */

    // Setup the CDB entries on the free list. The first "nCache" elements have data; the rest don't.
    ArcFree.clear();
    for (auto i=0;i<nCache;++i) ArcFree.emplace_back(&dataBase[i*(uLineSizeInWords+1)+1]);
    ArcFree.resize(2*nCache);
    }

ARC::~ARC() {
    }

void ARC::release(void *vp) {
    uint64_t *p = reinterpret_cast<uint64_t*>(vp);
    if (p>&dataBase.front() && p<=&dataBase.back()) { // Might have been a fast, read-only grab. If so ignore it.
	/* We will be given an element, but this needs to be turned into a cache line */
	p = &dataBase[(p - &dataBase.front()) / (uLineSizeInWords+1) * (uLineSizeInWords+1) + 1];
	uint64_t t = p[-1]-1;
	assert((t^_ARC_MAGIC_) < 0x00000000ffffffff); // Not an element or too many unlocks
	p[-1] = t;
	}
    }

ARC::CDB::CDB(uint64_t *data) : uId(0xdeadbeef), data(data) {}

bool ARC::CDB::dirty() { return (uId& _DIRTY_) != 0; }

uint32_t ARC::hash(uint32_t uLine,uint32_t uId) {
    return MurmurHash2(uLine,uId) & uHashMask;
    }

// Remove the CDB entry (given by a list iterator) from the hash table
ARC::CDBL::iterator ARC::remove_from_hash(const CDBL::iterator p) {
    uint32_t uHash = hash(p->uPage,p->uId&_IDMASK_);
    auto & Hash = HashChains[uHash];
    assert(!Hash.empty());
    // Move the matching hash table entry to the start of the collision chain, then move it to the free list
    std::partition(Hash.begin(), Hash.end(), [&](const CDBL::iterator & i) { return i == p; });
    assert(Hash.front() == p);
    HashFree.splice_after(HashFree.before_begin(),Hash,Hash.before_begin());
    return p;
    }

// Typically we want to remove the LRU (front) element from a list, in anticipation of reusing the entry
ARC::CDBL::iterator ARC::remove_from_hash(CDBL &list) {
    return remove_from_hash(list.begin());
    }

void ARC::destage(CDB &temp) {}

uint64_t *ARC::replace(bool bInB2) {
    uint64_t *data;
    uint32_t max = (target_T1 > 1)?(target_T1+(bInB2?0:1)):1;
    auto unlocked = [](CDB&i) {return i.data[-1] == _ARC_MAGIC_;}; // True if there are no locks
    CDBL::iterator tempX;
    if (T1.size() >= max) { // T1’s size exceeds target?
	tempX = std::find_if(T1.begin(),T1.end(),unlocked); // Grab the oldest unlocked entry
	if (tempX != T1.end()) goto replace_T1;
	tempX = std::find_if(T2.begin(),T2.end(),unlocked); // Try the same in T2 if all T1 entries are locked
	if (tempX != T2.end()) goto replace_T2;
	}
    else { // no: T1 is not too big
	tempX = std::find_if(T2.begin(),T2.end(),unlocked); // Try the same in T2 if all T1 entries are locked
	if (tempX != T2.end()) goto replace_T2;
	tempX = std::find_if(T1.begin(),T1.end(),unlocked); // Grab the oldest unlocked entry
	if (tempX != T1.end()) goto replace_T1;
	}
    fprintf(stderr,"ERROR: all ARC entries are locked, aborting\n");
    abort();
    if (0) {  /* using a Duff's device to handle the replacement */
    replace_T1:
	destage(*tempX);     /* if dirty, evict before overwrite */
	data = tempX->data;
	tempX->data = NULL; /*GHOST*/
	tempX->uId = (tempX->uId&_IDMASK_)|_B1_;  /* need to be careful here because it could have been _P1_ */
	B1.splice(B1.end(),T1,tempX); // Move element from T1 to B1
	}
    if (0) {
    replace_T2:
	destage(*tempX);     /* if dirty, evict before overwrite */
	data = tempX->data;
	tempX->data = NULL; /*GHOST*/
        tempX->uId |= _B2_;          /* note that fact */
	B2.splice(B2.end(),T2,tempX); // Move element from T2 to B2
	}
    return data;
    }

void ARC::RemoveAll() {
    for(auto &i : T1) { destage(i); } // Flush any dirty (modified) cache entries
    for(auto &i : T2) { destage(i); }
    for(auto &i : HashChains) { HashFree.splice_after(HashFree.before_begin(),i); } // Empty the hash table
    ArcFree.splice(ArcFree.begin(),T1); // Finally empty the lists making sure entries with data are at the start
    ArcFree.splice(ArcFree.begin(),T2);
    ArcFree.splice(ArcFree.end(),  B1); // B lists have no data so add them to the end
    ArcFree.splice(ArcFree.end(),  B2);
    }

void *ARC::fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual) {
    auto uLine = uIndex >> nLineBits;
    auto iInLine = uIndex & nLineMask;
    auto tuId = uId&_IDMASK_;
    uint32_t rat;
    bool inB2=false;

    /* First check our own cache */
    auto uHash = hash(uLine,tuId);
    CDBL::iterator tempX;
    auto &Hash = HashChains[uHash];
    auto iskey = [&](CDBL::iterator &i) {return i->uPage == uLine && (i->uId&_IDMASK_) == tuId;};
    auto match = std::find_if(Hash.begin(),Hash.end(),iskey);
    if (match != Hash.end()) {                       /* found in cache directory? */
	tempX = *match;
	switch (tempX->uId & _WHERE_) {                   /* yes, which list? */
	case _P1_:
	    tempX->uId = uId;     /* clears prefetch flag and sets WHERE = _T1_ (zero) and dirty bit */
	    T1.splice(T1.end(),T1,tempX);
	    goto cachehit;
	case _T1_:
	    tempX->uId |= _T2_ | uId;
	    T2.splice(T2.end(),T1,tempX);
	    goto cachehit;
	case _T2_:
	    tempX->uId |= uId;          /* if the dirty bit is now set we need to record this */
	    T2.splice(T2.end(),T2,tempX);
	cachehit:
	    if (bLock) {
		/*
		** We don't have to check if the lock counter rolls over, since it will increment the  
		** magic number first. This in turn will cause this page to be locked in a way that it 
		** cannot be unlocked without an error condition.
		*/
		++tempX->data[-1];       /* increase lock count */
		}
	    /*
	    ** Get me outa here.
	    */
	    return reinterpret_cast<char*>(tempX->data) + uDataSizeInBytes*iInLine;
	case _B1_:                            /* B1 hit: favor recency */
	    rat = B2.size()/B1.size();
	    if (rat < 1) rat = 1;
	    target_T1 += rat;
	    if (target_T1 > nCache) target_T1 = nCache;
	    /* adapt the target size */
	    T2.splice(T2.end(),B1,tempX);
	    goto doBcase;
	case _B2_:                            /* B2 hit: favor frequency */
	    rat = B1.size()/B2.size();
	    if (rat < 1) rat = 1;
	    if (rat > target_T1) target_T1 = 0;
	    else target_T1 = target_T1 - rat;
	    /* adapt the target size */
	    inB2=true;
	    T2.splice(T2.end(),B2,tempX);
	doBcase:
	    if (!bVirtual) invokeRequest(uLine,uId); // Request the element be fetched
	    tempX->data = replace(inB2);                                /* find a place to put new page */
	    tempX->uId = _T2_|uId;     /* temp->ARC_where = _T2_; and set the dirty bit for this page */
	    tempX->uPage = uLine;                          /* bookkeep */
	    finishRequest(uLine,uId,tempX->data,bVirtual);
	    break;
	    }
	}

    else {                                                              /* page is not in cache directory */
	/*
	** Can initiate the data request right here, and do the rest while waiting...
	*/
	if (!bVirtual) invokeRequest(uLine,uId);
	auto L1Length = T1.size() + B1.size();
	if (L1Length == nCache) {                                   /* B1 + T1 full? */
	    if (T1.size() < nCache) {                                           /* Still room in T1? */
		tempX = remove_from_hash(B1);        /* yes: take page off B1 */
		tempX->data = replace();                                /* find new place to put page */
		T1.splice(T1.end(),B1,tempX);
		}
	    else {                                                      /* no: B1 must be empty */
		tempX = remove_from_hash(T1);       /* take page off T1 */
		destage(*tempX);     /* if dirty, evict before overwrite */
		T1.splice(T1.end(),T1,tempX);
		}
	    }
	else {                                                          /* B1 + T1 have less than c pages */
	    uint32_t nInCache = T1.size() + T2.size() + B1.size() + B2.size();
	    if (nInCache >= nCache) {         /* cache full? */
		/* Yes, cache full: */
		if (nInCache == 2*nCache) {
		    /* directory is full: */
		    tempX = remove_from_hash(B2);
		    T1.splice(T1.end(),B2,tempX);
		    inB2=true;
		} else {                                                   /* cache directory not full, easy case */
		    T1.splice(T1.end(),ArcFree,ArcFree.begin());
		    assert(T1.back().data == NULL);
		}
		T1.back().data = replace(inB2);                                /* new place for page */
	    } else {                                                      /* cache not full, easy case */
		T1.splice(T1.end(),ArcFree,ArcFree.begin());
		assert(T1.back().data != NULL);               /* This CDB should have an unused page associated with it */
		T1.back().data[-1] = _ARC_MAGIC_; /* this also sets nLock to zero */
		}
	    }
	tempX = --T1.end();
	tempX->uId = uId;                  /* temp->dirty = dirty;  p->ARC_where = _T1_; as well! */
	tempX->uPage = uLine;
	finishRequest(uLine,uId,tempX->data,bVirtual);
	Hash.splice_after(Hash.before_begin(),HashFree,HashFree.before_begin());
	Hash.front() = tempX;
    }
    if (bLock) {
	/*
	** We don't have to check if the lock counter rolls over, since it will increment the  
	** magic number first. This in turn will cause this page to be locked in a way that it 
	** cannot be unlocked without an error condition.
	*/
	++tempX->data[-1];       /* increase lock count */
    }
    /* If we will modify the element, then it must eventually be flushed. */
    if (bModify) tempX->uId |= _DIRTY_;

    return reinterpret_cast<char*>(tempX->data) + uDataSizeInBytes*iInLine;
    }
