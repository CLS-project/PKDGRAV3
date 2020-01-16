#ifndef ARC_H
#define ARC_H

#include <stdint.h>
#include <vector>
#include <list>
#include <forward_list>

class ARC {
protected:
    /*
    ** The bit values of these flags are NOT arbitrary!
    ** The implementation depends on their values being set this way.
    */
    static const uint32_t _P1_     = 0x10000000;
    static const uint32_t _T1_     = 0x00000000;
    static const uint32_t _B1_     = 0x20000000;
    static const uint32_t _T2_     = 0x40000000;
    static const uint32_t _B2_     = 0x60000000;
    static const uint32_t _WHERE_  = 0x70000000;
    static const uint32_t _IDMASK_ = 0x0fffffff;
    static const uint32_t _DIRTY_  = 0x80000000;

    static const uint64_t _ARC_MAGIC_ = 0xa6c3a91c00000000;
    static const uint64_t _ARC_MASK_  = 0xa6c3a91cffffffff;

    class CDB {
	friend class ARC;
    public:
	uint64_t *data; /* page's location in cache */
	uint32_t uId;   /* upper 4 bits encode ARC_where and dirty bit */
	uint32_t uPage; /* page's ID number */
    public:
	explicit CDB(uint64_t *data=0);
	bool dirty();
    };
    typedef std::list<CDB> CDBL;
    typedef std::forward_list<CDBL::iterator> CDBIL; // Memory efficient linked-list
private:
    std::vector<uint64_t> dataBase; // Contains all cached data (continguous)
    CDBIL HashFree;                 // List of free hash entries (added to HashChains)
    std::vector<CDBIL> HashChains;  // List of hash entries for each value
    CDBL T1;
    CDBL B1;
    CDBL T2;
    CDBL B2;
    CDBL ArcFree;
    uint32_t nLineBits, nLineMask;
    uint32_t uHashMask;
    uint32_t target_T1;
    uint32_t nCache;
    uint32_t uLineSizeInWords;
    uint32_t uLineSizeInBytes;
    uint32_t uDataSizeInBytes;
protected:
    std::list<CDB>::iterator remove_from_hash(const std::list<CDB>::iterator p);
    std::list<CDB>::iterator remove_from_hash(std::list<CDB> &list);
    uint64_t *replace(bool iInB2=false);
    uint32_t hash(uint32_t uLine,uint32_t uId);
protected:
    virtual void destage(CDB &temp);
    virtual void invokeRequest(uint32_t uLine, uint32_t uId, bool bVirtual) = 0;
    virtual void finishRequest(uint32_t uLine, uint32_t uId, bool bVirtual, void *data) = 0;
public:
    explicit ARC();
    explicit ARC(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    virtual ~ARC();
    void initialize(uint32_t uCacheSizeInBytes,uint32_t uLineSizeInBytes,uint32_t nLineBits=0);
    void *fetch(uint32_t uIndex, int uId, int bLock,int bModify,bool bVirtual);
    void release(void *p);
    void RemoveAll();
    };
#endif
