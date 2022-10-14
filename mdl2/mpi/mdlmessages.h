#ifndef MDLMESSAGES_H
#define MDLMESSAGES_H
#include <stdint.h>
#include "mdl_config.h"
#include "mdlfft.h"
#include "basicmessage.h"
#include <vector>
namespace mdl {
enum class CacheMessageType : uint8_t {
    REQUEST = 0,
    REPLY = 1,
    FLUSH = 2,
    UNKNOWN,
};

#define MDL_FLUSH_DATA_SIZE 32000
#define MDL_MAX_KEY_SIZE 128

struct ServiceHeader {
    int32_t idFrom;      /* Global thread ID */
    int16_t replyTag;    /* Which MPI tag to send the reply */
    int16_t sid;
    int32_t nInBytes;
    int32_t nOutBytes;
};

struct CacheHeader {
    uint8_t cid;
    CacheMessageType mid;
    uint16_t nItems;
    int32_t idFrom;
    int32_t idTo;
    int32_t iLine;
};

class BufferTarget {
protected:
    uint32_t iRankTo;
public:
    uint32_t getRankTo() {return iRankTo;}
    void setRankTo(uint32_t iRank) { iRankTo=iRank; }
};

class FlushBuffer {
protected:
    uint32_t nBuffer;
    CacheMessageType mid;
    std::vector<char> Buffer;
    uint32_t mContains; // Bitmask of which CacheMessageType are contained
public:
    explicit FlushBuffer(uint32_t nSize=MDL_FLUSH_DATA_SIZE,CacheMessageType mid=CacheMessageType::FLUSH);
    char *getBuffer() {return &Buffer.front();}
    uint32_t getCount() {return nBuffer;}
    bool isEmpty() {return nBuffer==0;}
    void emptyBuffer() {nBuffer=0; mContains=0;}
    bool canBuffer(int nSize) { return nBuffer+nSize+sizeof(ServiceHeader) <= Buffer.size(); }
    void *getBuffer(int nSize);
    bool addBuffer(int nSize, const void *pData);
    bool addBuffer(uint8_t cid, int32_t idFrom, int32_t idTo, int32_t iLine, int nItems=1, int nSize=0, const void *pData=0);
    bool addBuffer(CacheMessageType mid,uint8_t cid, int32_t idFrom, int32_t idTo, int32_t iLine, int nItems=1, int nSize=0, const void *pData=0);
    bool addBuffer(int nSize, const CacheHeader *pData);
    bool contains(CacheMessageType mid) {return (mContains & (1<<static_cast<int>(mid))) != 0; }
};

class mdlMessage : public basicMessage {
public:
    mdlMessage();
    virtual void action(class mpiClass *mpi) {sendBack();}
    virtual void result(class mdlClass *mdl);
};
struct mdlMessageQueue : public messageQueue<mdlMessage> {
    mdlMessageQueue() : messageQueue<mdlMessage>() {}
};

class mdlMessageVote : public mdlMessage {
    int iVote;
public:
    explicit mdlMessageVote(int iVote=0) : iVote(iVote) {}
    int vote() {return iVote;}
    int vote(int i) { return (iVote=i); }
};

// Used to hold a sequence of cache lines to send to the MPI thread for processing
class mdlMessageFlushFromCore : public mdlMessage, public FlushBuffer {
public:
    virtual void action(class mpiClass *mdl);
};

// Used to hold a sequence of cache lines to send from the MPI thread to cores
class mdlMessageFlushToCore : public mdlMessage, public FlushBuffer {
public:
    virtual void result(class mdlClass *mdl);
};

class mdlMessageSTOP : public mdlMessage {
public:
    virtual void action(class mpiClass *mdl);
};

class mdlMessageCacheOpen : public mdlMessage {
public:
    virtual void action(class mpiClass *mdl);
};

class mdlMessageCacheClose : public mdlMessage {
public:
    virtual void action(class mpiClass *mdl);
};

class mdlMessageCacheFlushOut : public mdlMessage {
public:
    virtual void action(class mpiClass *mdl);
};

class mdlMessageCacheFlushLocal : public mdlMessage {
public:
    virtual void action(class mpiClass *mdl);
};

class mdlMessageGridShare : public mdlMessage {
    friend class mdlClass;
    friend class mpiClass;
protected:
    MDLGRID grid;
public:
    virtual void action(class mpiClass *mdl);
    explicit mdlMessageGridShare(MDLGRID grid);
};

#ifdef MDL_FFTW
class mdlMessageDFT_R2C : public mdlMessage {
    friend class mdlClass;
    friend class mpiClass;
protected:
    MDLFFT fft;
    FFTW3(real) *data;
    FFTW3(complex) *kdata;
public:
    virtual void action(class mpiClass *mdl);
    explicit mdlMessageDFT_R2C(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata);
};

class mdlMessageDFT_C2R : public mdlMessage {
    friend class mdlClass;
    friend class mpiClass;
protected:
    MDLFFT fft;
    FFTW3(real) *data;
    FFTW3(complex) *kdata;
public:
    virtual void action(class mpiClass *mdl);
    explicit mdlMessageDFT_C2R(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata);
};

class mdlMessageFFT_Sizes : public mdlMessage {
    friend class mdlClass;
    friend class mpiClass;
protected: // Input fields
    int n1,n2,n3;
protected: // Output fields
    ptrdiff_t nz, sz, ny, sy, nLocal;
public:
    virtual void action(class mpiClass *mdl);
    explicit mdlMessageFFT_Sizes(int n1, int n2, int n3);
};

class mdlMessageFFT_Plans : public mdlMessageFFT_Sizes {
    friend class mdlClass;
    friend class mpiClass;
protected: // Input fields
    FFTW3(real) *data;
    FFTW3(complex) *kdata;
protected: // Output fields
    FFTW3(plan) fplan, iplan;
public:
    virtual void action(class mpiClass *mdl);
    explicit mdlMessageFFT_Plans(int n1, int n2, int n3,FFTW3(real) *data=0,FFTW3(complex) *kdata=0);
};
#endif

} // namespace mdl
#endif
