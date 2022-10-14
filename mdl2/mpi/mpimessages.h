#ifndef MPIMESSAGES_H
#define MPIMESSAGES_H
#include "mdlmessages.h"
#include "mpi.h"
namespace mdl {

class mdlMessageMPI : public mdlMessage {
protected:
    friend class mdlClass;
    friend class mpiClass;
public:
    virtual void action(class mpiClass *mpi) override = 0;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status);
};

// Used to hold a sequence of cache lines to send from the MPI thread to cores
class mdlMessageFlushToRank : public mdlMessageMPI, public FlushBuffer {
protected:
    friend class mdlClass;
    friend class mpiClass;
public:
    virtual void action(class mpiClass *mpi) override;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
};

// Send a small reply message with a single cache line
class mdlMessageCacheReply : public mdlMessageMPI, public FlushBuffer {
protected:
    friend class mdlClass;
    friend class mpiClass;
public:
    mdlMessageCacheReply(uint32_t nSize) : FlushBuffer(nSize,CacheMessageType::REPLY) {}
    virtual void action(class mpiClass *mpi) override;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
};

class mdlMessageCacheReceive : public mdlMessageMPI, public FlushBuffer {
protected:
    friend class mdlClass;
    friend class mpiClass;
public:
    mdlMessageCacheReceive(uint32_t nSize) : FlushBuffer(nSize,CacheMessageType::UNKNOWN) {}
    virtual void action(class mpiClass *mpi) override;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
};

class mdlMessageAlltoallv : public mdlMessageMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    int dataSize;
    void *sbuff, *rbuff;
    int *scount, *sdisps, *rcount, *rdisps;
public:
    virtual void action(class mpiClass *mpi) override;
    mdlMessageAlltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps);
};

class mdlMessageBarrierMPI : public mdlMessageMPI {
public:
    virtual void action(class mpiClass *mpi) override;
};

class mdlMessageBufferedMPI : public mdlMessageMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    void *buf;
    int count;
    int target;
    int tag;
public:
    virtual void action(class mpiClass *mpi) override = 0;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
    explicit mdlMessageBufferedMPI(void *buf, int count, int target, int tag);
    int getCount() {return count;}
};

class mdlMessageSend : public mdlMessageBufferedMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
public:
    virtual void action(class mpiClass *mpi) override;
    explicit mdlMessageSend(void *buf,int32_t count,int source, int tag);
};

class mdlMessageReceive : public mdlMessageBufferedMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    int iCoreFrom;
public:
    virtual void action(class mpiClass *mpi) override;
    explicit mdlMessageReceive(void *buf,int32_t count, int source, int tag,int iCoreFrom);
};

class mdlMessageReceiveReply : public mdlMessageReceive {
protected:
    friend class mdlClass;
    friend class mpiClass;
    ServiceHeader header;
public:
    virtual void action(class mpiClass *mpi) override;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
    explicit mdlMessageReceiveReply(void *buf,int32_t count, int rID, int iCoreFrom);
};

class mdlMessageReceiveRequest : public mdlMessageMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    ServiceHeader header;
    std::vector<char> Buffer;
public:
//    virtual void action(class mpiClass *mpi) override;
    explicit mdlMessageReceiveRequest(int32_t count=0);
};

class mdlMessageSendRequest : public mdlMessageBufferedMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    ServiceHeader header;
public:
    virtual void action(class mpiClass *mpi) override;
    explicit mdlMessageSendRequest(int32_t idFrom,int16_t sid,int target,void *buf=0,int32_t count=0);
};

class mdlMessageSendReply : public mdlMessageMPI {
protected:
    friend class mdlClass;
    friend class mpiClass;
    ServiceHeader header;
    std::vector<char> Buffer;
    int iThreadTo;
public:
    virtual void action(class mpiClass *mpi) override;
//    explicit mdlMessageSendReply(int32_t idFrom,int16_t replyTag, int16_t sid,int target,void *buf=0,int32_t count=0);
    explicit mdlMessageSendReply(int32_t count=0);
    mdlMessageSendReply &makeReply(int32_t idFrom,int16_t replyTag,int16_t sid,int target,int32_t count);
};

struct mdlCacheRequestData {
    CacheHeader header;         // Request header
    char key[MDL_MAX_KEY_SIZE]; // Optional advanced key
};

class mdlMessageCacheRequest : public mdlMessageBufferedMPI, protected mdlCacheRequestData {
protected:
    friend class mdlClass;
    friend class mpiClass;
    void *pLine = nullptr;
    uint32_t key_size = 0;
public:
    virtual void action(class mpiClass *mpi) override;
    virtual void finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) override;
    explicit mdlMessageCacheRequest(uint8_t cid, int32_t idFrom);
    explicit mdlMessageCacheRequest(uint8_t cid, int32_t idFrom, uint16_t nItems, int32_t idTo, int32_t iLine, void *pLine);
    mdlMessageCacheRequest &makeCacheRequest(uint16_t nItems, int32_t idTo, int32_t iLine, uint32_t size, const void *pKey, void *pLine);
};

} // namespace mdl
#endif
