#include "mpimessages.h"
#include "mdl.h"

namespace mdl {
// The "action" is to pass along the message to the MPI thread
void mdlMessageFlushToRank::action(class mpiClass *mpi) { mpi->MessageFlushToRank(this); }
void mdlMessageCacheReply::action(class mpiClass *mpi)  { mpi->MessageCacheReply(this); }
void mdlMessageCacheReceive::action(class mpiClass *mpi) { mpi->MessageCacheReceive(this); }
void mdlMessageAlltoallv::action(class mpiClass *mpi)   { mpi->MessageAlltoallv(this); }
void mdlMessageBarrierMPI::action(class mpiClass *mpi)  { mpi->MessageBarrierMPI(this); }
void mdlMessageSend::action(class mpiClass *mpi)        { mpi->MessageSend(this); }
void mdlMessageReceive::action(class mpiClass *mpi)     { mpi->MessageReceive(this); }
void mdlMessageReceiveReply::action(class mpiClass *mpi) { mpi->MessageReceiveReply(this); }
void mdlMessageSendRequest::action(class mpiClass *mpi) { mpi->MessageSendRequest(this); }
void mdlMessageSendReply::action(class mpiClass *mpi)   { mpi->MessageSendReply(this); }
void mdlMessageCacheRequest::action(class mpiClass *mpi) { mpi->MessageCacheRequest(this); }

// What to do when the MPI request has completed. Default is to send it back to the worker.
void mdlMessageMPI::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) { sendBack(); }

// For MPI messages, what to do when the request completes
void mdlMessageFlushToRank::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    mpi->FinishFlushToRank(this);
}
void mdlMessageCacheReply::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    mpi->FinishCacheReply(this);
}
void mdlMessageCacheReceive::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    mpi->FinishCacheReceive(this,request,status);
}
void mdlMessageBufferedMPI::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    int bytes, source;
    MPI_Get_count(&status, MPI_BYTE, &bytes); // Relevant for Recv() only
    source = status.MPI_SOURCE; // Relevant for Recv() only
    count = bytes;
    target = source;
    sendBack();
}
void mdlMessageReceiveReply::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    int bytes;
    MPI_Get_count(&status, MPI_BYTE, &bytes); // Relevant for Recv() only
    count = bytes - sizeof(header);
    mpi->FinishReceiveReply(this,request,status);
}
void mdlMessageCacheRequest::finish(class mpiClass *mpi, MPI_Request request, MPI_Status status) {
    mpi->FinishCacheRequest(this,request,status);
}

// Constructors
mdlMessageAlltoallv::mdlMessageAlltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps)
    : dataSize(dataSize), sbuff(sbuff), rbuff(rbuff), scount(scount), sdisps(sdisps), rcount(rcount), rdisps(rdisps) {}
mdlMessageBufferedMPI::mdlMessageBufferedMPI(void *buf, int count, int target, int tag)
    : buf(buf), count(count), target(target), tag(tag) {}
mdlMessageSend::mdlMessageSend(void *buf,int32_t count, int source, int tag)
    : mdlMessageBufferedMPI(buf,count,source,tag) {}
mdlMessageReceive::mdlMessageReceive(void *buf,int32_t count, int source, int tag, int iCoreFrom)
    : mdlMessageBufferedMPI(buf,count,source,tag), iCoreFrom(iCoreFrom) {}
mdlMessageReceiveReply::mdlMessageReceiveReply(void *buf,int32_t count, int rID, int iCoreFrom)
    : mdlMessageReceive(buf,count+sizeof(header),rID,MDL_TAG_RPL,iCoreFrom) {}
mdlMessageReceiveRequest::mdlMessageReceiveRequest(int32_t count) : Buffer(count) {}
mdlMessageSendRequest::mdlMessageSendRequest(int32_t idFrom,int16_t sid,int target,void *buf,int32_t count)
    : mdlMessageBufferedMPI(buf,count,target,MDL_TAG_REQ) {
    header.idFrom = idFrom;
    header.replyTag = -1;
    header.sid = sid;
    header.nInBytes=count;
    header.nOutBytes = 0;
}
mdlMessageSendReply::mdlMessageSendReply(int32_t count) : Buffer(count) {}

// Update message with correct "reply" parameters
mdlMessageSendReply &mdlMessageSendReply::makeReply(int32_t idFrom,int16_t replyTag,int16_t sid,int target,int32_t count) {
    iThreadTo = target;
    header.idFrom = idFrom;
    header.replyTag = replyTag;
    header.sid = sid;
    header.nInBytes=count;
    header.nOutBytes = 0;
    return *this;
}

// Cache Request constructors and update
mdlMessageCacheRequest::mdlMessageCacheRequest(uint8_t cid, int32_t idFrom)
    : mdlMessageBufferedMPI(&header,sizeof(header),0,MDL_TAG_CACHECOM) {
    header.cid   = cid;
    header.mid   = CacheMessageType::REQUEST;
    header.nItems= 0;
    header.idFrom= idFrom;
    header.idTo  = 0;
    header.iLine = 0;
}
mdlMessageCacheRequest::mdlMessageCacheRequest(uint8_t cid, int32_t idFrom, uint16_t nItems, int32_t idTo, int32_t iLine, void *pLine)
    : mdlMessageBufferedMPI(&header,sizeof(header),idTo,MDL_TAG_CACHECOM), pLine(pLine) {
    header.cid   = cid;
    header.mid   = CacheMessageType::REQUEST;
    header.nItems= nItems;
    header.idFrom= idFrom;
    header.idTo  =  idTo;
    header.iLine = iLine;
}

mdlMessageCacheRequest &mdlMessageCacheRequest::makeCacheRequest(uint16_t nItems, int32_t idTo, int32_t iLine, uint32_t size, const void *pKey, void *pLine) {
    static_assert(offsetof(mdlCacheRequestData, header) + sizeof(header) == offsetof(mdlCacheRequestData, key),
                  "The header and key are not adjacent in memory. This shouldn't happen.");
    header.nItems= nItems;
    header.idTo  =  idTo;
    header.iLine = iLine;
    this->pLine = pLine;
    key_size = size;
    if (size) {
        assert(size <= sizeof(key));
        if (size <= sizeof(key)) memcpy(key,pKey,size);
        else abort();
    }
    return * this;
}

} // namespace mdl