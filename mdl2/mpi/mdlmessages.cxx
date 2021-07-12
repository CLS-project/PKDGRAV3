#include "mdlmessages.h"
#include "mdl.h"
#include <string.h>

namespace mdl {
mdlMessage::mdlMessage() : basicMessage() {}

// Some messages don't need a result action, so this is the default
void mdlMessage::result(class mdlClass *mdl) {}

FlushBuffer::FlushBuffer(uint32_t nSize,CacheMessageType mid) : nBuffer(0),mid(mid),Buffer(nSize) {}

void*FlushBuffer::getBuffer(int nSize) {
    assert(nBuffer+nSize <= Buffer.size());
    if (nBuffer+nSize > Buffer.size()) return nullptr;
    auto data = &Buffer[nBuffer];
    nBuffer += nSize;
    return data;
    }

bool FlushBuffer::addBuffer(int nSize, const void *pData) {
    if (nBuffer+nSize > Buffer.size()) return false;
    if (pData) memcpy(&Buffer[nBuffer],pData,nSize);
    else memset(&Buffer[nBuffer],0,nSize);
    nBuffer += nSize;
    return true;
    }

bool FlushBuffer::addBuffer(uint8_t cid, int32_t idFrom, int32_t idTo, int32_t iLine, int nItems, int nSize, const void *pData) {
    if (!canBuffer(nSize + sizeof(CacheHeader))) return false;
    CacheHeader *ca = reinterpret_cast<CacheHeader*>(&Buffer.front() + nBuffer);
    char *pBuffer = reinterpret_cast<char *>(ca+1);
    ca->cid = cid;
    ca->mid = mid;
    ca->nItems = nItems;
    ca->idFrom = idFrom;
    ca->idTo = idTo;
    ca->iLine = iLine;
    if (nSize) memcpy(pBuffer,pData,nSize);
    nBuffer += nSize + sizeof(CacheHeader);
    return true;
    }

bool FlushBuffer::addBuffer(int nSize, const CacheHeader *pData) {
    if (!canBuffer(nSize)) return false;
    CacheHeader *ca = reinterpret_cast<CacheHeader*>(&Buffer.front() + nBuffer);
    memcpy(ca,pData,nSize+sizeof(CacheHeader));
    nBuffer += nSize + sizeof(CacheHeader);
    return true;
    }

// What to do when the MPI request has completed. Default is to send it back to the worker.
void mdlMessageMPI::finish(class mpiClass *mdl, const MPI_Status &status) { sendBack(); }

void mdlMessageBufferedMPI::finish(class mpiClass *mdl, const MPI_Status &status) {
    MPI_Get_count(&status, MPI_BYTE, &count); // Relevant for Recv() only
    target = status.MPI_SOURCE; // Relevant for Recv() only
    sendBack();
    }

// The "result" is processed on the worker core
void mdlMessageFlushToCore::result(class mdlClass *mdl) { mdl->MessageFlushToCore(this); }

// The "action" is to pass along the message to the MPI thread
void mdlMessageSTOP::action(class mpiClass *mpi)        { mpi->MessageSTOP(this); }
void mdlMessageBarrierMPI::action(class mpiClass *mpi)  { mpi->MessageBarrierMPI(this); }
void mdlMessageFlushFromCore::action(class mpiClass *mpi){ mpi->MessageFlushFromCore(this); }
void mdlMessageFlushToRank::action(class mpiClass *mpi) { mpi->MessageFlushToRank(this); }
void mdlMessageCacheReply::action(class mpiClass *mpi)  { mpi->MessageCacheReply(this); }
void mdlMessageCacheReceive::action(class mpiClass *mpi){ mpi->MessageCacheReceive(this); }
void mdlMessageCacheOpen::action(class mpiClass *mpi)   { mpi->MessageCacheOpen(this); }
void mdlMessageCacheClose::action(class mpiClass *mpi)  { mpi->MessageCacheClose(this); }
void mdlMessageCacheFlushOut::action(class mpiClass *mpi)  { mpi->MessageCacheFlushOut(this); }
void mdlMessageCacheFlushLocal::action(class mpiClass *mpi){ mpi->MessageCacheFlushLocal(this); }
void mdlMessageGridShare::action(class mpiClass *mpi)   { mpi->MessageGridShare(this); }
void mdlMessageDFT_R2C::action(class mpiClass *mpi)     { mpi->MessageDFT_R2C(this); }
void mdlMessageDFT_C2R::action(class mpiClass *mpi)     { mpi->MessageDFT_C2R(this); }
void mdlMessageFFT_Sizes::action(class mpiClass *mpi)   { mpi->MessageFFT_Sizes(this); }
void mdlMessageFFT_Plans::action(class mpiClass *mpi)   { mpi->MessageFFT_Plans(this); }
void mdlMessageAlltoallv::action(class mpiClass *mpi)   { mpi->MessageAlltoallv(this); }
void mdlMessageSend::action(class mpiClass *mpi)        { mpi->MessageSend(this); }
void mdlMessageReceive::action(class mpiClass *mpi)     { mpi->MessageReceive(this); }
void mdlMessageReceiveReply::action(class mpiClass *mpi){ mpi->MessageReceiveReply(this); }
void mdlMessageSendRequest::action(class mpiClass *mpi) { mpi->MessageSendRequest(this); }
void mdlMessageSendReply::action(class mpiClass *mpi)   { mpi->MessageSendReply(this); }
void mdlMessageCacheRequest::action(class mpiClass *mpi) { mpi->MessageCacheRequest(this); }

// For MPI messages, what to do when the request completes
void mdlMessageReceiveReply::finish(class mpiClass *mpi, const MPI_Status &status) {
    MPI_Get_count(&status, MPI_BYTE, &count); // Relevant for Recv() only
    count -= sizeof(header);
    mpi->FinishReceiveReply(this);
    }
void mdlMessageFlushToRank::finish(class mpiClass *mpi, const MPI_Status &status) {
    mpi->FinishFlushToRank(this);
    }
void mdlMessageCacheReply::finish(class mpiClass *mpi, const MPI_Status &status) {
    mpi->FinishCacheReply(this);
    }
void mdlMessageCacheReceive::finish(class mpiClass *mpi, const MPI_Status &status) {
    mpi->FinishCacheReceive(this,status);
    }

// Constructors
mdlMessageGridShare::mdlMessageGridShare(MDLGRID grid) : grid(grid) {}
mdlMessageDFT_R2C::mdlMessageDFT_R2C(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata)
    : fft(fft), data(data), kdata(kdata)   {}
mdlMessageDFT_C2R::mdlMessageDFT_C2R(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata)
    : fft(fft), data(data), kdata(kdata)   {}
mdlMessageFFT_Sizes::mdlMessageFFT_Sizes(int n1, int n2, int n3)
    : n1(n1), n2(n2), n3(n3) {}
mdlMessageFFT_Plans::mdlMessageFFT_Plans(int n1, int n2, int n3,FFTW3(real) *data,FFTW3(complex) *kdata)
    : mdlMessageFFT_Sizes(n1,n2,n3), data(data), kdata(kdata) {}
mdlMessageAlltoallv::mdlMessageAlltoallv(int dataSize,void *sbuff,int *scount,int *sdisps,void *rbuff,int *rcount,int *rdisps)
    : dataSize(dataSize), sbuff(sbuff), rbuff(rbuff), scount(scount), sdisps(sdisps), rcount(rcount), rdisps(rdisps) {}
mdlMessageBufferedMPI::mdlMessageBufferedMPI(void *buf, int count, MPI_Datatype datatype, int target, int tag)
    : buf(buf), count(count), target(target), tag(tag), datatype(datatype) {}
mdlMessageSend::mdlMessageSend(void *buf,int32_t count,MPI_Datatype datatype, int source, int tag)
    : mdlMessageBufferedMPI(buf,count,datatype,source,tag) {}
mdlMessageReceive::mdlMessageReceive(void *buf,int32_t count,MPI_Datatype datatype, int source, int tag, int iCoreFrom)
    : mdlMessageBufferedMPI(buf,count,datatype,source,tag), iCoreFrom(iCoreFrom) {}
mdlMessageReceiveReply::mdlMessageReceiveReply(void *buf,int32_t count, int rID, int iCoreFrom)
    : mdlMessageReceive(buf,count+sizeof(header),MPI_BYTE,rID,MDL_TAG_RPL,iCoreFrom) {}

mdlMessageSendRequest::mdlMessageSendRequest(int32_t idFrom,int16_t sid,int target,void *buf,int32_t count)
    : mdlMessageBufferedMPI(buf,count,MPI_BYTE,target,MDL_TAG_REQ) {
    header.idFrom = idFrom;
    header.replyTag = -1;
    header.sid = sid;
    header.nInBytes=count;
    header.nOutBytes = 0;
    }

mdlMessageReceiveRequest::mdlMessageReceiveRequest(int32_t count) : Buffer(count) {}

// mdlMessageSendReply::mdlMessageSendReply(int32_t idFrom,int16_t replyTag,int16_t sid,int target,void *buf,int32_t count)
//     : mdlMessageBufferedMPI(buf,count,MPI_BYTE,target,MDL_TAG_RPL) {
mdlMessageSendReply::mdlMessageSendReply(int32_t count) : Buffer(count) {}

mdlMessageSendReply & mdlMessageSendReply::makeReply(int32_t idFrom,int16_t replyTag,int16_t sid,int target,int32_t count) {
    iThreadTo = target;
    header.idFrom = idFrom;
    header.replyTag = replyTag;
    header.sid = sid;
    header.nInBytes=count;
    header.nOutBytes = 0;
    return *this;
    }

mdlMessageCacheRequest::mdlMessageCacheRequest(uint8_t cid, int32_t idFrom) 
    : mdlMessageBufferedMPI(&header,sizeof(header),MPI_BYTE,0,MDL_TAG_CACHECOM) {
    header.cid   = cid;
    header.mid   = CacheMessageType::REQUEST;
    header.nItems= 0;
    header.idFrom= idFrom;
    header.idTo  = 0;
    header.iLine = 0;
    }

mdlMessageCacheRequest::mdlMessageCacheRequest(uint8_t cid, int32_t idFrom, uint16_t nItems, int32_t idTo, int32_t iLine, void *pLine)
    : mdlMessageBufferedMPI(&header,sizeof(header),MPI_BYTE,idTo,MDL_TAG_CACHECOM), pLine(pLine) {
    header.cid   = cid;
    header.mid   = CacheMessageType::REQUEST;
    header.nItems= nItems;
    header.idFrom= idFrom;
    header.idTo  =  idTo;
    header.iLine = iLine;
    }

mdlMessageCacheRequest & mdlMessageCacheRequest::makeCacheRequest(uint16_t nItems, int32_t idTo, int32_t iLine, uint32_t size, const void *pKey, void *pLine) {
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

// Normally, when an MPI request finishes, we send it back to the requesting thread. The process is
// different for cache requests. We do nothing because the result is actually sent back, not the request.
// You would think that the "request" MPI send would complete before the response message is received,
// but this is NOT ALWAYS THE CASE. Care must be take if new/delete is used on this type of message.
void mdlMessageCacheRequest::finish(class mpiClass *mdl, const MPI_Status &status) {}

} // namespace mdl