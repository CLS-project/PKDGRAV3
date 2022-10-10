#include "mdlmessages.h"
#include "mdl.h"
#include <string.h>

namespace mdl {
mdlMessage::mdlMessage() : basicMessage() {}

// Some messages don't need a result action, so this is the default
void mdlMessage::result(class mdlClass *mdl) {}

FlushBuffer::FlushBuffer(uint32_t nSize,CacheMessageType mid) : nBuffer(0),mid(mid),Buffer(nSize) {}

void *FlushBuffer::getBuffer(int nSize) {
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
    CacheHeader *ca = reinterpret_cast<CacheHeader *>(&Buffer.front() + nBuffer);
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
    CacheHeader *ca = reinterpret_cast<CacheHeader *>(&Buffer.front() + nBuffer);
    memcpy(ca,pData,nSize+sizeof(CacheHeader));
    nBuffer += nSize + sizeof(CacheHeader);
    return true;
}

// The "result" is processed on the worker core
void mdlMessageFlushToCore::result(class mdlClass *mdl) { mdl->MessageFlushToCore(this); }

// The "action" is to pass along the message to the MPI thread
void mdlMessageSTOP::action(class mpiClass *mpi)        { mpi->MessageSTOP(this); }
void mdlMessageFlushFromCore::action(class mpiClass *mpi) { mpi->MessageFlushFromCore(this); }
void mdlMessageCacheOpen::action(class mpiClass *mpi)   { mpi->MessageCacheOpen(this); }
void mdlMessageCacheClose::action(class mpiClass *mpi)  { mpi->MessageCacheClose(this); }
void mdlMessageCacheFlushOut::action(class mpiClass *mpi)  { mpi->MessageCacheFlushOut(this); }
void mdlMessageCacheFlushLocal::action(class mpiClass *mpi) { mpi->MessageCacheFlushLocal(this); }
void mdlMessageGridShare::action(class mpiClass *mpi)   { mpi->MessageGridShare(this); }
#ifdef MDL_FFTW
void mdlMessageDFT_R2C::action(class mpiClass *mpi)     { mpi->MessageDFT_R2C(this); }
void mdlMessageDFT_C2R::action(class mpiClass *mpi)     { mpi->MessageDFT_C2R(this); }
void mdlMessageFFT_Sizes::action(class mpiClass *mpi)   { mpi->MessageFFT_Sizes(this); }
void mdlMessageFFT_Plans::action(class mpiClass *mpi)   { mpi->MessageFFT_Plans(this); }
#endif

// Constructors
mdlMessageGridShare::mdlMessageGridShare(MDLGRID grid) : grid(grid) {}
#ifdef MDL_FFTW
mdlMessageDFT_R2C::mdlMessageDFT_R2C(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata)
    : fft(fft), data(data), kdata(kdata)   {}
mdlMessageDFT_C2R::mdlMessageDFT_C2R(MDLFFT fft, FFTW3(real) *data, FFTW3(complex) *kdata)
    : fft(fft), data(data), kdata(kdata)   {}
mdlMessageFFT_Sizes::mdlMessageFFT_Sizes(int n1, int n2, int n3)
    : n1(n1), n2(n2), n3(n3) {}
mdlMessageFFT_Plans::mdlMessageFFT_Plans(int n1, int n2, int n3,FFTW3(real) *data,FFTW3(complex) *kdata)
    : mdlMessageFFT_Sizes(n1,n2,n3), data(data), kdata(kdata) {}
#endif

} // namespace mdl
