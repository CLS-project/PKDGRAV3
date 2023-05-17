/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
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

/*
** This does DIRECT I/O to avoid swamping memory with I/O buffers. Ideally this
** would not be required, but buffered I/O causes it to CRASH during I/O on the Cray
** (or any Linux system probably) when the amount of available memory is quite low.
** The test was performed on Piz Daint and Piz Dora and the memory statistics at
** the time were:
** Piz Dora (Cray XC40):
**   free memory (GB): min=   2.522 @ 1195 avg=   2.797 of  2024 std-dev=   0.257
**      resident size: max=  57.671 @  129 avg=  57.415 of  2024 std-dev=   0.256
** Piz Daint (Cray XC30):
**   free memory (GB): min=   1.599 @21174 avg=   1.870 of 34300 std-dev=   0.029
**      resident size: max=  27.921 @17709 avg=  27.870 of 34300 std-dev=   0.013
*/

#include "pkd_config.h"
#include "iomodule.h"
#include <assert.h>
#ifdef HAVE_UNISTD_H
    #include <unistd.h>
#endif
#include <stdint.h>
#include <cinttypes>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <algorithm>
#ifdef _MSC_VER
    #define FILE_PROTECTION (_S_IREAD | _S_IWRITE)
    typedef int ssize_t;
    #define open _open
    #define write _write
    #define close _close
#else
    #define FILE_PROTECTION (S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)
#endif
#ifndef __linux__
    #define O_DIRECT 0
#endif

// Call this function to initialize the direct I/O context.
// The buffer size must be a multiple of the page size, or zero
// in which case the reader or writer must provide buffers that are
// properly aligned and will persist until the I/O completes.
void io_init(asyncFileInfo *info, size_t nBuffers,size_t nBufferSize,int method) {
    info->method = IO_REGULAR; // Default, but maybe we can do better
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    int i;
#ifdef HAVE_UNISTD_H
    info->nPageSize = sysconf(_SC_PAGESIZE);
#else
    info->nPageSize = 4096; /* A conservative guess */
#endif
    if (nBuffers > IO_MAX_ASYNC_COUNT) nBuffers = IO_MAX_ASYNC_COUNT;
    info->nBufferSize = nBufferSize;
    info->nBuffers = nBuffers;
    info->iBuffer = 0;
    info->nPending = 0;
    info->bWrite = 0;
#ifdef HAVE_AIO_H
    if ((method & IO_AIO) && info->method==IO_REGULAR) {
        info->method = IO_AIO;
        memset(&info->io.aio.cb,0,sizeof(info->io.aio.cb));
        for (i=0; i<info->nBuffers; ++i) {
            info->io.aio.pcb[i] = NULL;
            info->io.aio.cb[i].aio_fildes = info->fd;
            info->io.aio.cb[i].aio_offset = 0;
            info->io.aio.cb[i].aio_buf = NULL;
            info->io.aio.cb[i].aio_nbytes = 0;
            info->io.aio.cb[i].aio_sigevent.sigev_notify = SIGEV_NONE;
            info->io.aio.cb[i].aio_lio_opcode = LIO_NOP;
        }
    }
#endif
#ifdef HAVE_LIBAIO
    if ((method & IO_AIO) && info->method==IO_REGULAR) {
        info->method = IO_LIBAIO;
        info->io.libaio.ctx = 0;
        int rc = io_setup(info->nBuffers, &info->io.libaio.ctx);
        if (rc<0) { perror("io_setup"); abort(); }
    }
#endif
    if (info->method&(IO_AIO|IO_LIBAIO)) {
        for (i=0; i<info->nBuffers; ++i) {
            void *vBuffer;
            if (info->nBufferSize>0) {
                if (posix_memalign(&vBuffer,info->nPageSize,info->nBufferSize)) vBuffer = NULL;
                assert(vBuffer!=NULL);
            }
            else vBuffer = NULL;
            info->pBuffer[i] = reinterpret_cast<char *>(vBuffer);
        }
    }
#endif
}

void io_free(asyncFileInfo *info) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->nBufferSize>0 && (info->method&(IO_AIO|IO_LIBAIO))) {
        int i;
        for (i=0; i<info->nBuffers; ++i) {
            free(info->pBuffer[i]);
        }
    }
#endif
}

int io_create(asyncFileInfo *info, const char *pathname) {
    int flags = O_CREAT|O_WRONLY|O_TRUNC;
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) flags |= O_DIRECT;
#endif
    info->fd = open(pathname,flags,FILE_PROTECTION);
    info->iBuffer = 0; // Start in buffer zero
    info->iByte = 0; // Nothing in the buffer
    info->iFilePosition = 0;
    info->bWrite = 1;
    return info->fd;
}

int io_open(asyncFileInfo *info, const char *pathname) {
    int flags = O_RDONLY;
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) flags |= O_DIRECT;
#endif
    info->fd = open(pathname,flags,FILE_PROTECTION);
    info->iBuffer = 0; // Start in buffer zero
    info->iByte = 0; // Nothing in the buffer
    info->iFilePosition = 0;
    info->bWrite = 0;
    return info->fd;
}

#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
static void queue_dio(asyncFileInfo *info,int i,int bWrite) {
    size_t nBytes = info->iByte;
    size_t nTransfer;
    int rc;

    ++info->nPending;
    assert(info->nPending <= info->nBuffers);

    /* Align buffer size for direct I/O. File will be truncated before closing if writing */
    nTransfer = (nBytes+info->nPageSize-1) & ~(info->nPageSize-1);
    info->nExpected[i] = bWrite ? nTransfer : nBytes; // We need to ask for a page size, but could *read* less
    info->iOffset[i] = info->iFilePosition;
#ifdef HAVE_AIO_H
    if (info->method == IO_AIO) {
        info->io.aio.pcb[i] = info->io.aio.cb + i;
        info->io.aio.cb[i].aio_fildes = info->fd;
        info->io.aio.cb[i].aio_buf = info->pBuffer[i];
        info->io.aio.cb[i].aio_offset = info->iFilePosition;
        info->io.aio.cb[i].aio_nbytes = nTransfer;
        if (bWrite) rc = aio_write(&info->io.aio.cb[i]);
        else rc = aio_read(&info->io.aio.cb[i]);
        if (rc) { perror("aio_write/read"); abort(); }
    }
#endif
#ifdef HAVE_LIBAIO
    if (info->method == IO_LIBAIO) {
        struct iocb *pcb = &info->io.libaio.cb[i];
        if (bWrite) io_prep_pwrite(info->io.libaio.cb+i,info->fd,info->pBuffer[i],nTransfer,info->iFilePosition);
        else        io_prep_pread(info->io.libaio.cb+i,info->fd,info->pBuffer[i],nTransfer,info->iFilePosition);
        rc = io_submit(info->io.libaio.ctx,1,&pcb);
        if (rc<0) { perror("io_submit"); abort(); }
    }
#endif
    info->iFilePosition += nBytes;
}

static int wait_complete(asyncFileInfo *info, int nWait) {
#ifdef HAVE_AIO_H
    if (info->method == IO_AIO) {
        int iWait, rc, i;
        while (nWait) {
            rc = aio_suspend(info->io.aio.pcb,info->nBuffers,NULL);
            if (rc) { perror("aio_suspend"); abort(); }
            for (i=0; i<info->nBuffers && nWait; ++i) {
                char szError[100];
                if (info->io.aio.pcb[i] == NULL) continue;
                rc = aio_error(info->io.aio.pcb[i]);
                if (rc == EINPROGRESS) continue;
                else if (rc == 0) {
                    iWait = i;
                    info->io.aio.pcb[i] = NULL;
                    ssize_t nBytesTransferred = aio_return(&info->io.aio.cb[i]);
                    if (nBytesTransferred < info->nExpected[i]) {
                        snprintf(szError,sizeof(szError),"errno=%d nBytesExpected=%" PRIu64 " nBytesTransferred=%" PRIi64 "\n",
                                 errno,(uint64_t)info->nExpected[i],(int64_t)nBytesTransferred);
                        perror(szError);
                        abort();
                    }
                    --info->nPending;
                    --nWait;
                }
                else {
                    errno = rc;
                    snprintf(szError,sizeof(szError),"aio_error: rc=%d",rc);
                    perror(szError);
                    abort();
                }
            }
        }
        return iWait;
    }
#endif
#ifdef HAVE_LIBAIO
    if (info->method == IO_LIBAIO) {
        int nEvent = io_getevents(info->io.libaio.ctx,nWait,nWait,info->io.libaio.events,NULL);
        if (nEvent!=nWait) { perror("aio_getevents"); abort(); }
        info->nPending -= nWait;
        return info->io.libaio.events[0].obj - info->io.libaio.cb;
    }
#endif
    assert(0);
}
#endif

void io_write(asyncFileInfo *info, void *buf, size_t count) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) {
        // It is the responsibility of the caller to provide a buffer that
        // is properly aligned, and that will persist for the duration of the I/O.
        if (info->nBufferSize==0) {
            info->pBuffer[info->iBuffer] = reinterpret_cast<char *>(buf);
            info->iByte = count;
            queue_dio(info,info->iBuffer,info->bWrite);
            info->iByte = 0;
            if (info->nPending < info->nBuffers) info->iBuffer = info->nPending;
            else info->iBuffer = wait_complete(info,1);
            return;
        }
        auto pBuf = reinterpret_cast<char *>(buf);
        while (count) {
            size_t nBytes = info->nBufferSize - info->iByte;
            if (count < nBytes) nBytes = count;
            memcpy(info->pBuffer[info->iBuffer] + info->iByte,pBuf,nBytes);
            pBuf += nBytes;
            count -= nBytes;
            info->iByte += nBytes;
            if (info->iByte == info->nBufferSize) {
                queue_dio(info,info->iBuffer,info->bWrite);
                info->iByte = 0;
                if (info->nPending < info->nBuffers) info->iBuffer = info->nPending;
                else info->iBuffer = wait_complete(info,1);
            }
        }
        return;
    }
#endif
    if (write(info->fd,buf,count) != count) { perror("write"); abort(); }
    info->iFilePosition += count;
}

void io_read(asyncFileInfo *info, void *buf, size_t count) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) {
        assert(info->nBufferSize==0);
        if (info->nPending < info->nBuffers) info->iBuffer = info->nPending;
        else info->iBuffer = wait_complete(info,1);
        info->pBuffer[info->iBuffer] = reinterpret_cast<char *>(buf);
        info->iByte = count;
        queue_dio(info,info->iBuffer,0);
        info->iByte = 0;
        return;
    }
#endif
    if (read(info->fd,buf,count) != count) { perror("read"); abort(); }
    info->iFilePosition += count;
}

static void finish_read(asyncFileInfo *info, char *p, size_t offset, int i) {
    auto src = info->pBuffer[i];
    if (info->iOffset[i] >= offset) {
        memcpy(p+info->iOffset[i]-offset,src,info->nExpected[i]);
    }
    else {
        auto nSkip = offset - info->iOffset[i];
        memcpy(p,src+nSkip,info->nExpected[i]-nSkip);
    }
}


void io_read(asyncFileInfo *info, void *buf, size_t count,size_t offset) {
    auto p = static_cast<char *>(buf);
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) {
        assert(info->nBufferSize>0 && info->nBuffers>0);
        auto nExtra = offset % info->nPageSize;
        info->iFilePosition = offset - nExtra;
        count += nExtra;
        while (count) {
            auto iSize = std::min(count,info->nBufferSize);
            if (info->nPending < info->nBuffers) info->iBuffer = info->nPending;
            else {
                auto i = info->iBuffer = wait_complete(info,1);
                finish_read(info,p,offset,i);
            }
            info->iByte = iSize;
            queue_dio(info,info->iBuffer,0);
            info->iByte = 0;
            count -= iSize;
        }
        while (info->nPending) finish_read(info,p,offset,wait_complete(info,1));
        return;
    }
#endif
    info->iFilePosition = offset;
    if (lseek(info->fd, offset, SEEK_SET) != offset) {
        perror("Seek error");
        abort();
    }
    while (count) {
        auto iSize = std::min(count,info->nBufferSize);
        if (read(info->fd,p,iSize) != iSize) { perror("read"); abort(); }
        info->iFilePosition += iSize;
        count -= iSize;
        p += iSize;
    }
}

void io_close(asyncFileInfo *info) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->method&(IO_AIO|IO_LIBAIO)) {
        if (info->iByte) queue_dio(info,info->iBuffer,info->bWrite);
        if (info->nPending) wait_complete(info,info->nPending);
        assert(info->nPending==0);
        if (info->bWrite && ftruncate(info->fd, info->iFilePosition)) perror("ftruncate");
    }
#endif
    close(info->fd);
}
