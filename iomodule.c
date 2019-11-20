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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#include "iomodule.h"
#include <assert.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdint.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
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

void io_init(asyncFileInfo *info, size_t nBuffers,size_t nBufferSize) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    int i, rc;
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

#ifdef HAVE_LIBAIO
    info->ctx = 0;
    rc = io_setup(info->nBuffers, &info->ctx);
    if (rc<0) { perror("io_setup"); abort(); }
#else
    memset(&info->cb,0,sizeof(info->cb));
    for(i=0; i<info->nBuffers; ++i) {
	info->pcb[i] = NULL;
	info->cb[i].aio_fildes = info->fd;
	info->cb[i].aio_offset = 0;
	info->cb[i].aio_buf = NULL;
	info->cb[i].aio_nbytes = 0;
	info->cb[i].aio_sigevent.sigev_notify = SIGEV_NONE;
	info->cb[i].aio_lio_opcode = LIO_NOP;
	}
#endif
    for(i=0; i<info->nBuffers; ++i) {
	void *vBuffer;
	if (posix_memalign(&vBuffer,info->nPageSize,info->nBufferSize)) vBuffer = NULL;
	assert(vBuffer!=NULL);
	info->pBuffer[i] = vBuffer;
	}
#endif
    }

void io_free(asyncFileInfo *info) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    int i;
    for(i=0; i<info->nBuffers; ++i) {
	free(info->pBuffer[i]);
	}
#endif
    }


int io_create(asyncFileInfo *info, const char *pathname) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    info->fd = open(pathname,O_DIRECT|O_CREAT|O_WRONLY|O_TRUNC,FILE_PROTECTION);
    info->iBuffer = 0; // Start in buffer zero
    info->iByte = 0; // Nothing in the buffer
#else
    info->fd = open(pathname,O_CREAT|O_WRONLY|O_TRUNC,FILE_PROTECTION);
#endif
    info->iFilePosition = 0;
    return info->fd;
    }

#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
static void queue_dio(asyncFileInfo *info,int i,int bWrite) {
    size_t nBytes = info->iByte > info->nBufferSize ? info->nBufferSize : info->iByte;
    size_t nWrite;
    int rc;

    ++info->nPending;
    assert(info->nPending <= info->nBuffers);

    /* Align buffer size for direct I/O. File will be truncated before closing if writing */
    nWrite = (nBytes+info->nPageSize-1) & ~(info->nPageSize-1);
#ifdef HAVE_LIBAIO
    struct iocb *pcb = &info->cb[i];
    if (bWrite) io_prep_pwrite(info->cb+i,info->fd,info->pBuffer[i],nWrite,info->iFilePosition);
    else        io_prep_pread(info->cb+i,info->fd,info->pBuffer[i],nWrite,info->iFilePosition);
    rc = io_submit(info->ctx,1,&pcb);
    if (rc<0) { perror("io_submit"); abort(); }
#else
    info->pcb[i] = info->cb + i;
    info->cb[i].aio_fildes = info->fd;
    info->cb[i].aio_buf = info->pBuffer[i];
    info->cb[i].aio_offset = info->iFilePosition;
    info->cb[i].aio_nbytes = nWrite;
    if (bWrite) rc = aio_write(&info->cb[i]);
    else rc = aio_read(&info->cb[i]);
    if (rc) { perror("aio_write/read"); abort(); }
#endif
    info->iFilePosition += nBytes;
    }

static int wait_complete(asyncFileInfo *info, int nWait) {
#ifdef HAVE_LIBAIO
    int nEvent = io_getevents(info->ctx,nWait,nWait,info->events,NULL);
    if (nEvent!=nWait) { perror("aio_getevents"); abort(); }
    info->nPending -= nWait;
    return info->events[0].obj - info->cb;
#else
    int iWait, rc, i;
    while(nWait) {
	rc = aio_suspend(info->pcb,info->nBuffers,NULL);
	if (rc) { perror("aio_suspend"); abort(); }
	for(i=0; i<info->nBuffers && nWait; ++i) {
	    char szError[100];
	    if (info->pcb[i] == NULL) continue;
	    rc = aio_error(info->pcb[i]);
	    if (rc == EINPROGRESS) continue;
	    else if (rc == 0) {
		iWait = i;
		info->pcb[i] = NULL;
		ssize_t nWritten = aio_return(&info->cb[i]);
		if (nWritten != info->cb[i].aio_nbytes) {
		    sprintf(szError,"errno=%d nBytes=%"PRIu64" nBytesWritten=%"PRIi64"\n",
			errno,(uint64_t)info->cb[i].aio_nbytes,(int64_t)nWritten);
		    perror(szError);
		    abort();
		    }
		--info->nPending;
		--nWait;
		}
	    else {
		errno = rc;
		sprintf(szError,"aio_error: rc=%d",rc);
		perror(szError);
		abort();
		}
	    }
	}
    return iWait;
#endif
    }
#endif

void io_write(asyncFileInfo *info, void *buf, size_t count) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    char *pBuf = buf;
    while(count) {
	size_t nBytes = info->nBufferSize - info->iByte;
	if (count < nBytes) nBytes = count;
	memcpy(info->pBuffer[info->iBuffer] + info->iByte,pBuf,nBytes);
	pBuf += nBytes;
	count -= nBytes;
	info->iByte += nBytes;
	if (info->iByte == info->nBufferSize) {
	    queue_dio(info,info->iBuffer,1);
	    info->iByte = 0;
	    if (info->nPending < info->nBuffers) info->iBuffer = info->nPending;
	    else info->iBuffer = wait_complete(info,1);
	    }
	}
#else
    if (write(info->fd,buf,count) != count) { perror("write"); abort(); }
    info->iFilePosition += count;
#endif
    }

void io_close(asyncFileInfo *info) {
#if defined(HAVE_LIBAIO) || defined(HAVE_AIO_H)
    if (info->iByte) queue_dio(info,info->iBuffer,1);
    if (info->nPending) wait_complete(info,info->nPending);
    assert(info->nPending==0);
    if (ftruncate(info->fd, info->iFilePosition)) perror("ftruncate");
#endif
    close(info->fd);
    }
