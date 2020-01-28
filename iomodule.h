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

#ifndef IOMODULE_H
#define IOMODULE_H
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif

#include <fcntl.h>
#include <sys/types.h>

#if defined(HAVE_LIBAIO)
#include <libaio.h>
#elif defined(HAVE_AIO_H)
#include <aio.h>
#endif
#define IO_MAX_ASYNC_COUNT 8
typedef struct {
#if defined(HAVE_LIBAIO)
    struct iocb cb[IO_MAX_ASYNC_COUNT];
    struct io_event events[IO_MAX_ASYNC_COUNT];
    io_context_t ctx;
#elif defined(HAVE_AIO_H)
    struct aiocb cb[IO_MAX_ASYNC_COUNT];
    struct aiocb const * pcb[IO_MAX_ASYNC_COUNT];
#endif
    char *pBuffer[IO_MAX_ASYNC_COUNT];
    off_t iFilePosition;   /* File position */
    size_t nBufferSize;    /* Total size of the buffer */
    size_t iByte;          /* Index into current buffer */
    int nBuffers;          /* Total number of buffers */
    int iBuffer;           /* Current buffer */
    int nPending;
    int nPageSize;
    int fd;
    } asyncFileInfo;


void io_init(asyncFileInfo *info, size_t nBuffers,size_t nBufferSize);
void io_free(asyncFileInfo *info);
int io_create(asyncFileInfo *info, const char *pathname);
void io_write(asyncFileInfo *info, void *buf, size_t count);
void io_close(asyncFileInfo *info);

#endif
