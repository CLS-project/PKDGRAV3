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
#ifndef _FILE_OFFSET_BITS
    #define _FILE_OFFSET_BITS 64
#endif
#ifndef _LARGEFILE_SOURCE
    #define _LARGEFILE_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>
#include "iochunk.h"
#include "iomodule.h"

void io_chunk_read(const char *filename, void *buffer, uint64_t nBytes, uint64_t iOffset) {
    constexpr uint64_t chunk = 64*1024*1024; // Number of bytes to read at a time
    asyncFileInfo info;
    io_init(&info, 1, chunk, IO_AIO|IO_LIBAIO);
    auto fd = io_open(&info,filename);
    if (fd<0) {
        perror(filename);
        abort();
    }
    io_read(&info,buffer,nBytes,iOffset);
    io_close(&info);
    io_free(&info);
}

void io_chunk_write(const char *filename, const void *buffer, uint64_t nBytes, bool bAppend) {
    auto fd = open(filename,O_WRONLY | (bAppend?O_APPEND:O_CREAT | O_TRUNC),0666);
    if (fd<0) {
        perror(filename);
        abort();
    }
    constexpr uint64_t chunk = 1024*1024; // Number of bytes to read at a time
    auto p = static_cast<const char *>(buffer);
    for (auto i=0; i<nBytes; i+=chunk) {
        auto iSize = std::min(chunk,nBytes-i);
        auto nWrote = write(fd,p,iSize);
        if (nWrote != iSize) {
            fprintf(stderr,"Short write: ");
            perror(filename);
            abort();
        }
        p += nWrote;
    }
    close(fd);
}
