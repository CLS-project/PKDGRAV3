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
#include <unistd.h>
#include "iostaggered.h"
#include "iomodule.h"

void io_staggered_read(const char* filename, void* buffer1, void* buffer2, uint64_t nElements, uint64_t nBytes1, uint64_t nBytes2, uint64_t iOffset) {
    asyncFileInfo info;
    io_init(&info, 1, 0, IO_REGULAR);
    auto fd = io_open(&info, filename);
    if (fd < 0) {
        perror(filename);
        abort();
    }
    for (auto i = 0; i < nElements; i++) {
        io_read(&info, buffer1, nBytes1, iOffset);
        buffer1 = static_cast<void*>(static_cast<char*>(buffer1) + nBytes1);
        iOffset += nBytes1;

        io_read(&info, buffer2, nBytes2, iOffset);
        iOffset += nBytes2;
        buffer2 = static_cast<void*>(static_cast<char*>(buffer2) + nBytes2);
    }
    io_close(&info);
    io_free(&info);
}

void io_staggered_write(const char* filename, void* buffer1, void* buffer2, uint64_t nElements, uint64_t nBytes1, uint64_t nBytes2, bool bAppend) {
    auto fd = open(filename, O_WRONLY | (bAppend ? O_APPEND : O_CREAT | O_TRUNC), 0666);
    if (fd < 0) {
        perror(filename);
        abort();
    }
    auto p1 = static_cast<char*>(buffer1);
    auto p2 = static_cast<char*>(buffer2);
    for (auto i = 0; i < nElements; i++) {
        auto nWrote1 = write(fd, p1, nBytes1);
        if (nWrote1 != nBytes1) {
            fprintf(stderr, "Short write: ");
            perror(filename);
            abort();
        }
        p1 += nWrote1;
        auto nWrote2 = write(fd, p2, nBytes2);
        if (nWrote2 != nBytes2) {
            fprintf(stderr, "Short write: ");
            perror(filename);
            abort();
        }
        p2 += nWrote2;
    }
    close(fd);
}
