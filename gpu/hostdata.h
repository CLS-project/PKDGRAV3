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
#ifndef HOSTDATA_H
#define HOSTDATA_H
#include <cstdint>
#include <cstddef>
#include <cstdlib>
#ifdef __linux__
    #include <unistd.h>
#endif

#include "gravity/ilp.h"
#include "gravity/ilc.h"

namespace gpu {

/// Here we keep track of the input and output host buffers
class hostData {
protected:
    const std::size_t requestBufferSize;
    const std::size_t resultsBufferSize;
    void *pHostBufIn, *pHostBufOut;
public:
    explicit hostData(size_t requestBufferSize = 2*1024*1024,size_t resultsBufferSize = 2*1024*1024)
        : requestBufferSize(requestBufferSize), resultsBufferSize(resultsBufferSize) {
#ifdef __linux__
        std::uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
        std::uint64_t nPageSize = 4096;
#endif
        pHostBufIn  = std::aligned_alloc(nPageSize,requestBufferSize);
        pHostBufOut = std::aligned_alloc(nPageSize,resultsBufferSize);
    }
    virtual ~hostData() {
        std::free(pHostBufIn);
        std::free(pHostBufOut);
    }
};

} // namespace gpu
#endif
