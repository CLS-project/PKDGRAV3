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

/* -*- Mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <cuda.h>

/*#include <nvToolsExt.h>*/
#ifdef USE_NVML
    #include <nvidia/gdk/nvml.h>
#endif
#include <signal.h>

#include "cudautil.h"

#include <assert.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
    #include <unistd.h>
#endif
#include <pthread.h>
#ifdef __APPLE__
    #include "pthread_barrier.h"
#endif
#ifdef HAVE_SYS_PARAM_H
    #include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif

static void *CUDA_malloc(size_t nBytes) {
#ifdef __linux__
    uint64_t nPageSize = sysconf(_SC_PAGESIZE);
#else
    uint64_t nPageSize = 512;
#endif
    void *blk;
#ifdef _MSC_VER
    blk = _aligned_malloc(nBytes, nPageSize);
#else
    if (posix_memalign(&blk, nPageSize, nBytes)) blk = NULL;
#endif
    char *p = reinterpret_cast<char *>(blk);
    char *e = p + nBytes;
    for (; p<e; p+= nPageSize) *p = 0;

    return blk;
}

static void CUDA_free(void *data) {
    free(data);
}

void CUDA_Abort(cudaError_t rc, const char *fname, const char *file, int line) {
    fprintf(stderr,"%s error %d in %s(%d)\n%s\n", fname, rc, file, line, cudaGetErrorString(rc));
    exit(1);
}

#ifdef _MSC_VER
double CUDA_getTime() {
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    return clock / 10000000.0;
}
#else
double CUDA_getTime() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
}
#endif

/*****************************************************************************\
*   CudaClient interface (new!)
\*****************************************************************************/

extern "C"
void *CudaClientInitialize(MDL vmdl) {
    auto mdl = reinterpret_cast<mdl::mdlClass *>(vmdl);
    return new CudaClient(*mdl);
}

CudaClient::CudaClient(mdl::mdlClass &mdl) : mdl(mdl), ewald(nullptr), pp(nullptr), pc(nullptr) {
    if (mdl.isCudaActive()) {
        freeEwald.enqueue(new MessageEwald(*this));
        freeEwald.enqueue(new MessageEwald(*this));
        freePP.enqueue(new MessagePP(freePP));
        freePP.enqueue(new MessagePP(freePP));
        freePC.enqueue(new MessagePC(freePC));
        freePC.enqueue(new MessagePC(freePC));
    }
}

extern "C"
void CudaClientFlush(void *vcudaClient) {
    auto cuda = reinterpret_cast<CudaClient *>(vcudaClient);
    cuda->flushCUDA();
}

void CudaClient::flushCUDA() {
    if (ewald) { mdl.enqueue(*ewald);        ewald = nullptr; }
    flush(pp);
    flush(pc);
}

cudaDataMessage::cudaDataMessage() {
    pHostBufIn = CUDA_malloc(requestBufferSize);
    pHostBufOut= CUDA_malloc(resultsBufferSize);
    CUDA_CHECK(cudaHostRegister,(pHostBufIn, requestBufferSize, cudaHostRegisterPortable));
    CUDA_CHECK(cudaHostRegister,(pHostBufOut,resultsBufferSize, cudaHostRegisterPortable));
}

cudaDataMessage::~cudaDataMessage() {
    CUDA_free(pHostBufIn);
    CUDA_free(pHostBufOut);
}
