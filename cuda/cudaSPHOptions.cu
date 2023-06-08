/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2023 Joachim Stadel & Douglas Potter
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

#define OLD_CUDA

#ifdef HAVE_CONFIG_H
    #include "config.h"
#else
    #include "pkd_config.h"
#endif
#include <time.h>
#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif
#include <stdio.h>
#include "cudautil.h"
#include "SPH/SPHOptions.h"

__constant__ SPHOptionsGPU SPHoptions;

void CudaClient::setupSPHOptions(SPHOptionsGPU *const SPHoptions) {
    mdl::cudaMessageQueue wait;
    for (auto i=0; i<cuda.numDevices(); ++i) {
        auto m = new MessageSPHOptionsSetup(SPHoptions,i);
        cuda.enqueue(*m,wait);
    }
    for (auto i=0; i<cuda.numDevices(); ++i) {
        auto &m = wait.wait();
        delete &m;
    }
}

MessageSPHOptionsSetup::MessageSPHOptionsSetup(SPHOptionsGPU *const SPHoptions, int iDevice)
    : mdl::cudaMessage(iDevice), SPHoptionsIn(SPHoptions) {}

void MessageSPHOptionsSetup::launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) {
    CUDA_CHECK(cudaMemcpyToSymbolAsync, (SPHoptions, SPHoptionsIn, sizeof(SPHoptions), 0, cudaMemcpyHostToDevice, stream));
}