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

/*****************************************************************************\
*   CudaClient interface (new!)
\*****************************************************************************/

CudaClient::CudaClient( mdl::CUDA &cuda, mdl::gpu::Client &gpu) : cuda(cuda), gpu(gpu), ewald(nullptr), pp(nullptr), pc(nullptr), den(nullptr), denCorr(nullptr), sphForce(nullptr) {
    if (cuda.isActive()) {
        for (int i=0 ; i < 2 ; ++i) {
            freeEwald.enqueue(new MessageEwald(*this));
            freePP.enqueue(new MessagePP(freePP));
            freePC.enqueue(new MessagePC(freePC));
            freeDen.enqueue(new MessageDen(freeDen, &this->wps));
            freeDenCorr.enqueue(new MessageDenCorr(freeDenCorr));
            freeSPHForce.enqueue(new MessageSPHForce(freeSPHForce));
        }
    }
}

void CudaClient::flushCUDA() {
    if (ewald) { cuda.enqueue(*ewald,gpu);        ewald = nullptr; }
    flush(pp);
    flush(pc);
    flush(den);
    flush(denCorr);
    flush(sphForce);
}
