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

#ifndef CUDAPPPC_H
#define CUDAPPPC_H
#ifdef USE_CUDA
#include "cudautil.h"
#include "ilp.h"
#include "ilc.h"
#ifdef __cplusplus
extern "C" {
#endif
    int CUDA_queuePP(void *cudaCtx,workParticle *wp, ILPTILE tile, int bGravStep);
    int CUDA_queuePC(void *cudaCtx,workParticle *wp, ILCTILE tile, int bGravStep);
    void CUDA_sendWork(void *cudaCtx);
#ifdef __cplusplus
    }
#endif
#endif
#endif
