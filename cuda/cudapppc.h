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
#include "mdl.h"
#include "gpu/pppcdata.h"
#include "check.h"

template<class TILE,int WIDTH=32>
class MessagePPPC : public mdl::cudaMessage, public gpu::pppcData<TILE,WIDTH> {
protected:
    mdl::messageQueue<MessagePPPC> &freeQueue;
    bool bGravStep = false;
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
public:
    explicit MessagePPPC(mdl::messageQueue<MessagePPPC> &freeQueue)
        : freeQueue(freeQueue) {
        // For CUDA we want to "pin" the host memory for optimal performance
        CUDA_CHECK(cudaHostRegister,(this->pHostBufIn, this->requestBufferSize, cudaHostRegisterPortable));
        CUDA_CHECK(cudaHostRegister,(this->pHostBufOut,this->resultsBufferSize, cudaHostRegisterPortable));
    }
    virtual ~MessagePPPC() {
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufIn));
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufOut));
    }

};

typedef MessagePPPC<ilpTile> MessagePP;
typedef MessagePPPC<ilcTile> MessagePC;
#endif
#endif
