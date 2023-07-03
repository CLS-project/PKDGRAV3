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
#include "mdlcuda.h"
#include "gpu/pppcdata.h"
#include "gpu/dendata.h"
#include "gpu/dencorrdata.h"
#include "gpu/sphforcedata.h"
#include "check.h"
#include <queue>

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

template<int WIDTH=32>
class MessageDen : public mdl::cudaMessage, public gpu::denData<WIDTH> {
protected:
    mdl::messageQueue<MessageDen> &freeQueue;
    bool bGravStep = false;
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
public:
    std::queue<workParticle *> *wps;
    explicit MessageDen(mdl::messageQueue<MessageDen> &freeQueue, std::queue<workParticle *> *wps)
        : freeQueue(freeQueue), wps(wps) {
        // For CUDA we want to "pin" the host memory for optimal performance
        CUDA_CHECK(cudaHostRegister,(this->pHostBufIn, this->requestBufferSize, cudaHostRegisterPortable));
        CUDA_CHECK(cudaHostRegister,(this->pHostBufOut,this->resultsBufferSize, cudaHostRegisterPortable));
    }
    virtual ~MessageDen() {
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufIn));
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufOut));
    }

};

template<int WIDTH=32>
class MessageDenCorr : public mdl::cudaMessage, public gpu::denCorrData<WIDTH> {
protected:
    mdl::messageQueue<MessageDenCorr> &freeQueue;
    bool bGravStep = false;
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
public:
    explicit MessageDenCorr(mdl::messageQueue<MessageDenCorr> &freeQueue)
        : freeQueue(freeQueue) {
        // For CUDA we want to "pin" the host memory for optimal performance
        CUDA_CHECK(cudaHostRegister,(this->pHostBufIn, this->requestBufferSize, cudaHostRegisterPortable));
        CUDA_CHECK(cudaHostRegister,(this->pHostBufOut,this->resultsBufferSize, cudaHostRegisterPortable));
    }
    virtual ~MessageDenCorr() {
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufIn));
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufOut));
    }

};

template<int WIDTH=32>
class MessageSPHForce : public mdl::cudaMessage, public gpu::sphForceData<WIDTH> {
protected:
    mdl::messageQueue<MessageSPHForce> &freeQueue;
    bool bGravStep = false;
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
public:
    explicit MessageSPHForce(mdl::messageQueue<MessageSPHForce> &freeQueue)
        : freeQueue(freeQueue) {
        // For CUDA we want to "pin" the host memory for optimal performance
        CUDA_CHECK(cudaHostRegister,(this->pHostBufIn, this->requestBufferSize, cudaHostRegisterPortable));
        CUDA_CHECK(cudaHostRegister,(this->pHostBufOut,this->resultsBufferSize, cudaHostRegisterPortable));
    }
    virtual ~MessageSPHForce() {
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufIn));
        CUDA_CHECK(cudaHostUnregister,(this->pHostBufOut));
    }

};
#endif
#endif
