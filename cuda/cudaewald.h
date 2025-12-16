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

#ifndef CUDAEWALD_H
#define CUDAEWALD_H
#ifdef USE_CUDA
#include "mdlcuda.h"
#include "gpu/hostdata.h"


template <typename WORKPARTICLE>
class MessageEwaldBase : public mdl::cudaMessage,public gpu::hostData {
protected:
    class CudaClient &cuda;
    std::vector<WORKPARTICLE *> ppWP; // [CUDA_WP_MAX_BUFFERED]
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
    virtual void finish() override;
    int nParticles, nMaxParticles;
public:
    explicit MessageEwaldBase(class CudaClient &cuda);
    bool queue(WORKPARTICLE *work); // Add work to this message; return false if this is not possible
};

class MessageEwald : public MessageEwaldBase<workParticle> {
public:
    explicit MessageEwald(class CudaClient &cuda) : MessageEwaldBase(cuda) {};
};

class MessageEwaldnew : public MessageEwaldBase<GravityInteractionTraits::WORKPARTICLE> {
public:
    explicit MessageEwaldnew(class CudaClient &cuda) : MessageEwaldBase(cuda) {};
};
// typedef MessageEwaldBase<workParticle> MessageEwald ;
// typedef MessageEwaldBase<GravityInteractionTraits::WORKPARTICLE> MessageEwaldnew ;


class MessageEwaldSetup : public mdl::cudaMessage {
protected:
    virtual void launch(mdl::Stream &stream,void *pCudaBufIn, void *pCudaBufOut) override;
public:
    explicit MessageEwaldSetup(struct EwaldVariables *const ew, EwaldTable *const ewt,int iDevice=-1);
protected:
    struct EwaldVariables *const ewIn;
    EwaldTable *const ewt;
    std::vector<momFloat> dLx, dLy, dLz;
    std::vector<int> ibHole;
};
#endif
#endif
