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

#include "pkd_config.h"
#include "metal.h"
#include "metal/pkdgrav3_metallib.h"

MetalClient::MetalClient(mdl::mdlClass &mdl) : mdl(mdl), pp(nullptr), pc(nullptr) {
    if (mdl.isMetalActive()) {
        freePP.enqueue(new MessagePP(freePP));
        freePP.enqueue(new MessagePP(freePP));
        freePP.enqueue(new MessagePP(freePP));
        freePP.enqueue(new MessagePP(freePP));
        freePC.enqueue(new MessagePC(freePC));
        freePC.enqueue(new MessagePC(freePC));
        freePC.enqueue(new MessagePC(freePC));
        freePC.enqueue(new MessagePC(freePC));
    }
}

void MetalClient::flushMETAL() {
    flush(pp);
    flush(pc);
}

template <class TILE> struct kernel_name;
template<> struct kernel_name<ilpTile> { const char *operator()() const {return "pp";}};
template<> struct kernel_name<ilcTile> { const char *operator()() const {return "pc";}};

template<class TILE,int interactionWidth>
void MessagePPPC<TILE,interactionWidth>::launch(mdl::metal::Stream &stream,MTL::CommandBuffer *cbuf) {
    auto &device = stream.getDevice();
    auto input  = device.getBuffer(this->pHostBufIn, this->requestBufferSize);
    auto output = device.getBuffer(this->pHostBufOut,this->resultsBufferSize);
    kernel_name<TILE> name;
    auto pipeline = stream.getPipeline(pkdgrav3_library,sizeof(pkdgrav3_library),name());
    typedef gpu::Blk<interactionWidth,TILE> BLK;

    uint nBlocksInteractions = 2;
    // First assume that we can fit all the interactions into shared memory
    uint nBlocksParticles = pipeline->maxTotalThreadsPerThreadgroup() / interactionWidth / nBlocksInteractions;
    while ((nBlocksParticles&1) == 0 && nBlocksParticles>8) { nBlocksInteractions *= 2; nBlocksParticles /= 2; }
    // If necessary, reduce the number of interactions so that we have enough shared memory.
    auto maxThreadgroupMemoryLength = device.maxThreadgroupMemoryLength() - pipeline->staticThreadgroupMemoryLength();
    uint nReduceMemory = 0;
    assert(pipeline->threadExecutionWidth()==interactionWidth);
    // auto nReduceMemory = pipeline->threadExecutionWidth()==interactionWidth
    //                      ? 0 : pipeline->maxTotalThreadsPerThreadgroup()/pipeline->threadExecutionWidth() * sizeof(gpu::ppResult);
    NS::UInteger nBlocksMemory;
    while ( (nBlocksMemory = nBlocksInteractions*sizeof(BLK)) + nReduceMemory > maxThreadgroupMemoryLength) {
        nBlocksInteractions /= 2;
        nBlocksParticles *= 2;
    }

    // Not all devices support nonuniform threadgroups, so we adjust the size of the grid.
    // This means additional work units will be processed, but they were already set to zero.
    assert(nBlocksInteractions <= this->wp_unit_mask+1);
    auto remainder = this->nTotalInteractionBlocks % nBlocksInteractions;
    auto nGridZ = this->nTotalInteractionBlocks + (remainder ? nBlocksInteractions-remainder : 0);
    MTL::Size threadsPerThreadgroup(interactionWidth,nBlocksParticles,nBlocksInteractions);
    MTL::Size threadsPerGrid(interactionWidth,nBlocksParticles,nGridZ);

    input->didModifyRange(NS::Range::Make(0,this->requestBufferCount));

    // Zero the results buffer
    auto zero = cbuf->blitCommandEncoder();
    zero->fillBuffer(output,NS::Range::Make(0,this->resultsBufferCount),0);
    zero->endEncoding();

    auto cenc = cbuf->computeCommandEncoder();
    cenc->setComputePipelineState(pipeline);
    size_t blkOffset  = 0;
    size_t partOffset = blkOffset  + sizeof(BLK) * this->nTotalInteractionBlocks;
    size_t wuOffset   = partOffset + sizeof(gpu::ppInput) * this->nTotalParticles;
    cenc->setBuffer(input,blkOffset,2);
    cenc->setBuffer(input,partOffset,1);
    cenc->setBuffer(input,wuOffset,0);
    cenc->setBuffer(output,0,3);
    cenc->setThreadgroupMemoryLength(nBlocksMemory,0);
    cenc->setThreadgroupMemoryLength(nReduceMemory,1);
    cenc->dispatchThreads(threadsPerGrid,threadsPerThreadgroup);
    cenc->endEncoding();

    auto blit = cbuf->blitCommandEncoder();
    // Doesn't work?: output->didModifyRange(NS::Range::Make(0,4));
    blit->synchronizeResource(output);
    blit->endEncoding();
}

void pkdParticleWorkDone(workParticle *wp);

template<class TILE,int interactionWidth>
void MessagePPPC<TILE,interactionWidth>::finish() {
    auto *pR = reinterpret_cast<gpu::ppResult *>(this->pHostBufOut);

    for ( auto &w : this->work ) {
        auto nP = w.wp->nP;
        auto *pInfoOut = w.wp->pInfoOut;
        for (auto ip=0; ip<nP; ++ip) {
            pInfoOut[ip].a[0]    += pR[ip].ax;
            pInfoOut[ip].a[1]    += pR[ip].ay;
            pInfoOut[ip].a[2]    += pR[ip].az;
            pInfoOut[ip].fPot    += pR[ip].fPot;
            pInfoOut[ip].dirsum  += pR[ip].dirsum;
            pInfoOut[ip].normsum += pR[ip].normsum;
        }
        pR += this->align_nP(nP);
        pkdParticleWorkDone(w.wp);
    }
    this->clear();
    freeQueue.enqueue(*this);
}

/*****************************************************************************\
*   Explicit instantiation for methods we need elsewhere
\*****************************************************************************/

template void MessagePP::finish();
template void MessagePC::finish();
template void MessagePP::launch(mdl::metal::Stream &stream,MTL::CommandBuffer *cbuf);
template void MessagePC::launch(mdl::metal::Stream &stream,MTL::CommandBuffer *cbuf);
