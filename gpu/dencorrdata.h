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
#ifndef DENCORRDATA_H
#define DENCORRDATA_H
#include "hostdata.h"
#include "basetype.h"
#include "workunit.h"
#include "gravity/ilp.h"

namespace gpu {

inline double getFlopsDenCorr(workParticle *wp,ilpTile &tile) {
    return 0.0*wp->nP*tile.count();
}

inline int copyBLKs(denCorrInteract *out, ilpTile &in) {
    auto n = in.width;
    auto nIlp = in.count();
    int i, nBlk = (nIlp+n-1) / n;
    for (i=0; i<nBlk; ++i) {
        memcpy(&out[i].dx,      &in[i].dx,      sizeof(out[i].dx));
        memcpy(&out[i].dy,      &in[i].dy,      sizeof(out[i].dy));
        memcpy(&out[i].dz,      &in[i].dz,      sizeof(out[i].dz));
        memcpy(&out[i].T,       &in[i].T,       sizeof(out[i].T));
        memcpy(&out[i].P,       &in[i].P,       sizeof(out[i].P));
        memcpy(&out[i].expImb2, &in[i].expImb2, sizeof(out[i].expImb2));
        memcpy(&out[i].isGas,   &in[i].isGas,   sizeof(out[i].isGas));
    }
    return nBlk;
}

/// Extend hostData to keeps track of interaction lists and work packages.
/// We are still GPU agnostic.
template<int WIDTH=32>
class denCorrData : public hostData {
protected:
    bool bGravStep;
    std::size_t requestBufferCount, resultsBufferCount;
    int nTotalInteractionBlocks, nTotalParticles, nGrid;
    struct workInformation {
        workParticle *wp;
        std::size_t nInteractions;
        workInformation(workParticle *wp, size_t nInteractions) : wp(wp), nInteractions(nInteractions) {}
    };
    static constexpr size_t wp_max_buffer = 128;
    static constexpr size_t wp_unit_mask = 7; // Pad input buffer to this mask (so multiples of 8)
    std::vector<workInformation> work; // [wp_max_buffer]

public:
    denCorrData() : requestBufferCount(0), resultsBufferCount(0), nTotalInteractionBlocks(0), nTotalParticles(0) {
        work.reserve(wp_max_buffer);
    }

    void clear() {
        work.clear();
        requestBufferCount = resultsBufferCount = 0;
        nTotalInteractionBlocks = nTotalParticles = 0;
    }
    bool queue(workParticle *wp, ilpTile &tile, bool bGravStep) {
        typedef denCorrBlk<WIDTH> BLK;
        if (work.size() == wp_max_buffer) return false;  // Too many work packages so send the work
        this->bGravStep = bGravStep;
        const auto nP = wp->nP;                                 // Queue this many particles
        const auto nI = tile.count();                           // ... operating on this many interactions

        if (inputSize<BLK>(nP,nI) > requestBufferSize - 1024) return false;     // Refuse if this tile won't fit in this buffer
        if (outputSize<BLK>(nP,nI) > requestBufferSize) return false;           // Refuse if this response won't fit
        auto blk = reinterpret_cast<denCorrBlk<WIDTH> *>(pHostBufIn);             // Copy the blocks to the input buffer
        nTotalInteractionBlocks += copyBLKs(blk+nTotalInteractionBlocks,tile); // (the ILP tile can now be freed/reused)
        nTotalParticles += align_nP(nP);

        ++wp->nRefs;
        work.emplace_back(wp,tile.count());
        wp->dFlopSingleGPU += getFlopsDenCorr(wp,tile);
        return true;
    }
    void prepare() {
        typedef denCorrBlk<WIDTH> BLK;
        // The interation blocks -- already copied to the host memory
        auto *__restrict__ blkHost = reinterpret_cast<BLK *>(pHostBufIn);
        // The particle information
        auto *__restrict__ partHost = reinterpret_cast<denCorrInput *>(blkHost + nTotalInteractionBlocks);
        // The interaction block descriptors
        auto *__restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(partHost + nTotalParticles);
        uint32_t iP = 0;
        for ( auto &w : work ) {
            int nI = w.nInteractions;
            auto nP = w.wp->nP;
            auto *pInfoIn = w.wp->pInfoIn;
            auto nBlocks = (nI + WIDTH - 1) / WIDTH;

            // Generate a interaction block descriptor for each block
            for (auto j=0; j<nBlocks; ++j) {
                wuHost->nP = nP;
                wuHost->iP = iP;
                wuHost->nI = nI > WIDTH ? WIDTH : nI;
                nI -= wuHost->nI;
                ++wuHost;
            }

            // Copy in nP particles
            for (auto j=0; j<nP; ++j) {
                partHost[j].dx    = pInfoIn[j].r[0];
                partHost[j].dy    = pInfoIn[j].r[1];
                partHost[j].dz    = pInfoIn[j].r[2];
                partHost[j].fBall = pInfoIn[j].fBall;
                partHost[j].isGas = pInfoIn[j].isGas;
            }
            nP = align_nP(nP);
            partHost += nP;
            iP += nP;
        }
        // Pad the work unit
        for (auto i=nTotalInteractionBlocks; i&wp_unit_mask; ++i) {
            wuHost->nP = 0;
            wuHost->iP = 0;
            wuHost->nI = 0;
            ++wuHost;
        }
        requestBufferCount = reinterpret_cast<char *>(wuHost) - reinterpret_cast<char *>(pHostBufIn);
        resultsBufferCount = iP * sizeof(denCorrResult);
        nGrid = nTotalInteractionBlocks;
        assert(requestBufferCount <= requestBufferSize);
        assert(resultsBufferCount <= resultsBufferSize);
    }

    auto align_nP(int nP) {
        constexpr int mask = 32 * sizeof(float) / sizeof(denCorrInput) - 1;
        static_assert(32 * sizeof(float) == (mask+1)*sizeof(denCorrInput));
        return (nP + mask) & ~mask;
    }

    template<class BLK>
    int inputSize(int nP=0, int nI=0) {
        const auto nBlocks = (nI+WIDTH-1) / WIDTH; // number of interaction blocks needed
        return (nBlocks + nTotalInteractionBlocks) * (sizeof(BLK) + sizeof(ppWorkUnit))
               + (align_nP(nP) + nTotalParticles) * sizeof(denCorrInput)
               + wp_unit_mask * sizeof(ppWorkUnit);
    }

    template<class BLK>
    int outputSize(int nP=0, int nI=0) {
        //const auto nBlocks = (nI+BLK::width-1) / BLK::width; // number of interaction blocks needed
        return (align_nP(nP) + nTotalParticles) * sizeof(denCorrResult);
    }
};

} // namespace gpu

#endif
