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
#ifndef PPPCDATA_H
#define PPPCDATA_H
#include "hostdata.h"
#include "basetype.h"

namespace gpu {

template <int n,class BTILE> struct Blk;
template <int n> struct Blk<n,ilpTile> : BlockPP<n> {};
template <int n> struct Blk<n,ilcTile> : BlockPC<n> {};

// One of these entries for each interaction block
struct alignas(8) ppWorkUnit {
    std::uint32_t iP;   // Index of first particle
    std::uint16_t nP;   // Number of particles
    std::uint16_t nI;   // Number of interactions in the block
};
static_assert(sizeof(ppWorkUnit)==8);

struct alignas(32) ppInput {
    float dx, dy, dz;
    float ax, ay, az;
    float fSoft2;
    float dImaga;
};
static_assert(sizeof(ppInput)==32);

/* Each thread block outputs this for each particle */
struct alignas(32) ppResult {
    float ax;
    float ay;
    float az;
    float fPot;
    float dirsum;
    float normsum;
};
static_assert(sizeof(ppResult)==32);

inline double getFlops(workParticle *wp,ilpTile &tile) {
    return COST_FLOP_PP*wp->nP*tile.count();
}
inline double getFlops(workParticle *wp,ilcTile &tile) {
    return COST_FLOP_PC*wp->nP*tile.count();
}

/// Extend hostData to keeps track of interaction lists and work packages.
/// We are still GPU agnostic.
template<class TILE,int WIDTH=32>
class pppcData : public hostData {
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
    std::vector<workInformation> work; // [wp_max_buffer]

    template<int n,class BTILE>
    int copyBLKs2(Blk<n,BTILE> *out, BTILE &in) {
        assert(n==ILP_PART_PER_BLK);
        auto nIlp = in.count();
        int i, nBlk = (nIlp+n-1) / n;
        for (i=0; i<nBlk; ++i) memcpy(&out[i],&in[i],sizeof(out[i]));
        return nBlk;
    }

public:
    pppcData() : requestBufferCount(0), resultsBufferCount(0), nTotalInteractionBlocks(0), nTotalParticles(0) {
        work.reserve(wp_max_buffer);
    }

    void clear() {
        work.clear();
        requestBufferCount = resultsBufferCount = 0;
        nTotalInteractionBlocks = nTotalParticles = 0;
    }
    bool queue(workParticle *wp, TILE &tile, bool bGravStep) {
        typedef Blk<WIDTH,TILE> BLK;
        if (work.size() == wp_max_buffer) return false;  // Too many work packages so send the work
        this->bGravStep = bGravStep;
        const auto nP = wp->nP;                                 // Queue this many particles
        const auto nI = tile.count();                           // ... operating on this many interactions

        if (inputSize<BLK>(nP,nI) > requestBufferSize - 1024) return false;     // Refuse if this tile won't fit in this buffer
        if (outputSize<BLK>(nP,nI) > requestBufferSize) return false;           // Refuse if this response won't fit
        auto blk = reinterpret_cast<Blk<WIDTH,TILE> *>(pHostBufIn);             // Copy the blocks to the input buffer
        nTotalInteractionBlocks += copyBLKs2(blk+nTotalInteractionBlocks,tile); // (the ILP tile can now be freed/reused)
        nTotalParticles += align_nP(nP);

        ++wp->nRefs;
        work.emplace_back(wp,tile.count());
        wp->dFlopSingleGPU += getFlops(wp,tile);
        return true;
    }
    void prepare() {
        typedef Blk<WIDTH,TILE> BLK;
        // The interation blocks -- already copied to the host memory
        auto *__restrict__ blkHost = reinterpret_cast<BLK *>(pHostBufIn);
        // The particle information
        auto *__restrict__ partHost = reinterpret_cast<ppInput *>(blkHost + nTotalInteractionBlocks);
        // The interaction block descriptors
        auto *__restrict__ wuHost = reinterpret_cast<ppWorkUnit *>(partHost + nTotalParticles);
        uint32_t iP = 0;
        for ( auto &w : work ) {
            int nI = w.nInteractions;
            auto nP = w.wp->nP;
            auto *pInfoIn = w.wp->pInfoIn;
            auto nBlocks = (nI + BLK::width - 1) / BLK::width;

            // Generate a interaction block descriptor for each block
            for (auto j=0; j<nBlocks; ++j) {
                wuHost->nP = nP;
                wuHost->iP = iP;
                wuHost->nI = nI > BLK::width ? BLK::width : nI;
                nI -= wuHost->nI;
                ++wuHost;
            }

            // Copy in nP particles
            for (auto j=0; j<nP; ++j) {
                partHost[j].dx =  pInfoIn[j].r[0];
                partHost[j].dy =  pInfoIn[j].r[1];
                partHost[j].dz =  pInfoIn[j].r[2];
                partHost[j].ax =  pInfoIn[j].a[0];
                partHost[j].ay =  pInfoIn[j].a[1];
                partHost[j].az =  pInfoIn[j].a[2];
                partHost[j].fSoft2 = pInfoIn[j].fSmooth2;
                /*partHost[j].dImaga = 0;*/
            }
            nP = align_nP(nP);
            partHost += nP;
            iP += nP;
        }
        requestBufferCount = reinterpret_cast<char *>(wuHost) - reinterpret_cast<char *>(pHostBufIn);
        resultsBufferCount = iP * sizeof(ppResult);
        nGrid = nTotalInteractionBlocks;
        assert(requestBufferCount <= requestBufferSize);
        assert(resultsBufferCount <= resultsBufferSize);
    }

    auto align_nP(int nP) {
        constexpr int mask = 32 * sizeof(float) / sizeof(ppInput) - 1;
        static_assert(32 * sizeof(float) == (mask+1)*sizeof(ppInput));
        return (nP + mask) & ~mask;
    }

    template<class BLK>
    int inputSize(int nP=0, int nI=0) {
        const auto nBlocks = (nI+BLK::width-1) / BLK::width; // number of interaction blocks needed
        return (nBlocks + nTotalInteractionBlocks) * (sizeof(BLK) + sizeof(ppWorkUnit)) + (align_nP(nP) + nTotalParticles) * sizeof(ppInput);
    }

    template<class BLK>
    int outputSize(int nP=0, int nI=0) {
        //const auto nBlocks = (nI+BLK::width-1) / BLK::width; // number of interaction blocks needed
        return (align_nP(nP) + nTotalParticles) * sizeof(ppResult);
    }
};

} // namespace gpu

#endif
