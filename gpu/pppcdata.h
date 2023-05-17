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
#include "workunit.h"
#include "gravity/ilp.h"
#include "gravity/ilc.h"

namespace gpu {

template <int n,class BTILE> struct Blk;
template <int n> struct Blk<n,ilpTile> : ppBlk<n> {};
template <int n> struct Blk<n,ilcTile> : pcBlk<n> {};

inline double getFlops(workParticle *wp,ilpTile &tile) {
    return COST_FLOP_PP*wp->nP*tile.count();
}
inline double getFlops(workParticle *wp,ilcTile &tile) {
    return COST_FLOP_PC*wp->nP*tile.count();
}

// The interaction list blocks may contain more information than is necessary
// for the GPU kernel, so we copy just the fields that we need here.
inline int copyBLKs(ppInteract *out, ilpTile &in) {
    auto n = in.width;
    auto nIlp = in.count();
    int i, nBlk = (nIlp+n-1) / n;
    for (i=0; i<nBlk; ++i) {
        memcpy(&out[i].dx,    &in[i].dx,    sizeof(out[i].dx));
        memcpy(&out[i].dy,    &in[i].dy,    sizeof(out[i].dy));
        memcpy(&out[i].dz,    &in[i].dz,    sizeof(out[i].dz));
        memcpy(&out[i].m,     &in[i].m,     sizeof(out[i].m));
        memcpy(&out[i].fourh2,&in[i].fourh2,sizeof(out[i].fourh2));
    }
    return nBlk;
}

inline int copyBLKs(pcInteract *out, ilcTile &in) {
    auto n = in.width;
    auto nIlp = in.count();
    int i, nBlk = (nIlp+n-1) / n;
    for (i=0; i<nBlk; ++i) {
        memcpy(&out[i].dx,  &in[i].dx,  sizeof(out[i].dx));
        memcpy(&out[i].dy,  &in[i].dy,  sizeof(out[i].dy));
        memcpy(&out[i].dz,  &in[i].dz,  sizeof(out[i].dz));
        memcpy(&out[i].xxxx,&in[i].xxxx,sizeof(out[i].xxxx));
        memcpy(&out[i].xxxy,&in[i].xxxy,sizeof(out[i].xxxy));
        memcpy(&out[i].xxxz,&in[i].xxxz,sizeof(out[i].xxxz));
        memcpy(&out[i].xxyz,&in[i].xxyz,sizeof(out[i].xxyz));
        memcpy(&out[i].xxyy,&in[i].xxyy,sizeof(out[i].xxyy));
        memcpy(&out[i].yyyz,&in[i].yyyz,sizeof(out[i].yyyz));
        memcpy(&out[i].xyyz,&in[i].xyyz,sizeof(out[i].xyyz));
        memcpy(&out[i].xyyy,&in[i].xyyy,sizeof(out[i].xyyy));
        memcpy(&out[i].yyyy,&in[i].yyyy,sizeof(out[i].yyyy));
        memcpy(&out[i].xxx, &in[i].xxx, sizeof(out[i].xxx));
        memcpy(&out[i].xyy, &in[i].xyy, sizeof(out[i].xyy));
        memcpy(&out[i].xxy, &in[i].xxy, sizeof(out[i].xxy));
        memcpy(&out[i].yyy, &in[i].yyy, sizeof(out[i].yyy));
        memcpy(&out[i].xxz, &in[i].xxz, sizeof(out[i].xxz));
        memcpy(&out[i].yyz, &in[i].yyz, sizeof(out[i].yyz));
        memcpy(&out[i].xyz, &in[i].xyz, sizeof(out[i].xyz));
        memcpy(&out[i].xx,  &in[i].xx,  sizeof(out[i].xx));
        memcpy(&out[i].xy,  &in[i].xy,  sizeof(out[i].xy));
        memcpy(&out[i].xz,  &in[i].xz,  sizeof(out[i].xz));
        memcpy(&out[i].yy,  &in[i].yy,  sizeof(out[i].yy));
        memcpy(&out[i].yz,  &in[i].yz,  sizeof(out[i].yz));
#ifdef USE_DIAPOLE
        memcpy(&out[i].x,   &in[i].x,   sizeof(out[i].x));
        memcpy(&out[i].y,   &in[i].y,   sizeof(out[i].y));
        memcpy(&out[i].z,   &in[i].z,   sizeof(out[i].z));
#endif
        memcpy(&out[i].m,   &in[i].m,   sizeof(out[i].m));
        memcpy(&out[i].u,   &in[i].u,   sizeof(out[i].u));
    }
    return nBlk;
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
    static constexpr size_t wp_unit_mask = 7; // Pad input buffer to this mask (so multiples of 8)
    std::vector<workInformation> work; // [wp_max_buffer]

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
        nTotalInteractionBlocks += copyBLKs(blk+nTotalInteractionBlocks,tile); // (the ILP tile can now be freed/reused)
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
        // Pad the work unit
        for (auto i=nTotalInteractionBlocks; i&wp_unit_mask; ++i) {
            wuHost->nP = 0;
            wuHost->iP = 0;
            wuHost->nI = 0;
            ++wuHost;
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
        const auto nBlocks = (nI+WIDTH-1) / WIDTH; // number of interaction blocks needed
        return (nBlocks + nTotalInteractionBlocks) * (sizeof(BLK) + sizeof(ppWorkUnit))
               + (align_nP(nP) + nTotalParticles) * sizeof(ppInput)
               + wp_unit_mask * sizeof(ppWorkUnit);
    }

    template<class BLK>
    int outputSize(int nP=0, int nI=0) {
        //const auto nBlocks = (nI+BLK::width-1) / BLK::width; // number of interaction blocks needed
        return (align_nP(nP) + nTotalParticles) * sizeof(ppResult);
    }
};

} // namespace gpu

#endif
