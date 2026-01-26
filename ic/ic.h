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

#ifndef IC_H
#define IC_H
#include "TraversePST.h"

#include <cstddef>
#include <array>
#include <stdint.h>
#include <queue>
#include <vector>
#include "blitz/array.h"

struct basicParticle {
    blitz::TinyVector<float, 3> dr;
    blitz::TinyVector<float, 3> v;
    basicParticle() = default;
};

struct expandParticle {
    uint64_t ix : 21;
    uint64_t iy : 21;
    uint64_t iz : 21;
    blitz::TinyVector<float, 3> dr;
    blitz::TinyVector<float, 3> v;
    expandParticle() = default;
};

struct integerParticle {
    blitz::TinyVector<int32_t, 3> r;
    blitz::TinyVector<float, 3> v;
    integerParticle() = default;
};

typedef union {
    basicParticle b;
    expandParticle e;
    integerParticle i;
} overlayedParticle;

#ifdef MDL_FFTW

typedef struct {
    FFTW3(real) x, y, z;
} gridpos;

typedef struct {
    FFTW3(real) x, y, z;
    FFTW3(real) vx, vy, vz;
} gridpsc;

struct gridInfoLPT {
    int nGrids = 0;
    int indexPhi1 = -1; // 1LPT potential (persistent)
    int indexPhi2 = -1; // 2LPT potential (persistent)
    int indexPhi3 = -1; // 3LPT potentials
    int indexTmp0 = -1; // temporary grid
    int indexTmp1 = -1; // temporary grid

    /* Function for determining the total number of grids needed as well as
    ** the grid index mapping used during LPT computations.
    */
    gridInfoLPT(int iLPT) {
        bool bExtraTmpGridFor1LPT = true; // allow one additional temporary grid when doing 1LPT only?
        nGrids = 6;                       // particle positions and velocities
        if (iLPT >= 1) {
            /* For 1LPT we need 1 potential grid and one temporary grid */
            indexPhi1 = nGrids++;
            indexTmp0 = nGrids++;
            /* We can save time by having one additional temporary grid */
            indexTmp1 = indexTmp0;
            if (iLPT == 1 && bExtraTmpGridFor1LPT) {
                indexTmp1 = nGrids++;
            }
        }
        if (iLPT >= 2) {
            /* For 2LPT we need 1 additional potential grid
            ** and one additional temporary grid.
            */
            indexPhi2 = nGrids++;
            indexTmp1 = nGrids++;
        }
        if (iLPT >= 3) {
            /* For 3LPT we need 1 additional potential grid */
            indexPhi3 = nGrids++;
        }
        // It is better to adjust MAX_LPT_GRIDS above if it changes, but everything
        // should work with larger values except it will allocate on the heap.
        // assert(gridInfo.nGrids <= MAX_LPT_GRIDS);
    }
};
#endif

/* PLT_MOVEIC */
struct inMoveIC {
    overlayedParticle* pBase;
    uint64_t iStart;
    uint64_t nMove;
    float fMass;
    float fSoft;
    int nGrid;
    int bICgas;
    int nBucket;
    double dInitialT;
    double dInitialH;
    double dInitialHe;
    double dInitialC;
    double dInitialN;
    double dInitialO;
    double dInitialNe;
    double dInitialMg;
    double dInitialSi;
    double dInitialFe;
    double dInitialMetallicity;
    double dExpansion;
    double dBaryonFraction;
    double dTuFac;
};
int pltMoveIC(PST, void*, int, void*, int);
int pstMoveIC(PST, void*, int, void*, int);

class ServiceGenerateIC : public TraverseCount<ServiceGenerateIC> {
  public:
    static constexpr std::string_view name = "GenerateIC";
    static constexpr size_t max_tf = 4096;
    struct input {
        uint64_t nPerNode;
        double dBoxSize;
        double dBoxMass;
        double dBaryonFraction;
        double dExpansion;
        int iSeed;
        int bFixed;
        float fPhase;
        int nGrid;
        int iLPT;
        int bICgas;
        int nBucket;
        double dInitialT;
        double dInitialH;
        double dInitialHe;
        double dInitialC;
        double dInitialN;
        double dInitialO;
        double dInitialNe;
        double dInitialMg;
        double dInitialSi;
        double dInitialFe;
        double dInitialMetallicity;
        double dTuFac;
        int bComove;
        int nTf;
        std::array<double, max_tf> k;
        std::array<double, max_tf> tf;
    };
    struct output {
        uint64_t N;
        double noiseMean;
        double noiseCSQ;
        auto operator+=(const output& o) {
            N += o.N;
            noiseMean += o.noiseMean;
            noiseCSQ += o.noiseCSQ;
            return *this;
        }
    };
    explicit ServiceGenerateIC(PST pst) : TraverseCount(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
    auto AtNode(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

#endif
