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

#include <stdint.h>
#include <queue>
#include <vector>
#include "blitz/array.h"

typedef struct {
    blitz::TinyVector<float,3> dr;
    blitz::TinyVector<float,3> v;
} basicParticle;

typedef struct {
    uint64_t ix : 21;
    uint64_t iy : 21;
    uint64_t iz : 21;
    blitz::TinyVector<float,3> dr;
    blitz::TinyVector<float,3> v;
} expandParticle;

typedef struct {
    blitz::TinyVector<int32_t,3> r;
    blitz::TinyVector<float,3> v;
} integerParticle;

typedef union {
    basicParticle b;
    expandParticle e;
    integerParticle i;
} overlayedParticle;

#ifdef MDL_FFTW

typedef struct {
    FFTW3(real) x,y,z;
} gridpos;

typedef struct {
    FFTW3(real) x,y,z;
    FFTW3(real) vx,vy,vz;
} gridpsc;

typedef struct {
    int nGrids = 0;
    int indexPhi1 = -1;  // 1LPT potential (persistent)
    int indexPhi2 = -1;  // 2LPT potential (persistent)
    int indexPhi3 = -1;  // 3LPT potentials
    int indexTmp0 = -1;  // temporary grid
    int indexTmp1 = -1;  // temporary grid
} gridInfoLPT;

gridInfoLPT getGridInfoLPT(int iLPT);

int pkdGenerateIC(PKD pkd, MDLFFT fft, int iSeed, int bFixed, float fPhase, int nGrid, int iLPT, double dBoxSize,
                  double a, int nTf, double *tk, double *tf, double *noiseMean, double *noiseCSQ);
#endif
#endif
